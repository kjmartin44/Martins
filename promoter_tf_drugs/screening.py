"""Screen the whole approved-drug library against a receptor.

Instead of grading the handful of drugs already annotated against your TF, this
scores *every approved drug* with a structure, using the per-receptor
classifier from agonism.py, and ranks them.

About "effects not seen before"
-------------------------------
Be precise about what is and is not possible here, because the distinction
decides whether the output is useful or misleading.

**Not possible:** discovering a mechanism nobody has characterised. Every model
in this package is trained on published measurements. A classifier cannot know
about a binding mode that has never been observed. No amount of compute changes
this; it is a property of the training data, not of the method.

**Possible and genuinely useful:** finding drug/receptor *pairs* nobody has
tested. There are roughly 4,000 approved drugs and each has been assayed
against a tiny fraction of the proteome. A compound approved in 1974 for an
unrelated indication may well engage your receptor; if nobody ran the assay,
that pair is absent from every database. Surfacing those is real repurposing,
and it is what this module is for.

The catch is that these two live in tension. A prediction is reliable in
proportion to how similar the query is to the training set — and similar
molecules are exactly the ones already tested. So the interesting quadrant
(novel *and* predicted active) is also the least reliable one. Hiding that
tradeoff behind a single score would be dishonest, so every hit is placed in
one of four quadrants and labelled:

    known ligand          already annotated at this receptor. Not a discovery,
                          but a positive control: if your model ranks these
                          highly, it is working.
    confident prediction  close to known ligands, predicted active. Most likely
                          to replicate. Least likely to be new.
    NOVEL HYPOTHESIS      dissimilar to known ligands, predicted active. This is
                          the repurposing quadrant. It is an *extrapolation* —
                          the model is being asked about chemistry it was not
                          trained on. Treat every one as a hypothesis requiring
                          a binding assay, never as a finding.
    no call               too far outside the applicability domain to score.

The ranking within NOVEL HYPOTHESIS is not a probability of being right. Expect
most to be wrong; that is normal for virtual screening, where single-digit
percent hit rates are considered a success.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

import numpy as np

from .http import Client, HttpError

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

NOVELTY_CUTOFF = 0.55      # max Tanimoto below this = structurally novel here
CONFIDENT_P = 0.70         # predicted-active threshold
PROMISCUITY_FLAG = 0.25    # share of queries a compound may connect to
MIN_ACTIVES = 20           # below this a binding model is not worth training
MIN_BINDING_AUC = 0.70     # below this the model is not discriminating


@dataclass
class ScreenHit:
    drug: str
    chembl_id: str
    # The model's score. For a molecule the model was TRAINED on this is an
    # in-sample fit, so `probability_in_sample` says which it is. It is kept
    # rather than withheld because the enrichment check needs it: whether
    # known actives outscore novel candidates is the diagnostic that catches
    # an inverted ranking. It is not a number to rank drugs by.
    probability: float | None
    nearest_tanimoto: float
    # None for a known ligand. A training molecule is its own nearest
    # neighbour, so 1 - tanimoto is structurally 0.0 for every one of them -
    # an artifact of the lookup, not an assessment of novelty.
    novelty: float | None
    quadrant: str
    known_ligand: bool
    max_phase: float | None
    # The measured potency that made it a known ligand. This is the real
    # number for an established binder, and what its row is ranked by.
    measured_pchembl: float | None = None
    probability_in_sample: bool = False
    # What the molecule is, and what it is approved to treat. This table
    # is where a reader first meets a drug name, and a quadrant plus a
    # probability does not say whether a topical acne gel and an infused
    # cytotoxic are both on the list - which, for the receptors that
    # qualify for a screen, they usually are.
    indication_hint: str | None = None
    treats: str | None = None
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------ drug library
def approved_library(client: Client, limit: int = 4000,
                     page: int = 200) -> list[dict]:
    """All approved (max_phase 4) small molecules with a structure."""
    out: list[dict] = []
    offset = 0
    while len(out) < limit:
        try:
            res = client.get(f"{CHEMBL}/molecule", {
                "max_phase": 4, "molecule_type": "Small molecule",
                "limit": page, "offset": offset, "format": "json",
            })
        except HttpError as exc:
            log.warning("library fetch stopped at offset %d: %s", offset, exc)
            break
        mols = res.get("molecules", [])
        if not mols:
            break
        for m in mols:
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if not smi:
                continue
            out.append({
                "chembl_id": m["molecule_chembl_id"],
                "name": (m.get("pref_name") or m["molecule_chembl_id"]).title(),
                "smiles": smi,
                "max_phase": m.get("max_phase"),
            })
        offset += page
        if not (res.get("page_meta") or {}).get("next"):
            break
    log.info("approved library: %d molecules", len(out))
    return out


# ------------------------------------------------------------------ screen
def _quadrant(prob: float, tanimoto: float, in_domain: bool,
              known_binder: bool, linked_to_gene: bool = False) -> str:
    """Two independent axes, previously collapsed into one.

        drug <-> receptor : established binder, or predicted
        drug <-> gene     : already reported to affect it, or untested

    Calling an established binder "a positive control, not a finding" was
    wrong. A drug that reliably binds the receptor and has never been tested
    against YOUR gene is the lowest-risk candidate available: the
    drug-target step is already de-risked and only the target-gene step is a
    hypothesis. A structurally novel compound is a hypothesis on both counts.
    """
    if known_binder:
        return ("confirms a reported link" if linked_to_gene
                else "ESTABLISHED LIGAND, untested on this gene")
    if not in_domain:
        return "no call"
    if prob < CONFIDENT_P:
        return "predicted inactive"
    return ("NOVEL HYPOTHESIS" if tanimoto < NOVELTY_CUTOFF
            else "predicted binder, close to known chemistry")


def screen(model, library: list[dict], *, goal: str = "up",
           known_ids: set | None = None, domain_cutoff: float = 0.35,
           top_n: int = 50, gene_linked: set | None = None,
           known_potency: dict | None = None) -> list[ScreenHit]:
    """Score a drug library against a trained per-receptor model.

    `known_potency` maps a known ligand's ChEMBL id to its measured pChEMBL.
    Supplying it is what lets the established-ligand rows be ranked by
    measured affinity instead of by a model score that was fitted on those
    same molecules.
    """
    from . import agonism
    from . import indications as _ind

    if model is None:
        raise ValueError(
            "no trained model for this receptor - too few labelled ligands. "
            "Screening without a model would produce arbitrary rankings.")

    known_ids = known_ids or set()
    smiles = [d["smiles"] for d in library]
    X, bvs, keep = agonism.fingerprints(smiles)
    if len(keep) == 0:
        return []

    probs = model.clf.predict_proba(X)[:, 1]
    hits: list[ScreenHit] = []

    for row, idx in enumerate(keep):
        drug = library[idx]
        p_ag = float(probs[row])
        p_want = p_ag if goal == "up" else 1.0 - p_ag
        tan = model.nearest_tanimoto(bvs[row])
        in_domain = tan >= domain_cutoff
        known = drug["chembl_id"] in known_ids
        linked = (drug["name"] or "").strip().lower() in (gene_linked or set())
        quad = _quadrant(p_want, tan, in_domain, known, linked)
        potency = (known_potency or {}).get(drug["chembl_id"]) if known else None
        desc = _ind.describe(drug["name"])

        if quad == "predicted inactive":
            continue

        notes = []
        if quad == "NOVEL HYPOTHESIS":
            notes.append(
                f"extrapolation: nearest known ligand Tanimoto {tan:.2f}. "
                "The model has not seen chemistry like this at this receptor. "
                "Requires a binding assay before it means anything.")
        elif quad == "ESTABLISHED LIGAND, untested on this gene":
            notes.append(
                "binds this receptor with measured activity"
                + (f" (pChEMBL {potency:.1f})" if potency is not None else "")
                + ", and has no reported effect on your gene. The "
                "pharmacology is already de-risked; only the "
                "receptor-to-gene step is untested, which makes this the "
                "cheapest experiment on the list.")
            notes.append(
                "the model score on this row is an in-sample fit - the "
                "molecule is one the model was trained on - so the row is "
                "ranked by its measured potency instead")
        elif quad == "confirms a reported link":
            notes.append("binds this receptor AND is already reported to "
                         "affect your gene — a mechanism for a known effect, "
                         "not a new lead")
        elif quad == "no call":
            notes.append("outside applicability domain; score withheld")

        hits.append(ScreenHit(
            drug=drug["name"], chembl_id=drug["chembl_id"],
            probability=round(p_want, 3),
            probability_in_sample=known,
            nearest_tanimoto=round(tan, 3),
            novelty=(None if known else round(1.0 - tan, 3)),
            quadrant=quad, known_ligand=known,
            max_phase=drug.get("max_phase"),
            measured_pchembl=(round(potency, 2) if potency is not None
                              else None),
            indication_hint=desc["what_it_is"], treats=desc["treats"],
            notes=notes,
        ))

    # Established binders that are untested on the gene come first: they are
    # the most actionable, not a control to be skipped past.
    order = {"ESTABLISHED LIGAND, untested on this gene": 0,
             "confirms a reported link": 1,
             "predicted binder, close to known chemistry": 2,
             "NOVEL HYPOTHESIS": 3, "no call": 4}

    def rank(h):
        # Within a quadrant, rank by whichever number is real for it:
        # measured potency for a known ligand, model score for a prediction.
        # Sorting both by the model score put the best-memorised training
        # molecule on top, which is not a property of the drug.
        if h.measured_pchembl is not None:
            return (order[h.quadrant], -h.measured_pchembl)
        if h.probability is not None:
            return (order[h.quadrant], -h.probability)
        return (order[h.quadrant], 0.0)

    hits.sort(key=rank)
    return hits[:top_n]


def build_binding_model(actives: list[dict], decoys: list[dict], *,
                        seed: int = 0):
    """Train "binds this receptor or not" — the actual screening question.

    Separating this from the agonist/antagonist model fixes two problems at
    once. Binding labels come from /activity and are abundant, while direction
    labels come from /mechanism and are scarce, so a screen that refused to run
    for want of 8 direction labels can train on hundreds of actives. And the
    quality gate becomes meaningful: a model that cannot tell actives from
    decoys is broken, whereas a model whose known ligands score low on
    *agonism* may be entirely correct.

    Decoys are sampled from the approved library. Some will be undiscovered
    binders, which depresses the measured AUC — a conservative bias, so a model
    that passes this gate passes a hard test.
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score
    from . import agonism

    if len(actives) < MIN_ACTIVES or not decoys:
        return None, None

    smiles = [a["smiles"] for a in actives] + [d["smiles"] for d in decoys]
    X, bvs, keep = agonism.fingerprints(smiles)
    if len(keep) < MIN_ACTIVES:
        return None, None

    n_act = len(actives)
    y = np.array([1 if i < n_act else 0 for i in keep])
    if y.sum() < MIN_ACTIVES or (1 - y).sum() < MIN_ACTIVES:
        return None, None

    clf = RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                 class_weight="balanced", random_state=seed,
                                 n_jobs=-1)
    auc = None
    try:
        folds = min(5, int(y.sum()), int((1 - y).sum()))
        cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
        proba = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
        auc = float(roc_auc_score(y, proba))
    except Exception as exc:                              # noqa: BLE001
        log.debug("binding CV failed: %s", exc)

    clf.fit(X, y)
    model = agonism.AgonismModel(
        clf=clf, train_fps=X, train_bvs=bvs, labels=y, auc=auc,
        n_agonist=int(y.sum()), n_antagonist=int((1 - y).sum()),
        agonist_mws=[])
    return model, auc


def enrichment_check(hits: list[ScreenHit]) -> dict:
    """Do known actives outrank the novel candidates?

    Important correction. This must judge the screen on the task the model was
    trained for. For a BINDING model, known actives should score high. For a
    direction model, known ligands scoring low can be entirely correct — most
    annotated ESR1 ligands are antagonists, so a low P(agonist) is the right
    answer, and an absolute threshold rejected that working model.

    So the test is relative: known actives must rank at or above the novel
    candidates. The absolute quality gate for a binding screen is the
    cross-validated AUC, checked before this function is reached.

    The `threshold` parameter of the superseded absolute test is gone
    rather than left in place unread: a caller passing it would have
    changed nothing and had no way to find that out.
    """
    known = [h.probability for h in hits if h.known_ligand]
    novel = [h.probability for h in hits if h.quadrant == "NOVEL HYPOTHESIS"]

    if not known:
        return {"usable": False, "n_known_recovered": 0,
                "n_novel_hypotheses": len(novel),
                "reason": "no known actives appear among the top hits — the "
                          "model is not ranking them above the library, so the "
                          "novel list cannot be trusted"}

    med_known = float(np.median(known))
    med_novel = float(np.median(novel)) if novel else 0.0
    usable = med_known >= med_novel

    out = {
        "usable": usable,
        "median_probability_known": round(med_known, 3),
        "median_probability_novel": round(med_novel, 3),
        "n_known_recovered": len(known),
        "n_novel_hypotheses": len(novel),
    }
    if not usable:
        out["reason"] = (f"known actives ({med_known:.2f}) rank below the novel "
                         f"candidates ({med_novel:.2f}) — the ranking is "
                         f"inverted, so the novel list is noise")
    return out


def screen_receptor(client: Client, target_chembl_id: str, *, goal: str = "up",
                    library: list[dict] | None = None, top_n: int = 50,
                    decoy_ratio: int = 3, seed: int = 0,
                    gene_linked: set | None = None) -> dict:
    """Screen for BINDING, then annotate direction where it can be predicted.

    Two models, two questions, two data sources:

      binding    /activity, abundant  -> which approved drugs hit this receptor?
      direction  /mechanism, scarce   -> would they agonise or antagonise it?

    The screen runs on the first. The second is applied as an annotation when
    enough direction labels exist, and simply omitted when they do not —
    rather than blocking the whole screen, which is what used to happen.
    """
    import random

    from . import agonism

    if library is not None:
        lib = library
    else:
        lib = local_library()
        if not lib and client is not None:
            lib = approved_library(client)
    if not lib:
        # `fam` does not exist in this function. An earlier edit put this
        # guard here instead of in screen_family, where it was meant to
        # go, and it would have raised NameError the first time a run had
        # no library at all.
        return {"target": target_chembl_id, "hits": [],
                "error": "no screening library: no local drug table on disk and no client to fetch one. Merge one with `python -m promoter_tf_drugs.localdrugs --merge FILE...`."}
    actives = agonism.fetch_actives(client, target_chembl_id)
    active_ids = {a["chembl_id"] for a in actives}

    if len(actives) < MIN_ACTIVES:
        return {"target": target_chembl_id, "hits": [], "n_actives": len(actives),
                "error": f"only {len(actives)} measured actives (pChEMBL >= 6) "
                         f"at this receptor — too few to train a binding model. "
                         f"This is a data gap, not a judgement about the "
                         f"receptor."}

    pool = [d for d in lib if d["chembl_id"] not in active_ids]
    rng = random.Random(seed)
    decoys = rng.sample(pool, min(len(pool), len(actives) * decoy_ratio))

    model, auc = build_binding_model(actives, decoys, seed=seed)
    if model is None:
        return {"target": target_chembl_id, "hits": [], "n_actives": len(actives),
                "error": "binding model could not be trained (structures "
                         "unparseable or too few decoys)"}

    if auc is not None and auc < MIN_BINDING_AUC:
        return {
            "target": target_chembl_id, "hits": [], "n_actives": len(actives),
            "model_auc": round(auc, 3),
            "error": f"binding model CV AUC {auc:.2f} is below "
                     f"{MIN_BINDING_AUC}: it cannot separate known actives from "
                     f"decoys, so its predictions on untested drugs are noise. "
                     f"Screen not reported.",
        }

    hits = screen(model, lib, goal="up",      # P(binder), not P(agonist)
                  known_ids=active_ids, top_n=top_n,
                  gene_linked=gene_linked,
                  known_potency={a["chembl_id"]: a["pchembl"]
                                 for a in actives if a.get("pchembl")})

    # Direction is a separate, optional annotation.
    ligands = agonism.fetch_labelled_ligands(client, target_chembl_id)
    dir_model = agonism.train_model(ligands)
    direction_note = None
    if dir_model is None:
        direction_note = (
            f"Direction not annotated: only {len(ligands)} drugs at this "
            f"receptor carry an agonist/antagonist label. The hits below are "
            f"predicted BINDERS — whether each would activate or inhibit is "
            f"unknown and needs a functional assay.")
    else:
        for h in hits:
            pred = agonism.predict(dir_model, next(
                (d["smiles"] for d in lib if d["chembl_id"] == h.chembl_id), ""),
                drug_name=h.drug, target_symbol=target_chembl_id)
            if pred.predicted != "undetermined":
                h.notes.append(
                    f"predicted {pred.predicted} (p={pred.probability_agonist}, "
                    f"{pred.confidence} confidence)")
                if goal == "up" and pred.predicted == "antagonist":
                    h.notes.append("predicted direction is OPPOSITE to your goal")
                elif goal == "down" and pred.predicted == "agonist":
                    h.notes.append("predicted direction is OPPOSITE to your goal")

    return {
        "target": target_chembl_id,
        "task": "binding",
        "model_auc": round(auc, 3) if auc is not None else None,
        "n_actives": len(actives),
        "n_decoys": len(decoys),
        "library_size": len(lib),
        "direction_note": direction_note,
        "hits": [h.to_dict() for h in hits],
        "enrichment": enrichment_check(hits),
    }


# ================================================================= #
# Screening a receptor that has too few ligands to model on its own #
# ================================================================= #
# The single-target screen needs MIN_ACTIVES measured ligands, which almost
# nothing clears. That is why the screening section returns the same
# retinoid and steroid receptors whatever gene goes in: the filter selects
# for how much medicinal chemistry a receptor has attracted.
#
# Nuclear receptors share a ligand-binding domain fold, so a model trained
# on a receptor's RELATIVES may still recognise its ligands. That is the
# standard way to borrow strength, and it is also an easy way to fool
# yourself, because a family model will happily report the family's
# chemistry as if it were the target's.
#
# So nothing is reported on faith. The model is trained on the siblings
# ONLY, with the target held out entirely, and then asked to rank the
# target's own measured ligands against decoys it has never seen. That
# leave-one-target-out AUC is a direct measurement of whether pooling
# transfers to this particular receptor, and it is the gate. A family whose
# chemistry does not transfer is refused and says so.

# Enough of the target's own ligands to measure transfer against. Below
# this the family model cannot be validated, so it is not reported - an
# unvalidated family model is exactly the thing that would fill the report
# with a family's pharmacology dressed as the target's.
MIN_VALIDATION_ACTIVES = 5

# Transfer AUC below which the family's chemistry does not carry to this
# receptor. Same bar the single-target model has to clear.
MIN_TRANSFER_AUC = MIN_BINDING_AUC


def transfer_auc(model, own_actives: list, decoys: list) -> float | None:
    """How well a sibling-trained model ranks the target's OWN ligands.

    The target contributed nothing to training, so this is an honest
    out-of-target estimate rather than a cross-validation within one
    receptor's chemistry.
    """
    from . import agonism

    if len(own_actives) < MIN_VALIDATION_ACTIVES or not decoys:
        return None
    try:
        from sklearn.metrics import roc_auc_score
    except ImportError:
        return None

    smiles = ([a["smiles"] for a in own_actives]
              + [d["smiles"] for d in decoys])
    X, _, keep = agonism.fingerprints(smiles)
    if len(keep) < MIN_VALIDATION_ACTIVES:
        return None
    import numpy as np
    y = np.array([1] * len(own_actives) + [0] * len(decoys))[keep]
    if y.sum() < MIN_VALIDATION_ACTIVES or (1 - y).sum() < 2:
        return None
    probs = model.clf.predict_proba(X)[:, 1]
    try:
        return float(roc_auc_score(y, probs))
    except ValueError:
        return None


def uniprot_for(symbol: str) -> str | None:
    """Gene symbol to UniProt accession, from the local HGNC table.

    Delegates to `localannot.uniprot_for`, which is the single place that
    knows the accession is a dict KEY on TFAnnotation rather than an
    attribute. Returns None when no HGNC table is on disk: unknown, not a
    gene without a protein.
    """
    try:
        from . import localannot as la
        return la.uniprot_for(symbol)
    except Exception:                                       # noqa: BLE001
        return None


def local_actives_for(symbol: str, *, uniprot: str | None = None) -> list:
    """Measured ligands of one protein, from a local potency file.

    The same records the online path builds from `/activity`, in the same
    shape - name, smiles, chembl_id, pchembl - so a family screen reads
    identically offline. Returns [] when no potency file is on disk, which
    is not the same as a protein with no known ligands.
    """
    uni = uniprot or uniprot_for(symbol)
    if not uni:
        return []
    try:
        from . import localpotency as lp
        if not any(lp.available().values()):
            return []
        return lp.actives_for(uni)
    except Exception:                                       # noqa: BLE001
        return []


def local_library(path=None) -> list:
    """The approved screening library from a local table, or [] if none."""
    try:
        from . import localdrugs as ld
        return ld.library(path) if path else ld.library()
    except Exception:                                       # noqa: BLE001
        return []


def screen_family(client, symbol: str, *, goal: str = "up",
                  library: list[dict] | None = None, top_n: int = 25,
                  decoy_ratio: int = 3, seed: int = 0,
                  gene_linked: set | None = None) -> dict:
    """Screen `symbol` using a model trained on its receptor family.

    For a receptor with too few of its own ligands to train on. Returns the
    same shape as `screen_receptor`, with `family_model: True` and the
    measured transfer AUC, or an `error` saying which step refused.

    Runs offline. `client` may be None, and local potency files are read
    first in any case: on a blocked network the ChEMBL activity calls are
    what left this path with nothing to train on, and a family screen that
    cannot read ligands reports "the family is as under-studied as the
    receptor" when the truth is that nothing could be fetched. Which
    source answered is recorded in `source`.
    """
    import random

    from . import agonism, paralogs

    fam = paralogs.family_of(symbol)
    if fam is None:
        return {"target": symbol, "hits": [],
                "error": f"{symbol} is not in a curated receptor family, so "
                         f"there are no relatives to borrow ligands from."}

    siblings = sorted(m for m in fam.members if m.upper() != symbol.upper())
    if not siblings:
        return {"target": symbol, "hits": [],
                "error": f"{symbol} is the only member of {fam.name} here."}

    # The target's own ligands, which are the validation set and are never
    # trained on.
    sources: list = []
    own = local_actives_for(symbol)
    if own:
        sources.append("local potency file")
    elif client is not None:
        from .drugs import _chembl_targets
        own_ids = _chembl_targets(client, symbol, None)
        for tid in own_ids[:2]:
            own.extend(agonism.fetch_actives(client, tid))
        if own:
            sources.append("ChEMBL")

    if len(own) < MIN_VALIDATION_ACTIVES:
        return {
            "target": symbol, "hits": [], "family": fam.name,
            "n_own_actives": len(own),
            "source": sources[0] if sources else None,
            "error": (f"{symbol} has {len(own)} measured ligand(s) of its "
                      f"own, below the {MIN_VALIDATION_ACTIVES} needed to "
                      f"test whether {fam.name} chemistry transfers to it. A "
                      f"family model could be built and not checked, which "
                      f"would report the family's pharmacology as this "
                      f"receptor's. Not reported."
                      + ("" if sources else
                         " No ligand source was readable at all - no local "
                         "potency file and no client - so this is a missing "
                         "source, not a receptor without ligands.")),
        }

    # Siblings' ligands, which are the training set.
    pooled: list = []
    contributing: list = []
    own_smiles = {a["smiles"] for a in own}
    for sib in siblings:
        got = local_actives_for(sib)
        if not got and client is not None:
            from .drugs import _chembl_targets
            tids = _chembl_targets(client, sib, None)
            for tid in tids[:1]:
                got.extend(agonism.fetch_actives(client, tid))
        # A ligand shared with the target would leak the answer.
        got = [g for g in got if g["smiles"] not in own_smiles]
        if got:
            contributing.append({"symbol": sib, "n_actives": len(got)})
            pooled.extend(got)

    if len(pooled) < MIN_ACTIVES:
        return {
            "target": symbol, "hits": [], "family": fam.name,
            "n_pooled_actives": len(pooled),
            "error": (f"{fam.name} contributed {len(pooled)} ligand(s) "
                      f"across {len(siblings)} relative(s), below the "
                      f"{MIN_ACTIVES} needed to train. The family is as "
                      f"under-studied as the receptor."),
        }

    if library is not None:
        lib = library
    else:
        lib = local_library()
        if not lib and client is not None:
            lib = approved_library(client)
    if not lib:
        return {"target": symbol, "hits": [], "family": fam.name,
                "error": "no screening library: no local drug table on disk and no client to fetch one. Merge one with `python -m promoter_tf_drugs.localdrugs --merge FILE...`."}
    pooled_smiles = {a["smiles"] for a in pooled}
    pool = [d for d in lib if d["smiles"] not in pooled_smiles
            and d["smiles"] not in own_smiles]
    rng = random.Random(seed)
    decoys = rng.sample(pool, min(len(pool), len(pooled) * decoy_ratio))

    model, cv_auc = build_binding_model(pooled, decoys, seed=seed)
    if model is None:
        return {"target": symbol, "hits": [], "family": fam.name,
                "error": "family binding model could not be trained"}

    # The gate: does this family's chemistry actually carry to this
    # receptor? Measured against ligands the model has never seen.
    val_decoys = rng.sample(pool, min(len(pool), max(20, len(own) * 3)))
    t_auc = transfer_auc(model, own, val_decoys)
    if t_auc is None:
        return {"target": symbol, "hits": [], "family": fam.name,
                "error": "transfer to this receptor could not be measured"}
    if t_auc < MIN_TRANSFER_AUC:
        return {
            "target": symbol, "hits": [], "family": fam.name,
            "model_auc": round(cv_auc, 3) if cv_auc else None,
            "transfer_auc": round(t_auc, 3),
            "error": (f"a model trained on {fam.name} relatives ranks "
                      f"{symbol}'s own {len(own)} known ligands at AUC "
                      f"{t_auc:.2f}, below {MIN_TRANSFER_AUC}. The family "
                      f"shares a fold and not a usable preference, so its "
                      f"chemistry does not predict this receptor. Not "
                      f"reported."),
        }

    # `goal` was accepted and then overridden with "up" here, so a run
    # asking to LOWER a gene was screened for compounds that raise it.
    hits = screen(model, lib, goal=goal,
                  known_ids={a["chembl_id"] for a in own},
                  top_n=top_n, gene_linked=gene_linked,
                  known_potency={a["chembl_id"]: a["pchembl"]
                                 for a in own if a.get("pchembl")})
    for h in hits:
        h.notes.append(
            f"predicted by a model trained on {symbol}'s {fam.name} "
            f"relatives, not on {symbol} itself. That model ranks "
            f"{symbol}'s own known ligands at AUC {t_auc:.2f}, which is why "
            f"it is reported - but a family model predicts family binding, "
            f"and subtype selectivity is exactly what it cannot see")

    return {
        "target": symbol, "hits": [h.to_dict() for h in hits],
        "family": fam.name, "family_model": True,
        "contributing": contributing,
        "n_own_actives": len(own), "n_pooled_actives": len(pooled),
        "n_decoys": len(decoys), "library_size": len(lib),
        "model_auc": round(cv_auc, 3) if cv_auc else None,
        "transfer_auc": round(t_auc, 3),
        "enrichment": enrichment_check(hits),
        "family_note": fam.note,
        "source": sources[0] if sources else None,
        "offline": "local potency file" in sources,
    }
