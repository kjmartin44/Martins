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
    probability: float             # of the wanted direction
    nearest_tanimoto: float
    novelty: float                 # 1 - nearest_tanimoto
    quadrant: str
    known_ligand: bool
    max_phase: float | None
    indication_hint: str | None = None
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
           top_n: int = 50, gene_linked: set | None = None) -> list[ScreenHit]:
    """Score a drug library against a trained per-receptor model."""
    from . import agonism

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
                "binds this receptor with measured activity, and has no "
                "reported effect on your gene. The pharmacology is already "
                "de-risked; only the receptor-to-gene step is untested, which "
                "makes this the cheapest experiment on the list.")
        elif quad == "confirms a reported link":
            notes.append("binds this receptor AND is already reported to "
                         "affect your gene — a mechanism for a known effect, "
                         "not a new lead")
        elif quad == "no call":
            notes.append("outside applicability domain; score withheld")

        hits.append(ScreenHit(
            drug=drug["name"], chembl_id=drug["chembl_id"],
            probability=round(p_want, 3), nearest_tanimoto=round(tan, 3),
            novelty=round(1.0 - tan, 3), quadrant=quad, known_ligand=known,
            max_phase=drug.get("max_phase"), notes=notes,
        ))

    # Established binders that are untested on the gene come first: they are
    # the most actionable, not a control to be skipped past.
    order = {"ESTABLISHED LIGAND, untested on this gene": 0,
             "confirms a reported link": 1,
             "predicted binder, close to known chemistry": 2,
             "NOVEL HYPOTHESIS": 3, "no call": 4}
    hits.sort(key=lambda h: (order[h.quadrant], -h.probability))
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


def enrichment_check(hits: list[ScreenHit],
                     threshold: float = PROMISCUITY_FLAG) -> dict:
    """Do known actives outrank the novel candidates?

    Important correction. This must judge the screen on the task the model was
    trained for. For a BINDING model, known actives should score high. For a
    direction model, known ligands scoring low can be entirely correct — most
    annotated ESR1 ligands are antagonists, so a low P(agonist) is the right
    answer, and an absolute threshold rejected that working model.

    So the test is relative: known actives must rank at or above the novel
    candidates. The absolute quality gate for a binding screen is the
    cross-validated AUC, checked before this function is reached.
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

    lib = library if library is not None else approved_library(client)
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
                  gene_linked=gene_linked)

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
