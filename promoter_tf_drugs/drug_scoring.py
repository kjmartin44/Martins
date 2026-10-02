"""Grade each TF/drug pair by combining TF-level confidence with drug quality.

Why the two halves multiply rather than add
-------------------------------------------
Getting a protein's output to move via this route requires a chain of things to
all be true:

    the TF really binds this promoter   AND   the drug really moves that TF

An additive score lets a superb drug rescue a worthless TF prediction, which is
exactly backwards: a 1 nM selective agonist of a factor that does not bind your
gene does nothing at all. So the combined grade is a weighted **geometric
mean** of the two confidences. Either one near zero drags the product to zero,
which is the honest behaviour.

Drug quality components (each mapped to 0..1)
---------------------------------------------
  potency      median pChEMBL for this drug against this target. pChEMBL 9 is
               1 nM, 6 is 1 uM. Below ~5 the drug is unlikely to engage the
               target at clinical exposure.
  selectivity  how many distinct targets the drug has an annotated mechanism
               for. A drug hitting 12 targets will move a lot more than your
               gene.
  action       how unambiguous the pharmacology is. A full agonist is a
               cleaner instrument than "positive modulator" or an unlabelled
               interaction from an aggregator.
  maturity     approval phase. Gated to phase 4 by default.
  safety       withdrawn or black-box drugs are penalised, not hidden. You
               should see them and decide.
  tissue       optional. With --tissue cns, drugs unlikely to cross the
               blood-brain barrier are down-weighted on physicochemistry.

None of these weights are fitted to outcome data. They are a transparent
triage ordering, and every component is printed alongside the grade so you can
see which one is carrying a given row.
"""

from __future__ import annotations

import logging
import statistics
from dataclasses import dataclass, asdict, field

from .http import Client, HttpError

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Relative weights inside the drug-quality half. Must be positive.
DRUG_WEIGHTS = {
    "potency": 0.25,
    "selectivity": 0.15,        # how many PROTEINS the drug binds
    "regulon": 0.15,            # how many GENES that protein controls
    "action": 0.15,
    "maturity": 0.10,
    "tissue": 0.10,
    "evidence_depth": 0.10,     # how much was actually measured
}

ACTION_CONFIDENCE = {
    "AGONIST": 1.00,
    "ANTAGONIST": 1.00,
    "INVERSE AGONIST": 0.95,
    "DEGRADER": 0.95,
    "INHIBITOR": 0.90,
    "BLOCKER": 0.90,
    "PARTIAL AGONIST": 0.70,
    "POSITIVE ALLOSTERIC MODULATOR": 0.65,
    "NEGATIVE ALLOSTERIC MODULATOR": 0.65,
    "POSITIVE MODULATOR": 0.60,
    "NEGATIVE MODULATOR": 0.60,
    "ACTIVATOR": 0.75,
    "OPENER": 0.70,
    "UNKNOWN": 0.30,
}

MATURITY = {4.0: 1.00, 3.0: 0.60, 2.0: 0.35, 1.0: 0.20, 0.0: 0.10}


@dataclass
class DrugEvidence:
    pchembl_median: float | None = None
    pchembl_n: int = 0
    n_targets: int | None = None
    regulon_size: int | None = None
    regulon_score: float | None = None
    regulon_note: str | None = None
    withdrawn: bool = False
    black_box: bool = False
    properties: dict = field(default_factory=dict)
    source: str = "ChEMBL"


@dataclass
class DrugGrade:
    tf_symbol: str
    drug_name: str
    combined: float                 # 0..100
    grade: str                      # A..F
    tf_confidence: float            # 0..1
    drug_quality: float             # 0..1
    components: dict                # each 0..1, for inspection
    evidence: dict
    notes: list[str]

    def to_dict(self) -> dict:
        return asdict(self)


# ----------------------------------------------------------------- evidence
def fetch_evidence(client: Client, drug) -> DrugEvidence:
    """Pull potency, promiscuity and safety flags from ChEMBL for one drug."""
    ev = DrugEvidence(withdrawn=bool(drug.withdrawn),
                      properties=dict(getattr(drug, "properties", None) or {}))
    # A local potency file first. On a blocked network the ChEMBL activity
    # call is what leaves this axis empty, and an offline run then scored
    # every drug as potency-unknown - which the axis treats as neutral, so
    # a nanomolar ligand and an untested one ranked the same.
    uni = getattr(drug, "uniprot", None) or getattr(
        drug, "target_uniprot", None)
    if uni:
        try:
            from . import localpotency as _lp
            if any(_lp.available().values()):
                local = _lp.actives_for(uni)
                if local:
                    smi = (getattr(drug, "properties", None) or {}).get(
                        "canonical_smiles")
                    vals = [r["pchembl"] for r in local
                            if smi is None or r["smiles"] == smi]
                    if not vals and smi is None:
                        vals = [r["pchembl"] for r in local]
                    if vals:
                        import statistics as _st
                        ev.pchembl_median = round(_st.median(vals), 2)
                        ev.pchembl_n = len(vals)
                        ev.source = "local potency file"
                        return ev
        except Exception as exc:                            # noqa: BLE001
            log.debug("local potency lookup failed: %s", str(exc)[:80])

    mol_id = drug.chembl_id
    if not mol_id or not str(mol_id).startswith("CHEMBL"):
        ev.source = "unavailable"
        return ev

    # Potency against *this* target specifically.
    target_id = getattr(drug, "target_chembl_id", None)
    if target_id:
        try:
            res = client.get(f"{CHEMBL}/activity", {
                "molecule_chembl_id": mol_id,
                "target_chembl_id": target_id,
                "pchembl_value__isnull": "false",
                "limit": 200, "format": "json",
            })
            vals = [float(a["pchembl_value"]) for a in res.get("activities", [])
                    if a.get("pchembl_value")]
            if vals:
                ev.pchembl_median = statistics.median(vals)
                ev.pchembl_n = len(vals)
        except (HttpError, ValueError, KeyError) as exc:
            log.debug("activity lookup failed for %s: %s", mol_id, exc)

    # Promiscuity: distinct targets with an annotated mechanism.
    try:
        res = client.get(f"{CHEMBL}/mechanism", {
            "molecule_chembl_id": mol_id, "limit": 200, "format": "json",
        })
        targets = {m.get("target_chembl_id") for m in res.get("mechanisms", [])
                   if m.get("target_chembl_id")}
        ev.n_targets = len(targets) or None
    except HttpError:
        pass

    # Safety flags.
    try:
        res = client.get(f"{CHEMBL}/drug_warning", {
            "molecule_chembl_id": mol_id, "limit": 50, "format": "json",
        })
        for w in res.get("drug_warnings", []):
            kind = str(w.get("warning_type", "")).lower()
            if "withdraw" in kind:
                ev.withdrawn = True
            if "black" in kind or "box" in kind:
                ev.black_box = True
    except HttpError:
        pass

    return ev


# ------------------------------------------------------------- normalisation
def _potency(ev: DrugEvidence) -> tuple[float, str | None]:
    """pChEMBL 5 -> 0, 9 -> 1. Unknown gets a neutral-but-uncertain 0.5."""
    if ev.pchembl_median is None:
        return 0.50, "no potency data; scored neutral"
    v = (ev.pchembl_median - 5.0) / 4.0
    note = None
    if ev.pchembl_median < 5.0:
        note = f"weak affinity (pChEMBL {ev.pchembl_median:.1f})"
    if ev.pchembl_n < 3:
        note = (note or "") + f" (n={ev.pchembl_n})"
    return max(0.0, min(1.0, v)), note


def _evidence_depth(ev: DrugEvidence) -> tuple[float, str | None]:
    """Guards against rewarding ignorance.

    `selectivity` counts *annotated* targets, so a drug studied for forty years
    accumulates a long list and looks promiscuous, while an obscure compound
    with three papers looks pristine. A low target count can mean selective or
    can mean nobody looked. This component separates the two: it scores how
    much evidence exists, so "clean" and "unstudied" stop looking identical.
    """
    n_act = ev.pchembl_n
    n_tgt = ev.n_targets or 0
    if n_act == 0 and n_tgt == 0:
        return 0.25, "barely characterised: no potency or mechanism data found"
    if n_act < 3:
        return 0.50, f"thin evidence base (n={n_act} potency measurements)"
    if n_act < 10:
        return 0.75, None
    return 1.0, None


def _selectivity(ev: DrugEvidence) -> tuple[float, str | None]:
    n = ev.n_targets
    if not n:
        return 0.50, None
    if n == 1:
        return 1.0, None
    score = max(0.15, 1.0 / (1.0 + 0.35 * (n - 1)))
    note = f"promiscuous: {n} annotated targets" if n >= 5 else None
    return score, note


def _tissue(ev: DrugEvidence, tissue: str | None) -> tuple[float, str | None]:
    """Crude CNS-exposure filter on physicochemistry. Only applied on request."""
    if tissue != "cns":
        return 1.0, None
    p = ev.properties or {}
    try:
        mw = float(p.get("full_mwt") or p.get("mw_freebase") or 0) or None
        tpsa = float(p.get("psa")) if p.get("psa") is not None else None
        logp = float(p.get("cx_logp")) if p.get("cx_logp") is not None else None
    except (TypeError, ValueError):
        return 0.6, "CNS exposure not assessed (0 of 6 parameters available)"

    have = sum(1 for v in (mw, tpsa, logp) if v is not None)
    from . import delivery as _dlv
    if have < _dlv.MIN_MPO_PARAMETERS:
        # Neutral, not a fail: we do not know, and scoring it either way would
        # be inventing information.
        return 0.6, (f"CNS exposure not assessed ({have} of 6 parameters "
                     "available)")

    # CNS MPO (Wager et al.) rather than a three-way cutoff.
    from . import delivery as dlv
    score, n_params, comp = dlv.cns_mpo(properties={
        "mw": mw, "tpsa": tpsa, "logp": logp})
    if score is None or n_params < dlv.MIN_MPO_PARAMETERS:
        return 0.6, (f"CNS exposure not assessed ({n_params} of 6 parameters "
                     "available)")
    frac = score / n_params
    note = None
    if frac < 0.67:
        worst = sorted(comp.items(), key=lambda kv: kv[1])[:2]
        note = (f"CNS MPO {score}/{n_params}; unlikely to reach brain "
                f"(weakest: {', '.join(k for k, _ in worst)})")
    return round(max(0.15, frac), 3), note


def _safety_multiplier(ev: DrugEvidence) -> tuple[float, list[str]]:
    mult, notes = 1.0, []
    if ev.withdrawn:
        mult *= 0.20
        notes.append("WITHDRAWN from market")
    if ev.black_box:
        mult *= 0.65
        notes.append("black-box warning")
    return mult, notes


# ---------------------------------------------------------------- the grade
def _letter(score: float) -> str:
    for cutoff, letter in ((70, "A"), (55, "B"), (40, "C"), (25, "D")):
        if score >= cutoff:
            return letter
    return "F"


def apply_direction_check(grade: "DrugGrade", prediction, goal_direction: str) -> "DrugGrade":
    """Fold a structural/ligand-based direction prediction into an existing grade.

    This only ever *confirms or contradicts* the curated ChEMBL action_type; it
    never replaces it. A database label from a curator who read the assay beats
    a model prediction. The prediction earns its keep when it disagrees - that
    is a flag to go read the primary literature before trusting the row.
    """
    if prediction is None or prediction.predicted == "undetermined":
        grade.notes.append("direction not independently predicted")
        return grade

    want = "agonist" if goal_direction == "up" else "antagonist"
    agrees = prediction.predicted == want

    if agrees:
        boost = {"high": 1.10, "moderate": 1.05}.get(prediction.confidence, 1.0)
        grade.combined = round(min(100.0, grade.combined * boost), 1)
        grade.notes.append(
            f"structure supports {prediction.predicted} "
            f"(p={prediction.probability_agonist}, {prediction.confidence} confidence)")
    else:
        grade.combined = round(grade.combined * 0.75, 1)
        grade.notes.append(
            f"CONFLICT: database says this should {goal_direction}-regulate, but "
            f"the structure-based model predicts {prediction.predicted} "
            f"(p={prediction.probability_agonist}). Check the primary literature.")
    grade.grade = _letter(grade.combined)
    grade.components["direction_model"] = 1.0 if agrees else 0.0
    return grade


def grade_drug(drug, tf_confidence: float, evidence: DrugEvidence, *,
               tissue: str | None = None, tf_weight: float = 0.5) -> DrugGrade:
    """Combine TF confidence and drug quality into a 0-100 grade."""
    notes: list[str] = []

    potency, n = _potency(evidence)
    if n:
        notes.append(n.strip())
    selectivity, n = _selectivity(evidence)
    if n:
        notes.append(n)

    # Allosteric sites are not constrained to bind the endogenous ligand, so
    # they diverge across a receptor family where the orthosteric pocket does
    # not. That is a real selectivity advantage, and a real dependency.
    from . import modulators as mod_mod
    mech = mod_mod.classify(drug)
    bonus = mod_mod.profile_for(mech).selectivity_bonus
    if bonus != 1.0:
        selectivity = min(1.0, selectivity * bonus)
    if mech in mod_mod.ALLOSTERIC:
        notes.append(f"{mech}: subtype selectivity typically better than "
                     "orthosteric, but requires endogenous ligand present")
    tissue_score, n = _tissue(evidence, tissue)
    if n:
        notes.append(n)
    depth, n = _evidence_depth(evidence)
    if n:
        notes.append(n)

    # Transcriptional specificity: how many genes ride along with the target.
    regulon = evidence.regulon_score if evidence.regulon_score is not None else 0.55
    if evidence.regulon_note:
        notes.append(evidence.regulon_note)

    action = ACTION_CONFIDENCE.get((drug.action_type or "UNKNOWN").upper(), 0.40)
    if action <= 0.40:
        notes.append(f"pharmacology poorly defined ({drug.action_type})")
    maturity = MATURITY.get(float(drug.max_phase or 0), 0.10)

    components = {
        "potency": potency,
        "selectivity": selectivity,
        "regulon": regulon,
        "action": action,
        "maturity": maturity,
        "tissue": tissue_score,
        "evidence_depth": depth,
    }

    # Weighted geometric mean: any single near-zero component dominates.
    total_w = sum(DRUG_WEIGHTS.values())
    quality = 1.0
    for key, w in DRUG_WEIGHTS.items():
        quality *= max(components[key], 1e-3) ** (w / total_w)

    safety_mult, safety_notes = _safety_multiplier(evidence)
    notes.extend(safety_notes)
    quality *= safety_mult

    tf_c = max(0.0, min(1.0, tf_confidence))
    if tf_c < 0.01:
        combined = 0.0
    else:
        combined = 100.0 * (tf_c ** tf_weight) * (quality ** (1.0 - tf_weight))

    if tf_c < 0.4:
        notes.append("weak TF evidence is the limiting factor here")

    return DrugGrade(
        tf_symbol=drug.target_symbol,
        drug_name=drug.name,
        combined=round(combined, 1),
        grade=_letter(combined),
        tf_confidence=round(tf_c, 3),
        drug_quality=round(quality, 3),
        components={k: round(v, 3) for k, v in components.items()},
        evidence=asdict(evidence),
        notes=notes,
    )


def grade_all(client: Client, tf_results, *, tissue: str | None = None,
              tf_weight: float = 0.5, max_tf_score: float = 100.0,
              evidence_lookup=None, predict_direction: bool = False,
              goal: str = "up", regulons: dict | None = None) -> list[DrugGrade]:
    """Grade every TF/drug pair. `evidence_lookup` overrides ChEMBL (demo/tests)."""
    grades: list[DrugGrade] = []
    model_cache: dict = {}

    for res in tf_results:
        tf_conf = max(0.0, min(1.0, res.score / max_tf_score))
        for drug in res.drugs:
            ev = (evidence_lookup(drug) if evidence_lookup
                  else fetch_evidence(client, drug))

            from . import regulon as regulon_mod
            score, n_targets, note = regulon_mod.transcriptional_specificity(
                res.symbol, regulons)
            ev.regulon_score, ev.regulon_size, ev.regulon_note = score, n_targets, note

            g = grade_drug(drug, tf_conf, ev, tissue=tissue, tf_weight=tf_weight)

            if predict_direction:
                try:
                    from . import agonism
                    pred = agonism.predict_for_drug(client, drug,
                                                    _model_cache=model_cache)
                    g = apply_direction_check(g, pred, goal)
                except ImportError as exc:
                    g.notes.append(f"direction model unavailable: {exc}")
                except Exception as exc:              # noqa: BLE001
                    log.warning("direction prediction failed for %s: %s",
                                drug.name, exc)
            grades.append(g)

    grades.sort(key=lambda g: -g.combined)
    return grades
