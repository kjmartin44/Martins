"""From "the drug binds the receptor" to "my gene goes up": the chain, explicitly.

The single most common error in this kind of analysis is treating receptor-level
agonism as if it were a transcriptional outcome. It is not, and the reason is
worth stating precisely because it is not a technicality.

Tamoxifen is an ERalpha *antagonist* in breast tissue and an *agonist* in
endometrium and bone. Same drug, same receptor, same binding pose. The
difference is which coregulators the cell has: helix 12's position creates a
surface, and whether that surface recruits a coactivator (SRC/p160, LXXLL) or a
corepressor (NCoR/SMRT, CoRNR) depends on their relative abundance in that cell
type. Raloxifene, bazedoxifene and every other SERM exist in this grey zone by
design.

So structure gives you one link. This module makes all the links visible:

    P(gene moves) ~ P(engages) x P(direction) x P(activates here) x P(responsive)

      engages       drug reaches and binds the receptor at clinical exposure
                    -> potency, tissue access                  [drug_scoring]
      direction     binding produces the functional direction you want
                    -> ligand model, conformational docking    [agonism/docking]
      activates     that factor acts as an activator at THIS promoter
                    -> GO role, motif quality, ChIP evidence   [pipeline]
      responsive    the gene is actually rate-limited by this factor in your
                    cell type -> nothing in this package knows this

Multiplying honest probabilities gives small numbers. That is the correct
result, not a bug to be tuned away: a four-link inferential chain with 0.7 at
each link lands near 0.24. Anyone reporting 0.9 for this has hidden a link.

The output is a probability *band*, not a point estimate, and the widest
contributing uncertainty is named so you know which experiment would narrow it.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field

# Ligands that only partially reposition helix 12 are the ones whose direction
# flips between tissues. Flagging them is more useful than forcing a call.
SERM_LIKE_NOTE = (
    "partial/ambiguous H12 engagement: this is the SERM pattern, where the same "
    "ligand activates in one tissue and represses in another. Direction here is "
    "a property of your cell type, not of the molecule"
)


@dataclass
class OutcomeEstimate:
    drug: str
    tf_symbol: str
    gene: str
    goal: str
    p_engages: float
    p_direction: float
    p_activates: float
    p_responsive: float
    p_point: float
    p_low: float
    p_high: float
    limiting_step: str
    tissue_dependent: bool = False
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    @property
    def verdict(self) -> str:
        if self.p_point >= 0.35:
            return "plausible, worth testing"
        if self.p_point >= 0.15:
            return "speculative"
        return "unlikely via this route"


def _band(p: float, uncertainty: float) -> tuple[float, float]:
    """A crude multiplicative band. Honest about being crude."""
    lo = max(0.0, p * (1.0 - uncertainty))
    hi = min(1.0, p * (1.0 + uncertainty))
    return round(lo, 3), round(hi, 3)


def estimate(drug_name: str, tf_symbol: str, gene: str, *, goal: str = "up",
             grade=None, prediction=None, tf_result=None,
             p_responsive: float | None = None) -> OutcomeEstimate:
    """Assemble the chain from whatever evidence the pipeline gathered."""
    notes: list[str] = []

    # --- link 1: does the drug engage the receptor at clinical exposure? ---
    if grade is not None:
        c = grade.components
        p_engages = min(0.95, 0.35 + 0.45 * c.get("potency", 0.5)
                        + 0.20 * c.get("tissue", 1.0))
    else:
        p_engages = 0.6
        notes.append("no potency data; engagement assumed moderate")

    # --- link 2: is the functional direction the one you want? ---
    tissue_dependent = False
    if prediction is not None and prediction.predicted != "undetermined":
        p_dir = prediction.probability_agonist or 0.5
        if goal == "down":
            p_dir = 1.0 - p_dir
        # Confidence shrinks the estimate toward 0.5 rather than trusting it.
        shrink = {"high": 0.9, "moderate": 0.7, "low": 0.4}.get(
            prediction.confidence, 0.3)
        p_direction = 0.5 + shrink * (p_dir - 0.5)

        if 0.35 <= (prediction.probability_agonist or 0.5) <= 0.65:
            tissue_dependent = True
            notes.append(SERM_LIKE_NOTE)
        if prediction.protrusion.get("antagonist_like_bulk"):
            tissue_dependent = True
    else:
        # Fall back on the curated database label, which is evidence, just not
        # structural evidence.
        p_direction = 0.75
        notes.append("direction from database annotation only, not independently "
                     "predicted")

    # --- link 3: does that factor activate THIS promoter? ---
    if tf_result is not None:
        tf_conf = max(0.0, min(1.0, tf_result.score / 100.0))
        p_activates = 0.15 + 0.70 * tf_conf
        if getattr(tf_result, "chip_peaks", 0):
            p_activates = min(0.95, p_activates * 1.25)
            notes.append("ChIP evidence supports binding at this promoter")
        else:
            notes.append("no measured binding at this promoter; motif only")
    else:
        p_activates = 0.35

    # --- link 4: is the gene rate-limited by this factor in your cells? ---
    if p_responsive is None:
        p_responsive = 0.40
        notes.append("gene responsiveness in your cell type is unknown and "
                     "assumed 0.40; this is the link the package cannot see")

    p = p_engages * p_direction * p_activates * p_responsive

    links = {
        "engagement (potency/exposure)": p_engages,
        "functional direction": p_direction,
        "activation at this promoter": p_activates,
        "gene responsiveness in your cells": p_responsive,
    }
    limiting = min(links, key=links.get)

    # Uncertainty grows with the number of weakly-supported links.
    weak = sum(1 for v in links.values() if v < 0.55)
    lo, hi = _band(p, 0.35 + 0.15 * weak)

    return OutcomeEstimate(
        drug=drug_name, tf_symbol=tf_symbol, gene=gene, goal=goal,
        p_engages=round(p_engages, 3), p_direction=round(p_direction, 3),
        p_activates=round(p_activates, 3), p_responsive=round(p_responsive, 3),
        p_point=round(p, 3), p_low=lo, p_high=hi,
        limiting_step=limiting, tissue_dependent=tissue_dependent, notes=notes,
    )


def narrowing_experiment(est: OutcomeEstimate) -> str:
    """Which experiment would most reduce the uncertainty in this estimate."""
    return {
        "engagement (potency/exposure)":
            "Measure target engagement at achievable exposure: a cellular "
            "thermal shift or a reporter dose-response.",
        "functional direction":
            "Coactivator recruitment assay (FRET/TR-FRET with an LXXLL peptide) "
            "in your cell type. This settles agonist vs antagonist directly and "
            "is the right assay when the SERM flag is raised.",
        "activation at this promoter":
            "ChIP-qPCR for the factor at this promoter +/- drug, then a "
            "promoter-reporter with the predicted element mutated.",
        "gene responsiveness in your cells":
            "Knock the TF down (siRNA/CRISPRi) in your cell type and measure the "
            "target gene. If it does not move, no drug acting through that "
            "factor will help, and everything upstream is moot.",
    }[est.limiting_step]
