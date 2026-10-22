"""Design the experiment that settles what this whole package only predicts.

Why this module exists
----------------------
Read the outcome chain in `outcome.py` and notice which link is almost always
limiting: *gene responsiveness in your cell type*. The pipeline cannot see it,
and it routinely caps the predicted probability at ~0.15 regardless of how good
the motif and the drug are. Every other module is an argument about what to
test. This one designs the test.

A CRISPRa screen across the ~1,600 human TFs answers the question directly and
without a prior: which factors, when forced onto your promoter, actually move
your reporter? That is `p_activates` and `p_responsive` measured rather than
assumed. An arrayed approved-compound screen (Prestwick, Selleck FDA library,
NCATS NPC) answers the drug question the same way, and will surface exactly the
"anti-parasitic that massively activates your promoter" case — because it makes
no assumption about mechanism at all.

The honest framing: the computational work here does not replace that screen.
It prioritises it, interprets its hits, and tells you which of the thousands of
possible follow-ups is worth doing. A screen hit with a predicted mechanism is
worth far more than a screen hit alone, and a predicted mechanism with no
screen behind it is a hypothesis.

What this module produces
-------------------------
1. A reporter construct plan from the promoter window already fetched, with the
   controls that make the result interpretable.
2. CRISPRa sgRNAs in the empirically optimal window for each prioritised TF.
3. A compound cherry-pick list ordered by the pipeline's grades, so a 1,500-well
   library can be run as a focused 96-well plate first.
4. Realistic hit-rate and power expectations, because a screen read without
   them produces false confidence in both directions.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, asdict, field

log = logging.getLogger(__name__)

# dCas9-VP64/SAM activation is strongly position-dependent. Empirically the
# productive window is roughly -400 to -50 relative to the TSS; outside it,
# activation falls off sharply even with a perfectly good guide.
CRISPRA_WINDOW = (-400, -50)

POL3_TERMINATOR = "TTTT"


@dataclass
class Guide:
    sequence: str                  # 20 nt protospacer, no PAM
    pam: str
    strand: str
    tss_offset: int
    gc: float
    in_optimal_window: bool
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ScreenPlan:
    gene: str
    reporter: dict | None
    crispra: dict
    compounds: dict
    expectations: dict
    # What the primary screen measures. A promoter reporter reads
    # transcription from a cloned fragment; the protein endpoint reads the
    # thing the patient needs. See `design`.
    endpoint: str = "promoter reporter"
    protein: dict | None = None

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------- sgRNA design
def find_guides(sequence: str, offset_fn, *, window: tuple = CRISPRA_WINDOW,
                gc_range: tuple = (0.35, 0.80), max_guides: int = 12
                ) -> list[Guide]:
    """Find NGG-PAM protospacers, preferring the CRISPRa-productive window.

    Deliberately not an off-target predictor: that needs a genome index and a
    proper tool (CRISPOR, CHOPCHOP). Guides from here are candidates to run
    through one of those, not finished reagents.
    """
    seq = sequence.upper()
    guides: list[Guide] = []

    def _add(protospacer, pam, strand, idx):
        if "N" in protospacer:
            return
        gc = (protospacer.count("G") + protospacer.count("C")) / len(protospacer)
        off = offset_fn(idx)
        in_win = window[0] <= off <= window[1]
        g = Guide(protospacer, pam, strand, off, round(gc, 3), in_win)
        if POL3_TERMINATOR in protospacer:
            g.warnings.append("contains TTTT: terminates U6 transcription, "
                              "this guide will not be made")
        if not (gc_range[0] <= gc <= gc_range[1]):
            g.warnings.append(f"GC {gc:.0%} outside {gc_range[0]:.0%}-"
                              f"{gc_range[1]:.0%}")
        if not in_win:
            g.warnings.append(
                f"offset {off:+d} is outside the productive CRISPRa window "
                f"({window[0]} to {window[1]}); activation falls off sharply here")
        guides.append(g)

    for m in re.finditer(r"(?=([ACGT]{21})GG)", seq):
        block = m.group(1)
        _add(block[1:], seq[m.start() + 21:m.start() + 23], "+", m.start() + 1)

    rc = sequence.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]
    n = len(seq)
    for m in re.finditer(r"(?=([ACGT]{21})GG)", rc):
        start_rc = m.start() + 1
        fwd_idx = n - (start_rc + 20)
        _add(m.group(1)[1:], rc[m.start() + 21:m.start() + 23], "-", max(fwd_idx, 0))

    # Best first: in-window, no warnings, GC near 0.55
    guides.sort(key=lambda g: (not g.in_optimal_window, len(g.warnings),
                               abs(g.gc - 0.55)))
    return guides[:max_guides]


# --------------------------------------------------------------- the plan
def reporter_plan(promoter) -> dict:
    return {
        "insert": f"{promoter.chrom}:{promoter.start}-{promoter.end} "
                  f"({promoter.upstream} bp upstream / {promoter.downstream} bp "
                  f"downstream of the {promoter.transcript_id} TSS)",
        "sequence_file": "promoter.fa (written by this run)",
        "vector": "promoterless luciferase (pGL4.10-style) or GFP for FACS-based "
                  "screening; luciferase for arrayed compound screens, GFP for "
                  "pooled CRISPRa",
        "controls": [
            "Promoterless vector: establishes the floor. Without it you cannot "
            "tell activation from vector background.",
            "Constitutive promoter (CMV/EF1a): establishes the ceiling and "
            "reports transfection or transduction efficiency per well.",
            "Scrambled-motif variant of your top predicted element: the single "
            "most informative control here. If a drug induces the wild-type "
            "reporter but not the mutant, the effect runs through that element. "
            "If it induces both, the mechanism is elsewhere and the prediction "
            "was wrong.",
            "Dual-reporter normalisation (Renilla or similar) to separate real "
            "induction from changes in cell number or general translation — "
            "cytotoxic compounds otherwise read as repressors.",
        ],
        "caveat": "An episomal reporter lacks the chromatin context of the "
                  "endogenous locus. A hit should be confirmed at the native "
                  "gene by RT-qPCR before it is believed.",
    }


def crispra_plan(promoter, tf_symbols: list[str], guides: list[Guide]) -> dict:
    usable = [g for g in guides if not g.warnings]
    return {
        "approach": "dCas9-VP64 with SAM (MS2-p65-HSF1) or dCas9-VPR for "
                    "stronger activation",
        "library": {
            "focused": f"{len(tf_symbols)} TFs prioritised by this run "
                       f"({', '.join(tf_symbols[:10])}"
                       f"{'...' if len(tf_symbols) > 10 else ''})",
            "genome_wide": "~1,600 human TFs (Calabrese or a TF-focused sublibrary). "
                           "Run this if the focused set comes back empty — an "
                           "unbiased screen is the only thing that finds a factor "
                           "nobody has associated with your promoter.",
            "coverage": "3-5 sgRNAs per gene, >=500 cells per guide at "
                        "transduction, MOI ~0.3 to keep single integrations",
        },
        "guides_at_target_promoter": [g.to_dict() for g in guides],
        "n_usable_guides": len(usable),
        "note": "Guides listed target YOUR promoter directly, for a "
                "positive-control arm proving the reporter can be activated at "
                "all. The library guides target the TFs. Run these through "
                "CRISPOR or CHOPCHOP for off-target scoring before ordering — "
                "no off-target analysis is done here.",
        "readout": "FACS sort GFP-high vs GFP-low, then sequence guide "
                   "abundance; or arrayed luciferase for a focused set",
    }


def compound_plan(graded_drugs: list[dict], direct_hits: list[dict]) -> dict:
    """Which individual compounds to order, rather than screening a library.

    "Cherry-pick" is the physical act of having a vendor pull named compounds
    into a custom plate, instead of buying and running a whole library. It is
    a cost-and-time decision: ~48 wells you can run in a week, against ~1,520
    for Prestwick or ~3,000 for the Selleck set.

    What it buys and what it costs, which is the part worth being explicit
    about: this list is assembled **entirely from the predictions made
    earlier in this same run** — the direct-therapy lookup and the graded
    TF/drug table. So it tests the model's own hypotheses quickly and
    cheaply, and it is **structurally incapable of finding a mechanism the
    model did not propose.** The unbiased full-library screen is the one that
    surfaces the anti-parasitic nobody would have predicted. These are
    complementary experiments, not substitutes, and running the cherry-pick
    first is a way to fail fast, not a way to be thorough.
    """
    cherry = []
    for d in direct_hits[:10]:
        cherry.append({"compound": d.get("drug"),
                       "rationale": "direct gene-targeted agent",
                       "source": "direct-therapy lookup (ChEMBL/DGIdb + "
                                 "curated table)",
                       "priority": 1})
    for g in graded_drugs[:40]:
        cherry.append({"compound": g.get("drug_name"),
                       "rationale": f"{g.get('tf_symbol')} route, grade "
                                    f"{g.get('grade')}",
                       "source": f"graded TF route via {g.get('tf_symbol')}",
                       "priority": 2 if g.get("grade") in ("A", "B") else 3})
    return {
        "n_available": len(direct_hits) + len(graded_drugs),
        "pulled_from": [
            f"{len(direct_hits)} direct gene-targeted agent(s) — priority 1, "
            f"because they act on the gene or its product and skip every "
            f"inferential step",
            f"{len(graded_drugs)} graded TF/drug route(s) — priority 2 if "
            f"graded A/B, else 3",
        ],
        "priority_key": {
            1: "acts directly on your gene, transcript or protein",
            2: "TF route the grading believes (A or B)",
            3: "TF route with a weaker grade — included to fill the plate, "
               "and the first rows to cut if wells are scarce",
        },
        "cannot_do": (
            "Every compound here was proposed by this run. A cherry-pick "
            "therefore tests the model and nothing beyond it — it cannot "
            "discover a mechanism the pipeline never predicted. If the whole "
            "plate is negative, the honest reading is that the promoter/TF "
            "model is wrong for this gene, not that no drug moves it."),
        "libraries": [
            "Prestwick Chemical Library (~1,520 approved drugs, high chemical "
            "and pharmacological diversity)",
            "Selleck FDA-approved Drug Library (~3,000)",
            "NCATS Pharmaceutical Collection (NPC) — free, and the most complete "
            "approved set",
        ],
        "cherry_pick": cherry[:48],
        "design": "Arrayed, 384-well, 8-point dose-response on the actives. "
                  "Single-concentration primary screens at 10 uM produce a "
                  "cytotoxicity-dominated hit list; dose-response separates real "
                  "induction from cells dying.",
        "counter_screen": "Run the same plate against a constitutive-promoter "
                          "reporter. Anything that moves both is affecting "
                          "transcription or translation globally, not your "
                          "promoter, and that is the majority of primary hits.",
        "note": "A full-library screen is unbiased and will find mechanisms this "
                "package cannot predict — that is its advantage, not a criticism "
                "of it. Use the cherry-pick list to get an answer in a week; run "
                "the full library when you want the answer you did not expect.",
    }


# The approved-drug collections a core facility is likely to hold. Sizes
# are approximate and drift as approvals change.
APPROVED_LIBRARIES = (
    "Broad Drug Repurposing Hub, approved subset (~2,400 compounds)",
    "Prestwick Chemical Library (~1,520, all off-patent and approved)",
    "Selleck FDA-approved Drug Library (~2,900)",
    "NCATS Pharmaceutical Collection (NPC, ~2,800)",
)


def protein_abundance_plan(gene: str, *, why: str = "") -> dict:
    """A screen for compounds that raise `gene`'s PROTEIN in patient cells.

    Written against the failure mode this package spent its effort
    characterising: for a cell-cycle gene, anything that changes the share
    of dividing cells changes the measured protein without touching the
    gene. Every control below exists to separate those two, because
    without them the screen's top hits will be antiproliferatives.
    """
    g = (gene or "the gene").strip().upper()
    return {
        "endpoint": f"{g} protein per cell, in the patient's own cells",
        "why": why or (f"the declared problem is the {g} protein rather "
                       f"than the amount of transcript"),
        "cells": [
            "Patient-derived primary fibroblasts (skin biopsy) or EBV "
            "transformed lymphoblastoid lines. Primary fibroblasts are "
            "closer to the biology; LCLs scale better to 384-well.",
            "On the same plate: a carrier parent if available, and at "
            "least two unrelated healthy donors. The carrier is the "
            "informative control, because the target level is what a "
            "healthy heterozygote already demonstrates is tolerated.",
            "Keep every well sub-confluent and proliferating. This gene's "
            "protein tracks the cell cycle, so a density difference "
            "between plates is indistinguishable from a drug effect.",
        ],
        "readout": [
            "Primary: high-content immunofluorescence, per-cell protein "
            "intensity. Imaging is worth the extra effort here because it "
            "gives per-cell data, and per-cell data is what lets you gate "
            "on cell-cycle phase.",
            "Alternatives at 384-well: in-cell western (two-colour, "
            "normalised to a total-protein stain) or HTRF/AlphaLISA. Both "
            "give a well average, so neither can separate 'more protein "
            "per cell' from 'more cells in G2/M'.",
            "Confirmation: capillary or conventional western on the "
            "confirmed hits, which is the only readout that shows you the "
            "band is the right size - relevant when one allele may make a "
            "truncated product.",
            "Validate the antibody against a knockdown before screening. "
            "An antibody that does not drop with siRNA is measuring "
            "something else, and the whole screen then measures it too.",
        ],
        "normalisation": [
            "Total protein stain, or a housekeeping protein that is NOT "
            "cell-cycle regulated. Not PCNA, not a cyclin, not Ki-67.",
            "Co-stain DNA content or phospho-histone H3 in the same well. "
            "This is the control that matters most: it turns 'the signal "
            "went up' into 'the signal went up at matched cell-cycle "
            "phase', and without it the screen's strongest hits will be "
            "compounds that arrest cells in mitosis.",
        ],
        "library": {
            "what": "an FDA-approved collection, so anything that hits is "
                    "already through safety",
            "sources": list(APPROVED_LIBRARIES),
            "format": "single concentration primary, typically 10 uM, in "
                      "duplicate; two exposures (24 h and 48 h), because a "
                      "stabiliser and a transcriptional effect separate in "
                      "time",
        },
        "controls": [
            "Positive: a proteasome inhibitor (MG132, or bortezomib). It "
            "raises almost any protein, which makes it a poor therapy and "
            "an excellent control - it proves the assay can detect an "
            "increase at all, and it sets the ceiling.",
            "Negative: vehicle (DMSO) at the screening concentration, many "
            "wells, spread across the plate rather than in one column.",
            "Cytotoxicity in parallel (CellTiter-Glo or equivalent) on a "
            "sister plate. A dying well reads as whatever the "
            "normalisation does under stress, and most primary hits in "
            "screens like this are simply toxic.",
            "Cell-cycle counter-screen on the primary hits: anything that "
            "arrests cells in G2/M raises this protein for a reason that "
            "is not a treatment.",
        ],
        "thresholds": [
            "Z' >= 0.5 between the proteasome inhibitor and vehicle before "
            "screening anything. Below that the plate cannot tell a hit "
            "from noise.",
            "Primary hit: >= 1.5x vehicle, after excluding wells with more "
            "than 30% viability loss or a mitotic index more than twice "
            "vehicle.",
            "Confirm on an 8-point dose response in triplicate. A real "
            "stabiliser gives a dose-response; a plate artefact does not.",
        ],
        "triage": [
            "1. Dose-response confirms, cytotoxicity excluded.",
            "2. Western confirms the increase and the band size.",
            "3. The increase survives cell-cycle matching.",
            "4. RT-qPCR says whether the transcript moved too. Protein up "
            "with transcript flat points at stability or translation; both "
            "up points at transcription, and that route is testable "
            "against the rest of this report.",
            "5. Only then is it a candidate, and it is a candidate for a "
            "cell assay of function, not for a patient.",
        ],
        "caveat": "A screen in fibroblasts or LCLs is not the affected "
                  "tissue. A compound that raises the protein there may "
                  "not reach or work in the tissue that matters, and the "
                  "exposure check in `exposure.py` is the next question "
                  "for anything that survives.",
    }


def protein_expectations() -> dict:
    """What an approved-library protein screen realistically returns."""
    return {
        "compound_screen_hit_rate": "0.5-2% primary at a sensible "
                                    "threshold, and most of those are "
                                    "cytotoxic or cell-cycle effects. "
                                    "Expect to lose 60-80% at "
                                    "counter-screen and to finish with "
                                    "single digits of confirmed compounds, "
                                    "possibly none.",
        "replicates": "duplicate primary, triplicate for confirmation. A "
                      "single primary well cannot be distinguished from a "
                      "pipetting error.",
        "quality_metric": "Z' >= 0.5 against the proteasome-inhibitor "
                          "control, per plate, recorded for every plate.",
        "what_a_negative_result_means": "No approved compound raises this "
                                        "protein in this cell type at this "
                                        "exposure. That is an informative "
                                        "result and a cheap one: it says "
                                        "repurposing is unlikely to solve "
                                        "it, and redirects effort to "
                                        "chemistry aimed at the specific "
                                        "mechanism - the degradation route "
                                        "or the splice defect - rather "
                                        "than to more screening.",
        "why_this_screen_and_not_the_reporter": "A promoter reporter "
                                                "measures transcription "
                                                "from a cloned fragment. "
                                                "It cannot see a "
                                                "stabiliser, a splice "
                                                "modifier or a readthrough "
                                                "agent, which are the "
                                                "mechanisms in play when "
                                                "the protein is the "
                                                "problem.",
    }


def expectations() -> dict:
    return {
        "compound_screen_hit_rate": "0.5-2% at a sensible threshold. Most "
                                    "primary hits are cytotoxic or "
                                    "assay-interfering; expect to lose 60-80% at "
                                    "counter-screen.",
        "crispra_hit_rate": "typically a handful of TFs from ~1,600. A screen "
                            "returning hundreds of hits is reporting a global "
                            "effect on the reporter, not specific regulation.",
        "replicates": "3 biological replicates minimum; 2 gives you no estimate "
                      "of variance and a Z' you cannot trust",
        "quality_metric": "Z' >= 0.5 between positive and negative controls "
                          "before you screen anything. Below that the screen "
                          "cannot distinguish a real hit from noise and running "
                          "it wastes the library.",
        "what_a_negative_result_means": "If no TF in a genome-wide CRISPRa "
                                        "screen activates your reporter, the "
                                        "promoter fragment is probably not "
                                        "sufficient — the regulation likely lives "
                                        "in an enhancer outside your window. "
                                        "Widen it (--upstream 50000) or use the "
                                        "endogenous locus.",
    }


def design(promoter, tf_symbols: list[str], graded_drugs: list[dict],
           direct_hits: list[dict], *, protein_endpoint: bool = False,
           why: str = "") -> ScreenPlan:
    """The experiment this analysis prioritises.

    `protein_endpoint` swaps the promoter reporter for a screen that
    measures the gene's PROTEIN in patient cells. Set it when the
    declared problem is the protein rather than the transcript - a
    hypomorphic allele, a premature stop, a truncating second allele -
    because a reporter cannot see a stabiliser, a splice modifier or a
    readthrough agent, and those are the mechanisms in play.

    The CRISPRa arm is kept either way: it raises the endogenous gene, so
    it is a positive control for "would more of this protein help at all"
    whichever endpoint is read.
    """
    guides = find_guides(promoter.sequence, promoter.offset_to_tss)
    if protein_endpoint:
        return ScreenPlan(
            gene=promoter.gene_symbol,
            reporter=None,
            crispra=crispra_plan(promoter, tf_symbols, guides),
            compounds=compound_plan(graded_drugs, direct_hits),
            expectations=protein_expectations(),
            endpoint=f"{promoter.gene_symbol} protein in patient cells",
            protein=protein_abundance_plan(promoter.gene_symbol, why=why),
        )
    return ScreenPlan(
        gene=promoter.gene_symbol,
        reporter=reporter_plan(promoter),
        crispra=crispra_plan(promoter, tf_symbols, guides),
        compounds=compound_plan(graded_drugs, direct_hits),
        expectations=expectations(),
    )
