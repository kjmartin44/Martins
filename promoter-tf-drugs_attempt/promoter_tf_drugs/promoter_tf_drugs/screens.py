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
    reporter: dict
    crispra: dict
    compounds: dict
    expectations: dict

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
    cherry = []
    for d in direct_hits[:10]:
        cherry.append({"compound": d.get("drug"), "rationale": "direct gene-"
                       "targeted agent", "priority": 1})
    for g in graded_drugs[:40]:
        cherry.append({"compound": g.get("drug_name"),
                       "rationale": f"{g.get('tf_symbol')} route, grade "
                                    f"{g.get('grade')}",
                       "priority": 2 if g.get("grade") in ("A", "B") else 3})
    return {
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
           direct_hits: list[dict]) -> ScreenPlan:
    guides = find_guides(promoter.sequence, promoter.offset_to_tss)
    return ScreenPlan(
        gene=promoter.gene_symbol,
        reporter=reporter_plan(promoter),
        crispra=crispra_plan(promoter, tf_symbols, guides),
        compounds=compound_plan(graded_drugs, direct_hits),
        expectations=expectations(),
    )
