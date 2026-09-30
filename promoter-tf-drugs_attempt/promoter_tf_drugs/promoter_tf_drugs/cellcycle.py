"""Cell-cycle phase dependence of TF activity — a separate axis from everything else.

Not the same thing as a PAM preserving dynamics
-----------------------------------------------
A positive allosteric modulator preserves the *signalling* dynamics of the
endogenous ligand: pulsatile, circadian, tissue-restricted patterning survives
because the PAM scales what is already happening. That is a property of the
ligand-receptor relationship, on a timescale of seconds to hours.

Cell-cycle dependence is a different clock and a different failure mode. It
asks: **is this transcription factor even active in the cells you are
treating, right now?** A PAM at a G1/S-restricted factor still does nothing in
a cell sitting in G2.

Two consequences, and they are not subtle
-----------------------------------------
**1. Duty cycle dilutes potency in asynchronous populations.** E2F1 is active
across roughly the late-G1-to-early-S window. In an unsynchronised culture,
only that fraction of cells can respond at any moment. A drug acting through
E2F1 therefore has its effective population-level efficacy multiplied by the
duty cycle before anything else is considered — a real potency penalty that no
binding assay reveals, and one that makes a bulk RT-qPCR look like a failed
experiment when the per-cell effect was fine.

**2. Post-mitotic tissue is not merely "slower".** In neurons, cardiomyocytes
and other terminally differentiated cells, E2F-dependent promoters are held
down by the DREAM complex (DP, RB-like p130, E2F4, MuvB). The targets are not
quiescent-but-available; they are *actively repressed*. A route through a
G1/S factor in post-mitotic tissue is not weak, it is blocked. If your target
tissue is post-mitotic — which is exactly the case for most neurological and
cardiac indications — those routes should be removed from consideration, not
ranked lower.

The exception worth knowing: **mitotic bookmarking.** Most TFs are evicted from
chromatin during mitosis, but a minority stay bound and reseed transcription in
the daughter cells. Those are less sensitive to cycle phase.

Tell the pipeline what you are treating with `--proliferation`:
    cycling       tumour lines, proliferating primary cells
    post-mitotic  neurons, cardiomyocytes, mature muscle
    mixed         tissue containing both
    unknown       default; flags the issue without penalising
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

log = logging.getLogger(__name__)

# Approximate fraction of a ~24 h cycle during which the factor is active.
# Coarse by design: G1 ~46%, S ~33%, G2 ~17%, M ~4%.
PHASE_DUTY = {
    "G1": 0.46,
    "G1/S": 0.25,
    "S": 0.33,
    "G2": 0.17,
    "G2/M": 0.21,
    "M": 0.05,
    "constitutive": 1.0,
}


@dataclass
class PhaseProfile:
    symbol: str
    phase: str
    duty: float
    post_mitotic_active: bool
    bookmarking: bool
    note: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def _p(symbol, phase, post_mitotic, bookmarking=False, note=""):
    return PhaseProfile(symbol, phase, PHASE_DUTY[phase], post_mitotic,
                        bookmarking, note)


# Curated. Phase assignments are the well-established ones; treat as coarse.
PROFILES: dict[str, PhaseProfile] = {p.symbol: p for p in [
    # --- G1/S: the E2F/DREAM axis. Blocked in post-mitotic cells. -----------
    _p("E2F1", "G1/S", False, note="DREAM complex actively represses E2F targets "
                                   "in post-mitotic cells; not merely quiescent"),
    _p("E2F2", "G1/S", False), _p("E2F3", "G1/S", False),
    _p("E2F7", "S", False), _p("E2F8", "S", False),
    _p("MYBL2", "G2/M", False, note="B-MYB; MuvB core, proliferation-restricted"),
    _p("TFDP1", "G1/S", False),
    # E2F4 is the exception: it persists and acts as a repressor in quiescence.
    _p("E2F4", "G1", True, note="repressive arm of DREAM; active in quiescent and "
                                "post-mitotic cells, unlike activating E2Fs"),

    # --- G2/M ---------------------------------------------------------------
    _p("FOXM1", "G2/M", False, note="proliferation-restricted; essentially absent "
                                    "from post-mitotic tissue"),
    _p("MYBL1", "G2/M", False),

    # --- constitutive / cycle-independent, incl. mitotic bookmarking -------
    _p("SP1", "constitutive", True),
    _p("CREB1", "constitutive", True),
    _p("NFKB1", "constitutive", True), _p("RELA", "constitutive", True),
    _p("STAT3", "constitutive", True),
    _p("NR3C1", "constitutive", True), _p("VDR", "constitutive", True),
    _p("PPARG", "constitutive", True), _p("RARA", "constitutive", True),
    _p("RXRA", "constitutive", True), _p("THRB", "constitutive", True),
    _p("ESR1", "constitutive", True), _p("AR", "constitutive", True),
    _p("HNF4A", "constitutive", True),
    _p("CTCF", "constitutive", True, bookmarking=True),
    _p("GATA1", "constitutive", True, bookmarking=True,
       note="retained on mitotic chromatin; reseeds transcription in daughters"),
    _p("FOXA1", "constitutive", True, bookmarking=True,
       note="pioneer factor, mitotic bookmarking"),
    _p("ESRRB", "constitutive", True, bookmarking=True),
    _p("RUNX2", "constitutive", True, bookmarking=True),
    _p("REST", "constitutive", True),
    _p("MYC", "G1", True, note="cycle-associated but expressed in non-cycling "
                               "cells too; growth rather than phase"),
]}

# Factors that only make sense in dividing cells, by family prefix.
PROLIFERATION_PREFIXES = ("E2F", "CCN", "CDK", "FOXM", "MYBL", "MCM", "PLK", "AURK")


def profile_for(symbol: str) -> PhaseProfile:
    s = symbol.upper()
    if s in PROFILES:
        return PROFILES[s]
    if s.startswith(PROLIFERATION_PREFIXES):
        return _p(s, "G1/S", False,
                  note="assigned by family prefix; verify before relying on it")
    # Unknown factors are assumed cycle-independent, which is the common case
    # and the conservative one (no unearned penalty).
    return PhaseProfile(s, "unknown", 1.0, True, False,
                        "no cell-cycle annotation; assumed cycle-independent")


@dataclass
class CycleAssessment:
    symbol: str
    phase: str
    duty: float
    proliferation: str
    compatible: bool
    multiplier: float          # applied to efficacy/efficiency
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def assess(symbol: str, proliferation: str = "unknown") -> CycleAssessment:
    """How much does cell-cycle phase cost this TF route, in this tissue?"""
    prof = profile_for(symbol)
    notes: list[str] = []
    compatible = True
    mult = 1.0

    if proliferation == "post-mitotic":
        if not prof.post_mitotic_active:
            compatible = False
            mult = 0.0
            notes.append(
                f"{symbol} is proliferation-restricted ({prof.phase}) and your "
                "tissue is post-mitotic. This route is blocked, not merely weak: "
                "E2F/DREAM targets are actively repressed in terminally "
                "differentiated cells.")
        else:
            notes.append("cycle-independent; fine in post-mitotic tissue")

    elif proliferation == "cycling":
        if prof.duty < 1.0:
            mult = prof.duty
            notes.append(
                f"active only in {prof.phase} (~{prof.duty:.0%} of the cycle). In "
                f"an asynchronous population the population-level effect is "
                f"diluted accordingly — per-cell efficacy may be far higher than "
                f"a bulk measurement suggests.")
            if prof.duty <= 0.25:
                notes.append(
                    "consider synchronising, or a single-cell/phase-resolved "
                    "readout, before concluding the drug does not work")

    elif proliferation == "mixed":
        if not prof.post_mitotic_active:
            mult = 0.5 * prof.duty
            notes.append(
                f"{symbol} only responds in the dividing fraction of a mixed "
                "population; effect size will track the proliferative index")
        elif prof.duty < 1.0:
            mult = 0.5 + 0.5 * prof.duty

    else:   # unknown
        if prof.duty < 1.0 or not prof.post_mitotic_active:
            notes.append(
                f"{symbol} is cell-cycle restricted ({prof.phase}). Specify "
                "--proliferation for your tissue: this route is blocked in "
                "post-mitotic cells and diluted in asynchronous ones.")

    if prof.bookmarking:
        notes.append("mitotic bookmarking factor: stays chromatin-bound through "
                     "mitosis, so expression resumes promptly in daughter cells")
    if prof.note:
        notes.append(prof.note)

    return CycleAssessment(symbol, prof.phase, prof.duty, proliferation,
                           compatible, round(mult, 3), notes)


def filter_and_annotate(tf_results, proliferation: str = "unknown"):
    """Apply cycle multipliers to TF scores; return per-TF assessments.

    Incompatible routes are zeroed rather than merely demoted, because in
    post-mitotic tissue they do not work at all and ranking them below other
    options would still invite someone to try them.
    """
    assessments = {}
    for r in tf_results:
        a = assess(r.symbol, proliferation)
        assessments[r.symbol] = a.to_dict()
        if a.multiplier != 1.0:
            r.score = round(r.score * a.multiplier, 1)
    tf_results.sort(key=lambda r: -r.score)
    return assessments


def summarise(assessments: dict, proliferation: str) -> str | None:
    blocked = [s for s, a in assessments.items() if not a["compatible"]]
    diluted = [s for s, a in assessments.items()
               if a["compatible"] and a["multiplier"] < 1.0]
    if proliferation == "unknown":
        restricted = [s for s, a in assessments.items() if a["duty"] < 1.0]
        if restricted:
            return (f"{len(restricted)} factor(s) here are cell-cycle restricted "
                    f"({', '.join(restricted[:5])}). Set --proliferation to have "
                    "this scored rather than just flagged.")
        return None
    parts = []
    if blocked:
        parts.append(f"{len(blocked)} route(s) removed as incompatible with "
                     f"post-mitotic tissue: {', '.join(blocked[:5])}")
    if diluted:
        parts.append(f"{len(diluted)} route(s) down-weighted for partial duty "
                     f"cycle: {', '.join(diluted[:5])}")
    return ". ".join(parts) if parts else None
