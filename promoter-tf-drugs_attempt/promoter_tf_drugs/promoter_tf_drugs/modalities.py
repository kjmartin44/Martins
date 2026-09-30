"""Design a direct modality for YOUR gene, rather than looking up existing ones.

`direct.py` reports approved agents that already target a gene. For most genes
that returns nothing, which is not the same as "no direct route is possible".
This module designs one: siRNA, RNase H gapmer ASO, splice-switching oligo,
CRISPRi, or nuclease/base editing — and then tells you the thing that actually
decides whether any of it is feasible.

Delivery is the constraint, not design
--------------------------------------
Designing a potent siRNA is close to a solved problem: the sequence rules
below are well established and a competent design will usually knock down
80-95% in transfected cells. Programs do not fail there. They fail because the
molecule never reaches the cell type that matters.

What is actually solved, in vivo:

  hepatocytes      GalNAc conjugation, subcutaneous, quarterly-ish dosing.
                   Essentially routine — inclisiran, givosiran, lumasiran,
                   nedosiran, vutrisiran are all this.
  CNS neurons/glia Intrathecal delivery. Works, and requires repeated lumbar
                   puncture — nusinersen, tofersen.
  eye              Intravitreal injection.
  lung             Inhaled, for some chemistries.

What is not:

  T cells, macrophages, microglia — no routine systemic oligo delivery exists.
  This is the blocker for immune-cell targets, and it is why CCR5 work has gone
  the ex vivo route (edit HSCs or T cells outside the body, then reinfuse)
  rather than a systemic knockdown. Ex vivo editing is a real, demonstrated
  path, but it is a transplant programme, not a drug.

So this module scores deliverability first and design second. A beautiful siRNA
for a target in microglia is not a candidate; it is a cell-culture reagent.

Off-target
----------
Seed-region complementarity (positions 2-8 of the guide) drives most siRNA
off-target silencing, and it is a transcriptome-wide question. The checks here
are local sanity only. Every candidate must go through a proper homology search
(BLAST against RefSeq, or a seed-match scan) before it means anything.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, asdict, field

log = logging.getLogger(__name__)

# CRISPRi (dCas9-KRAB) works in a different window from CRISPRa: centred on the
# TSS rather than well upstream of it.
CRISPRI_WINDOW = (-50, 300)

DELIVERY = {
    "liver": {
        "status": "solved", "score": 1.0,
        "route": "GalNAc-conjugated siRNA or ASO, subcutaneous",
        "note": "hepatocyte delivery is routine; several approved agents",
    },
    "cns": {
        "status": "solved with a caveat", "score": 0.75,
        "route": "intrathecal ASO (lumbar puncture)",
        "note": "works and is approved (nusinersen, tofersen), but means "
                "repeated LPs; reaches neurons and glia, including microglia "
                "less efficiently",
    },
    "eye": {
        "status": "solved", "score": 0.85,
        "route": "intravitreal injection",
        "note": "local delivery, small doses",
    },
    "lung": {
        "status": "partial", "score": 0.55,
        "route": "inhaled oligonucleotide",
        "note": "chemistry-dependent; mucus and nuclease barriers",
    },
    "kidney": {
        "status": "partial", "score": 0.5,
        "route": "systemic PS-ASO accumulates in proximal tubule",
        "note": "accumulation is real but so is tubular toxicity",
    },
    "immune cells": {
        "status": "NOT solved in vivo", "score": 0.15,
        "route": "ex vivo editing of HSCs or T cells, then reinfusion",
        "note": "no routine systemic oligo delivery to T cells, macrophages or "
                "microglia. The demonstrated path is ex vivo — which is a "
                "transplant programme with conditioning, not a drug",
    },
    "muscle": {
        "status": "partial", "score": 0.45,
        "route": "high-dose systemic PMO, or AAV",
        "note": "the DMD exon-skippers work here but with famously modest "
                "dystrophin restoration; muscle is a large, poorly-penetrated "
                "compartment",
    },
    "tumour": {
        "status": "partial", "score": 0.4,
        "route": "LNP, intratumoural, or antibody conjugate",
        "note": "heterogeneous perfusion; delivery is the recurring failure",
    },
    "unknown": {
        "status": "unknown", "score": 0.5,
        "route": "—",
        "note": "specify a tissue: deliverability decides feasibility more than "
                "design quality does",
    },
}

# Map the pipeline's --tissue and common cell types onto delivery compartments.
TISSUE_ALIASES = {
    "cns": "cns", "brain": "cns", "neuron": "cns", "microglia": "immune cells",
    "liver": "liver", "hepatocyte": "liver",
    "t cell": "immune cells", "t-cell": "immune cells", "macrophage": "immune cells",
    "blood": "immune cells", "immune": "immune cells", "lymphocyte": "immune cells",
    "eye": "eye", "retina": "eye", "lung": "lung", "kidney": "kidney",
    "muscle": "muscle", "tumour": "tumour", "tumor": "tumour",
}


@dataclass
class Oligo:
    sequence: str                # sense/target-site sequence, 5'->3'
    antisense: str
    kind: str                    # "siRNA" | "gapmer ASO"
    position: int                # 0-based start in the transcript
    gc: float
    asymmetry: float             # >0 favours correct strand loading
    warnings: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


_COMP = str.maketrans("ACGTU", "TGCAA")


def revcomp(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


def _gc(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / max(len(s), 1)


# Nearest-neighbour-ish end stability proxy: AU-rich 5' antisense end favours
# loading of the correct strand into RISC.
def _asymmetry(target: str) -> float:
    """Positive = antisense 5' end is less stable than sense 5' end (good)."""
    def stab(s):
        return sum(2.0 if b in "GC" else 1.0 for b in s.upper())
    as_5 = revcomp(target)[:4]
    se_5 = target[:4]
    return round(stab(se_5) - stab(as_5), 2)


def design_sirna(transcript: str, *, n: int = 10, length: int = 21,
                 skip_5prime: int = 75) -> list[Oligo]:
    """Candidate siRNA target sites, by the standard sequence rules.

    Rules applied: GC 30-52%, no run of 4+ identical bases, thermodynamic
    asymmetry favouring antisense loading, avoidance of the first ~75 nt
    (ribosome-loaded, poorly accessible) and of immunostimulatory motifs.
    """
    seq = re.sub(r"[^ACGT]", "", transcript.upper())
    out: list[Oligo] = []

    for i in range(skip_5prime, len(seq) - length):
        site = seq[i:i + length]
        gc = _gc(site)
        o = Oligo(site, revcomp(site), "siRNA", i, round(gc, 3), _asymmetry(site))
        if not (0.30 <= gc <= 0.52):
            o.warnings.append(f"GC {gc:.0%} outside 30-52%")
        if re.search(r"(A{4,}|C{4,}|G{4,}|T{4,})", site):
            o.warnings.append("run of 4+ identical bases")
        if o.asymmetry <= 0:
            o.warnings.append("unfavourable thermodynamic asymmetry: the sense "
                              "strand may load into RISC instead")
        if "GTCCTTCAA" in site or "TGTGT" in site:
            o.warnings.append("contains a known immunostimulatory motif")
        out.append(o)

    out.sort(key=lambda o: (len(o.warnings), -o.asymmetry, abs(o.gc - 0.42)))
    return out[:n]


def design_gapmer(transcript: str, *, n: int = 10, length: int = 16,
                  skip_5prime: int = 75) -> list[Oligo]:
    """RNase H gapmer ASO sites (5-10-5 LNA/DNA/LNA is the usual architecture)."""
    seq = re.sub(r"[^ACGT]", "", transcript.upper())
    out: list[Oligo] = []

    for i in range(skip_5prime, len(seq) - length):
        site = seq[i:i + length]
        gc = _gc(site)
        o = Oligo(site, revcomp(site), "gapmer ASO", i, round(gc, 3), 0.0)
        if not (0.40 <= gc <= 0.60):
            o.warnings.append(f"GC {gc:.0%} outside 40-60%")
        if "CG" in site:
            o.warnings.append("contains CpG: immunostimulatory in "
                              "phosphorothioate backbones")
        if re.search(r"G{4,}", site):
            o.warnings.append("G-quadruplex risk (GGGG)")
        if re.search(r"(A{5,}|T{5,})", site):
            o.warnings.append("homopolymer run")
        out.append(o)

    out.sort(key=lambda o: (len(o.warnings), abs(o.gc - 0.50)))
    return out[:n]


def design_crispri(promoter_sequence: str, offset_fn, *, n: int = 8) -> list[dict]:
    """dCas9-KRAB guides. Different window from CRISPRa: centred on the TSS."""
    from .screens import find_guides
    guides = find_guides(promoter_sequence, offset_fn,
                         window=CRISPRI_WINDOW, max_guides=n * 3)
    out = []
    for g in guides:
        d = g.to_dict()
        d["in_optimal_window"] = CRISPRI_WINDOW[0] <= g.tss_offset <= CRISPRI_WINDOW[1]
        d["warnings"] = [w for w in g.warnings if "CRISPRa" not in w]
        if not d["in_optimal_window"]:
            d["warnings"].append(
                f"offset {g.tss_offset:+d} outside the CRISPRi window "
                f"({CRISPRI_WINDOW[0]} to {CRISPRI_WINDOW[1]}); KRAB silencing "
                "is most effective close to the TSS")
        out.append(d)
    out.sort(key=lambda d: (not d["in_optimal_window"], len(d["warnings"])))
    return out[:n]


def delivery_for(tissue: str | None) -> dict:
    key = TISSUE_ALIASES.get((tissue or "unknown").strip().lower(), "unknown")
    d = dict(DELIVERY[key])
    d["compartment"] = key
    return d


def design(gene: str, *, goal: str = "down", tissue: str | None = None,
           transcript: str | None = None, promoter=None) -> dict:
    """Propose direct modalities for this gene, gated on deliverability."""
    deliv = delivery_for(tissue)
    options: list[dict] = []

    if goal == "down":
        if transcript:
            options.append({
                "modality": "siRNA",
                "rationale": "catalytic, cytoplasmic; the best-validated "
                             "knockdown chemistry where delivery is solved",
                "candidates": [o.to_dict() for o in design_sirna(transcript)],
            })
            options.append({
                "modality": "gapmer ASO (RNase H)",
                "rationale": "works on nuclear and cytoplasmic RNA, so it can "
                             "hit pre-mRNA and nuclear-retained transcripts "
                             "that siRNA cannot",
                "candidates": [o.to_dict() for o in design_gapmer(transcript)],
            })
        if promoter is not None:
            options.append({
                "modality": "CRISPRi (dCas9-KRAB)",
                "rationale": "epigenetic silencing at the promoter; reversible, "
                             "and a fast way to test whether lowering this gene "
                             "does what you expect before committing to a "
                             "therapeutic modality",
                "candidates": design_crispri(promoter.sequence,
                                             promoter.offset_to_tss),
            })
        options.append({
            "modality": "nuclease knockout / base editing",
            "rationale": "permanent. For immune-cell targets this is the route "
                         "with human precedent — ex vivo editing of HSCs or T "
                         "cells followed by reinfusion, as done for CCR5.",
            "candidates": [],
            "note": "ex vivo only in practice; a transplant programme with "
                    "conditioning, not a drug",
        })
    else:
        options.append({
            "modality": "splice-switching oligonucleotide",
            "rationale": "the main way to RAISE output of a specific gene: "
                         "correct a splicing defect or force inclusion of an "
                         "exon, as nusinersen does for SMN2",
            "candidates": [],
            "note": "requires an exploitable splicing event — a pseudoexon, a "
                    "weak splice site, a skipped exon. Check the transcript "
                    "structure before assuming one exists.",
        })
        options.append({
            "modality": "CRISPRa (dCas9-VP64/SAM)",
            "rationale": "forces transcription from the endogenous promoter; "
                         "the cleanest test that raising this gene helps",
            "candidates": [],
            "note": "see the screen design section for guides",
        })
        options.append({
            "modality": "gene addition (AAV / lentiviral)",
            "rationale": "supply a working copy outright; the approved route "
                         "for monogenic loss of function",
            "candidates": [],
            "note": "AAV cargo limit ~4.7 kb; immunogenicity and pre-existing "
                    "neutralising antibodies exclude a substantial fraction of "
                    "patients",
        })

    feasible = deliv["score"] >= 0.5
    return {
        "gene": gene.upper(),
        "goal": goal,
        "delivery": deliv,
        "feasible_in_vivo": feasible,
        "options": options,
        "verdict": (
            f"Delivery to {deliv['compartment']} is {deliv['status']}. "
            + ("Design quality will be the limiting factor."
               if feasible else
               "DELIVERY, not design, is the blocker here. A potent oligo for "
               "this target is a cell-culture reagent until a delivery route "
               "exists — treat the ex vivo path as the realistic option.")),
        "off_target_warning":
            "Seed-region (guide positions 2-8) complementarity drives most "
            "siRNA off-target silencing and is a transcriptome-wide question. "
            "Run every candidate through BLAST against RefSeq, or a seed-match "
            "scan, before synthesising anything. The checks here are local "
            "sanity only.",
    }
