"""Raising a gene's output after transcription has already happened.

The gap this fills
------------------
The rest of this package is a transcriptional-initiation tool: promoter, then
transcription factor, then a drug acting on that factor. Everything between
"Pol II has started" and "functional protein exists" was invisible to it —
splicing, premature stop codons, transcript half-life, microRNA repression,
translational control. For some genes that interval is where the entire
problem lives, and no amount of extra transcription fixes it.

The clearest case: a nonsense allele. The promoter is fine, the factor is
bound, transcription is normal — and the message is destroyed by
nonsense-mediated decay before it is translated, or translated into a
truncated protein. Doubling transcription doubles the amount of RNA that gets
degraded.

Honesty about approval status
-----------------------------
Most of this class is **not FDA-approved**, and saying so is the point of
including it. A route listed here as `investigational` is not a treatment you
can use; it is a mechanism that exists, so that an empty transcriptional
result is not mistaken for "nothing can be done". Where the regulatory
picture differs between agencies — ataluren being the live example — both are
stated rather than the more favourable one.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class PostTranscriptionalRoute:
    mechanism: str
    what_it_does: str
    drugs: list
    status: str                  # "approved" | "approved elsewhere" | "investigational" | "none"
    applies_when: str
    raises: str                  # "transcript" | "functional protein" | "protein"
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


# --- Curated. Verify anything you intend to rely on; approvals move. --------
ROUTES: tuple = (
    PostTranscriptionalRoute(
        mechanism="nonsense readthrough",
        what_it_does=(
            "makes the ribosome read through a premature stop codon, so a "
            "full-length protein is produced from an allele that otherwise "
            "yields a truncated one"),
        drugs=["Gentamicin (FDA-approved as an antibiotic; readthrough use is "
               "off-label)",
               "Amikacin, tobramycin (same situation)",
               "Ataluren / Translarna (EMA conditional approval for nonsense "
               "DMD; the FDA did NOT approve it, and the EMA moved against "
               "renewal — do not record this as an approved option)",
               "ELX-02 (clinical, designed to separate readthrough from "
               "aminoglycoside toxicity)"],
        status="approved elsewhere",
        applies_when="the gene carries a premature termination codon (PTC)",
        raises="functional protein",
        caveats=[
            "Only helps a nonsense allele. A missense, frameshift or deletion "
            "allele gets nothing from readthrough.",
            "Aminoglycoside readthrough is weak and the doses implied are "
            "nephrotoxic and ototoxic; this is why the class has not become a "
            "therapy despite thirty years of interest.",
            "Readthrough inserts a near-cognate amino acid at the stop, so "
            "the product is full-length but not necessarily wild-type.",
            "Efficiency depends strongly on the stop codon identity and its "
            "sequence context — UGA with a C at +4 reads through best.",
        ]),
    PostTranscriptionalRoute(
        mechanism="nonsense-mediated decay (NMD) inhibition",
        what_it_does=(
            "stops the surveillance pathway destroying a PTC-bearing "
            "transcript, so more message survives to be translated"),
        drugs=["No approved NMD inhibitor exists.",
               "SMG1 and UPF1 inhibitors are preclinical."],
        status="investigational",
        applies_when="the gene carries a PTC and the transcript is being "
                     "degraded before translation",
        raises="transcript",
        caveats=[
            "NMD is a general surveillance pathway: blocking it stabilises "
            "every PTC-bearing transcript in the cell, not only yours, and "
            "roughly a third of normal transcripts are NMD substrates.",
            "Raising a PTC transcript without readthrough gives you more "
            "truncated protein, which may be inert or actively harmful. The "
            "two mechanisms are usually only interesting together.",
        ]),
    PostTranscriptionalRoute(
        mechanism="splice correction",
        what_it_does=(
            "redirects splicing so a productive transcript is made — exon "
            "inclusion, exon skipping, or suppression of a cryptic site"),
        drugs=["Risdiplam (FDA-approved, SMA — small molecule, oral)",
               "Nusinersen (FDA-approved, SMA — intrathecal ASO)",
               "Eteplirsen, golodirsen, viltolarsen, casimersen (FDA "
               "accelerated approval, DMD exon skipping)",
               "Branaplam (clinical; SMN2 splice modifier)"],
        status="approved",
        applies_when="a splice-site or splice-regulatory variant, a "
                     "poison exon, or a paralog whose splicing can be "
                     "redirected (SMN2)",
        raises="functional protein",
        caveats=[
            "Exon-skipping restores a partly functional, internally deleted "
            "protein, not a wild-type one — the DMD agents are judged on "
            "dystrophin that is shorter than normal.",
            "Accelerated approvals in this class rest on a surrogate "
            "endpoint (protein level), and the confirmatory functional "
            "benefit has been contested.",
            "Small-molecule splice modifiers are not fully gene-specific; "
            "risdiplam alters splicing of other transcripts too.",
        ]),
    PostTranscriptionalRoute(
        mechanism="transcript stabilisation",
        what_it_does=(
            "lengthens mRNA half-life, usually through AU-rich elements in "
            "the 3'UTR and the proteins that bind them. Doubling half-life "
            "doubles steady-state message with no change in transcription"),
        drugs=["No approved agent targets a specific transcript's stability.",
               "HuR/ELAVL1 modulators are preclinical."],
        status="investigational",
        applies_when="the transcript is short-lived and the promoter is "
                     "already active — the case where more transcription "
                     "achieves nothing",
        raises="transcript",
        caveats=[
            "Entirely invisible to a promoter-level analysis, which is why it "
            "is worth checking before concluding the transcriptional route is "
            "the only one.",
            "RNA-binding proteins are highly pleiotropic; HuR alone governs "
            "hundreds of transcripts.",
        ]),
    PostTranscriptionalRoute(
        mechanism="microRNA de-repression",
        what_it_does=(
            "blocks a microRNA that represses your transcript, raising both "
            "message and protein"),
        drugs=["No approved antagomir.",
               "Miravirsen and RG-101 (anti-miR-122) reached clinical trials."],
        status="investigational",
        applies_when="a specific miRNA is a documented repressor of your gene "
                     "in your cell type",
        raises="transcript",
        caveats=[
            "One miRNA has hundreds of targets, so this is a low-specificity "
            "intervention despite being sequence-directed.",
            "Needs the repressing miRNA identified in the right cell type; a "
            "TargetScan prediction is not sufficient evidence.",
        ]),
    PostTranscriptionalRoute(
        mechanism="translational enhancement",
        what_it_does=(
            "raises protein output from existing message, usually by "
            "relieving a uORF or modulating eIF2-alpha phosphorylation"),
        drugs=["No approved agent raises translation of a chosen transcript.",
               "ISRIB and related integrated-stress-response modulators are "
               "preclinical."],
        status="investigational",
        applies_when="protein is low while transcript is normal — the "
                     "signature of a translational, not transcriptional, "
                     "problem",
        raises="protein",
        caveats=[
            "Global translational modulators move the whole proteome.",
            "An RT-qPCR cannot see this problem OR this solution; you need "
            "protein and, ideally, ribosome profiling.",
        ]),
)


def routes_for(goal: str = "up", *, has_ptc: bool | None = None) -> list:
    """The post-transcriptional routes worth considering.

    `has_ptc` reorders rather than filters: with a known premature stop,
    readthrough and NMD move to the front, because for that allele they are
    the mechanisms that matter and extra transcription is not.
    """
    if goal != "up":
        # Lowering a gene post-transcriptionally is the ASO/siRNA story, and
        # `modalities` already designs those properly.
        return []
    # Default order: broadly applicable and approved first. Readthrough and
    # NMD are allele-specific - they do nothing for a missense, frameshift
    # or deletion allele - so leading with them would be wrong for most
    # genes. With a known PTC they become the mechanisms that matter, and
    # extra transcription becomes the thing that does not help.
    ALLELE_SPECIFIC = ("nonsense readthrough",
                       "nonsense-mediated decay (NMD) inhibition")
    APPROVED_FIRST = {"approved": 0, "approved elsewhere": 1,
                      "investigational": 2, "none": 3}

    def key(r):
        if has_ptc:
            return (r.mechanism not in ALLELE_SPECIFIC,
                    APPROVED_FIRST.get(r.status, 9))
        return (r.mechanism in ALLELE_SPECIFIC,
                APPROVED_FIRST.get(r.status, 9))

    return sorted(ROUTES, key=key)


def summarise(routes: list, *, has_ptc: bool | None = None) -> str:
    if not routes:
        return ""
    approved = [r for r in routes if r.status == "approved"]
    lines = [
        "A promoter-level analysis stops at transcription. These act after "
        "it, and for some genes that is where the whole problem is: if the "
        "transcript is destroyed or mistranslated, doubling transcription "
        "doubles the amount that is destroyed."]
    if has_ptc:
        lines.append(
            "**Your gene has a premature stop codon**, so readthrough and NMD "
            "inhibition are listed first — they address that allele directly, "
            "and no amount of extra transcription does.")
    lines.append(
        f"{len(approved)} of {len(routes)} mechanism(s) here have an "
        f"FDA-approved agent. The rest are listed so that an empty "
        f"transcriptional result is not mistaken for 'nothing can be done' — "
        f"they are mechanisms, not treatments.")
    return " ".join(lines)
