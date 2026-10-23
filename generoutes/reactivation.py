"""Reactivating a silenced allele, which is a different problem from a low gene.

The distinction this draws
-------------------------
"Raise this gene" collapses two situations that need opposite interventions:

**Low expression.** The promoter is active but weakly. More transcription
helps, and the rest of this package is about finding a way to get it.

**Silenced allele.** One copy is intact and switched off — by imprinting, by
a developmental switch, by X-inactivation. You do not need more
transcription; you need the brake released. Pushing harder on the active copy
is the wrong move, and often there is no active copy to push.

The second case has produced the field's clearest success and nothing in this
package could previously see it. Hydroxyurea does not touch HBB. It raises
HBG1/HBG2, the fetal globins, which were silenced after birth — and that is
an FDA-approved treatment for sickle cell disease, working entirely by
reactivating something the cell already knew how to make.

Why reactivation is easier than induction
-----------------------------------------
A silenced allele has an intact promoter, intact coding sequence and a cell
that has transcribed it before. The machinery exists; a repressive mark is
holding it shut. Releasing a brake is a smaller ask than building expression
where there was none, which is why this class has clinical successes while
generic "upregulate my gene" largely does not.

The cost, in every case
-----------------------
**These mechanisms are not gene-specific.** A DNA-methylation inhibitor
reactivates every methylation-silenced locus in the genome. A topoisomerase
inhibitor is a cytotoxic chemotherapeutic. The specificity has to come from
somewhere else — the dose, the tissue, or the fact that only your locus has
the relevant mark — and if it cannot, the route is not usable however
appealing the mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class ReactivationRoute:
    mechanism: str
    silencing_type: str
    what_it_does: str
    drugs: list
    status: str                  # "approved" | "clinical" | "preclinical"
    worked_example: str
    specificity: str             # where selectivity comes from, if anywhere
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


# --- Curated. Verify anything you intend to rely on. -----------------------
ROUTES: tuple = (
    ReactivationRoute(
        mechanism="fetal-gene reactivation (developmental switch reversal)",
        silencing_type="developmental silencing",
        what_it_does=(
            "re-expresses a paralog that was switched off after birth, by "
            "interfering with the repressor that closed it"),
        drugs=["Hydroxyurea (FDA-approved, sickle cell) — raises HbF",
               "Voxelotor, L-glutamine (approved for sickle cell by other "
               "mechanisms; not HbF inducers)",
               "Exagamglogene autotemcel / exa-cel (FDA-approved gene "
               "editing that disrupts the BCL11A erythroid enhancer — the "
               "switch itself)",
               "Decitabine, azacitidine (approved hypomethylating agents; "
               "raise HbF but are cytotoxic)"],
        status="approved",
        worked_example=(
            "sickle cell disease: HBB is mutated and undruggable, so the "
            "intervention raises HBG1/HBG2 instead. The clearest proof that "
            "reactivating a silenced paralog beats fixing the broken gene"),
        specificity=(
            "comes from the tissue and the repressor, not the drug: BCL11A's "
            "erythroid enhancer is lineage-restricted, which is what makes "
            "exa-cel targeted rather than global"),
        caveats=[
            "Hydroxyurea's HbF mechanism is still incompletely defined after "
            "decades of clinical use, and response varies widely between "
            "patients.",
            "Hypomethylating agents reactivate the whole methylome, not your "
            "locus; they are chemotherapeutics.",
        ]),
    ReactivationRoute(
        mechanism="imprinted-allele unsilencing",
        silencing_type="genomic imprinting",
        what_it_does=(
            "re-expresses the silent parental copy of an imprinted gene, "
            "leaving the already-expressed copy alone"),
        drugs=["Topotecan and other topoisomerase-I inhibitors "
               "(FDA-approved in oncology; this use is preclinical "
               "repurposing, not an approved indication)",
               "Antisense oligonucleotides against the silencing transcript "
               "(clinical for Angelman: rugonersen, ION582, GTX-102)"],
        status="clinical",
        worked_example=(
            "Angelman syndrome: maternal UBE3A is mutated and the paternal "
            "copy is imprinted silent. Topoisomerase inhibitors unsilence "
            "the paternal allele in mouse models by interfering with the "
            "long antisense transcript that represses it"),
        specificity=(
            "the ASO route is genuinely sequence-specific. The small-molecule "
            "route is not — a topoisomerase inhibitor is a cytotoxic drug "
            "and the doses are not plausibly tolerable for a neurological "
            "indication"),
        caveats=[
            "Only applies to imprinted loci — roughly 150 genes in humans. "
            "Check whether yours is one before reading further.",
            "The silent allele must be intact. If both copies carry the "
            "mutation, unsilencing achieves nothing.",
            "Most of the mouse data is developmental; the window in which "
            "reactivation helps may close.",
        ]),
    ReactivationRoute(
        mechanism="X-reactivation",
        silencing_type="X-chromosome inactivation",
        what_it_does=(
            "re-expresses the gene from the inactive X in a heterozygous "
            "female, where roughly half of cells silence the healthy copy"),
        drugs=["No approved agent.",
               "XIST-directed ASOs and combined DNMT/EZH2 inhibition are "
               "preclinical."],
        status="preclinical",
        worked_example=(
            "Rett syndrome: MECP2 is X-linked, and in heterozygous females "
            "the wild-type copy sits on the inactive X in half of neurons. "
            "Reactivating it would be curative in those cells"),
        specificity=(
            "essentially none at present. X-reactivation is genome-scale by "
            "definition, and reactivating the whole inactive X risks "
            "overexpressing hundreds of genes"),
        caveats=[
            "Only relevant for X-linked genes in heterozygous females.",
            "MECP2 in particular is dose-sensitive in both directions — "
            "duplication causes its own syndrome, so overshooting is harmful "
            "rather than merely wasteful.",
        ]),
    ReactivationRoute(
        mechanism="promoter de-methylation",
        silencing_type="DNA methylation",
        what_it_does=(
            "removes the methylation holding a CpG-island promoter shut, so "
            "the existing machinery can transcribe it again"),
        drugs=["Azacitidine, decitabine (FDA-approved, MDS/AML)",
               "Cedazuridine/decitabine (approved oral combination)"],
        status="approved",
        worked_example=(
            "tumour-suppressor reactivation in myelodysplastic syndrome, "
            "where the therapeutic effect is attributed partly to "
            "re-expressing methylation-silenced genes"),
        specificity=(
            "none at the locus level. These reactivate every "
            "methylation-silenced promoter in the genome, which is why they "
            "are useful in malignancy and hard to justify otherwise"),
        caveats=[
            "Myelosuppressive, and the global reactivation includes "
            "retroelements and genes you did not want.",
            "Only works if your promoter is actually methylated — worth "
            "measuring before assuming, and this package does not look at "
            "methylation at all.",
            "Already represented in the chromatin-level category, where its "
            "specificity is scored at ~0.05. Listed here because the "
            "*reason* to use it is different: a specific silenced allele, "
            "not a general nudge.",
        ]),
)


# Human imprinted loci where an allele-specific question is worth asking.
# Short and incomplete; the point is to prompt the check, not to be a registry.
IMPRINTED = frozenset({
    "UBE3A", "SNRPN", "MKRN3", "MAGEL2", "NDN", "SNORD116",
    "IGF2", "H19", "CDKN1C", "KCNQ1OT1", "PHLDA2", "SLC67A1",
    "GNAS", "NESP", "PLAGL1", "GRB10", "MEST", "PEG3", "PEG10",
    "DLK1", "MEG3", "RTL1", "ZIM2", "NNAT", "BLCAP", "INPP5F",
})

# Genes with a documented developmental-switch paralog, where the fetal-gene
# route is the first thing to ask about.
SWITCHED = {
    "HBB": "HBG1/HBG2 (fetal globin, silenced after birth by BCL11A)",
    "HBD": "HBG1/HBG2",
    "MYH7": "MYH6 (the cardiac alpha/beta myosin switch)",
    "DMD": "UTRN (utrophin, high in fetal and regenerating muscle)",
}


def assess(gene: str, *, goal: str = "up") -> dict:
    """Is this a silenced-allele problem rather than a low-expression one?"""
    if goal != "up":
        return {}

    sym = (gene or "").strip().upper()
    result = {"gene": sym, "routes": [], "notes": [], "applicable": False}

    if sym in IMPRINTED:
        result["applicable"] = True
        result["notes"].append(
            f"**{sym} is an imprinted locus.** One parental copy is silenced "
            f"by design, so if the expressed copy is the mutated one, the "
            f"question is not how to transcribe more — it is how to release "
            f"the intact silent copy. That is a different intervention, and "
            f"it is the Angelman/UBE3A situation.")

    if sym in SWITCHED:
        result["applicable"] = True
        result["notes"].append(
            f"**{sym} has a developmental-switch paralog: {SWITCHED[sym]}.** "
            f"Reactivating it may be easier than raising {sym}, because the "
            f"cell has transcribed it before and a repressor is holding it "
            f"shut. This is the hydroxyurea/HbF route, and it is "
            f"FDA-approved for sickle cell.")

    # The mechanism list is shown regardless: silencing is not something this
    # package can detect, so the user has to be told the question exists.
    result["routes"] = [r.to_dict() for r in ROUTES]
    if not result["applicable"]:
        result["notes"].append(
            f"{sym} is not in the curated imprinted or developmental-switch "
            f"lists, which are short. That is not evidence the gene is "
            f"unsilenced — **nothing in this package measures methylation, "
            f"imprinting or allele-specific expression.** If your gene is "
            f"silenced rather than merely low, every transcriptional route "
            f"above is aimed at the wrong problem, and the way to find out "
            f"is allele-specific RT-qPCR or bisulphite sequencing.")
    return result


def summarise(result: dict) -> str:
    if not result:
        return ""
    if result.get("applicable"):
        return (f"**{result['gene']} may be a silenced-allele problem rather "
                f"than a low-expression one.** If so, the intervention is to "
                f"release a brake on an intact copy, not to drive the "
                f"promoter harder — a different and historically more "
                f"successful class of fix.")
    return ("Whether this gene is silenced rather than merely low is "
            "unmeasured here, and it changes which routes make sense. "
            "Allele-specific expression is the experiment that settles it.")
