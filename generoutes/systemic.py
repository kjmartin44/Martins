"""Drugs and receptors whose action is systemic by design, not by accident.

The gap this closes
-------------------
Dexamethasone is a clean, potent, highly selective NR3C1 agonist. It passes
*every* filter this package already has: good potency, good protein
selectivity, no withdrawal flag, an unambiguous agonist annotation. It is also
a whole-body drug, because acting on nearly every tissue is what the
glucocorticoid receptor does.

`selectivity` counts how many **proteins** a drug binds. `regulon` counts how
many **genes** that protein controls. Neither counts **systems**, and a
compound can score perfectly on both while reorganising a patient's
metabolism, immune function and bone. This module is the third denominator.

Three independent signals
-------------------------
**Drug class.** A curated table of pharmacological classes that are systemic
by mechanism — glucocorticoids, retinoids, thyroid hormones, sex steroids.

**Target class.** Reused from `regulon.GLOBAL_REGULATORS` /
`BROAD_REGULATORS` rather than duplicated here: if the *receptor* is a global
regulator, the route is systemic even when the drug itself is pristine. This
is the signal that catches dexamethasone.

**Regulatory flags.** Withdrawal and black-box warnings, read from the *same*
note strings `shortlist` already filters on, so one drug cannot be withdrawn
for one filter and clean for another.

What this table is not
----------------------
It is hand-entered, short, and certainly incomplete. **Membership is evidence
of systemic burden; non-membership is not evidence of its absence.** A drug
absent from these sets is *unassessed*, and `ClinicalBurden.assessed` says so
rather than reporting it as clean.

**Route of administration is not in the underlying data.** Topical clobetasol
and systemic dexamethasone are both "glucocorticoid" here, and this module
cannot tell them apart. A flag is a prompt to check the label, never a
verdict — which is why flagged drugs are demoted and shown, never dropped.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .congeners import salt_free_name


@dataclass(frozen=True)
class SystemicClass:
    name: str
    burden: float          # 0..1, HIGHER = more systemic
    why: str               # the specific burden, not a generic warning
    members: frozenset


CLASSES: tuple[SystemicClass, ...] = (
    SystemicClass(
        "glucocorticoid", 0.90,
        "acts on nearly every tissue; chronic exposure brings adrenal "
        "suppression, hyperglycaemia, bone loss and immunosuppression",
        frozenset({
            "dexamethasone", "prednisone", "prednisolone",
            "methylprednisolone", "hydrocortisone", "cortisone",
            "betamethasone", "budesonide", "fluticasone", "triamcinolone",
            "beclomethasone", "beclometasone", "mometasone", "ciclesonide",
            "deflazacort", "clobetasol", "halobetasol", "desonide",
            "flunisolide", "halcinonide", "flurandrenolide", "fluocinolone",
            "fluocinonide", "diflorasone", "amcinonide",
        })),
    SystemicClass(
        "mineralocorticoid", 0.75,
        "sodium and volume handling; hypertension and hypokalaemia",
        frozenset({"fludrocortisone", "aldosterone", "desoxycorticosterone"})),
    SystemicClass(
        "retinoid", 0.85,
        "RAR/RXR drive differentiation programmes genome-wide; teratogenic, "
        "with mucocutaneous, hepatic and lipid toxicity across every trial",
        frozenset({
            "tretinoin", "isotretinoin", "alitretinoin", "acitretin",
            "etretinate", "bexarotene", "tazarotene", "adapalene",
            "trifarotene", "tamibarotene", "palovarotene",
        })),
    SystemicClass(
        "thyroid hormone", 0.85,
        "sets basal metabolic rate in every tissue; over-replacement causes "
        "atrial fibrillation and bone loss",
        frozenset({
            "levothyroxine", "liothyronine", "liotrix", "dextrothyroxine",
            "thyroid", "thyroglobulin",
        })),
    SystemicClass(
        "sex steroid", 0.80,
        "reproductive, skeletal, hepatic and cardiovascular effects well "
        "beyond the target tissue; thrombotic and hormone-sensitive-cancer risk",
        frozenset({
            "estradiol", "ethinylestradiol", "ethinyl estradiol", "estrone",
            "estriol", "mestranol", "diethylstilbestrol", "conjugated estrogens",
            "testosterone", "methyltestosterone", "oxandrolone", "nandrolone",
            "danazol", "progesterone", "medroxyprogesterone", "norethindrone",
            "norethisterone", "levonorgestrel", "drospirenone", "etonogestrel",
            "megestrol", "chlormadinone", "ulipristal", "mifepristone",
            "raloxifene", "tamoxifen", "fulvestrant", "elacestrant",
            "bicalutamide", "flutamide", "nilutamide", "enzalutamide",
            "abiraterone", "finasteride", "dutasteride",
        })),
    SystemicClass(
        "vitamin D receptor agonist", 0.65,
        "raises serum calcium body-wide; hypercalcaemia, nephrolithiasis and "
        "vascular calcification are the dose-limiting effects",
        frozenset({
            "calcitriol", "alfacalcidol", "paricalcitol", "doxercalciferol",
            "calcipotriene", "calcipotriol", "ergocalciferol",
            "cholecalciferol", "calcifediol", "dihydrotachysterol",
        })),
    SystemicClass(
        "PPAR agonist", 0.70,
        "systemic lipid and glucose handling; fluid retention, weight gain "
        "and heart-failure risk in the thiazolidinedione class",
        frozenset({
            "pioglitazone", "rosiglitazone", "fenofibrate", "bezafibrate",
            "gemfibrozil", "clofibrate", "elafibranor", "saroglitazar",
            "seladelpar", "pemafibrate",
        })),
    SystemicClass(
        "cytotoxic / antiproliferative", 0.95,
        "kills dividing cells wherever it reaches them; myelosuppression and "
        "mucositis are on-mechanism, not off-target",
        frozenset({
            # Anthracyclines and the anthracenedione. Mitoxantrone was
            # absent and is one of the compounds that crowds a real
            # shortlist through LINCS abundance signatures.
            "doxorubicin", "daunorubicin", "epirubicin", "idarubicin",
            "mitoxantrone",
            # Tubulin
            "paclitaxel", "docetaxel", "cabazitaxel", "vincristine",
            "vinblastine", "vinorelbine", "vindesine", "eribulin",
            "ixabepilone",
            # Topoisomerase
            "etoposide", "topotecan", "irinotecan",
            # Platinums
            "cisplatin", "carboplatin", "oxaliplatin",
            # Antimetabolites
            "methotrexate", "pralatrexate", "raltitrexed", "pemetrexed",
            "fluorouracil", "capecitabine", "gemcitabine", "cytarabine",
            "cladribine", "fludarabine", "clofarabine", "nelarabine",
            "mercaptopurine", "thioguanine", "hydroxyurea",
            "hydroxycarbamide",
            # Alkylating agents
            "cyclophosphamide", "ifosfamide", "melphalan", "busulfan",
            "chlorambucil", "bendamustine", "carmustine", "lomustine",
            "dacarbazine", "temozolomide", "procarbazine", "thiotepa",
            "streptozocin",
            # Antitumour antibiotics and other DNA-directed agents
            "dactinomycin", "bleomycin", "mitomycin", "plicamycin",
            "mithramycin", "trabectedin",
            # Proteasome inhibitors: not classical cytotoxics, but the
            # same on-mechanism marrow toxicity puts them here rather
            # than with the kinase inhibitors.
            "bortezomib", "carfilzomib", "ixazomib",
            # Other
            "arsenic trioxide", "asparaginase",
        })),
    SystemicClass(
        "antineoplastic kinase inhibitor", 0.80,
        "dosed to stop cells dividing, so rash, diarrhoea, myelosuppression "
        "and hepatotoxicity are on-mechanism; and because proliferation "
        "signalling runs through every tissue, an inhibitor moves the "
        "abundance of a great many transcription factors at once - which is "
        "why these crowd an expression-signature shortlist without being "
        "evidence of a route to any one gene",
        frozenset({
            # BCR-ABL
            "imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib",
            # EGFR and HER2
            "erlotinib", "gefitinib", "afatinib", "osimertinib",
            "dacomitinib", "lapatinib", "neratinib", "tucatinib",
            # ALK and ROS1
            "crizotinib", "ceritinib", "alectinib", "brigatinib",
            "lorlatinib",
            # BRAF and MEK
            "vemurafenib", "dabrafenib", "encorafenib", "trametinib",
            "cobimetinib", "binimetinib", "selumetinib",
            # VEGFR and multi-kinase
            "sunitinib", "sorafenib", "pazopanib", "axitinib",
            "regorafenib", "lenvatinib", "cabozantinib", "vandetanib",
            "tivozanib",
            # BTK
            "ibrutinib", "acalabrutinib", "zanubrutinib",
            # CDK4/6
            "palbociclib", "ribociclib", "abemaciclib",
            # PI3K
            "idelalisib", "duvelisib", "alpelisib", "copanlisib",
            # FLT3, KIT, RET, TRK, FGFR
            "midostaurin", "gilteritinib", "quizartinib", "ripretinib",
            "avapritinib", "selpercatinib", "pralsetinib", "entrectinib",
            "larotrectinib", "pemigatinib", "erdafitinib", "infigratinib",
            # JAK, in its oncology indications
            "ruxolitinib", "fedratinib",
            # Investigational and tool kinase inhibitors. These dominate
            # the LINCS perturbation library and therefore these
            # shortlists, and every name below appeared in a real run.
            # Membership is a claim about PHARMACOLOGY, not approval -
            # whether a compound is a drug at all is what
            # `catalogue.verify_approval` answers, separately.
            "dinaciclib", "seliciclib", "roscovitine", "alvocidib",
            "flavopiridol", "saracatinib", "bms-387032", "sns-032",
            "at-7519", "cgp-60474", "pd-0325901", "alw-ii-38-3",
            "as-605240", "gsk-1059615", "jnk-9l", "raf265",
            "nvp-bez235", "bi-2536", "tozasertib", "danusertib",
            "dorsomorphin",
        })),
    SystemicClass(
        "broad immunosuppressant", 0.85,
        "global suppression of adaptive immunity; infection and malignancy "
        "risk scale with exposure",
        frozenset({
            "cyclosporine", "ciclosporin", "tacrolimus", "sirolimus",
            "everolimus", "mycophenolate", "azathioprine", "leflunomide",
        })),
    SystemicClass(
        "chromatin-level agent", 0.95,
        "moves thousands of promoters at once; systemic by mechanism rather "
        "than by distribution",
        frozenset({
            "vorinostat", "romidepsin", "belinostat", "panobinostat",
            "azacitidine", "decitabine", "tazemetostat", "valproic acid",
            "valproate", "entinostat",
        })),
)

# Built once. Members are already lowercase and salt-free by construction.
_INDEX: dict[str, SystemicClass] = {m: c for c in CLASSES for m in c.members}
# The same table keyed without spaces or hyphens, so a code name
# written "RAF 265" in one source and "RAF265" in another resolves
# to the same class instead of silently missing.
_COMPACT: dict[str, SystemicClass] = {
    m.replace(" ", "").replace("-", ""): c
    for c in CLASSES for m in c.members}

BURDEN_LEVELS = ("low", "moderate", "high", "prohibitive")

# Compounds that disturb mitosis directly: they bind tubulin, or inhibit the
# kinases or motors that run the spindle. This is a TAG, not a class, because
# the same compound can also belong to a class above (paclitaxel is
# cytotoxic, BI-2536 a kinase inhibitor) and the classes must stay disjoint.
#
# It matters for one kind of target. For a gene whose job is mitosis, a
# compound that stalls or scrambles mitosis changes the cell's cycle state and
# its chromosome count, which is the wrong kind of evidence and the wrong kind
# of drug. Hand-entered and short: absence is not evidence that a compound
# leaves mitosis alone.
MITOTIC_DISRUPTORS = frozenset({
    # tubulin binders
    "paclitaxel", "docetaxel", "cabazitaxel", "vincristine", "vinblastine",
    "vinorelbine", "vinflunine", "eribulin", "ixabepilone", "colchicine",
    "colcemid", "demecolcine", "podophyllotoxin", "nocodazole",
    "mebendazole", "albendazole", "combretastatin",
    # polo-like kinase 1
    "bi-2536", "volasertib", "gsk-461364", "onvansertib", "rigosertib",
    # aurora kinases
    "alisertib", "mln8237", "barasertib", "azd1152", "danusertib",
    "tozasertib", "vx-680", "mk-0457", "zm-447439", "hesperadin",
    "mk-5108", "enmd-2076", "cenisertib",
    # kinesin Eg5, CDK1, MPS1
    "monastrol", "ispinesib", "filanesib", "ro-3306", "reversine",
    "nms-p715",
    # tubulin binders that CTD records raising BUB1B protein, added because a
    # report of that kind can be mitotic arrest and not the promoter
    "benomyl", "carbendazim", "rotenone",
    "2'-fluoro-6,7-methylenedioxy-2-phenyl-4-quinolone",
})
_MITOTIC_COMPACT = frozenset(m.replace(" ", "").replace("-", "").replace(",", "")
                             for m in MITOTIC_DISRUPTORS)


def _compact(name: str) -> str:
    """A name with spaces, hyphens and commas removed, for code-name matching."""
    return name.replace(" ", "").replace("-", "").replace(",", "")


def is_mitotic_disruptor(drug_name: str) -> bool:
    """Is this a compound that disturbs mitosis directly?

    Matched on the salt-free name with spaces and hyphens removed, the way
    `classify` matches code names, so "BI 2536" and "BI-2536" are one.
    """
    if not drug_name:
        return False
    return _compact(salt_free_name(drug_name)) in _MITOTIC_COMPACT


# Compounds whose therapeutic mechanism IS damage to DNA or to the
# machinery that segregates it. Separate from MITOTIC_DISRUPTORS, which is
# about a confound in the readout: this is about the patient.
#
# Not a hazard grading and not exhaustive. It is the approved classes
# whose labels carry a mutagenicity or secondary-malignancy warning, which
# is the bar that matters for someone whose cells already miss-segregate
# chromosomes. A compound absent from it is unexamined, not cleared.
GENOTOXIC = frozenset({
    # topoisomerase poisons
    "doxorubicin", "daunorubicin", "epirubicin", "idarubicin",
    "mitoxantrone", "etoposide", "teniposide", "topotecan", "irinotecan",
    # alkylators and platinum
    "cyclophosphamide", "ifosfamide", "chlorambucil", "melphalan",
    "busulfan", "lomustine", "carmustine", "temozolomide", "dacarbazine",
    "procarbazine", "bendamustine", "thiotepa", "mechlorethamine",
    "cisplatin", "carboplatin", "oxaliplatin",
    # antimetabolites incorporated into DNA
    "azacitidine", "decitabine", "cytarabine", "gemcitabine",
    "fluorouracil", "capecitabine", "mercaptopurine", "thioguanine",
    "azathioprine", "methotrexate", "hydroxyurea", "cladribine",
    "fludarabine", "clofarabine", "nelarabine",
    # radiomimetic and other direct DNA damage
    "bleomycin", "mitomycin", "dactinomycin", "plicamycin",
    "trabectedin", "lurbinectedin",
})

_GENOTOXIC_COMPACT = frozenset(_compact(m) for m in GENOTOXIC)


# Genes where moving the dose one way is a cancer risk in itself.
#
# Short and conservative on purpose. These are transcription factors and
# regulators whose role in human cancer is not in dispute, because the
# point is to catch an obviously bad idea and not to adjudicate the
# literature. Most genes are context-dependent and absent from both lists,
# and absence means unexamined.
#
# The upgrade, if this earns its keep, is the COSMIC Cancer Gene Census
# rather than a longer hand-written list.
ONCOGENES = frozenset("""
MYC MYCN MYCL MYB MYBL2 FOS FOSL1 JUN JUNB ETS1 ETS2 ERG FLI1 ETV1 ETV4
ETV6 SPI1 RELA RELB NFKB1 STAT3 STAT5A E2F1 E2F3 FOXM1 CCND1 CCNE1 CTNNB1
NOTCH1 GLI1 GLI2 YAP1 WWTR1 TWIST1 SNAI1 SNAI2 ZEB1 LMO2 TAL1 LYL1 HOXA9
MEIS1 PAX3 PAX7 SOX2 KLF4 NANOG MDM2 MDM4 BCL2 BCL6 EZH2 BMI1 SKP2 SRC
ABL1 KRAS HRAS NRAS PIK3CA AKT1 MTOR HIF1A
""".split())

# And the mirror: lowering one of these raises cancer risk the same way.
TUMOUR_SUPPRESSORS = frozenset("""
TP53 TP63 TP73 RB1 PTEN APC VHL NF1 NF2 BRCA1 BRCA2 ATM ATR CHEK2 CDKN1A
CDKN1B CDKN2A SMAD2 SMAD3 SMAD4 STK11 WT1 RUNX3 FOXO1 FOXO3 FOXP1 MXD1
MXI1 MNT HBP1 PML BAP1 ARID1A SMARCB1 KEAP1 TSC1 TSC2
""".split())


def cancer_risk_direction(symbol: str, want: str | None) -> str | None:
    """Why moving this gene `want` is a cancer risk, or None.

    `want` is the direction the route needs the gene to move, as
    `abundance.wanted_direction` computes it. None when the role is not
    recorded, and then no claim is made: guessing a direction here would
    invent the very thing the caution rests on.
    """
    if not symbol or want not in ("up", "down"):
        return None
    sym = symbol.strip().upper()
    if want == "up" and sym in ONCOGENES:
        return f"{sym} is a proto-oncogene and this route raises it"
    if want == "down" and sym in TUMOUR_SUPPRESSORS:
        return f"{sym} is a tumour suppressor and this route lowers it"
    return None


def is_genotoxic(drug_name: str) -> bool:
    """Does this compound work by damaging DNA or chromosome segregation?

    True also for the tubulin and spindle compounds in
    `MITOTIC_DISRUPTORS`: an aneugen gives a cell the wrong number of
    chromosomes, which is the very defect in a chromosomal-instability
    syndrome.
    """
    if not drug_name:
        return False
    key = _compact(salt_free_name(drug_name))
    return key in _GENOTOXIC_COMPACT or key in _MITOTIC_COMPACT

def classify(drug_name: str) -> SystemicClass | None:
    """Which systemic class this drug belongs to, or None if unlisted.

    Matches on the salt-free name so "Levothyroxine Sodium" resolves — the
    same normalisation `pipeline` uses when collapsing salt forms, so the two
    cannot disagree about what a drug is.
    """
    if not drug_name:
        return None
    name = salt_free_name(drug_name)
    hit = _INDEX.get(name)
    if hit:
        return hit

    # Code names are spelled inconsistently across sources: LINCS writes
    # "RAF 265" where ChEMBL writes "RAF265", and the head-token fallback
    # below is wrong for them - it would look up "raf" and miss, because
    # for a code name the number IS part of the identity. Compare with
    # spaces and hyphens removed before falling back.
    compact = _compact(name)
    hit = _COMPACT.get(compact)
    if hit:
        return hit

    # Multi-word brand-ish names: fall back to the head token, which is the
    # active moiety ("Fluticasone Propionate HFA" -> "fluticasone"). Only
    # when the head is a word rather than a code stem, so "raf" cannot
    # match something it is merely a prefix of.
    parts = name.split()
    head = parts[0] if parts else ""
    if len(parts) > 1 and not any(ch.isdigit() for ch in name):
        return _INDEX.get(head)
    return None


@dataclass
class ClinicalBurden:
    level: str = "low"
    score: float = 0.0
    reasons: list = field(default_factory=list)
    assessed: bool = False
    drug_class: str | None = None

    def to_dict(self) -> dict:
        return {"level": self.level, "score": round(self.score, 3),
                "reasons": self.reasons, "assessed": self.assessed,
                "drug_class": self.drug_class}


def _level_for(score: float, withdrawn: bool) -> str:
    if withdrawn:
        return "prohibitive"
    if score >= 0.80:
        return "high"
    if score >= 0.50:
        return "moderate"
    return "low"


def clinical_burden(drug: str, *, notes: str = "",
                    target: str | None = None) -> ClinicalBurden:
    """How much of the patient does this route move, beyond the target gene?

    `notes` is the joined route + grade note string, read for the SAME
    "WITHDRAWN" / "black-box" markers `shortlist` filters on.
    `target` is the route's `acts_on` symbol, checked against the regulator
    tiers already curated in `regulon`.
    """
    from . import regulon as regulon_mod

    b = ClinicalBurden()
    upper = (notes or "").upper()

    withdrawn = "WITHDRAWN" in upper
    if withdrawn:
        b.assessed = True
        b.score = max(b.score, 1.0)
        b.reasons.append("withdrawn from market")
    if "black-box" in (notes or ""):
        b.assessed = True
        b.score = max(b.score, 0.80)
        b.reasons.append("black-box warning")

    cls = classify(drug)
    if cls:
        b.assessed = True
        b.drug_class = cls.name
        b.score = max(b.score, cls.burden)
        b.reasons.append(f"{cls.name}: {cls.why}")

    if target:
        tier = regulon_mod.tier_for(target)
        if tier == "global":
            b.assessed = True
            b.score = max(b.score, 0.85)
            b.reasons.append(
                f"{target} is a global regulator — the receptor is systemic "
                f"even if the drug is selective for it")
        elif tier == "broad":
            b.assessed = True
            b.score = max(b.score, 0.60)
            b.reasons.append(
                f"{target} controls a broad programme; effects reach well "
                f"beyond your gene")

    b.level = _level_for(b.score, withdrawn)
    if not b.assessed:
        b.reasons.append(
            "not in the curated class table, no regulatory flag, and its "
            "target is untiered — UNASSESSED, which is not the same as clean")
    return b
