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
            "doxorubicin", "daunorubicin", "epirubicin", "idarubicin",
            "paclitaxel", "docetaxel", "vincristine", "vinblastine",
            "etoposide", "cisplatin", "carboplatin", "oxaliplatin",
            "methotrexate", "cyclophosphamide", "ifosfamide", "fluorouracil",
            "gemcitabine", "cytarabine", "dactinomycin", "bleomycin",
            "mitomycin", "plicamycin", "mithramycin",
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

BURDEN_LEVELS = ("low", "moderate", "high", "prohibitive")


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
    # Multi-word brand-ish names: fall back to the head token, which is the
    # active moiety ("Fluticasone Propionate HFA" -> "fluticasone").
    head = name.split()[0] if name.split() else ""
    return _INDEX.get(head)


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
