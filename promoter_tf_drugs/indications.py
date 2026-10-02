"""What each drug actually is, and what it is approved to treat.

Why the reports needed this
---------------------------
The tables name compounds and give a mechanism against one receptor, which
is the wrong unit for the decision being made. A reader looking at
tazarotene needs to know it is a topical acne and psoriasis retinoid before
the row means anything, and nothing in the report said so. Worse, the
compounds that dominate these tables are the ones with the largest ligand
series, which are oncology agents and hormones, so a list that looks
interchangeable is actually a mix of a topical gel and a cytotoxic
infusion.

Every route in this package is off-label by construction. The drug was
approved for one thing and is being considered for another, so the approved
indication is the baseline the proposal departs from, and the reader cannot
judge the departure without it.

Where the text comes from
-------------------------
Three sources, in order, and each row says which one it used.

  curated   a short hand-written entry for the compounds these reports
            actually return. Exact, and incomplete by construction.
  ChEMBL    `/drug_indication`, the approved indications ChEMBL records,
            taken at the highest phase reached. Needs the network, and is
            cached like every other request.
  class     the systemic class from `systemic.py`, which gives the kind of
            drug rather than its indication. A fallback, labelled as one.

None of these is a label. ChEMBL's indication list includes things studied
but never approved, the curated text is a summary written for orientation,
and the class line is a category. Check the actual prescribing information
before acting on any of it.
"""

from __future__ import annotations

import logging

from .congeners import salt_free_name

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


# --- What a systemic class is used for. The coarsest fallback, and the one
# --- that always works offline, because membership is already curated.
CLASS_USE = {
    "glucocorticoid":
        "inflammation and immune suppression - asthma, autoimmune disease, "
        "transplant rejection, some leukaemias",
    "mineralocorticoid":
        "adrenal insufficiency and salt-wasting states",
    "retinoid":
        "acne, psoriasis and other keratinisation disorders topically; "
        "acute promyelocytic leukaemia and cutaneous lymphoma systemically",
    "thyroid hormone":
        "hypothyroidism replacement",
    "sex steroid":
        "contraception, hormone replacement, and hormone-sensitive cancers",
    "vitamin D receptor agonist":
        "secondary hyperparathyroidism in renal disease, osteoporosis, and "
        "plaque psoriasis topically",
    "PPAR agonist":
        "type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia "
        "(fibrates)",
    "cytotoxic":
        "cancer chemotherapy",
    "broad immunosuppressant":
        "transplant rejection and autoimmune disease",
    "chromatin agent":
        "haematological malignancy, mostly myelodysplastic syndrome and "
        "T-cell lymphoma",
}


# --- Per-drug entries. Deliberately short: a phrase for what the molecule
# --- is, and a phrase for what it is approved to treat. Keyed salt-free
# --- and lowercased. This list exists because these are the compounds the
# --- reports actually return, not because it is a drug dictionary.
CURATED: dict = {
    # retinoids - the class that fills these tables
    "tretinoin": ("RAR agonist, all-trans retinoic acid",
                  "acute promyelocytic leukaemia; acne topically"),
    "isotretinoin": ("retinoid", "severe nodular acne"),
    "alitretinoin": ("RAR/RXR pan-agonist",
                     "Kaposi sarcoma skin lesions; chronic hand eczema"),
    "bexarotene": ("RXR-selective agonist (rexinoid)",
                   "cutaneous T-cell lymphoma"),
    "tamibarotene": ("RARalpha-selective agonist",
                     "relapsed acute promyelocytic leukaemia (Japan)"),
    "adapalene": ("RARbeta/gamma agonist",
                  "acne vulgaris, topical only"),
    "trifarotene": ("RARgamma-selective agonist",
                    "acne vulgaris, topical only"),
    "tazarotene": ("RARbeta/gamma agonist",
                   "plaque psoriasis and acne, topical only"),
    "acitretin": ("retinoid", "severe psoriasis"),
    "etretinate": ("retinoid", "severe psoriasis - WITHDRAWN, teratogenic"),
    # vitamin D
    "calcitriol": ("active vitamin D, VDR agonist",
                   "hypocalcaemia and secondary hyperparathyroidism in "
                   "renal failure"),
    "paricalcitol": ("selective VDR activator",
                     "secondary hyperparathyroidism in chronic kidney "
                     "disease"),
    "alfacalcidol": ("vitamin D prodrug", "renal osteodystrophy, rickets"),
    "calcifediol": ("25-hydroxyvitamin D", "vitamin D deficiency"),
    "doxercalciferol": ("vitamin D prodrug",
                        "secondary hyperparathyroidism in dialysis"),
    "tacalcitol": ("vitamin D analogue", "plaque psoriasis, topical"),
    "calcipotriol": ("vitamin D analogue", "plaque psoriasis, topical"),
    # corticosteroids
    "dexamethasone": ("potent glucocorticoid receptor agonist",
                      "inflammation, cerebral oedema, myeloma, and as an "
                      "antiemetic"),
    "prednisolone": ("glucocorticoid",
                     "asthma, autoimmune and inflammatory disease"),
    "prednisone": ("glucocorticoid prodrug",
                   "asthma, autoimmune and inflammatory disease"),
    "hydrocortisone": ("glucocorticoid and mineralocorticoid",
                       "adrenal insufficiency; topical dermatitis"),
    "betamethasone": ("potent glucocorticoid",
                      "inflammatory skin disease; fetal lung maturation"),
    "budesonide": ("glucocorticoid, high first-pass clearance",
                   "asthma, Crohn disease, ulcerative colitis"),
    "fluticasone": ("glucocorticoid, inhaled or intranasal",
                    "asthma and allergic rhinitis"),
    "clobetasol": ("very potent glucocorticoid",
                   "severe inflammatory skin disease, topical only"),
    "triamcinolone": ("glucocorticoid",
                      "inflammatory and allergic disease"),
    "mometasone": ("glucocorticoid", "asthma, rhinitis, dermatoses"),
    "fludrocortisone": ("mineralocorticoid",
                        "adrenal insufficiency, orthostatic hypotension"),
    # thyroid
    "levothyroxine": ("synthetic T4", "hypothyroidism"),
    "liothyronine": ("synthetic T3", "hypothyroidism, myxoedema coma"),
    "dextrothyroxine": ("T4 stereoisomer",
                        "formerly hyperlipidaemia - withdrawn for cardiac "
                        "mortality"),
    # PPAR
    "rosiglitazone": ("PPARgamma agonist (thiazolidinedione)",
                      "type 2 diabetes; restricted for cardiovascular risk"),
    "pioglitazone": ("PPARgamma agonist", "type 2 diabetes"),
    "fenofibrate": ("PPARalpha agonist", "hypertriglyceridaemia"),
    "gemfibrozil": ("PPARalpha agonist", "hypertriglyceridaemia"),
    # sex steroids and modulators
    "tamoxifen": ("selective oestrogen receptor modulator",
                  "ER-positive breast cancer"),
    "raloxifene": ("selective oestrogen receptor modulator",
                   "osteoporosis; breast cancer risk reduction"),
    "estradiol": ("oestrogen receptor agonist",
                  "menopausal hormone replacement"),
    "testosterone": ("androgen receptor agonist", "male hypogonadism"),
    "bicalutamide": ("androgen receptor antagonist", "prostate cancer"),
    "spironolactone": ("mineralocorticoid receptor antagonist",
                       "heart failure, ascites, hyperaldosteronism, acne"),
    # the kinase inhibitors the LINCS activity scan surfaces
    "trametinib": ("MEK1/2 inhibitor",
                   "BRAF-mutant melanoma and other BRAF-mutant tumours"),
    "selumetinib": ("MEK1/2 inhibitor",
                    "neurofibromatosis type 1 plexiform neurofibroma"),
    "gefitinib": ("EGFR tyrosine kinase inhibitor",
                  "EGFR-mutant non-small-cell lung cancer"),
    "erlotinib": ("EGFR tyrosine kinase inhibitor",
                  "non-small-cell lung cancer, pancreatic cancer"),
    # cytotoxics and others that recur
    "mitoxantrone": ("topoisomerase II inhibitor, DNA intercalator",
                     "acute myeloid leukaemia, prostate cancer, "
                     "multiple sclerosis"),
    "doxorubicin": ("anthracycline topoisomerase II inhibitor",
                    "many solid and haematological cancers"),
    "plicamycin": ("RNA synthesis inhibitor, binds GC-rich DNA",
                   "formerly hypercalcaemia of malignancy - withdrawn"),
    "vorinostat": ("histone deacetylase inhibitor",
                   "cutaneous T-cell lymphoma"),
    "romidepsin": ("histone deacetylase inhibitor",
                   "cutaneous and peripheral T-cell lymphoma"),
    "azacitidine": ("DNA methyltransferase inhibitor",
                    "myelodysplastic syndrome, acute myeloid leukaemia"),
    "decitabine": ("DNA methyltransferase inhibitor",
                   "myelodysplastic syndrome"),
    "maraviroc": ("CCR5 antagonist", "CCR5-tropic HIV-1 infection"),
    "celecoxib": ("COX-2 selective NSAID",
                  "osteoarthritis, rheumatoid arthritis, pain"),
    "metformin": ("biguanide, hepatic gluconeogenesis inhibitor",
                  "type 2 diabetes"),
    "atorvastatin": ("HMG-CoA reductase inhibitor",
                     "hypercholesterolaemia and cardiovascular prevention"),
    "simvastatin": ("HMG-CoA reductase inhibitor",
                    "hypercholesterolaemia"),
    "sirolimus": ("mTOR inhibitor",
                  "transplant rejection prophylaxis, lymphangioleiomyomatosis"),
    "thalidomide": ("cereblon modulator, immunomodulatory",
                    "multiple myeloma, erythema nodosum leprosum"),
    "lenalidomide": ("cereblon modulator", "multiple myeloma, MDS"),
    "hydroxyurea": ("ribonucleotide reductase inhibitor",
                    "sickle cell disease, chronic myeloid leukaemia"),
    "gentamicin": ("aminoglycoside antibiotic",
                   "serious Gram-negative infection"),
    "risdiplam": ("SMN2 splicing modifier", "spinal muscular atrophy"),
    "nusinersen": ("SMN2 splice-switching antisense oligonucleotide",
                   "spinal muscular atrophy"),
    "diazepam": ("GABA-A positive allosteric modulator",
                 "anxiety, seizures, muscle spasm"),
    "atenolol": ("beta-1 selective adrenergic antagonist",
                 "hypertension, angina"),
    "loperamide": ("peripheral mu-opioid agonist", "diarrhoea"),
    "caffeine": ("adenosine receptor antagonist",
                 "apnoea of prematurity; stimulant"),
    "cholecalciferol": ("vitamin D3, a precursor requiring activation",
                       "vitamin D deficiency and supplementation"),
    "pemafibrate": ("selective PPARalpha modulator",
                    "hypertriglyceridaemia (Japan)"),
    "bezafibrate": ("pan-PPAR agonist", "hypertriglyceridaemia"),
    "seocalcitol": ("vitamin D analogue",
                    "investigational, not approved"),
    "eldecalcitol": ("vitamin D analogue", "osteoporosis (Japan)"),
    "fluocinonide": ("potent glucocorticoid",
                     "inflammatory skin disease, topical only"),
}

# US and INN spellings of the same molecule. Kept separate from CURATED so
# the curated table stays one entry per drug: a report can name either
# spelling and `salt_free_name` will not bridge them, since they are
# different words rather than different salts.
SYNONYMS = {
    "calcipotriene": "calcipotriol",
    "tretinoin topical": "tretinoin",
    "all-trans retinoic acid": "tretinoin",
    "atra": "tretinoin",
    "13-cis-retinoic acid": "isotretinoin",
    "9-cis-retinoic acid": "alitretinoin",
    "dexamethasone acetate": "dexamethasone",
    "vitamin d3": "cholecalciferol",
    "1,25-dihydroxyvitamin d3": "calcitriol",
    "levothyroxine sodium": "levothyroxine",
    "mitoxantrone hydrochloride": "mitoxantrone",
}


def _from_class(drug_name: str) -> tuple | None:
    """(what it is, what it treats) from the systemic class table."""
    from . import systemic
    cls = systemic.classify(drug_name)
    if cls is None:
        return None
    use = CLASS_USE.get(cls.name)
    return (f"{cls.name} class", use) if use else None


def fetch_indications(client, chembl_id: str, *, limit: int = 6) -> list:
    """Approved indications for one molecule, highest phase first.

    ChEMBL's `/drug_indication` records what a compound was developed for,
    which includes indications studied and abandoned. `max_phase_for_ind`
    is what separates an approved use from a trial, so it is read and
    reported rather than treated as all one thing.
    """
    if not chembl_id or not str(chembl_id).startswith("CHEMBL"):
        return []
    try:
        res = client.get(f"{CHEMBL}/drug_indication", {
            "molecule_chembl_id": chembl_id, "format": "json",
            "limit": 50})
    except Exception as exc:                                # noqa: BLE001
        log.debug("indication lookup failed for %s: %s", chembl_id,
                  str(exc)[:90])
        return []

    seen: dict = {}
    for row in (res or {}).get("drug_indications", []):
        term = (row.get("mesh_heading") or row.get("efo_term") or "").strip()
        if not term:
            continue
        try:
            phase = float(row.get("max_phase_for_ind") or 0)
        except (TypeError, ValueError):
            phase = 0.0
        if term not in seen or phase > seen[term]:
            seen[term] = phase

    ranked = sorted(seen.items(), key=lambda kv: -kv[1])
    return [{"indication": t, "max_phase": p} for t, p in ranked[:limit]]


def route_of(drug_name: str) -> dict | None:
    """FDA route of administration, from the Orange Book if it is on disk.

    The curated entries below carry "topical only" for about seventy
    drugs, typed by hand. The Orange Book has the route for every approved
    ingredient from the agency that approved them, so it is consulted
    first and the hand list is the fallback.
    """
    try:
        from . import orangebook as ob
        return ob.lookup(drug_name)
    except Exception:                                       # noqa: BLE001
        return None


def describe(drug_name: str, *, mechanism: str | None = None,
             client=None, chembl_id: str | None = None) -> dict:
    """What this drug is and what it is approved to treat.

    Returns `what_it_is`, `treats` and `source`. `treats` is None when
    nothing could be established, which is not the same as a drug with no
    approved use - say so rather than leaving the cell looking empty.
    """
    stem = salt_free_name((drug_name or "").strip().lower())

    stem = SYNONYMS.get(stem, stem)

    fda = route_of(drug_name)

    def _with_route(d: dict) -> dict:
        if fda:
            d["route"] = fda["reach"]
            d["routes"] = fda["routes"]
            d["route_note"] = fda["reach_note"]
            d["route_source"] = "FDA Orange Book"
            d["generic_available"] = fda["generic_available"]
        else:
            # Absent is not systemic. A drug whose route is unknown has
            # not been shown to reach anything.
            d["route"] = None
            d["route_source"] = None
        return d

    hit = CURATED.get(stem)
    if hit:
        return _with_route({"what_it_is": hit[0], "treats": hit[1],
                            "source": "curated"})

    if client is not None and chembl_id:
        inds = fetch_indications(client, chembl_id)
        approved = [i["indication"] for i in inds if i["max_phase"] >= 4]
        shown = approved or [i["indication"] for i in inds]
        if shown:
            return _with_route({
                "what_it_is": mechanism or None,
                "treats": "; ".join(shown[:3]).lower()
                          + ("" if approved else
                             " (studied, not confirmed approved here)"),
                "source": "ChEMBL",
            })

    cls = _from_class(drug_name)
    if cls:
        return _with_route({"what_it_is": cls[0], "treats": cls[1],
                            "source": "class"})

    return _with_route({"what_it_is": mechanism or None, "treats": None,
                        "source": "unknown"})


def describe_all(drugs, *, client=None, chembl_ids: dict | None = None,
                 mechanisms: dict | None = None) -> dict:
    """`describe` over a set of names, deduplicated on the salt-free stem.

    Returns `{display name: description}`. Two salts of one drug share an
    entry, because repeating the same sentence for sodium and hydrochloride
    forms is the sort of thing that made these reports long.
    """
    out: dict = {}
    seen: set = set()
    for name in drugs:
        if not name:
            continue
        stem = salt_free_name(name.strip().lower())
        if stem in seen:
            continue
        seen.add(stem)
        out[name] = describe(
            name, client=client,
            chembl_id=(chembl_ids or {}).get(name),
            mechanism=(mechanisms or {}).get(name))
    return out


def summarise(described: dict) -> str:
    """One line on where the text came from, so it is not read as a label."""
    if not described:
        return ""
    by_source: dict = {}
    for d in described.values():
        by_source[d["source"]] = by_source.get(d["source"], 0) + 1
    bits = ", ".join(f"{n} {s}" for s, n in sorted(by_source.items()))
    return (f"Descriptions come from: {bits}. None of these is a label. The "
            f"curated text is a summary written for orientation, ChEMBL's "
            f"indication list includes uses studied but never approved, and "
            f"a class line gives the kind of drug rather than its "
            f"indication. Check the prescribing information before acting "
            f"on any row.")
