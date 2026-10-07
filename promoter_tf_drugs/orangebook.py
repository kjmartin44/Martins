"""The FDA's 20,000 approved products, used for what they actually add.

The number, and what it counts
------------------------------
The FDA says "there are over 20,000 prescription drug products approved
for marketing", and that is exact. The word doing the work is *products*.
A product is one strength of one dosage form from one applicant, so
metformin hydrochloride is seven products before you leave the first
manufacturer: three tablet strengths, an extended-release tablet, an oral
solution, and the same tablet again from each generic house. Those seven
are one molecule.

Collapse the Orange Book onto distinct active ingredients and about 20,000
becomes roughly 1,500 to 2,000. So the product list cannot enlarge a
screening library: docking metformin 500 mg and metformin 1 g returns the
same number twice, and running all 20,000 would be thirteen times the
compute for the same answers.

What it adds that no structure source has
-----------------------------------------
**Route of administration.** This is the reason to load it. Nothing in
ChEMBL, DrugCentral, the Repurposing Hub or ZINC says whether an approved
drug is a tablet or a gel, and the difference decides whether a candidate
can reach your target at all. A topical-only retinoid is not a route to a
gene in the liver however well it binds the receptor.

The package had this by hand for about seventy drugs, written into
`indications.CURATED`. The Orange Book has it for every approved
ingredient, from the agency that approved them.

Three more things come free:

  dosage form        extended release and depot formulations change the
                     exposure a dose achieves
  applicant count    several applicants for one ingredient means generics
                     exist, which is a practical fact about whether a
                     repurposing proposal is affordable
  approval date      the per-product FDA decision, which is what
                     "FDA-approved" actually refers to - ChEMBL's
                     max_phase cannot distinguish an FDA approval from
                     one granted elsewhere

What it does not have
---------------------
Structures. Not one SMILES anywhere in the Orange Book, so it can only be
joined to the structure sources on the ingredient name. That join is the
weak step: names carry salt forms, and "METFORMIN HYDROCHLORIDE" has to
reach "metformin". `congeners.salt_free_name` already does that work, and
a row that cannot be matched is reported rather than dropped quietly.

Nor does it cover biologics, which live in the Purple Book, or anything
approved outside the United States.
"""

from __future__ import annotations

import csv
import io
import logging
import zipfile
from pathlib import Path

from .congeners import salt_free_name

log = logging.getLogger(__name__)

DEFAULT_PATHS = ("data/products.txt", "data/Products.txt",
                 "data/orange_book_products.txt", "data/EOBZIP.zip")

SOURCE = """\
FDA Orange Book, "Approved Drug Products with Therapeutic Equivalence
Evaluations":

  https://www.fda.gov/drugs/drug-approvals-and-databases/approved-drug-products-therapeutic-equivalence-evaluations-orange-book

Download the current EOBZIP.zip and either unzip products.txt into data/
or drop the zip itself there - both are read. About 4 MB, no registration.
"""

# Routes that cannot deliver a drug to an internal target. A candidate
# available only this way is not a route to a gene anywhere but the
# tissue it is applied to.
LOCAL_ONLY_ROUTES = frozenset({
    "TOPICAL", "OPHTHALMIC", "OTIC", "NASAL", "DENTAL", "VAGINAL",
    "RECTAL", "INHALATION", "RESPIRATORY (INHALATION)", "BUCCAL",
    "SUBLINGUAL", "INTRAOCULAR", "INTRAVITREAL", "PERIODONTAL",
    "URETHRAL", "IRRIGATION", "INTRACAVITARY",
})
# Routes that reach the whole body.
SYSTEMIC_ROUTES = frozenset({
    "ORAL", "INTRAVENOUS", "INTRAMUSCULAR", "SUBCUTANEOUS",
    "TRANSDERMAL", "INTRAPERITONEAL", "INTRATHECAL", "PARENTERAL",
    "INJECTION", "INTRA-ARTERIAL", "EPIDURAL", "IMPLANTATION",
})

_PRODUCTS: dict = {}
_LOADED: str | None = None

# Whether the absence of the file has already been reported. `load` is
# called once per drug through `lookup`, and without this the same
# notice is logged hundreds of times in one run.
_ABSENT_REPORTED = False


def resolve(explicit=None) -> str | None:
    for c in ([explicit] if explicit else list(DEFAULT_PATHS)):
        if c and Path(c).exists():
            return str(c)
    return None


def _rows(path: Path):
    """Orange Book rows, from products.txt or straight out of EOBZIP.zip.

    The file is tilde-delimited rather than tab, which is unusual enough
    that a tab-based reader returns one enormous column and no error.
    """
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            name = next((n for n in zf.namelist()
                         if n.lower().endswith("products.txt")), None)
            if name is None:
                raise ValueError(
                    f"{path} has no products.txt in it. {SOURCE}")
            with zf.open(name) as fh:
                text = io.TextIOWrapper(fh, encoding="utf-8",
                                        errors="replace")
                yield from csv.DictReader(text, delimiter="~")
                return
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        yield from csv.DictReader(fh, delimiter="~")


def reset() -> None:
    """Forget what was loaded, and that the file was found missing.

    For tests, and for a process that drops a file into place and wants
    the next lookup to find it.
    """
    global _PRODUCTS, _LOADED, _ABSENT_REPORTED
    _PRODUCTS = {}
    _LOADED = None
    _ABSENT_REPORTED = False


def load(path=None) -> int:
    """Collapse the product table onto distinct active ingredients.

    Idempotent. Returns the number of ingredients, not products, because
    the collapse is the point.
    """
    global _PRODUCTS, _LOADED, _ABSENT_REPORTED
    if _PRODUCTS:
        return len(_PRODUCTS)
    # The lookup itself always runs: a file that appears after an earlier
    # miss - dropped into data/ mid-session, or written by a test - must
    # still be found. It is a handful of is_file() calls, so repeating it
    # costs nothing worth saving.
    #
    # Only the NOTICE is said once. `load` runs once per drug through
    # `lookup`, and a real run printed this line hundreds of times, which
    # buries the rest of the log and reads like a recurring fault rather
    # than a file that is simply absent.
    p = resolve(path)
    if not p:
        if not _ABSENT_REPORTED:
            _ABSENT_REPORTED = True
            log.info("no Orange Book on disk, so route of administration "
                     "is unknown for every row - which is not the same as "
                     "every row being systemic. To add it: %s",
                     " ".join(SOURCE.split()))
        return 0

    n_products = 0
    for row in _rows(Path(p)):
        ingredient = (row.get("Ingredient") or "").strip()
        if not ingredient:
            continue
        n_products += 1
        # A combination product lists its ingredients separated by ";".
        # Each is recorded, and the combination is noted, because a
        # combination is not evidence about the single agent.
        parts = [i.strip() for i in ingredient.split(";") if i.strip()]
        df_route = (row.get("DF;Route") or "").strip()
        form, _, route = df_route.partition(";")
        routes = [r.strip().upper() for r in route.split(",") if r.strip()]

        for part in parts:
            key = salt_free_name(part.lower())
            if not key:
                continue
            rec = _PRODUCTS.setdefault(key, {
                "ingredient": part.title(), "routes": set(), "forms": set(),
                "applicants": set(), "trade_names": set(),
                "n_products": 0, "approvals": set(),
                "is_combination_only": True, "marketing": set(),
            })
            rec["n_products"] += 1
            rec["routes"].update(routes)
            if form:
                rec["forms"].add(form.strip().upper())
            app = (row.get("Applicant") or "").strip()
            if app:
                rec["applicants"].add(app)
            tn = (row.get("Trade_Name") or "").strip()
            if tn:
                rec["trade_names"].add(tn)
            ad = (row.get("Approval_Date") or "").strip()
            if ad:
                rec["approvals"].add(ad)
            mk = (row.get("Type") or "").strip().upper()
            if mk:
                rec["marketing"].add(mk)
            if len(parts) == 1:
                rec["is_combination_only"] = False

    _LOADED = p
    log.info("Orange Book: %d product(s) collapse to %d distinct active "
             "ingredient(s) from %s", n_products, len(_PRODUCTS),
             Path(p).name)
    return len(_PRODUCTS)


def product_count() -> int:
    load()
    return sum(r["n_products"] for r in _PRODUCTS.values())


def lookup(name: str, path=None) -> dict | None:
    """FDA product facts for one drug name, or None if it is not listed."""
    load(path)
    if not _PRODUCTS:
        return None
    key = salt_free_name((name or "").strip().lower())
    rec = _PRODUCTS.get(key)
    if rec is None:
        return None

    routes = set(rec["routes"])
    systemic = sorted(routes & SYSTEMIC_ROUTES)
    local = sorted(routes & LOCAL_ONLY_ROUTES)
    unclassified = sorted(routes - SYSTEMIC_ROUTES - LOCAL_ONLY_ROUTES)

    out = {
        "ingredient": rec["ingredient"],
        "routes": sorted(routes),
        "systemic_routes": systemic,
        "local_routes": local,
        "unclassified_routes": unclassified,
        "forms": sorted(rec["forms"])[:8],
        "n_products": rec["n_products"],
        "n_applicants": len(rec["applicants"]),
        "generic_available": len(rec["applicants"]) > 1,
        "trade_names": sorted(rec["trade_names"])[:4],
        "combination_only": rec["is_combination_only"],
        "fda_approved": True,
    }

    if systemic:
        out["reach"] = "systemic"
        out["reach_note"] = (
            f"approved systemically ({', '.join(systemic).lower()}), so it "
            f"can in principle reach an internal target")
    elif local:
        out["reach"] = "local only"
        out["reach_note"] = (
            f"approved ONLY by local routes ({', '.join(local).lower()}). "
            f"It cannot reach a target outside the tissue it is applied "
            f"to, however well it binds - so for an internal gene this is "
            f"not a usable route without reformulation")
    else:
        out["reach"] = "unclear"
        out["reach_note"] = (
            f"route(s) {', '.join(unclassified) or 'not recorded'} are not "
            f"classified here, so whether it reaches an internal target is "
            f"unresolved rather than assumed")

    if rec["is_combination_only"]:
        out["notes"] = [
            "only ever approved as part of a combination product, so the "
            "approval is not evidence about this agent on its own"]
    return out


def enrich(library: list, path=None) -> dict:
    """Attach FDA product facts to a merged library, in place.

    Returns an accounting: how many molecules matched, how many are
    local-route-only, and which names could not be joined. The join is on
    the salt-free ingredient name, because the Orange Book has no
    structures at all - so unmatched rows are reported rather than
    silently left bare.
    """
    load(path)
    if not _PRODUCTS:
        return {"available": False, "matched": 0, "unmatched": [],
                "note": ("no Orange Book on disk, so route of "
                         "administration is unknown for every row. That is "
                         "different from every row being systemic.")}

    matched, unmatched, local_only = 0, [], []
    for row in library:
        fda = lookup(row.get("name") or "")
        if fda is None:
            unmatched.append(row.get("name"))
            continue
        matched += 1
        row["fda"] = fda
        row["route_reach"] = fda["reach"]
        row["generic_available"] = fda["generic_available"]
        if fda["reach"] == "local only":
            local_only.append(row.get("name"))

    return {
        "available": True, "source": _LOADED,
        "matched": matched, "n_library": len(library),
        "unmatched": unmatched[:40], "n_unmatched": len(unmatched),
        "local_only": local_only,
        "n_products_in_book": product_count(),
        "n_ingredients_in_book": len(_PRODUCTS),
    }


def summarise(acc: dict) -> str:
    if not acc.get("available"):
        return acc.get("note", "") + " " + SOURCE.splitlines()[0]

    bits = [
        f"The Orange Book lists {acc['n_products_in_book']:,} approved "
        f"products, which collapse to {acc['n_ingredients_in_book']:,} "
        f"distinct active ingredients. The FDA's 'over 20,000' counts "
        f"products, meaning each strength, dosage form and applicant "
        f"separately. So it cannot enlarge a screening library; what it "
        f"adds is what no structure source has.",
        f"{acc['matched']} of {acc['n_library']} library molecules were "
        f"matched by ingredient name and now carry route of "
        f"administration, dosage forms and generic availability.",
    ]
    if acc.get("local_only"):
        names = ", ".join(acc["local_only"][:6])
        bits.append(
            f"**{len(acc['local_only'])} "
            f"{'is' if len(acc['local_only']) == 1 else 'are'} "
            f"approved only by local routes** "
            f"({names}). Those cannot reach a target outside the tissue "
            f"they are applied to, which rules them out for an internal "
            f"gene no matter how they score here.")
    if acc.get("n_unmatched"):
        bits.append(
            f"{acc['n_unmatched']} could not be matched by name. The Orange "
            f"Book holds no structures, so the join is on the ingredient "
            f"name and a mismatch means unknown route rather than no "
            f"approval.")
    return " ".join(bits)
