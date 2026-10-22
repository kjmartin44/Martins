"""Can the drug reach the concentration at which the evidence was seen?

Why this exists
---------------
A signature in LINCS says a compound moved a gene at some concentration in a
dish, usually 10 uM. A binding assay says a compound engages a target at some
potency. Neither says a person taking the approved dose reaches that
concentration. A drug whose blood level is a hundred times below the one that
produced the effect is not a repurposing candidate, whatever its rank.

What this module does
---------------------
For each candidate drug in a saved `results.json` it sets two numbers against
each other:

  needed    the concentration the evidence was seen at: the lowest LINCS dose
            at which the drug moved this gene, else the drug's potency at the
            target the route goes through (from its pChEMBL)
  reached   a published human peak blood level (Cmax) for the approved use

and reports the ratio and one of four labels. It writes `exposure.json` and
`exposure.md`. It changes no rank and no score.

Where the exposure comes from
-----------------------------
In order of trust:

1. A local table you supply (`--cmax-table`). Values are taken as given and
   tagged with the source you name. This is the intended way to use a curated
   source such as a label compilation, because no open, machine-readable
   source of human Cmax values is reliable enough to read by program.
2. Human plasma records in ChEMBL, filtered strictly: the assay description
   must say human, and must not mention cerebrospinal fluid, tissue, or an
   animal. ChEMBL's coverage of approved drugs is thin, so most drugs will
   come back "not assessable". That is the honest answer, not a gap to fill
   with a guess.

`data/fda_label_cmax.tsv`, when it exists, is read as a table of this kind and
sits below one you name. It is built by `scripts/build_fda_cmax.py` from the
pharmacokinetics sections of FDA labels (see `fdalabel.py`). The first attempt
at label text was dropped because a Cmax can belong to an interaction study or
another dose and a wrong number read silently is worse than none. The reader
now refuses special populations, interactions, food and metabolites, and sets
aside any drug whose figures disagree, but it is a text reader: each row carries
the label, the dose and the sentence so a person can check it.

What it does not do
-------------------
**Total is not free.** Most values are total drug in plasma. Only the free
fraction acts, and for a highly protein bound drug it can be a few percent of
the total. So "far below" is firm (total already under the target, free is
lower still) and "reaches" is the weakest claim here. If your table gives the
unbound fraction, the free level is used and says so.

**A dish is not a body.** LINCS doses are nominal concentrations in medium. A
tissue the drug reaches poorly can sit far below the plasma level, and one it
concentrates in can sit above. This looks only at plasma.

**One dose is not the regimen.** A single-dose Cmax understates steady state
for a drug that accumulates and overstates it for one given at a lower dose.
The dose is kept when the source states it.

**Potency is not the gene's threshold.** Using a drug's potency at a
transcription factor as the needed concentration assumes the gene responds
when the factor does, which is the unmeasured step in a mechanism route. It is
used only when no LINCS dose exists for this gene.

**Not assessable is not "no".** A drug with no exposure data, or no known
effective concentration, is reported as such and never as unreachable.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import statistics
from pathlib import Path

from .signatures import breadth_key

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Read automatically when present, below any table named with --cmax-table.
DEFAULT_CMAX_PATHS = ("data/fda_label_cmax.tsv",)

# Ratio of reached to needed. Total drug in plasma is compared with the
# concentration the evidence was seen at.
# Records for one drug that differ by more than this are not one measurement.
# The median of values thirty-fold apart describes none of them, so the drug
# is left not assessable and the disagreement is reported.
MAX_SPREAD = 30.0

REACHES = 1.0            # at or above
WITHIN_TENFOLD = 0.1     # below, but within a factor of ten

VERDICTS = ("reaches", "within tenfold", "far below", "not assessable")

_HUMAN = ("human", "patient", "volunteer", "subject", "healthy", " men ",
          " women ", "adult", "child", "pediatric", "paediatric")
_NOT_PLASMA = ("csf", "cerebrospinal", "tissue", "urine", "bile", "milk",
               "saliva", "spinal", "intratumor")
# An organ counts as the matrix when the record is "in the brain" or "liver
# concentration". "Patients with solid tumor" and "lung cancer" name who was
# dosed, not what was measured; reading them as a tissue threw away bexarotene's
# only human record.
_ORGANS = r"(?:brain|tumou?r|liver|lung|skin|muscle|kidney)"
_ORGAN_AS_MATRIX = re.compile(
    r"\b(?:in|from|of)\s+(?:(?:the|human|healthy|normal)\s+)*" + _ORGANS
    + r"\b(?!\s+(?:cancer|carcinoma|"
    r"patient|disease|transplant|impairment|failure|function|metasta|malignan|"
    r"cirrhosis|resection|tumou?r))"
    r"|\b" + _ORGANS + r"\s+(?:tissue|homogenate|level|concentration|exposure|"
    r"penetration|uptake)")
# Studies of a drug given with another, or to people whose handling of it is
# changed, report levels that are not what the approved dose gives. A calcitriol
# value in ChEMBL that comes from a ketoconazole interaction study read as a
# blood level thousands of times above the real one.
_NOT_THE_APPROVED_USE = (
    "co-treat", "cotreat", "co-administ", "coadminist", "co-dos", "coadmin",
    "concomitant", "in combination", "combination with", "interaction",
    "ketoconazole", "ritonavir", "itraconazole", "rifampin", "rifampicin",
    "pretreat", "pre-treat", "impair", "overdose", "inhibitor", "inducer")
_ANIMAL = ("mouse", "mice", "rat ", "rats", "dog", "monkey", "rabbit", "pig",
           "cynomolgus", "beagle", "hamster", "guinea", "marmoset", "sheep",
           "macaque", "murine", "wistar", "sprague")

_MASS_PER_VOLUME = {          # factor to nM is factor / molecular weight
    "ng/ml": 1e3, "ug/l": 1e3,
    "ug/ml": 1e6, "mg/l": 1e6,
    "mg/ml": 1e9, "g/l": 1e9,
    "pg/ml": 1.0,
}
_MOLAR = {"nm": 1.0, "um": 1e3, "mm": 1e6, "pm": 1e-3, "m": 1e9}


# ------------------------------------------------------------------ units
def _unit_key(unit: str) -> str:
    u = (unit or "").strip().lower().replace("µ", "u").replace("μ", "u")
    u = u.replace("mcg", "ug").replace(".", "/").replace("-1", "")
    return u.replace(" ", "").replace("^", "")


def to_nm(value, unit: str, mw: float | None = None) -> float | None:
    """A concentration in nM, or None when it cannot be converted.

    Mass units need the molecular weight. Without it the value is left out,
    never guessed.
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if v <= 0:
        return None
    key = _unit_key(unit)
    if key in _MOLAR:
        return v * _MOLAR[key]
    if key in _MASS_PER_VOLUME:
        if not mw or mw <= 0:
            return None
        return v * _MASS_PER_VOLUME[key] / mw
    return None


def molecular_weight(smiles: str) -> float | None:
    """Molecular weight of the drug without its counter-ions, or None."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors
        from .docking import split_counter_ions
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(split_counter_ions(smiles)["smiles"])
        return float(Descriptors.MolWt(mol)) if mol is not None else None
    except Exception as exc:                               # noqa: BLE001
        # No weight means mass units cannot be converted, which leaves the
        # drug not assessable. Say why instead of leaving a silent gap.
        log.debug("no molecular weight for %s: %s", smiles[:40], exc)
        return None


# -------------------------------------------------------------- the sources
def _default_tables() -> dict:
    """Values from the tables read without being named (`DEFAULT_CMAX_PATHS`)."""
    table: dict = {}
    for path in DEFAULT_CMAX_PATHS:
        if Path(path).exists():
            table.update(read_cmax_table(path))
            log.info("human Cmax values from %s", path)
    return table


def read_cmax_table(path) -> dict:
    """Human exposures from a table: `{name key: record}`.

    Columns: `name`, then either `cmax_nM`, or `cmax` and `unit` (and `mw`
    for a mass unit). Optional: `fu` (unbound fraction, 0 to 1), `dose`,
    `source`, `note`. A row that cannot be read is skipped with a warning and
    never read as zero. The same drug twice keeps the larger value, since one
    row per dose is the usual shape.
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig")
    delim = "\t" if text.count("\t") >= text.count(",") else ","
    out: dict = {}
    for i, row in enumerate(csv.DictReader(text.splitlines(), delimiter=delim),
                            start=2):
        # A row with more fields than the header puts the extras under a None
        # key, as a list; they belong to no column and are dropped.
        row = {k.strip().lower(): (v or "").strip()
               for k, v in row.items() if isinstance(k, str)
               and not isinstance(v, list)}
        name = row.get("name") or row.get("drug")
        if not name:
            continue
        if row.get("cmax_nm"):
            nm = to_nm(row["cmax_nm"], "nM")
        else:
            nm = to_nm(row.get("cmax"), row.get("unit", ""),
                       _float(row.get("mw")))
        if nm is None:
            log.warning("%s line %d: no usable Cmax for %s", path.name, i, name)
            continue
        fu = _float(row.get("fu"))
        if fu is not None and not 0 < fu <= 1:
            log.warning("%s line %d: unbound fraction %s for %s is not in "
                        "(0, 1]; ignored", path.name, i, row.get("fu"), name)
            fu = None
        key = breadth_key(name)
        rec = {"cmax_nm": nm, "fu": fu, "dose": row.get("dose") or None,
               "source": row.get("source") or f"table {path.name}",
               "note": row.get("note") or None, "n_records": 1,
               "range_nm": [nm, nm]}
        if key not in out or nm > out[key]["cmax_nm"]:
            out[key] = rec
    return out


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def chembl_human_cmax(client, chembl_id: str, mw: float | None,
                      diag: dict | None = None) -> dict | None:
    """Human plasma Cmax records for one drug from ChEMBL, or None.

    Keeps a record only when its description says human and names no tissue,
    fluid other than plasma, or animal; converts it to nM; and takes the
    median across records, since they come from different studies and doses.
    Returns `{"cmax_nm", "range_nm", "n_records", "n_skipped", "dose",
    "source", "fu": None}`. None means ChEMBL has no usable human record,
    which is not evidence about the drug. `diag`, when given, is filled with
    `n_total` (Cmax records ChEMBL holds for the drug) and `n_human` (those
    that are human, plasma, and not an interaction study), so a caller can say
    whether there was nothing or only animal data.
    """
    if not chembl_id:
        return None
    parent = chembl_id
    n_human = 0
    try:
        mol = client.get(f"{CHEMBL}/molecule/{chembl_id}.json")
        if mw is None:
            try:
                mw = float(((mol or {}).get("molecule_properties") or {})
                           .get("full_mwt"))
            except (TypeError, ValueError):
                mw = None
        parent = ((mol or {}).get("molecule_hierarchy") or {}).get(
            "parent_chembl_id") or chembl_id
        res = client.get(f"{CHEMBL}/activity.json", {
            "standard_type": "Cmax", "parent_molecule_chembl_id": parent,
            "limit": 200})
    except Exception as exc:                                # noqa: BLE001
        raise RuntimeError(f"ChEMBL lookup failed for {chembl_id}: "
                           f"{str(exc)[:80]}") from exc
    kept, skipped, doses = [], 0, []
    if diag is not None:
        diag["n_total"] = len((res or {}).get("activities", []))
    for a in (res or {}).get("activities", []):
        desc = f" {(a.get('assay_description') or '').lower()} "
        if not any(t in desc for t in _HUMAN):
            continue
        if (any(t in desc for t in _NOT_PLASMA + _ANIMAL
                + _NOT_THE_APPROVED_USE)
                or _ORGAN_AS_MATRIX.search(desc)):
            continue
        n_human += 1
        nm = to_nm(a.get("standard_value"), a.get("standard_units") or "", mw)
        if nm is None:
            skipped += 1
            continue
        kept.append(nm)
        m = re.search(r"\bat\s+([\d.]+\s*(?:mg|ug|mcg|g)(?:/kg)?)", desc)
        if m:
            doses.append(m.group(1))
    if diag is not None:
        diag["n_human"] = n_human
    if not kept:
        return {"cmax_nm": None, "n_records": 0, "n_skipped": skipped} \
            if skipped else None
    if max(kept) / min(kept) > MAX_SPREAD:
        return {"cmax_nm": None, "range_nm": [min(kept), max(kept)],
                "n_records": len(kept), "n_skipped": skipped,
                "disagree": True}
    return {"cmax_nm": statistics.median(kept),
            "range_nm": [min(kept), max(kept)], "n_records": len(kept),
            "n_skipped": skipped, "fu": None,
            "dose": ", ".join(sorted(set(doses))) or None,
            "source": f"ChEMBL human Cmax records for {parent}",
            "note": None}


# ------------------------------------------------------------------- needed
def needed_concentrations(block: dict) -> dict:
    """The concentration each drug's evidence was seen at: `{key: record}`.

    The lowest LINCS dose at which the drug moved this gene wins, because it
    is an observed effect on the gene itself. Failing that, the drug's
    potency at the factor its route goes through (the best pChEMBL median
    among its grades).
    """
    out: dict = {}
    for hit in ((block.get("lincs") or {}).get("hits") or []):
        dose = _float(hit.get("dose"))
        if not dose or dose <= 0:
            continue
        nm = dose * 1e3                    # LINCS doses are micromolar
        key = breadth_key(hit.get("drug") or "")
        if key and (key not in out or nm < out[key]["needed_nm"]):
            out[key] = {"needed_nm": nm, "needed_basis": (
                f"lowest LINCS dose that moved this gene: {dose:g} uM"
                f"{', ' + hit['cell_line'] if hit.get('cell_line') else ''}")}
    for g in block.get("drug_grades") or []:
        key = breadth_key(g.get("drug_name") or "")
        p = ((g.get("evidence") or {}).get("pchembl_median"))
        if key and key not in out and isinstance(p, (int, float)) and p > 0:
            out[key] = {"needed_nm": 10 ** (9 - p), "needed_basis": (
                f"potency at {g.get('tf_symbol')}: pChEMBL {p:g} "
                f"(the factor, not the gene)")}
    # A receptor screen's established ligands carry a measured potency too,
    # and are not graded drugs, so without this every one of them had no
    # concentration to compare against. The most potent target wins.
    best: dict = {}
    for tf, scr in (block.get("screens") or {}).items():
        for h in (scr or {}).get("hits") or []:
            p = h.get("measured_pchembl")
            key = breadth_key(h.get("drug") or "")
            if (key and key not in out and isinstance(p, (int, float))
                    and p > 0 and (key not in best or p > best[key][0])):
                best[key] = (p, tf)
    for key, (p, tf) in best.items():
        out[key] = {"needed_nm": 10 ** (9 - p), "needed_basis": (
            f"binding potency at {tf}: pChEMBL {p:g} "
            f"(the factor, not the gene)")}
    return out


def chembl_ids_from_block(block: dict) -> dict:
    """`{key: ChEMBL id}` for drugs the run already identified, from the
    receptor screens, so a drug missing from the approved table is not looked
    up by name when its id is in hand."""
    out: dict = {}
    for scr in (block.get("screens") or {}).values():
        for h in (scr or {}).get("hits") or []:
            cid, key = h.get("chembl_id"), breadth_key(h.get("drug") or "")
            if key and cid and str(cid).startswith("CHEMBL"):
                out.setdefault(key, str(cid))
    return out


def candidate_drugs(block: dict, top: int, fda_only: bool = False,
                    hide_local_only: bool = False) -> list:
    """Distinct drugs in the block's routes, best first: `[(key, name)]`.

    Approved drugs come before the rest. A compound that is not approved has
    no human dose to compare with, and when LINCS research compounds fill the
    list first they take the places the check exists for.
    """
    routes = (block.get("all_routes") or {}).get("routes") or []
    if fda_only:
        from . import fda
        routes = [r for r in routes
                  if fda.shown(r.get("drug") or "", hide_local_only)]
    ordered = sorted(routes, key=lambda r: (0 if r.get("approved") else 1,
                                            r.get("tier", 9),
                                            -(r.get("score") or 0)))
    seen, out = set(), []
    for r in ordered:
        key = breadth_key(r.get("drug") or "")
        if key and key not in seen:
            seen.add(key)
            out.append((key, r["drug"]))
        if len(out) >= top:
            break
    return out


# ------------------------------------------------------------------ verdict
def verdict(reached_nm: float | None, needed_nm: float | None) -> tuple:
    """`(label, ratio)` for a reached and a needed concentration in nM."""
    if not reached_nm or not needed_nm:
        return "not assessable", None
    ratio = reached_nm / needed_nm
    if ratio >= REACHES:
        return "reaches", ratio
    if ratio >= WITHIN_TENFOLD:
        return "within tenfold", ratio
    return "far below", ratio


def assess(block: dict, *, client=None, table: dict | None = None,
           library: list | None = None, top: int = 60,
           fda_only: bool = False, hide_local_only: bool = False) -> dict:
    """Exposure against need for the top candidate drugs of one result block.

    Returns `{"drugs": {name: record}, "counts": {...}, "errors": [...]}`.
    Each record has `verdict`, `ratio`, `needed_nm`, `reached_nm`, and the
    wording behind both, or `reason` when it is not assessable. `client` is
    only used when ChEMBL is consulted; without it only the table is read.
    """
    table = table or {}
    by_key = {breadth_key(r["name"]): r for r in (library or [])}
    needed = needed_concentrations(block)
    known_ids = chembl_ids_from_block(block)
    drugs, errors = {}, []
    for key, name in candidate_drugs(block, top, fda_only, hide_local_only):
        rec: dict = {"drug": name}
        need = needed.get(key)
        if need:
            rec.update(need)
        exp = table.get(key)
        diag: dict = {}
        if exp is None and client is not None:
            row = by_key.get(key)
            cid = (row or {}).get("chembl_id") or known_ids.get(key)
            if not cid:
                from . import chemblid
                try:
                    cid, _ = chemblid.resolve(client, name, None)
                except Exception as exc:                    # noqa: BLE001
                    log.warning("ChEMBL id lookup failed for %s: %s", name,
                                str(exc)[:80])
                    cid = None
            rec["chembl_id"] = cid
            try:
                exp = chembl_human_cmax(
                    client, cid,
                    molecular_weight(row["smiles"]) if row else None, diag)
            except RuntimeError as exc:
                errors.append(str(exc))
                exp = None
                rec["error"] = str(exc)
        if exp and exp.get("disagree"):
            lo, hi = exp["range_nm"]
            rec["disagree"] = (f"{exp['n_records']} human records range from "
                               f"{_fmt_nm(lo)} to {_fmt_nm(hi)}, more than "
                               f"{MAX_SPREAD:g}-fold, so none is used")
        if exp and exp.get("cmax_nm"):
            fu = exp.get("fu")
            reached = exp["cmax_nm"] * (fu if fu else 1.0)
            rec.update({"cmax_nm": exp["cmax_nm"], "reached_nm": reached,
                        "level": "free" if fu else "total",
                        "range_nm": exp.get("range_nm"),
                        "n_records": exp.get("n_records"),
                        "dose": exp.get("dose"), "source": exp.get("source")})
        label, ratio = verdict(rec.get("reached_nm"), rec.get("needed_nm"))
        rec["verdict"] = label
        rec["ratio"] = None if ratio is None else round(ratio, 4)
        if label == "not assessable":
            missing = []
            if not rec.get("reached_nm"):
                kind, text = _exposure_gap(rec, diag, client is not None)
                rec["exposure_gap"] = kind
                missing.append(text)
            if not rec.get("needed_nm"):
                rec["need_gap"] = True
                missing.append("no effective concentration is known for this "
                               "gene")
            rec["reason"] = "; ".join(missing)
        drugs[name] = rec
    counts = {v: sum(1 for d in drugs.values() if d["verdict"] == v)
              for v in VERDICTS}
    return {"drugs": drugs, "counts": counts, "errors": errors,
            "n_assessed": len(drugs), "gaps": _gap_counts(drugs)}


GAP_TEXT = {
    "disagree": "records disagree by more than 30-fold",
    "error": "the ChEMBL lookup failed",
    "no_id": "no ChEMBL record under this name (a research compound, or a "
             "name that does not match)",
    "no_records": "ChEMBL holds no Cmax record for the drug",
    "no_human_records": "ChEMBL has Cmax records, none from people in plasma",
    "not_looked_up": "no exposure table row and ChEMBL was not consulted",
}


def _exposure_gap(rec: dict, diag: dict, consulted: bool) -> tuple:
    """`(kind, sentence)` for why a drug has no human exposure value."""
    if rec.get("disagree"):
        return "disagree", rec["disagree"]
    if rec.get("error"):
        return "error", rec["error"]
    if not consulted:
        return "not_looked_up", ("no human exposure value found: "
                                 + GAP_TEXT["not_looked_up"])
    if not rec.get("chembl_id"):
        return "no_id", ("no human exposure value: " + GAP_TEXT["no_id"])
    if not diag.get("n_total"):
        return "no_records", ("no human exposure value: "
                              + GAP_TEXT["no_records"])
    return "no_human_records", (
        f"no human exposure value: ChEMBL has {diag['n_total']} Cmax "
        f"record(s) for the drug, none from people in plasma (animal studies, "
        f"another tissue, or an interaction study)")


def _gap_counts(drugs: dict) -> dict:
    """How many not assessable drugs lack each kind of exposure data, and how
    many lack a concentration to compare with."""
    kinds: dict = {}
    for d in drugs.values():
        if d.get("verdict") == "not assessable" and d.get("exposure_gap"):
            kinds[d["exposure_gap"]] = kinds.get(d["exposure_gap"], 0) + 1
    return {"exposure": kinds,
            "need": sum(1 for d in drugs.values() if d.get("need_gap"))}


# ------------------------------------------------------------------- output
def _fmt_nm(v) -> str:
    if v is None:
        return "-"
    if v >= 1e3:
        return f"{v / 1e3:.3g} uM"
    return f"{v:.3g} nM"


def markdown(result: dict, gene: str) -> str:
    c = result["counts"]
    L = [f"# Exposure against need: {gene}", "",
         f"{result['n_assessed']} candidate drugs. Reaches: {c['reaches']}. "
         f"Within tenfold: {c['within tenfold']}. Far below: "
         f"{c['far below']}. Not assessable: {c['not assessable']}.", "",
         "This compares published human plasma levels with the concentration "
         "the evidence was seen at. It changes no rank. Total plasma drug is "
         "not free drug, so \"far below\" is firm and \"reaches\" is the "
         "weakest claim here. A drug with no data is not assessable, which "
         "is not the same as unreachable.", ""]
    gaps = result.get("gaps") or {}
    if c["not assessable"] and gaps:
        parts = [f"{GAP_TEXT[k]} ({n})"
                 for k, n in sorted((gaps.get("exposure") or {}).items(),
                                    key=lambda kv: -kv[1])]
        L += ["Why the not assessable drugs could not be compared: "
              + "; ".join(parts) + f". {gaps.get('need', 0)} drug(s) have no "
              "concentration to compare with, which is a separate gap from "
              "the exposure one.", ""]
    if result.get("errors"):
        L += [f"{len(result['errors'])} lookup(s) failed and are marked in "
              f"the table; those drugs were not assessed.", ""]
    L += ["| Drug | Verdict | Reached | Needed | Ratio | Basis |",
          "|---|---|---|---|---|---|"]
    order = {v: i for i, v in enumerate(VERDICTS)}
    for name, d in sorted(result["drugs"].items(),
                          key=lambda kv: order[kv[1]["verdict"]]):
        basis = d.get("needed_basis") or d.get("reason") or ""
        if d.get("reached_nm") and d.get("source"):
            basis += f" | exposure: {d['level']}, {d['source']}"
            if d.get("dose"):
                basis += f", {d['dose']}"
        L.append(f"| {name.replace('|', '/')} | {d['verdict']} | "
                 f"{_fmt_nm(d.get('reached_nm'))} | "
                 f"{_fmt_nm(d.get('needed_nm'))} | "
                 f"{'-' if d.get('ratio') is None else d['ratio']} | "
                 f"{basis.replace('|', '/')} |")
    L.append("")
    return "\n".join(L)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m generoutes.exposure",
        description=("Set each candidate drug's published human blood level "
                     "against the concentration its evidence was seen at."))
    p.add_argument("results", help="a saved results.json")
    p.add_argument("--cmax-table", default=None,
                   help="table of human Cmax values (name, cmax_nM or cmax "
                        "and unit, optional mw, fu, dose, source)")
    p.add_argument("--library", default=None,
                   help="drug table for ChEMBL ids and structures; the "
                        "bundled approved library by default")
    p.add_argument("--top", type=int, default=60,
                   help="how many distinct candidate drugs to assess")
    p.add_argument("--no-chembl", action="store_true",
                   help="use only the table; make no network request")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--include-local-only", action="store_true",
                   help="also assess drugs approved only for local use. They "
                        "are left out unless the run's target tissue is skin")
    p.add_argument("--include-non-fda", action="store_true",
                   help="also assess compounds that are not in the FDA "
                        "Orange Book. By default only FDA approved drugs "
                        "are assessed, since a research compound has no "
                        "human dose to compare with")
    p.add_argument("--out", default=None,
                   help="output folder; next to the results by default")
    args = p.parse_args(argv)
    if args.top < 1:
        p.error("--top must be at least 1")

    res_path = Path(args.results)
    if not res_path.exists():
        p.error(f"{res_path} not found")
    data = json.loads(res_path.read_text(encoding="utf-8"))
    gene = (data.get("query") or {}).get("gene", "?")

    table = _default_tables()
    if args.cmax_table:
        table.update(read_cmax_table(args.cmax_table))   # yours win
    from . import localdrugs as ld
    library = ld.library(args.library)
    client = None
    if not args.no_chembl and not args.offline:
        from .http import Client
        client = Client()

    from . import fda
    fda_on = not args.include_non_fda and fda.available()
    if not args.include_non_fda and not fda_on:
        log.warning("no Orange Book on disk, so the FDA filter was not "
                    "applied: %s", fda.NOTE_NO_TABLE)
    q = data.get("query") or {}
    hide_local = (fda_on and not args.include_local_only
                  and not fda.target_is_skin(q.get("tissue"),
                                             q.get("activity_tissue")))
    out_dir = Path(args.out) if args.out else res_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    blocks = []
    for blk in data.get("results") or []:
        blocks.append(assess(blk, client=client, table=table,
                             library=library, top=args.top,
                             fda_only=fda_on, hide_local_only=hide_local))
    merged = blocks[0] if len(blocks) == 1 else _merge(blocks)
    (out_dir / "exposure.json").write_text(
        json.dumps(merged, indent=2, default=str), encoding="utf-8")
    text = markdown(merged, gene)
    (out_dir / "exposure.md").write_text(text, encoding="utf-8")
    print(text)
    print(f"Wrote {out_dir / 'exposure.json'} and {out_dir / 'exposure.md'}")
    return 0


def _merge(blocks: list) -> dict:
    """One result for a run with several promoters: a drug keeps its best
    verdict across them, since the exposure does not depend on the promoter."""
    rank = {v: i for i, v in enumerate(VERDICTS)}
    drugs: dict = {}
    for b in blocks:
        for name, d in b["drugs"].items():
            if name not in drugs or rank[d["verdict"]] < rank[
                    drugs[name]["verdict"]]:
                drugs[name] = d
    counts = {v: sum(1 for d in drugs.values() if d["verdict"] == v)
              for v in VERDICTS}
    return {"drugs": drugs, "counts": counts,
            "errors": [e for b in blocks for e in b["errors"]],
            "n_assessed": len(drugs)}


if __name__ == "__main__":
    raise SystemExit(main())
