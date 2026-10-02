"""Approved drugs and their known targets, from a file instead of an API.

Why this exists
---------------
The docking and screening paths need three things from ChEMBL: the gene's
ChEMBL target id, the measured ligands at that target for the validation
gate, and the approved-drug library to screen. On a network that blocks or
breaks those calls, none of it runs - and a corporate network is the normal
case, not the exception.

The circuit breaker already stops a dead host costing more than one
request's worth of waiting. That makes a broken endpoint cheap; it does not
make the analysis possible. For that the data has to be on disk.

The one file that covers all three
----------------------------------
The Broad Drug Repurposing Hub publishes a single annotated table: around
6,800 compounds with SMILES, clinical phase, mechanism of action and target
gene symbols. From that one file you get the screening library, the clinical
status to filter it by, and the known ligands per target that the docking
validation gate needs - which is the whole ChEMBL dependency, offline, in a
few megabytes.

  https://repurposing.broadinstitute.org/download
  the file named repurposing_drugs_*.txt plus repurposing_samples_*.txt,
  or the merged export; any of them is accepted.

What a local table cannot give you
----------------------------------
Potency. The Hub records that a compound acts on a target, not how hard it
binds, so an `actives` set from it has no pChEMBL values. Everything that
reads potency - the established-ligand ranking in `screening`, the
`drug_scoring` potency axis - degrades to unranked rather than silently
scoring zero. The docking gate does not need potency, only the identity of
the known binders, so docking is unaffected.

It is also a curated snapshot. A compound approved after the file was
published is absent, and absence here is not evidence that nothing exists.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from .congeners import salt_free_name

log = logging.getLogger(__name__)

DEFAULT_LIBRARY_PATHS = (
    "data/repurposing_drugs.txt", "data/repurposing_hub.tsv",
    "data/approved_library.tsv", "data/approved_library.csv",
)

SOURCE = (
    "Broad Drug Repurposing Hub: https://repurposing.broadinstitute.org/"
    "download - save the drugs table as data/repurposing_drugs.txt. "
    "Any two-column name/SMILES table also works."
)

# Column names these tables actually use, in the order to try them.
_NAME_COLS = ("pert_iname", "name", "drug", "drug_name", "compound",
              "preferred_name", "title")
_SMILES_COLS = ("smiles", "canonical_smiles", "SMILES", "structure")
_TARGET_COLS = ("target", "targets", "gene", "gene_symbol", "target_symbol")
_PHASE_COLS = ("clinical_phase", "phase", "max_phase", "max_clinical_phase")

# Values in a clinical-phase column that mean "in people, approved".
_APPROVED = {"launched", "approved", "4", "4.0", "phase 4", "phase_4"}


def _sniff(path: Path) -> str:
    """Delimiter, from the header rather than the file extension."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = fh.readline()
    counts = {d: head.count(d) for d in ("\t", ",", ";", "|")}
    best = max(counts, key=counts.get)
    return best if counts[best] else "\t"


def _pick(fieldnames, candidates) -> str | None:
    lowered = {(f or "").strip().lower(): f for f in (fieldnames or [])}
    for c in candidates:
        if c.lower() in lowered:
            return lowered[c.lower()]
    return None


def resolve(explicit=None, candidates=DEFAULT_LIBRARY_PATHS) -> str | None:
    for c in ([explicit] if explicit else list(candidates)):
        if c and Path(c).exists():
            return str(c)
    return None


def read_table(path, *, approved_only: bool = True) -> list:
    """`[{name, smiles, targets, phase}]` from a local drug table.

    Rows with no structure are dropped: a name alone cannot be docked or
    fingerprinted, so keeping them would pad the library with compounds
    that silently fail later.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. {SOURCE}")

    delim = _sniff(path)
    out: list = []
    n_rows = n_no_smiles = n_filtered = 0

    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        name_c = _pick(reader.fieldnames, _NAME_COLS)
        smi_c = _pick(reader.fieldnames, _SMILES_COLS)
        tgt_c = _pick(reader.fieldnames, _TARGET_COLS)
        ph_c = _pick(reader.fieldnames, _PHASE_COLS)

        if smi_c is None:
            raise ValueError(
                f"{path} has no SMILES column. Looked for "
                f"{', '.join(_SMILES_COLS)}; found "
                f"{', '.join(reader.fieldnames or [])}. A table without "
                f"structures cannot be screened or docked.")

        for row in reader:
            n_rows += 1
            smi = (row.get(smi_c) or "").strip()
            if not smi:
                n_no_smiles += 1
                continue
            phase = (row.get(ph_c) or "").strip() if ph_c else ""
            if approved_only and phase:
                if phase.strip().lower() not in _APPROVED:
                    n_filtered += 1
                    continue
            targets = []
            if tgt_c:
                raw = (row.get(tgt_c) or "").strip()
                targets = [t.strip().upper() for t in
                           raw.replace(";", "|").replace(",", "|").split("|")
                           if t.strip()]
            name = (row.get(name_c) or "").strip() if name_c else ""
            out.append({"name": name or smi[:24], "smiles": smi,
                        "targets": targets, "phase": phase or None,
                        "chembl_id": None})

    log.info("local drug table %s: %d usable of %d row(s)"
             "%s%s", path.name, len(out), n_rows,
             f", {n_no_smiles} without a structure" if n_no_smiles else "",
             f", {n_filtered} not approved" if n_filtered else "")
    if approved_only and not any(r["phase"] for r in out):
        log.info("no clinical-phase column was found in %s, so nothing was "
                 "filtered on approval status - the whole table is being "
                 "treated as the library", path.name)
    return out


def library(path=None, *, approved_only: bool = True) -> list:
    """The screening library from a local table, or [] when none is on disk."""
    p = resolve(path)
    if not p:
        log.info("no local drug table found; the library needs the network. "
                 "%s", SOURCE)
        return []
    return read_table(p, approved_only=approved_only)


def actives_for(symbol: str, rows: list, *, uniprot: str | None = None,
                with_potency: bool = True) -> list:
    """Known ligands of one target, from the table's own annotation.

    This is what replaces `agonism.fetch_actives` offline. It carries no
    pChEMBL, because the table records that a compound acts on a target and
    not how hard it binds - so `pchembl` is None rather than a fabricated
    number, and anything ranking on potency degrades to unranked.
    """
    want = (symbol or "").strip().upper()
    if not want:
        return []

    # A local potency file, when there is one, answers this better: it
    # carries measured affinity rather than only the fact of an
    # interaction, so the potency axis has something to read.
    if with_potency and uniprot:
        try:
            from . import localpotency as lpot
            if any(lpot.available().values()):
                measured = lpot.actives_for(uniprot)
                if measured:
                    log.info("local potency file: %d measured ligand(s) for "
                             "%s (%s)", len(measured), want, uniprot)
                    return measured
        except Exception as exc:                            # noqa: BLE001
            log.debug("local potency unavailable: %s", str(exc)[:80])

    out = []
    for r in rows:
        if want in (r.get("targets") or []):
            out.append({"chembl_id": r.get("chembl_id"),
                        "name": r.get("name"),
                        "smiles": r["smiles"], "pchembl": None})
    log.info("local table: %d known ligand(s) for %s, without potency",
             len(out), want)
    return out


def summarise(rows: list, symbol: str | None = None) -> str:
    if not rows:
        return ("No local drug table is on disk, so the library and the "
                "known-ligand set both need the network. " + SOURCE)
    n_t = len({t for r in rows for t in (r.get("targets") or [])})
    bits = [f"{len(rows)} compound(s) with structures from a local table, "
            f"annotated against {n_t} target(s)."]
    if symbol:
        n = len(actives_for(symbol, rows))
        bits.append(
            f"{n} of them are recorded as acting on {symbol}"
            + (", which is the validation set the docking gate needs."
               if n else " - too few to validate a pocket against, so "
                         "docking will refuse rather than report."))
    bits.append(
        "A local table carries no potency, so anything that ranks on it is "
        "unranked rather than scored zero.")
    return " ".join(bits)


# ===================================================================== #
# Merging several sources into one library.                             #
# ===================================================================== #
# "There are 20,000 FDA-approved drugs and I only have 3,000" is the usual
# reaction to a library this size, and the premise is a counting difference
# rather than a gap.
#
# The FDA's often-quoted ~20,000 counts approved PRODUCTS: every strength,
# every formulation, every manufacturer of the same molecule listed
# separately. Metformin alone accounts for hundreds of entries. Distinct
# active moieties with a defined structure is a far smaller number -
# roughly 1,500 to 2,000 small molecules ever approved in the US, and about
# 4,000 in ChEMBL at max_phase 4, which also counts drugs approved
# elsewhere and counts salts separately.
#
# So a library of ~3,000 is not a seventh of the drugs. It is most of the
# distinct molecules, and the way to make it genuinely larger is to WIDEN
# THE DEFINITION deliberately - approved anywhere rather than in the US,
# compounds in trials, withdrawn drugs with human safety data - not to hunt
# for missing FDA approvals.
#
# This merges sources and reports what each contributed, so widening is a
# decision with a number attached.

# A withdrawn drug is not a lower phase: it REACHED approval and was then
# pulled. Putting it at -1 sorted it below preclinical, so asking for
# withdrawn drugs dragged in every lab compound. It sits at 4 with a
# separate flag, and `include_withdrawn` decides it.
PHASE_ORDER = {"launched": 4, "approved": 4, "4": 4, "4.0": 4,
               "withdrawn": 4,
               "phase 3": 3, "3": 3, "3.0": 3,
               "phase 2": 2, "2": 2, "2.0": 2,
               "phase 1": 1, "1": 1, "1.0": 1,
               "preclinical": 0, "": 0}
WITHDRAWN = {"withdrawn", "withdrawn from market", "suspended"}


def parent_key(smiles: str) -> str | None:
    """A deduplication key that survives salt and charge differences.

    An InChIKey alone does NOT merge salts: the counterion is part of the
    molecule, so metformin and metformin hydrochloride get entirely
    different keys. Measured, not assumed. So the parent is extracted
    first - largest fragment, uncharged, tautomer-canonicalised - and the
    key is the connectivity layer of that molecule's InChIKey.

    Residual limitation, also measured: a compound drawn in an explicit
    protonation state can still miss. Metformin written with a charged
    guanidinium does not merge with neutral metformin, because a
    delocalised charge cannot be canonically neutralised. Databases
    normally write the neutral parent or a simple salt, so the common cases
    are covered and this one is rare.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:
        return None
    RDLogger.DisableLog("rdApp.*")
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = rdMolStandardize.ChargeParent(mol)
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        key = Chem.MolToInchiKey(mol)
        return key[:14] if key else None
    except Exception:                                       # noqa: BLE001
        return None


def merge(paths, *, min_phase: int = 4, include_withdrawn: bool = False,
          fallback_name_key: bool = True) -> dict:
    """Combine several drug tables into one deduplicated library.

    `min_phase` is the widening dial: 4 keeps approved or launched only, 3
    adds phase 3, 1 adds everything given to a human. `include_withdrawn`
    is separate, because a withdrawn drug reached approval and was then
    pulled - it carries human safety data and a documented reason for
    removal, so it is an opt-in rather than a point on the phase scale.

    Returns the library plus an accounting of what each source contributed
    and what was collapsed, because more rows is only an improvement when
    the extra rows are distinct molecules.
    """
    merged: dict = {}
    per_source: dict = {}
    n_rows = n_no_key = n_below_phase = n_withdrawn = 0

    for path in paths:
        src = Path(path).name
        try:
            rows = read_table(path, approved_only=False)
        except Exception as exc:                            # noqa: BLE001
            # Keep the END of the message: it leads with a long path and
            # a head truncation preserved only that.
            per_source[src] = {
                "error": type(exc).__name__ + ": "
                         + str(exc).replace("\n", " ")[-220:]}
            continue
        added = already = 0
        for r in rows:
            n_rows += 1
            raw_phase = (r.get("phase") or "").strip().lower()
            phase = PHASE_ORDER.get(raw_phase, 0)
            if raw_phase in WITHDRAWN and not include_withdrawn:
                n_withdrawn += 1
                continue
            if phase < min_phase:
                n_below_phase += 1
                continue
            key = parent_key(r["smiles"])
            unparsed = False
            if key is None:
                if not fallback_name_key:
                    n_no_key += 1
                    continue
                # An unparseable structure still has a name, and dropping it
                # silently would lose a real drug. Keyed on the salt-free
                # name instead, and flagged so it is visible.
                stem = salt_free_name((r.get("name") or "").lower())
                if not stem:
                    n_no_key += 1
                    continue
                key = "name:" + stem
                unparsed = True
            if key in merged:
                already += 1
                prev = merged[key]
                if src not in prev["sources"]:
                    prev["sources"].append(src)
                prev["targets"] = sorted(set(prev["targets"])
                                         | set(r.get("targets") or []))
                prev["n_forms"] += 1
                if phase > PHASE_ORDER.get(
                        (prev.get("phase") or "").strip().lower(), 0):
                    prev["phase"] = r.get("phase")
            else:
                added += 1
                merged[key] = {
                    "name": r.get("name"), "smiles": r["smiles"],
                    "targets": sorted(set(r.get("targets") or [])),
                    "phase": r.get("phase"), "chembl_id": r.get("chembl_id"),
                    "sources": [src], "n_forms": 1, "key": key,
                    "structure_unparsed": unparsed,
                    "withdrawn": raw_phase in WITHDRAWN,
                }
        per_source[src] = {"rows": len(rows), "new": added,
                           "already_present": already}

    library = sorted(merged.values(), key=lambda r: (r["name"] or ""))
    collapsed = sum(r["n_forms"] - 1 for r in library)
    return {
        "library": library,
        "n_distinct": len(library),
        "n_rows_read": n_rows,
        "n_salt_forms_collapsed": collapsed,
        "n_below_phase": n_below_phase,
        "n_unkeyable": n_no_key,
        "n_withdrawn_excluded": n_withdrawn,
        "per_source": per_source,
        "min_phase": min_phase,
        "include_withdrawn": include_withdrawn,
    }


def explain_counts(result: dict) -> str:
    """Why the number is what it is."""
    L = [f"{result['n_distinct']} distinct molecules from "
         f"{result['n_rows_read']} rows across "
         f"{len(result['per_source'])} source(s)."]
    if result["n_salt_forms_collapsed"]:
        L.append(
            f"{result['n_salt_forms_collapsed']} row(s) were salt, hydrate "
            f"or duplicate forms of a molecule already present, collapsed "
            f"on the parent structure - screening the hydrochloride and the "
            f"free base as two candidates doubles the compute and finds the "
            f"same thing twice.")
    if result["n_below_phase"]:
        L.append(f"{result['n_below_phase']} row(s) were below the "
                 f"requested clinical phase.")
    if result["n_unkeyable"]:
        L.append(f"{result['n_unkeyable']} row(s) had neither a usable "
                 f"structure nor a name.")
    if result.get("n_withdrawn_excluded"):
        L.append(f"{result['n_withdrawn_excluded']} withdrawn drug(s) were "
                 f"excluded. They reached approval and were then pulled, so "
                 f"they carry human safety data and a documented reason for "
                 f"removal; include_withdrawn keeps them.")
    L.append(
        "On the counting: the FDA's often-quoted ~20,000 is approved "
        "PRODUCTS - every strength, formulation and manufacturer of the "
        "same molecule counted separately. Distinct active moieties with a "
        "structure is roughly 1,500 to 2,000 for the US, and about 4,000 in "
        "ChEMBL at max_phase 4, which also counts drugs approved elsewhere. "
        "A library of a few thousand is most of the molecules rather than a "
        "fraction of them.")
    if result["min_phase"] >= 4:
        L.append(
            "To make it genuinely larger, widen the definition on purpose: "
            "min_phase 3 adds phase 3 and min_phase 1 adds everything "
            "given to a human, while include_withdrawn keeps drugs that "
            "were approved and later pulled. Each step changes what a hit "
            "means, so it is a decision rather than a default.")
    return " ".join(L)


def write_table(library: list, path) -> Path:
    """Write a merged library back out as a table this package can read."""
    import csv as _csv

    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh, delimiter="\t")
        w.writerow(["name", "clinical_phase", "target", "smiles",
                    "sources", "n_forms"])
        for r in library:
            w.writerow([r.get("name") or "", r.get("phase") or "",
                        "|".join(r.get("targets") or []), r["smiles"],
                        "|".join(r.get("sources") or []), r.get("n_forms", 1)])
    return path


def main(argv=None) -> int:
    """`python -m promoter_tf_drugs.localdrugs` - inspect or merge tables."""
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.localdrugs",
        description="Inspect a local drug table, or merge several into one "
                    "deduplicated library.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=SOURCE)
    p.add_argument("--merge", nargs="+", metavar="FILE",
                   help="tables to combine. Salt, hydrate and duplicate "
                        "forms collapse onto the parent structure")
    p.add_argument("--out", metavar="FILE",
                   help="write the merged library here, in a format this "
                        "package reads back")
    p.add_argument("--min-phase", type=int, default=4, metavar="N",
                   help="4 approved or launched only, 3 adds phase 3, 1 "
                        "adds anything given to a human. This is the dial "
                        "that genuinely enlarges the library, because there "
                        "are no missing FDA approvals to find (default 4)")
    p.add_argument("--include-withdrawn", action="store_true",
                   help="keep drugs that were approved and later pulled. "
                        "They carry human safety data and a documented "
                        "reason for removal, which is why this is separate "
                        "from --min-phase")
    p.add_argument("--target", metavar="SYMBOL",
                   help="also report how many of the merged compounds are "
                        "annotated against this target")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s")

    if not args.merge:
        found = resolve()
        if not found:
            print("No local drug table found.")
            print()
            print(SOURCE)
            return 1
        rows = read_table(found)
        print(f"{found}: {len(rows)} approved compound(s) with structures")
        print()
        print(summarise(rows, args.target))
        return 0

    res = merge(args.merge, min_phase=args.min_phase,
                include_withdrawn=args.include_withdrawn)
    print(f"Merged library: {res['n_distinct']} distinct molecules")
    print()
    for src, info in res["per_source"].items():
        if "error" in info:
            print(f"  {src:34s} FAILED: {info['error']}")
        else:
            print(f"  {src:34s} {info['rows']:6d} rows -> "
                  f"{info['new']:5d} new, "
                  f"{info['already_present']:5d} already present")
    print()
    print(explain_counts(res))
    if args.target:
        n = sum(1 for r in res["library"]
                if args.target.strip().upper() in (r.get("targets") or []))
        print()
        print(f"{n} of them are annotated against {args.target.upper()}.")
    if args.out:
        write_table(res["library"], args.out)
        print()
        print(f"Wrote {args.out}. Use it with --library on the docking and "
              f"screening paths.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
