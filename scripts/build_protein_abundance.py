"""Turn a drug proteomics export into the table generoutes reads.

    python scripts/build_protein_abundance.py decryptE.csv
    python scripts/build_protein_abundance.py decryptE.csv ptm.csv
        --out data/drug_protein_abundance.tsv

Why this exists
---------------
`protein_abundance` is the only route that asks about the gene's own
PROTEIN. Every other abundance route reads the transcript, so a
stabiliser or a degradation blocker - which raises protein with the
transcript flat - is invisible without it. The route is on by default and
needs no flag; it needs the file, and nothing fetches that file, because
decryptE and ProTargetMiner are manual downloads with a column layout of
their own.

What it does and does not do
----------------------------
It renames columns and nothing else. Every input row becomes an output
row, because `proteomics.load` pools rows per compound and refuses a
compound with fewer than `MIN_OBSERVATIONS` of them: averaging replicates
here would turn two observations into one and silently drop the
compound.

A linear fold-change column is converted to log2. A non-positive one has
no logarithm, so it is counted and dropped rather than guessed at.

It will not invent a gene symbol. `proteomics.load` matches on the
symbol, so an export carrying only UniProt accessions is refused by name
rather than written out as a table that matches nothing.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# What the loader reads. `dose` and `time` are carried through because the
# format documents them, though only the first four are read today.
OUT_COLUMNS = ("drug", "protein", "log2fc", "cell_line", "dose", "time",
               "source")

# Column spellings seen in the two exports this is written for, plus the
# obvious neighbours. Matched lowercased with spaces and dots folded to
# underscores, so "Log2 Fold Change" and "log2_fold_change" are one name.
DRUG_COLUMNS = ("drug", "drug_name", "compound", "compound_name",
                "perturbagen", "treatment", "pert_iname", "inhibitor")
GENE_COLUMNS = ("gene", "gene_name", "gene_symbol", "symbol", "protein",
                "protein_name", "genes")
LOG2_COLUMNS = ("log2fc", "log2_fc", "logfc", "log2foldchange",
                "log2_fold_change", "log2_ratio", "effect_size")
LINEAR_COLUMNS = ("fold_change", "foldchange", "fc", "ratio")
CELL_COLUMNS = ("cell_line", "cell", "cellline", "cell_type", "cell_name")
DOSE_COLUMNS = ("dose", "concentration", "conc", "dose_um")
TIME_COLUMNS = ("time", "timepoint", "time_point", "duration")
SOURCE_COLUMNS = ("source", "dataset", "study", "experiment")

# Accession columns, named only so the refusal can say what it found when
# a gene symbol is missing. An accession is not a symbol and this script
# has no table to turn one into the other.
ACCESSION_COLUMNS = ("uniprot", "uniprot_id", "accession", "protein_id",
                     "majority_protein_ids")


def norm(name: str) -> str:
    """A column name, with spellings that differ only in punctuation folded."""
    out = (name or "").strip().lower()
    for ch in (" ", ".", "-", "/"):
        out = out.replace(ch, "_")
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def pick(header: list, names: tuple) -> str | None:
    """The first column in `header` whose normalised name is in `names`."""
    seen = {norm(h): h for h in header if h}
    for want in names:
        if want in seen:
            return seen[want]
    return None


def resolve_columns(header: list) -> dict:
    """Which input column feeds which output column.

    Raises `ValueError` naming what was wanted and what the file has, so
    a refusal says which column to add rather than that something is
    wrong somewhere.
    """
    found = {
        "drug": pick(header, DRUG_COLUMNS),
        "protein": pick(header, GENE_COLUMNS),
        "log2fc": pick(header, LOG2_COLUMNS),
        "linear": pick(header, LINEAR_COLUMNS),
        "cell_line": pick(header, CELL_COLUMNS),
        "dose": pick(header, DOSE_COLUMNS),
        "time": pick(header, TIME_COLUMNS),
        "source": pick(header, SOURCE_COLUMNS),
    }
    shown = ", ".join(str(h) for h in header[:12]) or "(no header)"
    if not found["drug"]:
        raise ValueError(
            f"no compound column. Looked for {', '.join(DRUG_COLUMNS)}; "
            f"found {shown}")
    if not found["protein"]:
        accession = pick(header, ACCESSION_COLUMNS)
        extra = (f" This file has {accession}, which is an accession and "
                 f"not a symbol; generoutes matches on the gene symbol, so "
                 f"add a symbol column." if accession else "")
        raise ValueError(
            f"no gene symbol column. Looked for "
            f"{', '.join(GENE_COLUMNS)}; found {shown}.{extra}")
    if not (found["log2fc"] or found["linear"]):
        raise ValueError(
            f"no effect column. Looked for {', '.join(LOG2_COLUMNS)} or "
            f"{', '.join(LINEAR_COLUMNS)}; found {shown}")
    return found


def _open(path: Path):
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_export(path: Path, source: str | None = None) -> tuple:
    """`(rows, counts)` for one export, in the output's own column names."""
    rows, counts = [], {"read": 0, "no_drug": 0, "no_gene": 0,
                        "no_number": 0, "not_positive": 0}
    with _open(path) as fh:
        sample = fh.read(8192)
        fh.seek(0)
        # Same sniff the loader uses: whichever of tab or comma is more
        # common, rather than guessing from the suffix - the Broad ships
        # .txt that is tab separated and people save it as .csv.
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(fh, delimiter=delim)
        cols = resolve_columns(list(reader.fieldnames or []))
        for raw in reader:
            counts["read"] += 1
            drug = (raw.get(cols["drug"]) or "").strip()
            gene = (raw.get(cols["protein"]) or "").strip()
            if not drug:
                counts["no_drug"] += 1
                continue
            if not gene:
                counts["no_gene"] += 1
                continue
            if cols["log2fc"]:
                text = raw.get(cols["log2fc"])
                try:
                    value = float(text)
                except (TypeError, ValueError):
                    counts["no_number"] += 1
                    continue
            else:
                text = raw.get(cols["linear"])
                try:
                    linear = float(text)
                except (TypeError, ValueError):
                    counts["no_number"] += 1
                    continue
                # log2 of zero or less is undefined. A ratio of 0 means
                # the protein was not detected after treatment, which is
                # not the same as a measured fall, so it is dropped and
                # counted rather than turned into a large negative number.
                if linear <= 0:
                    counts["not_positive"] += 1
                    continue
                value = math.log2(linear)
            if not math.isfinite(value):
                counts["no_number"] += 1
                continue
            rows.append({
                "drug": drug,
                "protein": gene.upper(),
                "log2fc": str(round(value, 4)),
                "cell_line": (raw.get(cols["cell_line"]) or "").strip()
                if cols["cell_line"] else "",
                "dose": (raw.get(cols["dose"]) or "").strip()
                if cols["dose"] else "",
                "time": (raw.get(cols["time"]) or "").strip()
                if cols["time"] else "",
                "source": ((raw.get(cols["source"]) or "").strip()
                           if cols["source"] else "") or source or path.stem,
            })
    return rows, counts


def write_table(rows: list, out: Path) -> Path:
    """Write the table, through a .part file so a crash leaves nothing.

    `proteomics.available()` is true for any non-empty file at this path,
    so a half-written one would be read as the table for ever.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".part")
    # csv.writer, because a compound name carrying a tab would otherwise
    # split the row and shift every column after it.
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n",
                       quoting=csv.QUOTE_MINIMAL)
        w.writerow(OUT_COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in OUT_COLUMNS])
    tmp.replace(out)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python scripts/build_protein_abundance.py",
        description="Convert a drug proteomics export (decryptE, "
                    "ProTargetMiner) into the table the protein abundance "
                    "route reads.",
        epilog="Replicates are kept: the route needs at least two "
               "observations of a compound before it reports it.")
    p.add_argument("exports", nargs="+",
                   help="one or more exports (.csv, .tsv, .txt, .gz)")
    p.add_argument("--out", default="data/drug_protein_abundance.tsv")
    p.add_argument("--source",
                   help="what to record as the source when the export has "
                        "no source column (default: the file name)")
    args = p.parse_args(argv)

    rows, failed = [], 0
    for name in args.exports:
        path = Path(name)
        if not path.is_file():
            print(f"  {path}: not found", file=sys.stderr)
            failed += 1
            continue
        try:
            got, counts = read_export(path, source=args.source)
        except ValueError as exc:
            print(f"  {path}: {exc}", file=sys.stderr)
            failed += 1
            continue
        dropped = {k: v for k, v in counts.items()
                   if k != "read" and v}
        print(f"  {path.name}: {len(got)} of {counts['read']} row(s)"
              + (", dropped " + ", ".join(f"{v} {k.replace('_', ' ')}"
                                          for k, v in dropped.items())
                 if dropped else ""))
        rows.extend(got)

    if not rows:
        print("nothing usable to write; no table written.", file=sys.stderr)
        return 1

    out = write_table(rows, Path(args.out))
    drugs = len({r["drug"].lower() for r in rows})
    genes = len({r["protein"] for r in rows})
    # The floor the route applies, said here rather than discovered later
    # when the run reports nothing.
    from generoutes.proteomics import MIN_OBSERVATIONS
    per_pair: dict = {}
    for r in rows:
        key = (r["drug"].lower(), r["protein"])
        per_pair[key] = per_pair.get(key, 0) + 1
    usable = sum(1 for n in per_pair.values() if n >= MIN_OBSERVATIONS)
    print(f"wrote {out} - {len(rows)} row(s), {drugs} compound(s), "
          f"{genes} protein(s)")
    print(f"  {usable} of {len(per_pair)} compound/protein pair(s) have the "
          f"{MIN_OBSERVATIONS} observations the route needs; the rest are "
          f"kept but will not be reported")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
