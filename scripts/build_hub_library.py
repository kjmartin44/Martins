"""Join the Broad Drug Repurposing Hub's two downloads into one table.

The Hub publishes drug annotations (phase, mechanism, targets) and structures
as separate files, and both start with comment lines that begin with "!".
The library loaders need one table with a SMILES column, so this joins them
on the drug name and writes ``data/repurposing_hub.tsv``.

    python scripts/build_hub_library.py \
        data/repurposing_drugs.txt data/repurposing_samples.txt \
        --out data/repurposing_hub.tsv

A drug with several samples takes the first sample that has a structure.
A drug with no structure is dropped, because it cannot be docked or screened.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_hub_table(path):
    """Rows of a Hub download as dicts, skipping the "!" comment lines."""
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        lines = (ln for ln in fh if not ln.startswith("!"))
        return list(csv.DictReader(lines, delimiter="\t"))


def build(drugs_path, samples_path):
    structures = {}
    for row in read_hub_table(samples_path):
        name = (row.get("pert_iname") or "").strip()
        smiles = (row.get("smiles") or "").strip()
        if name and smiles and smiles.lower() != "na":
            structures.setdefault(name.lower(), smiles)
    out, no_structure = [], 0
    for row in read_hub_table(drugs_path):
        name = (row.get("pert_iname") or "").strip()
        smiles = structures.get(name.lower())
        if not name:
            continue
        if not smiles:
            no_structure += 1
            continue
        out.append({"name": name, "smiles": smiles,
                    "clinical_phase": row.get("clinical_phase", ""),
                    "moa": row.get("moa", ""), "target": row.get("target", "")})
    return out, no_structure


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("drugs", help="repurposing_drugs_*.txt")
    p.add_argument("samples", help="repurposing_samples_*.txt")
    p.add_argument("--out", default="data/repurposing_hub.tsv")
    args = p.parse_args(argv)
    for path in (args.drugs, args.samples):
        if not Path(path).exists():
            print(f"{path} not found", file=sys.stderr)
            return 1
    rows, dropped = build(args.drugs, args.samples)
    if not rows:
        print("No drug matched a structure. Check that both files are the "
              "Hub downloads.", file=sys.stderr)
        return 1
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} drugs with structures to {args.out} "
          f"({dropped} dropped for having none).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
