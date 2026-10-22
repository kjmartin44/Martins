"""Build the drug-fame table from saved runs.

Why this table exists
---------------------
The retrospective check found that the ranking's apparent skill is mostly
drug popularity: a baseline that only counts how many OTHER genes record a
drug scores AUC 0.93 against CTD, while the routes score 0.53, and among
drugs of equal fame the routes fall to chance. That finding lived in
`docs/checks.md` and a reader of one gene's report never saw it.

This writes `data/ctd_fame.tsv`, which `generoutes.fame` reads so a single
run can say, per candidate, how many other genes the literature records it
as moving. A chemical recorded for every gene measured is telling you
about the literature, not about yours.

The counts come from runs you have already done, in `results/`, because
each one holds its gene's CTD record. No network, no bulk download.

    python scripts/build_fame_table.py
    python scripts/build_fame_table.py --results results --out data/ctd_fame.tsv

One row per chemical, with the number of DISTINCT GENES whose CTD record
names it. Genes are counted once however many runs they have, so eight
reruns of BUB1B do not make a chemical eight times as famous.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from generoutes.fame import MIN_GENES, write_table   # noqa: E402

LINK_KEYS = ("wanted_direction", "opposite_direction", "unclear_direction")


def chemicals_by_gene(results_dir: Path) -> dict:
    """`{gene: {chemical, ...}}` from every saved run that has a CTD record.

    Any direction counts. Fame is about how often the literature studied a
    chemical at all, not about which way it moved anything.
    """
    out: dict = {}
    for path in sorted(results_dir.glob("*/results.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            print(f"  skipped {path}: {type(exc).__name__}", file=sys.stderr)
            continue
        gene = ((data.get("query") or {}).get("gene") or "").upper()
        if not gene:
            continue
        names = set()
        for block in data.get("results") or []:
            links = ((block.get("improvements") or {})
                     .get("reported_links")) or {}
            for key in LINK_KEYS:
                for row in links.get(key) or []:
                    if row.get("chemical"):
                        names.add(row["chemical"])
        if names:
            out.setdefault(gene, set()).update(names)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="build_fame_table",
        description="Count how many genes each CTD chemical is recorded for.")
    p.add_argument("--results", default="results",
                   help="folder of saved run directories")
    p.add_argument("--out", default="data/ctd_fame.tsv")
    args = p.parse_args(argv)

    by_gene = chemicals_by_gene(Path(args.results))
    if len(by_gene) < MIN_GENES:
        print(f"only {len(by_gene)} gene(s) have a CTD record in "
              f"{args.results}; {MIN_GENES} are needed before a count means "
              f"anything. No table written.", file=sys.stderr)
        return 1

    counts: dict = {}
    for names in by_gene.values():
        for n in names:
            counts[n] = counts.get(n, 0) + 1

    out = Path(args.out)
    write_table(out, counts, n_genes=len(by_gene))
    top = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:5]
    print(f"wrote {out} - {len(counts)} chemical(s) over {len(by_gene)} gene(s)")
    print("  most recorded: "
          + ", ".join(f"{n} ({c}/{len(by_gene)})" for n, c in top))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
