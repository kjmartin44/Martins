"""Build a table of human peak plasma concentrations from FDA drug labels.

    python scripts/build_fda_cmax.py                       # drugs in results/*/results.json
    python scripts/build_fda_cmax.py --from-results results/ABCA4_up/results.json
    python scripts/build_fda_cmax.py --drugs tretinoin bexarotene
    python scripts/build_fda_cmax.py --all-orange-book     # every FDA ingredient, slow

It asks openFDA for each drug's label, reads the peak concentration sentences in
the pharmacokinetics section (see `generoutes/fdalabel.py` for what it
accepts and refuses), and writes two files:

    data/fda_label_cmax.tsv          drugs with a value it is willing to stand behind
    data/fda_label_cmax_review.tsv   drugs it set aside, with why and the sentence

The first is read by `python -m generoutes.exposure` when it is present.
Rows already in it are kept unless the drug is looked up again, so the table
grows run by run. Most labels give no usable number: a table of figures, a
sentence about time to peak with no value, or only an interaction study. Those
drugs get no row and are listed in the review file, which is the honest result
and not a failure.

openFDA allows 1,000 requests a day without a key (240 a minute) and 120,000
with one. Put a key in OPENFDA_API_KEY and it is sent.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from generoutes import exposure, fda, fdalabel, localdrugs  # noqa: E402
from generoutes.http import Client  # noqa: E402
from generoutes.signatures import breadth_key  # noqa: E402

TABLE = ROOT / "data" / "fda_label_cmax.tsv"
REVIEW = ROOT / "data" / "fda_label_cmax_review.tsv"
COLUMNS = ["name", "cmax_nM", "fu", "dose", "source", "note"]


def names_from_results(paths) -> list:
    """FDA approved drug names in saved results, best routes first."""
    seen, out = set(), []
    files = []
    for p in paths:
        p = Path(p)
        files += [p] if p.is_file() else sorted(p.rglob("results.json"))
    for f in files:
        try:
            data = json.loads(f.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            print(f"skipping {f}: {exc}", file=sys.stderr)
            continue
        for blk in data.get("results") or []:
            for _, name in exposure.candidate_drugs(blk, 10_000):
                k = breadth_key(name)
                if k not in seen and fda.is_approved(name):
                    seen.add(k)
                    out.append(name)
    return out


def read_rows(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as fh:
        return {breadth_key(r["name"]): r for r in csv.DictReader(fh, delimiter="\t")
                if r.get("name")}


def write_rows(path: Path, rows: dict, columns: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for k in sorted(rows):
            w.writerow({c: _one_line(rows[k].get(c, "")) for c in columns})


def _one_line(v) -> str:
    return " ".join(str(v).split())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--from-results", nargs="+", default=None, metavar="PATH",
                   help="results.json files, or folders searched for them "
                        "(default: results/)")
    p.add_argument("--drugs", nargs="+", default=None,
                   help="drug names to look up instead")
    p.add_argument("--all-orange-book", action="store_true",
                   help="every active ingredient in the FDA Orange Book")
    p.add_argument("--limit", type=int, default=None,
                   help="look up at most this many drugs")
    p.add_argument("--out", default=str(TABLE))
    p.add_argument("--review", default=str(REVIEW))
    p.add_argument("--library", default=None,
                   help="drug table with structures, for molecular weights; "
                        "the bundled approved library by default")
    p.add_argument("--offline", action="store_true",
                   help="use only what is already cached")
    args = p.parse_args(argv)
    if args.limit is not None and args.limit < 1:
        p.error("--limit must be at least 1")

    if args.drugs:
        names = list(args.drugs)
    elif args.all_orange_book:
        from generoutes import orangebook
        names = [i.title() for i in orangebook.ingredients()]
    else:
        names = names_from_results(args.from_results or [ROOT / "results"])
    if args.limit:
        names = names[:args.limit]
    if not names:
        print("no drug names to look up: give --drugs, --from-results, or run "
              "an analysis first", file=sys.stderr)
        return 1

    if not fda.available() and not args.drugs and not args.all_orange_book:
        print("no Orange Book on disk, so drugs cannot be limited to FDA "
              "approved ones; see orangebook.SOURCE", file=sys.stderr)
        return 1

    lib = {breadth_key(r["name"]): r for r in localdrugs.library(args.library)}
    client = Client(offline=args.offline)
    key = os.environ.get("OPENFDA_API_KEY")
    if key:
        # openFDA takes the key as a query parameter on every call.
        _get = client.get

        def get(url, params=None, **kw):
            return _get(url, dict(params or {}, api_key=key), **kw)
        client.get = get

    table, review = read_rows(Path(args.out)), read_rows(Path(args.review))
    counts = {"row": 0, "review": 0, "none": 0}
    for i, name in enumerate(names, 1):
        row = lib.get(breadth_key(name))
        mw = exposure.molecular_weight(row["smiles"]) if row else None
        res = fdalabel.read_drug(client, name, mw)
        counts[res["status"]] += 1
        k = breadth_key(name)
        if res["status"] == "row":
            table[k] = res["row"]
            review.pop(k, None)
            note = f"{res['row']['cmax_nM']} nM"
        else:
            table.pop(k, None)
            ex = (res["rejected"] or res["accepted"] or [{}])[0]
            review[k] = {"name": name, "status": res["status"],
                         "reason": res["reason"],
                         "example": ex.get("sentence", "")}
            note = res["reason"]
        print(f"[{i}/{len(names)}] {name}: {res['status']} - {note}")
    write_rows(Path(args.out), table, COLUMNS)
    write_rows(Path(args.review), review,
               ["name", "status", "reason", "example"])
    print(f"\n{counts['row']} with a value, {counts['review']} set aside for "
          f"review, {counts['none']} with nothing to read, of {len(names)}.")
    print(f"Wrote {args.out} ({len(table)} drugs) and {args.review}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
