"""List the tissues in an expression table and resolve a pick to one name.

The launcher uses this so the person choosing a tissue for the factor activity
test sees the names the table spells. A misspelt tissue used to log a warning
and carry on, and the activity route simply never produced anything.

    python scripts/tissues.py data/rna_tissue_consensus.tsv
    python scripts/tissues.py data/rna_tissue_consensus.tsv --pick 12
    python scripts/tissues.py data/rna_tissue_consensus.tsv --pick "lymph node"

With no `--pick` it prints a numbered list. With `--pick` it prints the one
resolved name and nothing else, so a batch file can read it with `for /f`. A
pick is a number from the list, a name in any case, or a part of a name that
matches exactly one tissue. When nothing or more than one tissue matches it
prints nothing and exits 1, with the reason on stderr.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


def read_tissues(path) -> list[str]:
    """Distinct tissue names, sorted, from the table's Tissue column."""
    p = Path(path)
    with open(p, "rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    opener = (lambda: gzip.open(p, "rt", encoding="utf-8", errors="replace")
              if gz else open(p, "r", encoding="utf-8-sig", errors="replace"))
    seen: set[str] = set()
    with opener() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None) or []
        names = [h.strip().lower() for h in header]
        if "tissue" not in names:
            raise ValueError(f"{p.name} has no Tissue column")
        i = names.index("tissue")
        for row in reader:
            if len(row) > i and row[i].strip():
                seen.add(row[i].strip())
    return sorted(seen, key=str.lower)


def pick(tissues: list[str], text: str) -> str:
    """Resolve `text` to one tissue, or raise ValueError saying why not."""
    text = (text or "").strip()
    if not text:
        raise ValueError("nothing entered")
    if text.isdigit():
        n = int(text)
        if 1 <= n <= len(tissues):
            return tissues[n - 1]
        raise ValueError(f"{n} is not a number between 1 and {len(tissues)}")
    want = text.lower()
    for t in tissues:
        if t.lower() == want:
            return t
    part = [t for t in tissues if want in t.lower()]
    if len(part) == 1:
        return part[0]
    if not part:
        raise ValueError(f"no tissue matches {text!r}")
    raise ValueError(f"{text!r} matches {len(part)} tissues: "
                     + ", ".join(part[:6]) + (", ..." if len(part) > 6 else ""))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("table")
    p.add_argument("--pick", default=None)
    args = p.parse_args(sys.argv[1:] if argv is None else argv)
    try:
        tissues = read_tissues(args.table)
    except (OSError, ValueError) as e:
        print(f"tissues: {e}", file=sys.stderr)
        return 2
    if not tissues:
        print("tissues: the table lists no tissues", file=sys.stderr)
        return 2
    if args.pick is None:
        width = len(str(len(tissues)))
        for n, t in enumerate(tissues, 1):
            print(f"  {n:>{width}}  {t}")
        return 0
    try:
        print(pick(tissues, args.pick))
    except ValueError as e:
        print(f"tissues: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
