"""Check that a downloaded .gmt file is a gene set library and not an error page.

The launcher downloads the Enrichr knockout libraries and runs this before it
keeps them. A .gmt line is a term, a description field (empty in Enrichr's
files) and then genes, all separated by tabs. A portal page saved in place of
the file has none of that.

    python scripts/check_gmt.py data/knockout_lincs_crispr_ko.gmt.part 1000

Exit status is 0 when the file has at least MIN_LINES lines and nearly all are
well formed, 1 otherwise, and 2 for bad arguments.
"""

from __future__ import annotations

import sys
from pathlib import Path

WELL_FORMED_SHARE = 0.95


def check(path, min_lines: int) -> tuple:
    """`(ok, n_lines, n_good)` for a file."""
    n = good = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            n += 1
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 3 and parts[0] and any(parts[2:]):
                good += 1
    return n >= min_lines and good >= WELL_FORMED_SHARE * n, n, good


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2 or not args[1].isdigit():
        print("usage: check_gmt.py FILE MIN_LINES", file=sys.stderr)
        return 2
    try:
        ok, n, good = check(args[0], int(args[1]))
    except OSError as e:
        print(f"check_gmt: {e}", file=sys.stderr)
        return 1
    if not ok:
        print(f"check_gmt: {Path(args[0]).name} has {n} line(s), {good} of "
              f"them gene sets; need at least {args[1]} and nearly all well "
              f"formed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
