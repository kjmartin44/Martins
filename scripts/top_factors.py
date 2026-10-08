"""Which transcription factors are worth docking against.

Docking one factor against a drug library is twenty minutes, so this
picks the few where it can produce something, rather than every factor
with a motif match - a promoter scan returns far more factors than are
really there, and most of a hit list is chance.

Three requirements, each for a reason:

**Measured occupancy.** A factor with no ChIP peak over this window has
a sequence match and nothing else. Docking against it answers "could
this drug bind a protein that may not be here", which is not a question
worth twenty minutes.

**A resolved role.** Activator or repressor decides which conformation
you want stabilised, and `docking` derives that from the role and your
direction. Without one it cannot choose, and docks whichever pocket
exists rather than the one that would move your gene.

**Not the gene itself.** The gene's own product is a separate and
better-posed run, with `--product`.

Prints one symbol per line, best first, so a shell loop can read it.

    python scripts/top_factors.py results/CCR5_down/results.json --top 3
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


# What a factor's role looks like when it was looked for and not found.
# annotate.py stores the STRING "unknown" for it, which is truthy, so the
# filter below used to test `not tf.get("role")` and let every unresolved
# factor straight through - on a real CCR5 result, 39 of 106 factors, and
# `require_role=True` returned exactly what `require_role=False` did.
UNRESOLVED_ROLES = {"", "unknown", "?", "none"}


def has_role(tf: dict) -> bool:
    """Is this factor's role resolved to something docking can use?"""
    role = tf.get("role")
    return role is not None and str(role).strip().lower() not in UNRESOLVED_ROLES


def pick(result: dict, *, top: int = 3, min_peaks: int = 1,
         require_role: bool = True) -> list[dict]:
    """The factors worth docking, best first."""
    gene = ((result.get("query") or {}).get("gene") or "").strip().upper()
    seen: dict[str, dict] = {}
    for block in result.get("results") or []:
        for tf in block.get("tfs") or []:
            sym = (tf.get("symbol") or "").strip().upper()
            if not sym or sym == gene:
                continue
            if (tf.get("chip_peaks") or 0) < min_peaks:
                continue
            if require_role and not has_role(tf):
                continue
            # The same factor can appear under several transcripts.
            best = seen.get(sym)
            if best is None or (tf.get("score") or 0) > (best.get("score") or 0):
                seen[sym] = tf
    ranked = sorted(seen.values(), key=lambda t: -(t.get("score") or 0))
    return ranked[:top]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python scripts/top_factors.py",
        description=__doc__.split("\n")[0])
    p.add_argument("results", help="results.json from a previous run")
    p.add_argument("--symbols-only", action="store_true",
                   help="print just the symbol, without its role")
    p.add_argument("--top", type=int, default=3,
                   help="how many factors to name (default 3)")
    p.add_argument("--min-peaks", type=int, default=1,
                   help="ChIP peaks required over the window. 0 includes "
                        "factors with a motif match and no measured "
                        "evidence that they are there (default 1)")
    p.add_argument("--any-role", action="store_true",
                   help="include factors whose role is unresolved. Their "
                        "conformation then cannot be chosen, so the dock "
                        "answers whether a drug binds rather than which "
                        "way it would push")
    p.add_argument("--out", metavar="FILE",
                   help="write the symbols here, one per line, instead of "
                        "to stdout. Preferred over a shell redirect: cmd "
                        "and PowerShell disagree about byte-order marks, "
                        "and a BOM on the first line turns the first "
                        "symbol into one no lookup resolves")
    p.add_argument("--verbose", action="store_true",
                   help="also print the role, peaks and score, to stderr")
    args = p.parse_args(argv)

    path = Path(args.results)
    if not path.is_file():
        print(f"no such file: {path}", file=sys.stderr)
        return 1
    try:
        result = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        print(f"could not read {path}: {exc}", file=sys.stderr)
        return 1

    chosen = pick(result, top=args.top, min_peaks=args.min_peaks,
                  require_role=not args.any_role)
    if not chosen:
        print("no factor met the bar for docking "
              f"(ChIP peaks >= {args.min_peaks}"
              + ("" if args.any_role else ", with a resolved role") + ")",
              file=sys.stderr)
        return 1

    for tf in chosen:
        if args.verbose:
            print(f"  {tf['symbol']:10} role={tf.get('role') or '?':10} "
                  f"peaks={tf.get('chip_peaks')} "
                  f"score={tf.get('score')}", file=sys.stderr)

    # "SYMBOL role" per line. The role travels with the symbol because
    # the docking CLI resolves it from local annotation, which often
    # comes back unknown for a factor the scan knows the role of - and
    # without it the dock is not filtered to a conformation at all,
    # which is the entire reason to dock a factor rather than a product.
    if args.symbols_only:
        symbols = [tf["symbol"] for tf in chosen]
    else:
        symbols = [f"{tf['symbol']} {tf.get('role') or 'unknown'}"
                   for tf in chosen]
    if args.out:
        # Plain UTF-8, no BOM, newline-terminated. A shell redirect is
        # not reliable here across cmd and PowerShell.
        text = "\n".join(symbols) + "\n"
        Path(args.out).write_text(text, encoding="utf-8")
    else:
        for sym in symbols:
            print(sym)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
