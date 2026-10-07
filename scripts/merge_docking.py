"""Combine several docking runs into one file the analysis can read.

`--docking-results` takes a single path, and `docking.load_results`
accepts either one run or a list of them. A workup that docks the gene's
own product and then two or three transcription factors has one file
each, so this concatenates them into the list form.

Runs that refused are kept, not dropped. A pocket that could not be
validated is a result - it says the scores for that target would carry
no information - and the report has a place to say so. Silently omitting
it would leave the reader unable to tell a target that was never tried
from one that was tried and refused.

    python scripts/merge_docking.py out.json a.json b.json c.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def expand(paths) -> list[Path]:
    r"""Resolve each argument, expanding any glob.

    cmd.exe does not expand wildcards for a program, so a batch file
    passing `dir\*.json` hands over that literal string. Expanding here
    means the caller does not need a shell that does it.

    Write the output somewhere the pattern does not match, or a second
    run folds the first run's merged file back into itself.
    """
    out: list[Path] = []
    for raw in paths:
        p = Path(raw)
        if any(ch in str(raw) for ch in "*?["):
            parent = p.parent if str(p.parent) else Path(".")
            out.extend(sorted(parent.glob(p.name)))
        else:
            out.append(p)
    return out


def merge(paths) -> list[dict]:
    """Every run in `paths`, flattened, in the order given."""
    runs: list[dict] = []
    for path in expand(paths):
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            print(f"skipping {path.name}: {exc}", file=sys.stderr)
            continue
        got = data if isinstance(data, list) else [data]
        for run in got:
            if isinstance(run, dict) and "usable" in run:
                runs.append(run)
            else:
                print(f"skipping a record in {path.name}: not a docking "
                      f"run", file=sys.stderr)
    return runs


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python scripts/merge_docking.py",
        description=__doc__.split("\n")[0])
    p.add_argument("out", help="where to write the combined file")
    p.add_argument("inputs", nargs="+", help="the runs to combine")
    args = p.parse_args(argv)

    runs = merge(args.inputs)
    if not runs:
        print("no docking runs to merge", file=sys.stderr)
        return 1

    Path(args.out).write_text(json.dumps(runs, indent=2), encoding="utf-8")
    usable = sum(1 for r in runs if r.get("usable"))
    targets = ", ".join(str(r.get("target") or "?") for r in runs)
    print(f"{len(runs)} run(s) -> {args.out}: {targets}")
    print(f"  {usable} validated, {len(runs) - usable} refused and kept "
          f"so the report can say which")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
