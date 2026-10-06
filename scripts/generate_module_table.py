"""Regenerate the module table in docs/architecture.md.

The table lists every module with its line count and the first line of
its own docstring. Written by a script rather than by hand because a
hand-written list of sixty-one modules is out of date the week after it
is written, and a map that is quietly wrong is worse than none.

    python scripts/generate_module_table.py          # rewrite the table
    python scripts/generate_module_table.py --check   # fail if it is stale

The --check form runs as a test.
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "promoter_tf_drugs"
DOC = ROOT / "docs" / "architecture.md"

START = "| Module | Lines | Purpose |"
END_MARKER = "Line counts are indicative"

SKIP = {"__init__", "__main__"}


def rows() -> list[tuple[str, int, str]]:
    out = []
    for f in sorted(PKG.glob("*.py")):
        if f.stem in SKIP:
            continue
        text = f.read_text(encoding="utf-8")
        try:
            doc = ast.get_docstring(ast.parse(text)) or ""
        except SyntaxError:
            doc = ""
        first = doc.strip().split("\n")[0].strip() or "(no docstring)"
        out.append((f.name, len(text.splitlines()), first))
    return out


def table() -> str:
    lines = [START, "| --- | --- | --- |"]
    for name, n, purpose in rows():
        # A pipe inside a description would break the table.
        lines.append(f"| `{name}` | {n} | {purpose.replace('|', '/')} |")
    return "\n".join(lines)


def splice(doc_text: str, new_table: str) -> str:
    i = doc_text.index(START)
    j = doc_text.index(END_MARKER, i)
    return doc_text[:i] + new_table + "\n\n" + doc_text[j:]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--check", action="store_true",
                   help="exit non-zero if the table is stale, changing "
                        "nothing")
    args = p.parse_args(argv)

    current = DOC.read_text(encoding="utf-8")
    updated = splice(current, table())

    if args.check:
        if current == updated:
            print("module table is current")
            return 0
        print("module table in docs/architecture.md is stale; run "
              "python scripts/generate_module_table.py", file=sys.stderr)
        return 1

    if current == updated:
        print("module table already current")
        return 0
    DOC.write_text(updated, encoding="utf-8")
    print(f"rewrote the module table: {len(rows())} modules")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
