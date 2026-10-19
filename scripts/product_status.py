"""Say whether a docking run of the gene's own protein produced a result.

The launcher asks this after the product dock so it can offer the exploratory
mode only when the pocket was refused. It prints one line, "usable UNIPROT" or
"refused UNIPROT" (UNIPROT is "-" when the record has none), and always exits
0, so a batch file can read the line with `for /f`.

A missing or unreadable file counts as refused. That is the safe reading: the
exploratory mode is offered, and the person running it decides.

    python scripts/product_status.py results/ABCA4_up/docking/product.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def status(path) -> tuple[str, str]:
    """Return ("usable" or "refused", UniProt accession or "-")."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (OSError, ValueError):
        return "refused", "-"
    if isinstance(data, list):
        data = data[0] if data else {}
    if not isinstance(data, dict):
        return "refused", "-"
    uniprot = str(data.get("uniprot") or "-").split()[0] or "-"
    return ("usable" if data.get("usable") else "refused"), uniprot


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        print("usage: product_status.py PRODUCT.json", file=sys.stderr)
        return 2
    state, uniprot = status(args[0])
    print(f"{state} {uniprot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
