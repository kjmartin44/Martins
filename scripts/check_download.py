"""Check a downloaded file against its published MD5 checksum.

The launcher downloads the DepMap and PRISM tables, which are hundreds of
megabytes, and a connection that drops partway or a portal page saved in place
of the file would otherwise look like data. It downloads to a ".part" name,
runs this, and renames the file only when the checksum matches.

    python scripts/check_download.py data/CRISPRGeneEffect.csv.part d10d7369...

Exit status is 0 when the checksum matches, 1 when it does not or the file
cannot be read, and 2 for bad arguments.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path


def md5_of(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("usage: check_download.py FILE MD5", file=sys.stderr)
        return 2
    path, want = args[0], args[1].strip().lower()
    try:
        got = md5_of(path)
    except OSError as e:
        print(f"check_download: {e}", file=sys.stderr)
        return 1
    if got != want:
        print(f"check_download: {Path(path).name} has checksum {got}, "
              f"expected {want}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
