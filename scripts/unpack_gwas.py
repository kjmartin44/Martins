"""Unpack the GWAS Catalog associations zip into the table the package reads.

EBI ships the associations as a zip, and the package reads a plain table, so
the launcher downloads the zip and runs this. It picks the table inside the
zip, checks that it has the gene and trait columns the loader needs, and
writes it to the destination. Nothing is written when the check fails, so a
captive portal page or a changed layout leaves no half-valid file behind.

    python scripts/unpack_gwas.py gwas-catalog-associations-full.zip \
        data/gwas_catalog_associations.tsv

Exit status is 0 on success and 1 when the zip holds no usable table.
"""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path

GENE_COLUMNS = ("mapped_gene", "reported gene(s)", "mapped gene")
TRAIT_COLUMNS = ("disease/trait", "mapped_trait", "trait")


def usable_header(line: str) -> bool:
    cols = {c.strip().lower() for c in line.rstrip("\r\n").split("\t")}
    return (any(c in cols for c in GENE_COLUMNS)
            and any(c in cols for c in TRAIT_COLUMNS))


def unpack(zip_path, dest) -> str:
    """Write the associations table to `dest` and return the member used."""
    with zipfile.ZipFile(zip_path) as zf:
        members = [i for i in zf.infolist()
                   if not i.is_dir()
                   and i.filename.lower().endswith((".tsv", ".txt"))]
        for info in sorted(members, key=lambda i: -i.file_size):
            with zf.open(info) as fh:
                first = fh.readline().decode("utf-8", "replace")
            if not usable_header(first):
                continue
            dest = Path(dest)
            tmp = dest.with_suffix(dest.suffix + ".part")
            with zf.open(info) as src, open(tmp, "wb") as out:
                while True:
                    chunk = src.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            tmp.replace(dest)
            return info.filename
    raise ValueError("no table with a gene column and a trait column in "
                     f"{Path(zip_path).name}")


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("usage: unpack_gwas.py ZIP DEST", file=sys.stderr)
        return 2
    try:
        member = unpack(*args)
    except (OSError, ValueError, zipfile.BadZipFile) as e:
        print(f"unpack_gwas: {e}", file=sys.stderr)
        return 1
    print(f"unpacked {member} to {args[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
