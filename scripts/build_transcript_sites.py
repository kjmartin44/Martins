"""Build `data/transcript_sites.tsv` from the public site catalogues.

Why this table exists
---------------------
`generoutes.utr` scans a transcript's UTRs for AU-rich elements, Pumilio
response elements, DRACH m6A consensus, upstream ORFs and G-quadruplexes.
A scan can say an element is in excess against shuffled sequence; it cannot
say which microRNA or RNA-binding protein uses it, and the name is what
makes the finding actionable - a microRNA has an antagomir, an RBP has a
knockdown. `utr.annotated_sites` reads this table to supply those names.

What goes in a row
------------------
    gene, factor, kind, region, evidence, source

`kind` groups the output ("microRNA", "RBP", "m6A"). `evidence` is what
separates a prediction from a measurement, and `utr.summarise_sites` counts
those apart, so every row carries the honest value for where it came from:
TargetScan rows are `predicted`, miRTarBase rows carry that entry's own
experiment list (reporter assay, western blot, CLIP), eCLIP rows say the
peak was called.

Inputs, each a bulk download under its own licence. This script does not
fetch anything; pass the files you downloaded.

    python scripts/build_transcript_sites.py \
        --targetscan Predicted_Targets_Info.txt \
        --mirtarbase hsa_MTI.csv \
        --eclip eclip_peaks/ --annotation data/gencode.basic.gtf.gz \
        --genes my_genes.txt

ENCODE eCLIP IDR peaks are genomic intervals, so they are mapped to gene
symbols against the GENCODE basic annotation passed as `--annotation`
(`data/gencode.basic.gtf.gz`, the same file the neighbour check uses), by
overlap on the same strand with the gene body, introns included - so those
rows say `gene body` rather than claiming the protein sits in a UTR. ENCODE
leaves the BED name column as ".", so save each download as
`<RBP>_<cell line>_<accession>.bed.gz` and the RBP is read from the file
name. A folder stands for the BED files in it, because naming all 252 peak
files is longer than a command line may be. `--m6a` takes an atlas of
methylation sites in the same shape and is mapped the same way.

Choosing which genes to keep
----------------------------
`utr.annotated_sites` reads the whole table to answer about one gene, so
its size is a cost every run pays. Two choices keep it near 115 MB rather
than 880 MB, and both are about what the rows mean, not only about size:

- miRTarBase's `Functional MTI (Weak)` class is a hit from a
  high-throughput screen, not a direct assay, and in release 11 it is
  7,309,069 of 7,339,072 human rows. It is left out unless you pass
  `--mirtarbase-weak`, which puts it back with every row labelled weak.
  `Non-Functional MTI` is always left out: it records that the pair does
  not act, so listing it as a regulator would say the opposite.
- `--genes` takes a file of gene symbols, one per line, and drops
  everything outside it. With no `--genes` the whole input is written.
  TargetScan's conserved-family file alone is millions of rows, so a
  table meant to live beside one project's runs usually wants this.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

COLUMNS = ("gene", "factor", "kind", "region", "evidence", "source")
BIN = 100_000


def _open(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace", newline="")


def read_genes(path) -> set:
    """Gene symbols to keep, one per line; blank and `#` lines ignored."""
    keep = set()
    with _open(path) as fh:
        for line in fh:
            name = line.split("\t")[0].strip().upper()
            if name and not name.startswith("#"):
                keep.add(name)
    return keep


def _pick(row: dict, *names):
    """First non-empty value among `names`, matched on lowercased headers."""
    for n in names:
        v = row.get(n)
        if v and str(v).strip():
            return str(v).strip()
    return ""


def _lower_keys(row: dict) -> dict:
    return {(k or "").strip().lower(): v for k, v in row.items()}


# ------------------------------------------------------------- TargetScan
def targetscan_rows(path):
    """Predicted conserved microRNA family sites, human rows only."""
    with _open(path) as fh:
        for raw in csv.DictReader(fh, delimiter="\t"):
            row = _lower_keys(raw)
            tax = _pick(row, "gene tax id", "species id")
            if tax and tax != "9606":
                continue
            gene = _pick(row, "gene symbol", "gene_symbol")
            factor = _pick(row, "mir family", "mirna", "mir id",
                           "representative mirna")
            if gene and factor:
                yield {"gene": gene.upper(), "factor": factor,
                       "kind": "microRNA", "region": "3'UTR",
                       "evidence": "predicted", "source": "TargetScan"}


# ------------------------------------------------------------- miRTarBase
def _xlsx_rows(path):
    from openpyxl import load_workbook

    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        rows = wb[wb.sheetnames[0]].iter_rows(values_only=True)
        header = [(c or "") for c in next(rows, [])]
        for values in rows:
            yield dict(zip(header, values))
    finally:
        wb.close()


def mirtarbase_rows(path, *, include_weak: bool = False):
    """microRNA targets with experimental support, human rows only.

    The entry's own experiment list becomes the evidence value, so a
    reporter assay is not reported as the same thing as a CLIP hit.

    Two classes are left out by default. A `Non-Functional MTI` is a
    measurement that the pair does NOT act, so listing it as a regulator
    would invert what it says. `Functional MTI (Weak)` is a hit from a
    high-throughput screen rather than a direct assay, and in release 11
    it is 7,309,069 of the 7,339,072 human rows against 17,377 direct
    ones - keeping it makes the table 880 MB, which every run then reads
    to answer about one gene. `--mirtarbase-weak` puts it back, labelled.
    """
    if str(path).lower().endswith((".xlsx", ".xlsm")):
        raw_rows = _xlsx_rows(path)
    else:
        with _open(path) as fh:
            sample = fh.read(8192)
            fh.seek(0)
            delim = "\t" if sample.count("\t") >= sample.count(",") else ","
            raw_rows = list(csv.DictReader(fh, delimiter=delim))
    for raw in raw_rows:
        row = _lower_keys(raw)
        # The xlsx spells the species out, the csv release writes "hsa".
        species = _pick(row, "species (target gene)").lower()
        if species and species not in ("homo sapiens", "hsa", "human"):
            continue
        gene = _pick(row, "target gene")
        factor = _pick(row, "mirna")
        if not (gene and factor):
            continue
        support = _pick(row, "support type").lower()
        if "non-functional" in support:
            continue
        weak = "weak" in support
        if weak and not include_weak:
            continue
        evidence = _pick(row, "experiments", "support type") or "reported"
        evidence = evidence.replace("//", ", ")
        yield {"gene": gene.upper(), "factor": factor, "kind": "microRNA",
               "region": "",
               "evidence": evidence + (" (weak)" if weak else ""),
               "source": "miRTarBase"}


# ------------------------------------------------- genomic intervals -> gene
def gene_index(annotation_path) -> dict:
    """`{(chrom, strand, bin): [(start, end, symbol), ...]}` from a GTF.

    The GTF parse is `neighbours.load_gtf`; only the lookup is built here,
    because its own `genes_in_window_local` scans every gene per query and
    an eCLIP file holds a hundred thousand peaks.
    """
    from generoutes import neighbours

    neighbours.load_gtf(annotation_path)
    index: dict = {}
    for chrom, start, end, strand, name, _biotype in neighbours._GTF_GENES:
        for b in range(start // BIN, end // BIN + 1):
            index.setdefault((chrom, strand, b), []).append(
                (start, end, name))
    return index


def genes_at(index: dict, chrom: str, start: int, end: int, strand: int):
    """Symbols of genes the interval overlaps on the same strand.

    `start`/`end` are BED coordinates, 0-based with an exclusive end; the
    GTF is 1-based inclusive, so the interval is compared as [start+1, end].
    """
    chrom = chrom.replace("chr", "")
    lo, hi = start + 1, end
    out = set()
    for b in range(start // BIN, end // BIN + 1):
        for g_start, g_end, name in index.get((chrom, strand, b), ()):
            if g_start <= hi and g_end >= lo:
                out.add(name.upper())
    return out


def bed_rows(path, index, *, kind, evidence, source, default_factor=""):
    """One row per gene a peak lands in, with the factor that bound it.

    The factor is the BED name column where there is one, then
    `default_factor`, then the first field of the file name - which is how
    an ENCODE download, whose name column is ".", carries its RBP.

    A peak is placed against the whole gene body, introns included, so the
    region is recorded as that rather than as a UTR site.
    """
    with _open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            name = f[3].strip() if len(f) > 3 else ""
            # ENCODE leaves the BED name column as "." on its eCLIP peaks,
            # so the RBP comes from the file name the download was saved
            # under: "<RBP>_<cell line>_<accession>.bed.gz".
            factor = (name.split("_")[0] if name and name != "." else "") \
                or default_factor \
                or Path(str(path)).name.split("_")[0].split(".")[0]
            strand = -1 if len(f) > 5 and f[5].strip() == "-" else 1
            for gene in genes_at(index, f[0], start, end, strand):
                yield {"gene": gene, "factor": factor, "kind": kind,
                       "region": "gene body", "evidence": evidence,
                       "source": source}


# -------------------------------------------------------------------- main
def bed_paths(given) -> list:
    """The BED files named, with a directory standing for the files in it.

    A set of eCLIP peaks is 252 files, and naming them all is longer than a
    command line is allowed to be.
    """
    out = []
    for path in given or ():
        p = Path(path)
        if p.is_dir():
            out.extend(sorted(q for q in p.iterdir()
                              if q.name.endswith((".bed", ".bed.gz"))))
        else:
            out.append(p)
    return out


def build(args) -> list:
    keep = read_genes(args.genes) if args.genes else None
    index = None
    if args.eclip or args.m6a:
        if not args.annotation:
            raise SystemExit("--eclip and --m6a need --annotation "
                             "(a GENCODE GTF) to map intervals to genes")
        index = gene_index(args.annotation)

    streams = []
    if args.targetscan:
        streams.append(targetscan_rows(args.targetscan))
    if args.mirtarbase:
        streams.append(mirtarbase_rows(args.mirtarbase,
                                       include_weak=args.mirtarbase_weak))
    for path in bed_paths(args.eclip):
        streams.append(bed_rows(path, index, kind="RBP",
                                evidence="eCLIP IDR peak",
                                source="ENCODE eCLIP"))
    for path in bed_paths(args.m6a):
        streams.append(bed_rows(path, index, kind="m6A",
                                evidence="atlas peak", source="m6A atlas",
                                default_factor="m6A"))

    seen, rows = set(), []
    for stream in streams:
        for row in stream:
            if keep is not None and row["gene"] not in keep:
                continue
            key = tuple(row[c] for c in COLUMNS)
            if key in seen:
                continue
            seen.add(key)
            rows.append(row)
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="build_transcript_sites",
        description=__doc__.split("\n\n")[0])
    p.add_argument("--targetscan", help="TargetScan Predicted_Targets_Info.txt")
    p.add_argument("--mirtarbase", help="miRTarBase hsa_MTI.xlsx (or a csv/tsv "
                                        "export of it)")
    p.add_argument("--mirtarbase-weak", action="store_true",
                   help="also keep miRTarBase's high-throughput screen hits "
                        "(Functional MTI (Weak)), which is 7.3M of its 7.34M "
                        "human rows")
    p.add_argument("--eclip", nargs="*", default=[],
                   help="ENCODE eCLIP IDR peak BED files, or a folder of them")
    p.add_argument("--m6a", nargs="*", default=[],
                   help="m6A atlas peaks, BED, or a folder of them")
    p.add_argument("--annotation", help="GENCODE GTF, for the BED inputs")
    p.add_argument("--genes", help="file of gene symbols to keep, one per line")
    p.add_argument("--out", default="data/transcript_sites.tsv")
    args = p.parse_args(argv)

    if not (args.targetscan or args.mirtarbase or args.eclip or args.m6a):
        print("nothing to read. Pass at least one of --targetscan, "
              "--mirtarbase, --eclip, --m6a.", file=sys.stderr)
        return 1
    for path in [args.targetscan, args.mirtarbase, args.annotation,
                 *args.eclip, *args.m6a]:
        if path and not Path(path).exists():
            print(f"{path} not found", file=sys.stderr)
            return 1

    rows = build(args)
    if not rows:
        print("no site survived. Check the input files, and the --genes "
              "list if you passed one.", file=sys.stderr)
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(COLUMNS), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    measured = sum(1 for r in rows if r["evidence"] != "predicted")
    by_kind: dict = {}
    for r in rows:
        by_kind[r["kind"]] = by_kind.get(r["kind"], 0) + 1
    print(f"wrote {out} - {len(rows)} site(s) over "
          f"{len({r['gene'] for r in rows})} gene(s), "
          f"{measured} measured and {len(rows) - measured} predicted")
    for kind, n in sorted(by_kind.items()):
        print(f"  {kind}: {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
