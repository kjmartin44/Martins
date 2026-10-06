"""Regulons built from ChIP peaks, for the factors curation never reached.

The gap this fills
------------------
A real run reported it plainly: 37 of the 200 factors with a motif hit in
the window had enough curated targets to test. The other 163 were not
absent from the result for any reason to do with biology. TRRUST holds
329 factors, 192 of them with ten or more targets, and JASPAR has motifs
for far more. A factor nobody has curated is untestable, and the report
correctly refused to read that as evidence about the factor.

ReMap is already on disk for the ChIP check, and it covers roughly a
thousand factors. Every peak names the factor that was pulled down, so
the set of genes whose promoter a factor binds can be read straight out
of it. That is a regulon, measured rather than curated, for factors the
curated tables never reach.

What a ChIP regulon is not
--------------------------
**Binding is not regulation.** A peak says the factor was at the promoter
in that experiment. Whether it drives the gene, represses it, or does
nothing there is a different question, and this file cannot answer it.
Curated regulons assert regulation; this asserts occupancy, and the two
should not be read as the same evidence.

**Unsigned.** There is no activation or repression direction in a peak
file, so every edge is written with mode `Unknown`. That is a shape
`load_signed_regulons` already handles - 46% of TRRUST edges are
unlabelled too - and it means the directional activity statistic cannot
use these edges, only the magnitude one. The honest consequence is that a
ChIP regulon tells you a factor's targets are moving, not which way.

**Pooled across cell types.** ReMap merges experiments from every cell
line assayed. A peak found only in K562 is in here with no marker saying
so, and your tissue may not have it.

**It does not fix the study-effort bias, it moves it.** A factor assayed
in two hundred experiments gets a large regulon and a factor assayed once
gets a small one, for the same reason TRRUST over-represents what people
have published on. The count reported per edge is how many experiments
support it, so a threshold is available, but no threshold turns an
unstudied factor into a studied one.
"""

from __future__ import annotations

import bisect
import gzip
import logging
from collections import defaultdict
from pathlib import Path

log = logging.getLogger(__name__)

# Promoter window a peak must fall in to count as targeting the gene.
# Asymmetric around the TSS because that is where promoter-proximal
# factors bind; a distal enhancer peak is a real regulatory event this
# deliberately does not capture, since assigning an enhancer to a gene
# needs contact data this package does not have.
UPSTREAM = 2000
DOWNSTREAM = 500

# Default output, picked up by `tfactivity.autoload_regulons` if it is
# added to DEFAULT_REGULON_PATHS, and usable with --regulon-file either
# way.
DEFAULT_OUT = "data/chip_regulon.tsv"


def _open(path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def _norm_chrom(c: str) -> str:
    c = c.strip()
    return c if c.startswith("chr") else "chr" + c


def promoter_index(gtf_path, *, upstream: int = UPSTREAM,
                   downstream: int = DOWNSTREAM,
                   protein_coding_only: bool = True) -> dict:
    """`{chrom: (starts, ends, genes)}` with each list sorted by start.

    One window per gene, around the gene's own 5' end. Sorted so a peak
    can be placed with a bisect rather than a scan: ReMap has of the order
    of a hundred million peaks, and a linear search per peak would not
    finish.

    Protein-coding only by default. A GENCODE GTF holds every biotype, so
    the first build of this indexed 77,114 genes - pseudogenes and lncRNAs
    among them - which inflated every regulon and pushed more factors over
    the promiscuity cap than it kept under it. Those genes are real and a
    factor may well bind them; they are simply not what a regulon in this
    package means, since nothing downstream can act on them.
    """
    spans: dict = defaultdict(list)
    seen: set = set()
    with _open(gtf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = f[8]
            if protein_coding_only and \
                    'gene_type "protein_coding"' not in attrs and \
                    'gene_biotype "protein_coding"' not in attrs:
                continue
            i = attrs.find('gene_name "')
            if i < 0:
                continue
            j = attrs.find('"', i + 11)
            name = attrs[i + 11:j].strip().upper()
            if not name or name in seen:
                continue
            seen.add(name)
            chrom = _norm_chrom(f[0])
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            if f[6] == "-":
                tss = end
                lo, hi = tss - downstream, tss + upstream
            else:
                tss = start
                lo, hi = tss - upstream, tss + downstream
            spans[chrom].append((max(0, lo), hi, name))

    out: dict = {}
    for chrom, rows in spans.items():
        rows.sort()
        out[chrom] = ([r[0] for r in rows], [r[1] for r in rows],
                      [r[2] for r in rows])
    log.info("promoter index: %d gene(s) across %d contig(s)",
             len(seen), len(out))
    return out


def _tf_cell(name: str) -> tuple:
    """`SMARCA4:HeLa` -> `("SMARCA4", "HELA")`.

    ReMap names the factor first and the cell line second, which is the
    only support signal in the file: how many distinct cell types show an
    edge. Counting peaks instead measures how finely the window was
    tiled, and came back as 1 for 92% of edges.
    """
    tf, _, cell = name.partition(":")
    return tf.strip().upper(), cell.strip().upper()


def build(bed_path, gtf_path, *, out_path=DEFAULT_OUT,
          upstream: int = UPSTREAM, downstream: int = DOWNSTREAM,
          min_cell_types: int = 1, tfs: set | None = None,
          max_targets: int | None = None,
          protein_coding_only: bool = True) -> dict:
    """Stream a ReMap BED into a regulon table, one pass.

    `min_cell_types` is how many distinct cell lines must show a
    factor-gene edge. One is permissive: raising it trades coverage, which
    is the thing this exists to add, for an edge that reproduced in a
    second cellular context.

    `max_targets` drops factors that bind more promoters than that. A
    factor with twenty thousand targets is not specific to anything and
    carrying it costs the activity statistic its power, so it is excluded
    and reported rather than silently kept.
    """
    index = promoter_index(gtf_path, upstream=upstream,
                           downstream=downstream,
                           protein_coding_only=protein_coding_only)
    want = {t.strip().upper() for t in tfs} if tfs else None

    # Each edge carries a bitmask of the cell types supporting it, so
    # distinct-cell-type support costs one integer per edge rather than a
    # set of strings per edge - there are 16 million placements and the
    # sets would not fit.
    cell_bit: dict = {}
    edges: dict = defaultdict(lambda: defaultdict(int))
    n_lines = n_placed = 0

    with _open(bed_path) as fh:
        for line in fh:
            if not line or line[0] in "#t":
                continue
            f = line.split("\t", 4)
            if len(f) < 4:
                continue
            n_lines += 1
            chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
            idx = index.get(chrom)
            if idx is None:
                continue
            tf, cell = _tf_cell(f[3])
            if not tf or (want is not None and tf not in want):
                continue
            bit = cell_bit.get(cell)
            if bit is None:
                bit = cell_bit[cell] = 1 << len(cell_bit)
            try:
                p_start, p_end = int(f[1]), int(f[2])
            except ValueError:
                continue

            starts, ends, genes = idx
            # Windows are sorted by start. Every window that could overlap
            # begins at or before p_end, and the ones that do also end at
            # or after p_start. Walking back from the insertion point is
            # bounded because promoter windows are short and few overlap.
            k = bisect.bisect_right(starts, p_end)
            hit = False
            for m in range(k - 1, max(-1, k - 60), -1):
                if ends[m] >= p_start:
                    edges[tf][genes[m]] |= bit
                    hit = True
            if hit:
                n_placed += 1

    log.info("%d distinct cell type(s) in the peak file", len(cell_bit))
    return _write(edges, out_path, n_lines=n_lines, n_placed=n_placed,
                  min_cell_types=min_cell_types, max_targets=max_targets,
                  upstream=upstream, downstream=downstream,
                  n_cell_types=len(cell_bit))


def _write(edges: dict, out_path, *, n_lines: int, n_placed: int,
           min_cell_types: int, max_targets: int | None,
           upstream: int, downstream: int, n_cell_types: int = 0) -> dict:
    """Write the table and report what was kept, dropped and why."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    too_broad: list = []
    n_edges = 0
    kept_tfs = 0
    sizes: list = []
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        # TRRUST's layout, which load_signed_regulons already reads. The
        # mode is Unknown for every edge because a peak carries no
        # direction, and the last column holds the experiment count so a
        # reader can threshold on support rather than take it on trust.
        for tf in sorted(edges):
            targets = {g: mask.bit_count() for g, mask in edges[tf].items()
                       if mask.bit_count() >= min_cell_types}
            if not targets:
                continue
            if max_targets is not None and len(targets) > max_targets:
                too_broad.append((tf, len(targets)))
                continue
            kept_tfs += 1
            sizes.append(len(targets))
            for gene in sorted(targets):
                # The last column is the number of distinct cell types
                # supporting the edge, so a reader can threshold on
                # reproducibility rather than take the edge on trust.
                fh.write(f"{tf}\t{gene}\tUnknown\tReMap_cells:"
                         f"{targets[gene]}\n")
                n_edges += 1

    acc = {
        "path": str(out_path),
        "n_factors": kept_tfs,
        "n_edges": n_edges,
        "n_peaks_read": n_lines,
        "n_peaks_in_a_promoter": n_placed,
        "median_targets": (sorted(sizes)[len(sizes) // 2] if sizes else 0),
        "excluded_too_broad": too_broad[:20],
        "n_excluded_too_broad": len(too_broad),
        "min_cell_types": min_cell_types,
        "n_cell_types_in_file": n_cell_types,
        "window": f"-{upstream}/+{downstream} around the TSS",
    }
    log.info("ChIP regulon: %d factor(s), %d edge(s) from %d peak(s)",
             kept_tfs, n_edges, n_lines)
    return acc


def summarise(acc: dict) -> str:
    """Prose for whoever runs this, including what it does not establish."""
    if not acc:
        return "no ChIP regulon was built."
    parts = [
        f"{acc['n_factors']} factor(s) and {acc['n_edges']} edge(s) from "
        f"{acc['n_peaks_read']} peak(s), of which "
        f"{acc['n_peaks_in_a_promoter']} fell in a promoter window "
        f"({acc['window']}). Median regulon size "
        f"{acc['median_targets']}.",
        "Every edge is occupancy, not regulation: the factor was at that "
        f"promoter in at least {acc['min_cell_types']} cell type(s), which "
        "is not the same as the factor driving the gene. Every edge is "
        "also unsigned, so the directional activity statistic cannot use "
        "these and the magnitude one can.",
        "ReMap pools cell lines. A peak seen only in one is in here "
        "unmarked, and your tissue may not have it.",
    ]
    if acc["n_excluded_too_broad"]:
        names = ", ".join(f"{t} ({n})"
                          for t, n in acc["excluded_too_broad"][:5])
        parts.append(
            f"{acc['n_excluded_too_broad']} factor(s) bound more promoters "
            f"than the cap and were excluded ({names}). A factor binding "
            f"most of the genome is not specific to your gene, and keeping "
            f"it costs the activity statistic its power.")
    parts.append(
        "This does not fix the study-effort bias in regulon size, it moves "
        "it: a factor assayed two hundred times gets a large regulon and "
        "one assayed once gets a small one, the same way curation "
        "over-represents what people have published on.")
    return " ".join(parts)


def main(argv=None) -> int:
    """`python -m promoter_tf_drugs.chipregulon --bed ... --gtf ...`"""
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.chipregulon",
        description="Build a regulon table from ChIP peaks, for the "
                    "factors curated tables never reached.")
    p.add_argument("--bed", required=True,
                   help="ReMap peak BED, optionally gzipped")
    p.add_argument("--gtf", required=True,
                   help="GENCODE/Ensembl GTF for promoter coordinates")
    p.add_argument("--out", default=DEFAULT_OUT)
    p.add_argument("--upstream", type=int, default=UPSTREAM)
    p.add_argument("--downstream", type=int, default=DOWNSTREAM)
    p.add_argument("--min-cell-types", type=int, default=1, metavar="N",
                   help="distinct cell lines that must show a factor-gene "
                        "edge. Raising it trades the coverage this exists "
                        "to add for an edge that reproduced in a second "
                        "cellular context (default 1)")
    p.add_argument("--all-biotypes", action="store_true",
                   help="index every GENCODE gene rather than protein "
                        "coding only. Inflates every regulon with "
                        "pseudogenes and lncRNAs that nothing downstream "
                        "can act on")
    p.add_argument("--max-targets", type=int, default=5000, metavar="N",
                   help="drop factors binding more promoters than this; "
                        "they are not specific to anything (default 5000)")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s")

    acc = build(args.bed, args.gtf, out_path=args.out,
                upstream=args.upstream, downstream=args.downstream,
                min_cell_types=args.min_cell_types,
                max_targets=args.max_targets,
                protein_coding_only=not args.all_biotypes)
    print(summarise(acc))
    print(f"\nWrote {acc['path']}. Use it with --regulon-file.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
