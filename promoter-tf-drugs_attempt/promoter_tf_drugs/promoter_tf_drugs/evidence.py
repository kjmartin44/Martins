"""Orthogonal evidence that a predicted site is real.

This is the single highest-value module in the package, and it is the one that
needs a file from you. A PWM hit on its own is close to worthless as evidence
of binding: a typical 8-12 bp motif matches by chance roughly every few
kilobases per strand, so a 1.2 kb promoter will "contain" several hundred
motifs for several hundred TFs in any species, including in shuffled sequence.

Intersecting with measured ChIP-seq peaks changes the claim from "this TF
*could* bind here" to "this TF *was observed* bound here in some cell type".

Get a peak catalogue (BED, one-time download, no API):

  ReMap 2022 (all human TF ChIP-seq, non-redundant peaks)
      https://remap.univ-amu.fr/  -> remap2022_nr_macs2_hg38_v1_0.bed.gz
  ChIP-Atlas (per-TF or all-TF BEDs, multiple assemblies)
      https://chip-atlas.org/

Then pass --chip-bed remap2022_nr_macs2_hg38_v1_0.bed.gz.

Or skip the download entirely: ReMap exposes a REST API that returns peaks for
one region, which is far lighter than scanning a multi-GB file per gene. Use
--chip-api instead of --chip-bed. The local file is still worth having if you
are analysing many genes or want reproducibility pinned to a release.

Open chromatin (optional, --dhs-bed): ENCODE cCREs or a DNase/ATAC peak set for
the tissue you actually care about. A site in closed chromatin in your cell
type is not available for binding no matter how good the motif is.
"""

from __future__ import annotations

import gzip
import io
import logging
import re
from collections import defaultdict
from pathlib import Path

log = logging.getLogger(__name__)


def _open(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _norm_chrom(c: str) -> str:
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def _tf_from_name(name: str) -> list[str]:
    """ReMap names look like 'GSE12345.CTCF.HeLa-S3'; ChIP-Atlas varies."""
    parts = re.split(r"[.:;|]", name)
    out = []
    for p in parts:
        p = p.strip().upper()
        if not p or p.startswith("GSE") or p.startswith("SRX") or p.startswith("ERX"):
            continue
        if re.fullmatch(r"[A-Z0-9\-]{2,15}", p):
            out.append(p)
    return out[:1] or ([parts[0].strip().upper()] if parts else [])


def peaks_in_region(bed_path: str | Path, chrom: str, start: int, end: int,
                    name_to_tf=True) -> dict[str, int]:
    """Stream a BED and count overlapping intervals, keyed by TF symbol.

    Streaming keeps memory flat on a 30 GB ReMap file at the cost of one full
    pass per region. For many genes, pre-filter once with tabix or bedtools.
    """
    chrom = _norm_chrom(chrom)
    counts: dict[str, int] = defaultdict(int)
    with _open(bed_path) as fh:
        for line in fh:
            if not line or line[0] in "#t":  # skip comments and 'track' lines
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            if _norm_chrom(f[0]) != chrom:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            if e < start or s > end:          # BED is 0-based half-open; region 1-based
                continue
            if not name_to_tf or len(f) < 4:
                counts["_ANY_"] += 1
                continue
            for tf in _tf_from_name(f[3]):
                counts[tf] += 1
    return dict(counts)


def region_is_open(bed_path: str | Path, chrom: str, start: int, end: int) -> bool:
    return bool(peaks_in_region(bed_path, chrom, start, end, name_to_tf=False))


def annotate_hits_with_chip(hits, promoter, chip_bed: str | Path | None):
    """Attach a `chip_peaks` count to each hit's TF. Returns {tf: count}."""
    if not chip_bed:
        return {}
    counts = peaks_in_region(chip_bed, promoter.chrom, promoter.start, promoter.end)
    log.info("%d TFs with ChIP peaks over %s", len(counts), promoter.region)
    return counts


# ----------------------------------------------------- ReMap REST API route
# Documented at https://remap-rest.univ-amu.fr/docs_page
#   /get_peaks/{release}/{assembly}/{datatype}/{chr:start-end}
# datatype is one of all | nr | crm. "nr" (non-redundant) matches the BED file
# recommended above and returns one merged region per target, which is what a
# presence/absence check wants.
REMAP_API = "https://remap-rest.univ-amu.fr/api/V1/get_peaks"

# The documented response nests the target three levels down:
#   {"peaks": [{"peakValues": {"name": {"data": {"TF": "CTCF", ...}}}}]}
# Reading row["target"] at the top level, as this module used to, finds
# nothing — which is why a whole promoter came back with one TF.
_TF_PATHS = (
    ("peakValues", "name", "data", "TF"),
    ("peakValues", "name", "data", "target"),
    ("name", "data", "TF"),
)


def _dig(row, path):
    """Follow a nested key path, returning None if any step is missing."""
    cur = row
    for key in path:
        if not isinstance(cur, dict):
            return None
        cur = cur.get(key)
    return cur if isinstance(cur, str) and cur.strip() else None


def _target_of(row) -> str | None:
    """Extract the TF name from one peak record, whatever shape it arrives in."""
    if isinstance(row, str):
        return row.strip() or None
    if not isinstance(row, dict):
        return None
    for path in _TF_PATHS:
        hit = _dig(row, path)
        if hit:
            return hit
    # Flat fallbacks, for other releases or a changed schema.
    for key in ("TF", "target", "target_name", "tf", "tf_name"):
        val = row.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return None


def peaks_from_api(client, chrom: str, start: int, end: int, *,
                   assembly: str = "hg38", version: str = "2022",
                   datatype: str = "nr") -> dict[str, int]:
    """Peak counts per TF over one region, straight from the ReMap API.

    Much faster than scanning the full catalogue: one small request instead of
    a pass over several GB. Requires network access; falls back to an empty
    result (never a crash) so a failed lookup degrades to "no ChIP evidence"
    rather than taking the run down.
    """
    c = chrom if str(chrom).lower().startswith("chr") else f"chr{chrom}"
    url = f"{REMAP_API}/{version}/{assembly}/{datatype}/{c}:{start}-{end}"
    try:
        res = client.get(url)
    except Exception as exc:                            # noqa: BLE001
        log.warning("ReMap API lookup failed for %s: %s", url, exc)
        return {}

    rows = _extract_rows(res)
    counts: dict[str, int] = defaultdict(int)
    unparsed = 0
    for row in rows:
        tf = _target_of(row)
        if not tf:
            unparsed += 1
            continue
        counts[tf.strip().upper()] += 1

    if unparsed and not counts:
        log.warning("ReMap returned %d peak record(s) but the target field "
                    "could not be located in any of them — the response schema "
                    "has probably changed. Use --chip-bed instead.", unparsed)

    region = url.split("/")[-1]
    log.info("ReMap API: %d TFs with peaks over %s", len(counts), region)

    # A promoter-sized region normally overlaps peaks for dozens of factors.
    # Two or fewer is far more likely a parsing or pagination problem than
    # biology, and silently treating every TF as "no measured binding" skews
    # the whole ranking. Say so rather than letting it pass.
    if len(counts) <= 2:
        log.warning(
            "Only %d TF(s) returned for %s. A promoter-sized region usually "
            "overlaps peaks for dozens of factors, so this is more likely an "
            "API or parsing problem than a real absence — ChIP evidence is "
            "effectively OFF for this run. Prefer --chip-bed with a local "
            "ReMap BED, which is parsed directly and not subject to this.",
            len(counts), region)
        log.debug("unparsed ReMap response shape: %s",
                  type(res).__name__ if not isinstance(res, dict)
                  else list(res)[:8])
    return dict(counts)


def _extract_rows(res):
    """Pull the peak records out of whatever shape the API returned.

    The response has changed shape between ReMap versions, and a wrong guess
    silently yields an empty list rather than an error.
    """
    if isinstance(res, list):
        return res
    if not isinstance(res, dict):
        return []
    for key in ("peaks", "results", "data", "features", "items", "hits"):
        val = res.get(key)
        if isinstance(val, list):
            return val
        if isinstance(val, dict):
            # e.g. {"data": {"chr15:1-2": [ ... ]}}
            out = []
            for v in val.values():
                if isinstance(v, list):
                    out.extend(v)
            if out:
                return out
    # Last resort: a dict keyed by region or by target name.
    out = []
    for v in res.values():
        if isinstance(v, list):
            out.extend(v)
    return out


# Below this many distinct factors over a promoter-sized region, the result is
# far likelier to be a parsing or pagination failure than a real absence of
# binding. `peaks_from_api` already warns at this threshold; here it decides
# whether the result may be used as scoring input at all.
MIN_CREDIBLE_TFS = 3


def chip_evidence(promoter, *, chip_bed=None, client=None, use_api: bool = False,
                  assembly: str = "hg38") -> dict[str, int]:
    """Unified entry point: local BED if given, else the API if requested."""
    counts, _ = chip_evidence_status(promoter, chip_bed=chip_bed, client=client,
                                     use_api=use_api, assembly=assembly)
    return counts


def chip_evidence_status(promoter, *, chip_bed=None, client=None,
                         use_api: bool = False, assembly: str = "hg38"):
    """As `chip_evidence`, but also report whether the source actually worked.

    The counts alone cannot answer that. An empty dict means "this promoter has
    no measured binding" when the BED parsed fine, and "the lookup broke" when
    the API changed shape — opposite conclusions with identical representation.
    Scoring must not treat the second as measured absence, so the caller is
    handed an explicit verdict instead of having to guess from a count.

    Returns (counts, status) where status has: requested, ok, detail, yielded.
    """
    if chip_bed:
        counts = peaks_in_region(chip_bed, promoter.chrom, promoter.start,
                                 promoter.end)
        # A local BED is parsed directly, so an empty result here is a real
        # statement about this region rather than a transport failure.
        return counts, {
            "requested": True, "ok": True, "yielded": len(counts),
            "detail": (f"local BED: peaks for {len(counts)} factor(s) over "
                       f"the window" if counts else
                       "local BED parsed, no peaks over this window — a real "
                       "absence of measured binding here"),
        }

    if use_api and client is not None:
        counts = peaks_from_api(client, promoter.chrom, promoter.start,
                                promoter.end, assembly=assembly)
        ok = len(counts) >= MIN_CREDIBLE_TFS
        return counts, {
            "requested": True, "ok": ok, "yielded": len(counts),
            "detail": (
                f"ReMap API: peaks for {len(counts)} factor(s)" if ok else
                f"ReMap API returned only {len(counts)} factor(s) for a "
                f"promoter-sized region, which normally overlaps peaks for "
                f"dozens. Treating this as a failed lookup rather than as "
                f"measured absence. Use --chip-bed with a local ReMap BED."),
        }

    return {}, {"requested": False, "ok": False, "yielded": 0,
                "detail": "no ChIP source given"}
