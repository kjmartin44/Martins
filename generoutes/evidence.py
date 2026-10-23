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

Then pass --chip-bed data/remap2022.bed.gz.

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
import hashlib
import io
import json
import logging
import re
import tempfile
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


def _memo_path(bed_path: Path, chrom: str, start: int, end: int,
               name_to_tf: bool, memo_dir: Path) -> Path | None:
    """Where the answer for this exact file and region is kept, or None.

    Keyed on the file's identity (absolute path, size, modification time)
    and the region, so replacing or editing the BED invalidates it without
    anyone remembering to. Reading a 1.5 GB gzip for a 1.2 kb window takes
    about 72 seconds, and the launcher asks for the same window in three
    stages.
    """
    try:
        st = bed_path.stat()
    except OSError:
        return None
    ident = "|".join([str(bed_path.resolve()), str(st.st_size),
                      str(st.st_mtime_ns), chrom, str(start), str(end),
                      str(bool(name_to_tf))])
    key = hashlib.sha256(ident.encode("utf-8")).hexdigest()[:32]
    return Path(memo_dir) / "chip_regions" / f"{key}.json"


def peaks_in_region(bed_path: str | Path, chrom: str, start: int, end: int,
                    name_to_tf=True, memo_dir: str | Path | None = None
                    ) -> dict[str, int]:
    """Stream a BED and count overlapping intervals, keyed by TF symbol.

    Streaming keeps memory flat on a 30 GB ReMap file at the cost of one full
    pass per region. For many genes, pre-filter once with tabix or bedtools.

    `memo_dir`, when given, remembers the answer for this exact file and
    region so a second request for the same window does not re-read the
    file. It is off by default: a caller that does not name a directory
    gets no files written.

    The count is per TRANSCRIPTION FACTOR over the WHOLE window. It does
    not say where in the window a peak lies, so it cannot tell a motif hit
    on top of a peak from one 900 bp away; every hit of a factor with any
    peak in the window carries the same support. That is a coarser claim
    than "ChIP overlap" suggests, and the documentation says so.
    """
    chrom = _norm_chrom(chrom)
    memo = None
    if memo_dir is not None:
        memo = _memo_path(Path(bed_path), chrom, start, end, name_to_tf,
                          Path(memo_dir))
        if memo is not None and memo.is_file():
            try:
                return {k: int(v) for k, v in
                        json.loads(memo.read_text(encoding="utf-8")).items()}
            except (OSError, ValueError):
                pass                      # unreadable memo: just recompute
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
            # BED is 0-based half-open [s, e); the region is 1-based
            # inclusive [start, end], i.e. 0-based half-open
            # [start - 1, end). They overlap iff s < end and e > start - 1.
            # The upper bound used to be `s > end`, which kept a peak whose
            # first base is the one just past the window.
            if e < start or s >= end:
                continue
            if not name_to_tf or len(f) < 4:
                counts["_ANY_"] += 1
                continue
            for tf in _tf_from_name(f[3]):
                counts[tf] += 1
    result = dict(counts)
    if memo is not None:
        _write_memo(memo, result)
    return result


def _write_memo(memo: Path, result: dict) -> None:
    """Save one region's counts, through a temp file in the same directory.

    The temp name must not be LONGER than the name it replaces. It used to
    be `<32 hex>.<pid>.tmp` against a final `<32 hex>.json`, five characters
    more, so on Windows a memo path that fit under MAX_PATH could still fail
    to write: at 256 characters the final name was fine and the temp one was
    261. `tempfile` gives a short name, in the shape `http._write_cache`
    already uses.
    """
    tmp = None
    try:
        memo.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile("w", dir=memo.parent, suffix=".tmp",
                                         delete=False,
                                         encoding="utf-8") as fh:
            tmp = Path(fh.name)
            json.dump(result, fh)
        tmp.replace(memo)
        tmp = None
    except OSError as exc:
        # Not fatal: the counts this returns are correct either way. Not
        # silent either, which it was at debug level - the cost is that
        # every later request for this window re-reads the whole BED, about
        # 72 seconds on ReMap, and nobody could see why.
        log.warning("could not memoise the ChIP counts at %s (%s). The "
                    "counts are correct; the next request for this window "
                    "reads the BED again instead of reusing them.",
                    memo, exc)
    finally:
        if tmp is not None:
            try:
                tmp.unlink(missing_ok=True)
            except OSError:
                log.debug("could not remove the memo temporary file %s", tmp)


def region_is_open(bed_path: str | Path, chrom: str, start: int, end: int) -> bool:
    return bool(peaks_in_region(bed_path, chrom, start, end, name_to_tf=False))


def annotate_hits_with_chip(hits, promoter, chip_bed: str | Path | None):
    """Attach a `chip_peaks` count to each hit's TF. Returns {tf: count}.

    The counts were computed and returned and the hits were never touched,
    so the docstring described something that did not happen. A hit now
    carries `chip_peaks` for the best-supported of its TF symbols, and
    `chip_checked` saying a BED was actually read - absent is not zero,
    because a motif whose factor nobody assayed is unexamined rather than
    unsupported.
    """
    if not chip_bed:
        for h in hits or []:
            _set(h, "chip_checked", False)
        return {}

    counts = peaks_in_region(chip_bed, promoter.chrom, promoter.start,
                             promoter.end)
    for h in hits or []:
        syms = _get(h, "tf_symbols") or []
        best = max((counts.get(sym.upper(), 0) for sym in syms), default=0)
        _set(h, "chip_peaks", best)
        _set(h, "chip_checked", True)
    log.info("%d TFs with ChIP peaks over %s", len(counts), promoter.region)
    return counts


def _get(obj, field: str):
    """Read a field from a dataclass hit or a plain dict, whichever came.

    `getattr` on a dict returns None for every key, which is how the first
    version of this read no TF symbols at all from the dict hits the
    report passes - the same mistake three `uniprot` callers made.
    """
    if isinstance(obj, dict):
        return obj.get(field)
    return getattr(obj, field, None)


def _set(obj, field: str, value) -> None:
    """Set a field on a dataclass hit or a plain dict, whichever arrived."""
    if isinstance(obj, dict):
        obj[field] = value
        return
    try:
        setattr(obj, field, value)
    except (AttributeError, TypeError):
        # A frozen or slotted hit cannot carry the annotation. Say so once
        # rather than failing the run over a report field.
        log.debug("could not annotate %r with %s", type(obj).__name__, field)


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
        # The client's cache directory doubles as the memo location, so
        # the 72 s pass over a ReMap file for this window is paid once and
        # not again in every stage of a workup that asks for it.
        counts = peaks_in_region(chip_bed, promoter.chrom, promoter.start,
                                 promoter.end,
                                 memo_dir=getattr(client, "cache_dir", None))
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
