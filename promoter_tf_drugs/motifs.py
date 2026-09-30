"""Load TF binding motifs and scan a promoter with them.

Motif source: JASPAR CORE (non-redundant, vertebrates). Prefer the bulk PFM
file over the REST API: it is one download instead of ~900 requests, and it
pins you to a named release, which matters for reproducibility.

    https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt

Pass it with --pfm-file to run the whole motif stage offline.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import pwm
from .http import Client

log = logging.getLogger(__name__)

JASPAR_API = "https://jaspar.elixir.no/api/v1"
JASPAR_BULK = (
    "https://jaspar.elixir.no/download/data/2024/CORE/"
    "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
)


@dataclass
class Motif:
    matrix_id: str
    name: str                       # may be a complex, e.g. "FOS::JUN"
    counts: np.ndarray = field(repr=False)   # L x 4, columns A C G T

    @property
    def length(self) -> int:
        return self.counts.shape[0]

    @property
    def tf_symbols(self) -> list[str]:
        """Split complexes and strip variant suffixes -> gene symbols."""
        parts = re.split(r"::|[+/]", self.name)
        out = []
        for p in parts:
            p = p.strip().split("(")[0].strip().upper()
            if p:
                out.append(p)
        return out


@dataclass
class Hit:
    matrix_id: str
    motif_name: str
    tf_symbols: list[str]
    start: int                      # 0-based index into promoter sequence
    end: int
    strand: str
    score: float                    # log2-odds bits
    rel_score: float                # 0..1, (score-min)/(max-min)
    pvalue_threshold: float
    tss_offset: int                 # negative = upstream of TSS
    matched_sequence: str


# --------------------------------------------------------------------- load
_HEADER = re.compile(r"^>\s*(\S+)\s*(.*)$")
_ROW = re.compile(r"^\s*([ACGT])?\s*\[?\s*([0-9.eE+\-\s]+?)\s*\]?\s*$")


def parse_jaspar_pfm_text(text: str) -> list[Motif]:
    """Parse the classic JASPAR .jaspar / .pfm text format."""
    motifs: list[Motif] = []
    mid = name = None
    rows: dict[str, list[float]] = {}
    order: list[str] = []

    def flush():
        if mid is None:
            return
        if len(order) != 4:
            log.warning("skipping %s: expected 4 rows, got %d", mid, len(order))
            return
        keys = order if set(order) == set("ACGT") else list("ACGT")
        mat = np.array([rows[k] for k in ["A", "C", "G", "T"]] if set(order) == set("ACGT")
                       else [rows[k] for k in keys], dtype=float).T
        motifs.append(Motif(mid, name or mid, mat))

    for line in text.splitlines():
        if not line.strip():
            continue
        h = _HEADER.match(line)
        if h:
            flush()
            mid, name = h.group(1), (h.group(2).strip() or h.group(1))
            rows, order = {}, []
            continue
        m = _ROW.match(line)
        if m and mid is not None:
            base = m.group(1) or "ACGT"[len(order)]
            rows[base] = [float(x) for x in m.group(2).split()]
            order.append(base)
    flush()
    return motifs


def load_motifs(client: Client, pfm_file: str | Path | None = None,
                matrix_ids: list[str] | None = None) -> list[Motif]:
    """Load motifs from a local file, an explicit ID list, or the bulk URL."""
    if pfm_file:
        text = Path(pfm_file).read_text()
        motifs = parse_jaspar_pfm_text(text)
        log.info("loaded %d motifs from %s", len(motifs), pfm_file)
    elif matrix_ids:
        motifs = []
        for mid in matrix_ids:
            data = client.get(f"{JASPAR_API}/matrix/{mid}/", {"format": "json"})
            pfm = data["pfm"]
            counts = np.array([pfm[b] for b in "ACGT"], dtype=float).T
            motifs.append(Motif(data["matrix_id"], data["name"], counts))
    else:
        text = client.request(JASPAR_BULK, parse_json=False)
        motifs = parse_jaspar_pfm_text(text)
        log.info("loaded %d motifs from JASPAR bulk download", len(motifs))

    if not motifs:
        raise ValueError("no motifs loaded")
    return motifs


# --------------------------------------------------------------------- scan
def scan_promoter(sequence: str, motifs: list[Motif], *, pvalue: float = 1e-4,
                  min_ic: float = 6.0, offset_fn=None,
                  background: np.ndarray | None = None,
                  skip_repeats: bool = False) -> list[Hit]:
    """Scan `sequence` with every motif at a fixed per-position p-value.

    `min_ic` drops very low-information matrices (< ~6 bits). Those produce
    hits every few dozen bases in any sequence and only add noise.
    """
    enc = pwm.encode(sequence, skip_repeats=skip_repeats)
    bg = pwm.background_from_encoding(enc) if background is None else background
    hits: list[Hit] = []

    for motif in motifs:
        ic = pwm.information_content(motif.counts, bg)
        if ic < min_ic:
            continue
        pssm = pwm.pfm_to_pssm(motif.counts, bg)
        lo, hi = pwm.min_max_score(pssm)
        if hi <= lo:
            continue
        thr = pwm.score_threshold(pssm, bg, pvalue)

        for start, end, strand, score in pwm.scan_both_strands(sequence, pssm, thr,
                                                               enc=enc):
            sub = sequence[start:end]
            hits.append(Hit(
                matrix_id=motif.matrix_id,
                motif_name=motif.name,
                tf_symbols=motif.tf_symbols,
                start=start,
                end=end,
                strand=strand,
                score=score,
                rel_score=(score - lo) / (hi - lo),
                pvalue_threshold=pvalue,
                tss_offset=offset_fn(start) if offset_fn else start,
                matched_sequence=sub if strand == "+" else pwm.reverse_complement(sub),
            ))

    hits.sort(key=lambda h: (-h.score, h.start))
    return hits


def collapse_overlapping(hits: list[Hit], max_overlap: float = 0.75) -> list[Hit]:
    """Keep the best hit per motif among mutually overlapping windows.

    Palindromic and repeat-like motifs otherwise report the same site three or
    four times at consecutive offsets, which inflates every per-TF count.
    """
    by_matrix: dict[str, list[Hit]] = {}
    for h in hits:
        by_matrix.setdefault(h.matrix_id, []).append(h)

    kept: list[Hit] = []
    for group in by_matrix.values():
        group.sort(key=lambda h: -h.score)
        chosen: list[Hit] = []
        for h in group:
            length = h.end - h.start
            clash = False
            for c in chosen:
                ov = min(h.end, c.end) - max(h.start, c.start)
                if ov > 0 and ov / length > max_overlap:
                    clash = True
                    break
            if not clash:
                chosen.append(h)
        kept.extend(chosen)

    kept.sort(key=lambda h: (-h.score, h.start))
    return kept


def expected_false_positives(n_motifs: int, seq_length: int,
                             pvalue: float) -> float:
    """How many hits this scan produces on random sequence.

    The single most useful sanity check in the package, and the one nobody
    runs. Per-position p-value x length x 2 strands x number of matrices is
    the expected count under the null. If your observed count is not
    comfortably above this, the scan found nothing: you are looking at a list
    of chance matches, and ranking them changes nothing.
    """
    return n_motifs * 2.0 * max(seq_length, 1) * pvalue


def scannable_length(sequence: str, skip_repeats: bool = True) -> int:
    """Positions actually scanned: excludes N and, when masking is on, repeats.

    Using the window length instead overstates the null expectation whenever
    masking removes a chunk of the window — which produced the nonsense line
    "293 sites is what chance alone predicts (1055)".
    """
    enc = pwm.encode(sequence, skip_repeats=skip_repeats)
    return int((enc < pwm.N_CODE).sum())


def scan_quality(n_hits: int, n_motifs: int, seq_length: int,
                 pvalue: float) -> dict:
    """Compare observed hits against the null expectation.

    `seq_length` must be the scannable length (see `scannable_length`), not
    the window width.
    """
    expected = expected_false_positives(n_motifs, seq_length, pvalue)
    ratio = n_hits / expected if expected > 0 else float("inf")

    if ratio < 0.8:
        verdict = "no signal"
        advice = (
            f"{n_hits} sites is FEWER than chance alone predicts "
            f"({expected:.0f}) over the {seq_length} scannable bases. There is "
            f"no enrichment to interpret. If the gap is large, much of the "
            f"window was repeat-masked or ambiguous — check the sequence, and "
            f"tighten to --pvalue 1e-6 with a narrower window so any real "
            f"sites stand out.")
    elif ratio < 1.2:
        verdict = "no signal"
        advice = (
            f"{n_hits} sites is what chance alone predicts ({expected:.0f} over "
            f"{seq_length} scannable bases). This scan found nothing. Tighten "
            f"to --pvalue 1e-6, or narrow the window (--upstream 1000 "
            f"--downstream 200), or both — and use --chip-bed so ranking rests "
            f"on measured binding rather than sequence.")
    elif ratio < 2.0:
        verdict = "weak"
        advice = (f"{n_hits} sites against {expected:.0f} expected by chance "
                  f"over {seq_length} scannable bases. "
                  f"Most of this list is noise; treat only the highest-scoring "
                  f"sites with ChIP support as real.")
    else:
        verdict = "enriched"
        advice = (f"{n_hits} sites against {expected:.0f} expected by chance "
                  f"over {seq_length} scannable bases ({ratio:.1f}x).")

    return {"observed": n_hits, "expected_by_chance": round(expected, 1),
            "scannable_bases": seq_length,
            "enrichment": round(ratio, 2), "verdict": verdict, "advice": advice}
