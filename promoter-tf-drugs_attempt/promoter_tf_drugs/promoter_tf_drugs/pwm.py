"""Position-weight-matrix mathematics.

Self-contained (numpy only). Implements:
  * 2-bit encoding of DNA with an N/ambiguity escape value
  * PFM -> PPM -> log-odds PSSM with a background-proportional pseudocount
  * an exact score distribution under a 0-order background, used to convert a
    desired per-position false-positive rate into a score threshold
  * both-strand scanning

The p-value machinery matters. Scanning at "80% of maximum score" is the most
common shortcut in the literature and it is not comparable across motifs:
an information-poor 6-mer and an information-rich 18-mer hit that cutoff at
wildly different false-positive rates. Thresholding on p-value puts every motif
on the same footing.
"""

from __future__ import annotations

import numpy as np

ALPHABET = "ACGT"
N_CODE = 4  # anything not ACGT

_LUT = np.full(256, N_CODE, dtype=np.int8)
for _i, _c in enumerate(ALPHABET):
    _LUT[ord(_c)] = _i
    _LUT[ord(_c.lower())] = _i

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


_LUT_HARD = np.full(256, N_CODE, dtype=np.int8)
for _i, _c in enumerate(ALPHABET):
    _LUT_HARD[ord(_c)] = _i          # lowercase deliberately left as N


def encode(seq: str, skip_repeats: bool = False) -> np.ndarray:
    """DNA string -> int8 array of 0..3, with 4 for N/ambiguous.

    With `skip_repeats=True`, soft-masked (lowercase) bases are treated as N.
    Ensembl returns repeat-soft-masked sequence when you request mask=soft;
    Alu and simple repeats are otherwise a major source of spurious hits,
    especially for GC-rich and homopolymer-like motifs.
    """
    table = _LUT_HARD if skip_repeats else _LUT
    return table[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def background_from_encoding(enc: np.ndarray, min_prob: float = 0.05) -> np.ndarray:
    """Base composition of an already-encoded sequence, ignoring N/masked."""
    counts = np.bincount(enc[enc < N_CODE], minlength=4).astype(float)
    if counts.sum() == 0:
        return np.full(4, 0.25)
    bg = np.clip(counts / counts.sum(), min_prob, None)
    return bg / bg.sum()


def background_from_sequence(seq: str, min_prob: float = 0.05) -> np.ndarray:
    """0-order base composition of the sequence itself.

    Human promoters are markedly GC-rich (CpG islands especially). Scoring a
    GC-rich motif such as SP1 against a flat 0.25 background inflates its score
    everywhere and produces a wall of false positives. Using the local
    composition is the cheapest correction available.
    """
    return background_from_encoding(encode(seq), min_prob)


def pfm_to_pssm(counts, bg: np.ndarray, pseudocount: float | None = None) -> np.ndarray:
    """Counts matrix (L x 4, columns A/C/G/T) -> log2-odds PSSM (L x 4)."""
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2 or counts.shape[1] != 4:
        raise ValueError(f"expected an L x 4 count matrix, got {counts.shape}")
    n_sites = counts.sum(axis=1, keepdims=True)
    if pseudocount is None:
        # Sqrt of the site count is the standard weakly-informative choice.
        pseudocount = float(np.sqrt(max(n_sites.mean(), 1.0)))
    probs = (counts + pseudocount * bg) / (n_sites + pseudocount)
    return np.log2(probs / bg)


def min_max_score(pssm: np.ndarray) -> tuple[float, float]:
    return float(pssm.min(axis=1).sum()), float(pssm.max(axis=1).sum())


def score_threshold(pssm: np.ndarray, bg: np.ndarray, pvalue: float,
                    scale: int = 100) -> float:
    """Lowest score whose right-tail probability under `bg` is <= `pvalue`.

    Exact dynamic programming over the discretised score distribution
    (the same construction MOODS and FIMO use).
    """
    ints = np.rint(np.asarray(pssm) * scale).astype(np.int64)
    row_min = ints.min(axis=1)
    offsets = ints - row_min[:, None]          # non-negative
    base = int(row_min.sum())

    dist = np.zeros(1, dtype=float)
    dist[0] = 1.0
    for row in offsets:
        new = np.zeros(len(dist) + int(row.max()), dtype=float)
        for b in range(4):
            o = int(row[b])
            new[o:o + len(dist)] += dist * bg[b]
        dist = new

    tail = np.cumsum(dist[::-1])[::-1]         # tail[s] = P(score >= s)
    hits = np.nonzero(tail <= pvalue)[0]
    offset = int(hits[0]) if hits.size else len(dist) - 1
    return (base + offset) / scale


def scan_strand(enc: np.ndarray, pssm: np.ndarray, threshold: float):
    """Score every window on one strand. Returns (start_indices, scores)."""
    L = pssm.shape[0]
    if len(enc) < L:
        return np.empty(0, dtype=int), np.empty(0, dtype=float)
    windows = np.lib.stride_tricks.sliding_window_view(enc, L)
    valid = (windows < N_CODE).all(axis=1)
    idx = np.clip(windows, 0, 3).astype(np.intp)
    scores = pssm[np.arange(L)[None, :], idx].sum(axis=1)
    keep = valid & (scores >= threshold)
    return np.nonzero(keep)[0], scores[keep]


def scan_both_strands(seq: str, pssm: np.ndarray, threshold: float,
                      enc: np.ndarray | None = None):
    """Yield (start, end, strand, score) for every window at or above threshold.

    Coordinates are 0-based half-open on the *given* (forward) sequence,
    regardless of which strand produced the hit.
    """
    if enc is None:
        enc = encode(seq)
    L = pssm.shape[0]

    starts, scores = scan_strand(enc, pssm, threshold)
    for s, sc in zip(starts.tolist(), scores.tolist()):
        yield int(s), int(s) + L, "+", float(sc)

    rc_pssm = pssm[::-1, ::-1]
    starts, scores = scan_strand(enc, rc_pssm, threshold)
    for s, sc in zip(starts.tolist(), scores.tolist()):
        yield int(s), int(s) + L, "-", float(sc)


def consensus(pssm_or_counts) -> str:
    m = np.asarray(pssm_or_counts, dtype=float)
    return "".join(ALPHABET[i] for i in m.argmax(axis=1))


def information_content(counts, bg: np.ndarray) -> float:
    """Total information content of a PFM in bits. Low-IC motifs are noisy."""
    counts = np.asarray(counts, dtype=float)
    n = counts.sum(axis=1, keepdims=True)
    p = (counts + 0.25) / (n + 1.0)
    return float(np.sum(p * np.log2(p / bg)))
