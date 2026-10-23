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
    # P(a random window of this motif scores at least this well), under the
    # window's own base composition. Every hit passes the same threshold, so
    # this is what ranks them against hits of other motifs.
    pvalue: float | None = None


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
        text = Path(pfm_file).read_text(encoding="utf-8")
        motifs = parse_jaspar_pfm_text(text)
        if matrix_ids:
            wanted = set(matrix_ids)
            motifs = [motif for motif in motifs if motif.matrix_id in wanted]
            missing = wanted - {motif.matrix_id for motif in motifs}
            if missing:
                raise ValueError("matrix IDs not found in PFM file: " +
                                 ", ".join(sorted(missing)))
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
                  skip_repeats: bool = False,
                  info: dict | None = None) -> list[Hit]:
    """Scan `sequence` with every motif at a fixed per-position p-value.

    `min_ic` drops very low-information matrices (< ~6 bits). Those produce
    hits every few dozen bases in any sequence and only add noise.

    Each strand has its own threshold, and a strand that cannot reach
    `pvalue` in this window (its best possible match is more common than
    that) is not scanned. `info`, when given, is filled with how many
    matrices were scanned, how many were below `min_ic`, and how many could
    not reach the p-value on either strand.
    """
    counts = {"n_scanned": 0, "n_low_information": 0, "n_unreachable": 0}
    enc = pwm.encode(sequence, skip_repeats=skip_repeats)
    bg = pwm.background_from_encoding(enc) if background is None else background
    hits: list[Hit] = []

    for motif in motifs:
        ic = pwm.information_content(motif.counts, bg)
        if ic < min_ic:
            counts["n_low_information"] += 1
            continue
        pssm = pwm.pfm_to_pssm(motif.counts, bg)
        lo, hi = pwm.min_max_score(pssm)
        if hi <= lo:
            continue
        base, tail = pwm.score_tail(pssm, bg)
        base_r, tail_r = pwm.score_tail(pssm[::-1, ::-1], bg)
        thr = pwm.threshold_from_tail(base, tail, pvalue, strict=True)
        thr_r = pwm.threshold_from_tail(base_r, tail_r, pvalue, strict=True)
        if thr is None and thr_r is None:
            counts["n_unreachable"] += 1
            continue
        counts["n_scanned"] += 1

        for start, end, strand, score in pwm.scan_strands(enc, pssm, thr, thr_r):
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
                pvalue=(pwm.tail_pvalue(score, base, tail) if strand == "+"
                        else pwm.tail_pvalue(score, base_r, tail_r)),
            ))

    hits.sort(key=lambda h: (-h.score, h.start))
    if info is not None:
        info.update(counts)
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


def dinucleotide_shuffle(seq: str, rng) -> str:
    """A random sequence with the same dinucleotide counts as `seq`.

    The Altschul and Erikson (1985) construction: every shuffle is an Eulerian
    path through the graph of observed A-to-C, C-to-G and so on steps, so the
    first and last bases and every two-base count are kept. A mononucleotide
    shuffle would also be a null, but it removes CpG depletion and runs of one
    base, which real promoters have and which change how many motifs match.
    """
    if len(seq) < 4:
        return seq
    edges: dict = {}
    for a, b in zip(seq, seq[1:]):
        edges.setdefault(a, []).append(b)
    last = seq[-1]
    nodes = [a for a in edges if a != last]
    # Pick each node's final exit so that the exits form a tree leading to the
    # last base; any other choice leaves the walk stranded before the end.
    for _ in range(1000):
        final = {a: rng.choice(edges[a]) for a in nodes}
        good = True
        for a in nodes:
            seen, cur = set(), a
            while cur != last:
                if cur in seen or cur not in final:
                    good = False
                    break
                seen.add(cur)
                cur = final[cur]
            if not good:
                break
        if good:
            break
    else:
        return seq
    pools = {}
    for a, outs in edges.items():
        rest = list(outs)
        if a in final:
            rest.remove(final[a])
            rng.shuffle(rest)
            rest.append(final[a])
        else:
            rng.shuffle(rest)
        pools[a] = rest
    used = {a: 0 for a in pools}
    out = [seq[0]]
    cur = seq[0]
    for _ in range(len(seq) - 1):
        nxt = pools[cur][used[cur]]
        used[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def shuffle_scannable(sequence: str, rng, *, skip_repeats: bool = True) -> str:
    """Shuffle each run of scanned bases, leaving the mask and Ns where they are.

    With repeat masking on, only upper-case runs are scanned, so only those
    are shuffled and every masked base stays at its position. Shuffling the
    whole sequence instead moves repeat bases into the scanned region and
    changes the composition the scan is judged against.
    """
    pattern = r"[ACGT]+" if skip_repeats else r"[ACGTacgt]+"
    out, last = [], 0
    for m in re.finditer(pattern, sequence):
        out.append(sequence[last:m.start()])
        run = m.group(0) if skip_repeats else m.group(0).upper()
        out.append(dinucleotide_shuffle(run, rng))
        last = m.end()
    out.append(sequence[last:])
    return "".join(out)


NULL_SHUFFLES = 30


def empirical_null(sequence: str, motifs: list, *, pvalue: float,
                   min_ic: float = 6.0, skip_repeats: bool = True,
                   n_shuffles: int = NULL_SHUFFLES, seed: int = 0) -> dict:
    """Run the same scan and the same collapsing on shuffled copies.

    The count and the scores of the sites each copy gives are the null this
    scan is judged against. Every setting is the one the real scan used, so
    nothing but the order of the bases differs.
    """
    import random
    rng = random.Random(seed)
    counts, scores = [], []
    for _ in range(max(0, n_shuffles)):
        shuffled = shuffle_scannable(sequence, rng, skip_repeats=skip_repeats)
        hits = collapse_overlapping(scan_promoter(
            shuffled, motifs, pvalue=pvalue, min_ic=min_ic,
            skip_repeats=skip_repeats))
        counts.append(len(hits))
        scores.append(sorted(h.pvalue for h in hits))
    return {"counts": counts, "scores": scores, "n_shuffles": len(counts),
            "seed": seed}


def _quantile(values: list, q: float) -> float:
    vs = sorted(values)
    if not vs:
        return 0.0
    pos = q * (len(vs) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(vs) - 1)
    return vs[lo] + (vs[hi] - vs[lo]) * (pos - lo)


def occupancy_partition(factors, quality: dict) -> dict:
    """Divide the scan's expected chance hits over the unmeasured sites.

    `factors` are the per-factor records the run built, each with a
    `chip_peaks` count. A factor with a measured peak at this window is
    one a shuffle of the window could not have produced, so the chance
    hits the null predicts are attributed to the factors without one.

    What this does and does not license
    -----------------------------------
    A peak lowers the risk that the SITE is an artefact of scanning. It
    does not show that THIS factor is what bound there - peaks are wide,
    and a window usually carries several motifs - and the experiment was
    run in whatever cell type the atlas used, not yours. So `supported`
    is "something was measured binding here", and nothing stronger.
    """
    rows = list(factors or [])
    if not rows or not quality:
        return {}
    supported = [f for f in rows if (f.get("chip_peaks") or 0) > 0]
    bare = [f for f in rows if not (f.get("chip_peaks") or 0)]
    n = len(rows)
    expected = quality.get("expected_by_chance")
    if expected is None or not n:
        return {}
    out = {
        "n_factors": n,
        "n_supported": len(supported),
        "n_unsupported": len(bare),
        "supported": sorted(f.get("symbol") or "?" for f in supported),
        "unsupported": sorted(f.get("symbol") or "?" for f in bare),
        "expected_by_chance": expected,
    }
    # Every chance hit is assumed to land among the unsupported sites.
    # That is the pessimistic reading for them and the generous one for
    # the measured sites, and it is the right way round: the cost of
    # over-trusting an unmeasured site is a wasted experiment.
    out["fdp_unsupported"] = (round(min(1.0, expected / len(bare)), 3)
                              if bare else 0.0)
    # Capped, because this is a proportion. The surplus spilling past
    # the unmeasured sites is real, but once it fills the measured ones
    # too the scan has no clean subset left and the number stops
    # carrying information; `saturated` below is what to read instead.
    out["fdp_supported"] = round(
        min(1.0, max(0.0, (expected - min(expected, len(bare)))
                     / len(supported))), 3
    ) if supported else None
    # True when chance alone accounts for every site, measured or not.
    out["saturated"] = bool(expected >= n)
    out["drugged_unsupported"] = sorted(
        f.get("symbol") or "?" for f in bare if f.get("drugs"))
    out["drugged_supported"] = sorted(
        f.get("symbol") or "?" for f in supported if f.get("drugs"))
    out["note"] = _occupancy_note(out)
    return out


def _occupancy_note(o: dict) -> str:
    n, s, b = o["n_factors"], o["n_supported"], o["n_unsupported"]
    if not b:
        return (f"All {n} factors with a site here also have a measured "
                f"peak over this window, so none of them rests on the scan "
                f"alone.")
    if not s:
        return (f"None of the {n} factors with a site here has a measured "
                f"peak over this window, so every one rests on the scan "
                f"alone and the window-wide false-discovery proportion "
                f"applies to all of them unchanged.")
    if o.get("saturated"):
        return (
            f"This scan is too loose for the split to help. Chance alone "
            f"accounts for {o['expected_by_chance']} sites against {n} "
            f"found, so there is no subset of it - measured or not - that "
            f"chance does not already explain, and {s} measured peaks "
            f"cannot rescue it. A stricter p-value is the fix; nothing "
            f"that can be done to this hit list is.")
    parts = [
        f"{s} of {n} factors have a measured peak over this window and {b} "
        f"do not. A shuffled copy of the window cannot produce a peak, so "
        f"the {o['expected_by_chance']} sites the null expects by chance "
        f"are attributed to those {b}: about "
        f"{o['fdp_unsupported']:.0%} of them are consistent with chance, "
        f"against roughly {o['fdp_supported']:.0%} of the measured ones."]
    if o["drugged_unsupported"]:
        parts.append(
            f"This matters for the ranking because "
            f"{', '.join(o['drugged_unsupported'])} "
            f"{'have' if len(o['drugged_unsupported']) > 1 else 'has'} a "
            f"ligand and no measured occupancy: every drug reaching the "
            f"gene through "
            f"{'them' if len(o['drugged_unsupported']) > 1 else 'it'} rests "
            f"on a predicted site.")
    if o["drugged_supported"]:
        parts.append(
            f"{', '.join(o['drugged_supported'])} "
            f"{'have' if len(o['drugged_supported']) > 1 else 'has'} both a "
            f"ligand and a measured peak, which is the stronger footing.")
    parts.append(
        "A peak says something bound here in the cell type the atlas used. "
        "It does not say this factor is what bound, and it is not a "
        "measurement in your cells.")
    return " ".join(parts)


def scan_quality(n_hits: int, n_motifs: int, seq_length: int,
                 pvalue: float, null: dict | None = None,
                 observed_scores: list | None = None) -> dict:
    """Compare observed hits against the null expectation.

    `seq_length` must be the scannable length (see `scannable_length`), not
    the window width.

    Without `null` the expectation is the formula `expected_false_positives`,
    which assumes every motif matches at exactly its stated p-value at every
    position. On a real repeat-masked promoter that undercounted: 13 expected
    against 24 to 27 (plus or minus 10) on shuffled copies, because score
    thresholds on a discrete distribution give more matches than the nominal
    p-value, and the formula ignores the mask. With `null` (see
    `empirical_null`) the expectation is the mean over shuffled copies, with
    a range from their spread, and the formula is kept beside it.
    """
    formula = expected_false_positives(n_motifs, seq_length, pvalue)
    counts = (null or {}).get("counts") or []
    empirical = len(counts) >= 10
    expected = (sum(counts) / len(counts)) if empirical else formula
    # No sites is no signal, even when chance gives none either; dividing 0 by
    # 0 would otherwise read as infinite enrichment.
    ratio = (n_hits / expected if expected > 0
             else float("inf") if n_hits > 0 else 0.0)

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

    # The share of the reported list expected to be chance. This is the
    # quantity the ranking actually needs, and it is not the same claim as the
    # enrichment ratio: 1.21x enrichment sounds like mild signal, but it means
    # roughly 83% of the rows below are noise. Stating it as a proportion makes
    # that unmissable, and it is a real Benjamini-style false discovery
    # PROPORTION — expected false positives over observed — not a fitted score.
    fdp = min(1.0, expected / n_hits) if n_hits > 0 else 0.0
    real = max(0, round(n_hits - expected))

    out = {"observed": n_hits, "expected_by_chance": round(expected, 1),
           "expected_by_formula": round(formula, 1),
           "scannable_bases": seq_length,
           "enrichment": round(ratio, 2), "verdict": verdict, "advice": advice,
           "false_discovery_proportion": round(fdp, 3),
           "expected_real_sites": real,
           "null": None}
    if not empirical:
        out["fdp_note"] = (
            f"about {fdp:.0%} of these {n_hits} sites are expected to be "
            f"chance matches, from a formula that assumes every motif "
            f"matches at exactly its stated p-value; roughly {real} are "
            f"real, but nothing here says which. Only ChIP evidence "
            f"separates them.")
        return out

    lo = _quantile(counts, 0.10)
    hi = _quantile(counts, 0.90)
    fdp_lo = min(1.0, lo / n_hits) if n_hits > 0 else 0.0
    fdp_hi = min(1.0, hi / n_hits) if n_hits > 0 else 0.0
    mean = sum(counts) / len(counts)
    sd = (sum((c - mean) ** 2 for c in counts) / len(counts)) ** 0.5
    out["null"] = {
        "method": ("shuffled copies of this window with the repeat mask and "
                   "every dinucleotide count kept, scanned and collapsed the "
                   "same way"),
        "n_shuffles": len(counts), "mean": round(mean, 1), "sd": round(sd, 1),
        "p10": round(lo, 1), "p90": round(hi, 1),
    }
    out["false_discovery_range"] = [round(fdp_lo, 3), round(fdp_hi, 3)]

    # How the proportion falls as the list is cut to its best sites. Observed
    # and null sites are ranked the same way, by the p-value of their score,
    # so the share of chance in the best k is the null's mean count at or
    # below the k-th p-value, over k.
    top = []
    obs = sorted(observed_scores or [])
    for k in (5, 10, 20):
        if k > len(obs) or k >= n_hits:
            continue
        cut = obs[k - 1]
        above = [sum(1 for s in sc if s <= cut) for sc in null["scores"]]
        top.append({"top": k, "fdp": round(min(1.0, sum(above) / len(above)
                                                / k), 3)})
    if top:
        out["top_site_fdp"] = top

    parts = [f"about {fdp:.0%} of these {n_hits} sites "
             f"(range {fdp_lo:.0%} to {fdp_hi:.0%}) would be expected on "
             f"shuffled sequence of this window, so roughly {real} are above "
             f"chance"]
    if top:
        parts.append("among the best-scoring sites the share is lower: "
                     + ", ".join(f"{t['fdp']:.0%} of the top {t['top']}"
                                 for t in top))
    out["fdp_note"] = (
        "; ".join(parts) + f". Estimated from {len(counts)} shuffled copies "
        f"that keep the repeat mask and every dinucleotide count. The "
        f"formula that assumes each motif matches at exactly its stated "
        f"p-value would expect {formula:.0f} and is not used. Nothing here "
        f"says which sites are real; only ChIP evidence separates them.")
    return out
