"""Which factors are demonstrably active, inferred from their targets.

The problem this exists for
---------------------------
Everything upstream of here guesses at transcription factors from
*sequence*. A PWM scan of a 6 kb promoter window returns dozens of sites,
and `motifs.scan_quality` is blunt about what that means: at a p-value of
1e-4 over a window that size, most calls are noise. A motif is evidence
that a factor *could* bind, under no constraint that it is present, active,
or doing anything in the cell you care about.

This asks the opposite question from the opposite data. A factor that is
actually active leaves a footprint in the expression of the genes it
regulates: activated targets up, repressed targets down. So instead of
asking "is there a binding site", ask "do this factor's known targets move
together, in the direction its known mode of regulation predicts". That is
the VIPER idea (Alvarez et al., Nat Genet 2016), and it is the closest this
package gets to measuring a factor rather than predicting one.

Two things it gives you
-----------------------
1. **A corrective on the motif scan.** A factor with a strong site and a
   flat regulon is probably a false positive of the scan. A factor with a
   weak site and a coherently perturbed regulon is probably real. The two
   lines of evidence are independent, which is what makes the combination
   worth more than either.

2. **A drug source that needs no annotated binding.** Run the same
   statistic on a drug's own perturbation signature and you learn which
   factors that drug actually moves in cells - whether or not anyone has
   recorded it as a ligand. This finds the ones `drugs.py` cannot, because
   ChEMBL only knows about compounds somebody assayed at the target.

Why the direction is only half-known
------------------------------------
Direction requires a *signed* regulon - you cannot tell activation from
repression without knowing which the edge was. TRRUST labels 46% of its
human edges `Unknown` (4,325 of 9,396), so for nearly half the data the
mode is simply not recorded. Those edges are not thrown away: they support
a two-tailed test ("is this regulon perturbed at all") while only signed
edges support a directional one. Both are reported, and an activity call
made from few signed edges says so.

The limit that matters most
---------------------------
**Inferred activity is not attributable to one factor.** Regulons overlap,
and co-regulated targets cannot distinguish a factor from its neighbour -
if two factors share 80% of their targets, they get near-identical scores
and the statistic cannot say which one moved. VIPER corrects for this
("pleiotropy"); here the overlap is measured and reported per call, so a
score that belongs to a neighbour is visible rather than silent.

Second: the p-values assume targets are independent, which is false.
Co-expression is exactly the structure being measured, so analytic
p-values are anticonservative. Use `permutations=` for an empirical null
when a number is going to be relied on.
"""

from __future__ import annotations

import logging
import math
import statistics as _statistics
from dataclasses import dataclass, field

import numpy as np

log = logging.getLogger(__name__)

# Below this many measured targets the score is noise regardless of how
# extreme it looks: a 3-target regulon can hit |NES| 2 by chance often
# enough to be worthless. VIPER's own implementations use a comparable
# floor.
MIN_TARGETS = 10

# An activity call is reported as directional only when enough of the
# contributing edges carry a mode. Below this it is magnitude-only - the
# regulon moved, but which way the factor went is not established.
MIN_SIGNED_FRACTION = 0.34

# Regulon overlap above which two factors cannot be told apart by this
# method at all, and the call is flagged as unattributable.
CONFOUNDED_JACCARD = 0.5

_SIGN = {
    "activation": 1, "activating": 1, "a": 1, "+": 1, "1": 1, "up": 1,
    "repression": -1, "repressing": -1, "inhibition": -1, "r": -1,
    "-": -1, "-1": -1, "down": -1,
}


@dataclass
class TFActivity:
    """One factor's inferred activity against one signature."""
    tf: str
    nes: float | None                 # directional, signed; None if unsigned
    p_value: float | None
    q_value: float | None = None
    magnitude: float | None = None    # two-tailed, uses unsigned edges too
    magnitude_p: float | None = None
    direction: str = "unknown"        # activated | repressed | unchanged
    n_measured: int = 0               # targets present in the signature
    n_signed: int = 0                 # of those, edges with a known mode
    n_regulon: int = 0                # the full regulon
    coverage: float = 0.0             # n_measured / n_regulon
    confounded_with: str | None = None
    confounded_jaccard: float | None = None
    method: str = ""
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


# --------------------------------------------------------------- regulons
# Where a regulon file is looked for when none is given. CollecTRI first:
# it is fully signed and covers far more factors than TRRUST, which is the
# single biggest limit on how many factors can be tested at all.
DEFAULT_REGULON_PATHS = (
    "data/collectri.tsv", "data/collectri.csv",
    "data/trrust.tsv",
)

COLLECTRI_SOURCE = (
    "CollecTRI: https://github.com/saezlab/CollecTRI, or in Python "
    "`import decoupler as dc; dc.get_collectri(organism='human')` then "
    "`.to_csv('data/collectri.tsv', sep='\\t', index=False)`"
)


def autoload_regulons(*, min_targets: int = 5) -> tuple:
    """(regulons, path) from the first default location that exists.

    Returns ({}, None) when none do, rather than raising: a missing
    regulon file costs you the activity section and nothing else.
    """
    from pathlib import Path as _P
    for cand in DEFAULT_REGULON_PATHS:
        if _P(cand).exists():
            return load_signed_regulons(cand, min_targets=min_targets), cand
    return {}, None


def load_signed_regulons(path, *, min_targets: int = 5) -> dict:
    """Signed regulons from a TRRUST, DoRothEA or CollecTRI table.

    Returns `{TF: {target: sign}}` with sign in {+1, -1, 0}; 0 is an edge
    that exists but whose mode was not recorded, which is 46% of TRRUST and
    is kept deliberately - see the module docstring.

    Accepts the two layouts these tables actually ship in:

      TRRUST     TF <tab> target <tab> Activation|Repression|Unknown <tab> pmid
      CollecTRI  source, target, weight, ...   with weight +1 or -1

    and either tab or comma separated, sniffed from the first line rather
    than from the file extension, because both are routinely saved under
    the wrong one.

    CollecTRI is the better input and worth fetching. It is fully signed
    where TRRUST leaves 46% of edges unlabelled, and it covers roughly
    three times as many factors, which is the single biggest limit on how
    many can be tested at all. See COLLECTRI_SOURCE.

    `regulon.load_regulons` already parses the same file but returns bare
    sets, discarding the mode column. That is the right shape for counting
    regulon size and the wrong one here, so this parses it again rather
    than widening that function's contract.
    """
    import csv
    from pathlib import Path

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")

    # Sniff the separator. Hardcoding a tab made a CollecTRI CSV load as
    # zero regulons and report no error, which is the worst way to fail.
    with open(path, "r", encoding="utf-8", errors="replace") as probe:
        head = probe.readline()
    delim = "," if head.count(",") > head.count("\t") else "\t"

    raw: dict = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.reader(fh, delimiter=delim):
            if len(row) < 2:
                continue
            tf, target = row[0].strip().upper(), row[1].strip().upper()
            if not tf or not target or tf == target:
                continue
            # Header rows and DoRothEA's `source/target/weight` layout.
            if tf in ("TF", "SOURCE", "REGULATOR"):
                continue
            mode = row[2].strip().lower() if len(row) > 2 else ""
            sign = _SIGN.get(mode, 0)
            if not sign and mode:
                # DoRothEA ships a numeric weight rather than a word.
                try:
                    w = float(mode)
                    sign = 1 if w > 0 else (-1 if w < 0 else 0)
                except ValueError:
                    sign = 0
            prev = raw.setdefault(tf, {}).get(target)
            # Contradictory edges for the same pair do occur in TRRUST
            # (activation in one paper, repression in another). Neither
                # wins: the pair is demoted to unsigned, which is the honest
            # reading and keeps it in the two-tailed test.
            raw[tf][target] = 0 if (prev is not None and prev != sign) else sign

    out = {tf: tg for tf, tg in raw.items() if len(tg) >= min_targets}
    n_signed = sum(1 for tg in out.values() for s in tg.values() if s)
    n_total = sum(len(tg) for tg in out.values())
    log.info("signed regulons: %d factor(s), %d edge(s), %d signed (%.0f%%)",
             len(out), n_total, n_signed,
             100.0 * n_signed / n_total if n_total else 0.0)
    return out


def signed_from_sets(regulons: dict) -> dict:
    """Adapt `regulon.load_regulons`-style `{TF: set}` to this module.

    Every edge becomes unsigned, so the result supports the two-tailed test
    and no directional one. Offered so a caller already holding unsigned
    regulons gets the magnitude answer rather than nothing.
    """
    return {tf: {t.upper(): 0 for t in targets}
            for tf, targets in (regulons or {}).items()}


# ------------------------------------------------------------- statistics
def _norm_sf(z: float) -> float:
    """Upper-tail standard normal probability, via erfc. No scipy."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _two_sided(z: float) -> float:
    return min(1.0, 2.0 * _norm_sf(abs(z)))


def rank_normalise(signature: dict) -> dict:
    """Signature scores replaced by their inverse-normal rank transform.

    The statistic below assumes each target's score is a standard normal
    draw under the null. A raw signature is not - a logFC distribution is
    heavy-tailed and its scale is arbitrary - so ranking and then mapping
    to normal quantiles makes the null exact and the result invariant to
    any monotone rescaling of the input. Ties share their mean rank.
    """
    items = [(g, float(v)) for g, v in (signature or {}).items()
             if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not items:
        return {}
    genes = [g for g, _ in items]
    vals = np.asarray([v for _, v in items], dtype=float)
    n = len(vals)

    order = np.argsort(vals, kind="mergesort")
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)
    # average ranks within ties, so a signature of mostly zeros does not
    # acquire an arbitrary ordering
    _, inv, counts = np.unique(vals, return_inverse=True, return_counts=True)
    if counts.max() > 1:
        sums = np.zeros(counts.size)
        np.add.at(sums, inv, ranks)
        ranks = (sums / counts)[inv]

    q = ranks / (n + 1.0)                      # open interval, no infinities
    z = np.sqrt(2.0) * _erfinv(2.0 * q - 1.0)
    return dict(zip(genes, z.tolist()))


def _erfinv(x):
    """Inverse error function, vectorised. numpy has no erfinv."""
    x = np.clip(np.asarray(x, dtype=float), -0.999999, 0.999999)
    # Giles' rational approximation; max abs error ~1e-7 over this range,
    # far below the resolution any of these p-values are read at.
    w = -np.log(1.0 - x * x)
    out = np.empty_like(x)

    lo = w < 5.0
    ww = w[lo] - 2.5
    p = 2.81022636e-08
    for c in (3.43273939e-07, -3.5233877e-06, -4.39150654e-06,
              0.00021858087, -0.00125372503, -0.00417768164,
              0.246640727, 1.50140941):
        p = p * ww + c
    out[lo] = p * x[lo]

    hi = ~lo
    if hi.any():
        ww = np.sqrt(w[hi]) - 3.0
        p = -0.000200214257
        for c in (0.000100950558, 0.00134934322, -0.00367342844,
                  0.00573950773, -0.0076224613, 0.00943887047,
                  1.00167406, 2.83297682):
            p = p * ww + c
        out[hi] = p * x[hi]
    return out


# A q-value that rounds to 0.0 is the most significant result possible, and
# `q or 1.0` would read it as the least - `0.0` is falsy. Every comparison
# goes through this instead, and `_bh` floors at Q_FLOOR so a reported
# number never claims to be exactly zero either.
Q_FLOOR = 1e-6


def _q(value, default: float = 1.0) -> float:
    """A q- or p-value for comparison, with None (not 0.0) as missing."""
    return default if value is None else float(value)


def _bh(pvals: list) -> list:
    """Benjamini-Hochberg q-values, order preserved."""
    vals = [(i, p) for i, p in enumerate(pvals) if p is not None]
    if not vals:
        return [None] * len(pvals)
    vals.sort(key=lambda t: t[1])
    m = len(vals)
    out = [None] * len(pvals)
    running = 1.0
    for rank in range(m, 0, -1):
        i, p = vals[rank - 1]
        running = min(running, p * m / rank)
        out[i] = max(Q_FLOOR, round(min(1.0, running), 6))
    return out


def _overlap(targets: set, regulons: dict, tf: str) -> tuple:
    """The factor whose regulon this one is least distinguishable from."""
    best, best_j = None, 0.0
    for other, tg in regulons.items():
        if other == tf:
            continue
        o = set(tg)
        if not o:
            continue
        j = len(targets & o) / float(len(targets | o))
        if j > best_j:
            best, best_j = other, j
    return best, round(best_j, 3)


def activity(signature: dict, regulons: dict, *,
             min_targets: int = MIN_TARGETS,
             permutations: int = 0,
             check_confounding: bool = True,
             seed: int = 0) -> list:
    """Inferred activity for every factor with enough measured targets.

    `signature` maps gene symbol to a continuous score - logFC, a t
    statistic, a z score, anything monotone in "how much this gene moved".
    Sign convention: positive means the gene went up.

    Two statistics per factor, because TRRUST only signs half its edges:

    * **directional** (`nes`) over signed edges only. Under the null that
      regulon membership is unrelated to the signature, the rank-normalised
      scores are standard normal, so `sum(sign_i * z_i) / sqrt(n)` is too.
      Positive means the factor looks activated.
    * **magnitude** (`magnitude`) over every edge including unsigned ones.
      `|z|` has mean sqrt(2/pi) and variance 1 - 2/pi, so the centred and
      scaled sum is again standard normal. This answers "did this regulon
      move at all", which is the only question the unsigned half of TRRUST
      can address.

    `permutations` replaces the analytic p-value with an empirical null
    built by resampling target sets of the same size. Slower, and worth it
    whenever a number will be relied on, because the analytic null assumes
    targets are independent and co-expression guarantees they are not.
    """
    z = rank_normalise(signature)
    if not z:
        log.warning("empty signature - no activity can be inferred")
        return []

    E_ABS = math.sqrt(2.0 / math.pi)
    SD_ABS = math.sqrt(1.0 - 2.0 / math.pi)
    zvals = np.asarray(list(z.values()), dtype=float)
    rng = np.random.default_rng(seed)

    out: list = []
    for tf, targets in (regulons or {}).items():
        measured = [(g, s) for g, s in targets.items() if g in z]
        if len(measured) < min_targets:
            continue
        n = len(measured)
        scores = np.asarray([z[g] for g, _ in measured], dtype=float)
        signs = np.asarray([s for _, s in measured], dtype=float)
        n_signed = int(np.count_nonzero(signs))

        # magnitude: every measured edge contributes
        mag = float((np.abs(scores).sum() - n * E_ABS) / (SD_ABS * math.sqrt(n)))
        mag_p = _two_sided(mag)

        nes = p = None
        if n_signed >= max(3, int(math.ceil(min_targets * MIN_SIGNED_FRACTION))):
            nes = float((signs * scores).sum() / math.sqrt(n_signed))
            p = _two_sided(nes)

        if permutations:
            # Empirical null: same-size random target sets from the same
            # signature, which keeps the signature's own structure and only
            # breaks the regulon assignment.
            idx = rng.integers(0, zvals.size, size=(permutations, n))
            draw = zvals[idx]
            null_mag = ((np.abs(draw).sum(axis=1) - n * E_ABS)
                        / (SD_ABS * math.sqrt(n)))
            mag_p = float((np.abs(null_mag) >= abs(mag)).mean())
            if nes is not None:
                s_signed = signs[signs != 0]
                d2 = draw[:, :s_signed.size] * s_signed
                null_nes = d2.sum(axis=1) / math.sqrt(n_signed)
                p = float((np.abs(null_nes) >= abs(nes)).mean())

        a = TFActivity(
            tf=tf, nes=(round(nes, 3) if nes is not None else None),
            p_value=(round(p, 6) if p is not None else None),
            magnitude=round(mag, 3), magnitude_p=round(mag_p, 6),
            n_measured=n, n_signed=n_signed, n_regulon=len(targets),
            coverage=round(n / float(len(targets)), 3),
            method=("empirical null, %d permutations" % permutations
                    if permutations else "analytic normal null"),
        )

        if nes is None:
            a.direction = "unknown"
            a.caveats.append(
                f"only {n_signed} of {n} measured edges carry a mode, so the "
                f"regulon can be tested for movement but the factor's "
                f"direction cannot be called")
        elif p is not None and p < 0.05:
            a.direction = "activated" if nes > 0 else "repressed"
        else:
            a.direction = "unchanged"

        if a.coverage < 0.2:
            a.caveats.append(
                f"only {a.coverage:.0%} of the regulon was measured in this "
                f"signature, so this is a sample of the factor's targets, "
                f"not its footprint")
        if permutations == 0:
            a.caveats.append(
                "analytic p-value: assumes targets are independent, which "
                "co-expression makes false, so it is anticonservative - "
                "re-run with permutations= before relying on it")

        if check_confounding:
            other, j = _overlap(set(targets), regulons, tf)
            a.confounded_with, a.confounded_jaccard = other, j
            if j >= CONFOUNDED_JACCARD:
                a.caveats.append(
                    f"regulon overlaps {other} at Jaccard {j} - this method "
                    f"cannot attribute the signal to one of the two, and the "
                    f"activity may be {other}'s")

        out.append(a)

    qs = _bh([a.p_value for a in out])
    for a, q in zip(out, qs):
        a.q_value = q

    out.sort(key=lambda a: -(abs(a.nes) if a.nes is not None
                             else abs(_q(a.magnitude, 0.0))))
    log.info("activity inferred for %d factor(s); %d directional, %d "
             "significant at q<0.05", len(out),
             sum(1 for a in out if a.nes is not None),
             sum(1 for a in out if _q(a.q_value) < 0.05))
    return out


def activity_from_sets(up: set, down: set, regulons: dict, *,
                       min_targets: int = MIN_TARGETS,
                       background: int = 12_000) -> list:
    """Activity from up/down gene SETS rather than a continuous signature.

    This is the shape LINCS data arrives in: Enrichr's
    `LINCS_L1000_Chem_Pert` .gmt files give the top up and top down genes
    per perturbation and no scores. A rank transform of a ternary vector
    would invent precision that is not there, so the statistic changes to
    match the data: among a factor's signed targets that appear in either
    set, count how many agree with the edge's mode (activation + up, or
    repression + down) and test against a coin flip.

    That is a binomial sign test. It is weaker than the continuous version
    - it uses only the targets that made the top-N cut - and it is the
    honest test for set data.
    """
    up = {g.strip().upper() for g in (up or set())}
    down = {g.strip().upper() for g in (down or set())}
    hit = up | down
    if not hit:
        return []

    out: list = []
    for tf, targets in (regulons or {}).items():
        pairs = [(g, s) for g, s in targets.items() if s and g in hit]
        if len(pairs) < min_targets:
            continue
        n = len(pairs)
        agree = sum(1 for g, s in pairs
                    if (s > 0 and g in up) or (s < 0 and g in down))
        # (k - n/2) / sqrt(n/4), the normal approximation to binomial(n, 1/2)
        z = (agree - n / 2.0) / math.sqrt(n / 4.0)
        p = _two_sided(z)
        a = TFActivity(
            tf=tf, nes=round(z, 3), p_value=round(p, 6),
            n_measured=n, n_signed=n, n_regulon=len(targets),
            coverage=round(n / float(len(targets)), 3),
            method=f"binomial sign test on {n} signed target(s) in the set",
            direction=("activated" if z > 0 else "repressed") if p < 0.05
                      else "unchanged",
        )
        a.caveats.append(
            f"set-based: only targets inside the perturbation's top-N lists "
            f"are used ({n} of {len(targets)}), so a factor whose targets "
            f"moved moderately is invisible here")
        if len(hit) < 50:
            a.caveats.append(
                f"the signature has only {len(hit)} genes, which makes any "
                f"regulon overlap small and the test underpowered")
        out.append(a)

    qs = _bh([a.p_value for a in out])
    for a, q in zip(out, qs):
        a.q_value = q
    out.sort(key=lambda a: -abs(_q(a.nes, 0.0)))
    return out


# ------------------------------------------------- the tissue-level answer
def tissue_signature(tissue: str, *, table: dict | None = None,
                     floor: float = 0.1) -> dict:
    """A contrast signature for one tissue: log2(tissue / median tissue).

    Built from the expression table `delivery.load_expression_file` already
    holds, so `--expression-file` is the only input needed.

    One subtlety that would otherwise corrupt every score: that table is
    *sparse*, because it is loaded with a `min_value` floor and a gene only
    carries the tissues where it cleared it. A gene missing from this
    tissue is therefore low, not unmeasured, so it is set to `floor` rather
    than dropped. Dropping it would bias every regulon toward the genes
    that happen to be broadly expressed.
    """
    if table is None:
        from . import delivery
        table = delivery._LOCAL_EXPRESSION
    if not table:
        log.warning("no local expression table loaded - tissue activity "
                    "needs --expression-file")
        return {}

    want = (tissue or "").strip().lower()
    sig: dict = {}
    for gene, per_tissue in table.items():
        if not per_tissue:
            continue
        here = None
        for t, v in per_tissue.items():
            if t.strip().lower() == want:
                here = float(v)
                break
        # statistics.median, not values[n // 2]: for an even number of
        # tissues the latter returns the upper of the two middle values, so
        # a gene measured in exactly two tissues got its own maximum as the
        # baseline and scored 0.0 in the tissue where it was highest.
        median = _statistics.median(float(v) for v in per_tissue.values())
        if median <= 0:
            continue
        sig[gene.upper()] = math.log2(max(here if here is not None else floor,
                                          floor) / max(median, floor))
    if not sig:
        log.warning("tissue %r not found in the expression table", tissue)
    return sig


def tissue_activity(tissue: str, regulons: dict, *, table: dict | None = None,
                    permutations: int = 0, **kw) -> list:
    """Which factors look active in one tissue, from their regulons.

    The corrective on the motif scan: a factor with a strong site here and
    a flat regulon is probably a false positive of the scan.

    What this is NOT: a measurement of the factor in *your* sample. It is
    the factor's footprint in a reference tissue average, which says
    nothing about a patient, a cell line, or a disease state. And a
    ligand-activated receptor can be transcribed everywhere while being
    active only where its ligand reaches, so for nuclear receptors in
    particular a quiet regulon here does not mean the receptor is absent -
    it means nothing is switching it on in healthy tissue.
    """
    sig = tissue_signature(tissue, table=table)
    if not sig:
        return []
    out = activity(sig, regulons, permutations=permutations, **kw)
    for a in out:
        a.caveats.append(
            "reference tissue average, not your sample: it cannot speak to a "
            "disease state, a cell line, or an individual")
    return out


def summarise(activities: list, *, scanned: list | None = None) -> str:
    """Prose for the report, including the comparison with the motif scan."""
    if not activities:
        return ("No factor had enough measured targets to infer activity. "
                "That is a statement about the regulon table and the "
                "signature, not about the factors: with TRRUST, most of the "
                "~1,600 human transcription factors have too few curated "
                "targets to test at all.")

    sig = [a for a in activities if _q(a.q_value) < 0.05
           and a.direction != "unknown"]
    directional = [a for a in activities if a.nes is not None]
    lines = [
        f"Activity was testable for {len(activities)} factor(s) - those with "
        f"at least {MIN_TARGETS} curated targets measured in this signature. "
        f"{len(directional)} had enough signed edges for a direction, and "
        f"{len(sig)} reached q<0.05."]

    if scanned:
        scanned_set = {s.upper() for s in scanned}
        both = [a for a in sig if a.tf.upper() in scanned_set]
        quiet = [a for a in activities if a.tf.upper() in scanned_set
                 and a.direction == "unchanged"]
        if both:
            lines.append(
                f"**Against the motif scan:** {len(both)} factor(s) have "
                f"both a binding site in the promoter and a coherently "
                f"perturbed regulon "
                f"({', '.join(a.tf for a in both[:6])}"
                f"{' and others' if len(both) > 6 else ''}). Two independent "
                f"lines of evidence agreeing is the strongest thing this "
                f"package can offer, so start there.")
        else:
            lines.append(
                "**Against the motif scan: no factor has both a site here "
                "and a perturbed regulon.** That is a real result, not a "
                "missing one. It says the factors this promoter has sites "
                "for show no functional footprint in this tissue, so the "
                "routes through them rest on sequence alone - and a route "
                "whose factor cannot be shown to be active is a weaker "
                "starting point than the drug tables' ranking suggests.")
        if quiet:
            lines.append(
                f"{len(quiet)} factor(s) with a site had a flat regulon "
                f"({', '.join(a.tf for a in quiet[:6])}"
                f"{' and others' if len(quiet) > 6 else ''}). A scan over a "
                f"6 kb window returns mostly noise, so a site with no "
                f"functional footprint is the expected failure mode rather "
                f"than a surprise.")
    unattributable = [a for a in sig
                      if _q(a.confounded_jaccard, 0.0) >= CONFOUNDED_JACCARD]
    if unattributable:
        lines.append(
            f"{len(unattributable)} significant call(s) sit on a regulon that "
            f"overlaps another factor's too heavily to attribute - the "
            f"signal is real and which factor produced it is not decided.")
    return " ".join(lines)


# ============================================================= drug source
# Which approved drugs move a factor, with no annotated binding required.
#
# `drugs.py` can only find compounds somebody assayed at the target, so a
# factor with no medicinal chemistry behind it returns nothing however
# important it is - the attrition that leaves every drug table full of
# nuclear receptors. This asks a different question of the LINCS data the
# package already has on disk: run the same regulon statistic on each
# compound's own perturbation signature, and a compound whose signature
# implies the factor's activity changed is a candidate, whether or not it
# has ever been recorded as a ligand.
#
# It is a weaker claim than a binding assay and a different kind of claim.
# A drug can move a factor's regulon by acting three steps upstream, or by
# killing the cell. What it establishes is that the factor's targets move
# when you give the drug, which is the thing a promoter analysis is trying
# to predict and never observes.


@dataclass
class DrugActivityHit:
    """One drug's effect on one factor, pooled over LINCS signatures."""
    drug: str
    tf: str
    direction: str                  # activated | repressed
    median_z: float
    n_signatures: int               # how many times this drug was profiled
    n_significant: int              # of those, how many reached q<0.05
    n_agreeing: int                 # how many moved it the same way
    consistency: float              # n_agreeing / n_signatures
    cell_lines: list = field(default_factory=list)
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


# A drug profiled once that looks significant is the weakest evidence LINCS
# can give: 33,000 signatures against one regulon will produce significant
# single hits by the thousand. Requiring repetition is the cheap defence,
# and LINCS supplies it for free because most compounds are profiled in
# several cell lines and doses.
MIN_SIGNATURES = 3
MIN_CONSISTENCY = 0.6

# Two more gates, both earned by looking at what got through without them.
#
# The scan tests every one of ~33,000 signatures against one regulon, so a
# nominal p<0.05 is met about 1,650 times by chance alone. Significance is
# counted after Benjamini-Hochberg across the signatures for that factor,
# not before.
#
# And consistency on its own is not enough. SP1 returned mitoxantrone at
# 72% consistency over 32 profiles with a median z of -0.30: agreement
# about an effect too small to matter. A floor on the median effect is what
# separates "moves the regulon" from "leans one way".
MIN_MEDIAN_Z = 1.0
MIN_SIGNIFICANT = 2


def _read_gmt_restricted(path, keep: set) -> dict:
    """{signature name: set of genes}, keeping only genes in `keep`.

    The .gmt files are 16 MB and 33,000 lines, and all that is ever needed
    from a line is its overlap with one regulon. Intersecting during the
    read keeps the whole thing small and makes a scan over every compound
    affordable.
    """
    from pathlib import Path as _P
    out: dict = {}
    with open(_P(path), "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            genes = {g.strip().upper() for g in fields[2:] if g.strip()}
            hit = genes & keep
            if hit:
                out[fields[0]] = hit
    return out


def drugs_moving_factor(tf: str, regulon: dict, *, up_path=None,
                        down_path=None, min_targets: int = MIN_TARGETS,
                        min_signatures: int = MIN_SIGNATURES,
                        min_consistency: float = MIN_CONSISTENCY,
                        min_median_z: float = MIN_MEDIAN_Z,
                        min_significant: int = MIN_SIGNIFICANT,
                        top_n: int = 25) -> list:
    """Approved compounds whose LINCS signature moves this factor's regulon.

    Returns `DrugActivityHit`s pooled per drug. A drug appears in LINCS many
    times over, at different cell lines, doses and durations, and that
    repetition is the whole value: a compound that moves the factor the same
    way in seven profiles out of eight is saying something, and one that
    does it in one profile out of nine is noise that a single-signature
    scan would have reported identically.
    """
    from . import signatures as _sig

    signed = {g: s for g, s in (regulon or {}).items() if s}
    if len(signed) < min_targets:
        log.info("%s has only %d signed targets; too few to read a drug "
                 "signature against", tf, len(signed))
        return []
    keep = set(signed)

    ups = _resolve_gmt(up_path, "up")
    downs = _resolve_gmt(down_path, "down")
    if not ups or not downs:
        return []

    up_sets = _read_gmt_restricted(ups, keep)
    down_sets = _read_gmt_restricted(downs, keep)

    # One pass over every signature, then a single multiplicity correction
    # across all of them - this factor was tested ~33,000 times and the
    # per-signature p-value does not know that.
    rows: list = []
    for name in set(up_sets) | set(down_sets):
        hits = activity_from_sets(up_sets.get(name, set()),
                                  down_sets.get(name, set()),
                                  {tf: signed}, min_targets=min_targets)
        if hits:
            rows.append((name, hits[0]))
    if not rows:
        return []

    qs = _bh([a.p_value for _, a in rows])

    per_drug: dict = {}
    for (name, a), q in zip(rows, qs):
        meta = _sig._parse_sig_name(name)
        rec = per_drug.setdefault(meta["drug"], {"z": [], "q": [],
                                                 "cells": set()})
        rec["z"].append(a.nes)
        rec["q"].append(q)
        if meta.get("cell_line"):
            rec["cells"].add(meta["cell_line"])

    out: list = []
    for drug, rec in per_drug.items():
        zs = rec["z"]
        n = len(zs)
        if n < min_signatures:
            continue
        med = float(_statistics.median(zs))
        if med == 0:
            continue
        if abs(med) < min_median_z:
            continue
        agreeing = sum(1 for z in zs if (z > 0) == (med > 0))
        consistency = agreeing / float(n)
        if consistency < min_consistency:
            continue
        n_sig = sum(1 for q in rec["q"] if q is not None and q < 0.05)
        if n_sig < min_significant:
            continue

        h = DrugActivityHit(
            drug=drug, tf=tf,
            direction="activated" if med > 0 else "repressed",
            median_z=round(med, 2), n_signatures=n, n_significant=n_sig,
            n_agreeing=agreeing, consistency=round(consistency, 2),
            cell_lines=sorted(rec["cells"])[:6],
        )
        h.caveats.append(
            "the signature says this factor's targets moved when the "
            "compound was given. It does NOT say the compound binds the "
            "factor - the effect can run through anything upstream, or "
            "through general toxicity")
        h.caveats.append(
            "the commonest confound in this scan is proliferation. A "
            "compound that arrests the cell cycle moves every regulon "
            "enriched for cell-cycle genes, so a growth inhibitor can look "
            "like a specific effect on any such factor. Check whether the "
            "agreeing compounds share an antiproliferative mechanism before "
            "reading the factor as the target.")
        if abs(med) < 2.0:
            h.caveats.append(
                f"the effect is modest (median z {med:+.2f}); it survived "
                f"because it repeated, not because any one profile was "
                f"striking")
        if len(rec["cells"]) < 2:
            h.caveats.append(
                f"all {n} profiles come from one cell line, so this is "
                f"repetition rather than replication across contexts")
        out.append(h)

    out.sort(key=lambda h: (-h.n_significant, -abs(h.median_z)))
    log.info("LINCS drug scan for %s: %d compound(s) move the regulon "
             "consistently", tf, len(out))
    return out[:top_n]


def _resolve_gmt(explicit, direction: str):
    """The .gmt path to use, or None, without raising."""
    from pathlib import Path as _P
    from . import signatures as _sig
    cands = [explicit] if explicit else list(
        _sig.DEFAULT_GMT_PATHS.get(direction, ()))
    for c in cands:
        if c and _P(c).exists():
            return c
    log.info("no local LINCS .gmt for direction %r; the drug scan needs "
             "both files", direction)
    return None


def summarise_drug_scan(results: dict, *, want: str = "up") -> str:
    """Prose for the report."""
    if not results:
        return ("No compound moved any testable factor's regulon "
                "consistently enough to report. That is the common outcome: "
                "the scan asks 33,000 signatures to agree with each other "
                "as well as with the regulon, and most compounds were "
                "profiled too few times to clear the repetition bar.")
    n = sum(len(v) for v in results.values())
    return (f"{n} compound-factor pair(s) across {len(results)} factor(s). "
            f"Each one is a compound whose own perturbation signature moves "
            f"that factor's targets, found without any record of the "
            f"compound binding the factor - which is what makes this "
            f"different from the drug tables above, and why it can reach "
            f"factors that have no medicinal chemistry at all.")
