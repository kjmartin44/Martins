"""Drugs that change how much of the factor exists, not what it does.

The gap this closes
-------------------
Every drug route in this package acts on a factor's *activity*: a ligand in
its pocket, an allosteric site, a kinase upstream, a degradation signal. All
of those need the factor to be druggable, or to sit downstream of something
druggable.

There is a cheaper lever that was never pulled. The perturbation corpora
record, for roughly 20,000 genes, which compounds changed each one's
transcript. The pipeline already queries that for YOUR GENE - the
"empirical signature" rows. It never queried it for the FACTOR. So a
compound that reproducibly doubles VDR transcript was invisible, even
though more VDR plausibly means more of every VDR target including yours.

That matters most exactly where the rest of the package fails. A factor
with no ligand, no structure and no curated upstream kinase still has an
expression level, and that level is measurable and perturbable. This is
the only route here that works for an entirely undruggable factor.

The direction is not obvious
----------------------------
Which way you want the factor to move depends on what it does to your
gene:

    goal   factor role   you want the factor
    ----   -----------   -------------------
    up     activator     higher
    up     repressor     lower
    down   activator     lower
    down   repressor     higher

A factor whose role is unrecorded gives no usable direction, and those are
reported separately rather than guessed into one column. Getting this
backwards would recommend the drug that moves your gene the wrong way, so
the mapping is explicit and tested in all four combinations.

The limit that matters most
---------------------------
**More factor is not more factor activity.** For a ligand-activated
receptor they are nearly unrelated: VDR transcript is near-ubiquitous and
1,25-D signalling is not, so doubling VDR without its ligand may do
nothing at all. Every nuclear receptor in these tables has that problem,
and they are most of the tables. So a route here is strongest for a factor
whose activity tracks its abundance - many bHLH, HOX and ETS factors - and
weakest for the receptors the report is otherwise full of. The caveat is
attached per factor, driven by the nuclear-receptor annotation the pipeline
already computes.

Second: a measured change in the factor's transcript says nothing about
protein, which is set by translation and degradation as well. And third,
this is one more inferential step than the "empirical signature" rows:
there the compound was measured against your gene, here it was measured
against the factor and the factor-to-gene link is an inference. The tier
reflects that.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

log = logging.getLogger(__name__)

# What to do to the factor, given what you want for the gene and what the
# factor does to it. The one place this mapping lives.
WANTED = {
    ("up", "activator"): "up",
    ("up", "repressor"): "down",
    ("down", "activator"): "down",
    ("down", "repressor"): "up",
}

# A compound profiled once is the weakest thing a corpus can say. The
# .gmt holds one consensus signature per cell line, dose and duration, so
# the same compound recurs, and that recurrence is the only cheap evidence
# available here. One appearance is not reported.
MIN_SIGNATURES = 2

# A compound that lowers most of the factors you asked about is not a
# regulator of any of them. Verified on the real corpora: the top rows for
# KLF4 were HSP90 inhibitors (geldanamycin, radicicol, NVP-AUY922) and CDK
# inhibitors (alvocidib, AT-7519, BMS-387032) - agents that shut down
# transcription and proteostasis globally, so they move nearly every
# factor's transcript. Same shape as the promiscuity filter
# `signatures.flag_non_specific` already applies to gene-level hits, keyed
# on factors instead of queries.
PROMISCUITY_SHARE = 0.5
MIN_FACTORS_TO_JUDGE = 4

# Roles that carry no direction. A dual-role factor is genuinely ambiguous
# rather than unknown, and both are unusable here for the same reason.
UNDIRECTED = ("unknown", "dual", "", None)


@dataclass
class AbundanceHit:
    """One compound that moved one factor's transcript."""
    drug: str
    tf: str
    tf_role: str
    moved_factor: str              # "up" | "down" - what the drug did
    wanted: str                    # what we needed it to do
    helps: bool
    tier: int
    n_signatures: int = 1          # distinct contexts it was seen in
    cell_lines: list = field(default_factory=list)
    feature_space: str | None = None
    cell_line: str | None = None
    nuclear_receptor: bool = False
    non_specific: bool = False
    n_factors_hit: int = 1
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def wanted_direction(goal: str, role: str | None) -> str | None:
    """Which way the factor must move, or None when the role cannot say."""
    r = (role or "").strip().lower()
    if r in UNDIRECTED:
        return None
    return WANTED.get((goal, r))


def _tier_for(feature_space: str | None) -> int:
    """One worse than the same measurement made on your own gene.

    `catalogue.from_lincs` calls a landmark measurement of YOUR GENE tier
    2, because the compound was measured against the thing you care about.
    Here it was measured against the factor and the factor-to-gene step is
    an inference, which is the same standing as a curated mechanism.
    """
    if feature_space == "landmark":
        return 3
    if feature_space == "best inferred":
        return 3
    return 4


def find(factors, *, goal: str = "up", up_path=None, down_path=None,
         limit_per_factor: int = 8) -> list:
    """Compounds that moved a factor's own transcript the way you need.

    `factors` is an iterable of dicts carrying at least `symbol` and
    `role`, as `TFResult.to_dict` produces. Both corpora directions are
    read, because a repressor needs the opposite move from an activator and
    both appear in the same run.
    """
    from . import signatures as sig

    from .congeners import salt_free_name

    out: list = []
    # {drug stem: every factor it moved}, over the untruncated lists.
    full_reach: dict = {}
    n_asked = 0
    for tf in factors:
        symbol = (tf.get("symbol") or "").strip()
        if not symbol:
            continue
        role = tf.get("role")
        want = wanted_direction(goal, role)
        if want is None:
            continue
        n_asked += 1

        try:
            hits = sig.load_enrichr_gmt(
                symbol, direction=want,
                up_path=up_path, down_path=down_path)
        except FileNotFoundError:
            log.info("no local corpus for direction %r; factor abundance "
                     "cannot be checked", want)
            return out
        except Exception as exc:                            # noqa: BLE001
            log.debug("abundance lookup failed for %s: %s", symbol,
                      str(exc)[:80])
            continue

        is_nr = bool(tf.get("nuclear_receptor"))

        # One row per compound, not one per signature. The .gmt carries a
        # consensus per cell line and dose, so a compound recurs; the
        # recurrence is evidence and the repetition in the table is not.
        from .congeners import salt_free_name
        grouped: dict = {}
        for h in hits:
            stem = salt_free_name((h.drug or "").lower())
            if not stem:
                continue
            g = grouped.setdefault(stem, {"hit": h, "n": 0, "cells": set()})
            g["n"] += 1
            if getattr(h, "cell_line", None):
                g["cells"].add(h.cell_line)
            if getattr(h, "feature_space", None) == "landmark":
                g["hit"] = h          # prefer a directly measured record

        ranked_groups = sorted(
            (g for g in grouped.values() if g["n"] >= MIN_SIGNATURES),
            key=lambda g: -g["n"])
        # Promiscuity has to be measured over the FULL list. Measuring it
        # after truncating to the top few per factor underestimates it:
        # geldanamycin read as 2 of 5 factors when it reaches more, and so
        # escaped a flag it plainly earns.
        for g in ranked_groups:
            full_reach.setdefault(
                salt_free_name((g["hit"].drug or "").lower()),
                set()).add(symbol)

        for g in ranked_groups[:limit_per_factor]:
            h = g["hit"]
            fs = getattr(h, "feature_space", None)
            hit = AbundanceHit(
                drug=h.drug, tf=symbol, tf_role=str(role),
                moved_factor=want, wanted=want, helps=True,
                tier=_tier_for(fs), feature_space=fs,
                cell_line=getattr(h, "cell_line", None),
                n_signatures=g["n"], cell_lines=sorted(g["cells"])[:6],
                nuclear_receptor=is_nr,
            )
            hit.caveats.append(
                f"measured against {symbol}'s own transcript, not against "
                f"your gene. That {symbol} controls your gene is the "
                f"inference this route rests on")
            if is_nr:
                hit.caveats.append(
                    f"{symbol} is a ligand-activated receptor, so its "
                    f"abundance and its activity are nearly unrelated - "
                    f"more receptor without more ligand may change nothing. "
                    f"Treat this as the weakest route to {symbol}")
            if fs and fs != "landmark":
                hit.caveats.append(
                    f"{symbol} was not measured directly on the array; its "
                    f"value is {fs}")
            hit.caveats.append(
                "transcript, not protein. Translation and degradation set "
                "protein independently")
            if hit.n_signatures < 3:
                hit.caveats.append(
                    f"seen in {hit.n_signatures} contexts, which is the "
                    f"floor for reporting at all")
            out.append(hit)

    # Flag compounds that move most of the factors asked about. With few
    # factors this cannot be judged, so it is not guessed at.
    n_factors = n_asked
    for h in out:
        seen = full_reach.get(salt_free_name(h.drug.lower()), {h.tf})
        h.n_factors_hit = len(seen)
        if (n_factors >= MIN_FACTORS_TO_JUDGE and len(seen) > 1
                and len(seen) / n_factors >= PROMISCUITY_SHARE):
            h.non_specific = True
            h.caveats.insert(0, (
                f"moves {len(seen)} of the {n_factors} factors asked about, "
                f"so it is not a regulator of any one of them. Compounds "
                f"that shut down transcription or proteostasis broadly - "
                f"HSP90 and CDK inhibitors above all - produce this "
                f"pattern"))

    # Tier, then specificity, then a receptor penalty (abundance says least
    # about those), then how many contexts it repeated in.
    out.sort(key=lambda h: (h.tier, h.non_specific, h.nuclear_receptor,
                            -h.n_signatures, h.tf))
    log.info("factor abundance: %d compound-factor pair(s) across %d "
             "factor(s)", len(out), len({h.tf for h in out}))
    return out


def undirected_factors(factors, *, goal: str = "up") -> list:
    """Factors skipped because their role gives no direction.

    Reported rather than silently dropped: for these the corpora may well
    hold a compound that moves the factor, and the reason it is not
    offered is a gap in the role annotation rather than in the data.
    """
    out = []
    for tf in factors:
        if wanted_direction(goal, tf.get("role")) is None:
            sym = (tf.get("symbol") or "").strip()
            if sym:
                out.append({"symbol": sym,
                            "role": tf.get("role") or "unknown"})
    return out


def summarise(hits: list, undirected: list | None = None, *,
              goal: str = "up") -> str:
    verb = "raise" if goal == "up" else "lower"
    if not hits:
        base = ("No compound reproducibly moved any of these factors' own "
                "transcripts in the direction needed. This route needs a "
                "perturbation corpus on disk and a factor whose role is "
                "recorded.")
    else:
        n_tf = len({h.tf for h in hits})
        nr = [h for h in hits if h.nuclear_receptor]
        base = (
            f"{len(hits)} compound-factor pair(s) across {n_tf} factor(s). "
            f"Each is a compound measured to move that factor's own "
            f"transcript the way you need in order to {verb} your gene. "
            f"This is the only route here that works for a factor with no "
            f"ligand, no structure and no curated upstream kinase - it "
            f"needs the factor to be expressed, not druggable.")
        ns = [h for h in hits if h.non_specific]
        if ns:
            base += (
                f" {len(ns)} move most of the factors asked about and are "
                f"flagged non-specific: a compound that lowers every "
                f"factor's transcript is shutting down transcription, not "
                f"regulating yours.")
        if nr:
            base += (
                f" {len(nr)} of them act on a ligand-activated receptor, "
                f"where abundance and activity are nearly unrelated; those "
                f"are the weakest rows here and are sorted last.")
    if undirected:
        base += (
            f" {len(undirected)} factor(s) were skipped because their role "
            f"is not recorded, so there is no way to say which direction "
            f"would help: {', '.join(d['symbol'] for d in undirected[:6])}"
            f"{' and others' if len(undirected) > 6 else ''}.")
    return base
