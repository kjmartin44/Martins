"""How many other genes the literature records a compound as moving.

The finding this exists to surface
----------------------------------
`docs/checks.md` records the hardest result in this package: the ranking's
apparent skill is mostly drug popularity. A baseline that counts only how
many OTHER genes record a drug scores an AUC of 0.93 against CTD; the
routes score 0.53. Among drugs recorded for three or more other genes,
where fame no longer separates the groups, the routes average 0.53 on the
reported-or-not test and 0.50 on the right-direction test. Both are chance.

That was measured across many runs and written in a document. A reader of
one gene's report never saw it, and could not tell whether the compound at
the top was there because of their gene or because it is famous. This puts
the number on the row.

What a count means, and what it does not
----------------------------------------
Tetrachlorodibenzodioxin, benzo(a)pyrene and bisphenol A are recorded for
every gene in the saved runs here. A chemical recorded for all of them is
telling you about the literature - what toxicologists dose cells with -
and not about your gene.

It is NOT a claim that a famous compound is wrong. Acetaminophen is
recorded for many genes and may still be right for one of them. The claim
is narrower: a compound that appears for everything is weak evidence about
anything in particular, so it loses a little score within its tier and
keeps its tier. The same shape as the LINCS breadth caution, from a
different source - breadth measures promiscuity in an assay, this measures
it in the literature, and it is the literature one that fooled the
ranking, because the ranking's inputs are literature-derived.

The table comes from runs you have already done
-----------------------------------------------
`scripts/build_fame_table.py` counts, over the saved results, how many
distinct genes name each chemical. So the resolution is however many genes
you have run: with thirteen, a count of 13 means "every gene measured
here", which is strong, and a count of 2 means very little. The table
records the denominator and every message states it.
"""

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)

DEFAULT_PATHS = ("data/ctd_fame.tsv",)

# Fewer genes than this and a count says nothing about how well known a
# compound is. Matches `retrospective.FAME_MIN_RUNS`, which is where the
# same judgement was first made.
MIN_GENES = 5
# A compound in this percentile of the table's counts is called
# promiscuous. A percentile rather than a fixed count, because the counts
# only go as high as the number of genes measured.
FAME_PERCENTILE = 90
# And it must be recorded for at least this share of them, so a thin
# table cannot make an unremarkable compound the 90th percentile of
# nothing.
MIN_SHARE = 0.4
# How much score a promiscuous compound loses, within its tier. Slight on
# purpose: being well studied is not evidence of being wrong, and this is
# a caution about what the row can tell you, not a verdict on the drug.
# It never changes a tier.
FAME_SCORE_FACTOR = 0.85

SOURCE = ("Build it from your saved runs with "
          "`python scripts/build_fame_table.py`.")


def _key(name: str) -> str:
    from .signatures import breadth_key
    return breadth_key(name)


def write_table(path, counts: dict, *, n_genes: int) -> int:
    """Write `{chemical: n_genes}` as a TSV, with the denominator."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    lines = [f"# genes_measured\t{n_genes}", "chemical\tn_genes"]
    lines += [f"{name}\t{n}" for name, n in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(rows)


def load(path=None) -> dict | None:
    """`{"n_genes", "counts", "cut"}`, or None when there is no table.

    `counts` is keyed the same way as the breadth and proliferation
    checks, so one compound's spellings meet. `cut` is the count at or
    above which a compound is called promiscuous.
    """
    cands = [path] if path else list(DEFAULT_PATHS)
    found = next((Path(c) for c in cands if c and Path(c).exists()), None)
    if found is None:
        return None
    n_genes, counts = 0, {}
    try:
        for line in found.read_text(encoding="utf-8").splitlines():
            if line.startswith("# genes_measured"):
                n_genes = int(line.split("\t")[1])
                continue
            if not line.strip() or line.startswith(("#", "chemical\t")):
                continue
            name, _, n = line.partition("\t")
            key = _key(name)
            if key:
                counts[key] = max(counts.get(key, 0), int(n))
    except (OSError, ValueError) as exc:
        log.warning("could not read the fame table %s (%s); no candidate "
                    "will be checked for being a generic responder",
                    found, type(exc).__name__)
        return None
    if n_genes < MIN_GENES or not counts:
        log.info("the fame table covers %d gene(s), fewer than the %d a "
                 "count needs to mean anything", n_genes, MIN_GENES)
        return None

    ranked = sorted(counts.values())
    i = min(len(ranked) - 1, int(len(ranked) * FAME_PERCENTILE / 100))
    cut = max(ranked[i], MIN_SHARE * n_genes)
    return {"n_genes": n_genes, "counts": counts, "cut": cut,
            "path": str(found)}


def assess(drugs, table: dict | None) -> dict:
    """`{drug key: record}` for the compounds the table has a count for."""
    if not table:
        return {}
    out = {}
    for d in drugs or []:
        key = _key(d)
        n = table["counts"].get(key)
        if not key or n is None:
            continue
        out[key] = {
            # The name as the caller spells it, so the report shows
            # "Valproic acid" and not the key "valproicacid".
            "name": d,
            "n_genes": n,
            "of": table["n_genes"],
            "share": round(n / table["n_genes"], 3),
            "cut": table["cut"],
            "percentile": FAME_PERCENTILE,
            "promiscuous": n >= table["cut"],
        }
    return out


def note_for(rec: dict) -> str:
    """The sentence that goes on a route through a generic responder."""
    return (
        f"recorded as moving {rec['n_genes']} of the {rec['of']} genes "
        f"measured here ({rec['share']:.0%}). A compound the literature "
        f"records for most genes is weak evidence about any one of them: "
        f"the ranking has been tested against this and orders such drugs "
        f"no better than chance. It keeps its tier - being well studied "
        f"is not evidence of being wrong - and loses a little score")


def summarise(assessed: dict, table: dict | None) -> str:
    """One paragraph for the report."""
    if not table:
        return (f"**No fame table, so no candidate was checked for being a "
                f"generic responder.** {SOURCE}")
    if not assessed:
        return ("**None of the candidates appears in the fame table**, so "
                "none could be checked for being a compound the literature "
                "records for every gene. That is a gap, not a pass.")
    n = sum(1 for v in assessed.values() if v["promiscuous"])
    head = (f"**{n} of {len(assessed)} checked candidate(s) are recorded "
            f"for most of the {table['n_genes']} genes measured here.**"
            if n else
            f"**No candidate is recorded for most of the "
            f"{table['n_genes']} genes measured here.**")
    return (
        f"{head} This is the control that matters most for reading the "
        f"order below. Tested across these runs, a baseline that only "
        f"counts how many other genes record a drug predicts CTD better "
        f"(AUC 0.93) than every route in this package combined (0.53), and "
        f"among drugs of equal fame the routes fall to chance. A compound "
        f"recorded for every gene is telling you about the literature, not "
        f"about yours. Promiscuous rows lose a little score within their "
        f"tier and say so; none changes tier.")
