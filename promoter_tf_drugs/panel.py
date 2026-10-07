"""Compare several candidate genes and narrow to the cleanest route.

The question this answers
-------------------------
The rest of the package asks "which drugs move THIS gene". If you have five
candidate targets and must pick one, that is the wrong unit: you need the
gene+drug pair with the least collateral effect on everything else, and that
comparison can only be made across genes.

**These genes are alternatives, not a panel to move together.** That inverts
the usual reading of a recurring hit: a drug appearing against four of your
five candidates is not four times validated, it is a drug that moves a lot of
things and one of them happens to be whichever gene you pick. `n_genes = 1`
is the good case.

A funnel, not a score
---------------------
Collapsing specificity into one number hides which criterion did the work.
Instead each stage applies one criterion and reports what it removed, so an
empty shortlist names the binding constraint and a stage that did nothing
says so. Three properties make the accounting honest:

**Missing data passes, and is counted.** Dropping a candidate for lacking a
field its route never had is backwards — the rule `shortlist._grade_ok`
already encodes. But a stage reporting "removed 0" when it had nothing to
work with reads like a clean bill of health, so `n_not_assessable` is
reported beside it.

**Correlated stages are exposed.** Global regulators tend to have large
regulons *and* broad expression *and* systemic drug classes. Three views of
one property, shown as three independent tests, overstate the evidence. Each
stage reports `n_unique_kills` — rows it removed that no later stage would
have removed anyway.

**Genes, not rows.** With `--all-tss` a gene with four TSSs contributes four
blocks of near-duplicate routes (CCR5's two windows overlapped by 71%), so
counting rows lets the most-annotated gene dominate. Counts are of distinct
(active moiety, target) pairs per gene.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

from . import progress as _progress
from . import (congeners as cong_mod, delivery as dlv_mod,
               paralogs as par_mod, pipeline,
               provenance as prov_mod, regulon as regulon_mod,
               shortlist as sl, systemic as sys_mod)

log = logging.getLogger(__name__)


# ------------------------------------------------------------------ genes
@dataclass
class PanelGene:
    symbol: str
    result: dict | None = None
    error: str | None = None
    excluded: str | None = None

    @property
    def verdict(self) -> dict:
        return (self.result or {}).get("evidence_verdict") or {}

    @property
    def level(self) -> str:
        return self.verdict.get("level", "unknown")

    @property
    def query(self) -> dict:
        return (self.result or {}).get("query") or {}

    def to_dict(self) -> dict:
        return {"gene": self.symbol, "level": self.level,
                "error": self.error, "excluded": self.excluded,
                "reasons": self.verdict.get("reasons", [])}


@dataclass
class Panel:
    genes: list = field(default_factory=list)

    @property
    def included(self) -> list:
        return [g for g in self.genes if g.result and not g.excluded]

    @property
    def dropped(self) -> list:
        return [g for g in self.genes if g.excluded or g.error]


def run_panel(genes, *, out_dir=None, require_interpretable: bool = False,
              **run_kwargs) -> Panel:
    """Run the full pipeline for each gene, tolerating individual failures.

    One gene whose Ensembl lookup dies must not cost you the other four —
    the same reasoning as `pipeline._optional`, one level up.
    """
    panel = Panel()
    genes = list(genes)
    for i, symbol in enumerate(genes, 1):
        pg = PanelGene(symbol=symbol)
        # Several genes run one after another, so each one's checkpoints
        # say which gene this is and how many are left.
        if _progress.enabled():
            _progress.enable(prefix=f"({i}/{len(genes)}) ")
        try:
            pg.result = pipeline.run(symbol, **run_kwargs)
        except Exception as exc:                            # noqa: BLE001
            pg.error = f"{type(exc).__name__}: {str(exc)[:200]}"
            log.warning("gene %s failed and was dropped from the panel: %s",
                        symbol, pg.error)
            panel.genes.append(pg)
            continue

        if out_dir:
            pipeline.write_outputs(pg.result, Path(out_dir) / symbol)

        level = pg.level
        if level == "not interpretable":
            pg.excluded = (
                "evidence verdict 'not interpretable' — its ranking is of "
                "chance motif matches, so letting it compete would let the "
                "weakest evidence supply the winning drug")
        elif require_interpretable and level == "degraded":
            pg.excluded = "evidence verdict 'degraded' and --require-interpretable set"
        panel.genes.append(pg)
    return panel


def parity(panel: Panel) -> dict:
    """Which run options differ across the panel, and is the comparison fair?

    The panel-level form of the bug `provenance` exists to prevent. A gene run
    without --chip-bed loses every chip-dependent stage for reasons that have
    nothing to do with its biology, and would look worse than a gene that had
    the evidence.
    """
    # Every argument that changes what evidence a gene can have. The list
    # was shorter than this because the result only recorded twenty-four
    # arguments; now that it records all of them, the check can see the
    # differences it was written to catch.
    fields = ("goal", "demo", "chip_bed", "chip_api", "cell_type_check",
              "use_lincs", "use_genetics", "regulon_file", "tissue",
              "approved_only", "pvalue", "upstream", "downstream",
              "tf_activity", "activity_tissue", "upstream_network",
              "docking_results", "promoter_dna", "screen_product",
              "sensitivity_scan", "expression_file", "lincs_file",
              "signature_file", "depmap_file", "prism_file",
              "predict_direction", "all_tss", "pfm_file")
    genes = panel.included
    seen: dict = {f: {} for f in fields}
    for g in genes:
        for f in fields:
            seen[f].setdefault(repr(g.query.get(f)), []).append(g.symbol)

    # Not comparable at all, rather than merely biased.
    for hard in ("goal", "demo"):
        if len(seen[hard]) > 1:
            raise ValueError(
                f"cannot compare genes run with different {hard!r}: "
                + "; ".join(f"{v} -> {', '.join(gs)}"
                            for v, gs in seen[hard].items())
                + ". Re-run them with the same setting.")

    differs = {f: v for f, v in seen.items()
               if len(v) > 1 and f not in ("goal", "demo")}
    warnings = []
    for f, variants in differs.items():
        warnings.append(
            f"**{f}** differs across the panel ("
            + "; ".join(f"{val}: {', '.join(gs)}" for val, gs in variants.items())
            + ") — genes without it lose evidence for reasons unrelated to "
              "their biology, so this comparison is not like-for-like.")
    return {"differs": {f: {k: v for k, v in d.items()}
                        for f, d in differs.items()},
            "warnings": warnings, "fair": not warnings}


def merge_results(panel: Panel) -> dict:
    """One synthetic results.json over the included genes.

    Blocks are concatenated so `shortlist.enrich_rows` works unchanged. The
    three top-level fields that are single-gene by construction are moved
    under `panel.per_gene` rather than blended, because merging them would
    silently misattribute one gene's data to another.
    """
    blocks, per_gene = [], {}
    for g in panel.included:
        blocks.extend(g.result.get("results") or [])
        per_gene[g.symbol] = {
            "query": g.query,
            "evidence_verdict": g.verdict,
            "evidence_sources": g.result.get("evidence_sources"),
            # Top-level, not per-block: merging would attribute one gene's
            # direct therapies to every gene in the panel.
            "direct_therapies": g.result.get("direct_therapies") or [],
        }
    goal = (panel.included[0].query.get("goal") if panel.included else "up")
    return {
        "results": blocks,
        "direct_therapies": [],      # deliberately empty; see per_gene
        "query": {"genes": [g.symbol for g in panel.included],
                  "goal": goal, "panel": True,
                  "demo": bool(panel.included
                               and panel.included[0].query.get("demo"))},
        "panel": {"per_gene": per_gene,
                  "dropped": [g.to_dict() for g in panel.dropped]},
    }


# ------------------------------------------------------------ extra axes
def axis_lookup(merged: dict) -> dict:
    """(gene, tf) -> the per-factor facts the panel-only stages need."""
    out: dict = {}
    for block in merged.get("results", []):
        gene = (block.get("promoter") or {}).get("gene_symbol")
        ct = block.get("cell_type") or {}
        n_total = ct.get("n_tissues_total")
        nr = {t["symbol"]: bool(t.get("nuclear_receptor"))
              for t in (block.get("tfs") or [])}
        for sym, ov in (ct.get("tfs") or {}).items():
            out[(gene, sym)] = {"n_tf": ov.get("n_tf"), "n_total": n_total,
                                "nuclear_receptor": nr.get(sym, False)}
        for sym, is_nr in nr.items():
            out.setdefault((gene, sym), {"n_tf": None, "n_total": n_total,
                                         "nuclear_receptor": is_nr})
    return out


def annotate_axes(rows: list, merged: dict, regulons: dict | None = None) -> None:
    """Stamp the three off-target axes onto each row, in place."""
    lookup = axis_lookup(merged)
    for r in rows:
        tf = r.get("acts_on")
        facts = lookup.get((r.get("gene"), tf), {})

        rev = regulon_mod.regulon_evidence(tf or "", regulons)
        r["_regulon_trustworthy"] = rev.trustworthy
        r["_regulon_understudied"] = rev.understudied
        r["_regulon_note"] = rev.note

        tb = dlv_mod.expression_breadth(
            tf or "", n_total=facts.get("n_total"),
            nuclear_receptor=facts.get("nuclear_receptor", False))
        # The block's own n_tf is authoritative when present — it came from
        # the run's expression source, which may not be loaded now.
        if facts.get("n_tf") is not None and tb["score"] is None \
                and not facts.get("nuclear_receptor") and facts.get("n_total"):
            tb = {"score": round(max(0.0, 1 - facts["n_tf"] / facts["n_total"]), 3),
                  "n_tf": facts["n_tf"], "n_total": facts["n_total"],
                  "basis": "expression table", "note": ""}
        r["_tissue_breadth"] = tb["score"]
        r["_tissue_basis"] = tb["basis"]

        # Can a selective ligand exist here at all? An orthosteric pocket is
        # conserved across a receptor family, so THRA/THRB and RARA/RARG
        # cross-react by construction, not by any drug's failing.
        par = par_mod.paralog_specificity(tf or "")
        r["_paralog_score"] = par["score"]
        r["_paralog_family"] = par["family"]
        r["_paralog_relatives"] = par["relatives"]
        r["_paralog_assessed"] = par["assessed"]
        r["_paralog_note"] = par["note"]
        r["_paralog_selective_agents"] = par.get("selective_agents") or []

        burden = sys_mod.clinical_burden(
            r["drug"], notes=r.get("_notes", ""), target=tf)
        r["_burden"] = burden.level
        r["_burden_score"] = burden.score
        r["_burden_assessed"] = burden.assessed
        r["_burden_class"] = burden.drug_class
        r["_burden_reasons"] = burden.reasons


# Panel-only predicates. Missing data passes — the stage counts it instead.
def _p_min_tissue_breadth(row, v):
    s = row.get("_tissue_breadth")
    return s is None or s >= v


_BURDEN_RANK = {lvl: i for i, lvl in enumerate(sys_mod.BURDEN_LEVELS)}


def _p_max_burden(row, v):
    return _BURDEN_RANK.get(row.get("_burden", "low"), 0) <= _BURDEN_RANK[v]


def _p_trustworthy_regulon(row, v):
    return not v or bool(row.get("_regulon_trustworthy"))


def _p_min_paralog_specificity(row, v):
    """Missing passes. A drug with a subtype-selective agent in its family is
    also exempt: penalising trifarotene for RARG's family size would punish
    the molecule developed precisely to solve that problem."""
    if row.get("_paralog_selective_agents"):
        return True
    s = row.get("_paralog_score")
    return s is None or s >= v


PREDICATES = dict(sl.PREDICATES)
PREDICATES.update({
    "min_tissue_breadth": _p_min_tissue_breadth,
    "max_burden": _p_max_burden,
    "trustworthy_regulon": _p_trustworthy_regulon,
    "min_paralog_specificity": _p_min_paralog_specificity,
})

# Which datum each panel axis depends on, so a stage can report that it had
# nothing to work with rather than reporting a reassuring zero.
_ASSESSABLE = {
    "min_tissue_breadth": lambda r: r.get("_tissue_breadth") is not None,
    "max_burden": lambda r: bool(r.get("_burden_assessed")),
    "max_regulon": lambda r: r.get("_regulon_size") is not None,
    "trustworthy_regulon": lambda r: r.get("_regulon_trustworthy") is not None,
    "min_paralog_specificity": lambda r: bool(r.get("_paralog_assessed")),
}


# ----------------------------------------------------------------- funnel
@dataclass(frozen=True)
class Stage:
    name: str
    axis: str
    filters: dict
    why: str


@dataclass
class StageResult:
    name: str
    axis: str
    why: str
    n_in: int = 0
    n_removed: int = 0
    n_remaining: int = 0
    n_not_assessable: int = 0
    n_unique_kills: int = 0
    genes_in: int = 0
    genes_out: int = 0
    genes_lost: list = field(default_factory=list)
    removed_examples: list = field(default_factory=list)
    requested: bool = True      # was any threshold actually set for it?

    @property
    def did_nothing(self) -> bool:
        """Applied, but had no data to apply. Distinct from 'not requested'."""
        return (self.requested and self.n_removed == 0
                and self.n_not_assessable > 0)

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["did_nothing"] = self.did_nothing
        return d


def _pairs(rows) -> dict:
    """Distinct (active moiety, target) pairs per gene — the honest unit.

    Counting rows lets a gene with four TSSs of near-identical windows
    dominate; counting salt-collapsed drug/target pairs does not.
    """
    out: dict = {}
    for r in rows:
        out.setdefault(r.get("gene"), set()).add(
            (cong_mod.salt_free_name(r["drug"]), r.get("acts_on")))
    return out


def _passes(row, stage: Stage) -> bool:
    return all(PREDICATES[k](row, v) for k, v in stage.filters.items()
               if v is not None)


def run_funnel(rows: list, stages: list) -> tuple[list, list]:
    """Apply each stage in turn, recording what it removed and why."""
    for st in stages:
        sl.validate_filters({k: v for k, v in st.filters.items()
                             if k in sl.PREDICATES})
        unknown = set(st.filters) - set(PREDICATES)
        if unknown:
            raise KeyError(f"stage {st.name!r} names unknown filter(s) {unknown}")

    current, results = list(rows), []
    for i, st in enumerate(stages):
        kept, removed = [], []
        for r in current:
            (kept if _passes(r, st) else removed).append(r)

        # How many survived only because the datum was missing.
        # Only ACTIVE filters count. A criterion the user never set is not a
        # stage that failed for want of data — reporting "did nothing, 20
        # unassessable" for a threshold nobody asked for would cry wolf and
        # devalue the warning where it is real.
        checks = [_ASSESSABLE[k] for k, v in st.filters.items()
                  if k in _ASSESSABLE and v is not None]
        n_unassessable = sum(
            1 for r in kept if checks and not all(c(r) for c in checks))

        # Redundancy with what comes after: pairs this stage removed that
        # no LATER stage would have removed anyway.
        #
        # Counted in PAIRS, like n_removed. Counting rows here while
        # n_removed counted pairs let unique exceed removed - CCR5's
        # burden stage reported 133 unique out of 124 removed - because a
        # drug reaching one gene through three factors is one pair and
        # three rows.
        later = stages[i + 1:]
        survives_later = [r for r in removed
                          if all(_passes(r, s2) for s2 in later)]
        unique = sum(len(v) for v in _pairs(survives_later).values())

        before, after = _pairs(current), _pairs(kept)
        lost = sorted(g for g in before if not after.get(g))

        results.append(StageResult(
            name=st.name, axis=st.axis, why=st.why,
            n_in=sum(len(v) for v in before.values()),
            n_remaining=sum(len(v) for v in after.values()),
            n_removed=(sum(len(v) for v in before.values())
                       - sum(len(v) for v in after.values())),
            n_not_assessable=n_unassessable,
            n_unique_kills=unique,
            genes_in=len(before), genes_out=len([g for g in after if after[g]]),
            genes_lost=lost,
            removed_examples=[f"{r['drug']} ({r.get('gene')} via "
                              f"{r.get('acts_on')})" for r in removed[:5]],
            requested=any(v is not None for v in st.filters.values()),
        ))
        current = kept
    return current, results


def default_stages(*, max_tier=4, min_grade=None, max_regulon=600,
                   min_tissue_breadth=0.0, max_burden="moderate",
                   require_chip=False, trustworthy_regulon=False,
                   min_paralog_specificity=0.0) -> list:
    return [
        Stage("Evidence floor", "evidence",
              {"max_tier": max_tier, "min_grade": min_grade,
               "exclude_conflicts": True,
               "require_chip": require_chip or None},
              "Removes routes the package itself does not believe. Runs first "
              "so the off-target counts below are not dominated by rows you "
              "were never going to use."),
        Stage("Not genome-wide", "regulon breadth",
              {"exclude_categories": ["chromatin-level"]},
              "Chromatin agents move thousands of promoters (specificity "
              "~0.05) and carry no letter grade, so they slip past every "
              "numeric gate on a technicality."),
        Stage("Narrow regulon", "regulon breadth",
              {"max_regulon": max_regulon,
               "trustworthy_regulon": trustworthy_regulon or None},
              "The factor must not run half the transcriptome. Gated here "
              "ONLY — axes.specificity already multiplies molecular "
              "selectivity by the regulon score, so filtering on both would "
              "count it twice."),
        Stage("Narrow TF tissue footprint", "tissue breadth",
              {"min_tissue_breadth": min_tissue_breadth or None},
              "How body-wide the factor is. A PROXY: it measures the "
              "transcription factor's expression, not the drug's "
              "distribution, and is withheld for ligand-gated receptors "
              "where mRNA presence does not imply activity."),
        Stage("Selective ligand achievable", "paralog family",
              {"min_paralog_specificity": min_paralog_specificity or None},
              "How many close relatives share the binding pocket. An "
              "orthosteric site is conserved across a receptor family, so "
              "THRA/THRB and RARA/RARG cross-react by construction rather "
              "than by any drug's failing. Drugs with a subtype-selective "
              "agent in their family are exempt — trifarotene exists "
              "because of this problem, not despite it. Left OFF unless you "
              "ask, because unlike the other axes this is a structural prior "
              "rather than a measurement: the test is a binding panel. The "
              "Family column is populated either way, so nothing is hidden "
              "by not filtering on it."),
        Stage("Tolerable clinical burden", "clinical burden",
              {"max_burden": max_burden},
              "Withdrawal and black-box flags, pleiotropic drug classes, and "
              "global-regulator targets. Catches the case every other filter "
              "misses: a clean, selective drug at a whole-body receptor."),
    ]


PANEL_PRESETS = {
    "narrow": {"max_tier": 4, "max_regulon": 600, "max_burden": "moderate",
               "min_paralog_specificity": 0.45},
    "narrow+clean": {"max_tier": 3, "min_grade": "B", "max_regulon": 300,
                     "max_burden": "low", "trustworthy_regulon": True,
                     "min_paralog_specificity": 0.6},
    "exploratory": {"max_tier": 5, "max_regulon": 2000, "max_burden": "high"},
}


# ------------------------------------------------------------- recurrence
# `catalogue` writes a literal em-dash into acts_on for routes with no
# protein target (ligand supplementation, for one). A truthiness check lets
# that through, and a target list of ["—"] then reads as "one shared
# regulator" when in fact there is no regulator at all.
_NO_TARGET = {"", "-", "—", "–", "?", "none", "unknown"}


def real_target(value) -> str | None:
    """The route's actual protein target, or None for a placeholder."""
    if not isinstance(value, str):
        return None
    v = value.strip()
    return None if v.lower() in _NO_TARGET else v


def recurrence(rows: list) -> dict:
    """Which candidate genes each drug shows up against.

    Counts GENES, not rows: one drug reaching one gene through three factors
    is recurrence 1. Keyed on the salt-free name so Levothyroxine and
    Levothyroxine Sodium are one drug, matching what the per-gene tables do.

    Diagnostic only. It is not in PREDICATES, no stage may name it, and this
    runs after the funnel has already returned.
    """
    out: dict = {}
    for r in rows:
        key = cong_mod.salt_free_name(r["drug"])
        e = out.setdefault(key, {"display": r["drug"], "genes": set(),
                                 "targets": set()})
        e["genes"].add(r.get("gene"))
        tgt = real_target(r.get("acts_on"))
        if tgt:
            e["targets"].add(tgt)
    return {k: {"display": v["display"], "genes": sorted(v["genes"]),
                "n_genes": len(v["genes"]), "targets": sorted(v["targets"])}
            for k, v in out.items()}


def annotate_recurrence(rows: list, rec: dict) -> None:
    """Stamp recurrence onto rows. Removes nothing, by construction."""
    for r in rows:
        e = rec.get(cong_mod.salt_free_name(r["drug"])) or {}
        r["n_genes"] = e.get("n_genes", 1)
        r["recurs_on"] = e.get("genes", [])
        r["recurs_targets"] = e.get("targets", [])


# ------------------------------------------------------------- the report
_BADGE = {"interpretable": "✓", "degraded": "⚠",
          "not interpretable": "⛔", "unknown": "?"}


def genetic_support(merged: dict) -> list:
    """Per-gene human genetic evidence, for the table printed BEFORE the funnel.

    Narrowing genes by how cleanly they can be drugged optimises for
    tractability, not therapeutic value — the easiest gene to drug is not
    necessarily one worth moving. This gates nothing; it sits above the funnel
    so tractability cannot quietly win the argument on its own.
    """
    by_gene = {sym: {"gene": sym, "summary": None, "n_phenotypes": 0, "top": []}
               for sym in (merged.get("panel", {}).get("per_gene") or {})}
    for block in merged.get("results", []):
        gene = (block.get("promoter") or {}).get("gene_symbol")
        gen = block.get("genetics") or {}
        row = by_gene.get(gene)
        if not row or not gen:
            continue
        phen = (gen.get("cis") or {}).get("phenotypes") or []
        if len(phen) > row["n_phenotypes"]:
            row["n_phenotypes"] = len(phen)
            row["summary"] = gen.get("summary")
            row["top"] = [p["phenotype"] for p in phen[:3]]
    return list(by_gene.values())


def cap_per_gene(rows: list, n: int | None) -> list:
    """At most `n` rows per gene, keeping the global order otherwise.

    `filter_rows` sorts by `(tier, -score)` across the whole panel, so a gene
    with more annotated pharmacology fills the table and a genuinely cleaner
    gene with three routes can fall off the bottom. For a tool whose job is
    comparing genes, that is the ranking answering the wrong question.
    """
    if not n or n <= 0:
        return rows
    seen: dict = {}
    out = []
    for r in rows:
        g = r.get("gene")
        if seen.get(g, 0) >= n:
            continue
        seen[g] = seen.get(g, 0) + 1
        out.append(r)
    return out


def markdown_report(panel: Panel, merged: dict, rows: list, stages: list,
                    rec: dict, criteria: dict | None = None) -> str:
    criteria = criteria or {}
    L: list[str] = []
    genes = [g.symbol for g in panel.included]
    goal = merged.get("query", {}).get("goal", "up")

    L.append(f"# Narrowing {len(genes)} candidate gene(s) — "
             f"{'up' if goal == 'up' else 'down'}regulation")
    L.append("")
    L.append("These genes are treated as **alternatives you are choosing "
             "between**, not a panel to move together. That inverts how a "
             "recurring drug reads: a compound hitting several of them is "
             "not corroborated, it is unselective.")
    L.append("")

    # -------------------------------------------- 1. genes in the comparison
    L.append("## Genes in the comparison")
    L.append("")
    L.append("| Gene | Evidence | Surviving pairs | Measured sources |")
    L.append("|---|---|---|---|")
    for g in panel.included:
        srcs = (g.result.get("evidence_sources") or {})
        n_left = len({(cong_mod.salt_free_name(r["drug"]), r.get("acts_on"))
                      for r in rows if r.get("gene") == g.symbol})
        L.append(f"| {g.symbol} | {_BADGE.get(g.level, '?')} {g.level} | "
                 f"{n_left} | {srcs.get('empirical_sources_ok', 0)} |")
    L.append("")
    if stages:
        L.append(f"*{stages[0].n_in} distinct (drug, target) pairs entered "
                 f"the funnel.*")
        L.append("")

    if panel.dropped:
        L.append("### Genes excluded before narrowing")
        L.append("")
        L.append("Never silently dropped — a gene removed here contributed "
                 "nothing to the comparison and you should know which.")
        L.append("")
        for g in panel.dropped:
            L.append(f"- **{g.symbol}** — {g.error or g.excluded}")
            for reason in g.verdict.get("reasons", [])[:3]:
                L.append(f"  - {reason}")
        L.append("")

    # --------------------------------------------------------- 2. parity
    par = (merged.get("panel") or {}).get("parity") or {}
    if par.get("warnings"):
        L.append("### ⚠ These genes were not run the same way")
        L.append("")
        L.append("A gene run without a source loses every stage depending on "
                 "it, for reasons unrelated to its biology.")
        L.append("")
        for w in par["warnings"]:
            L.append(f"- {w}")
        L.append("")

    # ------------------------------------------ 3. genetic support (no gate)
    gs = [r for r in genetic_support(merged) if r["n_phenotypes"]]
    if gs:
        L.append("## Is each gene worth moving at all?")
        L.append("")
        L.append("Human genetic support, printed **before** the funnel on "
                 "purpose. Everything below narrows on how cleanly a gene can "
                 "be drugged — that is tractability, not therapeutic value. "
                 "This table gates nothing.")
        L.append("")
        L.append("| Gene | Associated phenotypes | Top |")
        L.append("|---|---|---|")
        for r in sorted(gs, key=lambda x: -x["n_phenotypes"]):
            L.append(f"| {r['gene']} | {r['n_phenotypes']} | "
                     f"{', '.join(r['top']) or '—'} |")
        L.append("")

    # ----------------------------------------------------------- 4. funnel
    L.append("## The funnel")
    L.append("")
    L.append("Each stage applies one criterion. Counts are distinct "
             "(active moiety, target) pairs, not rows — with `--all-tss` a "
             "gene contributes several near-identical windows, and counting "
             "rows would let the most-annotated gene dominate.")
    L.append("")
    L.append("| Stage | Axis | In | Removed | Unique | Not assessable | Out |")
    L.append("|---|---|---|---|---|---|---|")
    for s in stages:
        L.append(f"| {s.name} | {s.axis} | {s.n_in} | {s.n_removed} | "
                 f"{s.n_unique_kills} | {s.n_not_assessable} | "
                 f"{s.n_remaining} |")
    L.append("")
    L.append("**Unique** = pairs this stage removed that no *later* stage "
             "would have removed anyway. Large Removed with small Unique "
             "means the stage is largely redundant with the next one — the "
             "three off-target axes are correlated, and this is where you see "
             "how much each actually adds.")
    L.append("")

    for s in stages:
        L.append(f"**{s.name}** — {s.why}")
        L.append("")
        if not s.requested:
            L.append("- *Not applied — no threshold was set for this axis, "
                     "so nothing was filtered on it. This is a choice, not a "
                     "missing measurement.*")
        elif s.did_nothing:
            L.append(f"- ⚠ **This stage did nothing.** It removed 0 and "
                     f"could not assess {s.n_not_assessable} of {s.n_in}. A "
                     f"zero here is missing data, not a clean result.")
        elif s.n_not_assessable:
            L.append(f"- {s.n_not_assessable} of {s.n_in} passed only because "
                     f"the datum was missing, not by clearing the bar.")
        if s.genes_lost:
            L.append(f"- **Genes eliminated here: {', '.join(s.genes_lost)}** "
                     f"— this is where those candidates ran out of routes.")
        if s.removed_examples:
            # Labelled, because the bullet above it is about rows that
            # SURVIVED unjudged and this one is about rows that were
            # judged and dropped. A bare "e.g." between them reads as
            # examples of the sentence before it, which is the opposite
            # set - and it has been read that way.
            L.append(f"- Removed here, for example: "
                     f"{'; '.join(s.removed_examples[:3])}")
        L.append("")

    # ------------------------------------------- 5. gene x stage survival
    if genes and stages:
        L.append("### Where each gene dropped out")
        L.append("")
        L.append("| Gene | " + " | ".join(s.name for s in stages) + " |")
        L.append("|---" * (len(stages) + 1) + "|")
        for g in genes:
            cells = ["**out**" if g in s.genes_lost else "·"
                     for s in stages]
            L.append("| " + " | ".join([g] + cells) + " |")
        L.append("")
        ctrl_sym = ((merged.get("panel") or {}).get("control_gene") or "")
        still = sorted({r.get("gene") for r in rows
                        if r.get("gene") != ctrl_sym})
        L.append(f"**Still standing: {', '.join(still) or 'none'}.**"
                 + (f" (control {ctrl_sym} excluded)" if ctrl_sym else ""))
        L.append("")

    # --------------------------------------------------- 6. surviving routes
    ctrl_sym = ((merged.get("panel") or {}).get("control_gene") or "")
    cand_rows = [r for r in rows if r.get("gene") != ctrl_sym]
    cand_rows = cap_per_gene(cand_rows, criteria.get("per_gene_top"))
    if cand_rows:
        L.append("## Surviving routes")
        L.append("")
        if ctrl_sym:
            L.append(f"*Candidates only — the control gene {ctrl_sym} is "
                     f"reported separately below.*")
            L.append("")
        L.append("| Gene | Tier | Grade | Drug | Acts on | Regulon | "
                 "TF tissues | Family | Burden | Other genes | Score |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for r in cand_rows[:criteria.get("top", 40)]:
            reg = r.get("_regulon_size")
            reg_s = "?" if reg is None else str(reg)
            if r.get("_regulon_understudied"):
                reg_s += " ⚠understudied"
            tb = r.get("_tissue_breadth")
            tb_s = "not assessed" if tb is None else f"{tb:.2f}"
            burden = r.get("_burden", "?")
            if not r.get("_burden_assessed"):
                burden += " (unassessed)"
            ng = r.get("n_genes", 1)
            other = "—" if ng <= 1 else f"⚠ {ng} of {len(genes)}"
            nrel = r.get("_paralog_relatives")
            if not r.get("_paralog_assessed"):
                fam_s = "—"
            elif r.get("_paralog_selective_agents"):
                fam_s = f"{len(nrel)} rel., selective agent exists"
            elif nrel:
                fam_s = f"⚠ {len(nrel)} rel. ({', '.join(nrel[:2])})"
            else:
                fam_s = "no paralogue"
            L.append(f"| {r.get('gene')} | {r['tier']} | "
                     f"{r.get('grade') or '—'} | {r['drug']} | "
                     f"{r.get('acts_on')} | {reg_s} | {tb_s} | {fam_s} | "
                     f"{burden} | {other} | {r['score']} |")
        L.append("")
        if len(cand_rows) >= 2:
            top = cand_rows[0]
            n_cmp = stages[0].n_in if stages else len(rows)
            L.append(f"**About the top row.** It is the maximum over ~{n_cmp} "
                     f"candidate pairs across {len(genes)} gene(s). You are "
                     f"looking at it *because* its predicted score came out "
                     f"high, so the same noise that surfaced it has already "
                     f"inflated it — and the more genes you compare, the more "
                     f"that is true.")
            L.append("")
            # The runner-up has to be at the SAME tier. Rows are ordered
            # by tier first, so taking cand_rows[1] reached into a weaker
            # tier and subtracted scores that do not share a scale - it
            # reported a NEGATIVE margin for CCR5 and called it a tie,
            # when being the only route at the strongest tier is the
            # opposite of a tie.
            peers = [r for r in cand_rows if r["tier"] == top["tier"]]
            if len(peers) >= 2:
                runner = peers[1]
                margin = round(top["score"] - runner["score"], 2)
                if margin < 1.5:
                    L.append(f"Its margin over {runner['drug']}, the next "
                             f"route at tier {top['tier']}, is {margin} "
                             f"points, which is **not separable from "
                             f"noise**. Treat the top few as a tied set, "
                             f"not a ranking.")
                else:
                    L.append(f"Margin over {runner['drug']}, the next route "
                             f"at tier {top['tier']}: {margin} points.")
            else:
                L.append(f"It is the **only route at tier {top['tier']}**, "
                         f"and everything below it rests on weaker "
                         f"evidence. There is no like-for-like runner-up to "
                         f"measure a margin against: scores are comparable "
                         f"within a tier, not across them, so the gap to "
                         f"the next row is not a quantity worth reporting.")
            L.append("")
    else:
        L.append("## Nothing survived")
        L.append("")
        if stages:
            binding = max(stages, key=lambda s: s.n_removed)
            L.append(f"The binding constraint is **{binding.name}**, which "
                     f"removed {binding.n_removed} of {binding.n_in}. Loosen "
                     f"that one first.")
            L.append("")
        L.append("An empty result is informative: it says no route clears "
                 "your bar for any of these genes on current evidence.")
        L.append("")

    # ------------------------------------------------- 6a. the control gene
    control = (merged.get("panel") or {}).get("control_gene")
    if control:
        ctrl = control.upper()
        surviving = _pairs(rows)
        n_ctrl = len(surviving.get(ctrl) or ())
        others = {g: len(v) for g, v in surviving.items() if g != ctrl}
        median = (sorted(others.values())[len(others) // 2]
                  if others else 0)
        L.append("### Negative control")
        L.append("")
        L.append(f"**{ctrl}** was run through the identical funnel as a "
                 f"control. It is not a candidate — it is the check that the "
                 f"funnel discriminates at all. Without one, an empty "
                 f"shortlist and a misconfigured filter look the same.")
        L.append("")
        L.append(f"- {ctrl}: **{n_ctrl}** surviving route(s)")
        L.append(f"- your candidates: median **{median}**"
                 + (f" (range {min(others.values())}–{max(others.values())})"
                    if others else ""))
        L.append("")
        if n_ctrl >= median and median > 0:
            L.append("> **⚠ The control survived as well as your candidates.** "
                     "The funnel is not discriminating — it is passing routes "
                     "on properties every gene shares, so the surviving list "
                     "above is not evidence about your genes specifically. "
                     "Tighten the stages that removed the least, or treat this "
                     "run as uninformative.")
        elif n_ctrl == 0:
            L.append("> ✓ The control was eliminated entirely while "
                     "candidates survived. The funnel is discriminating.")
        else:
            L.append(f"> The control retained fewer routes than the median "
                     f"candidate, so the funnel is discriminating — but it "
                     f"did not eliminate the control outright, so read the "
                     f"margin as modest rather than clean.")
        L.append("")

    # ------------------------------------ 6b. relatives among the targets
    pairs = par_mod.cross_reactive_pairs(
        [r.get("acts_on") for r in rows if r.get("acts_on")])
    if pairs:
        L.append("### Some of these targets are relatives")
        L.append("")
        L.append("Two members of one receptor family scoring well is **one "
                 "bet, not two**. The orthosteric pocket is under selection "
                 "to bind the same endogenous ligand across the family, so "
                 "it is conserved — a ligand at one will usually engage the "
                 "others, and that is a property of the family rather than a "
                 "defect in any particular drug.")
        L.append("")
        L.append("| Family | Members here | Shared pocket binds |")
        L.append("|---|---|---|")
        for p in pairs:
            L.append(f"| {p['family']} | {', '.join(p['members'])} | "
                     f"{p['ligand']} |")
        L.append("")
        for p in pairs:
            if p["selective_agents"]:
                L.append(f"- **{p['family']}** — subtype selectivity HAS been "
                         f"achieved: {'; '.join(p['selective_agents'])}. A "
                         f"large family is a risk, not a verdict.")
        L.append("")

    # ------------------------------------------------------ 7. recurrence
    L.append("## Cross-gene recurrence")
    L.append("")
    L.append("**A diagnostic, not a filter — it removed nothing above.** "
             "Because your genes are alternatives, a drug hitting several of "
             "them is *less* attractive, not more: it moves a lot of things "
             "and one of them happens to be whichever gene you pick. "
             "`1 gene` is the good case.")
    L.append("")
    multi = sorted([e for e in rec.values() if e["n_genes"] > 1],
                   key=lambda e: -e["n_genes"])
    if multi:
        L.append("| Drug | Genes hit | Via | Reading |")
        L.append("|---|---|---|---|")
        for e in multi[:20]:
            n_t = len(e["targets"])
            if n_t == 1:
                reading = "your genes share this regulator"
            elif n_t > 1:
                reading = "promiscuous — a different target per gene"
            else:
                # e.g. ligand supplementation, which acts on no protein here
                reading = "no protein target recorded — cannot tell which"
            L.append(f"| {e['display']} | {e['n_genes']}: "
                     f"{', '.join(e['genes'])} | "
                     f"{', '.join(e['targets'][:3]) or '—'} | "
                     f"{reading} |")
        L.append("")
        L.append("*Recurrence through the **same** factor means your genes "
                 "share a regulator — real biology, and a single point of "
                 "failure. Through **different** factors it means the drug is "
                 "unselective. Same column, opposite diagnoses.*")
    else:
        L.append("No drug appeared against more than one gene. On the "
                 "alternatives reading, that is the good outcome.")
    L.append("")

    # --------------------------------------- 8. what this cannot tell you
    L.append("## What this cannot tell you")
    L.append("")
    L.append("- **The drug's own tissue distribution.** Nothing here maps a "
             "compound to the tissues it reaches. The tissue axis is the "
             "*transcription factor's* expression breadth used as a proxy, "
             "and it is withheld entirely for ligand-gated receptors where "
             "mRNA presence does not imply activity.")
    L.append("- **Whether a low regulon count means specific or unstudied.** "
             "TRRUST is literature-curated, so a small count can mean nobody "
             "looked. Rows marked ⚠understudied are flagged for that.")
    L.append("- **Whether the drug-class table is complete.** It is "
             "hand-maintained and short. A row with no burden flag is "
             "*unassessed*, not cleared.")
    L.append("- **Whether these genes are worth moving.** That is the "
             "genetics table above, and it is deliberately not in the funnel.")
    L.append("")
    return "\n".join(L)


def write_outputs(panel: Panel, merged: dict, rows: list, stages: list,
                  rec: dict, outdir, criteria: dict | None = None) -> list:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []

    p = outdir / "panel.md"
    p.write_text(markdown_report(panel, merged, rows, stages, rec, criteria),
                 encoding="utf-8")
    written.append(p)

    payload = dict(merged)
    payload["funnel"] = [s.to_dict() for s in stages]
    payload["recurrence"] = rec
    payload["surviving"] = rows
    payload["criteria"] = criteria or {}
    p = outdir / "panel.json"
    p.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    written.append(p)

    cols = ["gene", "tier", "grade", "drug", "acts_on", "category",
            "_regulon_size", "_regulon_understudied", "_tissue_breadth",
            "_tissue_basis", "_burden", "_burden_assessed", "n_genes", "score"]
    lines = ["\t".join(c.lstrip("_") for c in cols)]
    for r in rows:
        lines.append("\t".join("" if r.get(c) is None else str(r.get(c))
                               for c in cols))
    p = outdir / "panel.tsv"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    written.append(p)
    return written


# -------------------------------------------------------------------- CLI
def narrow(panel: Panel, *, regulons=None, control: str | None = None,
           **criteria):
    """Merge, annotate, run the funnel, then annotate recurrence.

    Recurrence is applied only after `run_funnel` has returned, so it is
    structurally incapable of affecting what survived.
    """
    merged = merge_results(panel)
    merged.setdefault("panel", {})["parity"] = parity(panel)

    rows = sl.enrich_rows(merged)
    annotate_axes(rows, merged, regulons)

    stage_kwargs = {k: v for k, v in criteria.items()
                    if k in ("max_tier", "min_grade", "max_regulon",
                             "min_tissue_breadth", "max_burden",
                             "require_chip", "trustworthy_regulon",
                             "min_paralog_specificity")}
    stages = default_stages(**stage_kwargs)
    kept, results = run_funnel(rows, stages)
    kept.sort(key=lambda r: (r["tier"], -r["score"]))

    # Recurrence is a statement about the CANDIDATES. Counting the control
    # gene in it would inflate every drug's gene count by one and make the
    # control look like corroboration, which is the opposite of its job.
    ctrl = (control or "").upper()
    scored = [r for r in kept if not ctrl or r.get("gene") != ctrl]
    rec = recurrence(scored)
    annotate_recurrence(kept, rec)
    return merged, kept, results, rec


def build_parser():
    import argparse
    from .cli import add_analysis_args

    p = argparse.ArgumentParser(
        prog="promoter-tf-drugs.panel",
        description="Compare several CANDIDATE genes and narrow to the "
                    "route with the least collateral effect. The genes are "
                    "treated as alternatives, so a drug recurring across "
                    "them is a warning, not corroboration.")
    p.add_argument("--genes", help="comma-separated gene symbols")
    p.add_argument("--gene-file", help="one gene symbol per line")
    p.add_argument("--control-gene",
                   help="run this gene through the identical funnel as a "
                        "negative control. A funnel that treats a "
                        "housekeeping gene the same as your candidates is "
                        "not discriminating, and without one an empty "
                        "shortlist and a broken filter look identical")
    p.add_argument("--preset", choices=sorted(PANEL_PRESETS))
    p.add_argument("--max-tier", type=int)
    p.add_argument("--min-grade", choices=list("ABCDF"))
    p.add_argument("--max-regulon", type=int,
                   help="drop factors controlling more genes than this")
    p.add_argument("--trustworthy-regulon", action="store_true",
                   help="require a real regulon count, not a bundled tier or "
                        "a handful of curated targets that reflect study "
                        "effort rather than biology")
    p.add_argument("--min-tissue-breadth", type=float,
                   help="0..1, higher = narrower. A PROXY via the factor's "
                        "expression; withheld for nuclear receptors")
    p.add_argument("--max-burden", choices=list(sys_mod.BURDEN_LEVELS))
    p.add_argument("--min-paralog-specificity", type=float,
                   help="0..1, higher = fewer relatives share the binding "
                        "pocket, so a subtype-selective ligand is more "
                        "achievable. Drugs whose family already has a "
                        "selective agent are exempt")
    p.add_argument("--require-chip", action="store_true")
    p.add_argument("--require-interpretable", action="store_true",
                   help="also exclude genes whose evidence verdict is merely "
                        "'degraded', not just 'not interpretable'")
    p.add_argument("--renarrow", metavar="PANEL_JSON",
                   help="re-filter a previous panel.json with new cutoffs. "
                        "No network, no re-run — the point being that tuning "
                        "a threshold should take a second, not the tens of "
                        "minutes per gene a fresh run costs")
    p.add_argument("--per-gene-top", type=int,
                   help="at most this many rows per gene in the surviving "
                        "table. filter_rows sorts by (tier, -score) across "
                        "the whole panel, so without this a gene with more "
                        "annotated pharmacology fills the table and a "
                        "genuinely cleaner gene with three routes falls off "
                        "the bottom — the ranking answering the wrong "
                        "question for a tool that compares genes")
    p.add_argument("--top", type=int, default=40)
    p.add_argument("--out", default="results/panel")

    # Shared analysis flags, declared once in cli.add_analysis_args rather
    # than borrowed out of the other parser's private _actions.
    add_analysis_args(p)
    return p





def _check_flag_combinations(parser, args) -> None:
    """Reject flag combinations that cannot work, before anything runs.

    A full run is a quarter of an hour of scanning and lookups. Learning
    at the end that a requested source could never have produced
    anything wastes all of it, and the condition is knowable before the
    first request goes out.

    This is the CLI's check, not the pipeline's: `run()` still records a
    failed source and carries on, because a library caller may want the
    rest of the run regardless.
    """
    if getattr(args, "tf_activity", False) and not getattr(
            args, "activity_tissue", None):
        parser.error(
            "--tf-activity needs --activity-tissue: it infers which "
            "factors are active from their targets' expression in one "
            "tissue, and without a tissue there is no expression profile "
            "to infer from. Name a tissue as the expression table spells "
            "it, e.g. --activity-tissue liver.")
    if getattr(args, "activity_permutations", 0) and not getattr(
            args, "tf_activity", False):
        parser.error(
            "--activity-permutations only applies to --tf-activity, "
            "which was not requested. It would do nothing here.")


def main(argv=None) -> int:
    import sys
    parser = build_parser()
    args = parser.parse_args(argv)
    _check_flag_combinations(parser, args)
    # Checkpoints go to stderr. The panel report goes to stdout, so a
    # piped run is unaffected.
    if not getattr(args, "quiet", False):
        _progress.enable()
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    logging.basicConfig(
        level=[logging.WARNING, logging.INFO,
               logging.DEBUG][min(getattr(args, "verbose", 0), 2)],
        format="%(levelname)s %(name)s: %(message)s")

    # ---- re-narrow an existing panel, with no network and no re-run
    if args.renarrow:
        saved = json.loads(Path(args.renarrow).read_text(encoding="utf-8"))
        panel = Panel()
        for sym, info in (saved.get("panel", {}).get("per_gene") or {}).items():
            pg = PanelGene(symbol=sym)
            # Reconstruct just enough for the funnel and the report: the
            # blocks for this gene, plus its own query and verdict.
            pg.result = {
                "results": [b for b in saved.get("results", [])
                            if (b.get("promoter") or {}).get("gene_symbol") == sym],
                "query": info.get("query") or {},
                "evidence_verdict": info.get("evidence_verdict") or {},
                "evidence_sources": info.get("evidence_sources") or {},
                "direct_therapies": info.get("direct_therapies") or [],
            }
            panel.genes.append(pg)
        if not panel.included:
            print(f"{args.renarrow} holds no usable gene blocks.",
                  file=sys.stderr)
            return 1

        criteria = dict(PANEL_PRESETS.get(args.preset, {})) if args.preset else {}
        for key in ("max_tier", "min_grade", "max_regulon",
                    "min_tissue_breadth", "max_burden", "require_chip",
                    "trustworthy_regulon", "min_paralog_specificity"):
            v = getattr(args, key, None)
            if v not in (None, False):
                criteria[key] = v
        criteria["top"] = args.top
        control = (saved.get("panel", {}) or {}).get("control_gene")
        merged, kept, stages, rec = narrow(panel, control=control, **criteria)
        if control:
            merged.setdefault("panel", {})["control_gene"] = control
        written = write_outputs(panel, merged, kept, stages, rec, args.out,
                                criteria)
        print(markdown_report(panel, merged, kept, stages, rec, criteria))
        print("\nRe-narrowed with no network. Wrote: "
              + ", ".join(str(w) for w in written), file=sys.stderr)
        return 0

    genes = []
    if args.genes:
        genes += [g.strip().upper() for g in args.genes.split(",") if g.strip()]
    if args.gene_file:
        genes += [ln.strip().upper()
                  for ln in Path(args.gene_file).read_text().splitlines()
                  if ln.strip() and not ln.startswith("#")]
    if args.control_gene:
        genes.append(args.control_gene.strip().upper())
    if not genes:
        print("No genes given. Use --genes A,B,C or --gene-file.",
              file=sys.stderr)
        return 2

    # The same derivation the single-gene CLI uses. A hand-written list
    # here forwarded twenty-one of the fifty arguments run() takes, so
    # --genetics, --lincs, --upstream-network, --tf-activity and
    # --docking-results were accepted and silently dropped.
    from .cli import run_kwargs_from_args
    run_kwargs = run_kwargs_from_args(args)

    panel = run_panel(genes, out_dir=args.out,
                      require_interpretable=args.require_interpretable,
                      **run_kwargs)
    if not panel.included:
        print("Every gene was excluded or failed; nothing to compare.",
              file=sys.stderr)
        for g in panel.dropped:
            print(f"  {g.symbol}: {g.error or g.excluded}", file=sys.stderr)
        return 1

    criteria = dict(PANEL_PRESETS.get(args.preset, {})) if args.preset else {}
    for key in ("max_tier", "min_grade", "max_regulon", "min_tissue_breadth",
                "max_burden", "require_chip", "trustworthy_regulon",
                "min_paralog_specificity"):
        v = getattr(args, key, None)
        if v not in (None, False):
            criteria[key] = v
    criteria["top"] = args.top
    criteria["per_gene_top"] = args.per_gene_top

    regulons = (regulon_mod.load_regulons(args.regulon_file)
                if getattr(args, "regulon_file", None) else None)
    control = args.control_gene.strip().upper() if args.control_gene else None
    merged, kept, stages, rec = narrow(panel, regulons=regulons,
                                       control=control, **criteria)
    if control:
        merged.setdefault("panel", {})["control_gene"] = control

    written = write_outputs(panel, merged, kept, stages, rec, args.out,
                            criteria)
    print(markdown_report(panel, merged, kept, stages, rec, criteria))
    print("\nWrote: " + ", ".join(str(w) for w in written), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
