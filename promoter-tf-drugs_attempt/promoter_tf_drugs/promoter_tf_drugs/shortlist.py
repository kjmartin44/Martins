"""Narrow an existing results.json to a shortlist, without re-running anything.

A full run takes tens of minutes. Deciding what to actually order should take a
second. This reads the results.json a run already produced and applies filters
to it, so you can tighten and loosen criteria interactively.

    python -m promoter_tf_drugs.shortlist results/PAH_up/results.json --preset strict

The presets are opinionated on purpose. Each encodes a different question:

  strict      "what would I actually order next week?" Evidence tier 1-3,
              grade A/B, no withdrawn or black-box drugs, no direction
              conflicts, TF must have measured binding at the promoter.
              Usually returns single digits, and often zero — which is a real
              answer about your gene, not a failure of the filter.

  mechanism   "what do I understand well enough to design an experiment
              around?" Requires a named mechanism and a specific target;
              drops empirical hits with unknown mechanism and genome-wide
              chromatin drugs.

  empirical   "what has actually been measured?" Tier 1-2 only: approved for
              this gene, or observed moving it in a cell assay. Ignores every
              prediction in the package, including the good ones.

  specific    "what will not move half the transcriptome?" Hard cut on
              regulon breadth and protein selectivity; drops global regulators
              and chromatin drugs regardless of how well they score.

  novel       "what has nobody tried on this gene?" Excludes anything with a
              reported link, keeps tiers 3-5. The opposite risk profile from
              `strict` and should be read that way.

Filters compose, so a preset is a starting point you then adjust.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PRESETS = {
    "strict": {
        "max_tier": 3, "min_grade": "B", "require_chip": True,
        "exclude_withdrawn": True, "exclude_black_box": True,
        "exclude_conflicts": True, "max_regulon": 2000,
        # Genome-wide agents have no letter grade, so without this they slip
        # through the grade filter on a technicality despite specificity ~0.05.
        "exclude_categories": ["chromatin-level"], "min_specificity": 0.3,
    },
    "mechanism": {
        "max_tier": 4, "require_mechanism": True,
        "exclude_categories": ["chromatin-level", "empirical signature"],
        "exclude_withdrawn": True,
    },
    "empirical": {"max_tier": 2},
    "specific": {
        "max_tier": 4, "max_regulon": 600, "min_specificity": 0.45,
        "exclude_categories": ["chromatin-level"], "max_steps": 1,
    },
    "novel": {"max_tier": 5, "min_tier": 3, "only_novel": True},
}

GRADE_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "F": 4}


def _grade_ok(grade: str | None, minimum: str | None) -> bool:
    """Missing grade passes.

    Letter grades only exist for TF-mediated routes that went through
    drug_scoring. A direct siRNA or an approved gene therapy has no grade, and
    dropping a tier-1 clinical agent for lacking a field that does not apply to
    it is exactly backwards.
    """
    if not minimum or grade is None:
        return True
    return GRADE_ORDER.get(grade, 4) <= GRADE_ORDER.get(minimum, 4)


def enrich_rows(result: dict) -> list[dict]:
    """Every route from every promoter block, with its supporting evidence
    stamped on. No filtering happens here.

    Split out from `filter_routes` so the same rows can be fed to a staged
    funnel that applies one criterion at a time and counts what each removed.
    The underscore-prefixed keys exist so every predicate below can be a pure
    `(row, value) -> bool` rather than reaching back into the block.
    """
    out: list[dict] = []

    for block in result.get("results", []):
        # Index the supporting evidence so route rows can be enriched.
        grades = {g["drug_name"].lower(): g
                  for g in (block.get("drug_grades") or [])}
        tf_chip = {t["symbol"]: t.get("chip_peaks", 0)
                   for t in (block.get("tfs") or [])}
        linked = {c["name"].lower()
                  for c in ((block.get("improvements") or {})
                            .get("all_candidates") or [])
                  if c.get("previously_linked")}
        celltype = ((block.get("cell_type") or {}).get("tfs") or {})
        verdict = (block.get("evidence_verdict") or {}).get("level")

        for r in (block.get("all_routes") or {}).get("routes", []):
            g = grades.get(r["drug"].lower(), {})
            tf = r.get("acts_on")

            row = dict(r)
            row["grade"] = g.get("grade")
            row["gene"] = block.get("promoter", {}).get("gene_symbol")
            row["previously_linked"] = r["drug"].lower() in linked
            row["chip_peaks"] = tf_chip.get(tf)

            row["_notes"] = " ".join(r.get("notes", []) + g.get("notes", []))
            row["_regulon_size"] = (g.get("evidence") or {}).get("regulon_size")
            row["_gene_verdict"] = verdict
            # `require_chip` has only ever applied to TFs the block actually
            # listed. Preserving that as an explicit flag matters: without it,
            # a route whose factor is missing from the TF table would start
            # being dropped, which is a silent behaviour change.
            row["_chip_known"] = tf in tf_chip
            # Absent, NOT zero, when the expression lookup failed — see
            # delivery.expression_overlap's early return. Zero would read as
            # "expressed nowhere" and sail through a breadth filter.
            row["_tf_tissues"] = (celltype.get(tf) or {}).get("n_tf")
            out.append(row)

    return _dedupe_across_windows(out)


def _dedupe_across_windows(rows: list[dict]) -> list[dict]:
    """One row per distinct route per gene, not per promoter window.

    With `--all-tss` a gene contributes one block per TSS, and those windows
    can overlap heavily — CCR5's two shared 4,275 of 6,000 bp. Gene-level
    routes are folded into every block by `catalogue.from_direct`, so an
    approved drug against the gene product appeared once per window:
    Maraviroc twice in a two-row shortlist, with its caveat printed twice.

    Routes are keyed on (gene, active moiety, target, category) — the salt
    form collapses too, for the same reason it does in the graded tables.
    The best-scoring instance survives and records how many windows produced
    it, because that count is real information about the scan rather than
    something to hide.
    """
    from .congeners import salt_free_name

    best: dict = {}
    order: list = []
    for r in rows:
        key = (r.get("gene"), salt_free_name(r["drug"]), r.get("acts_on"),
               r.get("category"))
        prev = best.get(key)
        if prev is None:
            r["_n_windows"] = 1
            best[key] = r
            order.append(key)
        else:
            prev["_n_windows"] = prev.get("_n_windows", 1) + 1
            # Keep the stronger evidence, not whichever window came first.
            if (r["tier"], -r["score"]) < (prev["tier"], -prev["score"]):
                r["_n_windows"] = prev["_n_windows"]
                best[key] = r
    return [best[k] for k in order]


# --------------------------------------------------------------- predicates
# One per filter name. Keys are exactly the names PRESETS and the CLI already
# use, so neither had to change. Each returns True to KEEP the row.
def _p_max_tier(row, v):            return row["tier"] <= v
def _p_min_tier(row, v):            return row["tier"] >= v
def _p_min_grade(row, v):           return _grade_ok(row.get("grade"), v)
def _p_max_steps(row, v):           return row["directness"] <= v
def _p_min_specificity(row, v):     return row["axes"]["specificity"] >= v
def _p_exclude_categories(row, v):  return row["category"] not in (v or [])
def _p_require_approved(row, v):    return not v or bool(row.get("approved"))
def _p_require_mechanism(row, v):   return not v or bool(row.get("mechanism"))
def _p_only_novel(row, v):          return not v or not row["previously_linked"]
def _p_only_known(row, v):          return not v or row["previously_linked"]
def _p_exclude_withdrawn(row, v):
    return not v or "WITHDRAWN" not in row["_notes"].upper()
def _p_exclude_black_box(row, v):   return not v or "black-box" not in row["_notes"]
def _p_exclude_conflicts(row, v):   return not v or "CONFLICT" not in row["_notes"]


def _p_max_regulon(row, v):
    """Missing regulon size passes — see `_grade_ok` for the same principle."""
    size = row.get("_regulon_size")
    return size is None or size <= v


def _p_require_chip(row, v):
    if not v:
        return True
    # Only applies to TF-mediated routes; a siRNA has no TF.
    return not (row["directness"] >= 1 and row.get("_chip_known")
                and not row.get("chip_peaks"))


PREDICATES = {
    "max_tier": _p_max_tier,
    "min_tier": _p_min_tier,
    "min_grade": _p_min_grade,
    "max_steps": _p_max_steps,
    "min_specificity": _p_min_specificity,
    "exclude_categories": _p_exclude_categories,
    "require_approved": _p_require_approved,
    "require_mechanism": _p_require_mechanism,
    "require_chip": _p_require_chip,
    "only_novel": _p_only_novel,
    "only_known": _p_only_known,
    "exclude_withdrawn": _p_exclude_withdrawn,
    "exclude_black_box": _p_exclude_black_box,
    "exclude_conflicts": _p_exclude_conflicts,
    "max_regulon": _p_max_regulon,
}


def validate_filters(f: dict) -> None:
    """Reject unknown filter names before any row is examined.

    Checked up front rather than per row, for two reasons. An empty candidate
    list would otherwise never reach the lookup, so a typo'd criterion would
    pass silently on exactly the runs where you can least afford it. And a
    failure that depends on whether any data happened to arrive is not a
    failure anyone can debug.

    A silently-ignored filter quietly widens a shortlist that reads as though
    it had been narrowed, and it is what structurally prevents a
    diagnostic-only field (cross-gene recurrence) being used as a gate.
    """
    unknown = sorted(set(f) - set(PREDICATES))
    if unknown:
        raise KeyError(
            f"unknown filter(s) {', '.join(repr(u) for u in unknown)}; "
            f"known filters are {', '.join(sorted(PREDICATES))}")


def keep(row: dict, **f) -> bool:
    """Does this row survive every named filter?"""
    return all(PREDICATES[name](row, value)
               for name, value in f.items() if value is not None)


def filter_rows(rows: list[dict], **f) -> list[dict]:
    """Apply named filters to already-enriched rows, then rank."""
    validate_filters(f)
    out = [r for r in rows if keep(r, **f)]
    out.sort(key=lambda r: (r["tier"], -r["score"]))
    return out


def filter_routes(result: dict, **f) -> list[dict]:
    """Apply filters across every promoter block in a results.json."""
    return filter_rows(enrich_rows(result), **f)


def render(rows: list[dict], criteria: dict, total: int) -> str:
    L = ["# Shortlist", ""]
    L.append(f"{len(rows)} of {total} routes kept.")
    L.append("")
    L.append("Criteria: " + ", ".join(f"{k}={v}" for k, v in criteria.items()))
    L.append("")
    if not rows:
        L.append("**Nothing passed.** For `strict` that is a common and "
                 "informative result: it says no approved drug clears the bar "
                 "for this gene on current evidence. Loosen one criterion at a "
                 "time to see which one is binding — usually `require_chip` or "
                 "`min_grade`.")
        return "\n".join(L)

    L.append("| Tier | Grade | Drug | Category | Acts on | Steps | Spec | Score | New? |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        L.append(f"| {r['tier']} | {r.get('grade') or '—'} | {r['drug']} | "
                 f"{r['category']} | {r['acts_on']} | {r['directness']} | "
                 f"{r['axes']['specificity']} | {r['score']} | "
                 f"{'' if r['previously_linked'] else 'NEW'} |")
    L.append("")
    flagged = [r for r in rows if r.get("notes")]
    if flagged:
        L.append("Caveats on these:")
        L.append("")
        for r in flagged[:12]:
            L.append(f"- **{r['drug']}** — {r['notes'][0][:110]}")
    return "\n".join(L)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="promoter-tf-drugs.shortlist",
        description="Filter an existing results.json into an actionable "
                    "shortlist. No network, no re-run.")
    p.add_argument("results", help="path to results.json from a previous run")
    p.add_argument("--preset", choices=sorted(PRESETS),
                   help="starting point; individual flags override it")
    p.add_argument("--max-tier", type=int, help="1 clinical … 5 hypothesis")
    p.add_argument("--min-tier", type=int)
    p.add_argument("--min-grade", choices=list("ABCDF"))
    p.add_argument("--max-steps", type=int,
                   help="0 = acts on the gene, 1 = on a TF, 2+ = upstream")
    p.add_argument("--min-specificity", type=float)
    p.add_argument("--max-regulon", type=int,
                   help="drop TFs controlling more genes than this")
    p.add_argument("--require-chip", action="store_true",
                   help="TF must have measured binding at this promoter")
    p.add_argument("--require-approved", action="store_true")
    p.add_argument("--require-mechanism", action="store_true")
    p.add_argument("--exclude-withdrawn", action="store_true")
    p.add_argument("--exclude-black-box", action="store_true")
    p.add_argument("--exclude-conflicts", action="store_true",
                   help="drop rows where the structural model disagrees with "
                        "the database annotation")
    p.add_argument("--only-novel", action="store_true",
                   help="only drugs with no reported link to this gene")
    p.add_argument("--only-known", action="store_true")
    p.add_argument("--top", type=int, default=40)
    p.add_argument("--out", help="write the shortlist here as markdown")
    args = p.parse_args(argv)

    data = json.loads(Path(args.results).read_text(encoding="utf-8"))

    criteria = dict(PRESETS.get(args.preset, {})) if args.preset else {}
    for key in ("max_tier", "min_tier", "min_grade", "max_steps",
                "min_specificity", "max_regulon", "require_chip",
                "require_approved", "require_mechanism", "exclude_withdrawn",
                "exclude_black_box", "exclude_conflicts", "only_novel",
                "only_known"):
        v = getattr(args, key, None)
        if v not in (None, False):
            criteria[key] = v

    # Count distinct routes, not rows. Summing per-block route lists
    # double-counts every gene-level route once per TSS window.
    total = len(enrich_rows(data))
    rows = filter_routes(data, **criteria)[:args.top]
    text = render(rows, criteria, total)

    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
