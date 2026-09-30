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


def filter_routes(result: dict, **f) -> list[dict]:
    """Apply filters across every promoter block in a results.json."""
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

        for r in (block.get("all_routes") or {}).get("routes", []):
            g = grades.get(r["drug"].lower(), {})
            notes = " ".join(r.get("notes", []) + g.get("notes", []))
            a = r["axes"]

            if r["tier"] > f.get("max_tier", 5):
                continue
            if r["tier"] < f.get("min_tier", 0):
                continue
            if not _grade_ok(g.get("grade"), f.get("min_grade")):
                continue
            if r["category"] in (f.get("exclude_categories") or []):
                continue
            if r["directness"] > f.get("max_steps", 99):
                continue
            if a["specificity"] < f.get("min_specificity", 0.0):
                continue
            if f.get("require_approved") and not r.get("approved"):
                continue
            if f.get("exclude_withdrawn") and "WITHDRAWN" in notes.upper():
                continue
            if f.get("exclude_black_box") and "black-box" in notes:
                continue
            if f.get("exclude_conflicts") and "CONFLICT" in notes:
                continue
            if f.get("require_mechanism") and not r.get("mechanism"):
                continue
            if f.get("only_novel") and r["drug"].lower() in linked:
                continue
            if f.get("only_known") and r["drug"].lower() not in linked:
                continue

            max_reg = f.get("max_regulon")
            if max_reg is not None:
                ev = (g.get("evidence") or {})
                size = ev.get("regulon_size")
                if size is not None and size > max_reg:
                    continue

            if f.get("require_chip"):
                tf = r.get("acts_on")
                # Only applies to TF-mediated routes; a siRNA has no TF.
                if r["directness"] >= 1 and tf in tf_chip and not tf_chip[tf]:
                    continue

            row = dict(r)
            row["grade"] = g.get("grade")
            row["gene"] = block.get("promoter", {}).get("gene_symbol")
            row["previously_linked"] = r["drug"].lower() in linked
            row["chip_peaks"] = tf_chip.get(r.get("acts_on"))
            out.append(row)

    out.sort(key=lambda r: (r["tier"], -r["score"]))
    return out


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

    total = sum(len((b.get("all_routes") or {}).get("routes", []))
                for b in data.get("results", []))
    rows = filter_routes(data, **criteria)[:args.top]
    text = render(rows, criteria, total)

    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
