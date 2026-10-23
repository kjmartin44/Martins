"""Every compound worth putting on a plate, with nothing ranked.

Why this is separate from the report
------------------------------------
The report answers "what does the evidence say". This answers "what goes
in the wells", and the two want opposite things from the same routes.

The report folds a chemical family into one row, because five retinoids
making the same argument is one argument. A plate wants all five: they
differ in potency, selectivity and how much of each reaches a cell, and
the cost of an extra well is a well. The report demotes a compound the
literature records for every gene, because fame inflates the evidence.
A plate should not care: fame is a property of the literature, not of
what the compound does to your cells, and a screen measures the latter
directly. The report ranks. A screen of this size does not need a rank -
tested against CTD, the ordering these routes produce is no better than
chance once fame is held equal, so a plate ordered by score is a plate
in an arbitrary order with a misleading column attached.

So this takes the same routes and applies none of that: every distinct
compound, once, with what found it and what is known about it.

What it still leaves out
------------------------
Compounds excluded on safety grounds. Those are listed separately rather
than dropped, with the reason, because the reason is clinical and a
screen is not: a spindle poison is a perfectly good mechanistic control
for a chromosome-instability gene even though nobody would prescribe it.
Deciding to add one back is a judgement this module should inform and
not make.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)

# Which of the six steps a route's evidence actually places it at.
#
# Most routes do not place themselves anywhere. A compound recorded by
# CTD as raising a gene, or whose LINCS signature contains it, is a
# measurement of the OUTCOME with no mechanism attached - transcription,
# stability and degradation are indistinguishable in that record. Naming
# a step for those would invent the one thing the data does not say, so
# they are grouped as what they are.
STEP_BY_SOURCE = {
    "network": "transcription",
    "regulators": "transcription",
    "abundance": "transcription",
    "tfactivity": "transcription",
    "drug_scoring": "transcription",
    "screening": "transcription",
    "docking": "transcription",
    "catalogue": "transcription",
}

UNPLACED = "mechanism not established"

# Ordered for reading. "Mechanism not established" is last, not because
# those routes are weakest - several rest on a direct measurement of the
# gene - but because the group says nothing about where to look next.
STEP_ORDER = (
    "acts on the gene or its product",
    "transcription",
    "splicing and isoform choice",
    "mRNA stability",
    "translation",
    "folding",
    "degradation",
    UNPLACED,
)


def step_of(route: dict) -> str:
    """The step this route's evidence places it at, or `UNPLACED`.

    Deliberately returns the honest answer rather than a best guess: a
    route whose evidence is "the gene moved" is not evidence about which
    step moved it.
    """
    src = str(route.get("source") or "")
    if src == "direct":
        return "acts on the gene or its product"
    if src == "proteomics":
        # Measured on the protein, which bounds the answer no further:
        # a stabiliser, a translation effect and a transcriptional one
        # all show up as more protein.
        return UNPLACED
    cat = str(route.get("category") or "").lower()
    if "protein abundance" in cat:
        return UNPLACED
    return STEP_BY_SOURCE.get(src, UNPLACED)


def _key(name: str) -> str:
    from .signatures import breadth_key
    return breadth_key(name or "")


def build(result: dict) -> dict:
    """One row per distinct compound, grouped by step, nothing ranked."""
    from . import strata

    blocks = result.get("results") or []
    routes = [r for b in blocks
              for r in ((b.get("all_routes") or {}).get("routes") or [])]
    if not routes:
        return {"rows": [], "excluded": [], "n_compounds": 0}

    excluded = strata.exclusions(routes)
    out_keys = {_key(e.get("drug") or "") for e in excluded}

    occupancy = {}
    for b in blocks:
        occ = (b.get("scan_quality") or {}).get("occupancy") or {}
        for sym in occ.get("supported") or []:
            occupancy[sym.upper()] = True
        for sym in occ.get("unsupported") or []:
            occupancy.setdefault(sym.upper(), False)

    by_drug: dict = {}
    for r in routes:
        name = (r.get("drug") or "").strip()
        if not name:
            continue
        k = _key(name)
        if k in out_keys:
            continue
        rec = by_drug.setdefault(k, {
            "compound": name, "steps": set(), "targets": set(),
            "sources": set(), "categories": set(), "approved_for": set(),
            "n_routes": 0, "occupancy": set(),
        })
        rec["n_routes"] += 1
        rec["steps"].add(step_of(r))
        if r.get("acts_on"):
            rec["targets"].add(str(r["acts_on"]))
            sup = occupancy.get(str(r["acts_on"]).upper())
            if sup is not None:
                rec["occupancy"].add("measured" if sup else "predicted only")
        if r.get("source"):
            rec["sources"].add(str(r["source"]))
        if r.get("category"):
            rec["categories"].add(str(r["category"]))
        for field in ("approved_for", "indication"):
            if r.get(field):
                rec["approved_for"].add(str(r[field]))

    rows = []
    for rec in by_drug.values():
        steps = [s for s in STEP_ORDER if s in rec["steps"]]
        rows.append({
            "compound": rec["compound"],
            # A compound found by two routes at different steps keeps
            # both: which one is right is what the screen would settle.
            "step": steps[0] if steps else UNPLACED,
            "all_steps": steps,
            "n_routes": rec["n_routes"],
            "targets": sorted(rec["targets"])[:6],
            "sources": sorted(rec["sources"]),
            "site_evidence": ("measured" if "measured" in rec["occupancy"]
                              else ("predicted only" if rec["occupancy"]
                                    else "")),
            "approved_for": sorted(rec["approved_for"])[:2],
        })
    # Alphabetical inside a step, and said to be arbitrary. Sorting by
    # anything else would reintroduce the ranking this file exists to
    # avoid.
    rows.sort(key=lambda r: (STEP_ORDER.index(r["step"]),
                             r["compound"].lower()))
    log.info("screen list: %d compound(s) over %d route(s), %d excluded",
             len(rows), len(routes), len(excluded))
    return {"rows": rows, "excluded": excluded, "n_compounds": len(rows),
            "n_routes": len(routes)}


def to_tsv(plan: dict) -> str:
    """The plate list as a table, one compound per line.

    Through `csv.writer`, as `report.write_outputs` does for
    tf_summary.tsv and for the same reason: compound names and indication
    text come from external tables and APIs, and one carrying a tab or a
    newline split the row and shifted every column after it. This is the
    file someone orders compounds from, so it is the worst place for a
    row to say something other than what it means.
    """
    import csv
    import io

    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n",
                   quoting=csv.QUOTE_MINIMAL)
    w.writerow(["compound", "step", "site_evidence", "n_routes", "targets",
                "sources", "approved_for"])
    for r in plan.get("rows") or []:
        w.writerow([
            r["compound"], r["step"], r["site_evidence"] or "-",
            str(r["n_routes"]), ";".join(r["targets"]) or "-",
            ";".join(r["sources"]) or "-",
            ";".join(r["approved_for"]) or "-"])
    return buf.getvalue()


def summarise(plan: dict, gene: str) -> str:
    """What the list is, and the two things that would be misread."""
    n = plan.get("n_compounds") or 0
    if not n:
        return (f"**No compound reached the screen list for {gene}.** "
                f"Nothing was found that is not excluded on safety grounds.")
    by_step: dict = {}
    for r in plan["rows"]:
        by_step[r["step"]] = by_step.get(r["step"], 0) + 1
    placed = sum(v for k, v in by_step.items() if k != UNPLACED)
    return (
        f"**{n} distinct compounds, in no order.** Chemical families are "
        f"not folded together here and no compound is demoted for being "
        f"well studied, because a plate measures what a compound does to "
        f"cells and neither of those is about that. Rows are alphabetical "
        f"inside each group; the groups say where the evidence places a "
        f"compound, and for {n - placed} of the {n} it places them "
        f"nowhere - those rest on the gene having moved, with no mechanism "
        f"attached. {len(plan.get('excluded') or [])} further compound(s) "
        f"were found and left off on safety grounds, listed below: that "
        f"judgement is clinical, and some of them are reasonable "
        f"mechanistic controls on a plate.")
