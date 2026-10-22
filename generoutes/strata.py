"""What to test first, grouped by what was actually measured.

Why this replaces a ranked list
-------------------------------
The summary used to open with a numbered top ten, ordered by tier and
then score. That ordering has been tested against CTD and it carries far
less than it looks: the routes order drugs no better than chance once
drug fame is held equal (AUC 0.53, against 0.93 for a baseline that only
counts how many other genes record a drug), and within a tier most of the
top rows share an identical score anyway. A numbered list implies a
precision the package has measured itself not to have, and "#1" is the
one thing a reader takes away.

So the rows are grouped instead. The groups are ordered - what was
measured on your gene really does beat what was inferred about a factor -
and inside a group nothing is. That is the honest shape of the evidence.

How the groups are drawn
------------------------
Not by tier. Tier mixes two questions: how good the observation is, and
what it was an observation OF. A compound whose binding to a receptor is
measured and a compound measured to move your gene are both "tier 2", and
they are not the same kind of claim at all. These groups answer the
second question, which is the one that decides what an experiment would
test.
"""

from __future__ import annotations

# `source` strings as catalogue.py emits them. Read from there rather
# than guessed; an invented name would silently drop routes into the
# weakest group.
ON_THE_GENE = frozenset({"direct"})
MEASURED_ON_GENE = frozenset({"proteomics", "signatures", "ctd"})
MEASURED_ON_TARGET = frozenset({
    "drug_scoring", "screening", "docking", "regulators", "sensitivity",
    "abundance", "tfactivity",
})

STRATA = (
    ("on_the_gene", "Already acts on this gene",
     "A drug annotated against the gene or its product. Nothing is "
     "inferred between the compound and the target, which is why these "
     "sit above everything else whatever their score."),
    ("measured_on_gene", "Measured on this gene",
     "Somebody dosed cells and this gene, or its protein, moved. The "
     "mechanism is usually unknown, so these are strong about THAT and "
     "weak about WHY - the opposite of the group below."),
    ("measured_on_target", "Measured on the drug's target",
     "The compound's effect on its target is observed; the step from "
     "that target to your gene is inferred. Strong about mechanism, "
     "unverified about outcome."),
    ("predicted", "Predicted throughout",
     "Neither link is measured here. Worth testing only when something "
     "above has run out, or when a group above is empty for a reason "
     "the evidence table explains."),
)

_TITLES = dict((k, t) for k, t, _ in STRATA)


def stratum_of(route: dict) -> str:
    src = (route or {}).get("source") or ""
    if src in ON_THE_GENE:
        return "on_the_gene"
    if src in MEASURED_ON_GENE:
        return "measured_on_gene"
    if src in MEASURED_ON_TARGET:
        return "measured_on_target"
    return "predicted"


def n_methods(route: dict, by_drug: dict | None = None) -> int:
    """How many INDEPENDENT methods found this drug.

    The most useful number in the report, and the only thing used to
    order rows inside a group - because unlike the score it has not been
    measured to be chance. It is still not a vote: two weak agreements do
    not equal one measured result.
    """
    from .report import METHOD_OF_SOURCE

    if by_drug is None:
        return 1 if METHOD_OF_SOURCE.get(route.get("source")) else 0
    key = (route.get("drug") or "").strip().lower()
    return len({METHOD_OF_SOURCE.get(r.get("source"))
                for r in by_drug.get(key, [])
                if METHOD_OF_SOURCE.get(r.get("source"))})


def methods_by_drug(all_routes: list) -> dict:
    """`{drug lowercased: [route, ...]}` across every route in the run."""
    out: dict = {}
    for r in all_routes or []:
        key = (r.get("drug") or "").strip().lower()
        if key:
            out.setdefault(key, []).append(r)
    return out


def partition(rows: list, all_routes: list | None = None) -> list:
    """`[(key, title, blurb, rows)]` for the non-empty groups, in order.

    Inside a group the order is by how many independent methods found the
    drug, then by name. Never by score: the score's ordering is the thing
    this module exists to stop implying.
    """
    by_drug = methods_by_drug(all_routes or [])
    buckets: dict = {}
    for r in rows or []:
        buckets.setdefault(stratum_of(r), []).append(r)
    out = []
    for key, title, blurb in STRATA:
        got = buckets.get(key) or []
        if not got:
            continue
        got.sort(key=lambda r: (-n_methods(r, by_drug),
                                (r.get("drug") or "").lower()))
        out.append((key, title, blurb, got))
    return out


def exclusions(all_routes: list) -> list:
    """Compounds a gate pushed to the bottom, with the reason.

    Shown rather than dropped. A reader who expected azacitidine to be
    here should learn that it was found and rejected, not be left to
    wonder whether the search missed it.
    """
    seen: dict = {}
    for r in all_routes or []:
        why = None
        if r.get("cancer_risk"):
            why = r["cancer_risk"]
        elif r.get("genotoxic"):
            why = ("works by damaging DNA or the spindle, which is the "
                   "defect being treated")
        if not why:
            continue
        name = (r.get("drug") or "?").strip()
        seen.setdefault(name.lower(), {"drug": name, "why": why,
                                       "acts_on": r.get("acts_on")})
    return sorted(seen.values(), key=lambda d: d["drug"].lower())


def preamble(n_groups: int) -> str:
    """The sentence that has to be there if the grouping is to mean anything."""
    return (
        f"Grouped by what was actually measured, in {n_groups} group(s), "
        f"strongest kind of evidence first. **Within a group nothing is "
        f"ranked** - the rows are ordered by how many independent methods "
        f"found the drug, and beyond that alphabetically. That is "
        f"deliberate: tested against CTD, this package's scoring orders "
        f"drugs no better than chance once drug fame is held equal, so a "
        f"numbered list would imply a precision it does not have. Use the "
        f"groups, the method count and the exclusions; do not read "
        f"position within a group as priority.")
