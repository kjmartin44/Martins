"""Where CTD and LINCS speak about the same drug, and where either contradicts a route.

Why this exists
---------------
CTD records, from the literature, which way a chemical moved the gene. LINCS
records, from one set of cell lines, which way a compound moved it. A route
predicts a third thing: that a drug acting on a factor should move the gene the
wanted way. The three are kept apart in the results and were never compared, so
a drug could rank second while CTD recorded it moving the gene the other way.

What this does
--------------
Reads the CTD record already in the result (`improvements.reported_links`) and
the local LINCS consensus signatures, and reports

  * the drugs CTD records moving the gene the opposite way to the goal, so a
    ranked drug can carry that beside it, and
  * the drugs on which CTD and LINCS both give a direction and the directions
    differ.

What it does not do
-------------------
**It changes no tier and no score.** A caution is a reason to check, not a
verdict. CTD aggregates cell types, doses and assays, so a chemical can raise a
gene in one and lower it in another. Against LINCS, CTD's direction held for 24
of 37 pairs where both gave one (65%), and most CTD drugs have no LINCS call at
all, so a disagreement is common and a weak sign of which source is wrong.

**No LINCS file means no comparison.** The result then says LINCS was not
checked, which is not the same as the two agreeing.
"""

from __future__ import annotations

import logging

from .signatures import breadth_key

log = logging.getLogger(__name__)

# What a route's own notes say when no direction is recorded for the drug at
# the factor. The catalogue writes these; they are read here, not redefined.
UNKNOWN_DIRECTION_MARKERS = ("is not recorded", "unknown action")


def _names(entries) -> dict:
    """`{key: display name}` for the chemicals (not exposures or headings)."""
    out = {}
    for e in entries or []:
        if e.get("kind", "chemical") != "chemical":
            continue
        key = breadth_key(e.get("chemical") or "")
        if key:
            out.setdefault(key, e["chemical"])
    return out


def _label(wanted: bool, goal: str) -> str:
    """The direction a source records, as "up" or "down"."""
    if wanted:
        return goal
    return "down" if goal == "up" else "up"


def build(gene: str, goal: str, improvements: dict | None, *,
          use_lincs: bool = True) -> dict | None:
    """The cross-check for one gene, or None when CTD was not read.

    `improvements` is the block's `improvements` dict. LINCS comes from the
    local signature files (see `retrospective.lincs_keyset`); without them the
    CTD half is still returned and `lincs.checked` says so.
    """
    rl = (improvements or {}).get("reported_links")
    if not rl or goal not in ("up", "down"):
        return None
    wanted = _names(rl.get("wanted_direction"))
    opposite = _names(rl.get("opposite_direction"))
    out = {"gene": gene, "goal": goal,
           "ctd_opposite": sorted(opposite.values()),
           "ctd_wanted": sorted(wanted.values()),
           "lincs": {"checked": False, "reason": None,
                     "n_wanted": 0, "n_opposite": 0},
           "n_compared": 0, "n_agree": 0, "conflicts": []}
    if not use_lincs:
        out["lincs"]["reason"] = "LINCS was not requested for this run"
        return out
    try:
        from . import retrospective as rt
        ks, why = rt.lincs_keyset(gene, goal)
    except Exception as exc:                                # noqa: BLE001
        ks, why = None, f"LINCS could not be read ({str(exc)[:60]})"
    if ks is None:
        out["lincs"]["reason"] = why
        return out
    out["lincs"].update(checked=True, n_wanted=len(ks.wanted),
                        n_opposite=len(ks.opposite))
    for key, name in sorted({**opposite, **wanted}.items(),
                            key=lambda kv: kv[1].lower()):
        if key in ks.wanted:
            lincs_says = True
        elif key in ks.opposite:
            lincs_says = False
        else:
            continue
        ctd_says = key in wanted
        out["n_compared"] += 1
        if ctd_says == lincs_says:
            out["n_agree"] += 1
        else:
            out["conflicts"].append({"drug": name,
                                     "ctd": _label(ctd_says, goal),
                                     "lincs": _label(lincs_says, goal)})
    return out


def unknown_direction(route: dict) -> bool:
    """Does the route say that the drug's direction at its factor is not recorded?"""
    notes = " ".join(route.get("notes") or []).lower()
    return any(m in notes for m in UNKNOWN_DIRECTION_MARKERS)


def opposite_keys(xc: dict | None) -> set:
    """Salt-free keys of the drugs CTD records moving the gene the other way."""
    return {breadth_key(n) for n in (xc or {}).get("ctd_opposite") or []} - {""}


def cautions_for(route: dict, xc: dict | None) -> list:
    """Short cautions for one route, from the cross-check.

    Literature routes are CTD itself, so they carry none.
    """
    if not xc or (route.get("category") or "").startswith("literature:"):
        return []
    key = breadth_key(route.get("drug") or "")
    if not key or key not in opposite_keys(xc):
        return []
    if unknown_direction(route):
        return ["CTD records the opposite direction and this route's own "
                "direction is unknown"]
    return ["CTD records the opposite direction"]


def conflict_names(xc: dict | None) -> set:
    return {breadth_key(c["drug"]) for c in (xc or {}).get("conflicts") or []}


def markdown(xc: dict | None, ranked: list | None = None,
             limit: int = 20) -> list:
    """Report lines for the cross-check; empty when there is nothing to say.

    `ranked` is the routes to check for the opposite-direction caution, best
    first.
    """
    if not xc:
        return []
    gene, goal = xc["gene"], xc["goal"]
    L = ["## Where the sources disagree", ""]
    lincs = xc["lincs"]
    flagged = []
    seen = set()
    for r in ranked or []:
        key = breadth_key(r.get("drug") or "")
        if key in seen or not cautions_for(r, xc):
            continue
        seen.add(key)
        flagged.append(r)
    if flagged:
        L.append(f"**{len(flagged)} ranked drugs are recorded by CTD as moving "
                 f"{gene} the other way** (the goal is to move it {goal}). A "
                 f"route predicts the drug should help; CTD says someone "
                 f"measured the opposite. Nothing is demoted for this, since "
                 f"CTD mixes cell types and doses, but it is the first thing "
                 f"to check before testing the drug.")
        L.append("")
        L.append("| Drug | Tier | Acts on | Caution |")
        L.append("|---|---|---|---|")
        for r in flagged[:limit]:
            L.append(f"| {r.get('drug')} | {r.get('tier')} | "
                     f"{r.get('acts_on') or '-'} | "
                     f"{cautions_for(r, xc)[0]} |")
        if len(flagged) > limit:
            L.append(f"| ... {len(flagged) - limit} more | | | |")
        L.append("")
    if lincs["checked"]:
        n, a = xc["n_compared"], xc["n_agree"]
        if n:
            L.append(f"**CTD and LINCS** both give a direction for {n} "
                     f"drug{'s' if n != 1 else ''} and agree on {a} "
                     f"({a / n:.0%}).")
            L.append("")
        if xc["conflicts"]:
            L.append("| Drug | CTD records | LINCS records |")
            L.append("|---|---|---|")
            for c in xc["conflicts"][:limit]:
                L.append(f"| {c['drug']} | {c['ctd']} | {c['lincs']} |")
            if len(xc["conflicts"]) > limit:
                L.append(f"| ... {len(xc['conflicts']) - limit} more | | |")
            L.append("")
            L.append("One of the two is wrong for the cells you care about, "
                     "or the drug acts differently across cell types. LINCS "
                     "is a few cell lines at one dose and time, so neither "
                     "is a ruling.")
            L.append("")
    else:
        L.append(f"LINCS was not compared with CTD: {lincs['reason']}.")
        L.append("")
    return L if len(L) > 2 else []
