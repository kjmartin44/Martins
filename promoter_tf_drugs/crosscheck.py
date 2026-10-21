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
**It changes a tier in one case only.** A route whose drug has no recorded
direction at its factor, for a drug CTD records moving the gene the other way, is
lowered to tier 4 (`demotes`). Its tier rested on measured binding alone, so the
route was a coin flip on direction, and CTD's record tips it toward the wrong
way. Every other caution changes nothing. A caution is a reason to check, not a
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

# The tier a route falls to when its direction is unknown and CTD records the
# opposite. Tier 4 is "predicted": the other routes there have no measured step
# at the factor either.
DEMOTED_TIER = 4


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
           # What KIND of thing moves the gene your way. The list was
           # already here and only ever used for a novelty check.
           "wanted_profile": mover_profile(sorted(wanted.values()), goal),
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


def demotes(route: dict, xc: dict | None) -> bool:
    """Should this route fall to `DEMOTED_TIER`?

    True for a route that is not itself a CTD route, whose drug has no recorded
    direction at its factor, and which CTD records moving the gene the other way.
    A route with a recorded direction that CTD contradicts is a disagreement
    between two measurements and is flagged, not demoted.
    """
    return bool(cautions_for(route, xc)) and unknown_direction(route)


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
                 f"measured the opposite. A route whose own direction at its "
                 f"factor is not recorded is lowered to tier {DEMOTED_TIER}; "
                 f"any other is only flagged, since CTD mixes cell types and "
                 f"doses. Check it before testing the drug.")
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
# ===================================================================== #
# What KIND of chemical moves this gene the way you want?               #
# ===================================================================== #
# A keyword vocabulary, not a toxicology database. It is here to answer
# one coarse question - are the things that move this gene the way you
# want them to mostly medicines, or mostly poisons - and it answers that
# well enough to change a decision. It cannot grade a hazard, it does not
# know a dose, and a chemical it does not match is UNCLASSIFIED, which is
# not a clean bill.
#
# Matching is on lowercased substrings of the CTD chemical name.
MOVER_CLASSES = (
    ("carcinogen or genotoxin", (
        "dimethylhydrazine", "methylcholanthrene", "benzo(a)pyrene",
        "benzo[a]pyrene", "aflatoxin", "aristolochic", "acrylamide",
        "glycidol", "n-nitroso", "nitrosamine", "ochratoxin",
        "safrole", "urethane", "vinyl chloride", "diethylnitrosamine",
        "4-aminobiphenyl", "2-acetylaminofluorene", "thioacetamide",
        "carbon tetrachloride", "imidazo(4,5-b)pyridine", "dibenzo",
        "chromium", "cadmium", "arsen", "nickel",
    )),
    ("tumour promoter", (
        "tetradecanoylphorbol", "phorbol", "croton oil", "phenobarbital",
        "okadaic acid", "cuprizone", "clofibrate", "pirinixic",
        "wy-14643", "peroxisome prolifer",
    )),
    ("hormone or mitogen", (
        "estradiol", "estrogen", "dihydrotestosterone", "testosterone",
        "progesterone", "genistein", "coumestrol", "diethylstilbestrol",
        "bisphenol", "growth factor", "insulin",
    )),
    ("inflammogen or irritant", (
        "lipopolysaccharid", "oxazolone", "isocyanate", "trimellitic",
        "latex", "carrageenan", "turpentine",
    )),
    ("particulate or metal salt", (
        "nanotube", "titanium dioxide", "silicon dioxide", "silica",
        "silver nitrate", "copper sulfate", "asbestos", "carbon black",
        "particulate",
    )),
    ("pesticide or industrial chemical", (
        "chlorpyrifos", "imidacloprid", "acetamiprid", "decamethrin",
        "deltamethrin", "epoxiconazole", "picoxystrobin", "pcb ",
        "polychlorinated", "dioxin", "tetrafluoro", "perfluoro",
        "toluene", "dichloroacetic", "rotenone", "paraquat",
        "tributyltin", "bis(tri-n-butyltin)",
    )),
    ("nutrient or dietary", (
        "glucose", "fructose", "sucrose", "choline", "methionine",
        "dietary", "ethanol", "vitamin",
    )),
)


def classify_mover(name: str) -> str | None:
    """The coarse class of one CTD chemical, or None when unmatched."""
    low = (name or "").lower()
    for label, keys in MOVER_CLASSES:
        if any(k in low for k in keys):
            return label
    return None


def mover_profile(names: list, goal: str = "up") -> dict:
    """What kind of thing the literature says moves this gene your way.

    Returns the counts by class, how many are approved drugs, and a
    verdict sentence. `approved` is answered by the same Orange Book the
    report filter uses, so a chemical that is also a medicine - phenobarbital
    is both an approved barbiturate and a classic rodent tumour promoter -
    is counted in both columns rather than forced into one.
    """
    from .catalogue import verify_approval

    names = [n for n in (names or []) if n]
    if not names:
        return {}
    by_class: dict = {}
    for n in names:
        by_class.setdefault(classify_mover(n) or "unclassified",
                            []).append(n)
    approved = [n for n in names if verify_approval(n) is True]
    hazard = sum(len(v) for k, v in by_class.items()
                 if k not in ("unclassified", "nutrient or dietary"))
    share = hazard / len(names)
    word = "raise" if goal == "up" else "lower"

    if approved:
        drugs = (f"{len(approved)} of them are FDA approved drugs "
                 f"({', '.join(sorted(approved)[:6])}"
                 + (" …" if len(approved) > 6 else "") + ")")
    else:
        drugs = "none of them is an FDA approved drug"

    if share >= 0.5:
        verdict = (
            f"**Most of what is recorded as {word[:-1]}ing this gene is not "
            f"a medicine.** Of {len(names)} chemicals CTD records as "
            f"{word[:-1]}ing it, {hazard} fall into a carcinogen, tumour "
            f"promoter, mitogen, inflammogen or industrial class, and "
            f"{drugs}. That pattern is worth pausing on: when the published "
            f"way to {word} a gene is to damage or inflame the tissue, the "
            f"gene is most likely responding to that insult rather than "
            f"being regulated in a way a drug could copy. Check whether the "
            f"transcriptional objective is the right objective before "
            f"reading the ranking below.")
    else:
        verdict = (
            f"Of {len(names)} chemicals CTD records as {word[:-1]}ing this "
            f"gene, {hazard} fall into a hazard class and {drugs}. The "
            f"classes are a keyword vocabulary, so an unclassified chemical "
            f"is unexamined rather than safe.")

    return {
        "n": len(names),
        "goal": goal,
        "by_class": {k: sorted(v) for k, v in sorted(by_class.items())},
        "counts": {k: len(v) for k, v in sorted(by_class.items())},
        "n_hazard": hazard,
        "hazard_share": round(share, 3),
        "approved": sorted(approved),
        "n_approved": len(approved),
        "verdict": verdict,
    }
