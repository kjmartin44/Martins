"""Change the natural ligand, not the receptor.

A nuclear receptor's activity is set by how much of its endogenous ligand is
around. That concentration is controlled by enzymes — one set that synthesises
the ligand, another that degrades it, sometimes a carrier protein that
sequesters it. Those enzymes are frequently far more druggable than the
receptor, and some of the most successful drugs in the pharmacopoeia work this
way rather than by binding the receptor at all:

    anastrozole    inhibits aromatase  -> less estradiol -> less ER activity
    finasteride    inhibits 5a-reductase -> less DHT -> less AR activity
    metyrapone     inhibits CYP11B1 -> less cortisol -> less GR activity
    levothyroxine  supplies T4 -> more T3 -> more TR activity

None of these appear anywhere else in this package, because every other module
looks for ligands *of the receptor*. This one looks one step sideways.

Four strategies, and the direction each produces
------------------------------------------------
    supplement ligand      -> receptor activity UP
    inhibit catabolism     -> receptor activity UP   (ligand lingers)
    inhibit synthesis      -> receptor activity DOWN
    sequester ligand       -> receptor activity DOWN

Why this is worth trying
------------------------
**The enzyme is often druggable when the receptor is not.** A biosynthetic
cytochrome P450 has a well-formed active site; some receptors do not.

**Physiological patterning survives.** You are scaling an endogenous signal
rather than clamping the receptor on, so circadian and tissue-restricted
patterns persist — the same argument as for a PAM, applied upstream.

**Some routes are tissue-selective.** HSD11B2 converts cortisol to inert
cortisone only in tissues that express it, so inhibiting it raises cortisol
*locally* rather than systemically. Where such an enzyme has restricted
expression, you get selectivity the receptor itself cannot offer.

Why it often fails
------------------
**The ligand hits every receptor it binds, everywhere.** This is the dominant
risk and it is not subtle: aromatase inhibitors cause bone loss because
estradiol does much more than drive your gene. Changing a systemic hormone
level has a systemic phenotype, and that phenotype is usually the reason the
drug is hard to use.

**Feedback opposes you.** Endocrine axes are homeostatic. Lower cortisol and
ACTH rises to restore it; the effect you get is smaller and slower than the
enzyme inhibition suggests, and can rebound.

**Ligand supply may not be limiting.** If the receptor is already saturated,
raising the ligand does nothing. Check whether the axis has headroom before
assuming a supplement will move anything.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

log = logging.getLogger(__name__)

STRATEGY_DIRECTION = {
    "supplement ligand": "up",
    "inhibit catabolism": "up",
    "inhibit synthesis": "down",
    "sequester ligand": "down",
}


@dataclass
class LigandRoute:
    receptor: str
    ligand: str
    strategy: str
    acts_on: str                # the enzyme or carrier the drug binds
    drug: str
    direction: str              # effect on receptor activity
    status: str                 # approved | clinical | preclinical
    approval_year: int | None
    tissue_selective: bool = False
    note: str = ""
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _r(receptor, ligand, strategy, enzyme, drug, status, year,
       note="", tissue_selective=False, caveats=()):
    r = LigandRoute(receptor, ligand, strategy, enzyme, drug,
                    STRATEGY_DIRECTION[strategy], status, year,
                    tissue_selective, note, list(caveats))
    if status != "approved":
        r.caveats.append(f"NOT approved ({status})")
    return r


# --- Curated. Verify anything you act on; approvals and labels change. -----
ROUTES: dict[str, list[LigandRoute]] = {}
for _e in [
    # ---- ESR1 / estradiol ------------------------------------------------
    ("ESR1", "estradiol", "inhibit synthesis", "CYP19A1 (aromatase)",
     "Anastrozole", "approved", 1995,
     "Aromatase inhibition lowers systemic estradiol",
     False, ["estradiol maintains bone: aromatase inhibitors cause "
             "osteoporosis and arthralgia"]),
    ("ESR1", "estradiol", "inhibit synthesis", "CYP19A1 (aromatase)",
     "Letrozole", "approved", 1997, "", False, []),
    ("ESR1", "estradiol", "inhibit synthesis", "CYP19A1 (aromatase)",
     "Exemestane", "approved", 1999, "Irreversible (steroidal) inhibitor",
     False, []),
    ("ESR1", "estradiol", "supplement ligand", "—", "Estradiol",
     "approved", 1975, "Direct hormone replacement", False, []),

    # ---- AR / dihydrotestosterone ----------------------------------------
    ("AR", "dihydrotestosterone", "inhibit synthesis", "SRD5A2 (5a-reductase)",
     "Finasteride", "approved", 1992,
     "Blocks testosterone -> DHT conversion; DHT is the high-affinity AR ligand",
     True, ["5a-reductase is expressed mainly in prostate, skin and hair "
            "follicle, so the effect is partly tissue-restricted — but "
            "sexual and mood adverse effects are well documented"]),
    ("AR", "dihydrotestosterone", "inhibit synthesis", "SRD5A1/2",
     "Dutasteride", "approved", 2001, "Dual isoform inhibitor", True, []),
    ("AR", "testosterone", "inhibit synthesis", "CYP17A1",
     "Abiraterone", "approved", 2011,
     "Blocks androgen synthesis upstream, including adrenal and intratumoural",
     False, ["requires co-administered glucocorticoid because CYP17A1 "
             "blockade also raises mineralocorticoid precursors"]),
    ("AR", "testosterone", "supplement ligand", "—", "Testosterone",
     "approved", 1953, "", False, []),

    # ---- NR3C1 / cortisol ------------------------------------------------
    ("NR3C1", "cortisol", "inhibit synthesis", "CYP11B1 (11b-hydroxylase)",
     "Metyrapone", "approved", 1961, "Lowers cortisol synthesis", False,
     ["ACTH rises in response, partly overriding the block"]),
    ("NR3C1", "cortisol", "inhibit synthesis", "CYP11B1",
     "Osilodrostat", "approved", 2020, "Potent CYP11B1 inhibitor", False, []),
    ("NR3C1", "cortisol", "inhibit synthesis", "CYP17A1/CYP11B1",
     "Levoketoconazole", "approved", 2021, "", False, []),
    ("NR3C1", "cortisol", "inhibit catabolism", "HSD11B2",
     "Carbenoxolone", "approved", 1962,
     "Blocks cortisol -> cortisone inactivation, raising cortisol inside "
     "tissues that express HSD11B2",
     True, ["the classic tissue-selective example, but HSD11B2 inhibition in "
            "kidney causes sodium retention and hypertension"]),
    ("NR3C1", "cortisol", "supplement ligand", "—", "Hydrocortisone",
     "approved", 1952, "", False, []),

    # ---- THRA / THRB / thyroid hormone ------------------------------------
    ("THRB", "triiodothyronine", "supplement ligand", "—", "Levothyroxine",
     "approved", 1955, "T4, converted peripherally to active T3", False, []),
    ("THRB", "triiodothyronine", "supplement ligand", "—", "Liothyronine",
     "approved", 1956, "T3 directly", False, []),
    ("THRB", "triiodothyronine", "inhibit synthesis", "TPO",
     "Methimazole", "approved", 1950, "", False,
     ["TSH feedback rises; agranulocytosis is a rare but serious risk"]),
    ("THRB", "triiodothyronine", "inhibit synthesis", "TPO",
     "Propylthiouracil", "approved", 1947,
     "Also blocks peripheral T4 -> T3 conversion", False, []),

    # ---- VDR / calcitriol -------------------------------------------------
    ("VDR", "calcitriol", "supplement ligand", "—", "Cholecalciferol",
     "approved", 1940, "Precursor; requires intact 25- and 1a-hydroxylation",
     False, []),
    ("VDR", "calcitriol", "supplement ligand", "—", "Calcifediol",
     "approved", 2016, "25-OH-D, bypasses hepatic hydroxylation", False, []),
    ("VDR", "calcitriol", "supplement ligand", "—", "Alfacalcidol",
     "approved", 1978, "1a-OH-D, bypasses renal CYP27B1", False, []),
    ("VDR", "calcitriol", "inhibit catabolism", "CYP24A1",
     "Ketoconazole", "approved", 1981,
     "Non-selective CYP inhibition raises calcitriol; used in CYP24A1 "
     "deficiency and sarcoidosis-related hypercalcaemia",
     False, ["profoundly non-selective across CYPs; hepatotoxicity"]),

    # ---- RARA / retinoic acid ---------------------------------------------
    ("RARA", "all-trans retinoic acid", "supplement ligand", "—",
     "Tretinoin", "approved", 1971, "", False, []),
    ("RARA", "all-trans retinoic acid", "inhibit catabolism", "CYP26A1/B1",
     "Talarozole", "clinical", None,
     "RAMBA (retinoic acid metabolism blocking agent): raises endogenous ATRA "
     "rather than supplying exogenous retinoid",
     True, ["CYP26 expression is induced by retinoic acid itself, so the "
            "system compensates"]),
    ("RARA", "all-trans retinoic acid", "inhibit catabolism", "CYP26",
     "Liarozole", "clinical", None, "Orphan designation, not FDA approved",
     False, []),

    # ---- NR1H4 (FXR) / bile acids -----------------------------------------
    ("NR1H4", "bile acids", "sequester ligand", "intestinal bile acid pool",
     "Cholestyramine", "approved", 1973,
     "Binds bile acids in the gut, lowering FXR activation", False, []),
    ("NR1H4", "bile acids", "sequester ligand", "intestinal bile acid pool",
     "Colesevelam", "approved", 2000, "", False, []),

    # ---- AHR / kynurenine --------------------------------------------------
    ("AHR", "kynurenine", "inhibit synthesis", "IDO1",
     "Epacadostat", "clinical", None,
     "Lowers kynurenine, an endogenous AHR ligand", False,
     ["failed its phase 3 combination endpoint"]),
]:
    (rec, lig, strat, enz, drug, status, year, note, tsel, cav) = _e
    ROUTES.setdefault(rec, []).append(
        _r(rec, lig, strat, enz, drug, status, year, note, tsel, cav))


def _enzyme_key(name: str) -> str:
    """Normalise an enzyme label so 'CYP11B1' and 'CYP11B1 (11b-hydroxylase)'
    are not listed as two different enzymes."""
    return name.split("(")[0].strip()


def chain(route: "LigandRoute", gene: str = "your gene") -> str:
    """The full causal chain, written out.

    A route is only actionable if you can read the chain and judge each link.
    Printing just the drug name hides three inferential steps.
    """
    # Supplementation has one fewer link: no enzyme sits in the path.
    if route.strategy == "supplement ligand":
        return (f"{route.drug} supplies {route.ligand} → more "
                f"{route.receptor} activity → {route.direction} at {gene}")

    more_less = "more" if route.direction == "up" else "less"
    return (f"{route.drug} blocks {_enzyme_key(route.acts_on)} → "
            f"{more_less} {route.ligand} → {more_less} {route.receptor} "
            f"activity → {route.direction} at {gene}")


def covered_receptors() -> dict:
    """Which receptors have a curated ligand axis, and what the ligand is."""
    return {rec: routes[0].ligand for rec, routes in sorted(ROUTES.items())}


def routes_for(receptor: str, goal: str = "up", *,
               approved_only: bool = False) -> list[LigandRoute]:
    """Ligand-level routes that move this receptor in the wanted direction."""
    out = [r for r in ROUTES.get(receptor.upper(), []) if r.direction == goal]
    if approved_only:
        out = [r for r in out if r.status == "approved"]
    # Tissue-selective first: that is the whole advantage when it exists.
    out.sort(key=lambda r: (not r.tissue_selective, r.status != "approved",
                            -(r.approval_year or 0)))
    return out


def has_endogenous_ligand(receptor: str) -> bool:
    return receptor.upper() in ROUTES


def assess(receptor: str, goal: str = "up") -> dict:
    """Ligand-modulation options for one receptor, with the standing caveats."""
    hits = routes_for(receptor, goal)
    opposite = routes_for(receptor, "down" if goal == "up" else "up")
    result = {
        "receptor": receptor.upper(),
        "has_ligand_axis": has_endogenous_ligand(receptor),
        "ligand": (hits or opposite)[0].ligand if (hits or opposite) else None,
        "routes": [r.to_dict() for r in hits],
        "chains": [chain(r) for r in hits],
        "opposite_direction": [
            {"drug": r.drug, "strategy": r.strategy, "acts_on": r.acts_on,
             "status": r.status} for r in opposite],
        "enzymes": sorted({_enzyme_key(r.acts_on) for r in hits + opposite
                           if r.acts_on != "—"}),
        "n_approved": sum(1 for r in hits if r.status == "approved"),
        "notes": [],
    }
    if not result["has_ligand_axis"]:
        result["notes"].append(
            f"{receptor.upper()} has no curated endogenous-ligand axis here. "
            "Either it is ligand-independent (many TFs are), or its ligand and "
            "the enzymes controlling it are not in this table — check the "
            "literature before concluding there is no route.")
        return result
    if not hits:
        other = "down" if goal == "up" else "up"
        if routes_for(receptor, other):
            result["notes"].append(
                f"Ligand routes exist for {receptor.upper()} but only in the "
                f"{other} direction. Moving an endogenous ligand the other way "
                "usually means supplementation, which may be limited by "
                "receptor saturation.")
        return result

    if any(r.tissue_selective for r in hits):
        result["notes"].append(
            "At least one route is tissue-selective: the controlling enzyme is "
            "not expressed everywhere, so the ligand change is local rather "
            "than systemic. That is selectivity the receptor itself cannot "
            "offer, and the main reason to prefer this approach.")
    result["notes"].append(
        "Standing caveat: the ligand acts on every receptor it binds, in every "
        "tissue. Changing a systemic hormone level has a systemic phenotype — "
        "aromatase inhibitors cause bone loss because estradiol does far more "
        "than drive your gene.")
    result["notes"].append(
        "Endocrine axes are homeostatic. Feedback will partly oppose you, so "
        "the effect is smaller and slower than the enzyme inhibition suggests, "
        "and it can rebound on withdrawal.")
    if any(r.strategy == "supplement ligand" for r in hits):
        result["notes"].append(
            "Supplementation only works if ligand supply is limiting. If the "
            "receptor is already saturated in your tissue, adding more does "
            "nothing — check the axis has headroom first.")
    return result
