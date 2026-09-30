"""Orthosteric agonists vs allosteric modulators: a comparison, not a ranking.

A positive allosteric modulator (PAM) is not simply "another way to activate
the receptor". It differs on four axes that matter for what this pipeline is
trying to do, and the differences do not all point the same way.

Where PAMs win
--------------
**Subtype selectivity.** The orthosteric pocket is under selection to bind the
same endogenous ligand across a receptor family, so it is highly conserved and
hard to discriminate. Allosteric sites carry no such constraint and diverge
freely. This is why subtype-selective PAMs exist for targets where decades of
orthosteric medicinal chemistry produced nothing selective. For a pipeline
whose main failure mode is "the drug moves far more than your gene", this is
the single strongest argument for looking at PAMs.

**Preserved signalling dynamics.** An agonist clamps the receptor on wherever
the drug reaches. A PAM amplifies whatever the endogenous ligand is already
doing, so spatial and temporal patterning survives - tissue-restricted,
pulsatile or circadian signalling is scaled rather than overwritten. For a TF
under circadian or feedback control this can be the difference between a
physiological increase and a flat, dysregulated one.

**Ceiling effect.** PAM activity saturates once the allosteric site is
occupied; effect is bounded by available endogenous ligand. That caps overdose
risk and generally reduces receptor desensitisation and downregulation, which
chronic full agonists provoke.

Where PAMs lose
---------------
**They need endogenous tone.** No ligand, no effect. If the pathway is silent
in your cell type, or the patient's ligand is absent, a PAM does nothing at
all. An orthosteric agonist does not care. This is the failure mode to check
first, and nothing in this package can check it for you - it needs a
measurement in your system.

**Weaker maximal effect.** Bounded by definition. If you need a large fold
change rather than a modest, well-behaved one, an agonist is the honest choice.

**Probe dependence.** A PAM's effect can vary with which endogenous ligand is
present, so results in one assay may not transfer to another.

Nuclear-receptor reality check
------------------------------
Approved PAMs cluster on GPCRs and ion channels (benzodiazepines at GABA-A,
cinacalcet at CaSR, ivacaftor at CFTR). True allosteric modulators of nuclear
receptors are rare in the clinic, though allosteric sites are established -
the PPARgamma alternate site, RXR tetramer interface, Nurr1. So for a TF
target, expect most PAM hits to sit on *upstream signalling nodes* rather than
on the TF itself. That is still useful, and it is labelled as such rather than
being presented as direct TF pharmacology.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

from .drugs import Drug, chembl_drugs, dgidb_drugs
from .http import Client

log = logging.getLogger(__name__)

# ChEMBL action_type -> mechanism class
MECHANISM_CLASS = {
    "AGONIST": "orthosteric agonist",
    "PARTIAL AGONIST": "partial agonist",
    "ACTIVATOR": "activator",
    "OPENER": "activator",
    "POSITIVE ALLOSTERIC MODULATOR": "PAM",
    "POSITIVE MODULATOR": "PAM",
    "ANTAGONIST": "orthosteric antagonist",
    "INVERSE AGONIST": "inverse agonist",
    "NEGATIVE ALLOSTERIC MODULATOR": "NAM",
    "NEGATIVE MODULATOR": "NAM",
    "INHIBITOR": "inhibitor",
    "BLOCKER": "inhibitor",
    "DEGRADER": "degrader",
    "STABILISER": "stabiliser",
}

UP_CLASSES = {"orthosteric agonist", "partial agonist", "activator", "PAM",
              "stabiliser"}
DOWN_CLASSES = {"orthosteric antagonist", "inverse agonist", "NAM", "inhibitor",
                "degrader"}
ALLOSTERIC = {"PAM", "NAM"}


@dataclass
class MechanismProfile:
    """Qualitative properties of a mechanism class. Not fitted - these are
    textbook pharmacology, stated so the comparison is auditable."""
    mechanism: str
    max_effect: str              # bounded | full | partial
    subtype_selectivity: str     # typically higher/lower than orthosteric
    preserves_dynamics: bool
    needs_endogenous_ligand: bool
    desensitisation_risk: str
    selectivity_bonus: float     # applied to the specificity axis
    caveat: str
    surmountable: bool = False   # can rising endogenous ligand overcome it?

    def to_dict(self) -> dict:
        return asdict(self)


PROFILES = {
    "PAM": MechanismProfile(
        "PAM", "bounded (ceiling)", "typically higher", True, True, "lower", 1.30,
        "does nothing without endogenous ligand present - verify pathway tone "
        "in your cell type before anything else", surmountable=False),
    "NAM": MechanismProfile(
        "NAM", "bounded (ceiling)", "typically higher", True, True, "lower", 1.30,
        "only dampens existing signalling; ineffective where tone is already "
        "low. Not surmountable: a surge in endogenous ligand cannot displace "
        "it, because it does not occupy the ligand site",
        surmountable=False),
    "orthosteric agonist": MechanismProfile(
        "orthosteric agonist", "full", "typically lower", False, False, "higher", 1.0,
        "overrides endogenous patterning; conserved pocket makes subtype "
        "selectivity hard and chronic use invites downregulation"),
    "partial agonist": MechanismProfile(
        "partial agonist", "partial", "typically lower", False, False, "moderate", 1.05,
        "submaximal by design; can act as a functional antagonist where "
        "endogenous tone is high"),
    "activator": MechanismProfile(
        "activator", "full", "varies", False, False, "moderate", 1.0,
        "raises activity of existing protein, which is not the same as raising "
        "its abundance"),
    "orthosteric antagonist": MechanismProfile(
        "orthosteric antagonist", "full block", "typically lower", False, True,
        "n/a", 1.0,
        "competitive: occupies the ligand site, so a rise in endogenous ligand "
        "displaces it and the block fades. This matters more than it sounds - "
        "blocking a receptor often triggers feedback that RAISES its ligand",
        surmountable=True),
    "inverse agonist": MechanismProfile(
        "inverse agonist", "below baseline", "typically lower", False, False,
        "n/a", 1.0,
        "suppresses constitutive activity too, so the effect can exceed simple "
        "blockade", surmountable=True),
    "degrader": MechanismProfile(
        "degrader", "full", "varies", False, False, "n/a", 1.10,
        "removes the protein entirely, including any scaffolding role it has "
        "beyond the activity you targeted. Not surmountable by any amount of "
        "ligand - there is no receptor left to bind",
        surmountable=False),
}

DEFAULT_PROFILE = MechanismProfile(
    "unclassified", "unknown", "unknown", False, False, "unknown", 1.0,
    "pharmacology not classified in the source database")


@dataclass
class ModulatorHit:
    drug: str
    target: str
    mechanism: str
    action_type: str
    direction: str
    approved: bool
    withdrawn: bool
    first_approval: int | None
    is_allosteric: bool
    profile: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def classify(drug: Drug) -> str:
    return MECHANISM_CLASS.get((drug.action_type or "").upper(), "unclassified")


def profile_for(mechanism: str) -> MechanismProfile:
    return PROFILES.get(mechanism, DEFAULT_PROFILE)


def allosteric_selectivity_bonus(drug: Drug) -> float:
    """Multiplier applied to the specificity axis for allosteric mechanisms.

    Justified by the conservation argument, not by fitting: the orthosteric
    pocket is constrained to bind the endogenous ligand across a family, so
    allosteric sites are systematically easier to hit selectively.
    """
    return profile_for(classify(drug)).selectivity_bonus


def find_modulators(client: Client, symbol: str, uniprot: str | None = None, *,
                    goal: str = "up", approved_only: bool = True,
                    use_dgidb: bool = True) -> list[ModulatorHit]:
    """All approved drugs at this target, grouped by mechanism class."""
    found = chembl_drugs(client, symbol, uniprot, approved_only)
    if use_dgidb:
        found += dgidb_drugs(client, [symbol], approved_only=approved_only)

    wanted = UP_CLASSES if goal == "up" else DOWN_CLASSES
    hits, seen = [], set()

    for d in found:
        mech = classify(d)
        if mech not in wanted:
            continue
        key = (d.name.lower(), mech)
        if key in seen:
            continue
        seen.add(key)

        prof = profile_for(mech)
        notes = [prof.caveat] if prof.caveat else []
        hits.append(ModulatorHit(
            drug=d.name, target=symbol, mechanism=mech,
            action_type=d.action_type or "UNKNOWN", direction=d.direction,
            approved=d.approved, withdrawn=d.withdrawn,
            first_approval=d.first_approval,
            is_allosteric=mech in ALLOSTERIC,
            profile=prof.to_dict(), notes=notes,
        ))

    # Allosteric first: they are the ones worth comparing against the default.
    hits.sort(key=lambda h: (not h.is_allosteric, h.withdrawn,
                             -(h.first_approval or 0)))
    return hits


def surmountability_note(rows: list[dict], ligand_axis: bool = False
                         ) -> str | None:
    """Flag the competitive-block trap, which bites hardest on downregulation.

    A competitive antagonist occupies the ligand site, so it is displaced when
    the ligand rises — and blocking a receptor is precisely the thing that
    triggers homeostatic feedback raising its ligand. A NAM or a degrader does
    not have this problem: neither can be outcompeted.
    """
    surmountable = [r["mechanism"] for r in rows if r.get("surmountable")]
    insurmountable = [r["mechanism"] for r in rows if not r.get("surmountable")]
    if not surmountable:
        return None

    msg = (f"{', '.join(surmountable)} block(s) here are **competitive**: they "
           "occupy the ligand site and are displaced when endogenous ligand "
           "rises.")
    if ligand_axis:
        msg += (" This receptor has a curated endocrine axis, so feedback will "
                "push its ligand up when you block it — expect the effect to "
                "fade, and to need escalating doses.")
    if insurmountable:
        msg += (f" {', '.join(insurmountable)} cannot be outcompeted and is the "
                "more robust choice where ligand levels are high or variable.")
    else:
        msg += (" No insurmountable option was found here; consider a degrader "
                "or an allosteric modulator if one exists.")
    return msg


def compare_mechanisms(hits: list[ModulatorHit],
                       ligand_axis: bool = False) -> dict:
    """Side-by-side of the mechanism classes available at this target."""
    by_class: dict[str, list[ModulatorHit]] = {}
    for h in hits:
        by_class.setdefault(h.mechanism, []).append(h)

    rows = []
    for mech, group in by_class.items():
        prof = profile_for(mech)
        rows.append({
            "mechanism": mech,
            "n_drugs": len(group),
            "examples": [g.drug for g in group[:3]],
            "max_effect": prof.max_effect,
            "subtype_selectivity": prof.subtype_selectivity,
            "preserves_dynamics": prof.preserves_dynamics,
            "needs_endogenous_ligand": prof.needs_endogenous_ligand,
            "desensitisation_risk": prof.desensitisation_risk,
            "surmountable": prof.surmountable,
            "caveat": prof.caveat,
        })
    rows.sort(key=lambda r: (r["mechanism"] not in ALLOSTERIC, -r["n_drugs"]))

    allosteric = [r for r in rows if r["mechanism"] in ALLOSTERIC]
    summary = None
    if allosteric and len(rows) > len(allosteric):
        summary = (
            "Both allosteric and orthosteric options exist at this target. The "
            "allosteric route usually buys subtype selectivity and preserved "
            "signalling dynamics at the cost of a bounded maximal effect and a "
            "hard dependence on endogenous ligand. Which matters more depends on "
            "whether you need a large fold change or a clean one.")
    elif allosteric:
        summary = ("Only allosteric pharmacology is available here, so effect "
                   "size will be bounded by endogenous tone.")
    elif rows:
        summary = ("No approved allosteric modulator at this target; only "
                   "orthosteric pharmacology, with the usual subtype-selectivity "
                   "problem.")

    return {"by_mechanism": rows, "summary": summary,
            "has_allosteric_option": bool(allosteric),
            "surmountability": surmountability_note(rows, ligand_axis)}
