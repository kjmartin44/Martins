"""A plain sentence on how a candidate would act on the condition.

Why this exists
---------------
A ranking says a drug is a candidate and leaves the reader to work out why. The
route category is a code ("TF: abundance change") and the mechanism text is
written for the catalogue, not for someone deciding what to read first. This
turns a route into one sentence that says what the drug acts on, how that would
reach the gene, and which step is not established.

What it does
------------
Reads only fields the route already carries (its category, the target it acts
on, its mechanism text, its evidence label) and the factor's role at the
promoter, and fills a fixed sentence for each kind of route. It states a
mechanism the pipeline predicted. It adds no biology of its own.

`aim` gives the one paragraph that ties the list to the condition: what the
gene has to do (from the dosage scenario if one was declared) and the caveat
that applies to every row.

What it does not do
-------------------
**It is not evidence.** A sentence reads as though the mechanism were known. It
is the mechanism the route assumes, and for most routes the step from the drug's
target to the gene has not been measured. The sentence says so where the route
says so.

**It does not say the drug treats the condition.** It says how the drug would
move the gene. Whether moving the gene helps is a separate claim, which the
dosage scenario states and nothing here tests.
"""

from __future__ import annotations

from .drugs import desired_direction

_SCENARIO_HELP = {
    "haploinsufficiency": "one copy is lost, so the aim is to raise the working "
                          "copy until the total is back to normal",
    "hypomorph": "both copies work weakly, so the aim is to raise them toward "
                 "the normal total",
    "toxic_gain": "the harmful product comes from the gene, so the aim is to "
                  "lower it, not remove it",
    "null": "no working copy remains, so raising transcription cannot produce "
            "a working product and these routes are context",
    "silenced_allele": "an intact copy is switched off, so the aim is to switch "
                       "it back on",
    "modifier": "there is no mutation in this gene; the aim is to move a normal "
                "gene in a direction that helps the condition",
}


def _verb(goal: str) -> str:
    return "raise" if goal == "up" else "lower"


def _role_phrase(role: str | None) -> str:
    return {"activator": "switch on", "repressor": "hold down"}.get(
        role or "", "act on")


def _action_on_tf(role: str | None, goal: str) -> str:
    want = desired_direction(role or "unknown", goal)
    return {"up": "activating", "down": "blocking"}.get(
        want, "changing the activity of")


def aim(gene: str, goal: str, requirement: dict | None = None,
        condition: str | None = None, scenario: str | None = None) -> str:
    """The paragraph that ties the list to the condition.

    The dosage requirement from the run is used when there is one; failing
    that, the declared scenario is put in a sentence; failing that, the goal.
    """
    verb = _verb(goal)
    parts = []
    if requirement and requirement.get("requirement"):
        parts.append(
            f"For {gene}: {requirement.get('mechanism') or 'the gene'}. The "
            f"target is {requirement['requirement']}.")
    elif scenario in _SCENARIO_HELP:
        parts.append(f"For {gene}, {_SCENARIO_HELP[scenario]}.")
    else:
        what = f" to help with {condition}" if condition else ""
        parts.append(f"The aim is to {verb} {gene}{what}.")
    parts.append(
        "Each sentence below says how that drug would move the gene. That is "
        "the mechanism the route assumes, and for most routes the step from "
        "the drug's target to the gene has not been measured. Nothing here "
        "shows that moving the gene treats the condition.")
    return " ".join(parts)


def route_sentence(route: dict, gene: str, goal: str,
                   tf_roles: dict | None = None) -> str:
    """One sentence on how a route would act, from its own fields."""
    verb = _verb(goal)
    cat = route.get("category") or ""
    tf = route.get("acts_on") or "its target"
    role = (tf_roles or {}).get(tf)
    mech = (route.get("mechanism") or "").strip().rstrip(".")
    measured = (route.get("axes") or {}).get("evidence") == "measured binding"

    if cat in ("TF: orthosteric ligand", "TF: allosteric modulator"):
        untested = (f"Its binding is measured; the step from {tf} to {gene} "
                    f"is not." if measured else
                    f"Neither its effect on {tf} nor the step from {tf} to "
                    f"{gene} is measured here.")
        if role not in ("activator", "repressor", "dual"):
            return (f"Binds {tf}, a factor with a site at the {gene} "
                    f"promoter whose role there is not known, so whether to "
                    f"activate or block it to {verb} {gene} is not known. "
                    f"{untested}")
        return (f"Binds {tf}, a factor predicted to {_role_phrase(role)} "
                f"{gene}; {_action_on_tf(role, goal)} {tf} should {verb} "
                f"{gene}. {untested}")
    if cat == "TF: endogenous ligand":
        return (f"Changes the supply of {tf}'s natural ligand, which shifts "
                f"{tf} activity and so should {verb} {gene}.")
    if cat == "TF: abundance change":
        return (f"Changes how much {tf} the cell makes, which should {verb} "
                f"{gene} through the sites {tf} has in its promoter.")
    if cat == "TF: inferred activity change":
        return (f"Shifts the activity of {tf} in cells, inferred from "
                f"signatures, which should {verb} {gene}.")
    if cat in ("TF: degradation blocker", "TF: stabiliser or glue"):
        return (f"Changes how long {tf} lasts in the cell"
                f"{f' ({mech})' if mech else ''}, which should {verb} {gene}.")
    if cat == "upstream signalling":
        return (f"Acts on {tf}, upstream of factors that regulate {gene}; the "
                f"effect reaches {gene} only by inference through the "
                f"network.")
    if cat == "empirical signature":
        return (f"Was seen to {verb} {gene} in cell-line signatures at lab "
                f"doses; how it does so is not known.")
    if cat == "direct: gene product":
        return (f"Acts on the {gene} protein itself"
                f"{f' ({mech})' if mech else ''}, changing its activity and "
                f"not its amount.")
    if cat == "direct: gene or transcript":
        return (f"Acts on the {gene} gene or transcript"
                f"{f' ({mech})' if mech else ''}.")
    if cat == "direct: predicted product binder":
        return (f"Predicted to bind the {gene} protein; whether it blocks or "
                f"activates it is not known.")
    if cat == "chromatin-level":
        return (f"Changes chromatin across the genome, so it can move {gene} "
                f"together with thousands of other genes.")
    if cat == "promoter DNA binder":
        return (f"Predicted to bind the {gene} promoter DNA directly, which "
                f"could {verb} transcription; unproven.")
    if cat == "virtual screen hypothesis":
        return (f"Predicted by a model to bind {tf}; never tested at that "
                f"receptor.")
    return (f"{mech[:140]}." if mech else
            f"Acts through {tf}; the route category is {cat or 'unstated'}.")


def docking_sentence(target: str, gene: str) -> str:
    """For a drug found only by docking."""
    if (target or "").upper() == (gene or "").upper():
        return (f"Docked into the {gene} protein pocket and shaped to fit it; "
                f"whether it blocks or activates the protein is not known.")
    return (f"Docked into the pocket of {target}, a factor at the {gene} "
            f"promoter; whether it activates or blocks {target} is not known, "
            f"and the step from {target} to {gene} is not measured.")
