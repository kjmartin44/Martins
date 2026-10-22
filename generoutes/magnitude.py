"""How much the gene has to move, and which routes have ever been measured.

Every other module answers "would this move the gene". This one asks the
question that decides whether an answer is useful: **how much**, and is
any candidate known to deliver that much.

The report has always declined to predict a fold change, and that
remains right - deriving one needs the share of transcription the factor
drives, its occupancy change at the achieved dose, and the dose reached
in the tissue, none of which this package has. But two things ARE
available and were not being used.

**The requirement is arithmetic, not prediction.** If one allele is dead
and the other works, restoring the normal total means roughly doubling
what the working copy makes. If a toxic protein needs lowering by half,
that is a stated target. Neither is a guess; both come from the
mechanism, and stating them turns "here are two hundred routes" into
"here is what any of them would have to achieve".

**Some routes carry a measured effect and most do not.** A LINCS hit is
a compound that demonstrably moved this gene in a real experiment, with
a signed z-score. An annotation-derived route has no magnitude at all.
That difference is worth surfacing, because a route with a measured
effect of the right sign and size is a different proposition from one
inferred through three steps.

What this will NOT do
---------------------
Convert a z-score into a fold change. L1000's z is computed against the
plate population, so the conversion needs that population's spread in
log space - which varies by plate, cell line and time point and is not
in the data this package reads. A z of -2 is a real effect and is larger
than a z of -1; it is not "half as much protein". Anything claiming
otherwise would be inventing the number that matters most.

So magnitudes are reported, ordered and compared with each other, and
never against the requirement in the same units. Saying "this is the
largest measured effect available, and we cannot tell you whether it is
enough" is the honest shape of this answer.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict

# |z| at or above this is the conventional line for a real L1000 effect
# rather than plate noise. It orders effects; it does not size them.
Z_REAL = 2.0


@dataclass
class Ask:
    """What the gene has to do, stated from the mechanism."""
    gene: str
    goal: str                       # "up" | "down"
    mechanism: str = "unstated"
    requirement: str = ""           # the ask, in words
    target_fold: float | None = None     # multiple of current, when known
    acts_on: str = "transcript"     # "transcript" | "protein"
    # False when no route this package ranks can address the scenario.
    # Two of them cannot, and handing over two hundred routes for a
    # problem none of them solves is worse than saying so.
    feasible: bool = True
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def restore_one_allele(gene: str, *, broken_allele_makes: str = "unknown"
                       ) -> Ask:
    """One copy dead, raise the other to make up the shortfall.

    The arithmetic is simple and worth stating: the working allele
    currently supplies about half the normal total, so it has to roughly
    double. That is a small ask by the standards of transcriptional
    pharmacology, which is the useful part - it is far easier than
    asking for a level the gene has never had.
    """
    a = Ask(gene=gene, goal="up",
            mechanism="one allele lost, raising the other",
            requirement="about 2x from the working allele, to restore "
                        "the normal total",
            target_fold=2.0)
    a.notes.append(
        "The target is a return to the normal level, not an increase "
        "above it. Every healthy heterozygous carrier is evidence that "
        "the destination is tolerated; nothing here says the same about "
        "overshooting it.")
    if broken_allele_makes == "truncated":
        a.notes.append(
            "Raising transcription raises BOTH copies, so this also "
            "raises the truncated product. The requirement above is "
            "therefore not the whole specification: a route that "
            "doubles total transcript does not deliver a doubling of "
            "functional protein.")
    elif broken_allele_makes == "null":
        a.notes.append(
            "The broken copy makes nothing, so a doubling of total "
            "transcript is a doubling of functional protein.")
    return a


def lower_toxic_product(gene: str, *, target_fraction: float = 0.5,
                        acts_on: str = "protein",
                        wild_type_essential: bool = True) -> Ask:
    """Lower a harmful gene product to a stated fraction of current.

    `target_fraction` is the share you want left: 0.5 is a halving.
    Stated rather than derived, because the therapeutic target for a
    lowering strategy comes from the biology and the clinic, not from
    this package.

    `wild_type_essential` matters more than the number. Where the
    normal protein is needed, a non-selective lowering has a floor
    below which the treatment becomes the problem, and the useful
    window may be narrow.
    """
    pct = int(round((1 - target_fraction) * 100))
    a = Ask(gene=gene, goal="down",
            mechanism="lowering a harmful gene product",
            requirement=f"about a {pct}% reduction, leaving roughly "
                        f"{int(round(target_fraction * 100))}% of current",
            target_fold=target_fraction, acts_on=acts_on)
    if acts_on == "protein":
        a.notes.append(
            "The target is PROTEIN, and everything this package ranks "
            "acts on transcription. The two come apart: a long-lived "
            "protein lags a lowered transcript by its own half-life, "
            "and a measured transcript change is an upper bound on what "
            "you would see in a western blot at the same time point.")
    if wild_type_essential:
        a.notes.append(
            "The normal protein is needed, so this has a floor as well "
            "as a target. A route that lowers both the harmful and the "
            "normal product has a window rather than a direction, and "
            "nothing here predicts where the bottom of it is. An "
            "allele-selective approach removes that constraint; a "
            "promoter-level one cannot.")
    return a


def measured_effects(routes, gene: str) -> dict:
    """Split routes by whether anything measured their effect size.

    Returns the ones carrying a measured magnitude, ordered largest
    first, and a count of those that carry none. The split is the
    point: most routes are annotation, and an annotation has no size.
    """
    want = (gene or "").strip().upper()
    measured, unmeasured = [], 0
    for r in routes or []:
        # A measured effect records which gene it was measured on. A hit
        # for a different gene is not evidence about this one, and
        # counting it would inflate the one number this exists to
        # report. Routes that name no gene are taken at face value.
        on = (r.get("gene") or "").strip().upper()
        if want and on and on != want:
            continue
        score = r.get("lincs_score")
        if score is None:
            score = r.get("score_measured")
        if score is None and str(r.get("source") or "") == "lincs":
            score = r.get("score")
        if score is None:
            unmeasured += 1
            continue
        try:
            z = float(score)
        except (TypeError, ValueError):
            unmeasured += 1
            continue
        measured.append({
            "drug": r.get("drug"),
            "z": round(z, 2),
            "direction": "down" if z < 0 else "up",
            "cell_line": r.get("cell_line"),
            "real": abs(z) >= Z_REAL,
        })
    measured.sort(key=lambda m: -abs(m["z"]))
    return {"measured": measured, "n_measured": len(measured),
            "n_unmeasured": unmeasured}


def summarise(ask: Ask, effects: dict | None = None) -> str:
    """The requirement, and how far the evidence goes toward it."""
    if not ask.feasible:
        return " ".join(
            [f"**Changing transcription cannot solve this for "
             f"{ask.gene}: {ask.requirement}.**"] + ask.notes
            + ["Everything ranked below acts on how much the gene is "
               "transcribed, so none of it addresses the problem as "
               "stated. Read it as context rather than as candidates."])

    L = [f"**What {ask.gene} has to do: {ask.requirement}.**"]
    L.extend(ask.notes)

    e = effects or {}
    n_m, n_u = e.get("n_measured", 0), e.get("n_unmeasured", 0)
    if not (n_m or n_u):
        return " ".join(L)

    if n_m:
        best = e["measured"][0]
        wrong_way = [m for m in e["measured"]
                     if (m["direction"] != ask.goal)]
        L.append(
            f"Of {n_m + n_u} route(s), **{n_m} carry a measured effect** "
            f"on this gene and {n_u} carry none - those rest on "
            f"annotation, which has no size.")
        L.append(
            f"The largest measured effect is {best['drug']} at z="
            f"{best['z']}"
            + (f" in {best['cell_line']}" if best.get("cell_line") else "")
            + (", which clears the conventional line for a real effect."
               if best["real"] else
               ", which is below the conventional line for a real effect "
               "and may be plate noise."))
        if wrong_way:
            L.append(
                f"{len(wrong_way)} of the measured effects move the gene "
                f"the WRONG way for this goal.")
    else:
        L.append(
            f"**No route carries a measured effect on this gene.** All "
            f"{n_u} rest on annotation, so nothing here has a size to "
            f"compare against the requirement.")

    if n_m:
        # Only when there are numbers to caveat. "These numbers" reads
        # oddly under a sentence saying there are none.
        L.append(
            "A z-score cannot be converted to a fold change: it is "
            "computed against the plate population, and that "
            "population's spread is not in the data read here. These "
            "numbers order the measured effects against each other. "
            "They cannot tell you whether the largest of them is "
            "enough.")
    return " ".join(L)


# ===================================================================== #
# The shapes a dosage problem comes in.                                 #
# ===================================================================== #
# Declared by the caller, never inferred. The package cannot read a
# variant, and guessing the mechanism from a gene symbol would be the
# worst kind of help.
SCENARIOS = (
    "haploinsufficiency",   # one allele dead, raise the other
    "hypomorph",            # both alleles weak, raise both
    "biallelic_hypomorph",  # one weak allele, one dead - raise the weak one
    "toxic_gain",           # a harmful product to lower
    "null",                 # no functional transcript from either allele
    "silenced_allele",      # a working copy switched off
    "modifier",             # no mutation; move a normal gene to help
)

# The direction each scenario implies, and why. A scenario and --direction
# are declared separately, so nothing stopped a contradictory pair: the
# report opened "about 2x from the working allele" above a ranking built to
# LOWER the gene, led by a transcription inhibitor. `null` and `modifier`
# imply no direction (a null cannot be helped either way, and a modifier is
# moved whichever way helps), so they accept both.
SCENARIO_GOAL = {
    "haploinsufficiency": ("up", "one allele is lost, so the plan is to "
                                 "raise the working one"),
    "hypomorph": ("up", "both alleles are weak, so the plan is to raise "
                        "expression"),
    "biallelic_hypomorph": ("up", "the only working copy is a weak one, so "
                                  "the plan is to raise it"),
    "silenced_allele": ("up", "an intact copy is switched off, so the plan "
                              "is to switch it on"),
    "toxic_gain": ("down", "the product is harmful, so the plan is to "
                           "lower it"),
}


def goal_conflict(scenario: str | None, goal: str | None) -> str | None:
    """A sentence naming the contradiction, or None when the pair is fine.

    Checked before any request goes out: the combination is knowable from
    the two arguments alone, and learning of it from the report means a
    run spent ranking the opposite of what was asked.
    """
    need = SCENARIO_GOAL.get(scenario or "")
    if need is None or goal is None or goal == need[0]:
        return None
    return (f"scenario {scenario!r} needs direction {need[0]!r} - "
            f"{need[1]} - but direction {goal!r} was asked for. The "
            f"requirement and the ranking would point opposite ways.")


def ask_for(scenario: str, gene: str, **kw) -> Ask:
    """The requirement for one of the recognised scenarios.

    Raises for an unknown scenario rather than defaulting, because a
    silent default here would set the direction of the whole analysis.
    """
    if scenario not in SCENARIOS:
        raise ValueError(
            f"unknown scenario {scenario!r}. One of: "
            + ", ".join(SCENARIOS))
    return _SCENARIO_BUILDERS[scenario](gene, **kw)


def _biallelic_hypomorph(gene: str, *, target_fold: float = 2.0,
                         residual_function: float | None = None,
                         broken_allele_makes: str = "unknown", **_) -> Ask:
    """One allele hypomorphic, the other dead. The recessive arrangement.

    Distinct from `haploinsufficiency`, which assumes the surviving copy is
    WILD TYPE and only needs doubling. Here the surviving copy is itself
    impaired, so transcription is working against a lower ceiling: doubling
    it doubles a reduced-function product, and whether that clears the
    threshold depends on how much function the allele retains - which is a
    property of the variant and is not in any database this package reads.

    It is the usual arrangement in recessive disease, and declaring
    haploinsufficiency instead quietly promises a rescue the arithmetic
    does not support.
    """
    a = Ask(gene=gene, goal="up", mechanism="one weak allele, one dead",
            requirement=f"about {target_fold:g}x from the weak allele",
            target_fold=target_fold)
    a.notes.append(
        "Only one copy makes anything, and it is impaired. Raising "
        "transcription raises the amount of that impaired product; it does "
        "not improve it. Whether more of it is enough depends on how much "
        "function the allele keeps, which no database here holds - so the "
        "fold above is what transcription could deliver, not a dose that is "
        "known to be sufficient.")
    if residual_function is not None:
        reached = target_fold * residual_function
        a.notes.append(
            f"At {residual_function:.0%} residual function per transcript, "
            f"{target_fold:g}x of it reaches about {reached:.0%} of normal "
            f"activity."
            + ("" if reached >= 0.5 else
               " That is below the half a healthy carrier demonstrates is "
               "tolerated, so transcription alone would not close the gap."))
    if broken_allele_makes == "truncated":
        a.notes.append(
            "The dead allele is recorded as making a truncated product. "
            "Raising transcription raises BOTH copies, so that product "
            "rises alongside the one you want. Where it competes with the "
            "working protein, more transcription can make matters worse.")
    elif broken_allele_makes == "unknown":
        a.notes.append(
            "What the dead allele produces is not declared "
            "(--broken-allele-makes). If it makes a truncated protein "
            "rather than nothing, raising transcription amplifies it too.")
    a.notes.append(
        "If the problem is the quality of the protein rather than its "
        "amount, the routes that address it are stabilisation, readthrough "
        "and the protein-level sections, not the promoter.")
    return a


def _hypomorph(gene: str, *, target_fold: float = 2.0, **_) -> Ask:
    """Both copies present and weak. Raising both is the whole point."""
    a = Ask(gene=gene, goal="up", mechanism="both alleles hypomorphic",
            requirement=f"about {target_fold:g}x from both alleles",
            target_fold=target_fold)
    a.notes.append(
        "Both copies make a working, if reduced, product, so raising "
        "transcription raises the thing you want and nothing you do "
        "not. This is the most favourable shape for everything below: "
        "no allele-selectivity problem, and no bad product to amplify.")
    a.notes.append(
        "How much is enough depends on the residual function of the "
        "variant protein, which is a property of the variant and not of "
        "the gene. Nothing here can supply it.")
    return a


def _null(gene: str, **_) -> Ask:
    """Nothing to raise. The honest answer is that this will not work."""
    a = Ask(gene=gene, goal="up", mechanism="no functional product",
            requirement="not achievable by changing transcription",
            feasible=False)
    a.notes.append(
        "Neither allele makes a working product, so there is nothing "
        "for a promoter-level or factor-level route to amplify. "
        "Transcribing a null allele harder produces more of the same "
        "non-functional transcript.")
    a.notes.append(
        "The routes this package ranks cannot address this. What can: "
        "supplying the gene or protein from outside, reading through a "
        "premature stop where one exists, or correcting the sequence. "
        "`--design-modality` describes the first of those.")
    return a


def _silenced(gene: str, **_) -> Ask:
    """A working copy that is switched off, which is a different job."""
    a = Ask(gene=gene, goal="up", mechanism="a working allele silenced",
            requirement="re-express a silent but intact copy")
    a.notes.append(
        "The sequence is intact and the regulation is the problem, so "
        "this is not a dosage question at all - you are not asking an "
        "active promoter for more, you are asking a silent one to "
        "start. The machinery differs: chromatin and methylation "
        "rather than factor occupancy.")
    a.notes.append(
        "Agents that do this work genome-wide almost by definition, so "
        "the specificity axes matter more here than anywhere else in "
        "this package, not less.")
    return a


def _modifier(gene: str, *, goal: str = "down", condition: str = "",
              **_) -> Ask:
    """No mutation. A normal gene moved because moving it helps."""
    a = Ask(gene=gene, goal=goal,
            mechanism=f"modifier of {condition}" if condition
                      else "modifier, no mutation in this gene",
            requirement="no restoration target - the amount needed is a "
                        "clinical question, not an arithmetic one")
    a.notes.append(
        "There is no normal level to return to, because the gene is "
        "already at it. That removes the arithmetic that makes a "
        "rescue strategy checkable: nothing here can say how far is "
        "far enough, and the honest form of the answer is a direction "
        "with no number.")
    a.notes.append(
        "It also inverts the safety argument. A rescue aims at a level "
        "every healthy person demonstrates is tolerated; this moves a "
        "normal gene away from that level in someone whose copy works. "
        "The collateral axes - how many genes the factor controls, how "
        "body-wide it is, what else the drug moves - carry more weight "
        "here than in any rescue scenario.")
    return a


_SCENARIO_BUILDERS = {
    "haploinsufficiency": lambda g, **kw: restore_one_allele(
        g, broken_allele_makes=kw.get("broken_allele_makes", "unknown")),
    "hypomorph": _hypomorph,
    "biallelic_hypomorph": lambda g, **kw: _biallelic_hypomorph(
        g,
        target_fold=kw.get("target_fold", 2.0),
        residual_function=kw.get("residual_function"),
        broken_allele_makes=kw.get("broken_allele_makes", "unknown")),
    "toxic_gain": lambda g, **kw: lower_toxic_product(
        g,
        target_fraction=kw.get("target_fraction", 0.5),
        acts_on=kw.get("acts_on", "protein"),
        wild_type_essential=kw.get("wild_type_essential", True)),
    "null": _null,
    "silenced_allele": _silenced,
    "modifier": _modifier,
}
# Below this, losing one copy is not selected against in the population, so
# a dominant one-allele-lost model has to argue against the frequency data.
PLI_TOLERANT = 0.1


def constraint_conflict(scenario: str | None, pli: float | None,
                        loeuf: float | None = None) -> str | None:
    """Does the population data argue against the declared scenario?

    The scenario is declared and the constraint is measured, so only a
    check can notice when they disagree - the same reasoning as
    `goal_conflict`, against a different pair. This one cannot be run
    before the analysis, because the constraint tables are read during it,
    so it is reported rather than refused: the declaration is the
    clinician's and a gene can be haploinsufficient in one tissue and
    tolerated in the population.

    BUB1B is the case it was written for: pLI about 0, because MVA is
    recessive and heterozygous carriers are healthy, declared as
    haploinsufficiency - which promises that doubling a wild-type allele
    fixes it, when the second allele is usually a hypomorph.
    """
    if scenario != "haploinsufficiency" or pli is None or pli >= PLI_TOLERANT:
        return None
    extra = (f" (LOEUF {loeuf:.2f})" if isinstance(loeuf, (int, float))
             else "")
    return (
        f"**The declared scenario and the population data disagree.** "
        f"`--scenario haploinsufficiency` says one allele is lost and the "
        f"other is intact. But this gene's pLI is {pli:.3g}{extra}: losing "
        f"one copy is not selected against, so healthy people carry a "
        f"broken copy and a single loss is unlikely to cause disease on its "
        f"own. That usually means the condition is recessive and the second "
        f"allele is not wild type either - in which case the requirement "
        f"above understates the ask, because it assumes the surviving copy "
        f"works normally. Consider `--scenario biallelic_hypomorph`, and "
        f"declare `--broken-allele-makes`.")
