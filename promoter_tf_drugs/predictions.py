"""The headline answer: which medications the model predicts, and why.

Why this section exists
-----------------------
The evidence for any one drug was scattered across eight sections — the route
catalogue, the grade table and its components, the outcome chain, the
mechanism comparison, the ligand axis, the specificity axes and the
comparison against the standard of care. Reading it meant assembling the
argument yourself, from tables that each answer a different question.

That is the wrong default. The first question anyone brings to this tool is
"what should I give, and why do you think so", and the report answered it
only implicitly, by rank order.

So this states the chain in sentences, per drug:

    what the drug does  ->  to which factor  ->  why that factor is believed
    to act on your promoter  ->  which direction that moves your gene  ->
    how confident, and what would falsify it

What a "prediction" here is, and is not
---------------------------------------
Every row is a **hypothesis generated from annotation**, not a
recommendation. The tier tells you how far you are reaching: tier 1 is
documented in humans, tier 5 is a docking score. The predicted probability
that your gene actually moves is usually around 0.15, because four uncertain
steps multiply — and the section prints that number rather than hiding it
behind a rank.

Nothing here is medical advice. "Approved" means approved for its own
indication; using a drug to nudge an unrelated gene is off-label repurposing
carrying the whole systemic pharmacology with it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import delivery
from . import indications

# How much a tier should move you, in words rather than a number.
TIER_CONFIDENCE = {
    1: "documented in humans for this gene",
    2: "measured experimentally",
    3: "a human curator annotated this mechanism",
    4: "a model predicts it",
    5: "generated, never tested",
}


@dataclass
class Prediction:
    drug: str
    target: str | None = None
    category: str = ""
    mechanism: str = ""
    tier: int = 5
    grade: str | None = None
    score: float = 0.0
    p_move: float | None = None
    p_range: tuple | None = None
    limiting_step: str | None = None
    why: list = field(default_factory=list)
    caveats: list = field(default_factory=list)
    falsifier: str | None = None
    # "yes" / "no" / —. Screening-stage question, so it is carried
    # on every prediction rather than computed only under --tissue cns.
    bbb: str = "—"
    # What the molecule is, and what it is approved to treat. Every route
    # here is off-label by construction, so the approved indication is the
    # baseline the proposal departs from.
    what_it_is: str | None = None
    treats: str | None = None

    def to_dict(self) -> dict:
        return {"drug": self.drug, "target": self.target,
                "bbb": self.bbb, "what_it_is": self.what_it_is,
                "treats": self.treats,
                "category": self.category, "tier": self.tier,
                "grade": self.grade, "score": self.score,
                "p_move": self.p_move, "limiting_step": self.limiting_step,
                "why": self.why, "caveats": self.caveats,
                "falsifier": self.falsifier}


def _direction_clause(role: str | None, goal: str, is_direct: bool) -> str:
    """Why acting on this factor moves the gene the way you want."""
    if is_direct:
        return ("It acts on the gene, its transcript or its protein directly, "
                "so no transcription-factor reasoning is involved.")
    want = "raise" if goal == "up" else "lower"
    if role == "activator":
        act = "agonising" if goal == "up" else "blocking"
        return (f"The factor is annotated as an **activator**, and you want to "
                f"{want} the gene, so {act} it is the right direction.")
    if role == "repressor":
        act = "blocking" if goal == "up" else "agonising"
        return (f"The factor is annotated as a **repressor**, and you want to "
                f"{want} the gene, so {act} it is the right direction.")
    if role == "dual":
        return (f"The factor is annotated as **dual-role** — it activates at "
                f"some promoters and represses at others — so which direction "
                f"it moves your gene is not settled by the annotation.")
    return ("The factor has no activator/repressor annotation, so the "
            "direction of its effect on your gene is assumed, not known.")


def _promoter_clause(tf: dict | None, has_chip: bool) -> str:
    """Why this factor is believed to act on your promoter at all."""
    if not tf:
        return ("No binding-site evidence for this factor in the scanned "
                "window — this route came from elsewhere in the analysis.")
    n = tf.get("n_sites") or 0
    rel = tf.get("best_rel_score")
    off = tf.get("best_site_tss_offset")
    peaks = tf.get("chip_peaks") or 0
    bits = []
    if n:
        where = f"{off:+d} bp from the TSS" if isinstance(off, int) else "in the window"
        bits.append(f"{n} predicted site(s), the best scoring "
                    f"{rel} {where}")
    if peaks:
        bits.append(f"**{peaks} ChIP peak(s) confirm it actually binds here** "
                    f"— the strongest link in this chain")
    elif has_chip:
        bits.append("**no ChIP peak over this promoter**, although ChIP data "
                    "was available — the motif is unconfirmed")
    else:
        bits.append("no ChIP data in this run, so the site is a sequence "
                    "prediction only")
    return "Binding: " + "; ".join(bits) + "."


def _family_clause(target: str | None) -> str:
    """Can a ligand here even be subtype-selective?

    An orthosteric pocket is conserved across a receptor family, so THRA/THRB
    and RARA/RARG cross-react by construction rather than by any drug's
    failing. Worth saying next to the drug rather than only in a panel run.
    """
    if not target:
        return ""
    from . import paralogs as par_mod
    p = par_mod.paralog_specificity(target)
    if not p["assessed"] or not p["relatives"]:
        return ""
    sel = p.get("selective_agents") or []
    base = (f"Family: {target} shares its {p['ligand']} pocket with "
            f"{', '.join(p['relatives'])}, so a ligand here will tend to "
            f"engage them too")
    if sel:
        return base + f" — though subtype selectivity has been achieved ({sel[0]})."
    return base + ", and no subtype-selective agent is recorded for this family."


def _specificity_clause(grade: dict | None, burden: dict | None) -> str:
    parts = []
    ev = (grade or {}).get("evidence") or {}
    size = ev.get("regulon_size")
    if size:
        parts.append(f"the factor controls roughly {size} gene(s), which all "
                     f"move with yours")
    if burden and burden.get("assessed"):
        cls = burden.get("drug_class")
        parts.append(f"clinical burden **{burden.get('level')}**"
                     + (f" ({cls})" if cls else ""))
    elif burden is not None:
        parts.append("clinical burden not assessed — the class table is "
                     "incomplete, so this is not a clean bill of health")
    return "Collateral: " + "; ".join(parts) + "." if parts else ""


def build(block: dict, *, goal: str = "up", has_chip: bool = False,
          limit: int = 10, burdens: dict | None = None) -> list:
    """Assemble the top predicted medications with their reasoning chain."""
    routes = ((block.get("all_routes") or {}).get("routes") or [])
    grades = {g["drug_name"].lower(): g
              for g in (block.get("drug_grades") or [])}
    tfs = {t["symbol"]: t for t in (block.get("tfs") or [])}
    outcomes = {(o["drug"].lower(), o.get("tf_symbol")): o
                for o in (block.get("outcome_estimates") or [])}

    out: list = []
    seen: set = set()
    for r in routes:
        key = (r["drug"].lower(), r.get("acts_on"))
        if key in seen:
            continue
        seen.add(key)

        target = r.get("acts_on")
        tf = tfs.get(target)
        g = grades.get(r["drug"].lower())
        o = outcomes.get(key) or outcomes.get((r["drug"].lower(), target))
        is_direct = r["directness"] == 0

        p = Prediction(
            drug=r["drug"], target=target, category=r["category"],
            mechanism=r.get("mechanism") or "", tier=r["tier"],
            grade=(g or {}).get("grade"), score=r["score"],
        )
        desc = indications.describe(
            r["drug"], mechanism=(r.get("mechanism") or "").strip() or None)
        p.what_it_is, p.treats = desc["what_it_is"], desc["treats"]
        p.bbb = delivery.bbb_verdict(
            r["drug"], ((g or {}).get("evidence") or {}).get("properties")
        )["cell"]
        if o:
            p.p_move = o.get("p_point")
            p.p_range = (o.get("p_low"), o.get("p_high"))
            p.limiting_step = o.get("limiting_step")
            p.falsifier = o.get("next_experiment")

        # ---- the chain, in order
        mech = (r.get("mechanism") or "").strip().rstrip(".")
        if is_direct:
            p.why.append(f"**{r['drug']}** — {mech or 'acts on the gene '
                                                      'product directly'}.")
        elif mech:
            p.why.append(f"**{r['drug']}** acts on **{target}**: {mech}.")
        else:
            p.why.append(f"**{r['drug']}** acts on **{target}**.")
        p.why.append(_direction_clause(
            (tf or {}).get("role"), goal, is_direct))
        if not is_direct:
            p.why.append(_promoter_clause(tf, has_chip))
        spec = _specificity_clause(g, (burdens or {}).get(r["drug"].lower()))
        if spec:
            p.why.append(spec)
        if not is_direct:
            fam = _family_clause(target)
            if fam:
                p.why.append(fam)
        p.why.append(
            f"Evidence tier {r['tier']} — {TIER_CONFIDENCE.get(r['tier'], '')}"
            + (f"; grade {p.grade}" if p.grade else "") + ".")

        # ---- caveats worth carrying to the top of the report
        for note in (r.get("notes") or [])[:2]:
            p.caveats.append(note)
        for note in ((g or {}).get("notes") or [])[:2]:
            if "no potency data" in note or "barely characterised" in note:
                p.caveats.append(note)
            elif "WITHDRAWN" in note.upper() or "black-box" in note:
                p.caveats.append(note)
        if o and o.get("tissue_dependent"):
            p.caveats.append(
                "predicted direction is tissue-dependent — this is the SERM "
                "regime, where the cell type decides the sign, not the molecule")
        out.append(p)
        if len(out) >= limit:
            break
    return out


def render(preds: list, gene: str, goal: str, *, verdict: dict | None = None
           ) -> list:
    """Markdown for the headline section."""
    L: list = []
    verb = "raise" if goal == "up" else "lower"
    L.append(f"## Model-predicted medications for {gene}")
    L.append("")
    if not preds:
        L.append("No route was predicted for this gene. That is a result, not "
                 "an error — check the evidence table above before concluding "
                 "the pharmacology does not exist.")
        L.append("")
        return L

    L.append(f"Compounds this analysis predicts would **{verb}** {gene}, with "
             f"the reasoning for each. Everything here is a **hypothesis "
             f"generated from annotation**, not a recommendation — the tier "
             f"on each row tells you how far the claim is reaching.")
    L.append("")

    lvl = (verdict or {}).get("level")
    if lvl == "not interpretable":
        L.append("> **⛔ Do not act on this list.** The run's evidence verdict "
                 "is *not interpretable* — the candidates below are ranked "
                 "from chance motif matches. Fix the scan first.")
        L.append("")
    elif lvl == "degraded":
        L.append("> **⚠ Read as a reading order, not as evidence.** The run's "
                 "evidence verdict is *degraded*: at least one pillar of the "
                 "analysis is missing. See the evidence table above.")
        L.append("")

    # P(gene moves) is deliberately NOT a column here. Two of its four
    # factors are constants across candidates (direction defaults to 0.75
    # without --predict-direction; responsiveness is fixed at 0.40 because
    # nothing here can measure it), so the spread is ~0.05 and putting it in
    # the rightmost column of a numbered table read as the sort key — while
    # visibly disagreeing with the actual order. It conveys magnitude, not
    # rank, so it is stated once below and per-drug in the detail.
    L.append("| # | Medication | What it is | Approved to treat | "
             "Acts on | BBB? | Tier | Grade | Score |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for i, p in enumerate(preds, 1):
        what = p.what_it_is or "—"
        treats = p.treats or "not established here"
        L.append(f"| {i} | **{p.drug}** | {what} | {treats} | "
                 f"{p.target or '—'} | {p.bbb} | "
                 f"{p.tier} | {p.grade or '—'} | {p.score} |")
    L.append("")
    L.append("**Every row is off-label.** The drug was approved for what "
             "the fourth column says, and is listed here for something "
             "else. That column is the baseline this proposal departs "
             "from, and the distance between it and your gene is part of "
             "what you are judging - a topical acne gel and an infused "
             "cytotoxic can sit on adjacent rows with similar scores.")
    L.append("")
    L.append("*Ordered by evidence tier first, then score. A tier-1 route "
             "with a mediocre score beats a tier-5 route with a great one. "
             "**BBB?** is not part of the ordering — it is an "
             "applicability filter you apply yourself, since whether brain "
             "exposure is wanted or unwanted depends on where your gene needs "
             "to move. See the note under the graded table for what a yes "
             "does and does not establish.*")
    L.append("")

    probs = [p.p_move for p in preds if p.p_move is not None]
    if probs:
        lo, hi = min(probs), max(probs)
        L.append(f"**Every candidate here has a predicted P(gene moves) "
                 f"between {lo:.2f} and {hi:.2f}.** That narrow band is not a "
                 f"coincidence and not a ranking: two of the chain's four "
                 f"steps are the same constant for every row — direction "
                 f"defaults to 0.75 without `--predict-direction`, and "
                 f"responsiveness is fixed at 0.40 because nothing in this "
                 f"package can measure it. So the number tells you how "
                 f"uncertain this *class* of route is (marginal), not which "
                 f"of these is better. Order on tier and score; use P to "
                 f"decide whether to pursue a promoter route at all.")
        L.append("")

    # Grouped by TARGET, not by drug. Every drug acting through the same
    # factor shares its binding evidence, its collateral and its receptor
    # family word for word - in a ten-row section that repeated the same
    # falsifying experiment ten times, the tier line eight times and the
    # direction clause seven. Saying the shared thing once and the
    # distinguishing thing per drug loses nothing and removes about half
    # the section.
    SHARED_PREFIXES = ("The factor is annotated", "The factor has no",
                       "Binding:", "Collateral:", "Family:", "Evidence tier")

    def _m(q):
        """The drug's own mechanism clause, for deciding if it distinguishes."""
        c = next((x for x in q.why if not x.startswith(SHARED_PREFIXES)), "")
        return c.split(": ", 1)[-1] if ": " in c else c

    by_target: dict = {}
    for p in preds:
        by_target.setdefault(p.target or "—", []).append(p)

    L.append("### Why each one")
    L.append("")

    tiers = {p.tier for p in preds}
    if len(tiers) == 1:
        t = tiers.pop()
        L.append(f"Every route below is **evidence tier {t}** — "
                 f"{TIER_CONFIDENCE.get(t, '')}.")
        L.append("")
    # The direction rule is one sentence that applies to every group; only
    # each factor's ROLE differs, so the rule is stated once and the groups
    # carry the role.
    roles = {c.split("**")[1] for p in preds for c in p.why
             if c.startswith("The factor is annotated") and "**" in c}
    if len(roles) >= 1:
        want = "raise" if goal == "up" else "lower"
        L.append(f"**The direction rule**, which applies to every group "
                 f"below: you want to {want} the gene, so you agonise an "
                 f"activator and block a repressor"
                 + (" — and a dual-role factor settles neither way, because "
                    "it activates at some promoters and represses at others."
                    if "dual-role" in roles else ".") )
        L.append("")

    falsifiers = {p.falsifier for p in preds if p.falsifier}
    if len(falsifiers) == 1:
        L.append(f"**One experiment falsifies all of them:** "
                 f"{falsifiers.pop()}")
        L.append("")

    n = 0
    for target, group in by_target.items():
        group.sort(key=lambda p: -p.score)
        head = group[0]
        label = (f"Through **{target}**" if target != "—"
                 else "**Acting on the gene directly**")
        L.append(f"#### {label}"
                 + (f" — {len(group)} candidates" if len(group) > 1 else ""))
        L.append("")

        # the factor-level argument, stated once for the whole group
        for clause in head.why:
            if not clause.startswith(SHARED_PREFIXES):
                continue
            if clause.startswith("Evidence tier") and len(tiers) == 0:
                continue
            if clause.startswith("The factor is annotated") and roles:
                # the rule is stated above; only the role is news here
                role = clause.split("**")[1] if "**" in clause else "?"
                L.append(f"- Role: **{role}**")
                continue
            L.append(f"- {clause}")
        L.append("")

        # caveats shared by every drug in the group belong to the factor
        group_caveats = set()
        if len(group) > 1:
            common = [set(q.caveats or []) for q in group]
            group_caveats = set.intersection(*common) if common else set()
            for c in sorted(group_caveats):
                L.append(f"- ⚠ {c}")
            if group_caveats:
                L.append("")

        # then only what distinguishes each drug
        for p in group:
            n += 1
            bits = []
            if p.grade:
                bits.append(f"grade {p.grade}")
            if p.p_move is not None:
                bits.append(f"P={p.p_move:.2f}")
            # Only print the mechanism when it distinguishes this drug from
            # its group-mates; "ligand at the transcription factor" on every
            # row of a group is the same noise in a different place.
            mech = next((c for c in p.why
                         if not c.startswith(SHARED_PREFIXES)), "")
            mech = mech.split(": ", 1)[-1].rstrip(".") if ": " in mech else ""
            if len({_m(q) for q in group}) <= 1:
                mech = ""
            line = f"{n}. **{p.drug}**"
            if bits:
                line += f" — {', '.join(bits)}"
            if mech:
                line += f" — {mech}"
            L.append(line)
            # A caveat every drug in the group shares is a property of the
            # factor, and Collateral has already said it above.
            if p.caveats and p.caveats[0] not in group_caveats:
                L.append(f"   - ⚠ {p.caveats[0]}")
            if p.falsifier and len(falsifiers) > 1:
                L.append(f"   - *Would falsify it:* {p.falsifier}")
        L.append("")

    L.append("*Not medical advice. \"Approved\" means approved for its own "
             "indication — using one of these to move an unrelated gene is "
             "off-label repurposing with the full systemic pharmacology "
             "attached. The experimental sections at the end of this report "
             "are how you would actually test any of it.*")
    L.append("")
    return L
