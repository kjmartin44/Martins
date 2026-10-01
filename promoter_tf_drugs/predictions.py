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

    def to_dict(self) -> dict:
        return {"drug": self.drug, "target": self.target,
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
    L.append("| # | Medication | Acts on | Tier | Grade | Score |")
    L.append("|---|---|---|---|---|---|")
    for i, p in enumerate(preds, 1):
        L.append(f"| {i} | **{p.drug}** | {p.target or '—'} | {p.tier} | "
                 f"{p.grade or '—'} | {p.score} |")
    L.append("")
    L.append("*Ordered by evidence tier first, then score. A tier-1 route "
             "with a mediocre score beats a tier-5 route with a great one.*")
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

    L.append("### Why each one")
    L.append("")
    for i, p in enumerate(preds, 1):
        L.append(f"**{i}. {p.drug}**"
                 + (f" — via {p.target}" if p.target else ""))
        L.append("")
        for clause in p.why:
            L.append(f"- {clause}")
        if p.p_move is not None:
            rng = ""
            if p.p_range and p.p_range[0] is not None:
                rng = f" (range {p.p_range[0]}–{p.p_range[1]})"
            L.append(f"- Predicted probability the gene actually moves: "
                     f"**{p.p_move:.2f}**{rng}"
                     + (f", limited by *{p.limiting_step}*."
                        if p.limiting_step else "."))
        if p.caveats:
            L.append(f"- ⚠ {p.caveats[0]}")
        if p.falsifier:
            L.append(f"- *Would falsify it:* {p.falsifier}")
        L.append("")

    L.append("*Not medical advice. \"Approved\" means approved for its own "
             "indication — using one of these to move an unrelated gene is "
             "off-label repurposing with the full systemic pharmacology "
             "attached. The experimental sections at the end of this report "
             "are how you would actually test any of it.*")
    L.append("")
    return L
