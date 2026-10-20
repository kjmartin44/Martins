"""Find interventions that plausibly BEAT what already exists for this gene.

This module answers a different question from the rest of the package. Ranking
candidates tells you what to read first. "Stronger, more specific, or more
efficient" is a *comparative* claim, and a comparative claim needs three things
the ranking does not provide:

  1. a baseline    what currently controls this gene's expression at all
  2. axes          what "better" means, separately, so a win on one axis is not
                   laundered into a win overall
  3. a novelty test  whether the drug-gene link has already been reported

Why Pareto rather than one score
--------------------------------
Collapsing strength, specificity and efficiency into a weighted sum lets a
candidate that is slightly stronger and drastically less specific come out
ahead. That is precisely the trade this pipeline exists to avoid. So candidates
are compared by **dominance**: a candidate is interesting if it beats the
incumbent on at least one axis while not being meaningfully worse on any other.
The non-dominated set is reported, and the axis it wins on is named.

The asymmetry that matters most
-------------------------------
An incumbent drug's numbers are *measured*. A novel candidate's are *predicted*
by the models in this package. Comparing them directly is apples to oranges and
biased toward the challenger, because prediction error is symmetric but
selection is not: you are looking at the candidate precisely because its
predicted score came out high, so its score is inflated by the same noise that
made you notice it.

The correction is a required margin: a predicted candidate must beat a measured
incumbent by more than the prediction error before it is worth chasing. That is
implemented, not just described, and it is why some candidates are reported as
"within noise of the incumbent" rather than as improvements.

Novelty source
--------------
CTD (Comparative Toxicogenomics Database) curates chemical-gene interactions
with explicit "increases expression of" / "decreases expression of" actions -
exactly the relation this pipeline predicts. A candidate already in CTD for
your gene is not a new link; it is a confirmation, which is still useful and is
labelled as such.

  https://ctdbase.org/  (batch query API, no key required)
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, asdict, field

from .http import Client, HttpError

log = logging.getLogger(__name__)

CTD_BATCH = "https://ctdbase.org/tools/batchQuery.go"

# A predicted candidate must beat a measured incumbent by this much on an axis
# before the difference is treated as real rather than as model noise.
PREDICTION_MARGIN = 0.15
# Tolerance for "not meaningfully worse" on the other axes.
TOLERANCE = 0.05


@dataclass
class Axes:
    """The three things 'better' can mean, kept separate on purpose."""
    strength: float          # how hard it pushes the gene
    specificity: float       # how little else it touches
    efficiency: float        # effect per unit of dose, exposure and risk
    evidence: str = "predicted"     # "measured" | "predicted"

    def to_dict(self) -> dict:
        return asdict(self)

    def as_tuple(self) -> tuple:
        return (self.strength, self.specificity, self.efficiency)


@dataclass
class Candidate:
    name: str
    route: str               # "direct" | "via TF"
    tf_symbol: str | None
    axes: Axes
    previously_linked: bool = False
    # Whether the question could be asked at all. False means the CTD
    # lookup failed, so `previously_linked: False` is "unknown" rather
    # than "nobody has tried this".
    novelty_assessed: bool = True
    link_source: str | None = None
    wins_on: list = field(default_factory=list)
    loses_on: list = field(default_factory=list)
    verdict: str = ""
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["axes"] = self.axes.to_dict()
        return d


# ------------------------------------------------------------ axis scoring
def axes_from_grade(grade, tf_result=None, outcome=None) -> Axes:
    """Build the three axes from what the pipeline already computed."""
    c = grade.components

    # Strength: potency x unambiguous pharmacology x how well the TF drives
    # this promoter. A partial agonist is genuinely weaker, not just less certain.
    tf_conf = grade.tf_confidence
    strength = (0.45 * c.get("potency", 0.5)
                + 0.30 * c.get("action", 0.5)
                + 0.25 * tf_conf)
    if outcome is not None:
        strength = 0.6 * strength + 0.4 * min(1.0, outcome.p_point * 3)

    # Specificity: both layers, multiplied. One protein out of many is not
    # specific if that protein runs five thousand genes.
    specificity = c.get("selectivity", 0.5) * c.get("regulon", 0.55)

    # Efficiency: effect per unit dose and risk. A potent, oral, clean drug
    # achieves the same effect with less of everything.
    efficiency = (0.40 * c.get("potency", 0.5)
                  + 0.25 * c.get("maturity", 0.5)
                  + 0.20 * c.get("tissue", 1.0)
                  + 0.15 * c.get("evidence_depth", 0.5))
    ev = grade.evidence or {}
    if ev.get("withdrawn"):
        efficiency *= 0.3
    if ev.get("black_box"):
        efficiency *= 0.7

    return Axes(round(min(1.0, strength), 3),
                round(min(1.0, specificity), 3),
                round(min(1.0, efficiency), 3),
                evidence="predicted")


def axes_from_direct(therapy) -> Axes:
    """A sequence-directed agent is the specificity benchmark to beat."""
    strength = 0.85 if therapy.raises_protein else 0.55
    efficiency = 0.75
    if therapy.withdrawn:
        efficiency *= 0.3
    if therapy.modality in ("gene therapy",):
        efficiency *= 0.6      # one-time, but hardly a low-burden intervention
    return Axes(round(strength, 3), round(therapy.specificity, 3),
                round(efficiency, 3), evidence="measured")


# ------------------------------------------------------------- novelty test
# CTD's `InteractionActions` lists every action in a sentence, not the ones
# that concern this gene. "X results in increased activity of CYP3A4 protein
# which results in increased expression of CYP2B6 mRNA" carries both
# `increases^activity` and `increases^expression`, and reading the list alone
# files X as raising CYP3A4. "Y inhibits the reaction [X results in increased
# expression of CYP3A4 mRNA]" also lists `increases^expression` for Y, who
# blocks the increase. On CYP3A4, 27% of the rows with an expression action
# named another gene's expression, and the reaction forms are about one row in
# eight across the genes checked. So the sentence is read, not the list.
_EXPRESSION = re.compile(
    r"(?P<verb>results in (?P<dir>increased|decreased)|affects the|affects)"
    r"\s+expression of"
    r"(?:\s+and\s+(?:results in|affects)(?:\s+(?:increased|decreased))?"
    r"(?:\s+the)?[^\[\]]*?\bof)*"
    r"\s+(?P<obj>[^\s\[\]]+)(?:\s+(?P<form>mRNA|protein))?",
    re.IGNORECASE)
_ACTION_OF = {"up": "increases^expression", "down": "decreases^expression",
              "affects": "affects^expression"}
# Forms of the gene's expression that count as transcription. A sentence that
# does not say which form is read as the transcript, as CTD's own summary does.
TRANSCRIPT_FORMS = ("mRNA", None)
# The protein form: "results in increased expression of G protein". Protein
# level can move while the transcript does not (stability, translation,
# secretion), so these are kept apart from the transcript reports and are never
# read as evidence that the promoter is involved.
PROTEIN_FORMS = ("protein",)


def expression_effects(chemical: str, interaction: str, gene: str) -> list:
    """`[(direction, form, role)]` for what one CTD sentence says the
    chemical does to this gene's expression.

    `direction` is "up", "down" or "affects". `form` is "mRNA", "protein" or
    None when the sentence does not say. `role` is one of

      direct     the chemical is the subject of the expression clause
      modulated  the chemical changes another's effect, "inhibits the
                 reaction [X results in increased expression of G]". It is
                 read as "affects": blocking an inducer says nothing about
                 the gene's level without it, and the chemical is still a
                 known link to the gene
      indirect   the effect comes through another gene's change, "[A results
                 in ... of B] which results in ... expression of G"
      unparsed   the sentence names the gene but not in a shape read here

    A sentence that never mentions this gene's expression gives an empty
    list, which is how another gene's expression is kept out.
    """
    t = str(interaction or "")
    chem = str(chemical or "").strip().lower()
    sym = str(gene or "").strip().upper()
    if not t or not chem or not sym:
        return []
    low = t.lower()
    cpos = low.find(chem)
    out = []
    for m in _EXPRESSION.finditer(t):
        if m.group("obj").rstrip(".,;)").upper() != sym:
            continue
        word = (m.group("dir") or "").lower()
        direction = {"increased": "up", "decreased": "down"}.get(word,
                                                               "affects")
        form = m.group("form")
        form = {"mrna": "mRNA", "protein": "protein"}.get((form or "").lower())
        start = m.start()
        before = t[:start]
        depth = before.count("[") - before.count("]")
        role = "unparsed"
        if cpos < 0:
            pass
        elif depth <= 0:
            if re.search(r"\]\s+which\s*$", before, re.IGNORECASE):
                role = "indirect"
            elif cpos < start:
                role = "direct"
        else:
            # The innermost bracket still open where the clause begins.
            opened, stack = -1, []
            for i, ch in enumerate(before):
                if ch == "[":
                    stack.append(i)
                elif ch == "]" and stack:
                    stack.pop()
            opened = stack[-1] if stack else -1
            if opened < cpos < start:
                role = "direct"
            elif 0 <= cpos < opened:
                lead = low[:opened].rstrip()
                if lead.endswith(("inhibits the reaction",
                                  "promotes the reaction",
                                  "affects the reaction")):
                    role = "modulated"
                    direction = "affects"
        out.append((direction, form, role))
    return out


def strict_expression_links(rows: list, forms=TRANSCRIPT_FORMS) -> dict:
    """`{chemical_lower: [actions]}` from CTD rows, reading each sentence.

    Only `direct` and `modulated` clauses about the row's own gene count, and
    only in the requested forms. A row with no sentence falls back to its
    action list, so a table without the text still reads as before.
    """
    links: dict = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        chem = str(row.get("ChemicalName") or "").strip().lower()
        actions = str(row.get("InteractionActions") or "")
        if not chem or "expression" not in actions.lower():
            continue
        text = row.get("Interaction")
        if not text:
            links.setdefault(chem, []).extend(
                a for a in actions.split("|") if "expression" in a.lower())
            continue
        gene = row.get("GeneSymbol") or row.get("Input") or ""
        for direction, form, role in expression_effects(chem, text, gene):
            if role in ("direct", "modulated") and form in forms:
                links.setdefault(chem, []).append(_ACTION_OF[direction])
    return links


def expression_papers(rows: list, forms=TRANSCRIPT_FORMS) -> dict:
    """`{chemical_lower: {PubMed ids}}` for the rows `strict_expression_links`
    counts, so a link can say how many papers stand behind it."""
    out: dict = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        chem = str(row.get("ChemicalName") or "").strip().lower()
        actions = str(row.get("InteractionActions") or "")
        if not chem or "expression" not in actions.lower():
            continue
        text = row.get("Interaction")
        if text:
            gene = row.get("GeneSymbol") or row.get("Input") or ""
            ok = any(role in ("direct", "modulated") and form in forms
                     for _, form, role in expression_effects(chem, text, gene))
        else:
            ok = True
        if ok:
            ids = {p for p in str(row.get("PubMedIDs") or "").split("|") if p}
            out.setdefault(chem, set()).update(ids)
    return out


def ctd_expression_papers(client: Client, gene: str) -> dict:
    """`{chemical_lower: number of papers}`, or {} when the lookup fails.

    Asks the same question as `ctd_expression_links`, so the answer is already in
    the request cache.
    """
    try:
        res = client.get(CTD_BATCH, {
            "inputType": "gene", "inputTerms": gene.upper(),
            "report": "cgixns", "format": "json", "actionTypes": "ANY"})
    except HttpError:
        return {}
    rows = res if isinstance(res, list) else res.get("data", []) or []
    return {c: len(p) for c, p in expression_papers(rows).items()}


def ctd_protein_effects(client: Client, gene: str) -> dict | None:
    """What CTD records for this gene's protein level, read from the sentences.

    `{"links": {chemical_lower: [actions]}, "papers": {chemical_lower: n}}`, or
    `None` when the lookup failed (an empty `links` means CTD was searched and
    holds nothing). Only a sentence that names the protein form counts; the
    transcript reports are `ctd_expression_links`. It asks the same question as
    that function, so the answer is already in the request cache.
    """
    try:
        res = client.get(CTD_BATCH, {
            "inputType": "gene", "inputTerms": gene.upper(),
            "report": "cgixns", "format": "json", "actionTypes": "ANY"})
    except HttpError:
        return None
    rows = res if isinstance(res, list) else res.get("data", []) or []
    return {"links": strict_expression_links(rows, forms=PROTEIN_FORMS),
            "papers": {c: len(p) for c, p in
                       expression_papers(rows, forms=PROTEIN_FORMS).items()}}


def ctd_expression_links(client: Client, gene: str, *, strict: bool = True
                         ) -> dict[str, list[str]] | None:
    """{chemical_name_lower: [actions]} for curated expression effects.

    `None` when the lookup could not be made, `{}` when it was made and
    CTD knows of nothing. Those are different facts and the difference
    decides whether this package may call a candidate novel: an empty
    mapping says the literature was searched, and `None` says it was
    not.

    With `strict` (the default) each sentence is read for what it says about
    this gene's transcript; see `expression_effects`. `strict=False` is the
    earlier reading of the action list, kept to compare the two.
    """
    try:
        res = client.get(CTD_BATCH, {
            "inputType": "gene", "inputTerms": gene.upper(),
            "report": "cgixns", "format": "json", "actionTypes": "ANY",
        })
    except HttpError as exc:
        log.warning("CTD lookup failed for %s: %s - novelty cannot be "
                    "assessed for this run", gene, exc)
        return None

    rows = res if isinstance(res, list) else res.get("data", []) or []
    if strict:
        return strict_expression_links(rows)
    links: dict[str, list[str]] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        chem = str(row.get("ChemicalName") or "").strip().lower()
        actions = str(row.get("InteractionActions") or "")
        if chem and "expression" in actions.lower():
            links.setdefault(chem, []).extend(
                a for a in actions.split("|") if "expression" in a.lower())
    return links


def link_direction(actions) -> str:
    """Classify a CTD interaction as up / down / both / unclear.

    CTD encodes these as 'increases^expression' and 'decreases^expression'.
    The direction matters enormously: a chemical reported to RAISE your gene is
    a lead when you want it up and a warning when you want it down.
    """
    text = " ".join(actions).lower() if not isinstance(actions, str) else actions.lower()
    up = "increases^expression" in text or "increases expression" in text
    down = "decreases^expression" in text or "decreases expression" in text
    if up and down:
        return "both"          # usually different cell types or doses
    if up:
        return "up"
    if down:
        return "down"
    return "unclear"


# ----------------------------------------------------- what KIND of entry
# CTD indexes MeSH, which is a vocabulary for the toxicology and
# pharmacology literature rather than a drug list. A gene's chemical
# interactions therefore arrive as a mixture of three things that a single
# "predicted by this run?" column flattens into a wall of "no":
#
#   a drug        something somebody could be given
#   an exposure   formaldehyde, asbestos, bisphenol S, tobacco smoke
#   a category    "Anti-HIV Agents", "Plant Extracts", "Dietary Fats"
#
# Only the first can ever match a candidate, so reporting all three against
# one column makes a correct answer look like a broken join. On a real CCR5
# run, 33 of the 39 chemicals reported to lower the gene were exposures,
# research tools or categories.

# MeSH category and mixture headings. These are not compounds at all, so
# they can never match and are never leads.
_CATEGORY_WORDS = (
    "agents", "extracts", "compounds", "derivatives", "preparations",
    "substances", "pollution", "pollutants", "smoke", "dust", "fats",
    "oils", "tannins", "alkaloids", "flavonoids", "polyphenols",
    "mixtures", "particulate matter", "drugs, chinese herbal",
    "biological factors", "air pollutants", "hazardous substances",
)

# Hallmarks of an environmental or occupational exposure rather than a
# therapeutic. Deliberately conservative: a false "exposure" label would
# hide a real lead, so anything uncertain stays a candidate chemical.
_EXPOSURE_WORDS = (
    "asbestos", "benzo(a)pyrene", "benzopyrene", "formaldehyde", "arsenic",
    "cadmium", "chloroprene", "butadiene", "bisphenol", "phthalate",
    "dioxin", "dibenzofuran", "biphenyl", "atrazine", "diuron",
    "endosulfan", "pentachlorophenol", "methanesulfonate", "azoxymethane",
    "tetrabromodiphenyl", "dinitrotoluene", "naphthylamine",
    "sulfur mustard", "chloroethyl", "aflatoxin", "benzene", "crocidolite",
    "perchlorate", "dodecyl sulfate", "tobacco",
)


def entry_kind(name: str) -> str:
    """"category" | "exposure" | "chemical" for one CTD chemical name.

    "chemical" is the permissive default: a name this cannot place stays
    eligible to match a candidate, because a wrong "exposure" label would
    quietly hide a real lead.
    """
    n = (name or "").strip().lower()
    if not n:
        return "chemical"
    if any(w in n for w in _CATEGORY_WORDS):
        return "category"
    if any(w in n for w in _EXPOSURE_WORDS):
        return "exposure"
    return "chemical"


def summarise_links(links: dict, goal: str, candidate_names: set) -> dict:
    """Turn the raw CTD table into something readable.

    Two groups matter and they are different. Chemicals that also appear in our
    candidate list are confirmations of a route we predicted. Chemicals that do
    NOT appear are compounds with a reported effect on your gene that this
    pipeline missed entirely — often the more interesting group, because they
    are empirical and we have no mechanism for them.
    """
    from .congeners import salt_free_name

    wanted, opposed, unclear, ctd_only = [], [], [], []
    known = {n.lower() for n in candidate_names}
    # Salt-free too, the way every other name comparison in this package
    # works. CTD lists the parent and a drug table often lists a salt.
    known_stems = {salt_free_name(n) for n in known}

    for chem, actions in sorted(links.items()):
        d = link_direction(actions)
        stem = salt_free_name(chem)
        kind = entry_kind(chem)
        matched = (
            chem in known or stem in known_stems
            or any((len(chem) > 5 and chem in k) or (len(k) > 5 and k in chem)
                   for k in known))
        entry = {
            "chemical": chem.title(),
            "direction": d,
            "kind": kind,
            "actions": "; ".join(sorted(set(actions)))[:110],
            # A category heading is not a compound, so "did this run predict
            # it" has no answer rather than the answer "no".
            "in_candidates": False if kind == "category" else matched,
        }
        if d == goal:
            wanted.append(entry)
        elif d in ("up", "down"):
            opposed.append(entry)
        else:
            unclear.append(entry)
        # Only an actual compound counts as something this run missed. An
        # exposure or a MeSH heading was never a candidate, and listing
        # them as misses implied the scan had overlooked formaldehyde.
        if (not entry["in_candidates"] and d in (goal, "both")
                and kind == "chemical"):
            ctd_only.append(entry)

    return {
        "wanted_direction": wanted,
        "opposite_direction": opposed,
        "unclear_direction": unclear,
        "reported_but_not_predicted": ctd_only,
        "counts": {
            "total": len(links),
            "wanted": len(wanted),
            "opposite": len(opposed),
            "unclear": len(unclear),
            "not_predicted_by_us": len(ctd_only),
            # Across every direction, so a bucket with no match does not
            # read as "the matching is broken" when it matched elsewhere.
            "matched_any_direction": sum(
                1 for e in wanted + opposed + unclear if e["in_candidates"]),
            "n_category": sum(1 for e in wanted + opposed + unclear
                              if e["kind"] == "category"),
            "n_exposure": sum(1 for e in wanted + opposed + unclear
                              if e["kind"] == "exposure"),
            "n_chemical": sum(1 for e in wanted + opposed + unclear
                              if e["kind"] == "chemical"),
        },
    }


def check_novelty(candidate_name: str, links: dict | None
                  ) -> tuple[bool, str | None]:
    """(previously_linked, description). Substring match, deliberately generous.

    Over-calling a link as known is the safe error: it demotes a candidate you
    might have chased. Under-calling it means presenting a known relationship
    as a discovery, which is the one genuinely embarrassing failure mode.

    `links is None` means the lookup failed, and the description then
    says so rather than leaving a bare False that reads as "nobody has
    tried this". Use `novelty_assessed` to tell the two apart.
    """
    name = candidate_name.strip().lower()
    if links is None:
        return False, "not checked: the CTD lookup failed for this run"
    if not name or not links:
        return False, None
    if name in links:
        return True, "; ".join(sorted(set(links[name])))[:120]
    for chem, actions in links.items():
        if (len(name) > 5 and name in chem) or (len(chem) > 5 and chem in name):
            return True, f"matches CTD entry '{chem}': " + \
                         "; ".join(sorted(set(actions)))[:100]
    return False, None


# ------------------------------------------------------------- comparison
def dominates(a: Axes, b: Axes, margin: float = 0.0,
              tolerance: float = TOLERANCE) -> bool:
    """a beats b on >=1 axis by `margin`, and is not worse anywhere by > tolerance."""
    better = any(x > y + margin for x, y in zip(a.as_tuple(), b.as_tuple()))
    not_worse = all(x >= y - tolerance for x, y in zip(a.as_tuple(), b.as_tuple()))
    return better and not_worse


def required_margin(candidate: Axes, baseline: Axes) -> float:
    """Predicted challengers must clear a measured incumbent by a real margin."""
    if candidate.evidence == "predicted" and baseline.evidence == "measured":
        return PREDICTION_MARGIN
    return 0.0


AXIS_NAMES = ("stronger activation", "higher specificity", "superior efficiency")


def baseline_envelope(established: list) -> Axes | None:
    """Per-axis best across ALL established options, not one cherry-picked drug.

    Taking a single incumbent and summing its axes contradicts the Pareto
    design: a challenger could "beat the incumbent" on specificity while a
    different established drug is more specific still. Beating the standard of
    care means beating the best available value on that axis, whichever drug
    holds it. Anything weaker flatters challengers.
    """
    if not established:
        return None
    return Axes(
        strength=max(c.axes.strength for c in established),
        specificity=max(c.axes.specificity for c in established),
        efficiency=max(c.axes.efficiency for c in established),
        evidence=("measured" if any(c.axes.evidence == "measured"
                                    for c in established) else "predicted"),
    )


def compare(candidate: Candidate, baseline: Candidate | None,
            envelope: Axes | None = None) -> Candidate:
    """Annotate a candidate with how it stacks up against the standard of care.

    With no incumbent, every axis the candidate scores decently on is recorded
    as a win. Leaving `wins_on` empty there silently hid every candidate for
    genes with no existing treatment - which is the case you most want to see.
    """
    if baseline is None and envelope is None:
        candidate.verdict = "no incumbent for this gene; nothing to beat"
        candidate.notes.append(
            "No existing intervention was found that controls this gene's "
            "expression, so any credible candidate is worth pursuing - but "
            "check that the absence is real rather than a gap in the databases.")
        # Record axes the candidate is actually decent on, so it is not
        # filtered out of the novel-improvement list for lack of a comparator.
        candidate.wins_on = [n for n, v in zip(AXIS_NAMES, candidate.axes.as_tuple())
                             if v >= 0.5]
        if not candidate.wins_on:
            candidate.notes.append("scores below 0.5 on every axis; weak even "
                                   "with no competition")
        return candidate

    ref = envelope if envelope is not None else baseline.axes
    m = required_margin(candidate.axes, ref)
    for name, c, b in zip(AXIS_NAMES, candidate.axes.as_tuple(),
                          ref.as_tuple()):
        if c > b + m:
            candidate.wins_on.append(name)
        elif c < b - TOLERANCE:
            candidate.loses_on.append(name)

    if m > 0:
        candidate.notes.append(
            f"predicted scores compared against measured ones; a margin of "
            f"{m:.2f} was required before counting a win")
    if envelope is not None:
        candidate.notes.append(
            "compared against the best established value on each axis, not "
            "against a single incumbent drug")

    if candidate.wins_on and not candidate.loses_on:
        candidate.verdict = "beats the incumbent on " + ", ".join(candidate.wins_on)
    elif candidate.wins_on:
        candidate.verdict = ("trade-off: better on " + ", ".join(candidate.wins_on)
                             + "; worse on " + ", ".join(candidate.loses_on))
    elif candidate.loses_on:
        candidate.verdict = "worse than the incumbent; no reason to switch"
    else:
        candidate.verdict = "within noise of the incumbent"
    return candidate


def pareto_front(candidates: list[Candidate]) -> list[Candidate]:
    """Non-dominated set. Nothing here is beaten on every axis by another."""
    front = []
    for c in candidates:
        if not any(dominates(o.axes, c.axes) for o in candidates if o is not c):
            front.append(c)
    return front


# ------------------------------------------------------------------ driver
def discover(client: Client, gene: str, grades, tf_results, direct_hits,
             outcomes=None, *, goal: str = "up", check_ctd: bool = True,
             demo: bool = False) -> dict:
    """Assemble candidates, establish the baseline, and report what beats it."""
    links = {} if (demo or not check_ctd) else ctd_expression_links(client, gene)
    by_symbol = {r.symbol: r for r in tf_results}
    outcome_by = {(o["drug"], o["tf_symbol"]): o for o in (outcomes or [])}

    class _O:                                   # tiny adapter for p_point
        def __init__(self, d): self.p_point = d["p_point"]

    candidates: list[Candidate] = []

    for t in direct_hits:
        c = Candidate(name=t.drug, route="direct", tf_symbol=None,
                      axes=axes_from_direct(t))
        c.previously_linked = True
        c.link_source = f"established {t.modality} for {t.gene}"
        candidates.append(c)

    for g in grades:
        o = outcome_by.get((g.drug_name, g.tf_symbol))
        c = Candidate(
            name=g.drug_name, route="via TF", tf_symbol=g.tf_symbol,
            axes=axes_from_grade(g, by_symbol.get(g.tf_symbol),
                                 _O(o) if o else None))
        c.previously_linked, c.link_source = check_novelty(g.drug_name,
                                                           links)
        c.novelty_assessed = links is not None
        candidates.append(c)

    # Standard of care: the per-axis best across every established option.
    established = [c for c in candidates if c.previously_linked]
    envelope = baseline_envelope(established)
    baseline = (max(established, key=lambda c: sum(c.axes.as_tuple()))
                if established else None)

    for c in candidates:
        if c is not baseline:
            compare(c, baseline, envelope)
    if baseline:
        baseline.verdict = ("current best established option; the bar is the "
                            "per-axis best across all established options")

    front = pareto_front(candidates)
    novel_front = [c for c in front if not c.previously_linked and c.wins_on]

    # Nothing is silently dropped: account for every candidate.
    novel_all = [c for c in candidates if not c.previously_linked]
    novel_dominated = [c for c in novel_all
                       if c.wins_on and c not in front]
    novel_no_win = [c for c in novel_all if not c.wins_on]

    link_summary = summarise_links(links, goal, {c.name for c in candidates})
    papers = (ctd_expression_papers(client, gene)
              if (links and not demo) else {})
    protein = (None if (demo or not check_ctd)
               else ctd_protein_effects(client, gene))

    return {
        "gene": gene,
        "goal": goal,
        "reported_links": link_summary,
        "baseline": baseline.to_dict() if baseline else None,
        "baseline_envelope": envelope.to_dict() if envelope else None,
        "has_incumbent": envelope is not None,
        "n_candidates": len(candidates),
        "n_previously_linked": len(established),
        "n_novel": len(novel_all),
        "pareto_front": [c.to_dict() for c in front],
        "novel_improvements": [c.to_dict() for c in novel_front],
        "novel_dominated": [c.to_dict() for c in novel_dominated],
        "novel_no_advantage": [c.to_dict() for c in novel_no_win],
        "accounting": {
            "total": len(candidates),
            "previously_linked": len(established),
            "novel_and_on_front": len(novel_front),
            "novel_but_dominated": len(novel_dominated),
            "novel_no_advantage": len(novel_no_win),
        },
        "ctd_links_found": len(links or {}),
        # Papers behind each chemical's record, for the literature route.
        "ctd_papers": papers,
        # What CTD records for the protein level, kept apart from the
        # transcript. None when not looked up or the lookup failed.
        "reported_protein_links": (
            summarise_links(protein["links"], goal,
                            {c.name for c in candidates})
            if protein else None),
        "ctd_protein_papers": (protein or {}).get("papers") or {},
        # Whether the literature was searched at all. Without this a
        # failed lookup and a genuinely empty answer both arrive as
        # zero, and every candidate then reads as novel.
        "ctd_assessed": links is not None,
        "all_candidates": [c.to_dict() for c in candidates],
    }


def ctd_ledger_detail(n_links: int, n_wanted: int, assessed: bool) -> str:
    """What the evidence ledger should say about the CTD lookup.

    Three outcomes, and two of them used to share a sentence. "No
    curated effects returned" covered both a search that found nothing
    and a request that never completed, and only the first of those
    lets anything downstream be called novel.

    A function rather than three branches inside `run`, so the wording
    can be tested without a network round trip.
    """
    if n_links:
        return (f"{n_links} chemical(s) with curated expression effects "
                f"on this gene; {n_wanted} move it the way you want")
    if assessed:
        return ("CTD was searched and holds no curated expression effect "
                "for this gene")
    return ("the CTD lookup failed, so nothing here is known to be novel "
            "- absence of a retrieved link is not absence of a link")
