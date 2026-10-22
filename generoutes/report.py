"""Turning a finished run into something a person reads.

Split out of `pipeline.py`, which was four and a half thousand lines and
named for the orchestration rather than for this. Nothing here fetches,
scores or decides anything: every function takes the result dictionary a
run already produced and returns text.

The dependency runs one way. This module imports nothing from
`pipeline`, so the formatting can be read, tested and changed without
the orchestration in view, and `pipeline` re-exports `markdown_report`
and `write_outputs` so the public interface is unchanged.

The ordering of a report is deliberate and is the one thing worth
knowing before editing it. The evidence ledger comes first, because a
run whose lookups quietly failed is otherwise indistinguishable from one
that worked and found nothing. Routes are ordered by evidence tier and
only then by score, because scores are comparable within a tier and not
across one.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from . import (analogs as analog_mod, catalogue as cat_mod,
               delivery as dlv_mod, drugs as drugs_mod, fda as fda_mod,
               liabilities as liab_mod, ligands as lig_mod,
               knockout as ko_mod, outcome as outcome_mod,
               predictions as pred_mod,
               systemic as sys_mod,
               tfactivity as tfa_mod)

log = logging.getLogger(__name__)


# How many ranked factors the LINCS drug scan covers. Each is a pass over
# roughly 33,000 signatures, about a second, so this is a wall-clock cap
# rather than a claim that the rest are uninteresting.
# A blank table cell. Named rather than inlined because a
# backslash escape inside an f-string expression is PEP 701
# syntax: legal from 3.12 and a SyntaxError on 3.10 and 3.11,
# which pyproject declares support for. The module would not
# import at all on those versions.
EMDASH = "\u2014"


# ------------------------------------------------------------------ reports
def _fmt_offset(o) -> str:
    return f"{o:+d}" if isinstance(o, int) else "?"


def _dep_label(tf: dict) -> str:
    """The DepMap verdict, shortened for a table cell."""
    v = tf.get("dependency")
    if not v:
        return "not scored"
    return {
        "essential in nearly every line": "yes, nearly all",
        "essential in most lines": "yes, most",
        "essential in a substantial minority": "some lines",
        "essential in a minority of lines": "some lines",
        "essential in a handful of lines": "a few lines",
        "dispensable across the panel": "no",
    }.get(v, v)


def _protein_screen_block(pp: dict) -> list:
    """The patient-cell protein screen, when that is the endpoint.

    Rendered instead of the reporter construct. The ordering is the order
    somebody would set the experiment up in, and the normalisation
    section is deliberately before the library: for a cell-cycle gene it
    is the part that decides whether the screen measures anything.
    """
    L = [f"### Primary screen: {pp['endpoint']}", ""]
    L.append(f"The endpoint is the protein and not a reporter, because "
             f"{pp['why']}. A promoter reporter measures transcription "
             f"from a cloned fragment, so it cannot see a stabiliser, a "
             f"splice modifier or a readthrough agent - and those are the "
             f"mechanisms in play here. This screen is also "
             f"mechanism-agnostic, which means it does not depend on any "
             f"ranking in this report being right.")
    L.append("")
    for title, key in (("Cells", "cells"), ("Readout", "readout"),
                       ("Normalisation, and the control that matters most",
                        "normalisation"), ("Controls", "controls"),
                       ("Thresholds", "thresholds"),
                       ("Triage, in order", "triage")):
        L.append(f"**{title}**")
        L.append("")
        for item in pp.get(key) or []:
            L.append(f"- {item}")
        L.append("")
    lib = pp.get("library") or {}
    L.append("**Library**")
    L.append("")
    L.append(f"- {lib.get('what', '')}")
    for src in lib.get("sources") or []:
        L.append(f"  - {src}")
    L.append(f"- Format: {lib.get('format', '')}")
    L.append("")
    L.append(f"*{pp.get('caveat', '')}*")
    L.append("")
    return L


def _undrugged_regulators(res: dict) -> list:
    """Factors RECORDED as regulating the gene that no approved drug reaches.

    A different claim from the motif list above it. A motif hit says a
    sequence resembles a binding site; these factors are in a curated table
    or moved the gene in a knockout experiment. They produce no route,
    because there is no molecule to give - which is exactly why they used
    to be dropped, and exactly why they have to be named. They are the list
    of what a chemistry effort would aim at.
    """
    rows = res.get("regulators_without_drugs") or []
    if not rows:
        return []
    # Direct first, and all of them. With --regulator-steps 2 the two-step
    # factors outnumber the direct ones twenty to one, so a list sorted by
    # symbol showed ABL1 to CREB1 and hid FOXM1, MYC, TP53 and ZNF143 -
    # the factors this section exists to name.
    direct = [r for r in rows if not r.get("via")]
    indirect = [r for r in rows if r.get("via")]

    L = ["#### Recorded regulators with no approved drug", ""]
    L.append("These are not motif matches. Each is named by a curated table "
             "or by a knockout experiment as regulating this gene, and none "
             "has an approved drug that moves it the way this run needs. "
             "They generate no route and are listed so the regulator "
             "section is not read as the full set of regulators. This is "
             "the list a chemistry effort would aim at.")
    L.append("")
    if direct:
        L.append("| Factor | Effect on the gene | Would need to be | "
                 "Recorded by |")
        L.append("|---|---|---|---|")
        for r in sorted(direct, key=lambda r: r.get("symbol") or ""):
            L.append(
                f"| {r.get('symbol')} | "
                f"{'raises it' if r.get('sign') == 1 else 'lowers it'} | "
                f"{'activated' if r.get('need') == 'up' else 'inhibited'} | "
                f"{', '.join(r.get('sources') or []) or '?'} |")
        L.append("")
    else:
        L.append("No factor is recorded as regulating this gene directly "
                 "without already having a drug.")
        L.append("")
    if indirect:
        L.append(f"A further **{len(indirect)}** factor(s) reach the gene "
                 f"only through one of the regulators above, so each is a "
                 f"longer way round with a sign that is the product of two "
                 f"edges. They are in `results.json` rather than here, "
                 f"because at this distance the list is a property of the "
                 f"network and not of this gene.")
        L.append("")
    return L


def _corpus_note(res: dict) -> list:
    """Which corpora ran, which did not, and what agreed across them."""
    from . import corpora as cp

    scans = res.get("corpus_scans") or {}
    agree = res.get("corpus_agreement") or {}
    if not scans:
        return []

    present, absent = set(), set()
    for per_tf in scans.values():
        for key, r in per_tf.items():
            (present if r.get("available") else absent).add(key)

    L = []
    if len(present) > 1:
        L.append(f"**Scanned across {len(present)} independent corpora** "
                 f"({', '.join(sorted(present))}). They differ in species, "
                 f"tissue, platform and dose, so they share few systematic "
                 f"errors - which is the point of running more than one. "
                 f"Every signature in a single corpus shares that corpus's "
                 f"biases, so no amount of internal replication detects "
                 f"them.")
    else:
        only = ", ".join(sorted(present)) or "none"
        L.append(f"**Only one corpus was available ({only}).** Every "
                 f"signature in it shares its platform, species and dose "
                 f"regime, so internal replication cannot detect a "
                 f"systematic error common to all of them. A second corpus "
                 f"is what would.")
    if absent:
        L.append("")
        L.append(f"Not on disk, and therefore skipped rather than searched: "
                 f"{', '.join(sorted(absent))}. Absent is not the same as "
                 f"finding nothing. Run "
                 f"`python -m generoutes.corpora` for the downloads.")
    L.append("")

    flat = [a for per_tf in agree.values() for a in per_tf]
    if flat:
        consistent = [a for a in flat if a["consistent"]]
        conflicting = [a for a in flat if not a["consistent"]]
        if consistent:
            L.append("#### Agreeing across corpora")
            L.append("")
            L.append("The strongest empirical result this package can "
                     "produce. These compounds move the factor the same "
                     "way in systems that share almost no methodology.")
            L.append("")
            L.append("| Drug | Corpora | Direction |")
            L.append("|---|---|---|")
            for a in consistent[:15]:
                L.append(f"| {a['drug']} | "
                         f"{', '.join(sorted(a['corpora']))} | "
                         f"{a['direction']} |")
            L.append("")
        if conflicting:
            L.append(f"**{len(conflicting)} compound(s) appear in more than "
                     f"one corpus with opposite directions** "
                     f"({', '.join(a['drug'] for a in conflicting[:6])}). "
                     f"That usually means the systems are not comparable - "
                     f"rat liver in vivo against a human cell line, at "
                     f"different doses - rather than that either reading is "
                     f"wrong. It is shown because hiding it would make the "
                     f"agreeing rows look more uniform than the data is.")
            L.append("")
    elif len(present) > 1:
        L.append("No compound moved a regulon consistently in more than one "
                 "corpus. With this little overlap between their compound "
                 "sets that is the expected outcome, not a negative result.")
        L.append("")

    for key in sorted(present):
        c = cp.BY_KEY.get(key)
        if c and c.species != "human":
            L.append(f"*{c.name} is {c.species} ({c.system}). Gene symbols "
                     f"cross a species boundary and the regulon is not "
                     f"guaranteed to be conserved, so treat a hit found "
                     f"only there as a lead to check in a human system.*")
            L.append("")
    return L


def _docking_section(runs: list, consensus: dict) -> list:
    """Docking validation, and the compounds two independent methods both
    picked.

    The routes from a dock already appear in the tables. What does not
    appear anywhere is whether the pocket was validated at all, which is
    the only thing that separates a docking score from a random number.
    A screen that cannot separate known actives from decoys has not been
    shown to rank anything, and its hits are reported as unvalidated
    rather than withheld.
    """
    if not runs:
        return []
    from .docking import MIN_DOCKING_AUC, MIN_EF1

    L = ["### Structure-based docking", ""]
    L.append("Docking asks whether a molecule fits the pocket. It says "
             "nothing about direction: a compound that binds a factor may "
             "activate it, block it, or do neither to transcription of this "
             "gene.")
    L.append("")
    L.append("| Pocket | Docked | Failed | ROC AUC | EF 1% | Validated? |")
    L.append("|---|---|---|---|---|---|")
    for r in runs:
        auc = r.get("auc")
        ef1 = r.get("ef1")
        L.append(f"| {r.get('target') or '?'} | {r.get('n_docked', 0)} | "
                 f"{r.get('n_failed', 0)} | "
                 f"{f'{auc:.2f}' if auc is not None else 'not tested'} | "
                 f"{f'{ef1:.1f}x' if ef1 is not None else 'not tested'} | "
                 f"{'yes' if r.get('usable') else '**no**'} |")
    L.append("")
    L.append(f"Validated means the same scoring function separated known "
             f"actives from decoys in this pocket, at AUC "
             f"{MIN_DOCKING_AUC:.2f} or better and {MIN_EF1:.0f}x or better "
             f"enrichment in the top 1%. Below that the geometry may still "
             f"be right and the score is not evidence.")
    L.append("")

    failed = [r for r in runs if not r.get("usable")]
    if failed:
        n_f = len(failed)
        L.append(f"**{n_f} {'pocket' if n_f == 1 else 'pockets'} did not "
                 f"validate** "
                 f"({', '.join(str(r.get('target')) for r in failed[:5])}). "
                 f"Their hits are listed for the record and should not be "
                 f"ranked against hits from a pocket that did.")
        L.append("")
    for r in runs:
        for c in (r.get("caveats") or [])[:3]:
            L.append(f"- {r.get('target')}: {c}")
    if any(r.get("caveats") for r in runs):
        L.append("")

    if consensus:
        n = sum(len(v) for v in consensus.values())
        L.append(f"#### {n} {'compound' if n == 1 else 'compounds'} "
                 f"picked by both methods")
        L.append("")
        L.append("Docking reads pocket shape. The ligand screen reads "
                 "chemistry measured at this receptor. They share no input, "
                 "so agreement is worth more than a better score from "
                 "either - and a compound either method missed is not thereby "
                 "ruled out.")
        L.append("")
        L.append("| Factor | Drug | Dock score | Ligand screen quadrant |")
        L.append("|---|---|---|---|")
        for sym, hits in sorted(consensus.items()):
            for h in hits[:8]:
                sc = h.get("docking_score")
                L.append(f"| {sym} | **{h.get('drug')}** | "
                         f"{f'{sc:.1f}' if isinstance(sc, (int, float)) else EMDASH}"
                         f" | {h.get('ligand_quadrant') or EMDASH} |")
        L.append("")
        L.append("Still a binding hypothesis. Neither method addresses "
                 "direction, and neither was run in the tissue this gene "
                 "matters in.")
        L.append("")
    elif any(r.get("usable") for r in runs):
        L.append("No compound reached the top slice of both the dock and the "
                 "ligand screen. That is a real disagreement between two "
                 "independent arguments, not a missing result - treat the "
                 "single-method hits above accordingly.")
        L.append("")
    return L


def _sensitivity_section(sens: dict) -> list:
    """Drugs that behave as though they remove a factor.

    The only section here whose candidates needed no annotation of any
    kind, so it is the only one that can name a compound nobody has
    connected to the gene.
    """
    if not sens:
        return []
    from .sensitivity import MIN_ABS_R, MAX_Q, MIN_LINES

    assessed = {k: v for k, v in sens.items() if v.get("assessed")}
    if not assessed:
        reasons = {v.get("error") for v in sens.values() if v.get("error")}
        L = ["### Drugs that behave like losing the factor", ""]
        L.append("Nothing was scanned. " + (
            "; ".join(r for r in reasons if r)[:400] or "No reason recorded."))
        L.append("")
        return L

    L = ["### Drugs that behave like losing the factor", ""]
    L.append("Every other route in this report starts from something "
             "somebody wrote down - a curated ligand, a motif, a "
             "structure, a signed edge. None of them can surface a drug "
             "nobody has connected to the factor. This one reads a "
             "measured phenotype instead: if a compound hurts the same "
             "cell lines that knocking the factor out hurts, it behaves "
             "as though it removes the factor.")
    L.append("")

    rows = [(sym, h) for sym, v in sorted(assessed.items())
            for h in (v.get("hits") or [])]
    if not rows:
        n = sum(v.get("n_drugs_tested", 0) for v in assessed.values())
        L.append(f"{n} compound(s) were correlated across "
                 f"{len(assessed)} factor(s) and none cleared "
                 f"|r| >= {MIN_ABS_R} at q <= {MAX_Q}. That is the common "
                 f"outcome. It bounds what this method contributes here "
                 f"rather than saying no such drug exists: the panel is "
                 f"proliferating cancer lines, and a factor that matters "
                 f"elsewhere leaves no trace in it.")
        L.append("")
        return L

    L.append("| Factor | Drug | Partial r | Raw r | Lines | q | Behaves as |")
    L.append("|---|---|---|---|---|---|---|")
    for sym, h in sorted(rows, key=lambda t: -abs(t[1].get("r") or 0))[:20]:
        L.append(f"| {sym} | **{h.get('drug')}** | {h.get('r')} | "
                 f"{h.get('r_raw')} | {h.get('n_lines')} | "
                 f"{h.get('q_value')} | {h.get('direction')} |")
    L.append("")
    L.append(f"**Partial r** controls for how easily each cell line dies "
             f"under the whole panel. Some lines die under almost "
             f"anything and those same lines suffer most knockouts, which "
             f"puts a correlation between almost any drug and almost any "
             f"gene. The control moves a number either way - it removes a "
             f"shared component that can inflate a correlation or mask "
             f"one - so where the two columns differ a lot the raw value "
             f"was mostly about the cell lines. Estimates use at least "
             f"{MIN_LINES} shared lines.")
    L.append("")
    L.append("**This is not a target.** A correlated drug may act on the "
             "factor, on anything else the same lines need, or on a "
             "pathway the factor sits in, and the first is the least "
             "likely of the three for any one hit. Nothing here records "
             "any of them binding anything: the finding is a reason to "
             "run a binding assay, not a result from one.")
    L.append("")
    L.append("**A proliferation readout.** Both panels score growth of "
             "transformed cells in a dish, so a factor whose role is in a "
             "differentiated or post-mitotic cell looks silent here for "
             "reasons that have nothing to do with its importance.")
    L.append("")
    return L


def _depmap_section(deps: dict, res: dict) -> list:
    """What a cell loses without each factor, and what that does not mean."""
    from . import depmap as dm

    if not deps:
        return []
    n_scanned = len(res.get("tfs") or [])
    L = ["### What a cell loses without each factor (DepMap)", ""]
    L.append("The regulon column counts how many genes a factor controls, "
             "from literature curation, which measures how much anyone has "
             "published as much as what the factor does. This is the same "
             "collateral question answered by measurement: CRISPR knockout "
             "of each factor across roughly 1,100 cancer cell lines.")
    L.append("")
    L.append("It also reaches factors the curated tables cannot. A factor "
             "nobody has written a regulon for still has a knockout "
             "phenotype here.")
    L.append("")
    L.append(dm.summarise({k: _as_dep(v) for k, v in deps.items()},
                          requested=n_scanned))
    L.append("")
    L.append("| Factor | Cells need it? | Mean effect | Lines dependent | "
             "Most dependent lineage |")
    L.append("|---|---|---|---|---|")
    rows = sorted(deps.values(), key=lambda d: d.get("mean_effect") or 0.0)
    for d in rows[:25]:
        # Only when there is a dependency to localise. For a dispensable
        # factor the "most dependent" lineage is whichever one drifted
        # furthest below zero, which is noise dressed as a finding.
        lin = (d.get("top_lineages") or [{}])[0]
        lin_cell = "—"
        if d.get("n_dependent") and lin.get("lineage"):
            lin_cell = f"{lin['lineage']} ({lin['mean_effect']:+.2f})"
        L.append(f"| {d['gene']} | {d['verdict']} | "
                 f"{d['mean_effect']:+.2f} | "
                 f"{d['n_dependent']}/{d['n_lines']} "
                 f"({d['fraction_dependent']:.0%}) | {lin_cell} |")
    L.append("")
    L.append("**Mean effect** is on DepMap's own scale: 0 means knockout "
             "changed nothing, and -1 is the median of known "
             "common-essential genes. A line counts as dependent below "
             "-0.5.")
    L.append("")
    L.append("#### What this does not mean")
    L.append("")
    L.append("- **These are cancer cell lines dividing in a dish.** "
             "Essentiality here is a statement about proliferating "
             "transformed cells. A factor whose job is in a post-mitotic "
             "or differentiated cell looks dispensable for reasons that "
             "have nothing to do with its importance in your tissue.")
    L.append("- **Knockout is not inhibition.** The protein is removed "
             "entirely, including scaffolding roles no ligand would touch. "
             "Essential genes are drugged successfully at partial "
             "occupancy all the time, so a dependency is a reason for "
             "care rather than a veto.")
    L.append("- **Nothing here is about your gene.** This bounds the cost "
             "of acting on a factor. Whether acting on it moves your gene "
             "is the rest of the report.")
    L.append("- **A dispensable result is not permission.** It bounds the "
             "proliferative cost and nothing else.")
    L.append("")
    sel = [d for d in deps.values() if d.get("selective")]
    if sel:
        L.append(f"**The selective ones are the interesting ones.** "
                 f"{', '.join(d['gene'] for d in sel[:6])} matter in a "
                 f"minority of lines rather than all of them, so the "
                 f"requirement has a context. The lineage column is where "
                 f"to look for what that context is, and whether your "
                 f"tissue shares it is the question this panel cannot "
                 f"answer.")
        L.append("")
    return L


def _as_dep(d: dict):
    """A dict row back into the shape depmap.summarise reads."""
    from . import depmap as dm
    return dm.Dependency(**{k: v for k, v in d.items()
                            if k in dm.Dependency.__dataclass_fields__})


def _drug_glossary(routes: list, res: dict) -> list:
    """Every distinct drug in this block, with what it is and what it treats.

    One entry per drug rather than a column on each table. The tables are
    already wide, the same drug appears in several of them, and a sentence
    repeated five times is what made these reports long. Salts collapse
    onto one entry.
    """
    from . import indications as ind_mod

    mechanisms = {}
    for r in routes:
        name = r.get("drug")
        if name and name not in mechanisms:
            mechanisms[name] = (r.get("mechanism") or "").strip() or None

    described = ind_mod.describe_all(list(mechanisms), mechanisms=mechanisms)
    if not described:
        return []

    known = {k: v for k, v in described.items() if v.get("treats")}
    unknown = [k for k, v in described.items() if not v.get("treats")]

    L = ["### What these drugs are", ""]
    L.append("Every route in this report is off-label by construction: the "
             "drug was approved for one thing and is being considered for "
             "another. The approved use is the baseline that departure is "
             "measured from, so it is listed here once per drug rather "
             "than repeated in each table above.")
    L.append("")

    if known:
        L.append("| Drug | What it is | Approved to treat |")
        L.append("|---|---|---|")
        for name in sorted(known):
            d = known[name]
            L.append(f"| {name} | {d['what_it_is'] or '—'} | "
                     f"{d['treats']} |")
        L.append("")

    if unknown:
        L.append(f"**No description established for {len(unknown)}**: "
                 + ", ".join(sorted(unknown)[:25])
                 + ("…" if len(unknown) > 25 else "")
                 + ". That means this package could not find one offline, "
                 "not that the compound has no approved use. Several are "
                 "environmental agents and research chemicals from the CTD "
                 "tables, which have no indication because they are not "
                 "drugs.")
        L.append("")

    L.append(ind_mod.summarise(described))
    L.append("")
    return L


def _screen_scope_note(ran: list, skipped: list) -> list:
    """Why the screen covers the receptors it does, and no others.

    Written because the section looks nearly identical across unrelated
    genes, which invites the reading that these drugs are being selected
    for the gene. They are not. They are the receptors with enough
    measured ligands to train a model on.
    """
    from .agonism import PCHEMBL_CUTOFF
    from .screening import MIN_ACTIVES, MIN_BINDING_AUC

    out = ["## Virtual screening", ""]
    out.append(f"A screen needs a model, and a model needs data. This one "
               f"trains on measured binding at the receptor, so it runs "
               f"only where ChEMBL holds at least {MIN_ACTIVES} ligands at "
               f"pChEMBL {PCHEMBL_CUTOFF:.0f} or better, and only when the "
               f"trained model separates known actives from decoys at CV "
               f"AUC {MIN_BINDING_AUC} or above.")
    out.append("")
    if ran:
        out.append(f"**Screened:** {', '.join(sorted(ran))}.")
    if skipped:
        out.append(f"**Not screened:** {', '.join(sorted(skipped))} "
                   f"{chr(8212)} each one's reason is stated under its own "
                   f"heading below. A skip is always a limit of the "
                   f"available data, either too few measured ligands or a "
                   f"model that could not separate them from decoys, and "
                   f"never a judgement about the receptor.")
    out.append("")
    out.append("**This is why the same receptor families appear here "
               "whatever gene you run.** The filter selects for how much "
               "medicinal chemistry a receptor has attracted, not for "
               "anything about your promoter. Nuclear receptors, and the "
               "retinoid and steroid receptors above all, have ligand "
               "series running to thousands of compounds; most "
               "transcription factors have none at all. So a receptor "
               "reaching this section means it is well studied, and the "
               "drugs listed under it are that receptor's pharmacology, "
               "not a result about your gene.")
    out.append("")
    out.append("Read these rows as *candidate ligands for a factor that has "
               "a site in your promoter*. Whether moving that factor moves "
               "your gene is a separate question, and the activity and "
               "evidence sections above speak to it. A long list here is "
               "evidence about the receptor's chemistry; it is not evidence "
               "that the gene will respond.")
    out.append("")
    return out


def _regulon_column_note(any_regulon: bool, any_understudied: bool,
                         tfs: list) -> list:
    """What the regulon column means, and the two ways it misleads."""
    if not any_regulon:
        return ["*No regulon counts are available. Pass `--regulon-file` "
                "with a TRRUST-style table to fill the column; without one "
                "there is no basis for saying how much else a factor "
                "controls.*", ""]

    out = ["**Other genes it controls** is the collateral question. A drug "
           "acting through a factor moves everything that factor regulates, "
           "so this column says how wide the effect will be. Read it against "
           "the score: a strong route through a factor that controls two "
           "thousand genes is a worse starting point than a weaker route "
           "through one that controls thirty."]

    exact = [t for t in tfs if t.get("n_regulon") is not None
             and "tier" not in (t.get("regulon_provenance") or "")]
    tiered = [t for t in tfs if "tier" in (t.get("regulon_provenance") or "")]
    uncounted = [t for t in tfs if t.get("n_regulon") is None]

    if tiered:
        out.append("")
        out.append(f"A count written `~N` is a bundled estimate for a factor "
                   f"the regulon file does not list ({len(tiered)} here), not "
                   f"a count of curated edges. Treat it as an order of "
                   f"magnitude.")
    if uncounted:
        out.append("")
        out.append(f"`not counted` means the regulon table has no entry and "
                   f"there is no bundled estimate either ({len(uncounted)} "
                   f"here). It is not a small number. Most of the roughly "
                   f"1,600 human transcription factors have no curated "
                   f"regulon at all.")
    if any_understudied:
        out.append("")
        out.append("A ⚠ marks a count below the credible floor. These "
                   "are the trap on this column: a factor with eight curated "
                   "targets is far more likely to be understudied than "
                   "selective, and reading a small number as specificity "
                   "selects for obscurity. The same correction applies here "
                   "that `drug_scoring` applies to protein selectivity.")
    out.append("")
    out.append("Counts come from literature curation, so they measure what "
               "has been published as much as what the factor does. SP1 has "
               "thousands of recorded targets partly because it regulates "
               "many genes and partly because thousands of people assayed "
               "it.")
    out.append("")
    return out


def _sig_q(value) -> float:
    """A q-value for comparison, with None - not 0.0 - meaning missing.

    `value or 1` reads a q of exactly 0.0, the MOST significant result
    possible, as the least. `tfactivity._q` exists for this and this line
    reimplemented the broken version. Latent rather than live today only
    because `stats.benjamini_hochberg` floors its output at Q_FLOOR, an invariant
    held in another module.
    """
    return 1.0 if value is None else float(value)


def _route_cell(row: dict) -> str:
    """FDA route, compressed. "local only" is the one that rules a drug out."""
    r = row.get("route")
    if r == "local only":
        return "**local only**"
    if r in ("systemic", "unclear"):
        return r
    return "—"


def _bbb_footnote(counts: dict) -> list:
    """What the BBB column means, stated where it is read.

    `counts` maps the cells actually printed ("yes", "yes?", "no", "no?",
    "—") to how many rows carried them, so the note describes this table
    rather than warning about states that did not occur in it.
    """
    if not counts:
        return []
    total = sum(counts.values())
    n_unknown = counts.get("—", 0)
    n_yes = counts.get("yes", 0)
    n_borderline = counts.get("yes?", 0) + counts.get("no?", 0)

    out = ["**BBB?** is CNS MPO (Wager et al.) against the same scaled "
           "threshold the delivery section uses, so a drug cannot read "
           "\"yes\" here and \"unlikely\" there. It answers whether the "
           "molecule can reach brain at all: for a CNS target that decides "
           "whether a candidate is usable, and for a peripheral target it "
           "decides whether brain exposure is an unwanted side effect rather "
           "than the goal."]
    if n_yes or n_borderline:
        out.append("")
        out.append("A **yes** predicts *passive diffusion only*, and is the "
                   "weaker of the two answers. Active efflux by "
                   "P-glycoprotein at the endothelium pumps an otherwise "
                   "suitable molecule straight back out; known substrates "
                   "are checked here and a hit overrides a pass outright.")
    if n_borderline:
        out.append("")
        out.append(f"**A question mark means undecided** — the score landed "
                   f"within noise of the threshold, and {n_borderline} of "
                   f"{total} row(s) here did. This is not pedantry. "
                   f"**Atenolol clears a complete five-parameter MPO by 0.04 "
                   f"and is a textbook poor penetrant**, so a margin that "
                   f"thin carries no information. Read a question-marked "
                   f"cell as \"look this one up\", never as the answer it "
                   f"leans toward.")
    if n_unknown:
        out.append("")
        out.append(f"**— is not a no.** {n_unknown} of {total} row(s) "
                   f"returned no physicochemical properties — common for "
                   f"DGIdb- and Open Targets-sourced rows, which carry a "
                   f"name and a mechanism but no structure. Nothing was "
                   f"measured and nothing is being claimed; look the "
                   f"compound up before excluding it on this column.")
    out.append("")
    out.append("A **no** is the more reliable direction. A molecule failing "
               "on size, polarity or lipophilicity does not cross by "
               "diffusion, and no formulation fixes that short of bypassing "
               "the barrier — intrathecal delivery, convection-enhanced "
               "delivery, or focused ultrasound.")
    out.append("")
    return out


def escape_pipes(value):
    """Return `value` with each `|` in a string escaped, copying only what changes.

    Values reach the report from services and from notes this package builds
    ("A --> B |" arrows, "eQTLGen: ... | eQTL Catalogue"). Inside a table row
    an unescaped pipe starts a new cell, so one such value shifts every cell
    after it and the row no longer lines up with its header. Escaping once,
    at the entry point, covers every table without touching the thousand
    f-strings that build them. Outside a table, Markdown renders an escaped
    pipe as a plain one, so prose is unaffected. Dictionary keys are left
    alone.

    A container is rebuilt only when something inside it changed, so a
    result with no pipes (the usual case) is returned as it is, with no
    copy of its sequences and route tables.
    """
    if isinstance(value, str):
        return value.replace("|", "\\|") if "|" in value else value
    if isinstance(value, dict):
        changed = False
        out = {}
        for k, v in value.items():
            nv = escape_pipes(v)
            changed = changed or nv is not v
            out[k] = nv
        return out if changed else value
    if isinstance(value, (list, tuple)):
        items = [escape_pipes(v) for v in value]
        if all(a is b for a, b in zip(items, value)):
            return value
        return items if isinstance(value, list) else tuple(items)
    return value


def _standard_of_care_section(result: dict) -> list:
    """Approved treatments for the condition, beside the ranking and apart
    from it. Empty when the gene has no entry in `standard_of_care`."""
    soc = result.get("standard_of_care")
    if not soc or not soc.get("treatments"):
        return []
    gene = soc["gene"]
    L = [f"## Approved treatments for {soc['condition']}", "",
         f"Already prescribed for the condition {gene} causes, listed so the "
         f"ranking below can be read against what exists. They are context, "
         f"not candidates: they are not part of any score here, and they are "
         f"not routes to {gene}.", "",
         "| Drug | Brand | Approved | For | Acts on | Mechanism |",
         "|---|---|---|---|---|---|"]
    for t in soc["treatments"]:
        L.append(f"| {t['generic']} | {t['brand']} | {t['approved']} | "
                 f"{t['use']} | {t['target_gene']} | {t['mechanism']} |")
    L += ["", f"{soc['note']} *Source: {soc['source']}. The table is "
              f"hand-maintained; confirm dates and labels against Drugs@FDA.*",
          ""]
    return L


def _brake_exon_section(result: dict) -> list:
    """Exons only decay-prone isoforms use, as targets to look at for an
    antisense oligonucleotide. Empty unless the goal is to raise the gene."""
    brake = result.get("brake_exons") or {}
    if not brake:
        return []
    gene = result["query"]["gene"]
    L = [f"## Exons that may brake {gene} output", ""]
    if brake.get("error"):
        L += [f"Not assessed: {brake['error']}. That is not a finding that "
              f"{gene} has no such exon.", ""]
        return L
    L += [f"Some genes make part of their transcript with an extra exon that "
          f"puts a stop codon early, so the message is destroyed. An antisense "
          f"oligonucleotide that makes the cell skip that exon turns those "
          f"transcripts into working ones and raises the protein, from the "
          f"working copy only. {brake['note']}", ""]
    cands = brake.get("candidates") or []
    novel = [c for c in cands if c["kind"] == "novel exon"]
    other = [c for c in cands if c["kind"] != "novel exon"][:5]
    if novel or other:
        L += ["| Exon (genomic) | Length | Kind | Decay isoforms using it |",
              "|---|---|---|---|"]
        for c in novel + other:
            L.append(f"| {c['start']}\u2013{c['end']} | {c['length']} bp | "
                     f"{c['kind']} | {c['n_transcripts']} |")
        L.append("")
        if len(cands) - len(novel) > len(other):
            L += [f"{len(cands) - len(novel) - len(other)} more extended or "
                  f"shifted coding exons are not shown; most are first-exon "
                  f"start variants, which are not brakes.", ""]
    L += ["**This is a list of places to look, not a design.** An annotated "
          "exon is one Ensembl has seen in some tissue. Whether it is "
          "included in your cells, how often, and whether skipping it raises "
          "protein is a question for RNA sequencing and a splicing assay. "
          "Nothing here picks an oligonucleotide sequence or checks it "
          "against the rest of the transcriptome.", ""]
    return L


def _mover_profile_section(result: dict) -> list:
    """What KIND of chemical the literature records as moving this gene.

    First, because it frames everything under it. When the published way
    to move a gene is to damage the tissue, the ranking below is a list of
    ways to attempt something the literature has only ever achieved by
    insult - and that is worth knowing before reading a drug name.
    """
    for blk in result.get("results") or []:
        prof = (blk.get("crosscheck") or {}).get("wanted_profile") or {}
        if not prof.get("verdict"):
            continue
        L = ["## What the literature says moves this gene", "",
             prof["verdict"], "",
             "| Kind of chemical | How many |", "|---|---|"]
        for k, n in sorted(prof.get("counts", {}).items(),
                           key=lambda kv: (-kv[1], kv[0])):
            L.append(f"| {k} | {n} |")
        L.append("")
        L.append("*The classes come from a keyword vocabulary over the CTD "
                 "chemical names, not from a toxicology database. An "
                 "unclassified chemical was not examined and is not thereby "
                 "safe, and a chemical can be both a medicine and a hazard "
                 "- phenobarbital is an approved barbiturate and a classic "
                 "rodent tumour promoter.*")
        L.append("")
        return L
    return []


def _proliferation_section(result: dict) -> list:
    """Did the candidates move the gene, or move the cell cycle?"""
    from . import proliferation as pr_mod
    for blk in result.get("results") or []:
        assessed = blk.get("proliferation") or {}
        if not assessed:
            continue
        gene = (blk.get("promoter") or {}).get("gene_symbol") or "the gene"
        text = pr_mod.summarise(assessed, gene)
        if not text:
            continue
        L = ["### Moving the gene, or moving the cell cycle", "", text, ""]
        tracked = [(k, v) for k, v in sorted(assessed.items())
                   if v.get("tracks_proliferation")]
        if tracked:
            L.append("| Compound | Signatures moving the programme | "
                     "Background |")
            L.append("|---|---|---|")
            for k, v in tracked[:12]:
                L.append(f"| {k} | {v['n_with_programme']} of "
                         f"{v['n_signatures']} ({v['share']:.0%}) | "
                         f"{v['background_share']:.1%} |")
            L.append("")
        return L
    return []


def _coverage_section(result: dict) -> list:
    """Which of the six steps to a protein this run actually looked at.

    Placed before the route tables because it says how much of the
    problem they cover. A full table is not evidence that the other five
    steps were considered, and nothing else in the report would say so.
    """
    from . import coverage as cov

    rows = cov.audit(result)
    if not rows:
        return []
    L = ["## How much of the problem this run examined", "",
         cov.summarise(rows), ""]
    L.append("| Step | This run | Candidates | What is missing |")
    L.append("|---|---|---|---|")
    for r in rows:
        gap = r["gap"] or ("—" if r["state"] == cov.SEARCHED else "")
        L.append(f"| {r['step']} | {r['state']} | "
                 f"{r['found'] if r['state'] == cov.SEARCHED else '—'} | "
                 f"{_clip(gap, 190)} |")
    L.append("")
    missing = [r for r in rows if r["state"] != cov.SEARCHED and r["enable"]]
    if missing:
        L.append("What would close each gap:")
        L.append("")
        for r in missing:
            L.append(f"- **{r['step']}** — {r['enable']}.")
        L.append("")
    L.append("*\"Listed only\" means the mechanism and its drug classes "
             "are named in this report and nothing was looked up for this "
             "gene. It is a menu, not a search.*")
    L.append("")
    return L


def _protein_abundance_section(result: dict) -> list:
    """Compounds measured against the PROTEIN rather than the transcript."""
    from . import proteomics as pr

    for blk in result.get("results") or []:
        hits = blk.get("protein_abundance")
        if hits is None:
            continue
        gene = (blk.get("promoter") or {}).get("gene_symbol") or "the gene"
        L = ["### Measured on the protein", "",
             pr.summarise(hits, gene), ""]
        if hits:
            L.append("| Compound | Change | Observations | Cell lines |")
            L.append("|---|---|---|---|")
            for h in hits[:12]:
                L.append(f"| {h['drug']} | "
                         f"{2 ** h['median_log2fc']:.1f}x | "
                         f"{h['n_observations']} | "
                         f"{', '.join(h.get('cell_lines') or []) or '—'} |")
            L.append("")
        return L
    return []


def _protein_level_section(result: dict) -> list:
    """Routes that change the protein rather than how much of it there is.

    Fires on what the caller DECLARED about the alleles, never on a guess:
    a hypomorphic or biallelic scenario, a premature stop, or a broken
    allele recorded as making a truncated product. Each of those means the
    quality of the product is part of the problem, and transcription
    cannot reach quality.
    """
    q = result.get("query") or {}
    scenario = q.get("scenario")
    has_ptc = q.get("has_ptc")
    truncated = q.get("broken_allele_makes") == "truncated"
    gene = q.get("gene") or "this gene"
    if not (scenario in ("hypomorph", "biallelic_hypomorph")
            or has_ptc or truncated):
        return []

    why = []
    if scenario == "biallelic_hypomorph":
        why.append("the only working copy is itself impaired")
    elif scenario == "hypomorph":
        why.append("both copies are impaired")
    if has_ptc:
        why.append("the transcript carries a premature stop")
    if truncated:
        why.append("the broken copy is recorded as making a truncated "
                   "product")

    L = ["## The problem may be the protein, not the amount", ""]
    L.append(f"**You declared that {', and '.join(why)}.** Everything in "
             f"the promoter sections below changes how MUCH {gene} is "
             f"made. None of it changes what is made. Where the pathogenic "
             f"allele yields an impaired product, more transcription "
             f"yields more of the same impairment — and because "
             f"transcription raises both copies, it raises the broken one "
             f"too.")
    L.append("")
    L.append("The routes that address the product itself are computed in "
             "this run and sit lower down:")
    L.append("")

    blk = next(iter(result.get("results") or []), {}) or {}
    beyond = blk.get("beyond_transcription") or {}
    pointers = []
    if blk.get("stabilisation"):
        pointers.append(
            f"**Stabilisation** — {len(blk['stabilisation'])} factor(s) "
            f"assessed for a molecule that holds the protein in its "
            f"working shape instead of making more of it.")
    if (beyond.get("product_stability") or {}):
        pointers.append(
            "**Product stabilisation** — whether the gene's own protein "
            "can be kept from being degraded.")
    if (beyond.get("post_transcriptional") or {}).get("routes"):
        n = len(beyond["post_transcriptional"]["routes"])
        extra = (" Readthrough and NMD inhibition are the ones that "
                 "address a premature stop, and for a stop codon they "
                 "address the problem where extra transcription does not."
                 if has_ptc else "")
        pointers.append(
            f"**After transcription** — {n} route(s) acting on the "
            f"transcript or the protein rather than the promoter.{extra}")
    if result.get("brake_exons", {}).get("candidates"):
        pointers.append(
            f"**Brake exons** — "
            f"{len(result['brake_exons']['candidates'])} exon(s) used only "
            f"by decay-prone isoforms, candidates for an oligonucleotide "
            f"that raises useful output rather than total output.")
    if blk.get("product_screen", {}).get("hits"):
        pointers.append(
            f"**Binders of the product** — "
            f"{len(blk['product_screen']['hits'])} compound(s) predicted "
            f"to bind the protein itself.")
    if blk.get("modality_design"):
        pointers.append(
            "**Designed modality** — a direct construct for this gene.")

    if pointers:
        L.extend(f"- {p}" for p in pointers)
    else:
        L.append("- None ran. `--design-modality`, `--screen-product` and "
                 "the stabilisation section (on by default) are the ones "
                 "that would, and `--has-ptc` moves readthrough and NMD "
                 "inhibition to the front.")
    L.append("")
    L.append("*This is a pointer, not a ranking. Nothing here establishes "
             "that a protein-level route works either; it establishes that "
             "it is aimed at the right thing.*")
    L.append("")
    return L


def _fame_section(result: dict) -> list:
    """How many other genes the literature records each candidate for.

    Placed with the other specificity cautions rather than beside the
    ranking, because it is a statement about what the ranking can tell
    you rather than about any one drug.
    """
    from . import fame as fame_mod

    for blk in result.get("results") or []:
        assessed = blk.get("fame")
        if assessed is None:
            continue
        table = {"n_genes": next(iter(assessed.values()))["of"]} \
            if assessed else None
        L = ["### How much of this order is fame?", "",
             fame_mod.summarise(assessed, table), ""]
        loud = sorted((k, v) for k, v in assessed.items()
                      if v.get("promiscuous"))
        if loud:
            L.append("| Compound | Genes it is recorded for |")
            L.append("|---|---|")
            for k, v in sorted(loud, key=lambda kv: -kv[1]["n_genes"])[:12]:
                L.append(f"| {v.get('name') or k} | "
                         f"{v['n_genes']} of {v['of']} "
                         f"({v['share']:.0%}) |")
            L.append("")
        return L
    return []


def _crosscheck_section(result: dict) -> list:
    """Where CTD and LINCS disagree, and ranked drugs CTD records the other way."""
    from . import crosscheck as xc_mod
    for blk in result.get("results") or []:
        xc = blk.get("crosscheck")
        if xc:
            ranked = (blk.get("all_routes") or {}).get("routes") or []
            return xc_mod.markdown(xc, ranked)
    return []


def _mitotic_gene_section(result: dict) -> list:
    """A short caution for a gene whose job is mitosis."""
    if not result.get("mitotic_gene"):
        return []
    from .cellcycle import MITOTIC_GENE_CAUTION
    return ["## A mitotic gene", "",
            MITOTIC_GENE_CAUTION.format(gene=result["query"]["gene"]), ""]


_MITOTIC_REPORT_NOTE = (
    "⚠halts mitosis: the compound is suspected to hold cells in mitosis "
    "instead of letting them proceed, and {gene} is a mitotic gene. A report "
    "that it raises {gene} (mRNA or protein) can reflect cells held in "
    "mitosis, where such genes are at their highest, and not a switch at the "
    "promoter.")


def markdown_report(result: dict, top_n: int = 25) -> str:
    result = escape_pipes(fda_mod.apply(result))
    q = result["query"]
    goal = q.get("goal", "up")
    verb = "Upregulation" if goal == "up" else "Downregulation"
    L: list[str] = []
    L.append(f"# {verb} strategy for {q['gene']} — promoter / TF / drug analysis")
    L.append("")
    L.append(f"Objective: **{'increase' if goal == 'up' else 'decrease'}** "
             f"transcription of {q['gene']}.")
    L.append("")
    L.extend(fda_mod.notice(result))
    if q["demo"]:
        L.append("> **DEMO MODE.** Synthetic sequence, illustrative matrices and a "
                 "hand-entered drug table. Not usable for real conclusions.")
        L.append("")

    # ------------------------------------------------ what has to happen?
    # Above every table, because the ranking cannot tell you whether the
    # plan needs a doubling or a halving - and in two scenarios whether
    # anything below can work at all.
    req = result.get("requirement") or {}
    dose = result.get("dosage") or {}
    if req or (dose and dose.get("assessed")):
        L.append("## What would have to happen")
        L.append("")
        if req.get("summary"):
            if req.get("feasible") is False:
                L.append("> ⛔ **Changing transcription cannot address "
                         "this.** The routes below are context, not "
                         "candidates.")
                L.append("")
            L.append(req["summary"])
            L.append("")
        if req.get("constraint_conflict"):
            # The scenario is declared and the constraint is measured, so
            # only this says when the two disagree. It is above the dosage
            # line because it changes how that line should be read.
            L.append(f"> ⚠ {req['constraint_conflict']}")
            L.append("")
        if dose.get("summary") and dose.get("assessed"):
            L.append(dose["summary"])
            L.append("")
            # Which release said so. gnomAD v2.1.1 and v4.1 give
            # different numbers for the same gene and ClinGen is
            # re-dated weekly, so the verdict above is only comparable
            # between runs that name their sources.
            cite = []
            for src in dose.get("sources") or []:
                rel = src.get("release")
                cite.append(f"{src['name']} "
                            + (f"{rel}" if rel
                               else "(release not stated in the file)"))
            if cite:
                L.append(f"*Read from {'; '.join(cite)}.*")
                L.append("")
        elif dose and not dose.get("assessed"):
            # Nothing came back. Printing the requirement with no
            # dosage line beside it reads as though the question was
            # asked and the answer was reassuring. constraint.assess
            # already writes the correct sentence; it was never shown.
            for note in dose.get("notes") or []:
                L.append(note)
                L.append("")

    # -------------------------------------------- can this run be read at all?
    # Printed before any ranking. A reader who sees only the tables below
    # cannot tell a measured negative from a lookup that broke, and every
    # table in this report is structurally identical in both cases.
    ev = result.get("evidence_verdict") or {}
    srcs = (result.get("evidence_sources") or {}).get("sources") or {}
    if ev or srcs:
        L.append("## Evidence available in this run")
        L.append("")
        # Before the table, not after it. If the caller capped the run
        # and the cap bound, every FAILED row below has an explanation
        # that has nothing to do with the sources - and a reader who
        # meets the table first has already drawn the wrong conclusion.
        spend = result.get("spend") or {}
        if spend.get("refused"):
            L.append(f"> ⚠ **The request budget bound this run.** "
                     f"{spend['refused']} request(s) were not sent after "
                     f"the ceiling of {spend['budget']} was reached, so "
                     f"some sources below failed because they were never "
                     f"asked. Rerun with a larger `--max-requests`, or "
                     f"without one: everything already fetched is cached, "
                     f"so the second run resumes rather than repeating.")
            L.append("")
        badge = {"not interpretable": "⛔ NOT INTERPRETABLE",
                 "degraded": "⚠ DEGRADED",
                 "interpretable": "✓ INTERPRETABLE"}.get(ev.get("level"), "")
        if badge:
            L.append(f"**{badge}** — {ev.get('headline', '')}")
            L.append("")
        for reason in ev.get("reasons") or []:
            L.append(f"- {reason}")
        if ev.get("reasons"):
            L.append("")
        if srcs:
            L.append("| Source | State | What it returned |")
            L.append("|---|---|---|")
            for name, s in srcs.items():
                # "empty" is marked rather than left to the detail
                # column: a reader scanning this table for `ok` would
                # otherwise count a source that answered "nothing" as
                # one that delivered.
                state = {"ok": "ok",
                         "empty": "ok, nothing found",
                         "FAILED": "**FAILED**",
                         "not requested": "not requested"}.get(s["state"],
                                                               s["state"])
                L.append(f"| {name} | {state} | {s['detail']} |")
            L.append("")
        failed = [s for s in srcs.values() if s["state"] == "FAILED"]
        if failed:
            L.append("**These tests could not be completed.** Asking and "
                     "getting no answer is not the same as asking and "
                     "being told no, so nothing below rests on them. "
                     "What each one would have told you:")
            L.append("")
            for s in failed:
                # The source's own detail when no impact is written for
                # it. A bullet that stops after the dash promises an
                # explanation and then withholds it, which reads worse
                # than saying the plain thing.
                why = (s.get("impact") or s.get("detail")
                       or "no impact recorded for this source")
                L.append(f"- **{s['name']}** — {why}")
            L.append("")
        missing = [s for s in srcs.values() if s["state"] == "not requested"]
        if missing:
            L.append("Not requested (enable to strengthen the run): "
                     + ", ".join(f"{s['name']} (`{s['enable_with']}`)"
                                 for s in missing if s.get("enable_with"))
                     + ".")
            L.append("")

    direct_hits = result.get("direct_therapies") or []
    if direct_hits:
        L.append(f"## Drugs targeting {q['gene']} directly")
        L.append("")
        L.append("These act on the gene, transcript or protein itself, skipping "
                 "the transcription-factor route entirely. Where one of these "
                 "exists it is almost always the better starting point: "
                 "sequence-directed agents affect your gene and essentially "
                 "nothing else.")
        L.append("")
        L.append("| Drug | Modality | Specificity | Raises protein? | Mechanism | "
                 "Approved | Source |")
        L.append("|---|---|---|---|---|---|---|")
        for t in direct_hits:
            flag = " ⚠withdrawn" if t["withdrawn"] else ""
            raises = "yes" if t["raises_protein"] else "activity only"
            L.append(f"| {t['drug']}{flag} | {t['modality']} | {t['specificity']} | "
                     f"{raises} | {t['mechanism']} | {t['approval_year'] or '—'} | "
                     f"{t['source']} |")
        L.append("")
        if any(not t["raises_protein"] for t in direct_hits):
            L.append("*\"activity only\" means the drug raises the function of protein "
                     "already present without increasing its abundance. If your readout "
                     "is a western blot rather than an activity assay, that distinction "
                     "decides whether the experiment appears to work.*")
            L.append("")
    elif result["query"].get("include_direct", True):
        L.append(f"## Drugs targeting {q['gene']} directly")
        L.append("")
        L.append("None found in the curated table or in ChEMBL/DGIdb for this "
                 "direction. The curated list is hand-maintained and goes out of "
                 "date; check FDA/EMA and clinicaltrials.gov before concluding "
                 "nothing exists.")
        L.append("")

    L.extend(_standard_of_care_section(result))
    L.extend(_mitotic_gene_section(result))
    L.extend(_coverage_section(result))
    L.extend(_protein_level_section(result))
    L.extend(_protein_abundance_section(result))
    L.extend(_mover_profile_section(result))
    L.extend(_proliferation_section(result))
    L.extend(_fame_section(result))
    L.extend(_crosscheck_section(result))
    L.extend(_brake_exon_section(result))

    # ------------------------------------- the headline answer, up front
    # This is the question people actually bring to the tool. It used to be
    # answerable only by assembling eight tables yourself.
    if result["results"]:
        first = result["results"][0]
        burdens = {}
        for g in (first.get("drug_grades") or []):
            burdens[g["drug_name"].lower()] = sys_mod.clinical_burden(
                g["drug_name"], notes=" ".join(g.get("notes") or []),
                target=g.get("tf_symbol")).to_dict()
        preds = pred_mod.build(
            first, goal=goal,
            has_chip=bool(q.get("chip_bed") or q.get("chip_api")),
            burdens=burdens)
        L.extend(pred_mod.render(preds, q["gene"], goal,
                                 verdict=result.get("evidence_verdict")))

    ov = result.get("tss_window_overlap") or []
    if ov:
        L.append("> **⚠ These promoter windows overlap each other.** "
                 + " ".join(o["note"] for o in ov[:2])
                 + " Counts and route lists below are therefore not "
                 "independent across the windows.")
        L.append("")

    for res in result["results"]:
        p = res["promoter"]
        L.append(f"## Promoter of {p['gene_symbol']} ({p['transcript_id']})")
        L.append("")
        L.append(f"- Transcript chosen: **{p['transcript_selection']}**")
        L.append(f"- Assembly: {p['assembly']}  |  Strand: {'+' if p['strand'] == 1 else '-'}")
        L.append(f"- TSS: chr{p['chrom']}:{p['tss']}")
        L.append(f"- Window: chr{p['chrom']}:{p['start']}-{p['end']} "
                 f"(-{p['upstream']} / +{p['downstream']}, {p['sequence_length']} bp)")
        L.append(f"- Motifs scanned: {q['n_motifs']} at p ≤ {q['pvalue']:g} per position")
        qs = res.get("scan_quality") or {}
        ww = qs.get("window_widened")
        if ww:
            if ww.get("kept"):
                L.append(
                    f"- **The window was widened.** The requested window "
                    f"(-{ww['from'][0]}/+{ww['from'][1]}) showed no signal above "
                    f"chance ({ww['sites_before']} sites), so it was rescanned at "
                    f"-{ww['to'][0]}/+{ww['to'][1]}, which gave {ww['sites_after']} "
                    f"sites against about {ww['chance_after']} expected by chance. "
                    f"Sites further out say less about this gene than core "
                    f"promoter sites do. `--no-widen` keeps the requested window.")
            else:
                L.append(
                    f"- The requested window (-{ww['from'][0]}/+{ww['from'][1]}) "
                    f"showed no signal, and rescanning at "
                    f"-{ww['to'][0]}/+{ww['to'][1]} did not help either "
                    f"({ww['sites_after']} sites against about "
                    f"{ww['chance_after']} by chance), so the requested window "
                    f"was kept.")
        if qs.get("motifs_unreachable"):
            L.append(
                f"- {qs['motifs_unreachable']} of those motifs were not "
                f"scanned: in this window their best possible match is more "
                f"common than p = {q['pvalue']:g}, so any hit would be a "
                f"chance match at a higher rate than stated. A looser "
                f"`--pvalue` scans them, with more chance hits.")
        L.append(f"- Non-redundant sites: {res['n_sites']} across {res['n_tfs']} TFs")
        L.append("")
        q2 = res.get("scan_quality") or {}
        if q2:
            if q2["verdict"] == "no signal":
                L.append(f"> **⚠ This scan found no signal.** {q2['advice']}")
                L.append(">")
                L.append("> Everything below is ranked from a list of chance "
                         "matches. Fix the scan before reading any of it.")
            elif q2["verdict"] == "weak":
                L.append(f"> **⚠ Weak scan.** {q2['advice']}")
            else:
                L.append(f"*Scan quality: {q2['advice']}*")
            L.append("")
            if q2.get("fdp_note"):
                # The enrichment ratio understates this badly: 1.2x reads as
                # mild signal but means ~5 of every 6 rows are chance.
                L.append(f"*{q2['fdp_note']}*")
                L.append("")

        nb = res.get("neighbours") or {}
        if nb.get("checked") and nb.get("genes"):
            L.append("#### What else is in this window")
            L.append("")
            L.append("Every other specificity check here asks whether a drug "
                     "or factor affects other genes. This asks the prior "
                     "question: is the sequence that was scanned uniquely "
                     "this gene's?")
            L.append("")
            L.append(nb["summary"])
            L.append("")
            L.append("| Gene | Biotype | Location | Strand | Relation |")
            L.append("|---|---|---|---|---|")
            for g in nb["genes"][:12]:
                L.append(f"| {g['gene']} | {g['biotype']} | "
                         f"chr{g['chrom']}:{g['start']}-{g['end']} | "
                         f"{'+' if g['strand'] == 1 else '-'} | "
                         f"{g['relation']} |")
            L.append("")
            for g in nb["genes"][:4]:
                if g["relation"] != "nearby":
                    L.append(f"- **{g['gene']}** — {g['note']}")
            if nb.get("n_antisense"):
                L.append("")
                L.append("*An antisense overlap matters more than it looks. "
                         "The direction logic in this report assumes a factor "
                         "acts on YOUR promoter; where the window lies inside "
                         "an antisense transcript, a site may drive that "
                         "instead, and antisense RNA can stabilise or repress "
                         "the sense gene — so the sign does not carry.*")
            L.append("")
            L.append("*This checks annotation, not sequence uniqueness. BLAT "
                     "the window if you need to know it appears nowhere else "
                     "in the genome.*")
            L.append("")

        nrs = [t for t in res["tfs"] if t.get("nuclear_receptor")]
        if nrs:
            L.append("### Ligand-activated (nuclear) receptors with sites here")
            L.append("")
            L.append("These are the directly druggable ones.")
            L.append("")
            L.append("| TF | Role | Sites | Best rel. score | Best site (vs TSS) | ChIP peaks | Score |")
            L.append("|---|---|---|---|---|---|---|")
            for t in nrs[:top_n]:
                L.append(f"| {t['symbol']} | {t['role']} | {t['n_sites']} | "
                         f"{t['best_rel_score']} | {t['best_site_sequence']} "
                         f"({_fmt_offset(t['best_site_tss_offset'])}) | "
                         f"{t['chip_peaks']} | {t['score']} |")
            L.append("")

        at = res.get("attrition") or {}
        if at.get("n_tfs_with_sites"):
            L.append("#### Why the drug tables below look the way they do")
            L.append("")
            L.append(f"**{at['n_tfs_with_sites']} factor(s) had predicted "
                     f"sites here. {at['n_with_drug']} of them have an "
                     f"approved drug"
                     + (f", and all {at['n_drugged_nuclear_receptor']} are "
                        f"nuclear receptors.**" if at.get("all_drugged_are_nr")
                        else f" ({at['n_drugged_nuclear_receptor']} of those "
                             f"are nuclear receptors).**"))
            L.append("")
            L.append("Roughly 48 of the ~1,600 human transcription factors "
                     "have a ligand-binding pocket, and the steroid, "
                     "retinoid, thyroid and vitamin-D families are most of "
                     "them. So the drug tables are not a ranking that "
                     "hormone drugs won — they are what is left after every "
                     "factor without a ligand has been dropped.")
            L.append("")
            lost = at.get("top25_without_drug") or []
            if lost:
                L.append(f"**{at.get('n_top25_without_drug')} of the top 25 "
                         f"factors by score have no drug and appear in no "
                         f"drug table:** " + ", ".join(lost[:12])
                         + (" …" if len(lost) > 12 else "") + ".")
                L.append("")
                L.append("That is where the gene-specific signal is. A "
                         "half-site like AGGTCA occurs every few kilobases "
                         "by chance, so nuclear receptors surface for almost "
                         "any promoter; these factors surfaced for *this* "
                         "one. The routes that can reach them are the "
                         "upstream-network, stabilisation and "
                         "protein-level sections, not the ligand tables.")
                L.append("")

        L.extend(_undrugged_regulators(res))

        L.append("### Top-ranked transcription factors")
        L.append("")
        has_dep = any(t.get("dependency") for t in (res.get("tfs") or []))
        dep_col = " Cells need it? |" if has_dep else ""
        dep_sep = "---|" if has_dep else ""
        L.append("| Rank | TF | Role | NR | Other genes it controls |"
                 + dep_col +
                 " Sites | Best rel. score | Offset | ChIP | "
                 "Approved drugs | Score |")
        L.append("|---|---|---|---|---|" + dep_sep +
                 "---|---|---|---|---|---|")
        any_regulon = False
        any_understudied = False
        for i, t in enumerate(res["tfs"][:top_n], 1):
            n = t.get("n_regulon")
            if n is None:
                reg = "not counted"
            else:
                any_regulon = True
                reg = f"~{n}" if "tier" in (t.get("regulon_provenance") or "") \
                    else str(n)
                if t.get("regulon_understudied"):
                    any_understudied = True
                    reg += " \u26a0"
            dep_cell = ""
            if has_dep:
                dep_cell = f" {_dep_label(t)} |"
            L.append(f"| {i} | {t['symbol']} | {t['role']} | "
                     f"{'yes' if t['nuclear_receptor'] else ''} | {reg} |"
                     + dep_cell +
                     f" {t['n_sites']} | "
                     f"{t['best_rel_score']} | {_fmt_offset(t['best_site_tss_offset'])} | "
                     f"{t['chip_peaks']} | {len(t['drugs'])} | {t['score']} |")
        L.append("")
        L.extend(_regulon_column_note(any_regulon, any_understudied,
                                      res.get("tfs") or []))

        L.extend(_depmap_section(res.get("dependencies") or {}, res))

        L.extend(_sensitivity_section(res.get("sensitivity") or {}))

        L.extend(_docking_section(res.get("docking") or [],
                                  res.get("docking_consensus") or {}))

        act = res.get("tf_activity") or {}
        if act.get("n_testable"):
            L.append("### Is the factor actually active? (inferred from its "
                     "targets)")
            L.append("")
            L.append("Everything above this point is a *sequence* argument: a "
                     "PWM says a factor could bind. This is a different "
                     "question asked of different data - do the genes that "
                     "factor is known to regulate actually move, in the "
                     "direction its mode of regulation predicts? That is the "
                     "VIPER idea (Alvarez et al. 2016), and it is the closest "
                     "this package gets to measuring a factor rather than "
                     "predicting one.")
            L.append("")
            L.append("**Why it belongs next to the table above.** At p=1e-4 "
                     "over a 6 kb window most motif calls are noise. The two "
                     "lines of evidence are independent, so a factor with "
                     "both a site and a coherently perturbed regulon is "
                     "supported in a way that neither alone can match - and "
                     "a strong site with a flat regulon is the expected "
                     "failure mode of the scan, not a surprise.")
            L.append("")
            L.append(act.get("summary", ""))
            L.append("")

            scanned = act.get("scanned") or []
            if scanned:
                L.append(f"#### Factors with a site here that could be "
                         f"tested")
                L.append("")
                L.append(f"{len(scanned)} of the "
                         f"{len(res.get('tfs') or [])} factors with a motif "
                         f"hit in this window have enough curated targets to "
                         f"test. The rest are not absent from the result - "
                         f"they were never testable, which is a gap in the "
                         f"regulon table and not evidence about the factor.")
                L.append("")
                L.append("| TF | Activity | NES | q | Targets measured | "
                         "Coverage | Signed | Confounded with |")
                L.append("|---|---|---|---|---|---|---|---|")
                for a in scanned[:top_n]:
                    conf = a.get("confounded_with")
                    j = a.get("confounded_jaccard")
                    cstr = (f"{conf} (J={j})" if conf and (j or 0) >= 0.3
                            else "—")
                    L.append(
                        f"| {a['tf']} | {a['direction']} | "
                        f"{a['nes'] if a['nes'] is not None else '—'} | "
                        f"{a['q_value'] if a['q_value'] is not None else '—'} | "
                        f"{a['n_measured']}/{a['n_regulon']} | "
                        f"{a['coverage']:.0%} | "
                        f"{a['n_signed']} | {cstr} |")
                L.append("")
            else:
                L.append("*No factor with a site in this promoter had enough "
                         "curated targets to test. With TRRUST that is the "
                         "common case rather than a failure: most of the "
                         "~1,600 human transcription factors have too few "
                         "curated targets to support any activity estimate.*")
                L.append("")

            top = [a for a in (act.get("top") or [])
                   if a["tf"] not in {s["tf"] for s in scanned}]
            if top:
                L.append(f"#### The strongest calls overall, site or not")
                L.append("")
                L.append("These have no motif hit in your promoter window, so "
                         "they are not routes to your gene. They are here "
                         "because they say what this tissue's regulatory state "
                         "looks like, which is the context any route has to "
                         "work inside.")
                L.append("")
                L.append("| TF | Activity | NES | q | Targets measured |")
                L.append("|---|---|---|---|---|")
                for a in top[:10]:
                    L.append(
                        f"| {a['tf']} | {a['direction']} | "
                        f"{a['nes'] if a['nes'] is not None else '—'} | "
                        f"{a['q_value'] if a['q_value'] is not None else '—'} | "
                        f"{a['n_measured']}/{a['n_regulon']} |")
                L.append("")

            L.append("#### What this cannot tell you")
            L.append("")
            L.append(f"- **It is a reference tissue average, not your "
                     f"sample.** The signature is {act['tissue']} from a "
                     f"healthy-tissue expression table, so it says nothing "
                     f"about a patient, a cell line, or a disease state.")
            L.append("- **Activity is not attributable to one factor.** "
                     "Regulons overlap, and co-regulated targets cannot "
                     "distinguish a factor from its neighbour. Where the "
                     "overlap is large the table names the other factor; the "
                     "signal is real and which of the two produced it is "
                     "not decided.")
            L.append("- **Nuclear receptors are the hardest case here.** A "
                     "ligand-activated receptor can be transcribed "
                     "everywhere while being active only where its ligand "
                     "reaches, so a quiet regulon does not mean the receptor "
                     "is absent - it means nothing is switching it on in "
                     "healthy tissue. For those factors, read a flat result "
                     "as uninformative rather than negative.")
            L.append("- **Direction needs a signed edge, and half of them are "
                     "not.** TRRUST records no mode for 46% of its human "
                     "edges, so those support a two-tailed test (did the "
                     "regulon move) and no directional one. A row with a "
                     "dash under NES is that case.")
            src = (act.get("regulon_source") or "").lower()
            if "collectri" not in src:
                L.append(f"- **Coverage is the binding constraint, and it is "
                         f"fixable.** {act['n_regulons']} factors had a "
                         f"regulon here, against roughly 1,600 human "
                         f"transcription factors, so most were never "
                         f"testable. CollecTRI covers about three times as "
                         f"many and signs every edge, which would widen this "
                         f"section and let more of the dash rows resolve to "
                         f"a direction. Fetch it with "
                         f"`import decoupler as dc; "
                         f"dc.get_collectri(organism='human')`, save it to "
                         f"`data/collectri.tsv`, and it is picked up "
                         f"automatically.")
            if not act.get("permutations"):
                L.append("- **The p-values are analytic and "
                         "anticonservative.** They assume targets are "
                         "independent, and co-expression is exactly the "
                         "structure being measured. Re-run with "
                         "`--activity-permutations 2000` before relying on "
                         "one.")
            L.append("")

        ad = res.get("activity_drugs") or {}
        if ad:
            n_pairs = sum(len(v) for v in ad.values())
            L.append("### Drugs that move a factor, found without any "
                     "record that they bind it")
            L.append("")
            L.append("Every drug table above came from asking a database "
                     "which compounds were assayed at the target. That "
                     "question cannot return anything for a factor nobody "
                     "has done chemistry on, which is most of them, and it "
                     "is why those tables fill up with nuclear receptors.")
            L.append("")
            L.append("This asks the LINCS data instead. Each compound in "
                     "LINCS has its own transcriptional signature, so the "
                     "same regulon statistic used above can be run on it: "
                     "does this compound move the factor's known targets? "
                     "No annotated binding is needed, so a factor with no "
                     "pharmacology at all can still return candidates.")
            L.append("")
            L.append(f"**{n_pairs} compound-factor pair(s) across "
                     f"{len(ad)} factor(s).**")
            L.append("")
            L.append("| Factor | Drug | Corpus | Effect on factor | "
                     "Median z | Profiles agreeing | Significant |")
            L.append("|---|---|---|---|---|---|---|")
            for tf, hits in ad.items():
                for h in hits[:10]:
                    L.append(
                        f"| {tf} | {h['drug']} | "
                        f"{h.get('corpus', 'lincs')} | {h['direction']} | "
                        f"{h['median_z']:+.2f} | "
                        f"{h['n_agreeing']}/{h['n_signatures']} "
                        f"({h['consistency']:.0%}) | "
                        f"{h['n_significant']} at q<0.05 |")
            L.append("")
            L.extend(_corpus_note(res))
            L.append("**Repetition is the evidence here, not any single "
                     "profile.** The scan tests one regulon against roughly "
                     "33,000 signatures, so nominal significance is reached "
                     "about 1,650 times by chance. Significance is counted "
                     "after Benjamini-Hochberg across those signatures, and "
                     "a compound has to clear three further bars: profiled "
                     "at least "
                     f"{tfa_mod.MIN_SIGNATURES} times, moving the regulon "
                     f"the same way in at least "
                     f"{tfa_mod.MIN_CONSISTENCY:.0%} of them, with a median "
                     f"effect of at least z {tfa_mod.MIN_MEDIAN_Z}. A "
                     "compound profiled once is never reported.")
            L.append("")
            L.append("#### What a row here does not mean")
            L.append("")
            L.append("- **It is not a binding claim.** What was measured is "
                     "that the factor's targets moved when the compound was "
                     "given. The effect can run through anything upstream "
                     "of the factor, and the compound may never touch it.")
            L.append("- **The attribution to this factor is the weak "
                     "step.** Regulons overlap, so a compound that moves a "
                     "shared set of targets will score against every factor "
                     "that claims them. The activity section above names "
                     "the overlapping factor for each one.")
            L.append("- **Proliferation is the confound to rule out "
                     "first.** A compound that arrests the cell cycle moves "
                     "every regulon enriched for cell-cycle genes, so a "
                     "growth inhibitor can look like a specific effect on "
                     "any such factor. If the agreeing compounds share an "
                     "antiproliferative mechanism, suspect that before "
                     "crediting the factor.")
            L.append("- **LINCS cell lines are not your tissue.** These are "
                     "mostly cancer lines at micromolar doses for six to "
                     "twenty-four hours, which is a different regime from a "
                     "therapeutic exposure.")
            L.append("")
            L.append("These rows enter the route tables as **tier 3**, "
                     "alongside curated mechanisms rather than with the "
                     "tier-2 measured ones. The observation is real; the "
                     "attribution to the factor is an inference, and tier 3 "
                     "is what that combination is worth.")
            L.append("")

        ar = res.get("all_routes") or {}
        if ar.get("routes"):
            L.append("### Every route that acts on this gene")
            L.append("")
            L.append(f"{ar['n_routes']} routes across "
                     f"{len(ar['by_category'])} mechanism categories, pooled from "
                     "every method in this run and sorted by **evidence tier "
                     "first**, score second.")
            L.append("")
            L.append("A tier-1 route with a mediocre score beats a tier-5 route "
                     "with a great one. These are not the same kind of claim and "
                     "averaging them would destroy the information you need.")
            L.append("")
            L.append("**Why established binders sit at tier 2, just under the "
                     "direct agents.** A drug already measured to bind one of "
                     "your factors has half its chain settled: the "
                     "drug-to-receptor step has assay data behind it, so the "
                     "pharmacology is de-risked and only the receptor-to-gene "
                     "step is a hypothesis. A curated tier-3 ligand is "
                     "annotation on *both* halves. But they stay below tier 1, "
                     "because a direct agent acts on your gene itself and "
                     "skips the receptor-to-gene inference entirely — a "
                     "stronger claim however well-characterised the binder "
                     "is. The cost of the indirect route is in the row: "
                     "everything else that factor regulates moves with your "
                     "gene.")
            L.append("")
            for tier_label, n in ar["by_tier"].items():
                L.append(f"- {n} × {tier_label}")
            L.append("")
            L.append("| Tier | Route | Drug | Acts on | Steps | Mechanism | "
                     "Str | Spec | Eff | Score |")
            L.append("|---|---|---|---|---|---|---|---|---|---|")
            for r in ar["routes"][:max(top_n, 30)]:
                a = r["axes"]
                L.append(f"| {r['tier']} | {r['category']} | {r['drug']} | "
                         f"{r['acts_on']} | {r['directness']} | "
                         f"{r['mechanism'][:45]} | {a['strength']} | "
                         f"{a['specificity']} | {a['efficiency']} | {r['score']} |")
            L.append("")
            L.append("Best in each category:")
            L.append("")
            for cat, r in ar["best_per_category"].items():
                L.append(f"- **{cat}** — {r['drug']} (tier {r['tier']}, "
                         f"score {r['score']})")
            L.append("")
            if ar.get("categories_absent"):
                L.append("*No route found in: "
                         + ", ".join(ar["categories_absent"])
                         + ". An absent category is usually a real gap in the "
                         "pharmacology, but check you enabled the relevant flag "
                         "(--lincs, --screen, --upstream-network).*")
                L.append("")
            L.append("*Steps = how far the drug acts from your gene: 0 on the "
                     "gene or transcript, 1 on a TF that binds the promoter, 2+ "
                     "upstream of that TF. Each step costs specificity.*")
            L.append("")

            guide = cat_mod.CATEGORY_GUIDE
            present = [c for c in ar["by_category"] if c in guide]
            if present:
                L.append("#### What these categories mean")
                L.append("")
                L.append("Several of these sound alike and behave very "
                         "differently. The distinction that most often decides "
                         "an experiment is whether a route changes protein "
                         "**abundance** or protein **activity** — a western "
                         "blot only sees the first.")
                L.append("")

                # Which drugs are actually in each category, and what assay
                # would see the effect. Printing only "(21 found)" left the
                # reader unable to connect a category to their own candidates.
                by_cat: dict = {}
                for r in ar["routes"]:
                    by_cat.setdefault(r["category"], []).append(r["drug"])

                L.append("| Category | What it changes | Your drugs in it | "
                         "Assay that sees it |")
                L.append("|---|---|---|---|")
                for c in present:
                    g = guide[c]
                    drugs = by_cat.get(c) or []
                    seen, uniq = set(), []
                    for d in drugs:
                        if d.lower() not in seen:
                            seen.add(d.lower())
                            uniq.append(d)
                    shown = ", ".join(uniq[:4])
                    if len(uniq) > 4:
                        shown += f" … +{len(uniq) - 4}"
                    L.append(f"| {c} | **{g.get('changes', '?')}** | "
                             f"{shown or '—'} | {g.get('assay', '—')} |")
                L.append("")
                if "direct: gene product" in present:
                    L.append("**Only one of these changes activity rather "
                             "than abundance: *direct: gene product*.** "
                             "Everything else works through transcription, so "
                             "it moves how much of your gene's product exists "
                             "and an RT-qPCR will see it. A ligand or "
                             "inhibitor on the product itself changes what "
                             "the protein already present does while its "
                             "level barely moves — which is why ivacaftor "
                             "reads as a failure on a western blot and works "
                             "in patients.")
                else:
                    L.append("*Every category present here works through "
                             "transcription, so all of them move your gene's "
                             "product **abundance** and an RT-qPCR is the "
                             "right first readout. The one exception in the "
                             "taxonomy — a drug on the gene product itself, "
                             "which changes activity only — did not appear "
                             "in this run.*")
                L.append("")

                for c in present:
                    g = guide[c]
                    drugs = by_cat.get(c) or []
                    seen, uniq = set(), []
                    for d in drugs:
                        if d.lower() not in seen:
                            seen.add(d.lower())
                            uniq.append(d)
                    L.append(f"**{c}** ({ar['by_category'][c]} found)")
                    L.append("")
                    if uniq:
                        L.append(f"- *Yours:* {', '.join(uniq[:8])}"
                                 + (f" … +{len(uniq) - 8} more"
                                    if len(uniq) > 8 else ""))
                    L.append(f"- {g['what']}")
                    L.append(f"- *Why:* {g['why']}")
                    L.append(f"- *Cost:* {g['cost']}")
                    L.append(f"- *Differs:* {g['differs']}")
                    L.append("")
                absent = [c for c in (ar.get("categories_absent") or [])
                          if c in guide]
                if absent:
                    L.append("Not found here, and what you would be missing:")
                    L.append("")
                    for c in absent:
                        L.append(f"- **{c}** — {guide[c]['what']}")
                    L.append("")

        ar_g = res.get("all_routes") or {}
        if ar_g.get("routes"):
            L.extend(_drug_glossary(ar_g["routes"], res))

        if res.get("drug_grades"):
            L.append("### Graded candidates")
            L.append("")
            L.append("Combined = geometric mean of TF confidence and drug quality. "
                     "A weak TF prediction caps the grade no matter how good the drug.")
            L.append("")
            L.append("| Grade | Combined | Drug | BBB? | TF | TF conf. | "
                     "Drug qual. | Potency | Protein select. | Gene select. | "
                     "Evidence | Flags |")
            L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
            bbb_seen: dict = {}
            for g in res["drug_grades"][:top_n]:
                c = g["components"]
                flags = "; ".join(g["notes"])[:70] or "—"
                v = dlv_mod.bbb_verdict(
                    g["drug_name"],
                    (g.get("evidence") or {}).get("properties"))
                bbb_seen[v["cell"]] = bbb_seen.get(v["cell"], 0) + 1
                L.append(f"| **{g['grade']}** | {g['combined']} | {g['drug_name']} | "
                         f"{v['cell']} | "
                         f"{g['tf_symbol']} | {g['tf_confidence']} | {g['drug_quality']} | "
                         f"{c['potency']} | {c['selectivity']} | "
                         f"{c.get('regulon', '—')} | {c.get('evidence_depth', '—')} | "
                         f"{flags} |")
            L.append("")
            L.extend(_bbb_footnote(bbb_seen))

        if res.get("cell_cycle_summary") or res.get("cell_cycle_blocked"):
            L.append("### Cell-cycle compatibility")
            L.append("")
            L.append(f"Tissue proliferation state: **{q.get('proliferation')}**")
            L.append("")
            if res.get("cell_cycle_summary"):
                L.append(res["cell_cycle_summary"])
                L.append("")
            if res.get("cell_cycle_blocked"):
                L.append("**Removed as incompatible:** "
                         + ", ".join(res["cell_cycle_blocked"]))
                L.append("")
                L.append("These are not ranked lower — they are removed. In "
                         "post-mitotic tissue, E2F/DREAM target promoters are "
                         "actively repressed rather than merely quiet, so a route "
                         "through a G1/S factor does not work at any dose.")
                L.append("")
            restricted = [(s, a) for s, a in (res.get("cell_cycle") or {}).items()
                          if a["duty"] < 1.0 and a["compatible"]]
            if restricted:
                L.append("| TF | Phase | Duty cycle | Effect multiplier |")
                L.append("|---|---|---|---|")
                for sym, a in restricted[:top_n]:
                    L.append(f"| {sym} | {a['phase']} | {a['duty']:.0%} | "
                             f"{a['multiplier']} |")
                L.append("")
                L.append("*A partial duty cycle dilutes the population-level "
                         "effect in asynchronous cells: per-cell efficacy can be "
                         "much higher than a bulk RT-qPCR suggests. Synchronise or "
                         "use a phase-resolved readout before calling it a "
                         "failure.*")
                L.append("")

        ups = res.get("upstream") or {}
        up_search = res.get("upstream_search") or {}
        if up_search.get("factors_searched") and not ups:
            L.append("### Upstream routes to factors with no direct drug")
            L.append("")
            L.append(
                f"No route found. {up_search['nodes_checked']} of "
                f"{up_search['nodes_reachable']} upstream node(s) reachable "
                f"from {up_search['factors_searched']} factor(s) were checked "
                f"for an approved drug that moves them the needed way, in "
                f"order of how well each path explains the factor, and none "
                f"has one."
                + (f" {up_search['lookups_failed']} drug lookup(s) failed, so "
                   f"this is weaker than it reads."
                   if up_search.get("lookups_failed") else "")
                + (f" {up_search['factors_truncated']} factor(s) reached the "
                   f"limit of {up_search['max_checked']} nodes with more "
                   f"reachable, so a route further out is not ruled out."
                   if up_search.get("factors_truncated") else "")
                + (" Every reachable node was checked against an index of "
                   "approved-drug targets."
                   if up_search.get("drug_index") else ""))
            L.append("")
        if ups:
            L.append("### Upstream routes to factors with no direct drug")
            L.append("")
            L.append("Most TFs have no approved ligand. These paths walk upstream "
                     "through signed causal edges — kinases, phosphatases, "
                     "receptors — to proteins that *are* druggable. Signs "
                     "multiply along the path, so an inhibitor of an inhibitor "
                     "counts as activating.")
            L.append("")
            L.append("Distance buys druggability and spends selectivity: a node "
                     "three hops up reaches your gene through one branch and "
                     "everything else through the others. Read the fan-out column.")
            L.append("")
            for sym, hits in list(ups.items())[:6]:
                L.append(f"**{sym}** (no direct pharmacology)")
                L.append("")
                L.append("| Upstream target | Hops | Fan-out | Path | Act on it | "
                         "Drugs | Score |")
                L.append("|---|---|---|---|---|---|---|")
                for h in hits:
                    bp = h["best_path"]
                    act = "activate" if h["net_sign"] == 1 else "inhibit"
                    names = ", ".join(d["name"] for d in h["drugs"][:3]) or "—"
                    L.append(f"| {h['node']} | {bp['hops']} | {bp['fan_out']} | "
                             f"{bp['description']} | {act} | {names} | "
                             f"{h['score']} |")
                L.append("")
                warned = [w for h in hits for w in h["warnings"]]
                if warned:
                    for w in sorted(set(warned))[:3]:
                        L.append(f"- ⚠ {w}")
                    L.append("")

        mechs = res.get("mechanisms") or {}
        if mechs:
            L.append("### Mechanism options: orthosteric vs allosteric")
            L.append("")
            L.append("Allosteric sites are not constrained to bind the endogenous "
                     "ligand, so they diverge across a receptor family where the "
                     "orthosteric pocket cannot. That usually buys subtype "
                     "selectivity — at the cost of a bounded maximal effect and a "
                     "hard dependence on endogenous tone.")
            L.append("")
            for sym, cmp in list(mechs.items())[:5]:
                L.append(f"**{sym}**")
                L.append("")
                L.append("| Mechanism | Drugs | Max effect | Subtype selectivity | "
                         "Outcompeted by ligand? | Needs endogenous ligand | "
                         "Examples |")
                L.append("|---|---|---|---|---|---|---|")
                for r in cmp["by_mechanism"]:
                    star = " ⭐" if r["mechanism"] in ("PAM", "NAM") else ""
                    L.append(f"| {r['mechanism']}{star} | {r['n_drugs']} | "
                             f"{r['max_effect']} | {r['subtype_selectivity']} | "
                             f"{'**yes**' if r.get('surmountable') else 'no'} | "
                             f"{'YES' if r['needs_endogenous_ligand'] else 'no'} | "
                             f"{', '.join(r['examples'])} |")
                L.append("")
                if cmp.get("surmountability"):
                    L.append(f"- {cmp['surmountability']}")
                    L.append("")
                if cmp.get("summary"):
                    L.append(f"*{cmp['summary']}*")
                    L.append("")
            L.append("**Before pursuing any allosteric option, verify that the "
                     "pathway has endogenous tone in your cell type.** A PAM with "
                     "no ligand to modulate does nothing at all, and that is a "
                     "measurement, not a prediction — nothing in this package can "
                     "check it for you.")
            L.append("")

        lax = res.get("ligand_axes") or {}
        if lax:
            L.append("### Changing the natural ligand instead of the receptor")
            L.append("")
            L.append("A nuclear receptor's activity is set by how much of its "
                     "endogenous ligand is present, and the enzymes that make "
                     "or destroy that ligand are often far more druggable than "
                     "the receptor. Anastrozole, finasteride and metyrapone all "
                     "work this way — none of them bind a receptor at all.")
            L.append("")
            L.append("| Strategy | Effect on the receptor |")
            L.append("|---|---|")
            L.append("| supplement ligand | up |")
            L.append("| inhibit catabolism | up |")
            L.append("| inhibit synthesis | down |")
            L.append("| sequester ligand | down |")
            L.append("")

            # Work out which notes every receptor repeats, and say them once.
            _allnotes = [set((a.get("notes") or [])[:4]) for a in lax.values()]
            _shared_ligand_notes = (set.intersection(*_allnotes)
                                    if len(_allnotes) > 1 else set())
            for n2 in sorted(_shared_ligand_notes):
                L.append(f"- {n2}")
            if _shared_ligand_notes:
                L.append("")

            for sym, a in list(lax.items())[:6]:
                L.append(f"#### {sym} — ligand: {a.get('ligand') or 'unknown'}")
                L.append("")
                if a.get("enzymes"):
                    L.append(f"Controlling enzymes: "
                             f"{', '.join(a['enzymes'])}")
                    L.append("")
                L.append("| Drug | Strategy | Acts on | Local? | Status |")
                L.append("|---|---|---|---|---|")
                for r in a["routes"][:top_n]:
                    yr = f" ({r['approval_year']})" if r["approval_year"] else ""
                    L.append(f"| {r['drug']}{yr} | {r['strategy']} | "
                             f"{r['acts_on']} | "
                             f"{'**yes**' if r['tissue_selective'] else 'no'} | "
                             f"{r['status']} |")
                L.append("")

                if a.get("chains"):
                    L.append("The chain, link by link:")
                    L.append("")
                    for c in a["chains"][:top_n]:
                        L.append(f"- {c}")
                    L.append("")

                per_drug = [(r["drug"], r["caveats"]) for r in a["routes"]
                            if r.get("caveats")]
                if per_drug:
                    L.append("Drug-specific caveats:")
                    L.append("")
                    for drug, cav in per_drug[:top_n]:
                        L.append(f"- **{drug}** — {cav[0]}")
                    L.append("")

                if a.get("opposite_direction"):
                    opp = ", ".join(
                        f"{o['drug']} ({o['strategy']})"
                        for o in a["opposite_direction"][:6])
                    L.append(f"*Same axis, opposite direction — do not confuse "
                             f"these with the list above: {opp}.*")
                    L.append("")

                # Notes shared by every receptor on this axis are
                # properties of the axis, not of the receptor, and printing
                # them per receptor is the same text three times.
                _own = [n2 for n2 in a["notes"][:4]
                        if n2 not in _shared_ligand_notes]
                for n2 in _own:
                    L.append(f"- {n2}")
                if _own:
                    L.append("")

            covered = lig_mod.covered_receptors()
            scanned = [t["symbol"] for t in (res.get("tfs") or [])][:25]
            without = [t for t in scanned if t not in covered]
            if without:
                L.append(f"**No ligand axis for:** {', '.join(without[:20])}"
                         + (" …" if len(without) > 20 else ""))
                L.append("")
                L.append("Most transcription factors are ligand-independent, so "
                         "this is usually a real absence rather than a gap in "
                         "the table. Receptors with a curated axis here: "
                         + ", ".join(f"{k} ({v})" for k, v in covered.items())
                         + ".")
                L.append("")

        stab = res.get("stabilisation") or {}
        if stab:
            L.append("### Locking the factor into an active state")
            L.append("")
            L.append("Nuclear receptors have a ligand-binding pocket; almost no "
                     "other TF does, which is why the graded tables keep "
                     "returning NR agonists. Stabilisation sidesteps the pocket: "
                     "every factor has protein-protein interfaces and "
                     "degradation machinery, and both are druggable.")
            L.append("")
            # Split by whether anything actually exists for the factor. The
            # previous version printed the same paragraph and the same four
            # exemplars under every pocket-less TF, which read as though a
            # route were available for each of them. For most, none is.
            actionable = {s: a for s, a in stab.items() if a.get("actionable")}
            nothing = [s for s, a in stab.items() if not a.get("actionable")]

            for sym, a in list(actionable.items())[:6]:
                pocket = ("has a ligand pocket" if a["has_ligand_pocket"]
                          else "**no ligand pocket**")
                L.append(f"**{sym}** — {pocket} · **a curated route exists**")
                L.append("")
                L.append("| Drug | Mechanism | Binds | Raises | Status |")
                L.append("|---|---|---|---|---|")
                for st in a["stabilisers"][:top_n]:
                    yr = f" ({st['approval_year']})" if st["approval_year"] else ""
                    L.append(f"| {st['drug']}{yr} | {st['mechanism']} | "
                             f"{st['target_node']} | {st['raises']} | "
                             f"{st['status']} |")
                L.append("")
                for n2 in a["notes"][:3]:
                    L.append(f"- {n2}")
                L.append("")

            if nothing:
                L.append(f"#### No stabilisation route for "
                         f"{len(nothing)} factor(s)")
                L.append("")
                L.append(", ".join(sorted(nothing)))
                L.append("")
                L.append("These have no ligand pocket, no direct drug, and "
                         "**nothing curated against them or their degradation "
                         "machinery**. For these factors stabilisation is a "
                         "research programme, not a route you can take: "
                         "making it real needs the E3 ligase that degrades "
                         "the factor *and* a ligand for that ligase, and "
                         "neither is published here.")
                L.append("")
                L.append("The mechanism class is real — it has shipped "
                         "approved drugs, listed once below — but those act "
                         "on **other proteins entirely**, and none of them "
                         "is evidence that any factor above is tractable.")
                L.append("")

            if any(a.get("exemplars") for a in stab.values()):
                ex_list = next(a["exemplars"] for a in stab.values()
                               if a.get("exemplars"))
                L.append("<details><summary>Approved proof the mechanism "
                         "class works — on unrelated targets</summary>")
                L.append("")
                for ex in ex_list[:4]:
                    L.append(f"- **{ex['drug']}** ({ex['approval_year']}) — "
                             f"{ex['mechanism']}: {ex['description']}")
                L.append("")
                L.append("</details>")
                L.append("")
            L.append("*Conformational stabilisers raise activity with protein "
                     "levels unchanged; degradation blockers raise levels. A "
                     "western blot cannot see the first kind — pick the assay to "
                     "match the mechanism.*")
            L.append("")

        ct = res.get("cell_type") or {}
        if ct:
            L.append("### Cell-type reach")
            L.append("")
            L.append("Three things must overlap: where the drug goes, where the "
                     "gene is expressed, and where the transcription factor is. "
                     "Any empty intersection kills the route regardless of motif "
                     "quality or potency.")
            L.append("")
            tfs = ct.get("tfs") or {}
            if tfs:
                L.append("| TF | Co-expressed with gene | Shared tissues |")
                L.append("|---|---|---|")
                for sym, ov in list(tfs.items())[:top_n]:
                    if not ov.get("known"):
                        verdict = "unknown"
                    elif ov.get("compatible"):
                        verdict = f"yes ({ov.get('n_shared')})"
                    else:
                        verdict = "**NO**"
                    L.append(f"| {sym} | {verdict} | "
                             f"{', '.join(ov.get('shared') or [])[:60] or '—'} |")
                L.append("")
                if any(o.get("compatible") is False for o in tfs.values()):
                    L.append("*Factors with no shared tissue were heavily "
                             "down-weighted: a hepatocyte factor cannot regulate "
                             "your gene in neurons, whatever its motif score.*")
                    L.append("")
            delivery_rows = ct.get("delivery") or []
            assessed = [d for d in delivery_rows if d.get("cns_mpo") is not None]
            unassessed = [d for d in delivery_rows if d.get("cns_mpo") is None]

            for d in assessed[:top_n]:
                verdict = ("likely" if d.get("likely_penetrant")
                           else "unlikely")
                L.append(f"- **{d['drug']}** — CNS MPO {d['cns_mpo']}/"
                         f"{d['mpo_parameters']}, {verdict} to reach CNS")
                # Notes belong to THIS drug. Printing them outside the
                # per-drug block dumped every failed assessment under
                # whichever drug happened to be printed last.
                seen_notes = []
                for n2 in d.get("notes", []):
                    if n2 not in seen_notes:
                        seen_notes.append(n2)
                for n2 in seen_notes[:2]:
                    L.append(f"  - {n2}")
            if assessed:
                L.append("")

            if unassessed:
                L.append(f"*{len(unassessed)} drug(s) could not be assessed for "
                         "CNS exposure — no structure or physicochemical "
                         "properties were returned for them, which is common "
                         "for entries sourced from DGIdb rather than ChEMBL. "
                         "They are neither passed nor failed on this axis: "
                         + ", ".join(d["drug"] for d in unassessed[:15])
                         + (" …" if len(unassessed) > 15 else "") + ".*")
                L.append("")

        lincs = res.get("lincs") or {}
        if lincs.get("hits"):
            L.append("### Measured expression changes (LINCS L1000)")
            L.append("")
            L.append("**This is the only empirical evidence in the report.** "
                     "Everything else predicts; this measured what happened to "
                     "the transcriptome when cells were dosed. Where it agrees "
                     "with a predicted route, you have a real candidate.")
            L.append("")
            L.append("| Drug | Effect | Weighted | Cell line | Dose | Time | "
                     "Evidence |")
            L.append("|---|---|---|---|---|---|---|")
            for h in lincs["hits"][:top_n]:
                kind = ("measured" if h["is_landmark"] is True
                        else "INFERRED" if h["is_landmark"] is False else "unknown")
                flag = " ⚠non-specific" if h["non_specific"] else ""
                L.append(f"| {h['drug']}{flag} | {h['score']} | "
                         f"{h['weighted_score']} | {h['cell_line'] or '—'} | "
                         f"{h['dose'] or '—'} | {h['timepoint'] or '—'} | {kind} |")
            L.append("")
            if lincs.get("summary"):
                L.append(lincs["summary"])
                L.append("")
            L.append("*L1000 measures ~978 landmark transcripts and infers the "
                     "other ~11,350. If your gene shows INFERRED, that number is a "
                     "model prediction, not a measurement — confirm by RT-qPCR "
                     "before relying on it. Check too that your gene is expressed "
                     "in the cell line listed, or absence of signal means "
                     "nothing.*")
            L.append("")

        gen = res.get("genetics") or {}
        if gen.get("summary") or (gen.get("cis") or {}).get("phenotypes"):
            L.append("### Human genetic evidence")
            L.append("")
            if gen.get("summary"):
                L.append(gen["summary"])
                L.append("")
            phen = (gen.get("cis") or {}).get("phenotypes") or []
            if phen:
                L.append("| Phenotype | Genetic score | Overall |")
                L.append("|---|---|---|")
                for pz in phen[:top_n]:
                    L.append(f"| {pz['phenotype']} | {pz['genetic_score']} | "
                             f"{pz['overall_score']} |")
                L.append("")
            # The actual trans-eQTL answer, not advice to go and find it.
            # This used to print the same "look for trans-eQTLs of X in
            # eQTLGen" paragraph once per factor — identical wording, no
            # result. The significant-results file is 2 MB.
            teq = gen.get("trans_eqtl") or []
            if teq:
                L.append("#### Does variation in the factor actually move "
                         "this gene? (trans-eQTL)")
                L.append("")
                L.append("The strongest evidence available anywhere in this "
                         "report. A motif says a factor *could* bind; a "
                         "trans-eQTL says natural genetic variation near it "
                         "**does** change your gene's expression, in real "
                         "human tissue, with exposure randomised by meiosis.")
                L.append("")
                if gen.get("trans_eqtl_summary"):
                    L.append(gen["trans_eqtl_summary"])
                    L.append("")
                supported = [t for t in teq if t.get("supported")]
                if supported:
                    L.append("| Factor | Variants | Best SNP | Z | p | n |")
                    L.append("|---|---|---|---|---|---|")
                    for t in supported[:top_n]:
                        h = (t.get("hits") or [{}])[0]
                        L.append(f"| **{t['tf']}** | {len(t.get('hits') or [])} | "
                                 f"{h.get('snp', '—')} | {h.get('z', 0):.1f} | "
                                 f"{h.get('p', 0):.2g} | {h.get('n', '—')} |")
                    L.append("")
                else:
                    tested = ", ".join(t["tf"] for t in teq[:8])
                    L.append(f"*No factor tested ({tested}) has trans-eQTL "
                             f"support for this gene.*")
                    L.append("")
                L.append("*eQTLGen is whole blood at n=31,684, FDR<0.05. A "
                         "hit is strong causal evidence. **Absence is weak** "
                         "— trans effects are small and detecting them needs "
                         "enormous power, so a null overwhelmingly means "
                         "underpowered rather than unregulated, and a "
                         "relationship operating outside blood would not "
                         "appear at all.*")
                L.append("")
            else:
                for t in (gen.get("trans") or [])[:3]:
                    if t.get("notes"):
                        L.append(f"- **{t['tf_symbol']} → {t['gene']}**: "
                                 f"{t['notes'][0]}")
                L.append("")

        an = res.get("analogs") or {}
        if an:
            L.append("### Approved drugs that resemble these candidates")
            L.append("")
            L.append(analog_mod.summarise(an))
            L.append("")
            L.append("Everything already annotated at the factor has been "
                     "subtracted. That matters: if a compound looks like a "
                     "known ligand and anyone assayed it there, it is "
                     "already in the tables above — so a raw similarity "
                     "search rediscovers its own input and looks productive.")
            L.append("")
            L.append("| Factor | Candidate | Resembles | Similarity |")
            L.append("|---|---|---|---|")
            for target, blk in an.items():
                for a2 in blk["analogs"][:6]:
                    L.append(f"| {target} | **{a2['drug']}** | "
                             f"{a2['similar_to']} | {a2['similarity']}% |")
            L.append("")
            L.append("*Similarity is not activity. Two molecules at this "
                     "Tanimoto can differ in the one substituent the pocket "
                     "cares about, and for a nuclear receptor a small change "
                     "can flip agonist to antagonist — so even if binding "
                     "holds, the direction is unknown. Each row is a binding "
                     "assay, not a finding.*")
            L.append("")

        bt = res.get("beyond_transcription") or {}
        if bt:
            L.append("### Routes a promoter analysis cannot see")
            L.append("")
            L.append("Everything above acts on transcription. These act "
                     "after it, around it, or on a different gene entirely. "
                     "For some genes that is where the whole problem is: if "
                     "the transcript is destroyed, the allele is silenced, "
                     "or the protein is degraded too fast, more "
                     "transcription does not help and may not be possible.")
            L.append("")

            comp = bt.get("compensation") or {}
            if comp.get("summary"):
                L.append("#### Could a relative do the job instead?")
                L.append("")
                L.append(comp["summary"])
                L.append("")
                for n2 in (comp.get("notes") or [])[:3]:
                    L.append(f"- {n2}")
                L.append("")
                k = comp.get("known_redundancy") or {}
                if k and k.get("paralog") and k["paralog"] != "(none functional)":
                    L.append(f"| | |")
                    L.append(f"|---|---|")
                    L.append(f"| paralog | **{k['paralog']}** |")
                    L.append(f"| why it substitutes | {k['why_it_substitutes']} |")
                    L.append(f"| why it is low | {k['silencing']} |")
                    L.append(f"| precedent | {k['precedent']} |")
                    L.append("")
                    for c in (k.get("caveats") or [])[:2]:
                        L.append(f"- {c}")
                    L.append("")

            react = bt.get("reactivation") or {}
            if react.get("applicable"):
                L.append("#### Is this a silenced allele rather than a low gene?")
                L.append("")
                for n2 in (react.get("notes") or [])[:2]:
                    L.append(f"- {n2}")
                L.append("")
                L.append("| Mechanism | Silencing | Status | Where selectivity comes from |")
                L.append("|---|---|---|---|")
                for r2 in (react.get("routes") or [])[:4]:
                    L.append(f"| {r2['mechanism']} | {r2['silencing_type']} | "
                             f"{r2['status']} | {r2['specificity'][:70]} |")
                L.append("")

            ptx = bt.get("post_transcriptional") or {}
            if ptx.get("routes"):
                L.append("#### After transcription")
                L.append("")
                if ptx.get("summary"):
                    L.append(ptx["summary"])
                    L.append("")
                L.append("| Mechanism | Raises | Status | Applies when |")
                L.append("|---|---|---|---|")
                for r2 in ptx["routes"][:6]:
                    L.append(f"| {r2['mechanism']} | {r2['raises']} | "
                             f"{r2['status']} | {r2['applies_when'][:58]} |")
                L.append("")
                lead = ptx["routes"][0]
                L.append(f"**{lead['mechanism'].capitalize()}** — "
                         f"{lead['what_it_does']}.")
                L.append("")
                for c in (lead.get("caveats") or [])[:2]:
                    L.append(f"- {c}")
                L.append("")

            ps = bt.get("product_stability") or {}
            if ps.get("routes"):
                L.append("#### Raising the protein instead of the transcript")
                L.append("")
                for n2 in (ps.get("notes") or [])[:2]:
                    L.append(f"- {n2}")
                L.append("")
                L.append("| Mechanism | Specificity | Status | Why |")
                L.append("|---|---|---|---|")
                for r2 in ps["routes"]:
                    L.append(f"| {r2['mechanism']} | {r2['specificity']} | "
                             f"{r2['status']} | {r2['specificity_note'][:62]} |")
                L.append("")
                L.append("*Three of these four raise most of the proteome, "
                         "which makes them excellent diagnostics and poor "
                         "therapies. The pharmacological chaperone is the "
                         "exception — and it needs a chaperone to exist for "
                         "your protein.*")
                L.append("")
                # Which ligase, and which enzyme reverses it. Generic
                # "E3 ligase inhibition" is not actionable until both
                # halves of the balance are named.
                from . import stabilizers as _st
                for _txt in (_st.summarise_degradation(ps.get("degradation")),
                             _st.summarise_rescue(ps.get("rescue"))):
                    if _txt:
                        L.append(_txt)
                        L.append("")

        imp = res.get("improvements") or {}
        if imp:
            L.append("### Candidates that may beat the current standard")
            L.append("")
            env = imp.get("baseline_envelope")
            base = imp.get("baseline")
            if env:
                L.append(f"**The bar: strength {env['strength']}, specificity "
                         f"{env['specificity']}, efficiency {env['efficiency']}** "
                         f"({env['evidence']}) — the best established value on "
                         f"*each axis*, not one drug's scores. Representative "
                         f"incumbent: {base['name']} ({base['route']}).")
                L.append("")
                L.append(f"{imp['n_previously_linked']} of {imp['n_candidates']} "
                         f"candidates already have a reported link to this gene "
                         f"(CTD returned {imp['ctd_links_found']} chemicals with "
                         f"curated expression effects).")
            else:
                L.append("**No incumbent found.** Nothing in the curated tables or "
                         "CTD currently controls this gene's expression. There is "
                         "no bar to clear, so candidates below are ranked on their "
                         "own merits — check the gap is real rather than a hole in "
                         "the databases.")
            L.append("")
            acc = imp.get("accounting") or {}
            if acc:
                L.append(f"*Accounting for all {acc.get('total')} candidates: "
                         f"{acc.get('previously_linked')} already linked, "
                         f"{acc.get('novel_and_on_front')} novel and non-dominated, "
                         f"{acc.get('novel_but_dominated')} novel but beaten on "
                         f"every axis by another candidate, "
                         f"{acc.get('novel_no_advantage')} with no advantage. "
                         f"All of them are in results.json under "
                         f"`improvements.all_candidates`.*")
                L.append("")

            rl = imp.get("reported_links") or {}
            rc = rl.get("counts") or {}

            def _mit(e, _on=bool(result.get("mitotic_gene"))):
                return (" ⚠halts mitosis" if _on
                        and sys_mod.is_mitotic_disruptor(
                    e.get("chemical", "")) else "")
            if rc.get("total"):
                L.append("#### Compounds already reported to change this gene")
                L.append("")
                L.append("These come from CTD's curated chemical–gene "
                         "interactions: someone measured an expression change "
                         "and a curator recorded it. That makes them the most "
                         "concrete evidence in this report — but the direction "
                         "is what decides whether each is a lead or a warning.")
                L.append("")
                L.append(f"{rc['total']} chemicals total: **{rc['wanted']}** "
                         f"move it the way you want, {rc['opposite']} move it "
                         f"the opposite way, {rc['unclear']} unclear.")
                L.append("")
                L.append("**Expect the 'predicted by this run?' column to say "
                         "no, and do not read that as a failed match.** CTD "
                         "indexes MeSH, which is a vocabulary for the "
                         "toxicology and pharmacology literature rather than "
                         "a drug list, so what comes back is three different "
                         "kinds of thing:")
                L.append("")
                L.append(f"- **chemical** ({rc.get('n_chemical', 0)}) "
                         f"{chr(8212)} a compound. The only kind that can "
                         f"match a candidate.")
                L.append(f"- **exposure** ({rc.get('n_exposure', 0)}) "
                         f"{chr(8212)} formaldehyde, asbestos, bisphenol, "
                         f"tobacco smoke. Real measured effects on your "
                         f"gene, and not things anyone can be given.")
                L.append(f"- **category** ({rc.get('n_category', 0)}) "
                         f"{chr(8212)} MeSH headings such as "
                         f"\"Anti-HIV Agents\" or \"Plant Extracts\". Not "
                         f"compounds at all, so the column is blank rather "
                         f"than no.")
                L.append("")
                matched = rc.get("matched_any_direction", 0)
                L.append(f"Across all three directions, **{matched} of "
                         f"{rc['total']} matched a candidate from this run**. "
                         + ("A low number here is the normal result: this "
                            "report's drug tables are approved pharmacology, "
                            "and most of CTD is not. It means the two "
                            "sources are looking at different parts of "
                            "chemical space, not that either is wrong."
                            if matched <= 5 else
                            "That overlap is the part of CTD this analysis "
                            "independently arrived at."))
                L.append("")

                want = rl.get("wanted_direction") or []
                if want:
                    L.append(f"**Reported to move {q['gene']} "
                             f"{'UP' if goal == 'up' else 'DOWN'}** — leads:")
                    L.append("")
                    L.append("| Chemical | Kind | Predicted by this run? | "
                             "CTD actions |")
                    L.append("|---|---|---|---|")
                    # Compounds first: they are the only rows where the
                    # match column carries information.
                    _order = {"chemical": 0, "exposure": 1, "category": 2}
                    for e in sorted(want,
                                    key=lambda x: _order.get(x.get("kind"), 3)
                                    )[:top_n]:
                        kind = e.get("kind", "chemical")
                        if kind == "category":
                            mark = f"{chr(8212)} (not a compound)"
                        elif e["in_candidates"]:
                            mark = "**yes**"
                        else:
                            mark = "no"
                        L.append(f"| {e['chemical']}{_mit(e)} | {kind} | "
                                 f"{mark} | {e['actions']} |")
                    if len(want) > top_n:
                        L.append(f"| … {len(want) - top_n} more | | | |")
                    L.append("")
                    if any(_mit(e) for e in want[:top_n]):
                        L.append(_MITOTIC_REPORT_NOTE.format(gene=q["gene"]))
                        L.append("")

                pl = imp.get("reported_protein_links") or {}
                pwant = [e for e in pl.get("wanted_direction") or []
                         if e.get("kind") == "chemical"]
                if pwant:
                    pc = pl.get("counts") or {}
                    L.append(f"**Reported to move {q['gene']} protein "
                             f"{'UP' if goal == 'up' else 'DOWN'}** "
                             f"({len(pwant)} compounds of {pc.get('total', 0)} "
                             f"with a protein report). The protein level can "
                             f"change through stability, translation or "
                             f"secretion with no change at the promoter, so "
                             f"these are not evidence that the promoter is "
                             f"involved. A compound that also appears above "
                             f"is reported at both levels.")
                    L.append("")
                    L.append("| Chemical | Predicted by this run? | "
                             "CTD actions |")
                    L.append("|---|---|---|")
                    for e in pwant[:top_n]:
                        L.append(f"| {e['chemical']}{_mit(e)} | "
                                 f"{'**yes**' if e['in_candidates'] else 'no'}"
                                 f" | {e['actions']} |")
                    if len(pwant) > top_n:
                        L.append(f"| … {len(pwant) - top_n} more | | |")
                    L.append("")
                    if any(_mit(e) for e in pwant[:top_n]):
                        L.append(_MITOTIC_REPORT_NOTE.format(gene=q["gene"]))
                        L.append("")

                only = rl.get("reported_but_not_predicted") or []
                if only:
                    L.append(f"**{len(only)} compounds have a measured "
                             f"effect on this gene and no mechanism from "
                             f"this analysis.** Exposures and MeSH headings "
                             f"are excluded from that count, because the "
                             f"scan was never going to propose formaldehyde "
                             f"and listing it as a miss says nothing.")
                    L.append("")
                    L.append("Most of what remains is research chemistry "
                             "rather than medicine: tool compounds, probes "
                             "and natural-product fractions that are "
                             "measured because somebody had them on a "
                             "shelf. Read them as evidence the gene IS "
                             "movable, and check the few that are approved "
                             "drugs individually.")
                    L.append("")
                    L.append(", ".join(e["chemical"] for e in only[:25])
                             + (" …" if len(only) > 25 else ""))
                    L.append("")

                opp = rl.get("opposite_direction") or []
                if opp:
                    L.append(f"**Reported to move it the WRONG way** — avoid, "
                             f"or use as negative controls:")
                    L.append("")
                    L.append(", ".join(e["chemical"] for e in opp[:25])
                             + (" …" if len(opp) > 25 else ""))
                    L.append("")
                L.append("*CTD aggregates across cell types, doses and assays. "
                         "A 'both' or contradictory annotation usually means "
                         "the effect is context-dependent rather than that the "
                         "data are wrong — check the underlying study before "
                         "relying on any single row.*")
                L.append("")

            novel = imp.get("novel_improvements") or []
            if novel:
                heading = ("#### Not previously linked, and better on at least "
                           "one axis" if imp.get("has_incumbent")
                           else "#### Not previously linked (no incumbent to beat)")
                L.append(heading)
                L.append("")
                L.append("| Candidate | Route | TF | Strength | Specificity | "
                         "Efficiency | Wins on | Trade-offs |")
                L.append("|---|---|---|---|---|---|---|---|")
                for c in novel[:top_n]:
                    a = c["axes"]
                    L.append(f"| {c['name']} | {c['route']} | {c['tf_symbol'] or '—'} | "
                             f"{a['strength']} | {a['specificity']} | {a['efficiency']} | "
                             f"{', '.join(c['wins_on'])} | "
                             f"{', '.join(c['loses_on']) or 'none'} |")
                L.append("")
                if imp.get("has_incumbent"):
                    L.append("*These are predicted scores compared against "
                             "measured ones, so a margin was required before "
                             "counting a win. They are hypotheses to test, not "
                             "improvements demonstrated.*")
                else:
                    L.append("*No incumbent exists, so these cleared no bar — they "
                             "are simply the candidates that scored above 0.5 on "
                             "at least one axis. Hypotheses to test.*")
            else:
                if imp.get("has_incumbent"):
                    L.append("#### No unlinked candidate beat the standard of care")
                    L.append("")
                    L.append("Every candidate that scored well is already a "
                             "reported regulator of this gene, or failed to clear "
                             "the per-axis bar by more than the prediction error. "
                             "This is the common result and it is informative: it "
                             "says the obvious pharmacology has been explored.")
                else:
                    L.append("#### No candidate scored above 0.5 on any axis")
                    L.append("")
                    L.append("There is no incumbent here, so nothing was filtered "
                             "out for failing to beat one — the candidates "
                             "themselves are simply weak. Widen the promoter "
                             "window or check whether this gene is regulated "
                             "through an enhancer outside it.")
            L.append("")

            dominated = imp.get("novel_dominated") or []
            if dominated:
                L.append("Novel but beaten on every axis by another candidate "
                         "(shown so nothing is hidden): "
                         + ", ".join(c["name"] for c in dominated[:12])
                         + (" ..." if len(dominated) > 12 else ""))
                L.append("")

            front = [c for c in (imp.get("pareto_front") or [])
                     if c.get("wins_on") or c.get("verdict", "").startswith("current")]
            if front:
                L.append("#### Non-dominated set (nothing beats these on every axis)")
                L.append("")
                for c in front[:top_n]:
                    linked = "known link" if c["previously_linked"] else "NEW link"
                    L.append(f"- **{c['name']}** ({linked}) — {c['verdict']}")
                L.append("")

        if res.get("outcome_estimates"):
            L.append("### Strength of evidence for each factor route")
            L.append("")
            L.append("Two links between a drug and your gene come from "
                     "evidence and differ from row to row: whether the drug "
                     "engages its target at a dose that can be reached, and "
                     "whether the factor is likely to act on this promoter. "
                     "Each is a 0 to 1 heuristic and not a calibrated "
                     "probability. They are shown as they are, and they are "
                     "not multiplied into a chance that the gene moves.")
            L.append("")
            # Anything named above (the discovery winner, the non-dominated
            # set, the incumbent) must appear here regardless of truncation.
            imp2 = res.get("improvements") or {}
            highlighted = {c["name"].lower() for c in
                           (imp2.get("novel_improvements") or [])
                           + (imp2.get("pareto_front") or [])}
            if imp2.get("baseline"):
                highlighted.add(imp2["baseline"]["name"].lower())

            ests = res["outcome_estimates"]
            keys = outcome_mod.shown_links(ests)
            pinned = [e for e in ests if e["drug"].lower() in highlighted]
            rest = [e for e in ests if e["drug"].lower() not in highlighted]
            shown = pinned + rest[:max(top_n - len(pinned), 0)]

            has_ko = any(e.get("knockout") is not None for e in ests)
            heads = ["Drug", "TF"] + [outcome_mod.LINK_LABEL[k].capitalize()
                                      for k in keys] + ["Weaker link"]
            if has_ko:
                heads.append("Knockout evidence")
            L.append("| " + " | ".join(heads) + " |")
            L.append("|" + "---|" * len(heads))
            for e in shown:
                serm = " ⚠tissue-dependent" if e["tissue_dependent"] else ""
                star = " ★" if e["drug"].lower() in highlighted else ""
                cells = [f"{e['drug']}{star}{serm}", e["tf_symbol"]]
                cells += [str(e[k]) for k in keys]
                cells.append(outcome_mod.weaker_link(e, keys))
                if has_ko:
                    cells.append(ko_mod.cell(e.get("knockout")))
                L.append("| " + " | ".join(cells) + " |")
            L.append("")
            if highlighted:
                L.append("★ = named in the comparison against the standard of "
                         "care above, pinned here so the two sections cannot "
                         "disagree.")
                L.append("")
            ko_lines = outcome_mod.knockout_caveat(ests)
            if ko_lines:
                for c in ko_lines:
                    L.append(c)
                    L.append("")
                roles = {t["symbol"]: t.get("role") for t in
                         (res.get("tfs") or [])}
                clash = sorted({
                    f"{e['tf_symbol']}: "
                    + ko_mod.role_conflict(e.get("knockout"),
                                           roles.get(e["tf_symbol"]))
                    for e in ests if ko_mod.role_conflict(
                        e.get("knockout"), roles.get(e["tf_symbol"]))})
                if clash:
                    L.append("Where the knockout data and the annotated role "
                             "disagree, one of them is wrong for this gene, "
                             "and it matters whether the factor is to be "
                             "activated or blocked:")
                    L.append("")
                    for c in clash:
                        L.append(f"- {c}")
                    L.append("")
            L.append("Order candidates by tier and score, as the tables above "
                     "do. These columns say where the evidence behind a "
                     "route is thin.")
            L.append("")
            if any(e["tissue_dependent"] for e in res["outcome_estimates"]):
                L.append("*Rows marked tissue-dependent sit in the SERM regime: "
                         "tamoxifen antagonises ERalpha in breast and activates it "
                         "in endometrium. For those, direction is a property of "
                         "your cell type, not of the molecule, and no structural "
                         "model can settle it — a coactivator recruitment assay can.*")
                L.append("")

        all_screens = res.get("screens") or {}
        if all_screens:
            ran = [s for s, v in all_screens.items() if not v.get("error")]
            skipped = [s for s, v in all_screens.items() if v.get("error")]
            L.extend(_screen_scope_note(ran, skipped))

        for sym, sc in all_screens.items():
            L.append(f"### Virtual screen of approved drugs against {sym}")
            L.append("")
            if sc.get("error"):
                L.append(sc["error"])
                L.append("")
                continue
            enr = sc.get("enrichment", {})
            L.append(f"Binding model, CV AUC {sc.get('model_auc')} — trained on "
                     f"{sc.get('n_actives')} measured actives vs "
                     f"{sc.get('n_decoys')} decoys, screened "
                     f"{sc.get('library_size')} approved drugs.")
            if sc.get("direction_note"):
                L.append("")
                L.append(f"*{sc['direction_note']}*")
            if not enr.get("usable", False):
                L.append("")
                L.append(f"**Screen rejected:** {enr.get('reason')} "
                         "No candidates are reported from it.")
                L.append("")
                continue
            L.append("")
            L.append("| Quadrant | Drug | What it is | Approved to treat | "
                     "Measured pChEMBL | P(binds), model | Novelty |")
            L.append("|---|---|---|---|---|---|---|")
            for h in sc["hits"][:top_n]:
                pch = h.get("measured_pchembl")
                prob = h.get("probability")
                # A model score on a molecule the model was trained on is an
                # in-sample fit. It is shown, because hiding it would hide
                # the enrichment diagnostic, but marked so it is not read as
                # a prediction and not mistaken for the ranking key.
                prob_cell = (f"{prob} (in-sample)"
                             if h.get("probability_in_sample") else prob)
                L.append(f"| {h['quadrant']} | {h['drug']} | "
                         f"{h.get('indication_hint') or '—'} | "
                         f"{h.get('treats') or 'not established here'} | "
                         f"{pch if pch is not None else '—'} | "
                         f"{prob_cell if prob is not None else '—'} | "
                         f"{h['novelty'] if h.get('novelty') is not None else '—'} |")
            L.append("")
            L.append("Two independent axes here. **ESTABLISHED LIGAND, "
                     "untested on this gene** is the most actionable class: "
                     "the drug reliably binds this receptor and has no "
                     "reported effect on your gene, so the pharmacology is "
                     "already de-risked and only the receptor-to-gene step "
                     "is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses "
                     "on both counts, novel chemistry and an untested gene "
                     "link, so expect most to be wrong and treat each as "
                     "needing a binding assay first.")
            L.append("")
            L.append("**Read the indication column before the scores.** A "
                     "receptor only reaches this section when it has a "
                     "large ligand series, and such series are rarely "
                     "homogeneous: a topical gel, an oral oncology agent "
                     "and a withdrawn compound can sit on adjacent rows "
                     "with near-identical numbers. The scores rank binding "
                     "to the receptor and say nothing about whether the "
                     "drug could be given for your purpose.")
            L.append("")
            L.append("**The two score columns are not interchangeable, and "
                     "each quadrant has only one real number.** An "
                     "established ligand was a training molecule for this "
                     "model, so its model score is an in-sample fit and its "
                     "novelty is structurally zero, because the nearest "
                     "known ligand the lookup finds is itself. Those rows "
                     "are ranked by measured pChEMBL, which is data. "
                     "Predicted rows have no measurement, so they are ranked "
                     "by the model score, which is a prediction. Comparing a "
                     "pChEMBL to a probability compares a measurement to a "
                     "guess.")
            L.append("")

        du = res.get("direction_unknown") or {}
        if du:
            n = sum(len(v) for v in du.values())
            L.append("### Drug–gene links with no recorded direction")
            L.append("")
            L.append(f"{n} interaction(s) across {len(du)} factor(s) come from "
                     "curated sources that record a drug–gene relationship "
                     "**without saying which way it goes**. They are held back "
                     "from the recommendations above, because the entire "
                     "question here is direction and an undirected link cannot "
                     "answer it.")
            L.append("")
            L.append("They are still leads — someone curated a real "
                     "relationship — but they need a literature check, not a "
                     "prediction. Many are chemotherapeutics or broad stressors "
                     "whose association with a transcription factor reflects "
                     "cell-wide perturbation rather than a specific mechanism.")
            L.append("")
            L.append("| TF | Drug | Reported by | DGIdb score |")
            L.append("|---|---|---|---|")
            rows = [(sym, d) for sym, ds in du.items() for d in ds]
            rows.sort(key=lambda r: -(r[1].get("interaction_score") or 0))
            for sym, d in rows[:top_n]:
                score = d.get("interaction_score")
                L.append(f"| {sym} | {d['name']} | {d.get('sources') or '—'} | "
                         f"{round(score, 2) if score else '—'} |")
            if len(rows) > top_n:
                L.append(f"| … {len(rows) - top_n} more | | | |")
            L.append("")
            L.append("*Pass `--include-unknown-direction` to fold these into "
                     "the main tables anyway — but every such row will read "
                     "UNKNOWN / unclear, which is what they are.*")
            L.append("")

        L.append("### Approved drugs acting on these factors")
        L.append("")
        rows = []
        for t in res["tfs"]:
            for d in t["drugs"]:
                direction = drugs_mod.desired_direction(t["role"] or "unknown", goal)
                want = {"up": "increase", "down": "decrease"}.get(direction, "either")
                rows.append((t["score"], t["symbol"], t["role"], d, want))
        if not rows:
            kind = "activating" if goal == "up" else "blocking"
            L.append(f"No approved drug with {kind} pharmacology was found for any TF "
                     "matched in this promoter. That is the usual result outside the "
                     "nuclear-receptor family; see the README for indirect strategies.")
        else:
            L.append("| TF | TF role | Do what to this TF | Drug | Action | "
                     "Drug direction | First approval | Source |")
            L.append("|---|---|---|---|---|---|---|---|")
            # Untruncated, this ran to 850 rows per promoter — 2,550 lines of
            # a 3,943-line report, and unreadable at any length.
            ordered = sorted(rows, key=lambda r: -r[0])
            for _, sym, role, d, want in ordered[:top_n]:
                flag = " ⚠withdrawn" if d.get("withdrawn") else ""
                L.append(f"| {sym} | {role} | {want} | {d['name']}{flag} | "
                         f"{d['action_type']} | {d['direction']} | "
                         f"{d.get('first_approval') or '—'} | {d['source']} |")
            if len(ordered) > top_n:
                L.append(f"| … {len(ordered) - top_n} more, see "
                         f"tf_summary.tsv | | | | | | | |")
        L.append("")

    # Experimental design is deferred to the end of the report.
    # It is how you would TEST any of the above, and putting a
    # CRISPRa protocol between the reader and the drug list buried
    # the answer under the method.
    lab: list[str] = []
    # These two sections describe the GENE, not a particular promoter window,
    # so they are emitted once. Repeating them per TSS tripled the report.
    for res in result["results"][:1]:
        mo = res.get("modality_design") or {}
        if mo:
            lab.append(f"## Designing a direct modality for {mo['gene']}")
            lab.append("")
            d = mo["delivery"]
            flag = "" if mo["feasible_in_vivo"] else " ⚠"
            lab.append(f"**Delivery to {d['compartment']}: {d['status']}**{flag} — "
                     f"{d['route']}")
            lab.append("")
            lab.append(f"{d['note']}")
            lab.append("")
            lab.append(f"**{mo['verdict']}**")
            lab.append("")
            for opt in mo["options"]:
                lab.append(f"### {opt['modality']}")
                lab.append("")
                lab.append(opt["rationale"])
                lab.append("")
                if opt.get("note"):
                    lab.append(f"*{opt['note']}*")
                    lab.append("")
                cands = opt.get("candidates") or []
                # CRISPRi guides also carry a "sequence" key; the oligo branch
                # is distinguished by having an antisense strand.
                if cands and "antisense" in cands[0]:
                    lab.append("| Target site (5'→3') | Guide/antisense | Pos | "
                             "GC | Notes |")
                    lab.append("|---|---|---|---|---|")
                    for c in cands[:6]:
                        w = "; ".join(c["warnings"])[:45] or "clean"
                        lab.append(f"| `{c['sequence']}` | `{c['antisense']}` | "
                                 f"{c['position']} | {c['gc']:.0%} | {w} |")
                    lab.append("")
                elif cands:
                    lab.append("| Protospacer | PAM | Strand | vs TSS | Notes |")
                    lab.append("|---|---|---|---|---|")
                    for c in cands[:6]:
                        w = "; ".join(c["warnings"])[:40] or "clean"
                        lab.append(f"| `{c['sequence']}` | {c['pam']} | "
                                 f"{c['strand']} | {c['tss_offset']:+d} | {w} |")
                    lab.append("")
            lab.append(f"*{mo['off_target_warning']}*")
            lab.append("")

        sp = next((r.get("screen_plan") for r in result["results"]
                   if r.get("screen_plan")), None)
        if not sp:
            continue
        lab.append("## Screen design")
        lab.append("")
        if len(result["results"]) > 1:
            lab.append(f"*One plan for the gene. The reporter insert below uses "
                     f"the {result['results'][0]['promoter']['transcript_id']} "
                     f"window; clone a second construct for the other TSS if "
                     f"you want to compare them.*")
            lab.append("")
        lab.append("The outcome chain above is almost always limited by *gene "
                 "responsiveness in your cell type* — the link this package "
                 "cannot see. A screen measures it instead of assuming it. "
                 "Everything above prioritises this experiment; it does not "
                 "replace it.")
        lab.append("")
        if sp.get("protein"):
            lab.extend(_protein_screen_block(sp["protein"]))
        elif sp.get("reporter"):
            lab.append("### Reporter construct")
            lab.append("")
            lab.append(f"- Insert: {sp['reporter']['insert']}")
            lab.append(f"- Vector: {sp['reporter']['vector']}")
            lab.append("- Controls:")
            for c in sp["reporter"]["controls"]:
                lab.append(f"  - {c}")
            lab.append(f"- {sp['reporter']['caveat']}")
            lab.append("")
        cr = sp["crispra"]
        lab.append("### CRISPRa screen")
        lab.append("")
        lab.append(f"- Approach: {cr['approach']}")
        lab.append(f"- Focused library: {cr['library']['focused']}")
        lab.append(f"- Unbiased: {cr['library']['genome_wide']}")
        lab.append(f"- Coverage: {cr['library']['coverage']}")
        lab.append(f"- Readout: {cr['readout']}")
        lab.append("")
        gl = cr["guides_at_target_promoter"]
        if gl:
            lab.append(f"Positive-control sgRNAs at your promoter "
                     f"({cr['n_usable_guides']} clean of {len(gl)} shown):")
            lab.append("")
            lab.append("| Protospacer (20 nt) | PAM | Strand | vs TSS | GC | Note |")
            lab.append("|---|---|---|---|---|---|")
            for g in gl[:10]:
                note = g["warnings"][0][:45] if g["warnings"] else "ok"
                lab.append(f"| `{g['sequence']}` | {g['pam']} | {g['strand']} | "
                         f"{g['tss_offset']:+d} | {g['gc']:.0%} | {note} |")
            lab.append("")
        lab.append(f"*{cr['note']}*")
        lab.append("")
        cp = sp["compounds"]
        lab.append("### Compound screen")
        lab.append("")
        for lib in cp["libraries"]:
            lab.append(f"- {lib}")
        lab.append("")
        lab.append(f"- Design: {cp['design']}")
        lab.append(f"- Counter-screen: {cp['counter_screen']}")
        lab.append("")
        if cp["cherry_pick"]:
            n_shown = len(cp["cherry_pick"][:24])
            lab.append("#### Cherry-pick list")
            lab.append("")
            lab.append("**What this is.** \"Cherry-picking\" means having the "
                     "vendor pull these named compounds into a custom plate "
                     "rather than buying a whole library. It is a cost and "
                     "time decision: ~48 wells you can run in a week, against "
                     "~1,520 compounds for Prestwick or ~3,000 for the "
                     "Selleck set.")
            lab.append("")
            lab.append("**Where these came from.** Every compound below was "
                     "proposed earlier in *this run* — nothing external was "
                     "consulted to build the list:")
            lab.append("")
            for src in cp.get("pulled_from", []):
                lab.append(f"- {src}")
            lab.append("")
            if cp.get("priority_key"):
                lab.append("Priority means:")
                lab.append("")
                for p, meaning in sorted(cp["priority_key"].items()):
                    lab.append(f"- **{p}** — {meaning}")
                lab.append("")
            if cp.get("cannot_do"):
                lab.append(f"**What it cannot do.** {cp['cannot_do']}")
                lab.append("")
            total = len(cp["cherry_pick"])
            lab.append(f"Showing {n_shown} of {total} selected "
                     f"(from {cp.get('n_available', total)} candidates in "
                     f"this run), in priority order:")
            lab.append("")
            lab.append("| Priority | Compound | Why it is on the plate | Pulled from |")
            lab.append("|---|---|---|---|")
            for c in cp["cherry_pick"][:24]:
                lab.append(f"| {c.get('priority', '—')} | {c['compound']} | "
                         f"{c['rationale']} | {c.get('source', '—')} |")
            lab.append("")
        lab.append(f"*{cp['note']}*")
        lab.append("")
        lab.append("### What to expect")
        lab.append("")
        for k, v in sp["expectations"].items():
            lab.append(f"- **{k.replace('_', ' ')}**: {v}")
        lab.append("")

    L.append("---")
    L.append("")
    L.append("### How to read this")
    L.append("")
    if goal == "down":
        L.append("**For downregulation, check the direct modalities first.** "
                 "Approved siRNA (patisiran, givosiran, inclisiran, vutrisiran) and "
                 "antisense oligonucleotides (inotersen, nusinersen, tofersen) knock "
                 "down a named transcript with sequence-level specificity. A TF-level "
                 "intervention hits every other gene that factor controls. Use this "
                 "report when a direct modality is not an option, or to understand "
                 "*why* a drug you already have lowers your gene.")
        L.append("")
    L.append("A PWM hit means the sequence *resembles* a binding site. It is a "
             "hypothesis generator, not evidence of binding, and not evidence that "
             "binding would increase transcription of this gene in your cell type. "
             "Rank order here is a heuristic (weights printed in the JSON output), "
             "not a validated predictor. Every row needs wet-lab confirmation: "
             "promoter-reporter assay for the element, agonist dose-response with "
             "RT-qPCR/protein readout for the drug, and a site-directed mutant of the "
             "element to show the drug effect runs through it.")
    L.append("")

    # ------------------------------------- experimental design, at the end
    # Deferred deliberately. These sections are how you TEST the analysis,
    # and they are long: a CRISPRa protocol and a plate layout sitting
    # between the reader and the drug list buried the answer under the
    # method. Everything above is what the model concluded; everything
    # below is how you would find out whether it is right.
    if lab:
        L.append("---")
        L.append("")
        L.append("# Testing any of this in the lab")
        L.append("")
        L.append("Everything above is prediction. This part is experiment — "
                 "how you would find out whether any of it holds, and it is "
                 "placed last because it is long, not because it is "
                 "optional.")
        L.append("")
        L.append("**If you run one thing, run the reporter assay with the "
                 "scrambled-motif control.** A drug that induces the "
                 "wild-type promoter but not the mutant has an effect "
                 "running through the element this report predicted. A drug "
                 "that induces both means the prediction was wrong. That "
                 "single control turns a hit into a mechanism.")
        L.append("")
        L.extend(lab)
    return "\n".join(L)


# ===================================================================== #
# The short version.                                                    #
# ===================================================================== #
# report.md is long on purpose: every number in it is qualified, because
# an unqualified number from this pipeline would be read as stronger than
# it is. That makes it the wrong document to open when you want to know
# what to do next.
#
# This is the other document. It carries the shortlist, how much the
# evidence is worth, and what would change the answer - and it points at
# report.md for everything else rather than restating it in shorter and
# therefore less qualified words. The two rules it follows: say nothing
# here that report.md does not say, and never drop a caveat to make a row
# fit on one line.

SUMMARY_TOP = 10
# How many unique routes are read before folding, per row shown, so folding a
# family still leaves enough distinct drugs to fill the list.
SUMMARY_POOL = 8


def _clip(text: str, limit: int = 58) -> str:
    """Trim to a word boundary. A mid-word cut ("renal fai") reads as a
    rendering bug and costs the reader the one fact the column carries."""
    text = (text or "").strip()
    if len(text) <= limit:
        return text
    cut = text[:limit].rsplit(" ", 1)[0].rstrip(",;:")
    return (cut or text[:limit]) + "…"


def _verdict_line(verdict: dict) -> str:
    """How much the ranking is worth, in one sentence."""
    level = (verdict or {}).get("level", "unknown")
    return {
        "strong": "**Evidence: strong.** Measured evidence supports this "
                  "ranking from more than one independent source.",
        "usable": "**Evidence: usable.** At least one measured source "
                  "backs the ranking; read the tiers before acting.",
        "degraded": "**Evidence: degraded.** Read the ranking as a reading "
                    "order, not as evidence. A pillar of the analysis is "
                    "missing.",
        "not interpretable": "**Evidence: NOT INTERPRETABLE.** Too much is "
                             "missing for the ranking to mean anything. "
                             "Fix the inputs before reading the rows.",
    }.get(level, f"**Evidence: {level}.**")


# ===================================================================== #
# One row per drug: which methods found it, and what it should do.       #
# ===================================================================== #
# The route tables are one row per ROUTE, so a drug reached three ways
# appears three times and a reader cannot see that it was reached three
# ways. This inverts that: one row per compound, with the methods that
# found it counted.
#
# That count is the most useful number in the report. These methods share
# very little - a curated ligand annotation, a measured expression
# signature, a regulon statistic, a pocket shape, a chemical similarity -
# so they fail for unrelated reasons. A compound found by three is
# supported in a way no single method's score can express, and the route
# tables were hiding it.
#
# What it is NOT: a vote. Methods differ enormously in strength, so two
# weak agreements do not equal one measured result. The evidence tier
# still decides, and the method count breaks ties within a tier.

# Which methods are genuinely independent of each other. Two rows from
# sources sharing an input are not corroboration, so they collapse to one
# method here: from_grades and from_ligands both read the same ChEMBL
# annotation, and the two signature routes both read the same corpora.
# Keys are the `source` strings catalogue.py actually emits, read from it
# rather than guessed - the first version invented plausible names and
# matched none of them, so every row reported zero methods.
METHOD_OF_SOURCE = {
    "direct": "direct therapy for the gene",
    "drug_scoring": "annotated ligand (ChEMBL/DGIdb/Open Targets)",
    "ligands": "endogenous ligand availability",
    "signatures": "measured signature on your gene",
    "tfactivity": "regulon shift in perturbation corpora",
    "abundance": "measured change in the factor's own transcript",
    "screening": "ligand-based model (trained on this receptor)",
    "docking": "structure-based docking",
    "network": "upstream signalling network",
    "ctd": "literature record of an effect (CTD)",
    "regulators": "factor recorded as regulating the gene",
    "stabilizers": "protein stability / degradation",
    "catalogue": "chromatin-level agent",
}

# How much a method's agreement is worth, for ordering only. A measured
# result on your own gene is not the same kind of thing as a similarity
# search, and a flat count would say it was.
METHOD_WEIGHT = {
    "direct therapy for the gene": 1.0,
    "measured signature on your gene": 0.9,
    "annotated ligand (ChEMBL/DGIdb/Open Targets)": 0.7,
    "regulon shift in perturbation corpora": 0.6,
    "upstream signalling network": 0.5,
    "literature record of an effect (CTD)": 0.3,
    "factor recorded as regulating the gene": 0.4,
    "measured change in the factor's own transcript": 0.5,
    "endogenous ligand availability": 0.5,
    "ligand-based model (trained on this receptor)": 0.4,
    "protein stability / degradation": 0.4,
    "structure-based docking": 0.3,
    "chromatin-level agent": 0.1,
}


def unmapped_sources(result: dict) -> set:
    """Route sources with no method name, so a new source cannot be added
    without appearing in this table. Pinned by a test."""
    seen = set()
    for res in result.get("results", []):
        for r in (res.get("all_routes") or {}).get("routes") or []:
            src = r.get("source")
            if src and src not in METHOD_OF_SOURCE:
                seen.add(src)
    return seen


def _expected_change(routes: list, goal: str, outcomes: dict,
                     measured: dict) -> dict:
    """What should happen to the target gene, with every number it has.

    Three quantities are available and they answer different questions.

    `p_move` is the probability the gene moves at all, from outcome.py's
    four-step chain, with the interval it carries. It is a real computed
    number. It is also nearly constant across candidates, because two of
    the four steps are fixed - direction defaults to 0.75 without
    --predict-direction and responsiveness is pinned at 0.40 since nothing
    here can measure it - so it conveys how uncertain this CLASS of route
    is, not which candidate is better.

    `effect` is a measured signed effect size on your gene, present only
    where a corpus supplied one. A set-membership hit from a .gmt carries
    no magnitude at all and is reported as such rather than as its nominal
    placeholder score.

    What is deliberately absent is a predicted fold-change. Deriving one
    would need the share of the gene's transcription driven by that factor,
    the factor's occupancy change at the achieved dose, and the dose
    reached in your tissue. None of the three is available, so a number
    would be fabricated.
    """
    want = "increase" if goal == "up" else "decrease"
    drug = (routes[0].get("drug") or "").strip().lower()

    out = {"direction": want, "p_move": None, "p_low": None, "p_high": None,
           "effect": None, "effect_basis": None, "chain": None}

    o = outcomes.get(drug)
    if o:
        out["p_move"] = o.get("p_point")
        out["p_low"] = o.get("p_low")
        out["p_high"] = o.get("p_high")
        out["limiting_step"] = o.get("limiting_step")
        out["knockout"] = o.get("knockout")
        out["direction_measured"] = not any(
            "direction from database annotation only" in n
            for n in (o.get("notes") or []))
        out["chain"] = {
            "engages the target": o.get("p_engages"),
            "right direction": o.get("p_direction"),
            "acts on the promoter": o.get("p_activates"),
            "gene is responsive": o.get("p_responsive"),
        }

    m = measured.get(drug)
    if m is not None:
        out["effect"] = m["score"]
        out["effect_basis"] = m["basis"]

    if m is not None:
        out["basis"] = ("MEASURED: this compound moved your gene in a "
                        "cell-based assay, in the direction you want")
    elif any(r.get("source") == "signatures" for r in routes):
        out["basis"] = ("MEASURED: this compound moved your gene in a "
                        "cell-based assay, in the direction you want")
    elif any(r.get("source") == "direct" for r in routes):
        out["basis"] = ("established therapy acting on the gene or its "
                        "product")
    else:
        out["basis"] = ("PREDICTED from the factor's role and the drug's "
                        "action; no measurement of your gene supports it")
    return out


def _measured_effects(result: dict) -> dict:
    """{drug: signed effect on the gene} from whatever corpus supplied one.

    A `.gmt` hit is set membership and carries no magnitude - its score is
    a nominal placeholder - so it is reported as a direction without a
    number rather than as that placeholder.
    """
    from . import signatures as sig
    from .congeners import salt_free_name

    out: dict = {}
    for res in result.get("results", []):
        for h in ((res.get("lincs") or {}).get("hits") or []):
            name = salt_free_name((h.get("drug") or "").strip().lower())
            if not name:
                continue
            score = h.get("score")
            if score is None:
                continue
            nominal = abs(float(score)) == sig.GMT_SET_SCORE
            if nominal:
                out.setdefault(name, {
                    "score": None,
                    "basis": ("set membership in a consensus signature: the "
                              "compound moved this gene, with no magnitude "
                              "attached")})
            else:
                prev = out.get(name)
                if prev is None or prev["score"] is None or \
                        abs(score) > abs(prev["score"]):
                    out[name] = {
                        "score": round(float(score), 2),
                        "basis": ("measured signed effect on this gene in a "
                                  "perturbation corpus (z-like, not a "
                                  "fold-change)")}
    return out

def corroboration(result: dict, top_n: int = 20) -> list:
    """One row per drug across the top routes, with the methods that found it."""
    from .congeners import salt_free_name
    from . import indications as ind_mod, systemic as sys_mod

    goal = result["query"].get("goal", "up")
    routes = []
    outcomes: dict = {}
    for res in result["results"]:
        routes.extend((res.get("all_routes") or {}).get("routes") or [])
        for o in (res.get("outcome_estimates") or []):
            key = (o.get("drug") or "").strip().lower()
            # Keep the most favourable route's estimate for a drug reached
            # more than one way; the row reports the best available path.
            if key and (key not in outcomes
                        or (o.get("p_point") or 0)
                        > (outcomes[key].get("p_point") or 0)):
                outcomes[key] = o
    if not routes:
        return []
    measured = _measured_effects(result)

    # Properties for the BBB call, keyed salt-free so a salt form in the
    # route table finds the parent's structure.
    props: dict = {}
    for res in result["results"]:
        for gr in (res.get("drug_grades") or []):
            key = salt_free_name((gr.get("drug_name") or "").lower())
            pr = (gr.get("evidence") or {}).get("properties")
            if key and pr and key not in props:
                props[key] = pr

    # Rank by tier then score, and take the top routes to summarise.
    routes.sort(key=lambda r: (r.get("tier", 9), -(r.get("score") or 0)))

    grouped: dict = {}
    order: list = []
    for r in routes:
        stem = salt_free_name((r.get("drug") or "").lower())
        if not stem:
            continue
        if stem not in grouped:
            if len(order) >= top_n:
                continue          # the table is full; later routes only
                                  # add methods to drugs already in it
            order.append(stem)
            grouped[stem] = {"drug": r.get("drug"), "routes": []}
        grouped[stem]["routes"].append(r)

    out = []
    for stem in order:
        g = grouped[stem]
        rs = g["routes"]
        best = min(rs, key=lambda r: r.get("tier", 9))

        methods = []
        for r in rs:
            m = METHOD_OF_SOURCE.get(r.get("source"))
            if m and m not in methods:
                methods.append(m)

        desc = ind_mod.describe(
            g["drug"], mechanism=(best.get("mechanism") or "").strip() or None)
        cls = sys_mod.classify(g["drug"])
        _pr = props.get(salt_free_name((g["drug"] or "").lower()))
        bbb = dlv_mod.bbb_verdict(g["drug"], _pr)
        # Assay-interference structures undermine the EVIDENCE behind a
        # predicted route, which is what liabilities.py was written for
        # and was never connected to anything. It never demotes an
        # approved drug for being approved - only the confidence in a
        # prediction about it.
        liab = liab_mod.assess(g["drug"], (_pr or {}).get("canonical_smiles"))
        exp = _expected_change(rs, goal, outcomes, measured)
        from . import signatures as _sig
        _breadth = ((result.get("compound_breadth") or {})
                    .get("compounds") or {}).get(_sig.breadth_key(g["drug"]))

        out.append({
            "drug": g["drug"],
            "n_methods": len(methods),
            "methods": methods,
            "weight": round(sum(METHOD_WEIGHT.get(m, 0.2)
                                for m in methods), 2),
            "best_tier": best.get("tier"),
            "acts_on": sorted({r.get("acts_on") for r in rs
                               if r.get("acts_on")})[:3],
            "categories": sorted({r.get("category") for r in rs}),
            "mechanism": (best.get("mechanism") or "").strip(),
            "what_it_is": desc["what_it_is"],
            "treats": desc["treats"],
            "systems": (cls.name if cls else None),
            "systems_why": (cls.why if cls else None),
            "bbb": bbb["cell"],
            "bbb_why": bbb["why"],
            "route": desc.get("route"),
            "route_note": desc.get("route_note"),
            "generic_available": desc.get("generic_available"),
            "pains": bool(liab.pains),
            "evidence_multiplier": liab.evidence_multiplier,
            "expected_change": exp["direction"],
            "expected_basis": exp["basis"],
            "p_move": exp["p_move"],
            "p_low": exp["p_low"],
            "p_high": exp["p_high"],
            "effect_size": exp["effect"],
            "effect_basis": exp["effect_basis"],
            "chain": exp.get("chain"),
            "knockout": exp.get("knockout"),
            "direction_p": ((exp.get("chain") or {}).get("right direction")
                          if exp.get("direction_measured") else None),
            "limiting_step": exp.get("limiting_step"),
            "n_routes": len(rs),
            "breadth": _breadth,
            "mitotic_disruptor": sys_mod.is_mitotic_disruptor(g["drug"]),
        })

    out.sort(key=lambda r: (r["best_tier"], -r["weight"], -r["n_methods"]))
    return out


def _breadth_cell(b: dict | None) -> str:
    """How broadly a compound acts, for one table cell."""
    if not b:
        return "\u2014"
    if b.get("thin"):
        return "few signatures"
    pct = f"{b['percentile']:.0f}th pct, {b['median_genes']:.0f} genes"
    return f"\u26a0 broad ({pct})" if b.get("broad") else pct


def corroboration_table(result: dict, top_n: int = 20) -> list:
    """The cross-method table, as markdown lines."""
    rows = corroboration(result, top_n=top_n)
    if not rows:
        return []
    goal = result["query"].get("goal", "up")
    gene = result["query"]["gene"]

    L = [f"## All {len(rows)} candidates, and which methods found each", ""]
    L.append("One row per compound rather than per route, so a drug reached "
             "several ways is visible as one candidate with several "
             "supports.")
    L.append("")
    L.append(f"| Drug | Methods | Tier | Acts on | What it is | Systems it "
             f"acts on | Breadth | Reach | BBB? | {gene} change | "
             f"Engagement / promoter | "
             + ("Direction | " if any(r.get("direction_p") is not None
                                      for r in rows) else "")
             + ("Knockout evidence | " if any(r.get("knockout") is not None
                                              for r in rows) else "")
             + "Measured effect |")
    any_ko = any(r.get("knockout") is not None for r in rows)
    any_dir = any(r.get("direction_p") is not None for r in rows)
    L.append("|" + "---|" * (12 + any_ko + any_dir))
    any_p, any_effect = False, False
    for r in rows:
        sys_cell = r["systems"] or "no systemic class flagged"
        if r.get("mitotic_disruptor"):
            sys_cell += " + halts mitosis"
        what = r["what_it_is"] or r["mechanism"] or "\u2014"
        tag = ("measured" if "MEASURED" in r["expected_basis"]
               else "established" if "established" in r["expected_basis"]
               else "predicted")

        chain = r.get("chain") or {}
        eng, act = chain.get("engages the target"), chain.get(
            "acts on the promoter")
        if eng is None or act is None:
            p_cell = "not estimated"
        else:
            any_p = True
            p_cell = f"{eng:.2f} / {act:.2f}"

        eff = r.get("effect_size")
        if eff is not None:
            any_effect = True
            eff_cell = f"{eff:+.2f}"
        elif r.get("effect_basis"):
            eff_cell = "direction only"
        else:
            eff_cell = "\u2014"

        L.append(f"| **{r['drug']}** | {r['n_methods']} | "
                 f"{r['best_tier']} | "
                 f"{', '.join(r['acts_on']) or EMDASH} | "
                 f"{_clip(what, 34)} | {_clip(sys_cell, 28)} | "
                 f"{_breadth_cell(r.get('breadth'))} | "
                 f"{_route_cell(r)} | "
                 f"{r.get('bbb') or '—'} | "
                 f"{r['expected_change']} ({tag}) | {p_cell} | "
                 + (f"{r['direction_p']:.2f} | " if any_dir and r.get(
                     "direction_p") is not None else
                    "— | " if any_dir else "")
                 + (f"{ko_mod.cell(r.get('knockout'))} | " if any_ko else "")
                 + f"{eff_cell} |")
    L.append("")

    info = result.get("compound_breadth") or {}
    if any(r.get("breadth") for r in rows):
        L.append(
            f"**Breadth** is how many genes the compound moved across its "
            f"LINCS signatures (up and down together), against all "
            f"{info.get('n_signatures', 0):,} signatures. \u26a0 broad means "
            f"the top quartile. It is a caution, not a verdict: it rises "
            f"with dose and depends on the cell line, and the files hold "
            f"set membership without scores, so they cannot say whether "
            f"{gene} stands out among the genes a compound changed. A "
            f"compound moving hundreds of genes lowers {gene} among many "
            f"others. Compounds with fewer than "
            f"{info.get('min_signatures', 3)} signatures are not judged.")
        L.append("")

    # What the two numeric columns are, and what neither of them is.
    L.append(f"**There is no predicted fold-change here, and that is not an "
             f"omission.** Deriving one would need the share of {gene}'s "
             f"transcription driven by that factor, the factor's occupancy "
             f"change at the dose achieved, and the dose reached in your "
             f"tissue. None of the three is available to this package, so "
             f"any number in a fold-change column would be fabricated. The "
             f"two numbers below are the ones it can actually compute.")
    L.append("")
    if any_p:
        L.append("**Engagement / promoter** are the two links that come from "
                 "evidence: the drug engaging its target at a dose that can "
                 "be reached, and the factor acting on this promoter. Each "
                 "is a 0 to 1 heuristic and they differ from row to row. "
                 "**Direction** and **Knockout evidence** appear only when "
                 "a model predicted the direction or an experiment "
                 "perturbed the factor. Order on tier and method count.")
        L.append("")
    if any_effect:
        L.append("**Measured effect** is a signed effect size on your gene "
                 "from a perturbation corpus, z-like rather than a "
                 "fold-change. `direction only` means the corpus recorded "
                 "that the compound moved the gene without a magnitude, "
                 "which is what a consensus gene-set membership carries.")
        L.append("")

    multi = [r for r in rows if r["n_methods"] > 1]
    if multi:
        L.append(f"**{len(multi)} of {len(rows)} were found by more than one "
                 f"method.** These methods share little - a curated ligand "
                 f"annotation, a measured expression signature, a regulon "
                 f"statistic, a pocket shape - so they fail for unrelated "
                 f"reasons, and agreement between them is worth more than a "
                 f"better score from any one.")
        L.append("")
        for r in multi[:8]:
            L.append(f"- **{r['drug']}**: {'; '.join(r['methods'])}")
        L.append("")
    else:
        L.append(f"**No compound was found by more than one method.** With "
                 f"{len({m for r in rows for m in r['methods']})} method(s) "
                 f"contributing, that is expected rather than informative: "
                 f"corroboration needs two methods to have run and both to "
                 f"reach the same compound. `report.md` lists the flags that "
                 f"turn the others on.")
        L.append("")

    L.append("**The method count is not a vote.** These methods differ "
             "enormously in strength, so two weak agreements do not equal "
             "one measured result. Tier still decides the order and the "
             "count breaks ties within a tier.")
    L.append("")
    bbb_counts: dict = {}
    for r in rows:
        c = r.get("bbb") or "—"
        bbb_counts[c] = bbb_counts.get(c, 0) + 1
    L.extend(_bbb_footnote(bbb_counts))

    local = [r for r in rows if r.get("route") == "local only"]
    if local:
        n_l = len(local)
        L.append(f"**Reach** is the FDA route of administration, from the "
                 f"Orange Book. **{n_l} "
                 f"{'candidate here is' if n_l == 1 else 'candidates here are'}"
                 f" approved only by local routes** "
                 f"({', '.join(r['drug'] for r in local[:6])}) - topical, "
                 f"ophthalmic or inhaled. Those cannot reach a target "
                 f"outside the tissue they are applied to, however they "
                 f"score above, so for an internal gene they are not "
                 f"usable without reformulation.")
        L.append("")
    elif any(r.get("route") for r in rows):
        L.append("**Reach** is the FDA route of administration. Every "
                 "candidate here has a systemic route, so none is ruled "
                 "out for being topical-only.")
        L.append("")
    else:
        L.append("**Reach** is blank for every row because no Orange Book "
                 "is on disk. That is unknown route, not systemic - a "
                 "topical-only drug looks identical. See "
                 "`orangebook.SOURCE` for the download.")
        L.append("")

    pains = [r for r in rows if r.get("pains")]
    if pains:
        n_p = len(pains)
        L.append(f"**{n_p} "
                 f"{'candidate matches' if n_p == 1 else 'candidates match'} "
                 f"a pan-assay interference pattern** "
                 f"({', '.join(r['drug'] for r in pains[:5])}). They are "
                 f"approved drugs and this says nothing about their "
                 f"safety. What it undermines is the EVIDENCE: these "
                 f"structures read active across unrelated assays through "
                 f"aggregation, redox cycling or fluorescence "
                 f"interference, so a measured or predicted activity for "
                 f"them is likelier to be an artefact than the tier "
                 f"suggests.")
        L.append("")

    L.append(f"**The {gene} change column is a direction.** `measured` "
             f"means a compound moved {gene} in a cell-based assay; "
             f"`established` means it is an approved therapy acting on the "
             f"gene or its product; `predicted` means it follows from the "
             f"factor's role and the drug's action, with no measurement of "
             f"{gene} behind it.")
    L.append("")
    return L


def summary_report(result: dict, top_n: int = SUMMARY_TOP) -> str:
    """A one-page summary. Everything here is also in report.md."""
    result = escape_pipes(fda_mod.apply(result))
    q = result["query"]
    goal = q.get("goal", "up")
    verb = "raise" if goal == "up" else "lower"
    L: list[str] = []

    L.append(f"# {q['gene']}: how to {verb} it - summary")
    L.append("")
    L.append(f"Objective: **{verb} {q['gene']}**. "
             f"{'Demo data - not a real analysis. ' if q.get('demo') else ''}"
             f"Full detail, with every qualification, in `report.md`.")
    L.append("")
    L.extend(fda_mod.notice(result))

    soc = result.get("standard_of_care") or {}
    if soc.get("treatments"):
        L.append(f"**Already approved for {soc['condition']}:** "
                 + ", ".join(t["generic"] for t in soc["treatments"])
                 + f". {soc['note']}")
        L.append("")

    verdict = result.get("evidence_verdict") or {}
    L.append(_verdict_line(verdict))
    for r in (verdict.get("reasons") or [])[:3]:
        L.append(f"- {r}")
    # Only when the verdict's own reasons have not already said it. The
    # first draft printed the same 78% twice, two lines apart.
    fdp = verdict.get("false_discovery_proportion")
    said = any("chance" in r or "expected to be" in r
               for r in (verdict.get("reasons") or []))
    if fdp is not None and not said:
        n_sites = verdict.get("scan_sites")
        rng = verdict.get("false_discovery_range")
        of = f"the {n_sites} motif sites" if n_sites else "the motif sites"
        span = (f" (range {rng[0]:.0%} to {rng[1]:.0%})"
                if rng else "")
        L.append(f"- About **{fdp:.0%} of {of}**{span} would be expected on "
                 f"shuffled sequence of the same window. A site is a "
                 f"sequence match, not evidence of binding; ChIP peaks are "
                 f"what separate the real ones.")
    L.append("")

    # ---- the shortlist
    preds = []
    for res in result["results"]:
        preds.extend((res.get("all_routes") or {}).get("routes") or [])
    # One entry per drug, led by its best route, with every target it acts
    # through. The same drug reached through two factors is one candidate.
    from . import clusters as cl_mod
    from . import signatures as _sg
    by_drug: dict = {}
    for r in sorted(preds, key=lambda x: (x.get("tier", 9),
                                          -(x.get("score") or 0))):
        dk = _sg.breadth_key(r.get("drug", "")) or r.get("drug", "").lower()
        if dk in by_drug:
            if r.get("acts_on"):
                by_drug[dk]["targets"].add(r["acts_on"])
            by_drug[dk]["families"] |= cl_mod.family_keys([r])
            continue
        if len(by_drug) >= top_n * SUMMARY_POOL:
            continue
        by_drug[dk] = dict(r, name=r.get("drug", ""),
                           targets={r["acts_on"]} if r.get("acts_on")
                           else set(),
                           families=cl_mod.family_keys([r]))
    # Drugs of one mechanism family (the same factor, reached the same way)
    # are one idea to test first, so each family takes one row, led by its
    # best-ranked member. Close analogues acting through a shared target fold
    # in the same way. Every drug stays in the full table in report.md.
    clusters = cl_mod.group(
        list(by_drug.values()),
        same_route=lambda a, b: bool(a["targets"] & b["targets"]),
        same_family=lambda a, b: bool(a["families"] & b["families"]))[:top_n]
    rows = [c["leader"] for c in clusters]
    folded = {c["leader"]["drug"]: c for c in clusters if c["members"]}

    if rows:
        from . import indications as ind_mod
        from . import rationale as rat_mod
        _gene = q["gene"]
        _roles = {t["symbol"]: t.get("role")
                  for res in result["results"] for t in res.get("tfs", [])}
        from . import strata as st_mod
        _all_routes = [rt for res in result["results"]
                       for rt in ((res.get("all_routes") or {}).get("routes")
                                  or [])]
        _groups = st_mod.partition(rows, _all_routes)
        L.append("## What to test first")
        L.append("")
        L.append(st_mod.preamble(len(_groups)))
        L.append("")
        L.append("Drugs that reach the gene through the same factor by the "
                 "same kind of route, and close analogues acting through a "
                 "shared target, are folded into one row, so a family is "
                 "one idea to test rather than several.")
        L.append("")
        L.append(rat_mod.aim(_gene, goal, result.get("requirement"),
                             q.get("condition"), q.get("scenario")))
        L.append("")
        _line: dict = {}
        _any_broad = False
        _any_mitotic = False
        from . import crosscheck as xc_mod
        _xc = next((b.get("crosscheck") for b in result["results"]
                    if b.get("crosscheck")), None)
        _ctd_opposite: list = []
        for r in rows:
            d = ind_mod.describe(
                r.get("drug", ""),
                mechanism=(r.get("mechanism") or "").strip() or None)
            treats = d["treats"] or "not established here"
            from . import signatures as _sig
            _b = ((result.get("compound_breadth") or {}).get("compounds")
                  or {}).get(_sig.breadth_key(r.get("drug", "")))
            _broad = " ⚠broad" if (_b or {}).get("broad") else ""
            _any_broad = _any_broad or bool(_broad)
            if (result.get("mitotic_gene")
                    and sys_mod.is_mitotic_disruptor(r.get("drug", ""))):
                _broad += " ⚠halts mitosis"
                _any_mitotic = True
            _fold = folded.get(r.get("drug"))
            _more = f" (+{len(_fold['members'])} similar)" if _fold else ""
            for _m in [r] + (_fold["members"] if _fold else []):
                if xc_mod.cautions_for(_m, _xc):
                    _ctd_opposite.append(_m.get("drug"))
            if xc_mod.cautions_for(r, _xc):
                _broad += " ⚠CTD opposite"
            _targets = ", ".join(sorted(r["targets"])) or EMDASH
            _line[id(r)] = (
                f"| **{r.get('drug')}**{_broad}{_more} | "
                f"{_clip(treats)} | {_targets} | "
                f"{r.get('tier')} | "
                f"{rat_mod.route_sentence(r, _gene, goal, _roles)} |")

        for _key, _title, _blurb, _grp in _groups:
            L.append(f"### {_title} ({len(_grp)})")
            L.append("")
            L.append(_blurb)
            L.append("")
            L.append("| Drug | Approved to treat | Acts on | Tier | How it "
                     "would work |")
            L.append("|---|---|---|---|---|")
            for r in _grp:
                L.append(_line[id(r)])
            L.append("")

        _excluded = st_mod.exclusions(_all_routes)
        if _excluded:
            L.append(f"### Found and excluded ({len(_excluded)})")
            L.append("")
            L.append("These were reached by the search and then rejected. "
                     "They are shown because a reader who expected one of "
                     "them should learn it was considered, not be left "
                     "wondering whether the search missed it.")
            L.append("")
            L.append("| Compound | Why it is not a candidate |")
            L.append("|---|---|")
            for _e in _excluded[:20]:
                L.append(f"| {_e['drug']} | {_clip(_e['why'], 150)} |")
            if len(_excluded) > 20:
                L.append(f"| … | {len(_excluded) - 20} more, in "
                         f"`report.md` |")
            L.append("")
        if folded:
            L.append("**Folded into the rows above** (the same factor "
                     "reached the same way, or structure within Tanimoto "
                     "0.75 of the leader with a shared target). They are "
                     "not duplicates: members of a family can differ in "
                     "potency, selectivity, tissue reach and safety, and "
                     "each has its own row in `report.md`.")
            L.append("")
            for r in rows:
                c = folded.get(r.get("drug"))
                if c:
                    L.append(f"- {r.get('drug')}: "
                             f"{cl_mod.names(c, limit=8)}")
            L.append("")
        if _ctd_opposite:
            L.append(f"⚠CTD opposite: CTD records "
                     f"{', '.join(dict.fromkeys(_ctd_opposite))} moving "
                     f"{q['gene']} the other way. A route with no recorded "
                     f"direction at its factor is lowered to tier 4; the "
                     f"others are only flagged, since CTD mixes cell types "
                     f"and doses. `report.md` lists them all.")
            L.append("")
        if _xc and _xc.get("conflicts"):
            _cn = [c["drug"] for c in _xc["conflicts"]]
            L.append(f"**CTD and LINCS disagree on {len(_cn)} drugs** (they "
                     f"agree on {_xc['n_agree']} of {_xc['n_compared']} "
                     f"they both call): {', '.join(_cn[:8])}"
                     f"{' and more' if len(_cn) > 8 else ''}.")
            L.append("")
        if _any_mitotic:
            L.append("⚠halts mitosis: the compound is suspected to hold "
                     "cells in mitosis instead of letting them proceed, and "
                     "this gene is a mitotic one, so it moves the gene by "
                     "changing the cell's cycle. Its route is scored at "
                     "half.")
            L.append("")
        if _any_broad:
            L.append("⚠broad: the compound moves a large share of the "
                     "transcriptome (top quartile of LINCS signatures), so "
                     "lowering this gene is one of many changes it makes. "
                     "A caution about specificity, not a score; the full "
                     "table in `report.md` says how it is measured.")
            L.append("")
        L.append("**Every row is off-label.** The third column is what the "
                 "drug was approved for; using it to move your gene is a "
                 "different purpose, and the distance between the two is "
                 "part of what you are judging.")
        L.append("")
    else:
        L.append("## No route was found")
        L.append("")
        L.append("Nothing cleared the thresholds. `report.md` has the "
                 "per-stage accounting for where candidates were lost.")
        L.append("")

    L.extend(corroboration_table(result, top_n=20))

    # ---- what ran and what did not: the thing that decides how to read it
    sources = ((result.get("evidence_sources") or {}).get("sources") or {})
    if sources:
        # Four states, and they are four different pieces of news. The
        # one that must never be merged upward is `empty`: a source that
        # ran and measured nothing is a RESULT, and reporting it beside
        # one that broke throws away the distinction the ledger is for.
        _by = {}
        for k, v in sources.items():
            _by.setdefault(v.get("state"), []).append(k)
        ok = sorted(_by.get("ok") or [])
        empty = sorted(_by.get("empty") or [])
        failed = sorted(_by.get("FAILED") or [])
        absent = sorted(_by.get("not requested") or [])
        if ok or empty or failed or absent:
            L.append("## What the answer rests on")
            L.append("")
            if ok:
                L.append(f"- **Ran, and found something:** "
                         f"{', '.join(ok)}")
            if empty:
                L.append(f"- **Ran, and found nothing:** "
                         f"{', '.join(empty)}. That is a measurement, not "
                         f"a gap: these were asked and answered no.")
            if failed:
                L.append(f"- **Could not run:** {', '.join(failed)}. A "
                         f"source that broke is not a negative result, so "
                         f"nothing below should be read as evidence that "
                         f"it would have found nothing.")
            if absent:
                L.append(f"- **Not available here:** {', '.join(absent)}. "
                         f"No input for these on this machine, so they "
                         f"were never asked. `report.md` names the file or "
                         f"flag each one needs.")
            L.append("")
        elif not q.get("demo"):
            L.append("## What the answer rests on")
            L.append("")
            L.append("No optional evidence source was requested, so this "
                     "ranking is motif prediction and drug annotation "
                     "only. `report.md` lists the flags that add measured "
                     "evidence.")
            L.append("")

    # ---- the things most likely to mislead, from the block itself
    cautions = []
    for res in result["results"]:
        at = res.get("attrition") or {}
        if at.get("n_drugged_nuclear_receptor"):
            cautions.append(
                f"The drug tables are dominated by nuclear receptors "
                f"because {at.get('n_with_drug', 0)} of "
                f"{at.get('n_tfs_with_sites', 0)} factors with a site here "
                f"have any ligand at all. That is attrition, not a "
                f"ranking.")
        if at.get("n_top25_without_drug"):
            names = ", ".join((at.get("top25_without_drug") or [])[:4])
            cautions.append(
                f"{at['n_top25_without_drug']} of the best-scoring factors "
                f"have no approved ligand ({names}) and appear in no drug "
                f"table. The gene-specific signal is partly in them.")
        act = res.get("tf_activity") or {}
        if act.get("n_testable") and not [
                a for a in (act.get("scanned") or [])
                if _sig_q(a.get("q_value")) < 0.05]:
            cautions.append(
                "No factor with a site in this promoter has a detectably "
                "perturbed regulon in the tissue checked, so these routes "
                "rest on sequence alone.")
        ov = result.get("tss_window_overlap") or []
        if ov:
            cautions.append(
                "The promoter windows overlap each other, so counts across "
                "them are not independent.")
        break

    if cautions:
        L.append("## Read these before the table")
        L.append("")
        for c in cautions[:4]:
            L.append(f"- {c}")
        L.append("")

    # ---- what to do next
    L.append("## Next step")
    L.append("")
    L.append(f"Every row above is a hypothesis. The experiment that settles "
             f"one is the same in each case: a promoter-reporter assay for "
             f"the element, a dose-response with RT-qPCR or protein readout "
             f"for the drug, and a site-directed mutant of the element to "
             f"show the drug effect runs through it. Without the third, a "
             f"change in {q['gene']} does not establish the mechanism this "
             f"report proposes.")
    L.append("")
    L.append("*Not medical advice. Hypotheses for laboratory testing.*")
    L.append("")
    return "\n".join(L)


def write_outputs(result: dict, outdir: str | Path, *, top_n: int = 25) -> list[Path]:
    """Write every output; top_n limits report tables, not exported data."""
    if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
        raise ValueError("top_n must be a positive integer")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []

    p = outdir / "report.md"
    p.write_text(markdown_report(result, top_n=top_n), encoding="utf-8")
    written.append(p)

    # The short version. report.md is long because every number in it is
    # qualified; this is the page to open first.
    p = outdir / "summary.md"
    p.write_text(summary_report(result), encoding="utf-8")
    written.append(p)

    p = outdir / "results.json"
    p.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    written.append(p)

    p = outdir / "tf_summary.tsv"
    cols = ["symbol", "name", "role", "nuclear_receptor", "n_sites", "best_rel_score",
            "best_score_bits", "best_site_tss_offset", "best_site_sequence",
            "matrix_id", "chip_peaks", "score", "n_drugs", "drugs"]
    # csv.writer rather than "\t".join: names come from external tables
    # and APIs, and one carrying a tab or a newline used to split the row
    # and shift every column after it, silently. QUOTE_MINIMAL quotes
    # only the fields that need it, so an ordinary file is unchanged.
    import csv
    import io

    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n",
                   quoting=csv.QUOTE_MINIMAL)
    w.writerow(cols)
    for res in fda_mod.apply(result)["results"]:
        for t in res["tfs"]:
            names = "; ".join(d["name"] for d in t["drugs"])
            row = [t.get(c) for c in cols[:-2]] + [len(t["drugs"]), names]
            w.writerow(["" if v is None else str(v) for v in row])
    p.write_text(buf.getvalue(), encoding="utf-8")
    written.append(p)

    # BED of predicted sites, for genome-browser inspection.
    p = outdir / "sites.bed"
    lines = ['track name="predicted_TFBS" description="PWM hits" useScore=1']
    for res in result["results"]:
        prom = res["promoter"]
        for h in res["hits"]:
            strand = h["strand"]
            if prom["strand"] == 1:
                gs = prom["start"] - 1 + h["start"]
                ge = prom["start"] - 1 + h["end"]
            else:
                gs = prom["end"] - h["end"]
                ge = prom["end"] - h["start"]
                strand = {"+": "-", "-": "+"}.get(strand, strand)
            score = max(0, min(1000, int(h["rel_score"] * 1000)))
            # BED forbids whitespace in the name field and cannot be
            # quoted, so this is sanitised rather than escaped - lossy
            # by necessity. A motif name with a space in it would
            # otherwise be read as the start of the score column.
            name = f"{h['motif_name']}|{h['matrix_id']}"
            name = "_".join(str(name).split()) or "site"
            lines.append("\t".join([
                f"chr{prom['chrom']}", str(gs), str(ge),
                name, str(score), strand,
            ]))
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    written.append(p)

    p = outdir / "promoter.fa"
    fa = []
    for res in result["results"]:
        prom = res["promoter"]
        fa.append(f">{prom['gene_symbol']}|{prom['transcript_id']}|"
                  f"chr{prom['chrom']}:{prom['start']}-{prom['end']}|"
                  f"strand={prom['strand']}|-{prom['upstream']}/+{prom['downstream']}")
        seq = res["sequence"]
        fa += [seq[i:i + 60] for i in range(0, len(seq), 60)]
    p.write_text("\n".join(fa) + "\n", encoding="utf-8")
    written.append(p)

    return written
