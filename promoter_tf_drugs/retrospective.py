"""Does the pipeline recover drugs already known to move a gene?

Why this exists
---------------
A ranked list says nothing about whether the ranking works. This holds the
pipeline to the one question that can be asked of saved results without new
experiments: for a gene where earlier work already records which chemicals
raise or lower it, do the routes the pipeline predicted point at those
chemicals, and at the right direction?

Two answer keys
---------------
CTD. The Comparative Toxicogenomics Database curates chemical effects on gene
expression from papers. It is read from the saved result
(`improvements.reported_links`), so nothing is fetched. It plays no part in
route scores. It is wide and noisy: toxicology doses, cell lines and
exposures count the same as clinical findings, so most of what it records for
a gene is a toxicant, not a drug anyone would repurpose.

LINCS. The local L1000 signature files record which compounds moved the gene
in a dish. A drug counts as wanted when it moved the gene the wanted way in at
least two signatures and at least twice as often as the other way, as the
opposite when the reverse holds, and is set aside when the signs mix or it
appears once. Only drugs the files contain are compared, because a drug never
tested in LINCS has not been seen not moving the gene. Its weakness is that it
is a measurement in a dish at a high dose, not a finding about a patient.

The two keys fail in different ways, so agreement between them and the
pipeline carries more weight than either alone.

How the routes are split
------------------------
  predicted        routes through a transcription factor, an upstream
                   signalling node, a chromatin class, or promoter DNA. The
                   pipeline reasons its way to these. This is the group the
                   test is for.
  predicted, no    the same without the two TF categories built from LINCS
  LINCS            (a factor whose own transcript moves, or whose activity is
                   inferred from signatures). This is the only predicted
                   group scored against the LINCS key, because the others
                   would be graded on their own input.
  measured         routes from LINCS signatures: the drug was seen moving the
                   gene. Scored against CTD only, since against LINCS it
                   would be the key itself.

Left out, and why: routes that bind the gene's own product (they say nothing
about expression) and virtual screen hypotheses (CTD marks which screen hits
are already known, which shapes their labels).

Two tests per key, on the approved library:

  reported or not    among drugs in scope, do those the key records as moving
                     the gene in the wanted direction score higher than drugs
                     it has no record of? A drug with no route scores zero, so
                     this rewards reaching them at all as well as ranking them.
  right direction    among drugs the key records as moving the gene at all, do
                     those moving it the wanted way score higher than those
                     moving it the other way? This removes the effect of a
                     drug being well studied or broadly active, because both
                     groups were seen moving the gene.

Each is an area under the ROC curve: 0.5 is chance, 1.0 is perfect. A
permutation test says how often shuffled labels do as well. The first test
also reports a baseline that knows nothing about the gene, the number of routes
a drug has, so a pipeline that only favours drugs that appear everywhere does
not look good. Beside each score is how many of the positives have any route
at all, because a low AUC from a pipeline that never reaches the positives is
a coverage problem, not a ranking problem.

What it cannot show
-------------------
**The keys are noisy.** A perfect predictor would not score 1.0 against
either, and a drug missing from CTD is unstudied, not inert.

**It tests expression, not benefit.** Recovering drugs that move a gene says
the routes find real regulators. It does not say moving the gene treats a
disease.

**It tests the past.** Drugs the keys already record are the ones earlier work
could find. The value of the pipeline is in the ones it does not, which this
cannot score.

**A handful of genes is a handful of genes.** The mean over runs carries wide
uncertainty. The count of runs above chance is shown beside it so one strong
gene cannot carry the number, and a run with fewer than five positives or five
negatives is not scored.

**Independence is not total.** The routes that remain in "predicted, no LINCS"
come from target and network tables, not from signatures. Their scores may
still include terms (a direction check, a potency) that were fitted using
signature data somewhere upstream. That is not traced here.

**Saved results may predate the current code.** The report lists each file.

Drug fame
---------
Both keys record drugs the literature studied. A drug written up for many
genes is more likely to be recorded for this one too, whatever it does, so any
score that merely prefers well-known drugs looks good against them. When five
or more genes are scored together, the report adds a fame baseline for CTD: how
many of the OTHER genes record the drug. It also repeats the reported-or-not
test on only the drugs recorded for three or more other genes, where fame no
longer separates the groups. The second number is the one to trust: a score
that falls to 0.5 there was ranking on how well known a drug is.

Positive controls
-----------------
The two keys above say whether the routes agree with records. A different
question is whether the pipeline, given a gene whose regulation by an approved
drug is textbook pharmacology, puts that drug in front of the reader. A short
table of such genes (`POSITIVE_CONTROLS`) is checked against any run of those
genes: where the expected drugs rank among the run's candidates, by which
route, or that they are absent. The genes were picked for having a drug,
a receptor and a gene that every pharmacology text connects. The table is the
author's judgement of established pharmacology, written down before the runs
were made, and it carries no citations: read it critically and change it.
Passing a control shows the pipeline can find a known answer. Failing one is
the more informative result.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .signatures import breadth_key

log = logging.getLogger(__name__)

PREDICTED_PREFIXES = ("TF:", "upstream signalling", "chromatin-level",
                      "promoter DNA binder")
LINCS_DERIVED = ("TF: abundance change", "TF: inferred activity change")
MEASURED_CATEGORIES = ("empirical signature",)

GROUPS = {
    "predicted": "routes the pipeline reasoned its way to",
    "predicted_no_lincs": "the same without the two TF categories built "
                          "from LINCS",
    "measured": "LINCS signature routes (measured, not predicted)",
}
# Which route groups each key may grade. A group built from LINCS is never
# graded against LINCS.
KEY_GROUPS = {"ctd": ("predicted", "predicted_no_lincs", "measured"),
              "lincs": ("predicted_no_lincs",)}
KEY_LABELS = {"ctd": "CTD (curated literature)",
              "lincs": "LINCS (measured in a dish)"}

# Gene -> (goal, drug names accepted as the known answer, the mechanism). Every
# entry is a drug that acts on a receptor or factor that drives the gene, so a
# route through that factor is the expected finding. Names are matched after
# salt words and spelling are normalised.
POSITIVE_CONTROLS = {
    "CYP3A4": ("up", ("rifampin", "rifampicin"),
               "rifampicin activates PXR (NR1I2), which induces CYP3A4"),
    "CYP24A1": ("up", ("calcitriol", "paricalcitol", "calcipotriene",
                       "calcipotriol", "doxercalciferol", "alfacalcidol"),
                "vitamin D receptor ligands induce their own catabolic "
                "enzyme"),
    "FKBP5": ("up", ("dexamethasone", "prednisolone", "prednisone",
                     "hydrocortisone", "methylprednisolone", "cortisol",
                     "betamethasone", "triamcinolone", "fluticasone"),
              "glucocorticoids activate the glucocorticoid receptor "
              "(NR3C1), which induces FKBP5"),
    "ANGPTL4": ("up", ("fenofibrate", "bezafibrate", "pioglitazone",
                       "rosiglitazone", "gemfibrozil"),
                "PPAR agonists induce ANGPTL4"),
    "LDLR": ("up", ("atorvastatin", "simvastatin", "rosuvastatin",
                    "lovastatin", "pravastatin", "fluvastatin",
                    "pitavastatin"),
             "statins lower cell cholesterol, which activates SREBP2, which "
             "induces LDLR"),
    "CYP1A1": ("up", ("omeprazole", "lansoprazole", "leflunomide"),
               "these compounds activate the aryl hydrocarbon receptor, "
               "which induces CYP1A1"),
}

MIN_PER_CLASS = 5       # positives and negatives needed to score a run
MIN_SIGNATURES = 2      # LINCS signatures needed to call a direction
LINCS_MARGIN = 2.0      # times more often one way than the other
TOP_K = 20
FAME_MIN_RUNS = 5       # genes needed before fame means anything
FAME_FLOOR = 3          # a drug recorded for this many OTHER genes is "known"


@dataclass
class KeySet:
    """What an answer key says about one gene, as normalised drug names."""
    name: str
    wanted: set = field(default_factory=set)
    opposite: set = field(default_factory=set)
    unclear: set = field(default_factory=set)
    n_records: int = 0
    # Drugs the key could have recorded. None means every approved drug. For
    # LINCS it is the drugs the files contain: a drug never tested there has
    # not been seen not moving the gene.
    universe: set | None = None


@dataclass
class Run:
    gene: str
    goal: str
    path: str
    routes: list
    keys: dict                       # name -> KeySet
    skipped: dict                    # key name -> why it could not be built


# ------------------------------------------------------------------ loading
def _keys(links: list) -> set:
    return {breadth_key(x.get("chemical") or "") for x in links or []} - {""}


def ctd_keyset(blocks: list) -> tuple:
    """`(KeySet, None)` from the saved CTD links, or `(None, reason)`."""
    imp = {}
    for blk in blocks:
        cand = blk.get("improvements") or {}
        if cand.get("reported_links"):
            imp = cand
            break
    rl = imp.get("reported_links")
    if not rl:
        return None, "the saved result has no CTD record"
    # Older results lack the flag. The pipeline itself reads a missing flag as
    # assessed, so a present `reported_links` does too.
    if not imp.get("ctd_assessed", True):
        return None, "CTD was not read for this run"
    n = int(imp.get("ctd_links_found") or 0)
    if not n:
        return None, "CTD holds no expression record for this gene"
    return KeySet("ctd", wanted=_keys(rl.get("wanted_direction")),
                  opposite=_keys(rl.get("opposite_direction")),
                  unclear=_keys(rl.get("unclear_direction")),
                  n_records=n), None


def lincs_keyset(gene: str, goal: str, up_path=None, down_path=None) -> tuple:
    """`(KeySet, None)` from the local LINCS files, or `(None, reason)`."""
    from . import signatures as sg
    if goal not in ("up", "down"):
        return None, f"goal {goal!r} is not up or down"
    other = "down" if goal == "up" else "up"
    try:
        same = sg.load_enrichr_gmt(gene, direction=goal, up_path=up_path,
                                   down_path=down_path)
        opp = sg.load_enrichr_gmt(gene, direction=other, up_path=up_path,
                                  down_path=down_path)
        tested = sg.signature_breadth(up_path, down_path)
    except FileNotFoundError:
        return None, "no local LINCS .gmt files"
    if tested is None:
        return None, "no local LINCS .gmt files"
    n_w: dict = {}
    n_o: dict = {}
    for h in same:
        k = breadth_key(h.drug)
        n_w[k] = n_w.get(k, 0) + 1
    for h in opp:
        k = breadth_key(h.drug)
        n_o[k] = n_o.get(k, 0) + 1
    ks = KeySet("lincs", universe=set(tested["compounds"]))
    for k in (set(n_w) | set(n_o)) - {""}:
        w, o = n_w.get(k, 0), n_o.get(k, 0)
        if w >= MIN_SIGNATURES and w >= LINCS_MARGIN * o:
            ks.wanted.add(k)
        elif o >= MIN_SIGNATURES and o >= LINCS_MARGIN * w:
            ks.opposite.add(k)
        else:
            ks.unclear.add(k)
    ks.n_records = len(ks.wanted | ks.opposite | ks.unclear)
    if not ks.n_records:
        return None, f"no LINCS signature contains {gene}"
    return ks, None


def load_run(path, *, gmt_up=None, gmt_down=None) -> Run | None:
    """One saved result as a `Run`, or None when it cannot be read."""
    path = Path(path)
    if path.is_dir():
        path = path / "results.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        log.warning("%s: unreadable (%s)", path, str(exc)[:60])
        return None
    blocks = data.get("results") or []
    if not blocks:
        return None
    routes = []
    for blk in blocks:
        routes.extend((blk.get("all_routes") or {}).get("routes") or [])
    gene = (data.get("query") or {}).get("gene", "?")
    goal = (data.get("query") or {}).get("goal", "?")
    keys, skipped = {}, {}
    for name, got in (("ctd", ctd_keyset(blocks)),
                      ("lincs", lincs_keyset(gene, goal, gmt_up, gmt_down))):
        ks, why = got
        if ks is None:
            skipped[name] = why
        else:
            keys[name] = ks
    return Run(gene=gene.upper(), goal=goal, path=str(path), routes=routes,
               keys=keys, skipped=skipped)


def latest_per_gene_goal(runs: list) -> tuple:
    """One run per (gene, goal): the most recently written file.

    A folder often holds several reruns of the same gene. Counting each as an
    independent gene would multiply its weight, so the newest is kept and the
    rest are named.
    """
    best: dict = {}
    dropped = []
    for r in runs:
        k = (r.gene.upper(), r.goal)
        if k not in best:
            best[k] = r
        elif Path(r.path).stat().st_mtime > Path(best[k].path).stat().st_mtime:
            dropped.append(best[k].path)
            best[k] = r
        else:
            dropped.append(r.path)
    return list(best.values()), dropped


# ------------------------------------------------------------------- scores
def _in_group(category: str, group: str) -> bool:
    if group == "predicted":
        return category.startswith(PREDICTED_PREFIXES)
    if group == "predicted_no_lincs":
        return (category.startswith(PREDICTED_PREFIXES)
                and category not in LINCS_DERIVED)
    return category in MEASURED_CATEGORIES


def route_scores(routes: list, group: str) -> dict:
    """Each drug's best route in the group as one number: `{key: score}`.

    Tier decides first (a better tier always beats a worse one), then the
    route score. Drugs with no route in the group are absent and read as 0.
    """
    out: dict = {}
    for r in routes:
        if not _in_group(r.get("category") or "", group):
            continue
        key = breadth_key(r.get("drug") or "")
        if not key:
            continue
        s = (10 - int(r.get("tier") or 9)) * 1000 + float(r.get("score") or 0)
        if s > out.get(key, 0.0):
            out[key] = s
    return out


def drug_ranking(routes: list, group: str | None = None) -> list:
    """Distinct drugs, best first, each with the route that put it there.

    Ordered the way the report orders them: tier, then score. `group` limits
    the routes to one of GROUPS; None takes every category except the two
    left out of the tests (the gene's own product, and virtual screening).
    """
    best: dict = {}
    for r in routes:
        cat = r.get("category") or ""
        if group is not None and not _in_group(cat, group):
            continue
        if group is None and (cat.startswith(("direct:", "literature:"))
                              or cat == "virtual screen hypothesis"):
            continue
        key = breadth_key(r.get("drug") or "")
        if not key:
            continue
        order = (int(r.get("tier") or 9), -float(r.get("score") or 0))
        if key not in best or order < best[key][0]:
            best[key] = (order, r)
    ranked = sorted(best.items(), key=lambda kv: kv[1][0])
    return [{"key": k, "drug": r["drug"], "tier": r.get("tier"),
             "score": r.get("score"), "category": r.get("category"),
             "acts_on": r.get("acts_on")} for k, (_, r) in ranked]


def control_result(run: Run) -> dict | None:
    """Where the known drugs for a positive-control gene landed, or None when
    the run's gene is not a control or was run in the other direction."""
    spec = POSITIVE_CONTROLS.get(run.gene.upper())
    if spec is None or spec[0] != run.goal:
        return None
    goal, names, why = spec
    want = {breadth_key(n) for n in names}
    out = {"gene": run.gene, "goal": goal, "why": why,
           "expected": list(names), "path": run.path, "views": {}}
    for label, group in (("any route", None), ("predicted routes", "predicted")):
        order = drug_ranking(run.routes, group)
        hits = [(i + 1, d) for i, d in enumerate(order) if d["key"] in want]
        out["views"][label] = {
            "n_drugs": len(order),
            "found": bool(hits),
            "best_rank": hits[0][0] if hits else None,
            "best": hits[0][1] if hits else None,
            "n_expected_found": len(hits)}
    return out


def route_counts(routes: list) -> dict:
    """How many distinct routes each drug has, in every category: the
    baseline that favours drugs appearing everywhere."""
    seen: dict = {}
    for r in routes:
        key = breadth_key(r.get("drug") or "")
        if key:
            seen.setdefault(key, set()).add(
                (r.get("category"), r.get("acts_on"), r.get("mechanism")))
    return {k: float(len(v)) for k, v in seen.items()}


# ------------------------------------------------------------------ metrics
def _ranks(x: np.ndarray) -> np.ndarray:
    """Average ranks, ties shared."""
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    xs = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def auc(scores: np.ndarray, labels: np.ndarray) -> float | None:
    """Area under the ROC curve; ties count half. None without both classes."""
    pos = int(labels.sum())
    neg = len(labels) - pos
    if pos == 0 or neg == 0:
        return None
    r = _ranks(scores)
    return float((r[labels.astype(bool)].sum() - pos * (pos + 1) / 2)
                 / (pos * neg))


def permutation_p(scores: np.ndarray, labels: np.ndarray, observed: float,
                  n: int, seed: int) -> float | None:
    """One-sided: how often shuffled labels score at least `observed`."""
    if n <= 0:
        return None
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n):
        a = auc(scores, rng.permutation(labels))
        if a is not None and a >= observed - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


def _score_test(keys: list, labels: list, scores: dict, *, permutations: int,
                seed: int, baseline: dict | None = None,
                fame: dict | None = None) -> dict:
    y = np.array(labels, dtype=int)
    s = np.array([scores.get(k, 0.0) for k in keys], dtype=float)
    out = {"n_pos": int(y.sum()), "n_neg": int(len(y) - y.sum())}
    out["n_pos_with_route"] = int(sum(1 for k, v in zip(keys, y)
                                      if v and scores.get(k, 0.0) > 0))
    if out["n_pos"] < MIN_PER_CLASS or out["n_neg"] < MIN_PER_CLASS:
        out["auc"] = None
        out["reason"] = (f"needs {MIN_PER_CLASS} of each class; has "
                         f"{out['n_pos']} and {out['n_neg']}")
        return out
    a = auc(s, y)
    out["auc"] = round(a, 3)
    out["p"] = round(permutation_p(s, y, a, permutations, seed), 4) \
        if permutations else None
    if baseline is not None:
        b = auc(np.array([baseline.get(k, 0.0) for k in keys]), y)
        out["baseline_auc"] = None if b is None else round(b, 3)
    if fame is not None:
        f = np.array([fame.get(k, 0) for k in keys], dtype=float)
        fa = auc(f, y)
        out["fame_auc"] = None if fa is None else round(fa, 3)
        known = f >= FAME_FLOOR
        yk = y[known]
        if (int(yk.sum()) >= MIN_PER_CLASS
                and int(len(yk) - yk.sum()) >= MIN_PER_CLASS):
            out["matched_auc"] = round(auc(s[known], yk), 3)
            out["matched_n_pos"] = int(yk.sum())
            out["matched_n_neg"] = int(len(yk) - yk.sum())
        else:
            out["matched_auc"] = None
    return out


def top_k_breakdown(scores: dict, ks: KeySet, library: set,
                    k: int = TOP_K) -> dict:
    """What the best-scored drugs are, by what the key says of them."""
    ranked = sorted((key for key in scores if key in library),
                    key=lambda x: -scores[x])[:k]
    wanted = ks.wanted - ks.opposite
    opposite = ks.opposite - ks.wanted
    return {"k": len(ranked),
            "wanted": sum(1 for x in ranked if x in wanted),
            "opposite": sum(1 for x in ranked if x in opposite),
            "other": sum(1 for x in ranked
                         if x not in wanted and x not in opposite)}


def evaluate_key(run: Run, ks: KeySet, library: set, *, permutations: int,
                 seed: int, fame: dict | None = None) -> dict:
    """Both tests for one key, for each route group it may grade."""
    scope = library if ks.universe is None else library & ks.universe
    pos = (ks.wanted & scope) - ks.opposite
    reported = ks.wanted | ks.opposite | ks.unclear
    neg_unreported = scope - reported
    opp_only = (ks.opposite & scope) - ks.wanted
    counts = route_counts(run.routes)
    res = {"n_records": ks.n_records, "n_wanted_in_scope": len(pos),
           "n_opposite_in_scope": len(opp_only), "n_scope": len(scope),
           "groups": {}}
    for group in KEY_GROUPS[ks.name]:
        sc = route_scores(run.routes, group)
        t1_keys = sorted(pos | neg_unreported)
        t1 = _score_test(t1_keys, [1 if k in pos else 0 for k in t1_keys], sc,
                         permutations=permutations, seed=seed,
                         baseline=counts, fame=fame)
        t2_keys = sorted(pos | opp_only)
        t2 = _score_test(t2_keys, [1 if k in pos else 0 for k in t2_keys], sc,
                         permutations=permutations, seed=seed)
        res["groups"][group] = {
            "reported_or_not": t1, "right_direction": t2,
            "top": top_k_breakdown(sc, ks, library),
            "n_drugs_with_route": len([k for k in sc if k in library])}
    return res


def fame_counts(runs: list) -> dict:
    """`{drug: number of runs whose CTD record includes it}`, or {} when fewer
    than FAME_MIN_RUNS runs have a CTD record, since a count over a few genes
    says nothing about how well known a drug is."""
    sets = [r.keys["ctd"].wanted | r.keys["ctd"].opposite
            | r.keys["ctd"].unclear for r in runs if "ctd" in r.keys]
    if len(sets) < FAME_MIN_RUNS:
        return {}
    out: dict = {}
    for ks in sets:
        for k in ks:
            out[k] = out.get(k, 0) + 1
    return out


def evaluate_run(run: Run, library: set, *, permutations: int = 2000,
                 seed: int = 0, fame: dict | None = None) -> dict:
    """Every key that could be built for one run.

    `fame` is `fame_counts` over all the runs. A run's own record is taken out
    of its counts, so the baseline cannot see the answer it is graded on.
    """
    keys = {}
    for name, ks in run.keys.items():
        own = None
        if name == "ctd" and fame:
            mine = ks.wanted | ks.opposite | ks.unclear
            own = {k: c - (1 if k in mine else 0) for k, c in fame.items()}
        keys[name] = evaluate_key(run, ks, library, permutations=permutations,
                                  seed=seed, fame=own)
    return {"gene": run.gene, "goal": run.goal, "path": run.path,
            "skipped": dict(run.skipped), "keys": keys}


def summarise(per_run: list) -> dict:
    """Mean AUC over the runs that could be scored, and how many beat chance."""
    out: dict = {}
    for key, groups in KEY_GROUPS.items():
        out[key] = {}
        for group in groups:
            out[key][group] = {}
            for test in ("reported_or_not", "right_direction"):
                vals = [r["keys"][key]["groups"][group][test]
                        for r in per_run if key in r["keys"]
                        and r["keys"][key]["groups"][group][test].get("auc")
                        is not None]
                aucs = [v["auc"] for v in vals]
                entry = {"n_runs": len(aucs),
                         "mean_auc": round(float(np.mean(aucs)), 3)
                         if aucs else None,
                         "min_auc": min(aucs) if aucs else None,
                         "max_auc": max(aucs) if aucs else None,
                         "n_above_chance": sum(1 for a in aucs if a > 0.5),
                         "n_positives": sum(v["n_pos"] for v in vals),
                         "n_positives_with_route": sum(
                             v["n_pos_with_route"] for v in vals)}
                if test == "reported_or_not":
                    base = [v.get("baseline_auc") for v in vals
                            if v.get("baseline_auc") is not None]
                    entry["mean_baseline_auc"] = (
                        round(float(np.mean(base)), 3) if base else None)
                    fa = [v["fame_auc"] for v in vals
                          if v.get("fame_auc") is not None]
                    entry["mean_fame_auc"] = (round(float(np.mean(fa)), 3)
                                              if fa else None)
                    # Matched runs are scored on their own, so they are
                    # taken from every run that had a value, not only those
                    # in `vals`.
                    ma = [r["keys"][key]["groups"][group][test]
                          for r in per_run if key in r["keys"]]
                    ma = [v["matched_auc"] for v in ma
                          if v.get("matched_auc") is not None]
                    entry["mean_matched_auc"] = (round(float(np.mean(ma)), 3)
                                                 if ma else None)
                    entry["n_matched_runs"] = len(ma)
                out[key][group][test] = entry
    return out


# ------------------------------------------------------------------- output
def _cell(t: dict) -> str:
    if t.get("auc") is None:
        return f"not scored ({t.get('reason', 'no data')})"
    p = "" if t.get("p") is None else f", p {t['p']}"
    b = ("" if t.get("baseline_auc") is None
         else f", routes-per-drug baseline {t['baseline_auc']}")
    return (f"{t['auc']} ({t['n_pos']} vs {t['n_neg']}; "
            f"{t['n_pos_with_route']} of {t['n_pos']} positives have a "
            f"route{p}{b})")


def markdown(report: dict) -> str:
    L = ["# Retrospective check against known effects", "",
         "Do the routes point at drugs already recorded as moving the gene? "
         "AUC 0.5 is chance. Both answer keys are noisy and the check tests "
         "expression, not benefit. Everything below is read from saved "
         "results and local files.", ""]
    for r in report["runs"]:
        L += [f"## {r['gene']} ({r['goal']})", "", f"File: `{r['path']}`", ""]
        for key, why in r["skipped"].items():
            L.append(f"- {KEY_LABELS[key]}: not built, {why}.")
        if r["skipped"]:
            L.append("")
        for key, e in r["keys"].items():
            L += [f"### Key: {KEY_LABELS[key]}", "",
                  f"{e['n_records']} drugs recorded; {e['n_wanted_in_scope']} "
                  f"moving it the wanted way and {e['n_opposite_in_scope']} "
                  f"the other way are among the {e['n_scope']} approved "
                  f"drugs in scope.", "",
                  "| Routes | Reported or not | Right direction | Top 20: "
                  "wanted / opposite / other |", "|---|---|---|---|"]
            for g, d in e["groups"].items():
                top = d["top"]
                L.append(f"| {g}: {GROUPS[g]} | {_cell(d['reported_or_not'])}"
                         f" | {_cell(d['right_direction'])} | "
                         f"{top['wanted']} / {top['opposite']} / "
                         f"{top['other']} of {top['k']} |")
            L.append("")
    s = report["summary"]
    L += ["## Across genes", "",
          "| Key | Routes | Test | Genes scored | Mean AUC | Range | Above "
          "chance | Positives with a route |", "|---|---|---|---|---|---|---|---|"]
    for key, groups in KEY_GROUPS.items():
        for g in groups:
            for test, label in (("reported_or_not", "reported or not"),
                                ("right_direction", "right direction")):
                e = s[key][g][test]
                rng = ("-" if e["min_auc"] is None
                       else f"{e['min_auc']} to {e['max_auc']}")
                L.append(
                    f"| {key} | {g} | {label} | {e['n_runs']} | "
                    f"{'-' if e['mean_auc'] is None else e['mean_auc']} | "
                    f"{rng} | {e['n_above_chance']} of {e['n_runs']} | "
                    f"{e['n_positives_with_route']} of {e['n_positives']} |")
    if report.get("fame_used"):
        e = s["ctd"]["predicted"]["reported_or_not"]
        f = e.get("mean_fame_auc")
        m = e.get("mean_matched_auc")
        L += ["", "### Drug fame", "",
              f"A baseline that scores a drug by how many OTHER genes record "
              f"it in CTD, and knows nothing about this gene, averages "
              f"{'-' if f is None else f} on the reported-or-not test, "
              f"against {'-' if e['mean_auc'] is None else e['mean_auc']} "
              f"for the predicted routes."]
        if m is not None:
            L.append(f"Among only the drugs recorded for at least "
                     f"{FAME_FLOOR} other genes, where fame no longer "
                     f"separates the groups, the predicted routes average "
                     f"{m} over {e['n_matched_runs']} genes. That is the "
                     f"figure to trust: 0.5 there means the routes rank "
                     f"known regulators no better than chance once drug "
                     f"fame is held equal.")
        else:
            L.append(f"Too few drugs recorded for at least {FAME_FLOOR} other "
                     f"genes to repeat the test with fame held equal.")
    for key in KEY_GROUPS:
        base = s[key]["predicted_no_lincs"]["reported_or_not"].get(
            "mean_baseline_auc")
        if base is not None:
            L += ["", f"Against {key}, a baseline that scores a drug by how "
                  f"many routes it has, and knows nothing about the gene, "
                  f"averages {base} on the reported-or-not test. A predicted "
                  f"AUC that does not beat it may reflect drugs that appear "
                  f"everywhere."]
    if report.get("controls"):
        L += ["", "## Positive controls", "",
              "Genes whose regulation by an approved drug is textbook "
              "pharmacology. The table was written before the runs and is the "
              "author's judgement, not a citation. Rank is among the distinct "
              "drugs the run proposes, best first.", "",
              "| Gene | Known answer | Any route | Predicted routes | Route "
              "that found it |", "|---|---|---|---|---|"]
        for c in report["controls"]:
            cells = []
            for view in ("any route", "predicted routes"):
                v = c["views"][view]
                cells.append(
                    f"rank {v['best_rank']} of {v['n_drugs']}" if v["found"]
                    else f"absent from {v['n_drugs']}")
            best = c["views"]["any route"]["best"]
            how = (f"{best['drug']}: {best['category']} via {best['acts_on']}"
                   if best else "-")
            L.append(f"| {c['gene']} | {', '.join(c['expected'][:3])}"
                     f"{'...' if len(c['expected']) > 3 else ''} | {cells[0]}"
                     f" | {cells[1]} | {how} |")
        n = len(report["controls"])
        k = sum(1 for c in report["controls"]
                if c["views"]["any route"]["found"])
        L += ["", f"{k} of {n} controls have an expected drug among the "
              f"candidates at all."]
    L += ["", "## How to read this", "",
          "- Few genes means wide uncertainty. Look at the range and the "
          "count above chance, not only the mean.",
          "- A low AUC with few positives that have a route is a coverage "
          "result: the pipeline does not reach those drugs at all, which is "
          "different from reaching them and ranking them badly.",
          "- \"Right direction\" is the stricter test: both groups were seen "
          "moving the gene, so being well studied does not help.",
          "- A drug the key does not record is unstudied, not inert, so some "
          "drugs counted as misses may be real regulators.",
          "- Passing says the routes find known regulators. It does not say "
          "the unrecorded candidates are right, or that moving the gene "
          "treats a disease.", ""]
    if report.get("dropped"):
        L += ["Older results for the same gene and direction, not used:",
              *[f"- `{p}`" for p in report["dropped"]], ""]
    return "\n".join(L)


def run_all(paths: list, library_keys: set, *, permutations: int = 2000,
            seed: int = 0, gmt_up=None, gmt_down=None) -> dict:
    runs = [r for r in (load_run(p, gmt_up=gmt_up, gmt_down=gmt_down)
                        for p in paths) if r is not None]
    kept, dropped = latest_per_gene_goal(runs)
    fame = fame_counts(kept)
    per_run = [evaluate_run(r, library_keys, permutations=permutations,
                            seed=seed, fame=fame) for r in kept]
    controls = [c for c in (control_result(r) for r in kept) if c]
    return {"runs": per_run, "summary": summarise(per_run),
            "controls": controls,
            "dropped": dropped, "permutations": permutations, "seed": seed,
            "library_size": len(library_keys),
            "fame_used": bool(fame)}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.retrospective",
        description=("Check saved results against known effects: do the "
                     "predicted routes point at drugs known to move the "
                     "gene?"))
    p.add_argument("results", nargs="+",
                   help="results.json files or the folders that hold them")
    p.add_argument("--library", default=None,
                   help="approved drug table; the bundled one by default")
    p.add_argument("--gmt-up", default=None,
                   help="LINCS up signatures (.gmt); the local one by default")
    p.add_argument("--gmt-down", default=None,
                   help="LINCS down signatures (.gmt)")
    p.add_argument("--permutations", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="results/retrospective")
    args = p.parse_args(argv)
    if args.permutations < 0:
        p.error("--permutations cannot be negative")

    from . import localdrugs as ld
    rows = ld.library(args.library)
    if not rows:
        p.error("no approved drug table was found; pass --library")
    keys = {breadth_key(r["name"]) for r in rows} - {""}
    report = run_all(args.results, keys, permutations=args.permutations,
                     seed=args.seed, gmt_up=args.gmt_up,
                     gmt_down=args.gmt_down)
    if not report["runs"]:
        p.error("none of the given paths held a readable result")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "retrospective.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8")
    text = markdown(report)
    (out / "retrospective.md").write_text(text, encoding="utf-8")
    print(text)
    print(f"Wrote {out / 'retrospective.json'} and {out / 'retrospective.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
