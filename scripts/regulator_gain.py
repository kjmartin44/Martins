"""What the start-from-regulators route adds, checked against CTD.

For each saved run this builds the regulator route from CollecTRI and the
knockout files, adds it to the run's saved routes, and scores the result the
way `retrospective.py` does:

    python scripts/regulator_gain.py <results.json files>

Variants, so the sources can be told apart:

  collectri        CollecTRI signed edges only
  knockout, TFs    knockout data only, factors that are CollecTRI regulators
  knockout, any    knockout data only, any perturbed gene
  both             CollecTRI and knockout (TFs only) together

CollecTRI and CTD draw on the same literature, so a gain from the first
variant is not an independent confirmation. The knockout variants use
measured perturbations and are the cleaner test. The LINCS key is not used,
because the knockout files and the LINCS compound files come from the same
platform.

It reads the OmniPath-free parts of the pipeline only: the saved routes, the
CollecTRI table, the knockout files and the approved-drug index from the
cache. It sends no request that is not already cached.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from promoter_tf_drugs import catalogue as cat  # noqa: E402
from promoter_tf_drugs import drugs as drugs_mod  # noqa: E402
from promoter_tf_drugs import localdrugs, regulators  # noqa: E402
from promoter_tf_drugs import retrospective as rt  # noqa: E402
from promoter_tf_drugs import tfactivity as tfa  # noqa: E402
from promoter_tf_drugs.http import Client  # noqa: E402
from promoter_tf_drugs.signatures import breadth_key  # noqa: E402

VARIANTS = ["baseline", "collectri", "knockout, TFs", "knockout, any",
            "both", "knockout, TFs, 1 experiment",
            "knockout, TFs, 3 experiments", "collectri, two steps",
            "both, two steps", "both, two steps, measured"]
MIN_RESPONSES = {"knockout, TFs, 1 experiment": 1,
                 "knockout, TFs, 3 experiments": 3}
BOOTSTRAPS = 5000
# Set by main() for the variant that reads measured activities.
BINDERS = None


def bootstrap_share(per_run: list, chance: list, seed: int = 0) -> tuple:
    """`(share, low, high, gain_low, gain_high)`.

    `per_run` is one `(wanted, opposite)` pair per gene and `chance` the
    matching pair for every approved drug CTD classes for that gene. Genes are
    resampled, not drugs, because the drugs of one gene share its biology. The
    gain is the share minus the chance share in the same resample.
    """
    import random
    rng = random.Random(seed)
    n = len(per_run)

    def share(idx, table):
        w = sum(table[i][0] for i in idx)
        o = sum(table[i][1] for i in idx)
        return w / (w + o) if (w + o) else None
    full = share(range(n), per_run)
    shares, gains = [], []
    for _ in range(BOOTSTRAPS):
        idx = [rng.randrange(n) for _ in range(n)]
        s, c = share(idx, per_run), share(idx, chance)
        if s is None or c is None:
            continue
        shares.append(s)
        gains.append(s - c)
    if not shares:
        return full, None, None, None, None
    shares.sort()
    gains.sort()
    lo, hi = int(0.025 * len(shares)), int(0.975 * len(shares)) - 1
    return full, shares[lo], shares[hi], gains[lo], gains[hi]


def regulator_routes(variant, gene, goal, skip, signed, index, ko_paths):
    """Saved-route dicts for one variant of the route."""
    if variant == "baseline":
        return []
    regs = regulators.candidates(
        gene, signed, knockout_paths=ko_paths, skip=skip,
        tf_only=(variant != "knockout, any"),
        min_responses=MIN_RESPONSES.get(variant, regulators.MIN_RESPONSES),
        use_collectri=variant in ("collectri", "both", "collectri, two steps",
                                  "both, two steps",
                                  "both, two steps, measured"),
        use_knockout=variant not in ("collectri", "collectri, two steps"),
        steps=2 if "two steps" in variant else 1)
    found = regulators.find_drugs(
        regs, goal, index,
        binders=BINDERS if variant.endswith("measured") else None)
    routes = cat.from_regulators(found)
    out = []
    for r in routes:
        r.score = cat._score(r)
        out.append(r.to_dict())
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("results", nargs="+")
    p.add_argument("--library", default=None)
    p.add_argument("--permutations", type=int, default=500)
    p.add_argument("--out", default="results/retrospective_fate")
    args = p.parse_args(argv)

    signed, _ = tfa.autoload_regulons()
    if not signed:
        p.error("no signed regulon table (CollecTRI) found in data/")
    ko_paths = regulators.knockout_mod.resolve_paths()
    client = Client()
    rows = localdrugs.library(args.library)
    index = drugs_mod.approved_target_index(client, rows)
    library = {breadth_key(r["name"]) for r in rows} - {""}
    cache: dict = {}

    def binders(sym):
        if sym not in cache:
            try:
                cache[sym] = drugs_mod.measured_binders(client, sym, rows)
            except Exception as exc:                          # noqa: BLE001
                print(f"  measured binders for {sym} failed: {str(exc)[:60]}")
                cache[sym] = []
        return cache[sym]
    global BINDERS
    BINDERS = binders

    runs = [r for r in (rt.load_run(x) for x in args.results) if r]
    runs, _ = rt.latest_per_gene_goal(runs)
    runs = [r for r in runs if "ctd" in r.keys]
    fame = rt.fame_counts(runs)

    summary = {}
    tally = {}
    chance = [0, 0]
    chance_runs: list = []
    gene_pairs: dict = {}
    for variant in VARIANTS:
        per_run = []
        gene_pairs[variant] = []
        new_wanted = new_opp = new_other = n_new = n_with = 0
        for run in runs:
            data = json.loads(Path(run.path).read_text(encoding="utf-8"))
            tfs = data["results"][0].get("tfs") or []
            skip = {t["symbol"].upper() for t in tfs if t.get("drugs")}
            extra = regulator_routes(variant, run.gene, run.goal, skip,
                                     signed, index, ko_paths)
            aug = rt.Run(run.gene, run.goal, run.path,
                         list(run.routes) + extra, run.keys, run.skipped)
            per_run.append(rt.evaluate_run(
                aug, library, permutations=args.permutations, fame=fame))
            ks = run.keys["ctd"]
            wanted, opp = ks.wanted - ks.opposite, ks.opposite - ks.wanted
            if variant == "baseline":
                # The routes the run already had, for comparison.
                keys = set(rt.route_scores(run.routes, "predicted")) & library
                chance[0] += len(wanted & library)
                chance[1] += len(opp & library)
                chance_runs.append((len(wanted & library),
                                    len(opp & library)))
            else:
                keys = {breadth_key(r["drug"]) for r in extra} & library
            gene_pairs[variant].append((len(keys & wanted), len(keys & opp)))
            n_with += bool(keys)
            n_new += len(keys)
            new_wanted += len(keys & wanted)
            new_opp += len(keys & opp)
            new_other += len(keys - wanted - opp)
        summary[variant] = rt.summarise(per_run)["ctd"]["predicted"]
        tally[variant] = (n_with, n_new, new_wanted, new_opp, new_other)

    L = ["# What the regulator route adds", "",
         f"{len(runs)} runs with a CTD record. Scored on the 'predicted' "
         f"route group against CTD. `found` counts approved drugs CTD records "
         f"as moving the gene the wanted way that have a route in that group.",
         "",
         "| Variant | Runs with new routes | New drugs | Wanted | Opposite | "
         "Not in CTD | Right-direction share | Found | Mean AUC, right "
         "direction | Mean AUC, reported or not | Matched to fame |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for v in VARIANTS:
        s = summary[v]
        n_with, n_new, w, o, ot = tally[v]
        share = f"{w / (w + o):.0%}" if (w + o) else "-"
        rd, ro = s["right_direction"], s["reported_or_not"]
        L.append(
            f"| {v} | {n_with} | {n_new} | {w} | {o} | {ot} | {share} | "
            f"{ro['n_positives_with_route']} of {ro['n_positives']} | "
            f"{rd['mean_auc']} | {ro['mean_auc']} | "
            f"{ro.get('mean_matched_auc')} |")
    L += ["", "### How sure the direction rates are", "",
          "Genes are resampled 5,000 times, not drugs, because the drugs of "
          "one gene share its biology. The gain is the right-direction share "
          "minus the share an approved drug picked at random would have in "
          "the same resample. An interval that includes zero does not "
          "separate the route from chance.", "",
          "| Variant | Right-direction share | 95% interval | Gain over "
          "chance, 95% interval |", "|---|---|---|---|"]
    boots = {}
    for v in VARIANTS:
        full, lo, hi, glo, ghi = bootstrap_share(gene_pairs[v], chance_runs)
        boots[v] = [full, lo, hi, glo, ghi]
        if full is None or lo is None:
            L.append(f"| {v} | - | - | - |")
        else:
            L.append(f"| {v} | {full:.0%} | {lo:.0%} to {hi:.0%} | "
                     f"{glo:+.0%} to {ghi:+.0%} |")
    L += ["", "The baseline row counts the drugs the run's own predicted "
          "routes already had, so its drug columns are what the new routes "
          "should be compared with. Right-direction share is wanted / "
          "(wanted + opposite) among the drugs CTD has an opinion on. With "
          f"no information it would be {chance[0] / (chance[0] + chance[1]):.0%}"
          f", the share among every approved drug CTD classes for these "
          f"genes ({chance[0]} wanted, {chance[1]} opposite)."]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "regulator_gain.md").write_text("\n".join(L) + "\n",
                                           encoding="utf-8")
    (out / "regulator_gain.json").write_text(
        json.dumps({"summary": summary, "tally": tally,
                    "bootstrap": boots}, indent=2,
                   default=str), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
