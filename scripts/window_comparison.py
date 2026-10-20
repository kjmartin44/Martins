"""Run a set of genes at two promoter windows, for comparing drug recovery.

    python scripts/window_comparison.py run <results.json files> --out results/window_cmp
    python scripts/window_comparison.py score --out results/window_cmp

`run` reads the gene and direction of each saved result and runs the pipeline for
it twice with the options that drive drug recovery (the motif scan, ChIP, the
factor drugs, the upstream network, the known regulators route and LINCS) and
without the slow extras (virtual screens, designed modalities, DNA-directed agents,
genetics, the sensitivity scan). The two arms are:

  default   -1000/+200, the window the launcher uses, with no widening
  wide      -5000/+1000, with no widening

`score` scores each arm against CTD with `retrospective.py`, and also an adaptive
arm: the default window, with the wide one used for a gene whose default scan shows
no signal when the wide scan shows some. That is what `--widen` does inside a run,
worked out here from the two arms so that each gene is run twice and not three
times. Runs are skipped if their folder already holds a `results.json`.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ARMS = {"default": (1000, 200), "wide": (5000, 1000)}

OPTIONS = [
    "--pfm-file", "data/jaspar2024.txt", "--regulon-file", "data/chip_regulon.tsv",
    "--chip-bed", "data/remap2022.bed.gz", "--gtf", "data/gencode.basic.gtf.gz",
    "--pvalue", "1e-5", "--all-tss", "--tf-weight", "0.7", "--predict-direction",
    "--upstream-network", "--upstream-depth", "3", "--regulator-route",
    "--regulator-steps", "2", "--lincs", "--lincs-gmt",
    "data/lincs_chem_pert_up.gmt", "data/lincs_chem_pert_down.gmt",
    "--lincs-geneinfo", "data/lincs_geneinfo.txt", "--no-widen",
]


def genes_from(paths) -> list:
    """`[(gene, goal)]`, one per gene and direction, from saved results."""
    seen, out = set(), []
    for p in paths:
        p = Path(p)
        if p.is_dir():
            p = p / "results.json"
        try:
            q = json.loads(p.read_text(encoding="utf-8"))["query"]
        except (OSError, ValueError, KeyError):
            continue
        key = (q["gene"].upper(), q["goal"])
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def one_run(gene: str, goal: str, arm: str, out: Path) -> str:
    dest = out / f"{gene}_{goal}_{arm}"
    if (dest / "results.json").is_file():
        return f"{gene} {arm}: already done"
    dest.mkdir(parents=True, exist_ok=True)
    up, down = ARMS[arm]
    cmd = [sys.executable, "-m", "promoter_tf_drugs", "--gene", gene,
           "--direction", goal, "--upstream", str(up), "--downstream", str(down),
           *OPTIONS, "--out", str(dest), "--quiet"]
    with open(dest / "run.log", "w", encoding="utf-8") as log:
        rc = subprocess.run(cmd, cwd=ROOT, stdout=log,
                            stderr=subprocess.STDOUT).returncode
    return f"{gene} {arm}: exit {rc}"


def run_all(paths, out: Path, workers: int) -> None:
    jobs = [(g, goal, arm) for g, goal in genes_from(paths) for arm in ARMS]
    print(f"{len(jobs)} runs, {workers} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for line in pool.map(lambda j: one_run(*j, out), jobs):
            print(line, flush=True)


def score(out: Path, permutations: int = 500) -> str:
    import numpy as np  # noqa: F401
    from promoter_tf_drugs import localdrugs
    from promoter_tf_drugs import retrospective as rt
    from promoter_tf_drugs.signatures import breadth_key

    library = {breadth_key(r["name"]) for r in localdrugs.library()} - {""}

    def verdict(folder):
        data = json.loads((folder / "results.json").read_text(encoding="utf-8"))
        qs = [b.get("scan_quality") or {} for b in data["results"]]
        order = {"no signal": 0, "weak": 1, "enriched": 2}
        return min((order.get(q.get("verdict"), 0) for q in qs), default=0), qs

    pairs = []
    for d in sorted(out.glob("*_default")):
        w = out / (d.name[:-len("default")] + "wide")
        if (d / "results.json").is_file() and (w / "results.json").is_file():
            pairs.append((d, w))
    adaptive, switched = [], []
    for d, w in pairs:
        vd, _ = verdict(d)
        vw, _ = verdict(w)
        use_wide = vd == 0 and vw > 0
        adaptive.append(w if use_wide else d)
        if use_wide:
            switched.append(d.name[:-len("_default")])

    def evaluate(folders):
        runs = [r for r in (rt.load_run(f) for f in folders) if r]
        runs, _ = rt.latest_per_gene_goal(runs)
        runs = [r for r in runs if "ctd" in r.keys]
        fame = rt.fame_counts(runs)
        per = [rt.evaluate_run(r, library, permutations=permutations, fame=fame)
               for r in runs]
        return rt.summarise(per)["ctd"]["predicted"], len(runs)

    rows = []
    for name, folders in (("default", [d for d, _ in pairs]),
                          ("wide", [w for _, w in pairs]),
                          ("adaptive", adaptive)):
        s, n = evaluate(folders)
        ro, rd = s["reported_or_not"], s["right_direction"]
        rows.append((name, n, ro["n_positives_with_route"], ro["n_positives"],
                     rd["mean_auc"], ro["mean_auc"], ro.get("mean_matched_auc")))
    counts = {"default": {}, "wide": {}}
    for d, w in pairs:
        for arm, f in (("default", d), ("wide", w)):
            v, _ = verdict(f)
            k = ["no signal", "weak", "enriched"][v]
            counts[arm][k] = counts[arm].get(k, 0) + 1
    L = ["# Promoter window comparison", "",
         f"{len(pairs)} genes run at both windows. Scored against CTD as in "
         f"`retrospective.py`. `adaptive` is the default window, with the wide one "
         f"for a gene whose default scan shows no signal when the wide one shows "
         f"some ({len(switched)} genes: {', '.join(switched) or 'none'}).", "",
         "| Arm | Genes | Wanted drugs with a route | Of | Mean AUC, right direction "
         "| Mean AUC, reported or not | Matched to fame |",
         "|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append("| " + " | ".join(str(x) for x in r) + " |")
    L += ["", "Scan verdict by arm (worst verdict across a gene's windows):", "",
          "| Arm | no signal | weak | enriched |", "|---|---|---|---|"]
    for arm in ("default", "wide"):
        c = counts[arm]
        L.append(f"| {arm} | {c.get('no signal', 0)} | {c.get('weak', 0)} | "
                 f"{c.get('enriched', 0)} |")
    text = "\n".join(L) + "\n"
    (out / "window_comparison.md").write_text(text, encoding="utf-8")
    return text


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("results", nargs="+")
    r.add_argument("--out", default="results/window_cmp")
    r.add_argument("--workers", type=int, default=6)
    s = sub.add_parser("score")
    s.add_argument("--out", default="results/window_cmp")
    s.add_argument("--permutations", type=int, default=500)
    args = p.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.cmd == "run":
        run_all(args.results, out, args.workers)
    else:
        print(score(out, args.permutations))
    return 0


if __name__ == "__main__":
    sys.exit(main())
