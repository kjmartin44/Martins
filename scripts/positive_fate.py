"""Where do the drugs the literature says move a gene get lost?

The retrospective check finds that the routes recover few of the approved drugs
CTD records as raising or lowering a gene. This follows each recorded drug that
has no route and says the first reason it was missed, so the fix can aim at the
stage that loses the most.

    python scripts/positive_fate.py results/validation --out results/retrospective_fate

For every run it takes the drugs CTD records moving the gene the wanted way
(approved drugs in the local library only, as the retrospective does), then
splits them:

  found               a predicted route already reaches the drug
  no target recorded  ChEMBL lists no mechanism target for it, so no route
                      that starts from a target can reach it
  acts on a found factor
                      its target is one of the factors this run ranked, yet
                      the drug has no route: it was filtered (direction, role,
                      approval) or graded out
  factor not found    its target is a transcription factor (one with targets
                      in CollecTRI) that the scan did not rank for this
                      promoter: the motif, the window or the score missed it
  in reach upstream   its target sits within three signalling steps upstream
                      of a ranked factor in OmniPath, yet the search did not
                      return it
  too far upstream    in the network, but more than three steps away
  not in the network  a target OmniPath has no edge for

A drug with several targets takes the stage closest to being found. The first
reason is a classification of what the data say about the drug's targets, not a
test of why the pipeline missed it: a drug can act through a target that is
not why CTD records it, and CTD includes toxicology doses and stress
responses that no regulator route would predict.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
import time
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from generoutes import chemblid, localdrugs, network  # noqa: E402
from generoutes import retrospective as retro  # noqa: E402
from generoutes.http import Client, HttpError  # noqa: E402
from generoutes.signatures import breadth_key  # noqa: E402

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
MAX_HOPS = 3

FOUND = "found"
NO_TARGET = "no target recorded"
FOUND_FACTOR = "acts on a found factor"
FACTOR_MISSED = "factor not found"
IN_REACH = "in reach upstream"
TOO_FAR = "too far upstream"
OFF_NETWORK = "not in the network"
ORDER = [FOUND, FOUND_FACTOR, IN_REACH, FACTOR_MISSED, TOO_FAR, OFF_NETWORK,
         NO_TARGET]
# The stage closest to being found wins when a drug has several targets.
CLOSENESS = {FOUND_FACTOR: 0, IN_REACH: 1, FACTOR_MISSED: 2, TOO_FAR: 3,
             OFF_NETWORK: 4, NO_TARGET: 5}


def target_symbols(client, chembl_ids, cache: dict, failures: list) -> set:
    """Gene symbols of the targets in a drug's ChEMBL mechanism records.

    `chembl_ids` is every library row with the drug's name once salts are
    ignored, because a mechanism is recorded on one form of a drug (lapatinib
    ditosylate and not lapatinib) and the package's own index reads each row.
    A failed request is retried, never cached as "no target", and counted.
    """
    key = tuple(sorted(chembl_ids))
    if key in cache:
        return cache[key]
    out: set = set()
    failed = False
    for mid in key:
        res = None
        for attempt in range(3):
            try:
                res = client.get(f"{CHEMBL}/mechanism.json",
                                 {"molecule_chembl_id": mid, "limit": 100})
                break
            except HttpError:
                time.sleep(2 * (attempt + 1))
        if res is None:
            failed = True
            continue
        for m in (res or {}).get("mechanisms", []):
            tid = m.get("target_chembl_id")
            if tid:
                out |= component_symbols(client, tid, cache, failures)
    if failed and not out:
        failures.append(key)
        return out                      # not cached: it is "could not ask"
    cache[key] = out
    return out


def component_symbols(client, target_id: str, cache: dict,
                      failures: list) -> set:
    key = ("t", target_id)
    if key in cache:
        return cache[key]
    syms: set = set()
    t = None
    for attempt in range(3):
        try:
            t = client.get(f"{CHEMBL}/target/{target_id}.json")
            break
        except HttpError:
            time.sleep(2 * (attempt + 1))
    if t is None:
        failures.append(("target", target_id))
        return syms
    for comp in (t or {}).get("target_components", []) or []:
        for syn in comp.get("target_component_synonyms", []) or []:
            if syn.get("syn_type") == "GENE_SYMBOL":
                syms.add(str(syn.get("component_synonym", "")).upper())
    cache[key] = syms
    return syms


def upstream_reach(incoming: dict, sources, max_hops: int = MAX_HOPS) -> dict:
    """`{node: hops}` for every node within `max_hops` steps upstream of any of
    `sources`, walking edges backwards."""
    dist = {s: 0 for s in sources}
    queue = deque(sources)
    while queue:
        node = queue.popleft()
        if dist[node] >= max_hops:
            continue
        for e in incoming.get(node, []):
            if e.source not in dist:
                dist[e.source] = dist[node] + 1
                queue.append(e.source)
    return dist


def classify(targets: set, ranked: set, tfs: set, reach: dict,
             in_network: set) -> str:
    """The stage closest to being found over a drug's targets."""
    if not targets:
        return NO_TARGET
    best = None
    for t in targets:
        if t in ranked:
            c = FOUND_FACTOR
        elif t in reach:
            c = IN_REACH
        elif t in tfs:
            c = FACTOR_MISSED
        elif t in in_network:
            c = TOO_FAR
        else:
            c = OFF_NETWORK
        if best is None or CLOSENESS[c] < CLOSENESS[best]:
            best = c
    return best


def load_tfs(path: Path) -> set:
    if not path.exists():
        return set()
    with path.open(encoding="utf-8") as fh:
        return {r["source"].upper() for r in csv.DictReader(fh, delimiter="\t")}


def analyse_run(run, lib: dict, client, cache, incoming, in_network, tfs,
                failures):
    """Fate of each approved drug CTD records as moving the gene the right way."""
    ks = run.keys.get("ctd")
    if ks is None:
        return None
    blocks = json.loads(Path(run.path).read_text(encoding="utf-8"))["results"]
    ranked = {t["symbol"].upper() for b in blocks for t in b.get("tfs") or []}
    reach = upstream_reach(incoming, ranked)
    have = retro.route_scores(run.routes, "predicted")
    positives = sorted(k for k in ks.wanted if k in lib)
    rows = []
    for k in positives:
        if k in have:
            rows.append({"drug": lib[k]["name"], "fate": FOUND,
                         "targets": []})
            continue
        targets = target_symbols(client, lib[k]["ids"], cache, failures)
        rows.append({"drug": lib[k]["name"], "targets": sorted(targets),
                     "fate": classify(targets, ranked, tfs, reach,
                                      in_network)})
    return {"gene": run.gene, "goal": run.goal, "n_ranked_factors": len(ranked),
            "n_positives": len(positives), "drugs": rows}


def summarise(per_run: list) -> dict:
    total = collections.Counter()
    by_target = collections.defaultdict(collections.Counter)
    for r in per_run:
        for d in r["drugs"]:
            total[d["fate"]] += 1
            for t in d["targets"]:
                by_target[d["fate"]][t] += 1
    return {"total": dict(total), "n": sum(total.values()),
            "top_targets": {f: by_target[f].most_common(8) for f in ORDER
                            if f != FOUND and by_target[f]}}


def markdown(per_run: list, summ: dict) -> str:
    n = summ["n"]
    L = ["# Where the literature's drugs get lost", "",
         f"{n} approved drugs that CTD records moving a gene the wanted way, "
         f"across {len(per_run)} runs. Each is placed at the first stage that "
         f"explains why it has no route, or marked found.", "",
         "| Stage | Drugs | Share |", "|---|---|---|"]
    for f in ORDER:
        c = summ["total"].get(f, 0)
        L.append(f"| {f} | {c} | {c / n:.1%} |")
    L += ["", "A drug with several targets takes the stage closest to being "
          "found. This classifies what the data say about each drug's "
          "targets. It does not prove why a route was missed, and CTD "
          "includes toxicology and stress effects that no regulator route "
          "would predict.", ""]
    L += ["## Targets that account for the misses", ""]
    for f, items in summ["top_targets"].items():
        L.append(f"- **{f}**: " + ", ".join(f"{t} ({c})" for t, c in items))
    L += ["", "## By gene", "", "| Gene | Positives | Found | "
          + " | ".join(ORDER[1:]) + " |", "|---|---|---|"
          + "---|" * (len(ORDER) - 1)]
    for r in per_run:
        c = collections.Counter(d["fate"] for d in r["drugs"])
        L.append(f"| {r['gene']} ({r['goal']}) | {r['n_positives']} | "
                 f"{c.get(FOUND, 0)} | "
                 + " | ".join(str(c.get(f, 0)) for f in ORDER[1:]) + " |")
    L.append("")
    return "\n".join(L)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("results", nargs="+")
    p.add_argument("--library", default=None)
    p.add_argument("--collectri", default=str(ROOT / "data" / "collectri.tsv"))
    p.add_argument("--out", default="results/retrospective_fate")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)

    rows = localdrugs.library(args.library)
    if not rows:
        p.error("no approved drug table found; pass --library")
    lib = {}
    for r in rows:
        k = breadth_key(r["name"])
        if not k or not r.get("chembl_id"):
            continue
        row = lib.setdefault(k, dict(r, ids=set()))
        row["ids"].add(r["chembl_id"])

    runs = []
    for path in args.results:
        pth = Path(path)
        files = [pth] if pth.is_file() else sorted(pth.glob("*/results.json"))
        runs += [x for x in (retro.load_run(f) for f in files) if x]
    kept, _ = retro.latest_per_gene_goal(runs)
    if not kept:
        p.error("no readable results")

    client = Client(offline=args.offline)
    edges = network.fetch_network(client)
    if not edges:
        p.error("could not load the OmniPath network (needs a cached copy or "
                "network access)")
    incoming, outgoing, _ = network.build_graph(edges)
    in_network = set(incoming) | set(outgoing)
    tfs = load_tfs(Path(args.collectri))
    cache: dict = {}
    failures: list = []
    per_run = []
    for i, run in enumerate(kept, 1):
        got = analyse_run(run, lib, client, cache, incoming, in_network, tfs,
                          failures)
        if got:
            per_run.append(got)
            c = collections.Counter(d["fate"] for d in got["drugs"])
            print(f"[{i}/{len(kept)}] {run.gene} {run.goal}: "
                  f"{got['n_positives']} positives, "
                  f"{c.get(FOUND, 0)} found", flush=True)
    summ = summarise(per_run)
    summ["lookup_failures"] = len(failures)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "fate.json").write_text(
        json.dumps({"summary": summ, "runs": per_run}, indent=2),
        encoding="utf-8")
    text = markdown(per_run, summ)
    (out / "fate.md").write_text(text, encoding="utf-8")
    print()
    print(text)
    print(f"Wrote {out / 'fate.json'} and {out / 'fate.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
