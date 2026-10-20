"""Why drugs a target away from a ranked factor are not returned by the search.

`positive_fate.py` puts many drugs at "in reach upstream": their target sits
within three signalling steps of a ranked factor, yet no upstream route lists
the drug. This replays the filters of `network.find_druggable_upstream` on those
drugs and says which one removes each:

    python scripts/upstream_misses.py <results.json files> \
        --fate results/retrospective_fate/fate.json

The stages, in the order the search applies them:

  outside the factors searched   the search covers the first eight ranked
                                 factors that have no direct drug of their own
  wrong sign                     no path of three steps or fewer whose signs
                                 multiply to the effect the goal needs
  not in the approved-drug index the drug has no recorded mechanism at that
                                 target in the index the search reads
  moves the node the wrong way   the drug's action at the node is opposite to
                                 what the path needs
  passes every test              nothing above removes it, so the search's
                                 limits do (ten druggable nodes per factor and
                                 five drugs per node)

It reads the same OmniPath network and approved-drug index the pipeline does,
from the cache. It does not run the random walk that orders nodes, so it cannot
say which cap removed a drug that passes every test, only that the filters did
not.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import positive_fate as pf  # noqa: E402
from promoter_tf_drugs import drugs as drugs_mod  # noqa: E402
from promoter_tf_drugs import localdrugs, network  # noqa: E402
from promoter_tf_drugs.http import Client  # noqa: E402
from promoter_tf_drugs.signatures import breadth_key  # noqa: E402

N_SEARCHED = network.UPSTREAM_FACTORS

OUTSIDE = "reaches only a factor outside the eight searched"
UNREACHED = "no searched factor within three steps"
WRONG_SIGN = "no path with the net effect the goal needs"
NOT_INDEXED = "drug not in the approved-drug index for that target"
WRONG_WAY = "drug moves the node the wrong way for that path"
SURVIVES = "passes every filter, so the limits left it out"
ORDER = [OUTSIDE, UNREACHED, WRONG_SIGN, NOT_INDEXED, WRONG_WAY, SURVIVES]


def signed_reach(incoming: dict, factor: str, max_hops: int = 3) -> dict:
    """`{node: {net signs}}` over signed paths of up to `max_hops` edges into
    `factor`, walking edges backwards."""
    out: dict = collections.defaultdict(set)
    queue = deque([(factor, 1, 0, (factor,))])
    seen: set = set()
    while queue:
        node, sign, hops, chain = queue.popleft()
        if hops >= max_hops:
            continue
        for e in incoming.get(node, []):
            if e.source in chain:
                continue
            net = sign * e.sign
            state = (e.source, net, hops + 1)
            if state in seen:
                continue
            seen.add(state)
            out[e.source].add(net)
            queue.append((e.source, net, hops + 1, chain + (e.source,)))
    return out


def classify(drug: str, targets: set, searched: dict, ranked_reach: dict,
             index_dirs: dict) -> str:
    """The first filter that removes a drug, or SURVIVES."""
    key = breadth_key(drug)
    if not any(t in sr for _, sr in searched.values() for t in targets):
        return OUTSIDE if any(t in ranked_reach for t in targets) else UNREACHED
    verdict = WRONG_SIGN
    for want, sr in searched.values():
        for t in targets:
            for net in sr.get(t, set()):
                if net != want:
                    continue
                dirs = index_dirs.get(key, {}).get(t, set())
                need = "up" if want * net == 1 else "down"
                if not dirs:
                    verdict = NOT_INDEXED
                elif need in dirs or "unclear" in dirs:
                    return SURVIVES
                elif verdict == WRONG_SIGN:
                    verdict = WRONG_WAY
    return verdict


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("results", nargs="+", help="results.json files")
    p.add_argument("--fate", default="results/retrospective_fate/fate.json")
    p.add_argument("--library", default=None)
    p.add_argument("--out", default="results/retrospective_fate")
    args = p.parse_args(argv)

    fate = json.loads(Path(args.fate).read_text(encoding="utf-8"))
    byrun = {(r["gene"], r["goal"]): r for r in fate["runs"]}
    client = Client()
    edges = network.fetch_network(client)
    if not edges:
        p.error("could not load the OmniPath network")
    incoming, _, _ = network.build_graph(edges)
    rows = localdrugs.library(args.library)
    index = drugs_mod.approved_target_index(client, rows)
    index_dirs: dict = collections.defaultdict(dict)
    for sym, ds in index.items():
        for d in ds:
            index_dirs[breadth_key(d.name)].setdefault(sym, set()).add(
                d.direction)

    counts = collections.Counter()
    wide = collections.Counter()
    targets = collections.defaultdict(collections.Counter)
    for path in args.results:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        q = data["query"]
        run = byrun.get((q["gene"].upper(), q["goal"]))
        if not run:
            continue
        tfs = data["results"][0].get("tfs") or []
        ranked = {t["symbol"].upper() for t in tfs}
        def wanted(t):
            role = t.get("role") or "activator"
            return 1 if (role in ("activator", "dual")) == (
                q["goal"] == "up") else -1

        searched = {}
        for t in [t for t in tfs if not t.get("drugs")][:N_SEARCHED]:
            sym = t["symbol"].upper()
            searched[sym] = (wanted(t), signed_reach(incoming, sym))
        # The wider search: every ranked factor, with or without a drug of its
        # own, and no cap on how many.
        everything = dict(searched)
        for t in tfs:
            sym = t["symbol"].upper()
            if sym not in everything:
                everything[sym] = (wanted(t), signed_reach(incoming, sym))
        ranked_reach = pf.upstream_reach(incoming, ranked)
        for d in run["drugs"]:
            if d["fate"] != pf.IN_REACH:
                continue
            tg = set(d["targets"])
            why = classify(d["drug"], tg, searched, ranked_reach, index_dirs)
            counts[why] += 1
            for t in tg:
                targets[why][t] += 1
            wide[classify(d["drug"], tg, everything, ranked_reach,
                          index_dirs)] += 1
    n = sum(counts.values())
    L = ["# Why drugs within reach upstream are not returned", "",
         f"{n} drugs sit at \"in reach upstream\" in the funnel. Each is placed "
         f"at the first filter of the upstream search that removes it.", "",
         "| Filter | Eight factors searched | Every ranked factor |",
         "|---|---|---|"]
    for k in ORDER:
        L.append(f"| {k} | {counts.get(k, 0)} ({counts.get(k, 0) / n:.1%}) | "
                 f"{wide.get(k, 0)} ({wide.get(k, 0) / n:.1%}) |")
    L += ["", "The second column asks every ranked factor, with or without a "
          "drug of its own and with no cap on how many. Drugs that pass every "
          "filter there are the most a wider search could add before the "
          "per-factor limits apply. They are an upper bound, because the "
          "limits and the node ordering are not replayed."]
    L += ["", "Targets most often involved:", ""]
    for k in ORDER:
        if targets[k]:
            L.append(f"- **{k}**: " + ", ".join(
                f"{t} ({c})" for t, c in targets[k].most_common(8)))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "upstream_misses.md").write_text("\n".join(L) + "\n",
                                            encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
