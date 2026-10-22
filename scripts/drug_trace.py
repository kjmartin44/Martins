"""Why is a drug that should move this gene not in the list?

    python scripts/drug_trace.py FXN up omaveloxolone
    python scripts/drug_trace.py TERT up danazol "valproic acid" --results results/TERT_up/results.json

For each named drug this follows the chain the routes need and stops at the first
link that is missing:

  1  the drug is in the approved library
  2  it is on the FDA Orange Book, which the reports require
  3  ChEMBL records a mechanism for it at some protein (the index every route reads)
  4  one of those proteins is a recorded regulator of the gene, in CollecTRI, TRRUST,
     the knockout experiments or ChIP occupancy
  5  the drug moves that protein the way the goal needs
  6  if no protein regulates the gene, whether one is within three signed steps of a
     regulator (the upstream search's reach)

With `--results` it also says whether the drug is in a finished run, and by which
route. Nothing here changes a run. It reads the same cached tables the pipeline reads
and sends no request that is not already cached, except the one CTD lookup for the gene.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from generoutes import discovery as disc  # noqa: E402
from generoutes import drugs as drugs_mod  # noqa: E402
from generoutes import fda, knockout, localdrugs, network  # noqa: E402
from generoutes import regulon as regulon_mod  # noqa: E402
from generoutes import tfactivity  # noqa: E402
from generoutes.http import Client  # noqa: E402
from generoutes.signatures import breadth_key  # noqa: E402

MAX_HOPS = 3


def regulator_sources(gene: str, data: Path = ROOT / "data") -> dict:
    """`{factor: {"sources": {table: detail}, "sign": +1, -1 or 0}}`.

    The sign is the first one given, in the order CollecTRI, TRRUST, knockout
    data. ChIP occupancy gives none, so a factor named only there has sign 0.
    """
    gene = gene.upper()
    out: dict = {}

    def add(tf, src, detail, sign=0):
        rec = out.setdefault(tf.upper(), {"sources": {}, "sign": 0})
        rec["sources"][src] = detail
        if not rec["sign"] and sign:
            rec["sign"] = sign

    for name, label in (("collectri.tsv", "CollecTRI"), ("trrust.tsv", "TRRUST")):
        path = data / name
        if path.is_file():
            for tf, targets in tfactivity.load_signed_regulons(
                    path, min_targets=1).items():
                if gene in targets:
                    sg = targets[gene]
                    add(tf, label, {1: "raises", -1: "lowers",
                                    0: "sign not recorded"}[sg], sg)
    for tf, rec in knockout.regulators_of(gene).items():
        add(tf, "knockout",
            f"{rec['activates']} raise, {rec['represses']} lower", rec["sign"])
    chip = data / "chip_regulon.tsv"
    if chip.is_file():
        for tf, targets in regulon_mod.load_regulons(chip).items():
            if gene in targets:
                add(tf, "ChIP occupancy", "near the gene, no sign")
    return out


def reach_upstream(incoming: dict, sources, max_hops: int = MAX_HOPS) -> dict:
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


def trace(gene: str, goal: str, drug: str, *, library, index, regs, incoming,
          ctd, run_routes) -> list:
    """`[(step, ok, text)]` for one drug, ending at the first missing link."""
    lines = []
    key = breadth_key(drug)
    rows = [r for r in library if breadth_key(r["name"]) == key]
    if not rows:
        rows = [r for r in library if key and key in breadth_key(r["name"])]
    if not rows:
        lines.append((1, False, "not in the approved library, so no route can name it"))
        return lines
    lines.append((1, True, f"in the library as {rows[0]['name']} ({rows[0].get('chembl_id')})"))

    ok = fda.is_approved(drug)
    lines.append((2, ok, "on the FDA Orange Book" if ok else
                  "NOT on the Orange Book, so reports hide it"))

    names = {breadth_key(r["name"]) for r in rows}
    mech = []
    for sym, ds in index.items():
        for d in ds:
            if breadth_key(d.name) in names:
                mech.append((sym, d.direction, d.action_type,
                             "family" if d.target_type else ""))
    if not mech:
        lines.append((3, False, "ChEMBL records no mechanism for it at any protein, "
                      "so no route that reads mechanisms can reach it"))
    else:
        lines.append((3, True, f"mechanism at {len({m[0] for m in mech})} protein(s): "
                      + ", ".join(sorted({f"{m[0]} ({m[1]})" for m in mech})[:10])))
        hits = [(s, dr, a, fam) for s, dr, a, fam in mech if s in regs]
        if hits:
            for sym, dr, act, fam in hits:
                sign = regs[sym]["sign"]
                srcs = "; ".join(f"{k}: {v}" for k, v in
                                 regs[sym]["sources"].items())
                lines.append((4, True, f"{sym} is a recorded regulator ({srcs})"))
                if not sign:
                    lines.append((5, False, f"no table gives the sign of {sym}'s "
                                  f"effect, so the route cannot say which way "
                                  f"{drug} (which moves it {dr}) should go"))
                    continue
                need = "up" if (goal == "up") == (sign == 1) else "down"
                good = dr in (need, "unclear")
                lines.append((5, good, f"{drug} is recorded as moving {sym} {dr}; the "
                              f"goal needs {need}"
                              + ("" if good else " - the wrong way")))
        else:
            lines.append((4, False, "none of those proteins is a recorded regulator of "
                          "the gene in CollecTRI, TRRUST, the knockout data or ChIP "
                          "occupancy"))
            near = reach_upstream(incoming, set(regs))
            reach = [(s, near[s]) for s, *_ in mech if s in near]
            if reach:
                lines.append((6, True, "but within " + ", ".join(
                    f"{s} ({h} step{'s' if h != 1 else ''})" for s, h in reach[:5])
                    + " of a regulator, which the upstream search could reach if "
                      "that regulator is ranked and has no drug of its own"))
            else:
                lines.append((6, False, f"and none is within {MAX_HOPS} signed steps "
                              "of a regulator in the OmniPath network"))
    if ctd is not None:
        d = disc.link_direction(ctd.get(key, []) or [])
        hit = [c for c in ctd if breadth_key(c) == key]
        lines.append((7, bool(hit), f"CTD records it as moving the gene: {d}"
                      if hit else "CTD has no expression record for it on this gene"))
    if run_routes is not None:
        found = [r for r in run_routes if breadth_key(r["drug"]) == key]
        lines.append((8, bool(found), "in the finished run via " + ", ".join(
            sorted({f"{r['category']} (tier {r['tier']})" for r in found}))
            if found else "not in the finished run"))
    return lines


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("gene")
    p.add_argument("goal", choices=["up", "down"])
    p.add_argument("drugs", nargs="+")
    p.add_argument("--results", default=None, help="a finished results.json")
    p.add_argument("--no-ctd", action="store_true")
    args = p.parse_args(argv)

    client = Client()
    library = localdrugs.library()
    index = drugs_mod.approved_target_index(client, library)
    regs = regulator_sources(args.gene)
    edges = network.fetch_network(client)
    incoming, _, _ = network.build_graph(edges) if edges else ({}, {}, {})
    ctd = None if args.no_ctd else disc.ctd_expression_links(client, args.gene)
    run_routes = None
    if args.results:
        data = json.loads(Path(args.results).read_text(encoding="utf-8"))
        run_routes = []
        for blk in data.get("results") or []:
            run_routes += (blk.get("all_routes") or {}).get("routes") or []

    print(f"{args.gene.upper()} {args.goal}: {len(regs)} recorded regulator(s) in "
          f"{len({s for v in regs.values() for s in v})} table(s)\n")
    for drug in args.drugs:
        print(f"## {drug}")
        for step, ok, text in trace(args.gene, args.goal, drug, library=library,
                                    index=index, regs=regs, incoming=incoming,
                                    ctd=ctd, run_routes=run_routes):
            print(f"  {step} {'ok  ' if ok else 'MISS'} {text}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
