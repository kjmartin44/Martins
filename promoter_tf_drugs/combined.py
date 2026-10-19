"""One ranked list from docking and the rest of the analysis.

Why this exists
---------------
The run ranks routes by evidence tier, then score. Docking hits against the
gene's own protein enter that ranking as tier 4 ("a model says so"), so a
pocket whose docking recovered its own known ligands (CCR5: AUC 0.93, EF1 9.7)
sits below every curated transcription factor route, however weak. Those two
kinds of evidence answer different questions and are measured differently, and
no single tier scale compares them. This keeps them as separate rankings and
fuses them.

How it fuses
------------
Each source gives its own ordered list of drugs:

  routes     the panel's surviving routes if a panel file is given (the
             funnel has already applied its off-target criteria), else the
             run's routes. Ordered by tier, then score. Docking-derived routes
             are left out here so docking is not counted twice.
  docking    one list per docking file whose validation gate passed, in the
             order the docking run ranked them. A docking file whose gate did
             not pass is not used and is listed with the reason.

A drug's combined score is the sum over the sources that list it of
1 / (60 + rank), reciprocal rank fusion. It uses only ranks, so the different
scales never meet; it has one constant, the conventional 60, which is not
tuned here; and a drug in several lists beats one in a single list at the same
rank, which is what agreement between independent methods should do. Tied
drugs share the best rank of their tie group, so a tie is not rewarded with a
lucky position.

What it does not do
-------------------
**It is not a probability.** The score orders drugs. It says nothing about how
likely any of them is to work.

**The lists answer different questions.** A docking hit is a compound shaped
to fit the protein; a route through a factor is a hypothesis about changing
the gene's expression. They are put in one list so none is buried, and the
columns say which questions each drug answered. Docking says binding, not
whether the compound activates or blocks, unless the run recovered the
conformation; that is flagged on every docking row.

**Fusion does not add evidence.** A drug high in both lists is high in two
rankings, each as weak as it was alone. Rankings from methods that share an
input are not independent, and nothing here checks for that.

**No source is weighted.** A pocket that validated well and one that barely
passed count the same. The validation numbers are shown beside each docking
source so the reader can discount one.

Near-identical drugs
--------------------
Close analogues score alike and would fill the top with one idea, so they are
folded into the best-ranked one (see `clusters`): structures within Tanimoto 0.75
of it, or the same drug under another ester name. The folded drugs are named
beside their leader with their own ranks, and the full list is kept in the
JSON. Similar structure is not equivalent behaviour, so a folded drug is a
neighbour to check, not a duplicate to discard.

Each row also carries one sentence on how that drug would act on the gene
(`rationale`). It states the mechanism the route assumes, which for most routes
is not measured.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path

from .signatures import breadth_key

log = logging.getLogger(__name__)

RRF_K = 60


# ------------------------------------------------------------------ sources
def _competition_ranks(keys_in_order: list, sort_keys: list) -> dict:
    """`{key: rank}` where equal sort keys share the best rank of the group."""
    out, first = {}, {}
    for i, (k, sk) in enumerate(zip(keys_in_order, sort_keys), start=1):
        first.setdefault(sk, i)
        out[k] = first[sk]
    return out


def routes_source(routes: list, gene: str | None = None) -> dict:
    """Ranking from route dicts: `{"ranks": {key: rank}, "best": {key: route}}`.

    Routes made by docking are skipped (docking has its own source). When
    `gene` is given, routes carrying another gene's name are skipped too,
    since a panel file holds the control gene's routes as well.
    """
    best: dict = {}
    for r in routes:
        if r.get("source") == "docking":
            continue
        if gene and r.get("gene") and str(r["gene"]).upper() != gene.upper():
            continue
        key = breadth_key(r.get("drug") or "")
        if not key:
            continue
        order = (int(r.get("tier") or 9), -float(r.get("score") or 0))
        if key not in best or order < best[key][0]:
            best[key] = (order, r)
    ordered = sorted(best.items(), key=lambda kv: kv[1][0])
    keys = [k for k, _ in ordered]
    ranks = _competition_ranks(keys, [v[0] for _, v in ordered])
    return {"ranks": ranks, "best": {k: v[1] for k, v in best.items()},
            "n": len(keys)}


def docking_source(doc: dict, name: str, keep=None) -> dict:
    """Ranking from one docking result, or `{"unused": reason}`.

    Used only when the run says its gate passed (`usable`). A run that failed
    its gate has hits that were scored by a pocket that could not recover its
    own known ligands, so they carry no information.
    """
    if not doc.get("usable"):
        return {"name": name, "unused": doc.get("error") or
                "the docking validation gate did not pass"}
    ranks, hits = {}, {}
    # `keep` drops compounds before ranking, so a hit the reports do not show
    # does not take a rank. Ranks are then counted among the kept ones.
    listed = [h for h in doc.get("hits") or []
              if keep is None or keep(h.get("drug") or "")]
    for i, h in enumerate(listed, start=1):
        key = breadth_key(h.get("drug") or "")
        if key and key not in ranks:
            ranks[key] = i
            hits[key] = h
    if not ranks:
        return {"name": name, "unused": "the docking run listed no hits"}
    return {"name": name, "ranks": ranks, "hits": hits,
            "auc": doc.get("auc"), "ef1": doc.get("ef1"),
            "n_docked": doc.get("n_docked"),
            "conformation_verified": bool(doc.get("conformation_verified")),
            "pocket": (doc.get("pocket") or {}).get("pdb_id")
            or (doc.get("pocket") or {}).get("structure")}


# ------------------------------------------------------------------- fusion
def fuse(routes: dict | None, dockings: list) -> list:
    """Combined ranking: a list of rows, best first.

    Every row has `rank`, `score`, the per-source ranks, and the route and
    docking records behind them.
    """
    sources = []
    if routes and routes["ranks"]:
        sources.append(("routes", routes["ranks"]))
    seen: dict = {}
    for d in dockings:
        if "ranks" in d:
            # Two files can name the same target. Each is still its own list,
            # so the name is made unique rather than letting one overwrite
            # the other.
            label = f"docking: {d['name']}"
            seen[label] = seen.get(label, 0) + 1
            if seen[label] > 1:
                label += f" ({seen[label]})"
            sources.append((label, d["ranks"]))
    keys = {k for _, rk in sources for k in rk}
    rows = []
    for k in keys:
        per = {name: rk[k] for name, rk in sources if k in rk}
        score = sum(1.0 / (RRF_K + r) for r in per.values())
        row = {"key": k, "score": round(score, 6), "source_ranks": per,
               "n_sources": len(per)}
        if routes and k in routes["best"]:
            row["route"] = routes["best"][k]
        for d in dockings:
            if "hits" in d and k in d["hits"]:
                row.setdefault("docking", []).append(
                    {"source": d["name"], "hit": d["hits"][k],
                     "conformation_verified": d["conformation_verified"]})
        rows.append(row)
    # most sources first at an equal score would hide the score; sort on the
    # score, then on the better single rank, then on the name, so the order
    # is reproducible
    rows.sort(key=lambda r: (-r["score"], min(r["source_ranks"].values()),
                             r["key"]))
    rank = 0
    last = None
    for i, r in enumerate(rows, start=1):
        if r["score"] != last:
            rank, last = i, r["score"]
        r["rank"] = rank
    return rows


# ------------------------------------------------------------------- output
def _name(row: dict) -> str:
    if "route" in row:
        return row["route"]["drug"]
    hit = row["docking"][0]["hit"]
    return str(hit.get("drug") or row["key"]).title()


def flags(row: dict) -> list:
    """Short cautions that belong beside a drug in the list."""
    out = []
    r = row.get("route") or {}
    if r.get("breadth_penalty"):
        out.append("moves many genes")
    if r.get("mitotic_penalty"):
        out.append("disturbs mitosis")
    if str(r.get("burden") or "") in ("high", "prohibitive"):
        out.append(f"{r['burden']} clinical burden")
    # A drug that arrived by docking alone has no route record, so the
    # route's own cautions never reach it. Its class does.
    from . import systemic
    cls = systemic.classify(_name(row))
    if cls is not None and cls.burden >= 0.5:
        out.append(f"{cls.name}: acts across the body")
    for d in row.get("docking") or []:
        if not d["conformation_verified"]:
            out.append("docking: binding, not direction")
            break
    if row["n_sources"] == 1:
        out.append("one source only")
    return out


def cluster_rows(rows: list, threshold: float | None = None) -> list:
    """Cluster the fused rows, best first (see `clusters.group`)."""
    from . import clusters
    items = [{"name": _name(r), "row": r} for r in rows]
    return clusters.group(items, threshold=clusters.SIMILARITY
                          if threshold is None else threshold)


def _how(row: dict, gene: str, goal: str, roles: dict) -> str:
    from . import rationale
    if row.get("route"):
        return rationale.route_sentence(row["route"], gene, goal, roles)
    dock = (row.get("docking") or [{}])[0]
    return rationale.docking_sentence(dock.get("source") or "?", gene)


def table_rows(rows: list, *, clusters: list | None = None, gene: str = "?",
               goal: str = "up", roles: dict | None = None) -> list:
    """Flat dicts for the TSV and the report.

    With `clusters` only each cluster's leader is a row, and `similar` names
    the drugs folded into it.
    """
    from . import clusters as cl_mod
    folded = {}
    if clusters is not None:
        folded = {id(c["leader"]["row"]): c for c in clusters}
        rows = [c["leader"]["row"] for c in clusters]
    out = []
    for r in rows:
        route = r.get("route") or {}
        c = folded.get(id(r))
        dock = (r.get("docking") or [{}])[0]
        hit = dock.get("hit") or {}
        out.append({
            "rank": r["rank"], "drug": _name(r), "score": r["score"],
            "sources": "; ".join(f"{k} #{v}" for k, v in
                                 r["source_ranks"].items()),
            "route": (f"{route['category']} via {route['acts_on']} "
                      f"(tier {route['tier']}, {route['score']})"
                      if route else ""),
            "docking_score": hit.get("score"),
            "docking_cnn": hit.get("cnn_score"),
            "how": _how(r, gene, goal, roles or {}),
            "similar": (cl_mod.names(c, limit=8) if c and c["members"]
                        else ""),
            "n_similar": len(c["members"]) if c else 0,
            "cautions": "; ".join(flags(r))})
    return out


def markdown(report: dict, top: int = 40) -> str:
    L = [f"# Combined ranking: {report['gene']}", "",
         "One list from the routes and from docking, fused by reciprocal rank "
         f"(score = sum of 1 / ({RRF_K} + rank) over the lists a drug is in). "
         "The score orders drugs and is not a probability. The lists answer "
         "different questions, and the columns say which each drug answered.",
         ""]
    L += ["## Sources", ""]
    r = report["routes_source"]
    L.append(f"- Routes: {r['label']}, {r['n']} drugs.")
    for d in report["docking_sources"]:
        if "unused" in d:
            L.append(f"- Docking `{d['name']}`: **not used**, {d['unused']}.")
        else:
            L.append(f"- Docking `{d['name']}`: used. AUC {d['auc']}, EF1 "
                     f"{d['ef1']}, {d['n_docked']} compounds docked, "
                     f"{len(d['ranks'])} hits listed. Pocket "
                     f"{d['pocket'] or 'not recorded'}."
                     + ("" if d["conformation_verified"] else
                        " The conformation could not be read, so a hit is a "
                        "binder and not known to activate or block."))
    L += ["", "## Ranking", ""]
    if report.get("aim"):
        L += [report["aim"], ""]
    if report.get("folded"):
        L += [f"Near-identical drugs are folded into the best-ranked one "
              f"(structure within Tanimoto {report['similarity']} of it): "
              f"{report['n_folded']} drugs are folded into "
              f"{report['n_clusters_folding']} rows below.", ""]
    L += ["| Rank | Drug | Score | In lists | Best route | Docking score "
          "(CNN) | How it would work | Cautions |",
          "|---|---|---|---|---|---|---|---|"]
    for t in report["table"][:top]:
        dock = ("" if t["docking_score"] is None else
                f"{t['docking_score']} ({t['docking_cnn']})")
        more = f" (+{t['n_similar']} similar)" if t["n_similar"] else ""
        L.append(f"| {t['rank']} | {t['drug'].replace('|', '/')}{more} | "
                 f"{t['score']} | {t['sources'].replace('|', '/')} | "
                 f"{t['route'].replace('|', '/')} | {dock} | "
                 f"{t['how'].replace('|', '/')} | {t['cautions']} |")
    shown = [t for t in report["table"][:top] if t["n_similar"]]
    if shown:
        L += ["", "**Folded into the rows above.** Analogues can differ in "
              "potency, selectivity, tissue reach and safety, so these are "
              "neighbours to check, not duplicates.", ""]
        for t in shown:
            L.append(f"- {t['drug'].replace('|', '/')}: "
                     f"{t['similar'].replace('|', '/')}")
    L += ["", "## How to read this", "",
          "- A drug in two lists is high in two rankings. That is agreement, "
          "not extra evidence, and lists that share an input are not "
          "independent.",
          "- No list is weighted. A pocket that barely passed its gate counts "
          "the same as one that passed well; read the AUC and EF1 above.",
          "- Tied drugs share the best rank of their group.",
          "- A docking row says the compound is shaped to fit the pocket. "
          "It does not say it changes the gene, or in which direction.", ""]
    return "\n".join(L)


def build(*, results: dict | None = None, panel: dict | None = None,
          dockings: list | None = None, gene: str | None = None,
          fold: bool = True, similarity: float | None = None) -> dict:
    """The whole report from loaded files. `dockings` is `[(name, dict)]`.

    `fold` folds near-identical drugs into one row (see `clusters`), at
    `similarity` if given.
    """
    from . import clusters as cl_mod, fda as fda_mod, rationale
    fda_on = bool(((results or {}).get("fda_filter") or {}).get("applied")) or (
        results is None and bool(((panel or {}).get("query") or {})
                                 .get("fda_only")) and fda_mod.available())
    if results is not None:
        results = fda_mod.apply(results)
    gene = gene or ((results or {}).get("query") or {}).get("gene") or (
        ((panel or {}).get("query") or {}).get("genes") or ["?"])[0]
    if panel is not None:
        routes = routes_source(panel.get("surviving") or [], gene)
        routes["label"] = (f"the panel's surviving routes for {gene} after "
                           f"its funnel")
    elif results is not None:
        pool = [r for blk in results.get("results") or []
                for r in (blk.get("all_routes") or {}).get("routes") or []]
        routes = routes_source(pool)
        routes["label"] = f"every route in the run for {gene}"
    else:
        routes = {"ranks": {}, "best": {}, "n": 0,
                  "label": "none given"}
    docs = [docking_source(d, n, fda_mod.is_approved if fda_on else None)
            for n, d in (dockings or [])]
    rows = fuse(routes, docs)
    goal = (((results or {}).get("query") or {}).get("goal")
            or ((panel or {}).get("query") or {}).get("goal") or "up")
    roles = {t["symbol"]: t.get("role")
             for blk in (results or {}).get("results") or []
             for t in blk.get("tfs") or []}
    threshold = cl_mod.SIMILARITY if similarity is None else similarity
    clusters = cluster_rows(rows, threshold) if fold else None
    out = {"gene": gene, "goal": goal,
           "routes_source": {"label": routes["label"], "n": routes["n"]},
           "docking_sources": [{k: v for k, v in d.items() if k != "hits"}
                               for d in docs],
           "n_ranked": len(rows),
           "aim": rationale.aim(gene, goal,
                                (results or {}).get("requirement"),
                                ((results or {}).get("query") or {}).get(
                                    "condition"),
                                ((results or {}).get("query") or {}).get(
                                    "scenario")),
           "similarity": threshold if fold else None,
           "table": table_rows(rows, clusters=clusters, gene=gene, goal=goal,
                               roles=roles)}
    folded = [c for c in (clusters or []) if c["members"]]
    out["folded"] = [{"leader": c["leader"]["name"],
                      "leader_rank": c["leader"]["row"]["rank"],
                      "members": [{"name": m["name"],
                                   "rank": m["row"]["rank"],
                                   "similarity": c["similarity"][m["name"]]}
                                  for m in c["members"]]}
                     for c in folded]
    out["n_folded"] = sum(len(c["members"]) for c in folded)
    out["n_clusters_folding"] = len(folded)
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.combined",
        description=("One ranked list from docking and the rest of the "
                     "analysis, fused by reciprocal rank."))
    p.add_argument("--results", default=None, help="a saved results.json")
    p.add_argument("--panel", default=None, help="a saved panel.json")
    p.add_argument("--docking", action="append", default=[],
                   help="a docking result file (repeat for several)")
    p.add_argument("--gene", default=None)
    p.add_argument("--similarity", type=float, default=None,
                   help="Tanimoto similarity (0 to 1) at which a drug is "
                        "folded into a better-ranked one; 0.75 by default")
    p.add_argument("--no-fold", action="store_true",
                   help="list every drug, folding none")
    p.add_argument("--top", type=int, default=40)
    p.add_argument("--out", default="results/combined")
    args = p.parse_args(argv)
    if args.top < 1:
        p.error("--top must be at least 1")
    if args.similarity is not None and not 0.0 < args.similarity <= 1.0:
        p.error("--similarity must be above 0 and at most 1")
    if not (args.results or args.panel or args.docking):
        p.error("give at least one of --results, --panel or --docking")

    def load(path):
        path = Path(path)
        if not path.exists():
            p.error(f"{path} not found")
        return json.loads(path.read_text(encoding="utf-8"))

    results = load(args.results) if args.results else None
    panel = load(args.panel) if args.panel else None
    dockings = [(load(d).get("target") or Path(d).stem, load(d))
                for d in args.docking]
    report = build(results=results, panel=panel, dockings=dockings,
                   gene=args.gene, fold=not args.no_fold,
                   similarity=args.similarity)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "combined.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8")
    text = markdown(report, args.top)
    (out / "combined.md").write_text(text, encoding="utf-8")
    cols = ["rank", "drug", "score", "sources", "route", "docking_score",
            "docking_cnn", "how", "similar", "n_similar", "cautions"]
    with open(out / "combined.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(report["table"])
    print(text)
    print(f"Wrote {out / 'combined.md'}, combined.json and combined.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
