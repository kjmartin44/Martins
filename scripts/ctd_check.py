"""How much of CTD's expression record is about the gene's own transcript.

    python scripts/ctd_check.py <results.json files> [--out results/retrospective_fate]

For each saved run this reads the CTD rows for the gene (from the request
cache, so it sends nothing that was not asked before) and sorts every row that
carries an expression action by what its sentence says:

  mRNA of this gene       the chemical is the subject and the transcript is named
  protein of this gene    expression of the protein, which can move without the
                          transcript
  modulated               the chemical changes another's effect on the gene
  indirect                the effect comes through another gene
  another gene            the expression named is not this gene's
  unparsed                the gene is named in a shape that is not read

It then scores the routes in each saved run against the CTD key built the old
way (every `expression` action in the list) and the new way (the sentence),
using the same tests as `retrospective.py`. Nothing is validated by this. It
says how far the answer key itself can be trusted, and whether the routes look
better or worse against the cleaner one.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from generoutes import discovery as disc  # noqa: E402
from generoutes import localdrugs  # noqa: E402
from generoutes import retrospective as rt  # noqa: E402
from generoutes.http import Client  # noqa: E402
from generoutes.signatures import breadth_key  # noqa: E402

LABELS = ["mRNA of this gene", "protein of this gene", "modulated",
          "indirect", "another gene", "unparsed"]


def rows_for(client, gene):
    res = client.get(disc.CTD_BATCH, {
        "inputType": "gene", "inputTerms": gene.upper(), "report": "cgixns",
        "format": "json", "actionTypes": "ANY"})
    return res if isinstance(res, list) else res.get("data", []) or []


def classify_row(row) -> str:
    gene = row.get("GeneSymbol") or row.get("Input") or ""
    eff = disc.expression_effects(row.get("ChemicalName"),
                                  row.get("Interaction"), gene)
    if not eff:
        return "another gene"
    roles = {r for _, _, r in eff}
    if "direct" in roles:
        forms = {f for _, f, r in eff if r == "direct"}
        if forms & {"mRNA", None}:
            return "mRNA of this gene"
        return "protein of this gene"
    if "modulated" in roles:
        return "modulated"
    if "indirect" in roles:
        return "indirect"
    return "unparsed"


def keyset(links: dict, goal: str):
    ks = rt.KeySet("ctd")
    for chem, actions in links.items():
        key = breadth_key(chem)
        if not key:
            continue
        d = disc.link_direction(actions)
        if d == goal:
            ks.wanted.add(key)
        elif d in ("up", "down"):
            ks.opposite.add(key)
        else:
            ks.unclear.add(key)
    ks.n_records = len(links)
    return ks


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("results", nargs="+")
    p.add_argument("--library", default=None)
    p.add_argument("--permutations", type=int, default=300)
    p.add_argument("--out", default="results/retrospective_fate")
    args = p.parse_args(argv)

    client = Client()
    rows_lib = localdrugs.library(args.library)
    library = {breadth_key(r["name"]) for r in rows_lib} - {""}
    runs = [r for r in (rt.load_run(x) for x in args.results) if r]
    runs, _ = rt.latest_per_gene_goal(runs)

    kinds = collections.Counter()
    per_gene = []
    old_runs, new_runs = [], []
    flips = collections.Counter()
    for run in runs:
        try:
            rows = rows_for(client, run.gene)
        except Exception as exc:                              # noqa: BLE001
            print(f"{run.gene}: CTD not readable ({str(exc)[:60]})")
            continue
        expr = [r for r in rows
                if "expression" in str(r.get("InteractionActions", "")).lower()]
        c = collections.Counter(classify_row(r) for r in expr)
        kinds.update(c)
        old = {}
        for r in expr:
            chem = str(r.get("ChemicalName") or "").strip().lower()
            old.setdefault(chem, []).extend(
                a for a in str(r["InteractionActions"]).split("|")
                if "expression" in a.lower())
        new = disc.strict_expression_links(rows)
        ko, kn = keyset(old, run.goal), keyset(new, run.goal)
        for k in (ko.wanted | ko.opposite) & (kn.wanted | kn.opposite):
            if (k in ko.wanted) != (k in kn.wanted):
                flips[run.gene] += 1
        per_gene.append((run.gene, run.goal, len(expr), len(old), len(new),
                         len(ko.wanted & library), len(kn.wanted & library)))
        for bucket, ks in ((old_runs, ko), (new_runs, kn)):
            bucket.append(rt.Run(run.gene, run.goal, run.path, run.routes,
                                 {"ctd": ks}, run.skipped))

    out = {}
    for label, bucket in (("old", old_runs), ("new", new_runs)):
        fame = rt.fame_counts(bucket)
        per_run = [rt.evaluate_run(r, library, permutations=args.permutations,
                                   fame=fame) for r in bucket]
        out[label] = rt.summarise(per_run)["ctd"]["predicted"]

    total = sum(kinds.values())
    L = ["# How much of CTD's expression record concerns the gene's transcript",
         "",
         f"{len(per_gene)} genes, {total} CTD rows with an expression action.",
         "", "| What the sentence says | Rows | Share |", "|---|---|---|"]
    for k in LABELS:
        L.append(f"| {k} | {kinds[k]} | {kinds[k] / total:.0%} |")
    L += ["", "Only the first row is read as an effect of the chemical on the "
          "gene's transcript. A sentence that names no form is read as the "
          "transcript and counted there.", "",
          "| Gene | Rows | Chemicals, old | Chemicals, new | Approved wanted, "
          "old | Approved wanted, new |", "|---|---|---|---|---|---|"]
    for g, goal, n, a, b, wo, wn in per_gene:
        L.append(f"| {g} ({goal}) | {n} | {a} | {b} | {wo} | {wn} |")
    L += ["", "### The routes scored against each key", "",
          "| Key | Approved drugs wanted | With a route | Mean AUC, right "
          "direction | Mean AUC, reported or not | Matched to fame |",
          "|---|---|---|---|---|---|"]
    for label in ("old", "new"):
        s = out[label]
        ro, rd = s["reported_or_not"], s["right_direction"]
        L.append(f"| {label} | {ro['n_positives']} | "
                 f"{ro['n_positives_with_route']} | {rd['mean_auc']} | "
                 f"{ro['mean_auc']} | {ro.get('mean_matched_auc')} |")
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "ctd_check.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
