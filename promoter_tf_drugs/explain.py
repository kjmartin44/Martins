"""Everything the pipeline knows about one compound, in one place.

A run scatters evidence about a single drug across eight sections: the route
catalogue, the grade table with its components, the outcome chain, the
mechanism comparison, the ligand axis, the LINCS signature, the discovery
comparison against the standard of care, and the screen cherry-pick list. To
judge one candidate you currently have to find it in each of them.

    python -m promoter_tf_drugs.explain results/PAH_up/results.json calcitriol

Matching is case-insensitive and substring-based, so a partial name works.
Where several compounds match, the ambiguity is reported rather than the first
one being picked silently.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _matches(name: str, query: str) -> bool:
    n, q = (name or "").lower(), query.lower()
    return q in n or n in q


def collect(result: dict, query: str) -> dict:
    """Pull every mention of `query` out of a results.json."""
    found: dict = {"query": query, "names": set(), "blocks": []}

    # direct_therapies is stored once at the top level, not per promoter block
    top_direct = [t for t in (result.get("direct_therapies") or [])
                  if _matches(t.get("drug", ""), query)]

    for block in result.get("results", []):
        gene = (block.get("promoter") or {}).get("gene_symbol")
        entry: dict = {"gene": gene}
        if top_direct:
            entry["direct"] = list(top_direct)
            found["names"].update(t["drug"] for t in top_direct)

        for r in (block.get("all_routes") or {}).get("routes", []):
            if _matches(r["drug"], query):
                found["names"].add(r["drug"])
                entry.setdefault("routes", []).append(r)

        for g in block.get("drug_grades") or []:
            if _matches(g["drug_name"], query):
                found["names"].add(g["drug_name"])
                entry.setdefault("grades", []).append(g)

        for o in block.get("outcome_estimates") or []:
            if _matches(o["drug"], query):
                entry.setdefault("outcomes", []).append(o)

        for sym, cmp in (block.get("mechanisms") or {}).items():
            for h in cmp.get("hits", []):
                if _matches(h["drug"], query):
                    entry.setdefault("mechanisms", []).append({**h, "tf": sym})

        for sym, a in (block.get("ligand_axes") or {}).items():
            for r in a.get("routes", []):
                if _matches(r["drug"], query):
                    entry.setdefault("ligand", []).append({**r, "tf": sym})

        for sym, a in (block.get("stabilisation") or {}).items():
            for st in a.get("stabilisers", []):
                if _matches(st["drug"], query):
                    entry.setdefault("stabilisers", []).append({**st, "tf": sym})

        for t in (block.get("direct_therapies") or []):
            if _matches(t["drug"], query):
                entry.setdefault("direct", []).append(t)
                found["names"].add(t["drug"])

        for h in (block.get("lincs") or {}).get("hits", []):
            if _matches(h["drug"], query):
                entry.setdefault("lincs", []).append(h)

        imp = block.get("improvements") or {}
        for c in imp.get("all_candidates") or []:
            if _matches(c["name"], query):
                entry.setdefault("comparison", []).append(c)
        base = imp.get("baseline")
        if base and _matches(base["name"], query):
            entry["is_baseline"] = True
        entry["bar"] = imp.get("baseline_envelope")

        for sym, hits in (block.get("upstream") or {}).items():
            for h in hits:
                for d in h.get("drugs", []):
                    if _matches(d.get("name", ""), query):
                        entry.setdefault("upstream", []).append(
                            {**d, "tf": sym, "node": h["node"],
                             "path": h["best_path"]["description"]})

        if len(entry) > 2:
            found["blocks"].append(entry)

    found["names"] = sorted(found["names"])
    return found


def render(found: dict) -> str:
    q = found["query"]
    names = found["names"]
    L: list[str] = []

    if not found["blocks"]:
        return (f"# {q}\n\nNot found in this results file. Check the spelling, "
                f"or list what is present with:\n\n"
                f"    python -m promoter_tf_drugs.shortlist <results.json> "
                f"--preset mechanism --top 100\n")

    if len(names) > 1:
        L.append(f"> **{len(names)} compounds match '{q}'**: "
                 + ", ".join(names) + ". Everything below is pooled across "
                 "them — narrow the query if that is not what you want.")
        L.append("")

    L.append(f"# {names[0] if names else q}")
    L.append("")

    for b in found["blocks"]:
        L.append(f"## For {b.get('gene')}")
        L.append("")

        for r in b.get("routes", []):
            a = r["axes"]
            L.append(f"**Route**: {r['category']} — {r['mechanism']}")
            L.append("")
            L.append(f"- Acts on **{r['acts_on']}**, {r['directness']} step(s) "
                     f"from the gene")
            L.append(f"- Evidence tier {r['tier']}: {r['tier_label']}")
            L.append(f"- Strength {a['strength']} · specificity "
                     f"{a['specificity']} · efficiency {a['efficiency']} "
                     f"· score {r['score']}")
            for n in r.get("notes", [])[:3]:
                L.append(f"- {n}")
            L.append("")

        for d in b.get("direct", []):
            L.append(f"**Direct agent** ({d['modality']}, {d['approval_year']}): "
                     f"{d['mechanism']}")
            L.append("")
            L.append(f"- Raises protein: "
                     f"{'yes' if d['raises_protein'] else '**activity only**'}")
            L.append("")

        for g in b.get("grades", []):
            c = g["components"]
            L.append(f"**Grade {g['grade']}** ({g['combined']}/100) via "
                     f"{g['tf_symbol']}")
            L.append("")
            L.append("| Component | Value |")
            L.append("|---|---|")
            for k, v in c.items():
                L.append(f"| {k.replace('_', ' ')} | {v} |")
            L.append(f"| TF confidence | {g['tf_confidence']} |")
            L.append(f"| drug quality | {g['drug_quality']} |")
            L.append("")
            ev = g.get("evidence") or {}
            bits = []
            if ev.get("pchembl_median"):
                bits.append(f"median pChEMBL {ev['pchembl_median']} "
                            f"(n={ev.get('pchembl_n')})")
            if ev.get("n_targets"):
                bits.append(f"{ev['n_targets']} annotated targets")
            if ev.get("regulon_size"):
                bits.append(f"TF controls ~{ev['regulon_size']} genes")
            if ev.get("withdrawn"):
                bits.append("**WITHDRAWN**")
            if ev.get("black_box"):
                bits.append("**black-box warning**")
            if bits:
                L.append("Evidence: " + "; ".join(bits))
                L.append("")
            for n in g.get("notes", [])[:4]:
                L.append(f"- {n}")
            L.append("")

        for m in b.get("mechanisms", []):
            p = m.get("profile") or {}
            L.append(f"**Mechanism class**: {m['mechanism']} at {m['tf']} "
                     f"({m['action_type']})")
            L.append("")
            if p:
                L.append(f"- Max effect {p.get('max_effect')}; subtype "
                         f"selectivity {p.get('subtype_selectivity')}")
                L.append(f"- Preserves signalling dynamics: "
                         f"{'yes' if p.get('preserves_dynamics') else 'no'}")
                L.append(f"- Needs endogenous ligand: "
                         f"{'**YES**' if p.get('needs_endogenous_ligand') else 'no'}")
                if p.get("caveat"):
                    L.append(f"- {p['caveat']}")
            L.append("")

        for r in b.get("ligand", []):
            L.append(f"**Ligand axis** at {r['tf']}: {r['strategy']} "
                     f"({r['ligand']}) via {r['acts_on']}")
            L.append("")
            L.append(f"- Moves the receptor **{r['direction']}**"
                     + (" · tissue-selective" if r["tissue_selective"] else ""))
            for c in r.get("caveats", [])[:3]:
                L.append(f"- {c}")
            L.append("")

        for st in b.get("stabilisers", []):
            L.append(f"**Stabilisation** at {st['tf']}: {st['mechanism']}, "
                     f"binds {st['target_node']}, raises {st['raises']} "
                     f"({st['status']})")
            L.append("")

        for u in b.get("upstream", []):
            L.append(f"**Upstream route** to {u['tf']}: {u['path']}")
            L.append("")

        for h in b.get("lincs", []):
            kind = ("measured" if h.get("is_landmark") is True
                    else "INFERRED, not measured")
            L.append(f"**Measured signature**: score {h['score']} in "
                     f"{h.get('cell_line') or 'unknown line'} — {kind}")
            L.append("")

        for o in b.get("outcomes", []):
            L.append(f"**Predicted outcome**: P(gene moves) = "
                     f"**{o['p_point']}** ({o['p_low']}–{o['p_high']}) — "
                     f"{o['verdict']}")
            L.append("")
            L.append(f"- engages {o['p_engages']} × direction "
                     f"{o['p_direction']} × promoter {o['p_activates']} × "
                     f"responsive {o['p_responsive']}")
            L.append(f"- **Limiting step**: {o['limiting_step']}")
            if o.get("next_experiment"):
                L.append(f"- Next experiment: {o['next_experiment']}")
            if o.get("tissue_dependent"):
                L.append("- ⚠ tissue-dependent direction (SERM regime)")
            L.append("")

        for c in b.get("comparison", []):
            a = c["axes"]
            L.append(f"**Against the standard of care**: {c['verdict']}")
            L.append("")
            L.append(f"- Axes: strength {a['strength']} · specificity "
                     f"{a['specificity']} · efficiency {a['efficiency']} "
                     f"({a['evidence']})")
            if b.get("bar"):
                bar = b["bar"]
                L.append(f"- The bar: {bar['strength']} · {bar['specificity']} "
                         f"· {bar['efficiency']} ({bar['evidence']})")
            L.append(f"- Previously linked to this gene: "
                     f"{'yes — ' + (c.get('link_source') or 'reported') if c['previously_linked'] else '**no**'}")
            if c.get("wins_on"):
                L.append(f"- Wins on: {', '.join(c['wins_on'])}")
            if c.get("loses_on"):
                L.append(f"- Worse on: {', '.join(c['loses_on'])}")
            for n in c.get("notes", [])[:3]:
                L.append(f"- {n}")
            L.append("")

        if b.get("is_baseline"):
            L.append("*This is the incumbent — the option everything else is "
                     "measured against.*")
            L.append("")

    return "\n".join(L)


def check_schema(data: dict, tool: str) -> str | None:
    """Warn when a saved result was written by a newer package.

    Returns the warning text, and prints it to stderr. None when the
    file is one this version understands - including a file with no
    schema at all, which predates the field and is read quietly rather
    than warned about every time.
    """
    import sys

    from .provenance import RESULT_SCHEMA

    found = (data.get("provenance") or {}).get("schema")
    if found is None or found <= RESULT_SCHEMA:
        return None
    version = (data.get("provenance") or {}).get("package_version")
    msg = (f"{tool}: this result uses schema {found} and this package "
           f"understands {RESULT_SCHEMA}"
           + (f" (written by version {version})" if version else "")
           + ". Reading it anyway; fields added since are ignored, and "
             f"any whose meaning changed are read with the older one. "
             f"Upgrade the package to read it as it was written.")
    print(msg, file=sys.stderr)
    return msg


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="promoter-tf-drugs.explain",
        description="Show everything a run found about one compound.")
    p.add_argument("results", help="results.json from a previous run")
    p.add_argument("drug", nargs="+", help="compound name (partial is fine)")
    p.add_argument("--out", help="write the write-up here")
    args = p.parse_args(argv)

    data = json.loads(Path(args.results).read_text(encoding="utf-8"))
    check_schema(data, "explain")
    text = render(collect(data, " ".join(args.drug)))
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
