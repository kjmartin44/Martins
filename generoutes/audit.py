"""Are the hand-curated tables still true?

The problem
-----------
Several tables here are hand-entered: the systemic drug classes, the receptor
families, the curated stabilisers, the P-gp substrates. Each says in its own
docstring that it is incomplete and will go out of date — which is honest,
and useless on its own, because nothing tells you *when* it has drifted.

Approvals happen constantly, drugs get withdrawn, and a name that resolved
two years ago may now be a synonym of something else. A table that is quietly
wrong is worse than one known to be incomplete, because the report presents
both with the same confidence.

Two levels of audit
-------------------
**Offline, and run as a test.** Internal consistency: no drug in two systemic
classes, no receptor in two families, no empty or too-thin fields, and every
name normalising to itself under `congeners.salt_free_name` so the tables and
the salt-collapsing logic agree about what a drug is called.

**Online, opt-in.** Resolve every curated drug name against ChEMBL and report
the ones that no longer resolve, plus any that ChEMBL now flags withdrawn.
One cached request per name. It also compares the nuclear receptor list with
the HGNC gene group, and the docking residue tables with the PDB chemical
component dictionary. Those are the tables a public list defines, so for
them a gap can be named. `docs/curated-tables.md` says which tables have no
such list.

    python -m generoutes.audit             # offline only
    python -m generoutes.audit --online    # also resolve in ChEMBL

Neither level can tell you what is **missing** from a curated table — only a
human reading the current literature can. What this catches is drift in what
is already there, which is the failure that would otherwise be silent.
"""

from __future__ import annotations

import csv
import io
import logging
import urllib.parse

from .congeners import salt_free_name
from . import provenance as _prov

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Group 71 is HGNC's nuclear hormone receptors; type=branch takes every subgroup.
HGNC_NR_GROUP = ("https://www.genenames.org/cgi-bin/genegroup/download"
                 "?id=71&type=branch")
HGNC_SEARCH = "https://rest.genenames.org/search/"
# Entries in the gene lists that are not approved symbols, and stay.
_NOT_SYMBOLS = {
    "SNORD116": "the name of a cluster of snoRNA genes",
    "NESP": "a transcript of GNAS, which is listed beside it",
}
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
_CHEMCOMP_QUERY = ("query($ids:[String!]!){chem_comps(comp_ids:$ids){"
                   "chem_comp{id type}}}")
# Codes the dictionary has retired, which older structures still write: WAT is
# the old water name (it now holds only HOH), and the rest are disaccharides
# that were split into their sugars.
_LEGACY_CODES = frozenset({"WAT", "CBK", "LAK", "MAL", "SUC", "TRE"})


# ---------------------------------------------------------------- offline
def check_internal() -> dict:
    """Consistency checks needing no network. Returns {check: [problems]}."""
    from . import systemic, paralogs, stabilizers, delivery

    problems: dict = {}

    # A drug in two systemic classes makes `classify` order-dependent, so the
    # burden it reports would depend on table declaration order.
    seen: dict = {}
    dupes = []
    for c in systemic.CLASSES:
        for m in c.members:
            if m in seen:
                dupes.append(f"{m}: {seen[m]} and {c.name}")
            seen[m] = c.name
    problems["systemic: drug in two classes"] = dupes

    # Same hazard for `family_of`.
    seen = {}
    fam_dupes = []
    for f in paralogs.FAMILIES:
        for m in f.members:
            if m in seen:
                fam_dupes.append(f"{m}: {seen[m]} and {f.name}")
            seen[m] = f.name
    problems["paralogs: receptor in two families"] = fam_dupes

    # Names must survive salt-stripping unchanged, or `classify` will miss
    # them: the lookup is done on the salt-free form.
    unstable = [m for c in systemic.CLASSES for m in c.members
                if salt_free_name(m) != m]
    problems["systemic: name changes under salt_free_name"] = unstable

    # A generic warning is not actionable; the point of `why` is the specific
    # systemic phenotype.
    problems["systemic: why too thin to act on"] = [
        c.name for c in systemic.CLASSES if len(c.why) < 40]
    problems["paralogs: note too thin"] = [
        f.name for f in paralogs.FAMILIES if len(f.note) < 30]

    problems["systemic: burden out of range"] = [
        c.name for c in systemic.CLASSES if not 0.0 < c.burden <= 1.0]

    # A stabiliser with no target node cannot be acted on.
    problems["stabilizers: no target node"] = [
        f"{tf}/{s.drug}" for tf, lst in stabilizers.CURATED.items()
        for s in lst if not s.target_node]

    # P-gp names are matched lowercase.
    problems["delivery: PGP name not lowercased"] = [
        n for n in delivery.PGP_SUBSTRATES if n != n.lower()]

    return problems


def summarise(problems: dict) -> str:
    bad = {k: v for k, v in problems.items() if v}
    if not bad:
        return ("All curated tables are internally consistent. This does NOT "
                "mean they are complete or current — only a human reading the "
                "literature can say what is missing.")
    lines = [f"{sum(len(v) for v in bad.values())} problem(s) found:"]
    for check, items in sorted(bad.items()):
        lines.append(f"  {check}:")
        for i in items[:8]:
            lines.append(f"    - {i}")
        if len(items) > 8:
            lines.append(f"    ... {len(items) - 8} more")
    return "\n".join(lines)


# ----------------------------------------------------------------- online
def resolve_names(client, names, *, limit: int | None = None) -> dict:
    """Which curated drug names still resolve in ChEMBL?

    A name that no longer resolves is the clearest signal a table has
    drifted: renamed, merged, or wrong from the start. Resolution proves
    existence, never correctness.
    """
    names = sorted(set(names))
    if limit:
        names = names[:limit]
    found, missing, withdrawn, errors = {}, [], [], 0

    for name in names:
        try:
            res = client.get(f"{CHEMBL}/molecule/search",
                             {"q": name, "format": "json", "limit": 1})
        except Exception as exc:                            # noqa: BLE001
            errors += 1
            log.warning("ChEMBL lookup failed for %r: %s", name,
                        str(exc)[:100])
            continue
        mols = (res or {}).get("molecules") or []
        if not mols:
            missing.append(name)
            continue
        m = mols[0]
        found[name] = m.get("molecule_chembl_id")
        if m.get("withdrawn_flag"):
            withdrawn.append(name)

    return {"checked": len(names), "resolved": found, "unresolved": missing,
            "withdrawn_in_chembl": withdrawn, "lookup_errors": errors}


def all_curated_names() -> dict:
    """Every hand-entered drug name in the package, by the table it lives in.

    A "full" audit means all of them. Auditing only the systemic classes
    would leave the direct-therapy list — the one the report prints FIRST and
    calls the better starting point — unchecked.
    """
    from . import (catalogue, systemic, stabilizers, delivery, ligands, direct,
                   indications, promoterdna, retrospective, standard_of_care)

    out: dict = {}
    out["systemic.CLASSES"] = sorted(
        {m for c in systemic.CLASSES for m in c.members})
    out["systemic.MITOTIC_DISRUPTORS"] = sorted(systemic.MITOTIC_DISRUPTORS)
    out["systemic.GENOTOXIC"] = sorted(systemic.GENOTOXIC)
    out["indications.CURATED"] = sorted(indications.CURATED)
    out["standard_of_care.CONDITIONS"] = sorted(
        {t.generic for c in standard_of_care.CONDITIONS
         for t in c.treatments})
    out["promoterdna.AGENTS"] = sorted({a.name for a in promoterdna.AGENTS})
    out["catalogue.CHROMATIN_DRUGS"] = sorted(
        {row[0] for row in catalogue.CHROMATIN_DRUGS})
    out["retrospective.POSITIVE_CONTROLS"] = sorted(
        {d for v in retrospective.POSITIVE_CONTROLS.values() for d in v[1]})
    out["stabilizers.CURATED"] = sorted(
        {s.drug for lst in stabilizers.CURATED.values() for s in lst})
    out["stabilizers.EXEMPLARS"] = sorted(
        {s.drug for s in stabilizers.EXEMPLARS})
    out["delivery.PGP_SUBSTRATES"] = sorted(delivery.PGP_SUBSTRATES)
    try:
        out["ligands.ROUTES"] = sorted(
            {r.drug for lst in ligands.ROUTES.values() for r in lst})
    except Exception as exc:                                # noqa: BLE001
        log.warning("audit could not read ligands.ROUTES (%s: %s); its "
                    "names are not checked", type(exc).__name__,
                    str(exc)[:100])
    try:
        rows = direct.CURATED
        vals = rows.values() if isinstance(rows, dict) else rows
        names = set()
        for v in vals:
            for item in (v if isinstance(v, (list, tuple)) else [v]):
                d = getattr(item, "drug", None) or (
                    item.get("drug") if isinstance(item, dict) else None)
                if d:
                    names.add(d)
        if names:
            out["direct.CURATED"] = sorted(names)
    except Exception as exc:                                # noqa: BLE001
        log.warning("audit could not read direct.CURATED (%s: %s); its "
                    "names are not checked", type(exc).__name__,
                    str(exc)[:100])
    return {k: v for k, v in out.items() if v}


def check_withdrawn_agreement(client, *, limit: int | None = None) -> dict:
    """Does our own withdrawn flag agree with ChEMBL's?

    Stronger than asking whether a name resolves. `direct.CURATED` carries a
    `withdrawn` field and the report prints that table FIRST, calling a
    direct agent "almost always the better starting point". An entry we mark
    available that ChEMBL now marks withdrawn is the single most consequential
    way these tables can drift — it would put a withdrawn drug at the top of
    the report with no flag on it.

    Disagreement in the other direction matters less but is still reported:
    approvals are regional and ChEMBL is not the regulator, so we may be
    right and it may be stale.
    """
    from . import direct

    rows = direct.CURATED
    vals = list(rows.values()) if isinstance(rows, dict) else list(rows)
    entries = [it for v in vals
               for it in (v if isinstance(v, (list, tuple)) else [v])
               if getattr(it, "drug", None)]
    if limit:
        entries = entries[:limit]

    we_say_ok_chembl_says_withdrawn = []
    we_say_withdrawn_chembl_says_ok = []
    unchecked = []

    for e in entries:
        try:
            res = client.get(f"{CHEMBL}/molecule/search",
                             {"q": e.drug, "format": "json", "limit": 1})
        except Exception:                                   # noqa: BLE001
            unchecked.append(e.drug)
            continue
        mols = (res or {}).get("molecules") or []
        if not mols:
            unchecked.append(e.drug)
            continue
        chembl_withdrawn = bool(mols[0].get("withdrawn_flag"))
        ours = bool(getattr(e, "withdrawn", False))
        if chembl_withdrawn and not ours:
            we_say_ok_chembl_says_withdrawn.append(e.drug)
        elif ours and not chembl_withdrawn:
            we_say_withdrawn_chembl_says_ok.append(e.drug)

    return {
        "checked": len(entries),
        "we_say_available_chembl_says_withdrawn":
            we_say_ok_chembl_says_withdrawn,
        "we_say_withdrawn_chembl_says_available":
            we_say_withdrawn_chembl_says_ok,
        "unchecked": unchecked,
    }


def all_curated_symbols() -> dict:
    """Every gene symbol in the hand-entered gene lists, by table."""
    from . import (annotate, cellcycle, direct, ligands, paralogs,
                   proliferation, reactivation, regulon, retrospective,
                   signatures, stabilizers, standard_of_care, systemic)

    tables = {
        "annotate.NUCLEAR_RECEPTORS": annotate.NUCLEAR_RECEPTORS,
        "cellcycle.MITOTIC_GENES": cellcycle.MITOTIC_GENES,
        "cellcycle.PROFILES": cellcycle.PROFILES,
        "proliferation.MARKERS": proliferation.MARKERS,
        "reactivation.IMPRINTED": reactivation.IMPRINTED,
        "regulon.GLOBAL_REGULATORS": regulon.GLOBAL_REGULATORS,
        "regulon.BROAD_REGULATORS": regulon.BROAD_REGULATORS,
        "regulon.MODERATE_REGULATORS": regulon.MODERATE_REGULATORS,
        "signatures.KNOWN_LANDMARKS": signatures.KNOWN_LANDMARKS,
        "systemic.ONCOGENES": systemic.ONCOGENES,
        "systemic.TUMOUR_SUPPRESSORS": systemic.TUMOUR_SUPPRESSORS,
        "paralogs.FAMILIES": [m for f in paralogs.FAMILIES
                              for m in f.members],
        "stabilizers.CURATED": stabilizers.CURATED,
        "stabilizers.DEGRON": stabilizers.DEGRON,
        "stabilizers.DUBS": list(stabilizers.DUBS) + [
            d for v in stabilizers.DUBS.values() for d in v["dubs"]],
        "direct.CURATED": direct.CURATED,
        "ligands.ROUTES": ligands.ROUTES,
        "retrospective.POSITIVE_CONTROLS": retrospective.POSITIVE_CONTROLS,
        "standard_of_care.CONDITIONS": (
            [g for c in standard_of_care.CONDITIONS for g in c.genes]
            + [t.target_gene for c in standard_of_care.CONDITIONS
               for t in c.treatments]),
    }
    return {k: {str(s).strip() for s in v if s} for k, v in tables.items()}


def check_gene_symbols(client, tables: dict, *, chunk: int = 40) -> dict:
    """Which curated gene symbols are not HGNC approved symbols?

    A symbol goes stale when HGNC renames the gene (SLC22A18 is now SLC67A1),
    and an alias can name a different gene (TAZ is an alias of WWTR1 and the
    old symbol of TAFAZZIN). A lookup by the stale symbol then finds nothing,
    or the wrong gene. Returns `{table: [symbol, ...]}` for the tables with any.
    """
    every = sorted(set().union(*tables.values()))
    approved = set()
    for i in range(0, len(every), chunk):
        query = "symbol:(" + " OR ".join(every[i:i + chunk]) + ")"
        res = client.get(HGNC_SEARCH + urllib.parse.quote(query, safe=":()"),
                         headers={"Accept": "application/json"})
        docs = ((res or {}).get("response") or {}).get("docs") or []
        approved |= {d["symbol"] for d in docs}
    bad = {t: sorted(s - approved - set(_NOT_SYMBOLS))
           for t, s in tables.items()}
    return {t: v for t, v in bad.items() if v}


def check_local_references() -> dict:
    """The checks that need a downloaded table. Each is None when it is absent.

    `landmarks`: genes in `signatures.KNOWN_LANDMARKS` that the CLUE gene
    table does not list as landmarks. `orangebook_routes`: routes the Book
    uses that neither route table classifies (their drugs stay "unclear"),
    and table entries the Book no longer uses (older releases did).
    """
    from . import orangebook, signatures

    out = {"landmarks": None, "orangebook_routes": None}
    if signatures.load_geneinfo():
        out["landmarks"] = sorted(
            g for g in signatures.KNOWN_LANDMARKS
            if signatures.feature_space(g) != "landmark")
    vocab = orangebook.route_vocabulary()
    if vocab:
        known = orangebook.SYSTEMIC_ROUTES | orangebook.LOCAL_ONLY_ROUTES
        out["orangebook_routes"] = {
            "unclassified": sorted(vocab - known),
            "not_in_this_release": sorted(known - vocab)}
    return out


def check_nuclear_receptors(text: str) -> dict:
    """Compare `annotate.NUCLEAR_RECEPTORS` with the HGNC receptor group.

    `text` is the group's TSV download. Only protein-coding rows count: the
    group also lists the pseudogene NR1H5P, which is not a receptor anyone can
    dose. An empty reference raises, because a portal page that parses to
    nothing would otherwise read as every receptor in our list being extra.
    """
    from .annotate import NUCLEAR_RECEPTORS

    rows = csv.DictReader(io.StringIO(text), delimiter="\t")
    ref = {r.get("Approved symbol") for r in rows
           if r.get("Locus type") == "gene with protein product"} - {None}
    if not ref:
        raise ValueError("the HGNC download held no protein-coding genes")
    return {"reference": len(ref),
            "missing_from_ours": sorted(ref - NUCLEAR_RECEPTORS),
            "not_in_reference": sorted(NUCLEAR_RECEPTORS - ref)}


def check_residue_codes(client) -> dict:
    """Do the docking residue tables agree with the PDB component dictionary?

    A code in `MODIFIED_RESIDUES` must be a peptide-linking component, or a
    structure's real ligand could be skipped. A solvent or crystallisation
    additive must not be one, or part of the protein would be ignored, and a
    sugar must be a saccharide. A code the dictionary does not know is
    reported, because it can never match.
    """
    from . import docking, exploratory

    out = {}
    for label, codes, needle, expect in (
            ("MODIFIED_RESIDUES", docking.MODIFIED_RESIDUES,
             "PEPTIDE LINKING", True),
            ("SOLVENTS", docking.SOLVENTS, "PEPTIDE LINKING", False),
            ("CRYSTALLISATION_ADDITIVES", docking.CRYSTALLISATION_ADDITIVES,
             "PEPTIDE LINKING", False),
            ("exploratory.SUGARS", exploratory.SUGARS, "SACCHARIDE", True)):
        res = client.post(RCSB_GRAPHQL, {
            "query": _CHEMCOMP_QUERY, "variables": {"ids": sorted(codes)}})
        comps = ((res or {}).get("data") or {}).get("chem_comps") or []
        types = {c["chem_comp"]["id"]: (c["chem_comp"].get("type") or "")
                 for c in comps if c and c.get("chem_comp")}
        out[label] = {
            "codes": len(codes),
            "not_found": sorted(set(codes) - set(types) - _LEGACY_CODES),
            "wrong_type": sorted(
                k for k, t in types.items()
                if (needle in t.upper()) != expect)}
    return out


def main(argv=None) -> int:
    import argparse
    import sys

    ap = argparse.ArgumentParser(
        prog="generoutes.audit",
        description="Check the hand-curated tables for internal consistency "
                    "and, where a download is on disk, against it. With "
                    "--online, also check them against ChEMBL, HGNC and the "
                    "PDB, to detect drift.")
    ap.add_argument("--online", action="store_true",
                    help="also check curated drug names against ChEMBL "
                         "(one cached request per name), gene symbols and "
                         "the nuclear receptor list against HGNC, and "
                         "residue codes against the PDB")
    ap.add_argument("--limit", type=int,
                    help="only check this many names online, for a quick look")
    args = ap.parse_args(argv)

    _prov.use_utf8_console()
    logging.basicConfig(level=logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")

    problems = check_internal()
    print(summarise(problems))
    failed = any(problems.values())

    local = check_local_references()
    print("\nChecks against tables on disk:")
    if local["landmarks"] is None:
        print("  KNOWN_LANDMARKS: not checked, data/lincs_geneinfo.txt absent")
    elif local["landmarks"]:
        print("  ** NOT LANDMARKS in the CLUE gene table: "
              + ", ".join(local["landmarks"]))
        failed = True
    else:
        print("  KNOWN_LANDMARKS: every gene is a landmark in the CLUE table")
    ob = local["orangebook_routes"]
    if ob is None:
        print("  Orange Book routes: not checked, no Orange Book on disk")
    else:
        print(f"  Orange Book routes in neither table: "
              f"{', '.join(ob['unclassified']) or 'none'}")
        print(f"  Table routes this release does not use: "
              f"{', '.join(ob['not_in_this_release']) or 'none'}")

    if args.online:
        from .http import Client
        tables = all_curated_names()
        client = Client()
        total = sum(len(v) for v in tables.values())
        print(f"\nResolving {total} curated drug name(s) across "
              f"{len(tables)} table(s) against ChEMBL...\n")

        agg = {"checked": 0, "resolved": 0, "unresolved": [],
               "withdrawn": [], "errors": 0}
        for table, names in tables.items():
            r = resolve_names(client, names, limit=args.limit)
            agg["checked"] += r["checked"]
            agg["resolved"] += len(r["resolved"])
            agg["errors"] += r["lookup_errors"]
            agg["unresolved"] += [f"{table}: {n}" for n in r["unresolved"]]
            agg["withdrawn"] += [f"{table}: {n}"
                                 for n in r["withdrawn_in_chembl"]]
            print(f"  {table:26} {r['checked']:4} checked, "
                  f"{len(r['resolved']):4} resolved, "
                  f"{len(r['unresolved']):3} not, "
                  f"{r['lookup_errors']:3} error(s)")

        print(f"\nTOTAL {agg['checked']} checked: {agg['resolved']} resolved, "
              f"{len(agg['unresolved'])} did not, "
              f"{agg['errors']} lookup error(s)")
        if agg["unresolved"]:
            print("\nDid NOT resolve - renamed, merged, or wrong:")
            for n in agg["unresolved"]:
                print(f"  - {n}")
        if agg["withdrawn"]:
            print("\nChEMBL flags these withdrawn. Harmless in a burden "
                  "table, but read carefully in direct.CURATED, which the "
                  "report presents as therapies to consider:")
            for n in agg["withdrawn"]:
                print(f"  - {n}")
        if agg["errors"]:
            print(f"\n{agg['errors']} name(s) could not be checked because "
                  f"ChEMBL did not answer - those are UNKNOWN, not confirmed.")
        # The consequential check: our own withdrawn flags vs ChEMBL's.
        w = check_withdrawn_agreement(client, limit=args.limit)
        print(f"\nWithdrawn-flag agreement on direct.CURATED "
              f"({w['checked']} entries):")
        bad = w["we_say_available_chembl_says_withdrawn"]
        if bad:
            print("  ** WE SAY AVAILABLE, ChEMBL SAYS WITHDRAWN **")
            print("  These would sit at the top of the report unflagged:")
            for n in bad:
                print(f"    - {n}")
            failed = True
        else:
            print("  no entry we present as available is withdrawn in ChEMBL")
        if w["we_say_withdrawn_chembl_says_available"]:
            print("  we flag withdrawn, ChEMBL does not (approvals are "
                  "regional; we may be right):")
            for n in w["we_say_withdrawn_chembl_says_available"]:
                print(f"    - {n}")
        if w["unchecked"]:
            print(f"  {len(w['unchecked'])} could not be checked - UNKNOWN, "
                  f"not confirmed")

        print("\nNuclear receptor list against the HGNC group "
              "(annotate.NUCLEAR_RECEPTORS):")
        try:
            nr = check_nuclear_receptors(
                client.request(HGNC_NR_GROUP, parse_json=False))
        except Exception as exc:                            # noqa: BLE001
            print(f"  could not be checked - UNKNOWN, not confirmed "
                  f"({type(exc).__name__}: {str(exc)[:100]})")
        else:
            print(f"  HGNC lists {nr['reference']} protein-coding receptors")
            if nr["missing_from_ours"]:
                print("  ** HGNC LISTS THESE, WE DO NOT **")
                for n in nr["missing_from_ours"]:
                    print(f"    - {n}")
                failed = True
            else:
                print("  none missing from our list")
            if nr["not_in_reference"]:
                print("  we list these, HGNC does not:")
                for n in nr["not_in_reference"]:
                    print(f"    - {n}")

        print("\nGene symbols in the hand-entered lists against HGNC:")
        try:
            syms = all_curated_symbols()
            bad = check_gene_symbols(client, syms)
        except Exception as exc:                            # noqa: BLE001
            print(f"  could not be checked - UNKNOWN, not confirmed "
                  f"({type(exc).__name__}: {str(exc)[:100]})")
        else:
            total = len(set().union(*syms.values()))
            print(f"  {total} distinct symbols in {len(syms)} lists, "
                  f"{sum(len(v) for v in bad.values())} not approved")
            for table, items in sorted(bad.items()):
                print(f"    {table}: {', '.join(items)}")
            failed = failed or bool(bad)

        print("\nDocking residue codes against the PDB component dictionary:")
        try:
            rc = check_residue_codes(client)
        except Exception as exc:                            # noqa: BLE001
            print(f"  could not be checked - UNKNOWN, not confirmed "
                  f"({type(exc).__name__}: {str(exc)[:100]})")
        else:
            for label, r in rc.items():
                print(f"  {label:26} {r['codes']:3} codes, "
                      f"{len(r['not_found'])} unknown, "
                      f"{len(r['wrong_type'])} of the wrong type")
                for kind in ("not_found", "wrong_type"):
                    for n in r[kind]:
                        print(f"    - {n} ({kind.replace('_', ' ')})")
                if r["not_found"] or r["wrong_type"]:
                    failed = True

        print("\nResolution proves existence, not correctness, and no check "
              "here can tell you what is MISSING from a table.")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
