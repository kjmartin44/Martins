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
One cached request per name.

    python -m promoter_tf_drugs.audit             # offline only
    python -m promoter_tf_drugs.audit --online    # also resolve in ChEMBL

Neither level can tell you what is **missing** from a curated table — only a
human reading the current literature can. What this catches is drift in what
is already there, which is the failure that would otherwise be silent.
"""

from __future__ import annotations

import logging

from .congeners import salt_free_name

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


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
    from . import systemic, stabilizers, delivery, ligands, direct

    out: dict = {}
    out["systemic.CLASSES"] = sorted(
        {m for c in systemic.CLASSES for m in c.members})
    out["stabilizers.CURATED"] = sorted(
        {s.drug for lst in stabilizers.CURATED.values() for s in lst})
    out["stabilizers.EXEMPLARS"] = sorted(
        {s.drug for s in stabilizers.EXEMPLARS})
    out["delivery.PGP_SUBSTRATES"] = sorted(delivery.PGP_SUBSTRATES)
    try:
        out["ligands.ROUTES"] = sorted(
            {r.drug for lst in ligands.ROUTES.values() for r in lst})
    except Exception:                                       # noqa: BLE001
        pass
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
    except Exception:                                       # noqa: BLE001
        pass
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


def main(argv=None) -> int:
    import argparse
    import sys

    ap = argparse.ArgumentParser(
        prog="promoter-tf-drugs.audit",
        description="Check the hand-curated tables for internal consistency, "
                    "and optionally resolve every curated drug name against "
                    "ChEMBL to detect drift.")
    ap.add_argument("--online", action="store_true",
                    help="also resolve curated names against ChEMBL "
                         "(one cached request per name)")
    ap.add_argument("--limit", type=int,
                    help="only check this many names online, for a quick look")
    args = ap.parse_args(argv)

    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    logging.basicConfig(level=logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")

    problems = check_internal()
    print(summarise(problems))
    failed = any(problems.values())

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

        print("\nResolution proves existence, not correctness, and no check "
              "here can tell you what is MISSING from a table.")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
