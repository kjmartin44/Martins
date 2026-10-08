"""Do the web services still return the shape the parsers expect?

The failure this catches
------------------------
Every parser here reaches into a response by name, defensively:

    hits = (((res.get("data") or {}).get("search") or {}).get("hits")) or []
    return hits[0]["id"] if hits else None

That style cannot raise, which is the point — one renamed field must
not take a thirty-minute run down. But it means a renamed field returns
`None`, and the caller writes "could not resolve TP53 in Open Targets".
A reader sees *this gene has no genetic evidence*. What happened is
*the API changed*.

So the package's central rule — a failed lookup is not a measured
negative — holds at the level of the evidence ledger and leaks one
layer below it. The ledger can only report what the parser hands back,
and a parser defeated by a rename hands back an honest-looking zero.

`genes._is_mane` already tests two spellings of the same field, and
`localdrugs` sniffs several column names per table. Those are scars.
Nothing systematically checks that a response still looks like itself.

What this does
--------------
One probe per service, against an entity that will not disappear —
TP53, aspirin — and then a check that the paths the parsers actually
walk are present and populated. It reports which module would break and
what the symptom would look like, because "ChEMBL changed" is not
actionable and "the screen would report no known ligands" is.

Three outcomes, never two
-------------------------
`ok`, `CONTRACT CHANGED`, `unreachable`. Conflating the last two is the
same mistake one level up: an outage is a reason to try later, a schema
change is a reason to edit the parser, and a check that cannot tell
them apart sends you to do the wrong one. A probe that times out is
reported as unreachable and never as a contract break.

Opt-in, and not free
--------------------
Nothing calls this during a run. Each probe is a live request, which is
why it is a command rather than a stage. It is deliberately NOT served
from the response cache: a probe answered from a copy kept for up to 30
days reports last month's schema and cannot see a change made last
week. The fresh response is still written to the cache, so a normal run
afterwards benefits.

    python -m promoter_tf_drugs.services
    python -m promoter_tf_drugs.services --key opentargets --verbose
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass, field
from typing import Callable

from .http import Client, HttpError

log = logging.getLogger(__name__)

OK = "ok"
CHANGED = "CONTRACT CHANGED"
UNREACHABLE = "unreachable"
SKIPPED = "skipped"
# The probe itself is broken - a bad import, a renamed constant. It
# says NOTHING about the service, and reporting it as an outage sends
# someone to check a host that is answering perfectly. This module
# exists to keep "it broke" and "it changed" apart; it has to keep its
# own failures apart from both.
BROKEN_PROBE = "PROBE BROKEN"


@dataclass(frozen=True)
class Contract:
    """One service, and what the code that reads it depends on."""
    key: str
    name: str
    why: str                    # what the package uses it for
    reads: str                  # which module parses it
    probe: Callable             # (client) -> parsed response
    requires: tuple             # dotted paths that must exist and be filled
    symptom: str                # what a silent break looks like in a report
    notes: tuple = field(default_factory=tuple)


def _dig(obj, path: str):
    """Walk a dotted path, with `[]` meaning "first element of a list".

    Returns a (found, value) pair rather than raising or returning None,
    because `None` is itself a value a service can legitimately send and
    the whole point here is telling absence from content.
    """
    cur = obj
    for part in path.split("."):
        if part.endswith("[]"):
            part = part[:-2]
            if part:
                if not isinstance(cur, dict) or part not in cur:
                    return False, None
                cur = cur[part]
            if not isinstance(cur, list) or not cur:
                return False, None
            cur = cur[0]
            continue
        if not isinstance(cur, dict) or part not in cur:
            return False, None
        cur = cur[part]
    return True, cur


# --------------------------------------------------------------- probes
# Each uses an entity chosen to outlive schema churn: TP53 is the most
# studied gene there is and aspirin the most studied drug. A probe that
# failed because its example was retired would be indistinguishable
# from one that failed because the schema moved.
def _ensembl(client):
    # The package's own function, not a request built here. A probe
    # that composes its own parameters checks a fiction: this exact
    # mistake hid the MANE bug, because `mane=1` was missing from the
    # real call and the first version of this probe did not send it
    # either, so the two agreed about the wrong thing.
    from .genes import _lookup_symbol
    return _lookup_symbol(client, "TP53", "homo_sapiens")


def _opentargets(client):
    from .genetics import OPENTARGETS
    q = """query ($q: String!) {
        search(queryString: $q, entityNames: ["target"]) {
            hits { id name entity }
        } }"""
    return client.post(OPENTARGETS, {"query": q, "variables": {"q": "TP53"}})


def _opentargets_assoc(client):
    from .genetics import OPENTARGETS, _ASSOC_QUERY
    return client.post(OPENTARGETS, {
        "query": _ASSOC_QUERY,
        "variables": {"ensemblId": "ENSG00000141510", "size": 2}})


def _chembl(client):
    from .drugs import CHEMBL
    return client.get(f"{CHEMBL}/molecule/search",
                      {"q": "aspirin", "format": "json", "limit": 2})


def _mygene(client):
    from .annotate import MYGENE
    return client.get(f"{MYGENE}/query",
                      {"q": "TP53", "species": "human", "size": 2})


def _jaspar(client):
    from .motifs import JASPAR_API
    return client.get(f"{JASPAR_API}/matrix/",
                      {"collection": "CORE", "tax_group": "vertebrates",
                       "page_size": 2, "format": "json"})


def registry() -> tuple:
    """Every contract, in the order a run would notice them breaking."""
    return (
        Contract(
            key="ensembl",
            name="Ensembl REST",
            why="resolves a symbol to a transcript and a promoter window",
            reads="genes.py",
            probe=_ensembl,
            # genes.py reads every one of these off the lookup.
            requires=("id", "seq_region_name", "start", "end", "strand",
                      "display_name", "Transcript[]", "Transcript[].id",
                      "Transcript[].start", "Transcript[].end"),
            symptom="the run stops at the first step with 'could not "
                    "resolve', which reads as a bad gene symbol",
            notes=("`genes._is_mane` already tests both `mane_select` and "
                   "`is_mane_select`, and both `tag` and `tags`. A third "
                   "spelling would silently cost MANE selection, not the "
                   "run - so it is checked separately below.",),
        ),
        Contract(
            key="ensembl_mane",
            name="Ensembl MANE flags",
            why="picks the representative transcript, so the promoter "
                "window is the right one",
            reads="genes._is_mane",
            probe=_ensembl,
            requires=(),      # checked by the predicate instead
            symptom="a non-representative TSS is scanned and nothing says "
                    "so; the window is simply the wrong piece of DNA",
        ),
        Contract(
            key="opentargets",
            name="Open Targets search",
            why="resolves a symbol to an Ensembl id before any genetics",
            reads="genetics._resolve_ensembl",
            probe=_opentargets,
            requires=("data.search.hits[]", "data.search.hits[].id"),
            symptom="'could not resolve X in Open Targets', so --genetics "
                    "reports no phenotypes - which reads as a gene with no "
                    "genetic evidence rather than a lookup that broke",
        ),
        Contract(
            key="opentargets_assoc",
            name="Open Targets associations",
            why="the human genetic support for moving the gene at all",
            reads="genetics.target_consequences",
            probe=_opentargets_assoc,
            requires=("data.target.associatedDiseases.rows[]",
                      "data.target.associatedDiseases.rows[].score",
                      "data.target.associatedDiseases.rows[].disease",
                      "data.target.associatedDiseases.rows[]"
                      ".datatypeScores[].id",
                      "data.target.associatedDiseases.rows[]"
                      ".datatypeScores[].score"),
            symptom="zero phenotypes, reported as a measured absence of "
                    "genetic support",
            notes=("The parser keeps only rows whose datatypeScores carry "
                   "`genetic_association`. Renaming that ID alone empties "
                   "the section while every path above still resolves, so "
                   "it is checked by value and not by shape.",),
        ),
        Contract(
            key="chembl",
            name="ChEMBL molecule search",
            why="the drug library and the known ligands a docking gate needs",
            reads="drugs.py, agonism.py",
            probe=_chembl,
            requires=("molecules[]", "molecules[].molecule_chembl_id",
                      "molecules[].pref_name", "page_meta.total_count"),
            symptom="an empty screening library, or a docking gate that "
                    "refuses for want of known ligands - both of which "
                    "look like a target nobody has studied",
        ),
        Contract(
            key="mygene",
            name="MyGene.info",
            why="the fallback symbol resolver when Ensembl cannot",
            reads="genes.py",
            probe=_mygene,
            requires=("hits[]", "hits[].symbol", "total"),
            symptom="a symbol Ensembl missed stays unresolved, so the gene "
                    "looks unrecognised",
        ),
        Contract(
            key="jaspar",
            name="JASPAR API",
            why="motif matrices when no local PFM file is given",
            reads="motifs.py",
            probe=_jaspar,
            requires=("results[]", "results[].matrix_id", "count"),
            symptom="the scan falls back to the small bundled set, so most "
                    "real factors are never tested",
        ),
    )


# Checks that a shape cannot express. Keyed by contract.
def _check_mane(res) -> str | None:
    """Can at least one transcript still be recognised as MANE Select?

    `genes._is_mane` accepts `mane_select`, `is_mane_select`, `tag` and
    `tags`. If TP53 — which has a MANE Select transcript — returns none
    that the predicate recognises, the spelling has moved again and
    every promoter window silently shifts to whatever sorts first.
    """
    from .genes import _is_mane
    trs = res.get("Transcript") or []
    if not trs:
        return "no Transcript list to look in"
    if not any(_is_mane(t) for t in trs):
        keys = sorted({k for t in trs for k in t if "mane" in k.lower()
                       or k in ("tag", "tags")})
        return ("no transcript of TP53 is recognisable as MANE Select; "
                f"fields present that look related: {keys or 'none'}")
    return None


def _check_genetic_association(res) -> str | None:
    """Is `genetic_association` still a datatype id?

    The parser keeps only rows carrying it, so renaming that one string
    empties the genetics section while every structural path still
    resolves. A shape check cannot see this.
    """
    rows = (((res.get("data") or {}).get("target") or {})
            .get("associatedDiseases") or {}).get("rows") or []
    seen = {d.get("id") for r in rows for d in (r.get("datatypeScores") or [])}
    if not seen:
        return "no datatypeScores on any row"
    if "genetic_association" not in seen:
        return (f"no row carries datatype 'genetic_association'; "
                f"ids seen: {sorted(x for x in seen if x)[:6]}")
    return None


VALUE_CHECKS = {
    "ensembl_mane": _check_mane,
    "opentargets_assoc": _check_genetic_association,
}


def check(contract: Contract, client: Client) -> dict:
    """Probe one service and compare what came back against the contract."""
    row = {"key": contract.key, "name": contract.name, "state": OK,
           "detail": "", "missing": [], "reads": contract.reads}

    if client.offline:
        row["state"] = SKIPPED
        row["detail"] = "offline; a contract check has to reach the service"
        return row

    # Ask the service, not the cache. Client keeps a response for 30 days,
    # so a probe answered from it reports last month's schema and a change
    # made last week goes unseen, which defeats the whole point of the
    # check. Restored afterwards so a caller's client is left as it was.
    previous = getattr(client, "bypass_cache", False)
    client.bypass_cache = True
    try:
        res = contract.probe(client)
    except HttpError as exc:
        # Unreachable, NOT changed. An outage is a reason to try later;
        # a schema change is a reason to edit a parser. A check that
        # cannot tell them apart sends you to do the wrong one.
        row["state"] = UNREACHABLE
        row["detail"] = str(exc)[:200]
        return row
    except (ImportError, AttributeError, NameError, TypeError) as exc:
        # Raised before any request went out: the probe is wrong, not
        # the service. Keeping this separate is the whole point.
        row["state"] = BROKEN_PROBE
        row["detail"] = (f"{type(exc).__name__}: {str(exc)[:160]} "
                         f"— this is a fault in the check, not in "
                         f"{contract.name}")
        return row
    except Exception as exc:                                # noqa: BLE001
        row["state"] = UNREACHABLE
        row["detail"] = f"{type(exc).__name__}: {str(exc)[:160]}"
        return row
    finally:
        client.bypass_cache = previous

    missing = []
    for path in contract.requires:
        found, value = _dig(res, path)
        if not found:
            missing.append(path)
        elif value in ((), [], {}, "", None):
            # Present but empty is still a break for a probe chosen to
            # return something: TP53 has transcripts, aspirin exists.
            missing.append(f"{path} (present but empty)")

    if missing:
        row["state"] = CHANGED
        row["missing"] = missing
        row["detail"] = f"{len(missing)} expected field(s) absent or empty"
        return row

    value_check = VALUE_CHECKS.get(contract.key)
    if value_check is not None:
        try:
            problem = value_check(res)
        except Exception as exc:                          # noqa: BLE001
            # The value checks read the response with .get(). A service
            # that starts returning a list, a string or null raises here,
            # and that is exactly the event this command exists to
            # report: it must not end as a traceback.
            problem = (f"the response is {type(res).__name__}, not the "
                       f"object the parser reads "
                       f"({type(exc).__name__}: {str(exc)[:80]})")
        if problem:
            row["state"] = CHANGED
            row["detail"] = problem
            return row

    row["detail"] = "every field the parsers read is present"
    return row


def check_all(client: Client | None = None, *, keys=None) -> list:
    """Every contract, or the ones named."""
    client = client or Client()
    wanted = registry()
    if keys:
        wanted = tuple(c for c in wanted if c.key in set(keys))
    return [check(c, client) for c in wanted]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.services",
        description="Check that each web service still returns the shape "
                    "this package's parsers read. A renamed field does "
                    "not raise here - it returns an honest-looking zero - "
                    "so nothing else catches it.")
    p.add_argument("--key", action="append",
                   help="check only this service; repeatable")
    p.add_argument("--offline", action="store_true",
                   help="do not reach the network. Every check is then "
                        "reported as skipped, which is the honest result "
                        "rather than a pass")
    p.add_argument("--cache-dir")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="list every field checked, not only the missing")
    args = p.parse_args(argv)

    known = {c.key for c in registry()}
    unknown = sorted(set(args.key or []) - known)
    if unknown:
        print(f"no such service: {', '.join(unknown)}. "
              f"Known: {', '.join(sorted(known))}", file=sys.stderr)
        return 2

    client = Client(cache_dir=args.cache_dir, offline=args.offline) \
        if args.cache_dir else Client(offline=args.offline)
    rows = check_all(client, keys=args.key)
    by_key = {c.key: c for c in registry()}

    width = max(len(r["key"]) for r in rows)
    changed = [r for r in rows if r["state"] == CHANGED]
    down = [r for r in rows if r["state"] == UNREACHABLE]
    broken = [r for r in rows if r["state"] == BROKEN_PROBE]

    print(f"{len(rows)} contract(s) checked, {len(changed)} changed, "
          f"{len(down)} unreachable, {len(broken)} probe(s) broken\n")
    for r in rows:
        c = by_key[r["key"]]
        mark = {OK: "ok  ", CHANGED: "****", UNREACHABLE: "down",
                SKIPPED: "--  ", BROKEN_PROBE: "??  "}[r["state"]]
        print(f"  [{mark}] {r['key']:<{width}}  {c.name}")
        print(f"         {c.why}")
        print(f"         {r['state']}: {r['detail']}")
        for path in r["missing"]:
            print(f"           - {path}")
        if r["state"] == CHANGED:
            print(f"         parsed in {c.reads}")
            print(f"         if unfixed: {c.symptom}")
        if args.verbose:
            for note in c.notes:
                print(f"         note: {note}")
            if r["state"] == OK and c.requires:
                print(f"         checked: {', '.join(c.requires)}")
        print()

    if changed:
        print("A changed contract does not raise during a run. It returns "
              "an empty result, which the report shows as a source that "
              "found nothing — so fix the parser rather than waiting for "
              "a failure that will not come.")
    if broken:
        print("A probe that cannot run says nothing about its service. "
              "Fix the check.")
    # Unreachable is not a failure of this command: the service being
    # down says nothing about whether its schema moved. A broken probe
    # IS a failure, because it means a contract is going unchecked.
    return 1 if (changed or broken) else 0


if __name__ == "__main__":        # pragma: no cover
    logging.basicConfig(level=logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    raise SystemExit(main())
