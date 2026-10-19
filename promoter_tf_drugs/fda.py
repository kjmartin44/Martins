"""Show only drugs the FDA has approved.

Why this exists
---------------
`approved` in the rest of the package means ChEMBL's max phase is 4 or a source
called the drug approved. That is approval by any regulator. Pemafibrate is
approved in Japan and not by the FDA. The LINCS perturbation library adds
research tool compounds that no regulator has approved, and ten route
collectors mark everything they find as approved. A report that lists them
beside Orange Book drugs mixes candidates a person could be prescribed with
ones nobody can.

What counts as FDA approved here
--------------------------------
A drug is FDA approved when its name matches an active ingredient in the FDA
Orange Book (`orangebook.py`), with salt forms ignored, or when it is one of the
gene therapies in the curated direct therapy table (`direct.CURATED`), which
the Orange Book does not cover because it lists drugs and not biologics.

That rule is narrow on purpose and has a known cost. Antibodies and other
biologics are in the FDA's Purple Book, which this package does not read, so
an approved antibody is hidden along with the non-approved compounds. The
curated gene therapy list is written by hand and can be out of date. A name
that differs from the Orange Book's, such as a trade name or an older synonym,
is hidden too. `--include-non-fda` shows everything.

What is changed
---------------
Nothing is deleted. `results.json` keeps every route and every compound. This
module builds a view of a result with the non-approved drugs removed, and the
reports render that view. Counts inside the view, such as "271 routes", are
counts after the removal, and the report says how many entries were hidden.

If there is no Orange Book on disk the filter cannot tell approved from not,
and it is not applied. It never empties a report because a file is missing,
and the report says so.
"""

from __future__ import annotations

import copy
import logging

log = logging.getLogger(__name__)

NOTE_NO_TABLE = ("the FDA filter was requested but there is no Orange Book on "
                 "disk, so it could not tell approved drugs from the rest and "
                 "was not applied. See `orangebook.SOURCE` for the download.")

# Lists of drug entries inside one result block, as (path, key holding the
# drug's name). "*" is any dict key and "[]" is every list element. The keys
# differ because the lists were written at different times; "name" appears on
# transcription factors elsewhere, which is why this is a list of paths and not
# a search for any dict with a name.
PATHS = (
    ("drug_grades", "drug_name"),
    ("outcome_estimates", "drug"),
    ("activity_drugs.*", "drug"),
    ("corpus_scans.*.lincs.hits", "drug"),
    ("lincs.hits", "drug"),
    ("factor_abundance", "drug"),
    ("docking.[].hits", "drug"),
    ("mechanisms.*.hits", "drug"),
    ("mechanisms.*.by_mechanism.[].examples", None),
    ("screens.*.hits", "drug"),
    ("sensitivity.*.hits", "drug"),
    ("product_screen.hits", "drug"),
    ("improvements.*", "name"),
    ("improvements.reported_links.*", "chemical"),
    ("screen_plan.compounds.cherry_pick", "compound"),
    ("tfs.[].drugs", "name"),
    ("upstream.*.[].drugs", "name"),
    ("direction_unknown.*", "name"),
    ("analogs.*.analogs", "drug"),
    ("cell_type.delivery", "drug"),
    ("promoter_dna.candidates", "drug"),
    ("stabilisation.*.stabilisers", "drug"),
)

_CACHE: dict = {}
_INDEX: dict = {}

# Words after the active moiety that the Orange Book spells out and a candidate
# list may leave off, so "lapatinib" has to find "LAPATINIB DITOSYLATE". The
# counterions and esters are the package's own vocabulary; the rest are
# counterions it lacks.
_EXTRA_COUNTERIONS = frozenset({
    "ditosylate", "dimesylate", "disodium", "dimaleate", "ditartrate",
    "difumarate", "hemihydrate", "sesquihydrate", "hyclate", "tosilate",
    "dihydrobromide", "hemifumarate", "hemitartrate", "oxalate",
})


def reset() -> None:
    _CACHE.clear()
    _INDEX.clear()


def base_name(name: str) -> str:
    """The name with trailing counterion, ester and hydration words removed.

    Looser than `congeners.salt_free_name`, which keeps esters because they are
    different molecules. Approval is asked about the active ingredient as the
    Orange Book spells it, and it spells "REZAFUNGIN ACETATE" where a candidate
    list says "rezafungin". The cost is that an ester prodrug is read as its
    parent, which is also an approved drug in nearly every case.
    """
    from . import congeners
    drop = congeners.SALTS | congeners.ESTERS | congeners.NOISE | _EXTRA_COUNTERIONS
    toks = [t for t in congeners._tokens(name)]
    while len(toks) > 1 and toks[-1] in drop:
        head = toks[:-1]
        if all(t in drop for t in head):
            break
        toks = head
    return " ".join(toks)


def _index() -> frozenset:
    if "ob" not in _INDEX:
        from . import orangebook
        _INDEX["ob"] = frozenset(base_name(k) for k in orangebook.ingredients())
    return _INDEX["ob"]


def available() -> bool:
    """True when an Orange Book is loaded, so approval can be told."""
    from . import orangebook
    return orangebook.load() > 0


def _curated_gene_therapies() -> frozenset:
    from .congeners import salt_free_name
    from . import direct
    out = set()
    for rows in direct.CURATED.values():
        for t in rows:
            if t.modality == "gene therapy":
                out.add(salt_free_name(t.drug.lower()))
    return frozenset(out)


def is_approved(name: str) -> bool:
    """Whether the name is an FDA approved drug by the rule above."""
    key = (name or "").strip().lower()
    if not key:
        return False
    if key in _CACHE:
        return _CACHE[key]
    from . import orangebook
    from .congeners import salt_free_name
    ok = (orangebook.lookup(key) is not None
          or base_name(key) in _index()
          or salt_free_name(key) in _curated_gene_therapies())
    _CACHE[key] = ok
    return ok


def status(*, demo: bool, requested: bool = True) -> dict:
    """What a run records about the filter: asked for, applied, and why not.

    Demo data is synthetic and its drug names are not real products, so the
    filter is never applied to it.
    """
    if not requested:
        return {"requested": False, "applied": False,
                "note": "shown with every compound, as requested"}
    if demo:
        return {"requested": True, "applied": False,
                "note": "not applied to demo data, whose drugs are synthetic"}
    if not available():
        return {"requested": True, "applied": False, "note": NOTE_NO_TABLE}
    return {"requested": True, "applied": True, "note": None}


def _walk(node, tokens: list, key: str, drop: set) -> int:
    """Remove non-approved entries at `tokens` below `node`; return how many."""
    if not tokens:
        if not isinstance(node, list):
            return 0
        keep, gone = [], 0
        for item in node:
            # `key` None means the list holds the names themselves.
            nm = (item if key is None else
                  item.get(key) if isinstance(item, dict) else None)
            if isinstance(nm, str) and nm and not is_approved(nm):
                drop.add(nm)
                gone += 1
            else:
                keep.append(item)
        node[:] = keep
        return gone
    head, rest = tokens[0], tokens[1:]
    if head == "[]":
        return sum(_walk(x, rest, key, drop) for x in node) \
            if isinstance(node, list) else 0
    if not isinstance(node, dict):
        return 0
    if head == "*":
        return sum(_walk(v, rest, key, drop) for v in list(node.values()))
    return _walk(node.get(head), rest, key, drop) if head in node else 0


def _rebuild_routes(block: dict) -> None:
    """Recompute the route summaries after routes were removed."""
    from . import catalogue
    ar = block.get("all_routes")
    if not isinstance(ar, dict) or "routes" not in ar:
        return
    routes = ar["routes"]
    by_category: dict = {}
    by_tier: dict = {}
    best: dict = {}
    for r in routes:
        by_category[r["category"]] = by_category.get(r["category"], 0) + 1
        by_tier[r["tier"]] = by_tier.get(r["tier"], 0) + 1
        best.setdefault(r["category"], r)
    ar["n_routes"] = len(routes)
    ar["by_category"] = by_category
    ar["by_tier"] = {catalogue.TIERS[k]: v for k, v in sorted(by_tier.items())}
    ar["best_per_category"] = best
    ar["categories_absent"] = [c for c in catalogue.CATEGORIES
                               if c not in by_category]


def apply(result: dict) -> dict:
    """The result with non-approved drugs removed, or `result` itself.

    Returned unchanged unless the run recorded `fda_filter.applied`, so a
    hand built result, a demo, and a run without an Orange Book all render as
    they did. The returned copy carries `fda_filter.hidden`, the count of
    entries removed and the names, for the report to say.
    """
    flt = (result or {}).get("fda_filter") or {}
    if not flt.get("applied"):
        return result
    out = copy.deepcopy(result)
    dropped: set = set()
    n = 0
    for block in out.get("results") or []:
        routes_before = len((block.get("all_routes") or {}).get("routes") or [])
        for path, key in PATHS:
            n += _walk(block, path.split("."), key, dropped)
        # all_routes is not in PATHS: its removals feed the summaries below.
        ar = block.get("all_routes")
        if isinstance(ar, dict) and isinstance(ar.get("routes"), list):
            keep = []
            for r in ar["routes"]:
                if is_approved(r.get("drug") or ""):
                    keep.append(r)
                else:
                    dropped.add(r.get("drug") or "?")
            ar["routes"] = keep
            n += routes_before - len(keep)
            _rebuild_routes(block)
    dt = out.get("direct_therapies")
    if isinstance(dt, list):
        keep = [t for t in dt if not isinstance(t, dict)
                or is_approved(t.get("drug") or "")]
        n += len(dt) - len(keep)
        out["direct_therapies"] = keep
    out["fda_filter"] = dict(flt, hidden={"entries": n,
                                          "drugs": sorted(dropped, key=str.lower)})
    return out


def notice(result: dict) -> list:
    """Lines for the top of a report saying what the filter did."""
    flt = (result or {}).get("fda_filter") or {}
    if not flt.get("requested"):
        return []
    if not flt.get("applied"):
        return [f"> **FDA filter not applied:** {flt.get('note')}", ""]
    h = flt.get("hidden") or {}
    return [f"> **Only FDA approved drugs are shown.** {h.get('entries', 0)} "
            f"entries for {len(h.get('drugs', []))} compounds that are not in "
            f"the FDA Orange Book were left out: research compounds, drugs "
            f"approved only outside the United States, and biologics, which "
            f"the Orange Book does not list. `results.json` keeps all of "
            f"them. `--include-non-fda` shows everything.", ""]
