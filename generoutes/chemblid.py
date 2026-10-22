"""Find the ChEMBL id of a drug that arrived without one.

Why this exists
---------------
Drug-gene tables name the same drug by different identifiers. DGIdb returns a
ChEMBL id for some drugs and an RxNorm ("rxcui:") or NCI Thesaurus ("ncit:")
id for others. The potency, selectivity and safety lookups in `drug_scoring`
all query ChEMBL by molecule id, and each began by giving up on anything that
did not start with "CHEMBL". In saved runs about half of the transcription
factor drugs were in that group (42 of 83 grades in the CYP3A4 run, 38 of 60 in
CYP1A1).

A drug given up on was not marked unknown. It was scored as average on potency
and selectivity and had its withdrawn and boxed-warning lookups skipped without
a note. So a nanomolar ligand and an untested one tied, a long list of drugs
shared one score, and a drug withdrawn for safety could pass with nothing said.
Rifampicin, the textbook ligand of PXR, was one of eight drugs at NR1I2 with
an identical grade for this reason.

A second gap sat beside it. Potency is looked up against the drug's target,
by ChEMBL target id, and DGIdb supplies none (all 42 of its drugs in the CYP3A4
run), so none of them got a potency even when the drug id was fine. The target
id is found here from the factor's gene symbol.

What this does
--------------
Resolves a name to a ChEMBL id in two steps, and says how it did:

1. The local approved drug table, which carries a ChEMBL id beside each name.
   Names are matched after case, hyphens and salt words are ignored.
2. Failing that, a ChEMBL name search, accepted only on an exact preferred
   name match. A near match is a different drug as often as the same one.

The id is then moved to the parent molecule. A salt form such as a mesylate has
a molecule id of its own, with few assay records, and the activity, mechanism
and warning records sit on the parent.

What it does not do
-------------------
It does not guess. A name with no exact match stays unresolved and the caller
records that, which is different from the drug having no data. It does not
change an id ChEMBL already supplied.
"""

from __future__ import annotations

import logging

from .http import HttpError
from .signatures import breadth_key

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

_LIBRARY: dict = {"index": None}
_TARGETS: dict = {}


def library_index(rows: list | None = None) -> dict:
    """`{normalised name: chembl id}` from the approved drug table.

    Read once and remembered. Pass `rows` to use another table (tests do);
    that replaces what is remembered.
    """
    if rows is not None:
        _LIBRARY["index"] = _index(rows)
    elif _LIBRARY["index"] is None:
        from . import localdrugs
        _LIBRARY["index"] = _index(localdrugs.library(None))
    return _LIBRARY["index"]


def _index(rows: list) -> dict:
    out: dict = {}
    for r in rows:
        key = breadth_key(r.get("name") or "")
        cid = str(r.get("chembl_id") or "")
        if key and cid.startswith("CHEMBL") and key not in out:
            out[key] = cid
    return out


def reset() -> None:
    """Forget the remembered table and target ids, so the next lookup reads
    them again."""
    _LIBRARY["index"] = None
    _TARGETS.clear()


def target_id(client, symbol: str | None) -> str | None:
    """The ChEMBL target id of a human protein from its gene symbol, or None.

    Asks ChEMBL for the human single-protein target that has this gene symbol
    among its synonyms. Complexes and protein-protein interaction targets that
    merely contain the protein are not single proteins and are not returned,
    since a potency measured against a complex is not one against the factor.
    No local table is needed. None when ChEMBL has no such target, or when the
    symbol matches more than one (a guess between them would be a guess); a
    failed request is logged and not remembered, so a later call can try
    again.
    """
    key = (symbol or "").strip().upper()
    if not key:
        return None
    if key in _TARGETS:
        return _TARGETS[key]
    try:
        res = client.get(f"{CHEMBL}/target.json", {
            "target_components__target_component_synonyms__component_synonym"
            "__iexact": key,
            "target_type": "SINGLE PROTEIN", "organism": "Homo sapiens",
            "limit": 5})
    except HttpError as exc:
        log.warning("ChEMBL target lookup failed for %s (%s)", key,
                    str(exc)[:80])
        return None
    ids = [t.get("target_chembl_id") for t in (res or {}).get("targets", [])
           if t.get("target_type") == "SINGLE PROTEIN"
           and t.get("target_chembl_id")]
    found = ids[0] if len(ids) == 1 else None
    if len(ids) > 1:
        log.info("%s matches %d ChEMBL targets (%s); none is used", key,
                 len(ids), ", ".join(ids))
    _TARGETS[key] = found
    return found


def parent_id(client, chembl_id: str) -> str:
    """The parent molecule's id, or `chembl_id` itself when it has none or the
    lookup fails (the failure is logged, not hidden)."""
    try:
        mol = client.get(f"{CHEMBL}/molecule/{chembl_id}.json")
    except HttpError as exc:
        log.warning("ChEMBL could not give the parent of %s (%s); using it as "
                    "it is", chembl_id, str(exc)[:80])
        return chembl_id
    return ((mol or {}).get("molecule_hierarchy") or {}).get(
        "parent_chembl_id") or chembl_id


def _search(client, name: str) -> str | None:
    """A ChEMBL id whose preferred name equals `name` exactly, else None."""
    try:
        res = client.get(f"{CHEMBL}/molecule/search.json",
                         {"q": name, "limit": 10})
    except HttpError as exc:
        log.warning("ChEMBL name search failed for %s (%s)", name,
                    str(exc)[:80])
        return None
    want = breadth_key(name)
    for m in (res or {}).get("molecules", []):
        if breadth_key(m.get("pref_name") or "") == want:
            return m.get("molecule_chembl_id")
    return None


def resolve(client, name: str, given: str | None) -> tuple:
    """`(chembl id or None, how)` for a drug.

    `how` is "given" when `given` was already a ChEMBL id (returned as it
    is), "approved table" or "name search" when it was found, and
    "unresolved" when it was not. A resolved id is the parent molecule's.
    """
    if given and str(given).startswith("CHEMBL"):
        return str(given), "given"
    if not (name or "").strip():
        return None, "unresolved"
    found, how = library_index().get(breadth_key(name)), "approved table"
    if found is None:
        found, how = _search(client, name), "name search"
    if found is None:
        return None, "unresolved"
    return parent_id(client, found), how
