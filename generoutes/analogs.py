"""Approved drugs that resemble a known modulator but are not annotated at it.

The question this asks
----------------------
Something is already known to move your factor — an approved ligand, a CTD
entry, a LINCS hit. Chemical similarity says: what else that is already
approved looks like it?

Why the naive version is nearly useless
---------------------------------------
Run a similarity search on calcitriol and most of what comes back is already
in your drug table. That is not a coincidence. If a compound is structurally
similar to a known ligand **and** anyone has assayed it at that target, it is
already annotated there — so similarity rediscovers what ChEMBL already told
you. The rest is salt forms and esters of the query itself, which
`congeners` exists to collapse.

Measured on real queries: calcitriol returns alfacalcidol, calcifediol and
tacalcitol; tretinoin returns isotretinoin and alitretinoin; dexamethasone
returns betamethasone, triamcinolone and four salts. Nearly all of those are
already in the graded tables.

What is actually interesting
----------------------------
The leftover: a drug that **looks like a known modulator and is not recorded
at that target**. That is the repurposing quadrant — the pair nobody assayed
— and it is the only part of a similarity search worth reporting. Tacalcitol
is the worked example: an approved vitamin D analogue that did not appear in
a real CCR5 run's VDR list.

So this subtracts rather than adds. Everything already annotated at the
target is removed, every salt form of the query is collapsed, and what
survives is a short list of untested pairs.

The honest limit
----------------
Similarity is not activity. Two molecules at 76% Tanimoto can differ in the
one substituent the pocket cares about, and the whole premise fails for
targets where small changes flip agonist to antagonist — which is most
nuclear receptors. Every row here is a binding assay, not a finding, and the
hit rate for this kind of inference is low by construction.
"""

from __future__ import annotations

import logging
import urllib.parse
from dataclasses import dataclass, field

from .congeners import salt_free_name

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

# Tanimoto floor. Below ~70% the hits stop resembling the query in any way a
# binding pocket would notice; above ~85% almost everything returned is a
# salt or ester of the query itself.
DEFAULT_THRESHOLD = 70


@dataclass
class Analog:
    drug: str
    chembl_id: str | None
    similarity: float
    similar_to: str
    target: str
    note: str = ""
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def _parent_molecule(client, name: str) -> dict | None:
    """The parent compound for a name, not one of its esters or salts.

    A plain search for "tretinoin" returns TRETINOIN TOCOFERIL first, whose
    structure is the ester and whose similarity neighbourhood is therefore
    the wrong one. An exact name match wins; otherwise the shortest name,
    which is the parent by naming convention.
    """
    try:
        res = client.get(f"{CHEMBL}/molecule/search",
                         {"q": name, "format": "json", "limit": 10})
    except Exception as exc:                                # noqa: BLE001
        log.warning("ChEMBL lookup failed for %r: %s", name, str(exc)[:90])
        return None

    mols = [m for m in (res or {}).get("molecules", [])
            if (m.get("molecule_structures") or {}).get("canonical_smiles")]
    if not mols:
        return None
    want = name.strip().lower()
    exact = [m for m in mols
             if (m.get("pref_name") or "").strip().lower() == want]
    pool = exact or mols
    return min(pool, key=lambda m: len(m.get("pref_name") or "z" * 99))


def find_analogs(client, known_drug: str, target: str, *,
                 already_known: set | None = None,
                 threshold: int = DEFAULT_THRESHOLD,
                 limit: int = 100) -> list:
    """Approved drugs similar to `known_drug` and NOT annotated at `target`.

    `already_known` is the set of drug names this run has already found for
    the target, in any form; they are subtracted, because rediscovering them
    is what makes a naive similarity search look productive and be useless.
    """
    parent = _parent_molecule(client, known_drug)
    if not parent:
        return []
    smiles = parent["molecule_structures"]["canonical_smiles"]
    parent_name = parent.get("pref_name") or known_drug
    parent_stem = salt_free_name(parent_name)

    try:
        res = client.get(
            f"{CHEMBL}/similarity/{urllib.parse.quote(smiles)}/{threshold}",
            {"format": "json", "limit": limit, "max_phase": 4})
    except Exception as exc:                                # noqa: BLE001
        log.warning("similarity search failed for %s: %s",
                    known_drug, str(exc)[:90])
        return []

    known = {salt_free_name(n) for n in (already_known or set())}
    known.add(parent_stem)

    out: list = []
    seen: set = set()
    for h in (res or {}).get("molecules", []):
        name = (h.get("pref_name") or "").strip()
        if not name:
            continue
        stem = salt_free_name(name.lower())
        # Subtract: the query itself and its salts, and anything this run
        # already found at the target.
        if stem in known or stem in seen:
            continue
        seen.add(stem)
        try:
            sim = float(h.get("similarity"))
        except (TypeError, ValueError):
            continue
        out.append(Analog(
            drug=name.title() if name.isupper() else name,
            chembl_id=h.get("molecule_chembl_id"),
            similarity=round(sim, 1),
            similar_to=parent_name,
            target=target,
            note=(f"{round(sim)}% similar to {parent_name}, which is a known "
                  f"{target} ligand — and not itself recorded at {target}"),
            caveats=[
                "Similarity is not activity. Two molecules at this Tanimoto "
                "can differ in the one substituent the pocket cares about.",
                "For a nuclear receptor in particular, a small change can "
                "flip agonist to antagonist, so the DIRECTION is unknown "
                "even if the binding holds.",
            ]))
    out.sort(key=lambda a: -a.similarity)
    log.info("analogs of %s at %s: %d approved drug(s) not already annotated "
             "there", known_drug, target, len(out))
    return out


def search(client, graded, *, limit_seeds: int = 4,
           threshold: int = DEFAULT_THRESHOLD) -> dict:
    """Run the analog search from this report's own best candidates.

    Seeds are the top graded drugs, because they are the compounds with the
    strongest evidence of acting at the factor — so "what looks like this"
    is being asked of the best-supported starting point rather than an
    arbitrary one.
    """
    by_target: dict = {}
    for g in graded:
        by_target.setdefault(g.get("tf_symbol"), []).append(g)

    results: dict = {}
    for target, rows in list(by_target.items())[:limit_seeds]:
        if not target:
            continue
        known = {r.get("drug_name", "") for r in rows}
        seed = rows[0].get("drug_name")
        if not seed:
            continue
        found = find_analogs(client, seed, target,
                             already_known=known, threshold=threshold)
        if found:
            results[target] = {
                "seed": seed,
                "analogs": [a.to_dict() for a in found[:10]],
            }
    return results


def summarise(results: dict) -> str:
    if not results:
        return ("No approved drug resembles this report's candidates without "
                "already being annotated at the same factor. That is the "
                "common outcome and it is informative: if a compound looks "
                "like a known ligand and anyone has assayed it there, it is "
                "already in the tables above.")
    n = sum(len(v["analogs"]) for v in results.values())
    return (f"{n} approved drug(s) across {len(results)} factor(s) resemble a "
            f"known ligand without being recorded at that factor. Each is an "
            f"untested pair — a binding assay, not a finding — and the hit "
            f"rate for this kind of inference is low by construction.")
