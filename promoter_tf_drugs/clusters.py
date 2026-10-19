"""Fold near-identical drugs into one row of a ranking.

Why this exists
---------------
Drugs come in families. A factor that has one good ligand usually has a dozen
close relatives that bind it equally well, and they score the same. Listed
separately they fill the top of a ranking with one idea repeated, and push out
the second and third ideas the reader came for. In a CCR5 list the top ten was
five vitamin D analogues and two retinoids.

What this does
--------------
Walks a ranking from the best row down and folds each drug into the best-ranked
row it is nearly identical to, so the best-ranked member leads its cluster and
the others are named beside it. Two drugs are nearly identical when either:

  * they are the same drug under another ester or prodrug name
    (`congeners.ester_free_name` matches), or
  * their structures, with counter-ions removed, have a Tanimoto similarity of
    at least 0.75 on Morgan fingerprints of radius 2 (2,048 bits). That is
    close, as an analogue series is, not merely related.

A caller may add a condition on top. The summary folds two drugs only when they
act through the same target, since the same structure acting through different
factors is two routes. The combined list folds on structure alone, because a
drug found only by docking has no route to compare.

What it does not do
-------------------
**Similar structure is not equivalent behaviour.** Close analogues can differ in
potency, selectivity, tissue reach and safety. Folding them hides those
differences from the top of the list, so the folded members are always named
and their own ranks kept, and the full list is unchanged underneath.

**It is greedy, by rank.** A cluster is every drug close to its leader, not a
chain of drugs each close to the last, so a long series does not merge into one
blob. A drug close to two leaders joins the better-ranked one.

**It needs a structure.** A compound with none in the approved drug table (a
research compound from a signature file, for instance) stays on its own unless
its name shows it is an ester of a listed drug. Without RDKit only the name
test runs and the report says so.

**The threshold is a convention.** 0.75 is a common cut for "highly similar" on
this fingerprint. It was not tuned on anything here.
"""

from __future__ import annotations

import logging

from .congeners import ester_free_name
from .signatures import breadth_key

log = logging.getLogger(__name__)

SIMILARITY = 0.75
RADIUS = 2
N_BITS = 2048

_SMILES: dict = {"index": None}
_FP: dict = {}
_WARNED: dict = {"rdkit": False}


# ---------------------------------------------------------------- structures
def smiles_index(rows: list | None = None) -> dict:
    """`{normalised name: smiles}` from the approved drug table.

    Read once and remembered. Pass `rows` to use another table (tests do);
    that replaces what is remembered. Empty when no table is on disk.
    """
    if rows is not None:
        _SMILES["index"] = _index(rows)
        _FP.clear()
    elif _SMILES["index"] is None:
        from . import localdrugs
        _SMILES["index"] = _index(localdrugs.library(None))
    return _SMILES["index"]


def _index(rows: list) -> dict:
    out: dict = {}
    for r in rows:
        key = breadth_key(r.get("name") or "")
        smi = r.get("smiles")
        if key and smi and key not in out:
            out[key] = smi
    return out


def reset() -> None:
    """Forget the remembered table and fingerprints."""
    _SMILES["index"] = None
    _FP.clear()


def fingerprint(smiles: str | None):
    """A Morgan bit vector of the drug without its counter-ions, or None."""
    if not smiles:
        return None
    if smiles in _FP:
        return _FP[smiles]
    fp = None
    try:
        from rdkit import Chem, DataStructs        # noqa: F401
        from rdkit.Chem import rdFingerprintGenerator
        from .docking import split_counter_ions
        parent = split_counter_ions(smiles)["smiles"]
        mol = Chem.MolFromSmiles(parent)
        if mol is not None:
            gen = rdFingerprintGenerator.GetMorganGenerator(
                radius=RADIUS, fpSize=N_BITS)
            fp = gen.GetFingerprint(mol)
    except ImportError:
        if not _WARNED["rdkit"]:
            _WARNED["rdkit"] = True
            log.warning("RDKit is not installed, so drugs are folded by name "
                        "only (ester and prodrug forms), not by structure")
    except Exception as exc:                                # noqa: BLE001
        log.debug("no fingerprint for %s: %s", smiles[:40], exc)
    _FP[smiles] = fp
    return fp


def tanimoto(a, b) -> float:
    from rdkit import DataStructs
    return float(DataStructs.TanimotoSimilarity(a, b))


# ------------------------------------------------------------------ grouping
def group(items: list, *, smiles_of=None, same_route=None,
          threshold: float = SIMILARITY) -> list:
    """Cluster a best-first list of drugs.

    `items` are dicts with at least `name`. `smiles_of(item)` gives a SMILES or
    None (default: the approved drug table, by name). `same_route(leader,
    item)` is an extra condition; return False to keep two drugs apart.

    Returns clusters in the order of their leaders, each
    `{"leader": item, "members": [item, ...], "similarity": {name: value}}`,
    where members exclude the leader and `similarity` is each member's
    similarity to it (1.0 for a name match).
    """
    if smiles_of is None:
        index = smiles_index()

        def smiles_of(item, _index=index):                  # noqa: E306
            return _index.get(breadth_key(item.get("name") or ""))
    clusters: list = []
    for item in items:
        name = item.get("name") or ""
        fam = ester_free_name(name)
        fp = fingerprint(smiles_of(item))
        placed = False
        for cl in clusters:
            leader = cl["leader"]
            if same_route is not None and not same_route(leader, item):
                continue
            sim = None
            if ester_free_name(leader.get("name") or "") == fam and fam:
                sim = 1.0
            elif fp is not None and cl["fp"] is not None:
                s = tanimoto(fp, cl["fp"])
                if s >= threshold:
                    sim = s
            if sim is not None:
                cl["members"].append(item)
                cl["similarity"][name] = round(sim, 3)
                placed = True
                break
        if not placed:
            clusters.append({"leader": item, "members": [], "similarity": {},
                             "fp": fp})
    for cl in clusters:
        del cl["fp"]
    return clusters


def names(cluster: dict, limit: int = 6) -> str:
    """The folded drugs as text: "A, B, C and 4 more"."""
    ns = [m.get("name") for m in cluster["members"]]
    if len(ns) <= limit:
        return ", ".join(ns)
    return ", ".join(ns[:limit]) + f" and {len(ns) - limit} more"
