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

A caller may also pass `same_family`, which folds two drugs whatever their
structures when it returns True. The summary uses it for one mechanism family:
drugs that reach the gene through the same factor by the same kind of route
(see `family_keys`). Six unrelated chemotypes that all push one receptor the
right way are one idea for choosing what to test first, and take one row.

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

# Routes where the drug is recorded as acting on the factor in one way. Other
# categories do not qualify: an abundance or activity change is inferred from
# signatures, so unrelated drugs land there for unrelated reasons, and upstream
# signalling names a node, not a mechanism. Two drugs there stay two rows.
FAMILY_CATEGORIES = ("TF: orthosteric ligand", "TF: allosteric modulator",
                     "TF: endogenous ligand", "TF: stabiliser or glue",
                     "TF: degradation blocker", "TF: known regulator")
NO_FACTOR = {"", "not stated", "(unknown / whole cell)", "promoter DNA",
             "chromatin (genome-wide)"}


def family_keys(routes: list) -> set:
    """The (route category, factor) pairs a drug reaches the gene through.

    Two drugs with a pair in common are in one mechanism family. Routes whose
    category does not name a factor give no pair, so they never fold this way.
    """
    out = set()
    for r in routes or []:
        cat = r.get("category") or ""
        tgt = r.get("acts_on") or ""
        if cat in FAMILY_CATEGORIES and tgt not in NO_FACTOR:
            out.add((cat, tgt))
    return out

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
def group(items: list, *, smiles_of=None, same_route=None, same_family=None,
          threshold: float = SIMILARITY) -> list:
    """Cluster a best-first list of drugs.

    `items` are dicts with at least `name`. `smiles_of(item)` gives a SMILES or
    None (default: the approved drug table, by name). `same_route(leader,
    item)` is an extra condition; return False to keep two drugs apart.
    `same_family(leader, item)` folds two drugs without a structure test when
    it returns True (after `same_route` has passed).

    Returns clusters in the order of their leaders, each
    `{"leader": item, "members": [item, ...], "similarity": {name: value}}`,
    where members exclude the leader and `similarity` is each member's
    similarity to it (1.0 for a name match, None for a family fold of drugs
    with no comparable structures). `basis` says why each member was folded:
    "name", "structure" or "family".
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
            sim, basis = None, None
            s = (tanimoto(fp, cl["fp"])
                 if fp is not None and cl["fp"] is not None else None)
            if ester_free_name(leader.get("name") or "") == fam and fam:
                sim, basis = 1.0, "name"
            elif s is not None and s >= threshold:
                sim, basis = s, "structure"
            elif same_family is not None and same_family(leader, item):
                sim, basis = s, "family"
            if basis is not None:
                cl["members"].append(item)
                cl["similarity"][name] = None if sim is None else round(sim, 3)
                cl["basis"][name] = basis
                placed = True
                break
        if not placed:
            clusters.append({"leader": item, "members": [], "similarity": {},
                             "basis": {}, "fp": fp})
    for cl in clusters:
        del cl["fp"]
    return clusters


def names(cluster: dict, limit: int = 6) -> str:
    """The folded drugs as text: "A, B, C and 4 more"."""
    ns = [m.get("name") for m in cluster["members"]]
    if len(ns) <= limit:
        return ", ".join(ns)
    return ", ".join(ns[:limit]) + f" and {len(ns) - limit} more"
