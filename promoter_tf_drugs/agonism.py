"""Predict the *functional direction* of a drug at a receptor from its structure.

Read this before using it
-------------------------
The obvious approach - dock the drug into the receptor and look at the score -
does not work for this question, and it is worth being explicit about why.

Agonism in a nuclear receptor is not decided by binding affinity. It is decided
by whether helix 12 (the AF-2 helix) folds into its active position, closing
the hydrophobic groove that recruits LXXLL coactivator motifs. Antagonists
overwhelmingly bind the *same* pocket as agonists, with similar affinity, but
carry a bulky substituent that protrudes past the pocket and physically blocks
H12 from docking - tamoxifen's dimethylaminoethoxy arm, raloxifene's
piperidine, mifepristone's 11-beta dimethylaminophenyl. Estradiol and tamoxifen
both dock well into ERalpha. A docking score cannot tell them apart because
affinity and direction are nearly orthogonal properties.

So this module uses two things that *are* informative:

1. A supervised classifier trained on ChEMBL's own AGONIST/ANTAGONIST
   mechanism labels for that specific receptor, over Morgan fingerprints.
   This works because the structural determinants of antagonism are largely
   scaffold-local and recur across a receptor's chemotypes. It is honest about
   its limits: it reports cross-validated AUC, refuses to predict when the
   training set is too small, and refuses when the query is outside the
   applicability domain (max Tanimoto to training set below a threshold).

2. An interpretable "H12 protrusion" descriptor: bulk in excess of the known
   agonist set, plus a basic amine on a flexible tether - the recurring
   structural signature of NR antagonists. Reported alongside, never on its
   own, because it is a heuristic with obvious counterexamples.

For structure-based work that *can* address direction, see docking.py, which
does conformational-selection docking against agonist- and antagonist-bound
structures separately rather than scoring a single pocket.

What none of this does: tell you whether the drug will change your target
gene's expression in your cells. Direction at the receptor is one link in the
chain, not the answer.
"""

from __future__ import annotations

import logging
import statistics
from dataclasses import dataclass, asdict, field

import numpy as np

from .http import Client, HttpError

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

AGONIST_LABELS = {"AGONIST", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR",
                  "POSITIVE MODULATOR", "ACTIVATOR", "OPENER"}
ANTAGONIST_LABELS = {"ANTAGONIST", "INVERSE AGONIST", "INHIBITOR", "BLOCKER",
                     "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR",
                     "DEGRADER"}

MIN_PER_CLASS = 8          # below this, don't pretend to have a model
MIN_AUC = 0.70             # below this, report the model as unreliable
DOMAIN_TANIMOTO = 0.35     # below this, the query is out of applicability domain


def _rdkit():
    try:
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
        RDLogger.DisableLog("rdApp.*")
        return Chem, DataStructs, AllChem, Descriptors, rdMolDescriptors
    except ImportError as exc:
        raise ImportError(
            "agonism prediction needs RDKit and scikit-learn:\n"
            "    pip install rdkit scikit-learn"
        ) from exc


@dataclass
class AgonismPrediction:
    drug_name: str
    target_symbol: str
    predicted: str                  # "agonist" | "antagonist" | "undetermined"
    probability_agonist: float | None = None
    confidence: str = "none"        # high | moderate | low | none
    model_auc: float | None = None
    n_train_agonist: int = 0
    n_train_antagonist: int = 0
    nearest_tanimoto: float | None = None
    in_domain: bool = False
    protrusion: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------- bioactivity (binding) data
# The potency floor for calling a molecule a measured binder. pChEMBL 6 is
# 1 uM, the conventional line between a real ligand and an assay artefact.
# Named so the report can cite the same number the fetch applies.
PCHEMBL_CUTOFF = 6.0


def fetch_actives(client: Client, target_chembl_id: str, *,
                  pchembl_cutoff: float = PCHEMBL_CUTOFF,
                  limit: int = 2000) -> list[dict]:
    """Molecules with measured binding at this target, from /activity.

    This is a different and far larger table than /mechanism. Mechanism
    annotations exist only for drugs whose pharmacology a curator wrote up —
    a few dozen per target. Activity records exist for everything anyone ever
    assayed, which for a well-studied nuclear receptor is thousands.

    Using /mechanism to answer "what binds this receptor" was the reason every
    screen reported 1-27 ligands and refused to run.
    """
    actives: dict[str, dict] = {}
    offset = 0
    page = 250

    while len(actives) < limit:
        try:
            res = client.get(f"{CHEMBL}/activity", {
                "target_chembl_id": target_chembl_id,
                "pchembl_value__gte": pchembl_cutoff,
                "limit": page, "offset": offset, "format": "json",
            })
        except HttpError as exc:
            log.debug("activity fetch stopped at offset %d: %s", offset, exc)
            break
        rows = res.get("activities", [])
        if not rows:
            break
        for a in rows:
            mol = a.get("molecule_chembl_id")
            smiles = a.get("canonical_smiles")
            if not mol:
                continue
            try:
                pch = float(a.get("pchembl_value"))
            except (TypeError, ValueError):
                continue
            prev = actives.get(mol)
            if prev is None or pch > prev["pchembl"]:
                actives[mol] = {"chembl_id": mol, "smiles": smiles,
                                "pchembl": pch}
        offset += page
        if not (res.get("page_meta") or {}).get("next"):
            break

    # Activity rows do not always carry the structure; fill the gaps.
    missing = [m for m, v in actives.items() if not v["smiles"]]
    for i in range(0, len(missing), 50):
        chunk = missing[i:i + 50]
        try:
            res = client.get(f"{CHEMBL}/molecule", {
                "molecule_chembl_id__in": ",".join(chunk),
                "limit": 50, "format": "json"})
        except HttpError:
            continue
        for m in res.get("molecules", []):
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if smi:
                actives[m["molecule_chembl_id"]]["smiles"] = smi

    out = [v for v in actives.values() if v["smiles"]]
    log.info("%s: %d actives at pChEMBL >= %.1f", target_chembl_id,
             len(out), pchembl_cutoff)
    return out


# ------------------------------------------------------------ training data
def fetch_labelled_ligands(client: Client, target_chembl_id: str) -> list[dict]:
    """All ChEMBL molecules with an action_type at this target, with SMILES."""
    try:
        res = client.get(f"{CHEMBL}/mechanism", {
            "target_chembl_id": target_chembl_id, "limit": 500, "format": "json",
        })
    except HttpError as exc:
        log.warning("mechanism fetch failed for %s: %s", target_chembl_id, exc)
        return []

    labels: dict[str, str] = {}
    for m in res.get("mechanisms", []):
        mol = m.get("molecule_chembl_id")
        action = (m.get("action_type") or "").upper()
        if not mol:
            continue
        if action in AGONIST_LABELS:
            labels[mol] = "agonist"
        elif action in ANTAGONIST_LABELS:
            labels.setdefault(mol, "antagonist")

    out: list[dict] = []
    ids = list(labels)
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        try:
            res = client.get(f"{CHEMBL}/molecule", {
                "molecule_chembl_id__in": ",".join(chunk),
                "limit": 50, "format": "json",
            })
        except HttpError:
            continue
        for mol in res.get("molecules", []):
            smiles = (mol.get("molecule_structures") or {}).get("canonical_smiles")
            if not smiles:
                continue
            out.append({
                "chembl_id": mol["molecule_chembl_id"],
                "smiles": smiles,
                "label": labels[mol["molecule_chembl_id"]],
                "name": mol.get("pref_name"),
                "mw": (mol.get("molecule_properties") or {}).get("full_mwt"),
            })
    return out


# ------------------------------------------------------------- descriptors
# Free-text mechanism phrases, as a drug table writes them. ChEMBL gives
# an action_type from a controlled vocabulary; a table like the Broad
# Repurposing Hub gives prose. These are the words that decide direction
# in that prose, longest first so "inverse agonist" is not read as
# "agonist".
_MOA_ANTAGONIST = (
    "inverse agonist", "negative allosteric modulator", "antagonist",
    "inhibitor", "blocker", "degrader", "downregulator", "suppressor",
    "negative modulator",
)
_MOA_AGONIST = (
    "positive allosteric modulator", "partial agonist", "agonist",
    "activator", "opener", "stimulant", "positive modulator",
    "upregulator",
)


def label_from_text(text: str | None) -> str | None:
    """`"CCR5 antagonist"` -> `"antagonist"`. None when it does not say.

    Checked against the antagonist vocabulary first and longest phrase
    first, because the words nest: "inverse agonist" contains "agonist"
    and means the opposite of it, and a substring test in the other
    order labels every inverse agonist as an agonist.

    Returns None rather than guessing. A wrong direction label is worse
    than a missing one - it trains the classifier to predict backwards,
    and nothing downstream could tell.
    """
    if not text:
        return None
    low = str(text).strip().lower()
    if not low:
        return None
    for phrase in _MOA_ANTAGONIST:
        if phrase in low:
            return "antagonist"
    for phrase in _MOA_AGONIST:
        if phrase in low:
            return "agonist"
    return None


def labelled_ligands_local(symbol: str, rows: list) -> list[dict]:
    """Direction-labelled ligands of one target, from a local drug table.

    The offline counterpart of `fetch_labelled_ligands`. Rows need a
    mechanism string and a SMILES; rows whose mechanism names no
    direction are dropped rather than defaulted.
    """
    want = (symbol or "").strip().upper()
    if not want:
        return []
    out, seen = [], set()
    for r in rows or []:
        smi = r.get("smiles")
        if not smi or smi in seen:
            continue
        targets = r.get("targets") or []
        if isinstance(targets, str):
            targets = [targets]
        if want not in {str(t).strip().upper() for t in targets}:
            continue
        label = label_from_text(r.get("moa") or r.get("mechanism"))
        if label is None:
            continue
        seen.add(smi)
        out.append({"smiles": smi, "label": label,
                    "name": r.get("name") or "?", "mw": r.get("mw")})
    return out


def fingerprints(smiles_list: list[str], radius: int = 2, n_bits: int = 2048):
    """Morgan fingerprints as a dense array, plus the RDKit bit vectors."""
    Chem, DataStructs, AllChem, _, _ = _rdkit()
    from rdkit.Chem import rdFingerprintGenerator

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    arrs, bvs, keep = [], [], []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        bv = gen.GetFingerprint(mol)
        arr = np.zeros(n_bits, dtype=np.int8)
        DataStructs.ConvertToNumpyArray(bv, arr)
        arrs.append(arr)
        bvs.append(bv)
        keep.append(i)
    return (np.array(arrs) if arrs else np.empty((0, n_bits))), bvs, keep


def protrusion_features(smiles: str, agonist_mws: list[float] | None = None) -> dict:
    """The recurring structural signature of an NR antagonist.

    Bulk beyond the agonist range, plus a basic amine on a rotatable tether.
    Interpretable, and wrong often enough that it is reported as context
    rather than used as a decision on its own.
    """
    Chem, _, _, Descriptors, rdMolDescriptors = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"parsed": False}

    basic_n = Chem.MolFromSmarts("[NX3;H0,H1,H2;!$(N[#6]=[O,N,S]);!$(N[a])]")
    n_basic = len(mol.GetSubstructMatches(basic_n)) if basic_n else 0

    mw = Descriptors.MolWt(mol)
    rot = Descriptors.NumRotatableBonds(mol)
    feats = {
        "parsed": True,
        "mw": round(mw, 1),
        "rotatable_bonds": rot,
        "basic_amines": n_basic,
        "tpsa": round(Descriptors.TPSA(mol), 1),
        "logp": round(Descriptors.MolLogP(mol), 2),
    }

    if agonist_mws:
        ref = statistics.median(agonist_mws)
        feats["mw_excess_vs_agonists"] = round(mw - ref, 1)
        # The classic antagonist pattern: bigger than the agonists, with a
        # flexible basic tail long enough to reach past the pocket.
        feats["antagonist_like_bulk"] = bool(
            mw - ref > 60 and n_basic >= 1 and rot >= 4)
    return feats


# ------------------------------------------------------------------- model
@dataclass
class AgonismModel:
    clf: object
    train_fps: np.ndarray
    train_bvs: list
    labels: np.ndarray
    auc: float | None
    n_agonist: int
    n_antagonist: int
    agonist_mws: list[float]

    def nearest_tanimoto(self, bv) -> float:
        from rdkit import DataStructs
        if not self.train_bvs:
            return 0.0
        return float(max(DataStructs.BulkTanimotoSimilarity(bv, self.train_bvs)))


def train_model(ligands: list[dict], seed: int = 0) -> AgonismModel | None:
    """Train a per-receptor agonist/antagonist classifier. None if infeasible."""
    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import StratifiedKFold, cross_val_predict
        from sklearn.metrics import roc_auc_score
    except ImportError as exc:
        raise ImportError("pip install scikit-learn") from exc

    smiles = [l["smiles"] for l in ligands]
    X, bvs, keep = fingerprints(smiles)
    if len(keep) == 0:
        return None
    y = np.array([1 if ligands[i]["label"] == "agonist" else 0 for i in keep])

    n_ag, n_ant = int(y.sum()), int((1 - y).sum())
    if n_ag < MIN_PER_CLASS or n_ant < MIN_PER_CLASS:
        # Expected for most receptors now that direction is an optional
        # annotation rather than a gate, so this is debug-level noise.
        log.debug("too few direction labels (%d agonist / %d antagonist); "
                  "direction will not be annotated", n_ag, n_ant)
        return None

    clf = RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                 class_weight="balanced", random_state=seed, n_jobs=-1)

    auc = None
    try:
        folds = min(5, n_ag, n_ant)
        cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
        proba = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
        auc = float(roc_auc_score(y, proba))
    except Exception as exc:                      # noqa: BLE001
        log.debug("cross-validation failed: %s", exc)

    clf.fit(X, y)
    ag_mws = [float(ligands[i]["mw"]) for i in keep
              if y[list(keep).index(i)] == 1 and ligands[i].get("mw")] if keep else []

    return AgonismModel(clf=clf, train_fps=X, train_bvs=bvs, labels=y, auc=auc,
                        n_agonist=n_ag, n_antagonist=n_ant, agonist_mws=ag_mws)


def predict(model: AgonismModel | None, smiles: str, *, drug_name: str = "?",
            target_symbol: str = "?") -> AgonismPrediction:
    """Classify one molecule, with applicability-domain and AUC gating."""
    pred = AgonismPrediction(drug_name=drug_name, target_symbol=target_symbol,
                             predicted="undetermined")
    if model is None:
        pred.notes.append("no model: too few labelled ligands at this receptor")
        return pred

    pred.model_auc = round(model.auc, 3) if model.auc is not None else None
    pred.n_train_agonist = model.n_agonist
    pred.n_train_antagonist = model.n_antagonist

    X, bvs, keep = fingerprints([smiles])
    if len(keep) == 0:
        pred.notes.append("could not parse SMILES")
        return pred

    sim = model.nearest_tanimoto(bvs[0])
    pred.nearest_tanimoto = round(sim, 3)
    pred.in_domain = sim >= DOMAIN_TANIMOTO
    pred.protrusion = protrusion_features(smiles, model.agonist_mws)

    p_ag = float(model.clf.predict_proba(X)[0, 1])
    pred.probability_agonist = round(p_ag, 3)

    if not pred.in_domain:
        pred.notes.append(
            f"OUT OF APPLICABILITY DOMAIN (nearest training ligand Tanimoto "
            f"{sim:.2f} < {DOMAIN_TANIMOTO}); prediction withheld")
        return pred

    if model.auc is not None and model.auc < MIN_AUC:
        pred.notes.append(f"model discriminates poorly (CV AUC {model.auc:.2f}); "
                          "treat the call as unreliable")
        pred.confidence = "low"
    else:
        margin = abs(p_ag - 0.5)
        pred.confidence = "high" if margin > 0.30 else "moderate" if margin > 0.15 else "low"

    pred.predicted = "agonist" if p_ag >= 0.5 else "antagonist"

    if pred.protrusion.get("antagonist_like_bulk") and pred.predicted == "agonist":
        pred.notes.append("structural flag disagrees: bulk + basic amine is the "
                          "classic H12-displacing antagonist pattern")
        pred.confidence = "low"
    return pred


def predict_for_drug(client: Client, drug, smiles: str | None = None,
                     _model_cache: dict | None = None) -> AgonismPrediction:
    """Convenience wrapper: build/reuse the model for the drug's target."""
    cache = _model_cache if _model_cache is not None else {}
    target_id = getattr(drug, "target_chembl_id", None)
    if not target_id:
        p = AgonismPrediction(drug.name, drug.target_symbol, "undetermined")
        p.notes.append("no ChEMBL target id; cannot build a receptor-specific model")
        return p

    if target_id not in cache:
        cache[target_id] = train_model(fetch_labelled_ligands(client, target_id))

    if smiles is None:
        smiles = (getattr(drug, "properties", None) or {}).get("canonical_smiles")
    if not smiles and drug.chembl_id:
        try:
            mol = client.get(f"{CHEMBL}/molecule/{drug.chembl_id}", {"format": "json"})
            smiles = (mol.get("molecule_structures") or {}).get("canonical_smiles")
        except HttpError:
            pass
    if not smiles:
        p = AgonismPrediction(drug.name, drug.target_symbol, "undetermined")
        p.notes.append("no structure available")
        return p

    return predict(cache[target_id], smiles,
                   drug_name=drug.name, target_symbol=drug.target_symbol)
