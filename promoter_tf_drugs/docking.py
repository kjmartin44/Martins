"""Structure-based prediction of functional direction, done the way that can
actually work.

The mistake to avoid
--------------------
Docking a drug into one receptor structure and reading the score tells you
about *affinity*, not direction. Agonists and antagonists of a nuclear receptor
bind the same pocket with comparable scores. Ranking by docking score will
happily hand you tamoxifen when you asked for an ERalpha agonist.

Conformational selection
------------------------
What does carry directional information is *which receptor conformation the
ligand is compatible with*. Nuclear receptor LBDs crystallise in two
distinguishable states:

  active (agonist-bound)     helix 12 packs against the LBD, completing the
                             AF-2 groove. Co-crystals very often include an
                             LXXLL coactivator peptide occupying that groove.
  inactive (antagonist)      H12 is displaced into the coactivator groove
                             itself, often with a corepressor peptide or
                             nothing bound there.

So: dock the ligand into an ensemble of *both* kinds of structure and compare.

    delta = best_score_in_active - best_score_in_inactive

A ligand that fits the active conformation better than the inactive one is
agonist-like; the reverse is antagonist-like. This is a real and much weaker
signal than people expect, so the output is reported with the raw scores, never
as a bare verdict.

Detecting which is which
------------------------
`classify_conformation` uses the presence of a short peptide chain carrying an
LXXLL motif (the NR box) as evidence of an active, coactivator-bound structure,
and a CoRNR box (LXXI/HIXXXI/L) as evidence of a corepressor-bound one. This is
a structural fact you can read straight out of the PDB file rather than a
literature lookup, which makes it automatable.

Requirements (none bundled - this is opt-in)
--------------------------------------------
    pip install rdkit meeko
    smina:  https://github.com/mwaskom/smina  (single static binary)
      or AutoDock Vina: https://github.com/ccsb-scripps/AutoDock-Vina/releases
    Structures are fetched from RCSB: https://files.rcsb.org/download/XXXX.pdb

Honest expectation: on retrospective NR benchmarks this kind of ensemble
approach separates agonists from antagonists meaningfully better than chance
and far worse than a ligand-based classifier trained on the same receptor. Use
agonism.py first. Use this when the receptor has too few labelled ligands to
train on, or when you want a mechanistic picture of *why* a compound reads as
an antagonist.
"""

from __future__ import annotations

import logging
import math
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field, asdict, field
from pathlib import Path

log = logging.getLogger(__name__)

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILE = "https://files.rcsb.org/download"

# LXXLL, the coactivator NR box. L-x-x-L-L in one-letter code.
NR_BOX = re.compile(r"L..LL")
# CoRNR box, corepressor: L-x-x-I/H-I-x-x-x-I/L (NCoR/SMRT interaction motif)
CORNR_BOX = re.compile(r"L..[IH]I...[IL]")

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


@dataclass
class Conformation:
    pdb_id: str
    state: str                      # "active" | "inactive" | "unknown"
    evidence: str
    peptide_chains: list = field(default_factory=list)
    path: str | None = None


@dataclass
class DockingResult:
    drug_name: str
    target_symbol: str
    best_active: float | None = None      # kcal/mol, more negative = better
    best_inactive: float | None = None
    delta: float | None = None            # active - inactive
    predicted: str = "undetermined"
    n_active_structures: int = 0
    n_inactive_structures: int = 0
    per_structure: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# -------------------------------------------------------------- PDB parsing
def chain_sequences(pdb_text: str) -> dict[str, str]:
    """One-letter sequence per chain, from ATOM/HETATM CA records."""
    seqs: dict[str, list] = {}
    seen: set = set()
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[12:16].strip() != "CA":
            continue
        resn = line[17:20].strip().upper()
        chain = line[21].strip() or "_"
        resi = line[22:27].strip()
        key = (chain, resi)
        if key in seen:
            continue
        seen.add(key)
        seqs.setdefault(chain, []).append(AA3_TO_1.get(resn, "X"))
    return {c: "".join(s) for c, s in seqs.items()}


def classify_conformation(pdb_text: str, pdb_id: str = "?",
                          max_peptide_len: int = 35) -> Conformation:
    """Active vs inactive, from bound coactivator/corepressor peptides.

    A short chain carrying an LXXLL motif means a coactivator is sitting in the
    AF-2 groove, which is only possible when helix 12 is in its active position.
    That makes it agonist-conformation evidence you can read off the file.
    """
    seqs = chain_sequences(pdb_text)
    if not seqs:
        return Conformation(pdb_id, "unknown", "no CA records parsed")

    longest = max(len(s) for s in seqs.values())
    peptides = {c: s for c, s in seqs.items()
                if len(s) <= max_peptide_len and len(s) < longest}

    coact = [c for c, s in peptides.items() if NR_BOX.search(s)]
    corep = [c for c, s in peptides.items() if CORNR_BOX.search(s)]

    if coact and not corep:
        return Conformation(pdb_id, "active",
                            f"LXXLL coactivator peptide in chain(s) {','.join(coact)}",
                            list(peptides))
    if corep and not coact:
        return Conformation(pdb_id, "inactive",
                            f"CoRNR corepressor peptide in chain(s) {','.join(corep)}",
                            list(peptides))
    if coact and corep:
        return Conformation(pdb_id, "unknown", "both motifs matched", list(peptides))
    return Conformation(pdb_id, "unknown",
                        "no coregulator peptide; conformation must be assigned "
                        "another way (H12 geometry or literature)", list(peptides))


# ------------------------------------------------------------ RCSB fetching
def find_structures(client, uniprot: str, limit: int = 25) -> list[str]:
    """PDB IDs for a UniProt accession, via the RCSB search API."""
    query = {
        "query": {
            "type": "terminal", "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers."
                             "reference_sequence_identifiers.database_accession",
                "operator": "exact_match", "value": uniprot,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    try:
        res = client.post(RCSB_SEARCH, query)
    except Exception as exc:                      # noqa: BLE001
        log.warning("RCSB search failed for %s: %s", uniprot, exc)
        return []
    return [r["identifier"] for r in res.get("result_set", [])]


def fetch_pdb(client, pdb_id: str, out_dir: Path) -> Path | None:
    out = Path(out_dir) / f"{pdb_id}.pdb"
    if out.exists():
        return out
    try:
        text = client.request(f"{RCSB_FILE}/{pdb_id}.pdb", parse_json=False)
    except Exception as exc:                      # noqa: BLE001
        log.debug("could not fetch %s: %s", pdb_id, exc)
        return None
    out.write_text(text)
    return out


# ------------------------------------------------------------------ docking
def _docking_binary() -> str | None:
    for name in ("smina", "smina.static", "vina", "gnina"):
        path = shutil.which(name)
        if path:
            return path
    return None


def dock_one(receptor_pdb: Path, ligand_sdf: Path, *, autobox_ligand: Path | None = None,
             exhaustiveness: int = 16, seed: int = 0) -> float | None:
    """Run one docking job; return the best score in kcal/mol."""
    binary = _docking_binary()
    if binary is None:
        raise RuntimeError(
            "no docking binary found on PATH. Install smina "
            "(https://github.com/mwaskom/smina) or AutoDock Vina and make sure "
            "it is callable as 'smina' or 'vina'."
        )

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "docked.sdf"
        cmd = [binary, "-r", str(receptor_pdb), "-l", str(ligand_sdf),
               "-o", str(out), "--exhaustiveness", str(exhaustiveness),
               "--seed", str(seed), "--cpu", "1"]
        # Box: reuse the crystallographic ligand when we have one, which is the
        # only reliable way to place the box on the orthosteric pocket.
        if autobox_ligand:
            cmd += ["--autobox_ligand", str(autobox_ligand), "--autobox_add", "4"]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        except subprocess.TimeoutExpired:
            log.warning("docking timed out: %s", receptor_pdb.name)
            return None
        if proc.returncode != 0:
            log.warning("docking failed (%s): %s", receptor_pdb.name,
                        proc.stderr.strip()[:200])
            return None

        scores = [float(m.group(1)) for m in
                  re.finditer(r"minimizedAffinity[^\d\-]*(-?\d+\.\d+)", out.read_text())]
        if not scores:
            scores = [float(m.group(1)) for m in
                      re.finditer(r"^\s+1\s+(-?\d+\.\d+)", proc.stdout, re.M)]
        return min(scores) if scores else None


def smiles_to_sdf(smiles: str, path: Path, n_confs: int = 1, seed: int = 0) -> Path:
    """3D-embed and minimise a ligand ready for docking."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"could not parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    with Chem.SDWriter(str(path)) as w:
        w.write(mol)
    return path


def predict_direction(client, drug, smiles: str, uniprot: str, *,
                      work_dir: str | Path = ".docking", max_structures: int = 8,
                      exhaustiveness: int = 16) -> DockingResult:
    """Ensemble-dock into active and inactive conformations and compare."""
    res = DockingResult(drug_name=drug.name, target_symbol=drug.target_symbol)
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    pdb_ids = find_structures(client, uniprot)
    if not pdb_ids:
        res.notes.append(f"no PDB structures found for {uniprot}")
        return res

    active, inactive = [], []
    for pid in pdb_ids[:max_structures * 2]:
        path = fetch_pdb(client, pid, work)
        if path is None:
            continue
        conf = classify_conformation(path.read_text(), pid)
        conf.path = str(path)
        if conf.state == "active" and len(active) < max_structures:
            active.append(conf)
        elif conf.state == "inactive" and len(inactive) < max_structures:
            inactive.append(conf)

    res.n_active_structures, res.n_inactive_structures = len(active), len(inactive)
    if not active or not inactive:
        res.notes.append(
            f"need both conformations; found {len(active)} active / "
            f"{len(inactive)} inactive. Without both, a docking score says "
            "nothing about direction.")
        return res

    ligand = smiles_to_sdf(smiles, work / f"{drug.name.replace(' ', '_')}.sdf")

    for group, confs in (("active", active), ("inactive", inactive)):
        for conf in confs:
            score = dock_one(Path(conf.path), ligand, exhaustiveness=exhaustiveness)
            if score is not None:
                res.per_structure[f"{conf.pdb_id} ({group})"] = score

    act = [s for k, s in res.per_structure.items() if "(active)" in k]
    ina = [s for k, s in res.per_structure.items() if "(inactive)" in k]
    if not act or not ina:
        res.notes.append("docking produced no usable scores")
        return res

    res.best_active, res.best_inactive = min(act), min(ina)
    res.delta = res.best_active - res.best_inactive

    # More negative delta = fits the active conformation better = agonist-like.
    if abs(res.delta) < 0.5:
        res.predicted = "undetermined"
        res.notes.append(f"delta {res.delta:+.2f} kcal/mol is within docking noise "
                         "(~1 kcal/mol); no directional call")
    else:
        res.predicted = "agonist" if res.delta < 0 else "antagonist"
        res.notes.append(f"delta {res.delta:+.2f} kcal/mol. Docking scores carry "
                         "errors near 1-2 kcal/mol, so treat this as weak evidence "
                         "to be checked against agonism.py and the literature.")
    return res


# ====================================================================
# AlphaFold structures — useful, and the easiest way to fool yourself
# ====================================================================
"""
AlphaFold gives you a structure for every TF, including the ones with no PDB
entry. That is genuinely enabling, and it is also why docking into AlphaFold
models has a poor reputation: benchmarks repeatedly show enrichment far below
what holo crystal structures give, for two reasons that both apply here.

**Apo side chains.** Predicted models are unliganded. Side chains sit in
rotamers that have not reorganised around a ligand, and cryptic pockets are
closed. Rigid-receptor docking into that under-scores real binders.

**Low pLDDT is a disorder prediction, not a pocket.** This is the trap specific
to transcription factors. Their activation domains are frequently intrinsically
disordered; AlphaFold reports that as pLDDT below ~50-70. Those regions look
like open, greasy cavities to a docking program, which will happily score
compounds into them with excellent affinities. Every one of those hits is an
artefact. The pLDDT filter below is not optional hygiene — it is the difference
between a hypothesis and a fabrication.

Use AlphaFold docking to triage which pockets are worth a real experiment,
never as evidence of binding.
"""

ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"

# Below this mean pLDDT a region is predicted-disordered, not a binding site.
PLDDT_ORDERED = 70.0
PLDDT_CONFIDENT = 90.0


@dataclass
class PocketAssessment:
    uniprot: str
    mean_plddt: float
    ordered_fraction: float
    n_residues: int
    dockable: bool
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def residue_plddt(pdb_text: str) -> dict:
    """AlphaFold stores per-residue pLDDT in the B-factor column."""
    out: dict = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            resi = int(line[22:26])
            out[resi] = float(line[60:66])
        except ValueError:
            continue
    return out


def fetch_alphafold(client, uniprot: str, out_dir) -> "Path | None":
    from pathlib import Path
    out = Path(out_dir) / f"AF-{uniprot}.pdb"
    if out.exists():
        return out
    try:
        text = client.request(ALPHAFOLD_FILE.format(uniprot=uniprot),
                              parse_json=False)
    except Exception as exc:                              # noqa: BLE001
        log.warning("AlphaFold fetch failed for %s: %s", uniprot, exc)
        return None
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    return out


def assess_model(pdb_text: str, uniprot: str = "?") -> PocketAssessment:
    """Decide whether this predicted structure can be docked into at all."""
    plddt = residue_plddt(pdb_text)
    if not plddt:
        return PocketAssessment(uniprot, 0.0, 0.0, 0, False,
                                ["no pLDDT values parsed"])
    vals = list(plddt.values())
    mean = sum(vals) / len(vals)
    ordered = sum(1 for v in vals if v >= PLDDT_ORDERED) / len(vals)

    a = PocketAssessment(uniprot, round(mean, 1), round(ordered, 3), len(vals),
                         dockable=ordered >= 0.5)
    if mean >= PLDDT_CONFIDENT:
        a.notes.append("high-confidence model; still apo, so side chains are not "
                       "in ligand-bound rotamers and cryptic pockets stay closed")
    elif ordered < 0.5:
        a.notes.append(
            f"only {ordered:.0%} of residues are confidently folded (pLDDT>="
            f"{PLDDT_ORDERED:.0f}). This is a disorder prediction, not a "
            "structure with pockets. Docking here produces confident artefacts — "
            "low-pLDDT regions look like open greasy cavities to a scoring "
            "function. Do not dock this model.")
    else:
        a.notes.append(
            f"{ordered:.0%} of residues confidently folded; restrict docking to "
            "those regions and discard any pose contacting low-pLDDT residues")
    return a


def ordered_residues(pdb_text: str, cutoff: float = PLDDT_ORDERED) -> set:
    return {r for r, v in residue_plddt(pdb_text).items() if v >= cutoff}


def pose_is_in_ordered_region(contact_residues, ordered: set,
                              min_fraction: float = 0.8) -> bool:
    """Reject poses that lean on predicted-disordered residues."""
    contacts = set(contact_residues)
    if not contacts:
        return False
    return len(contacts & ordered) / len(contacts) >= min_fraction


# ===================================================================== #
# Docking as a way to FIND compounds, not just to explain one.          #
# ===================================================================== #
# predict_direction above takes a drug you already have and asks which
# conformation it prefers. This asks the other question: dock the whole
# approved library into the pocket and see what comes out. It needs no
# prior ligands for the receptor, so unlike the ligand-based screen in
# screening.py it can propose a scaffold nothing was trained on.
#
# That independence is the reason to run it. A compound a ligand model and
# a dock both pick is supported by two arguments sharing no input.
#
# Two things it will not do, both worth saying plainly:
#
# It will not reach the factors you cannot currently drug. Docking needs a
# pocket and a structure, and most transcription factors have neither -
# they work through DNA-binding surfaces and intrinsically disordered
# activation domains with nothing concave to dock into. The receptors that
# dock well are the nuclear-receptor LBDs already dominating the report.
# `assess_model` exists to refuse the rest rather than produce numbers for
# them.
#
# And a docking score is not an affinity. These functions enrich for
# binders and order them badly. The output is a shortlist to assay, not a
# ranking to trust.

# The gate. Docking is the method here most able to produce confident
# nonsense, so before any score on an untested drug is reported, the same
# pocket must rank the receptor's KNOWN ligands above decoys.
MIN_DOCKING_AUC = 0.70
# Fold enrichment in the top 1%, which is what virtual screening is
# actually judged on and a harsher test than AUC: a method can hold a
# respectable AUC while putting nothing useful at the top, and the top is
# the only part anyone looks at. 1.0 is random.
MIN_EF1 = 2.0
MIN_VALIDATION_ACTIVES = 5


def roc_auc(scores: list, labels: list, *, lower_is_better: bool = True
            ) -> float | None:
    """AUC via Mann-Whitney U, ties at half. None entries are failed docks
    and are dropped rather than scored worst, which would flatter the
    method by burying every compound it could not handle."""
    pairs = [(s, y) for s, y in zip(scores, labels) if s is not None]
    if not pairs:
        return None
    if lower_is_better:
        pairs = [(-s, y) for s, y in pairs]
    pos = [s for s, y in pairs if y]
    neg = [s for s, y in pairs if not y]
    if not pos or not neg:
        return None
    wins = sum((1.0 if p > n else 0.5 if p == n else 0.0)
               for p in pos for n in neg)
    return wins / float(len(pos) * len(neg))


def enrichment_factor(scores: list, labels: list, *, fraction: float = 0.01,
                      lower_is_better: bool = True) -> float | None:
    """Fold enrichment of actives in the top `fraction` of the ranking."""
    pairs = [(s, y) for s, y in zip(scores, labels) if s is not None]
    if not pairs:
        return None
    n = len(pairs)
    n_act = sum(1 for _, y in pairs if y)
    if not n_act or n_act == n:
        return None
    k = max(1, int(math.ceil(n * fraction)))
    pairs.sort(key=lambda t: t[0] if lower_is_better else -t[0])
    hits = sum(1 for _, y in pairs[:k] if y)
    return round((hits / float(k)) / (n_act / float(n)), 2)


def _dock_smiles(smiles_list: list, receptor: Path, autobox: Path | None,
                 *, exhaustiveness: int, seed: int, score_fn=None) -> list:
    """Score each SMILES, None where docking or embedding failed.

    `score_fn` exists so the gate and the ranking can be exercised without
    a docking binary; it takes a SMILES and returns a score.
    """
    out: list = []
    for smi in smiles_list:
        if score_fn is not None:
            try:
                out.append(score_fn(smi))
            except Exception:                               # noqa: BLE001
                out.append(None)
            continue
        with tempfile.TemporaryDirectory() as tmp:
            sdf = Path(tmp) / "lig.sdf"
            try:
                smiles_to_sdf(smi, sdf, seed=seed)
            except Exception as exc:                        # noqa: BLE001
                log.debug("embedding failed for %s: %s", smi[:40], exc)
                out.append(None)
                continue
            try:
                out.append(dock_one(receptor, sdf, autobox_ligand=autobox,
                                    exhaustiveness=exhaustiveness, seed=seed))
            except Exception as exc:                        # noqa: BLE001
                log.debug("docking failed: %s", str(exc)[:80])
                out.append(None)
    return out


def validate_pocket(receptor: Path, actives: list, decoys: list, *,
                    autobox: Path | None = None, exhaustiveness: int = 8,
                    seed: int = 0, score_fn=None) -> dict:
    """Does docking into THIS pocket recover ligands known to bind it?

    `actives` and `decoys` are `[{"smiles": ...}]`. Returns the measured
    AUC and top-1% enrichment and whether they clear the bar. Nothing this
    module says about an untested compound should be reported when this
    says no.
    """
    if len(actives) < MIN_VALIDATION_ACTIVES:
        return {"usable": False, "auc": None, "ef1": None,
                "reason": (f"only {len(actives)} known ligand(s) for this "
                           f"pocket, below the {MIN_VALIDATION_ACTIVES} "
                           f"needed to test whether docking works here. An "
                           f"unvalidated dock is exactly the thing that "
                           f"produces confident nonsense")}
    if not decoys:
        return {"usable": False, "auc": None, "ef1": None,
                "reason": "no decoys to validate against"}

    smiles = [a["smiles"] for a in actives] + [d["smiles"] for d in decoys]
    labels = [1] * len(actives) + [0] * len(decoys)
    scores = _dock_smiles(smiles, receptor, autobox,
                          exhaustiveness=exhaustiveness, seed=seed,
                          score_fn=score_fn)

    auc = roc_auc(scores, labels)
    ef1 = enrichment_factor(scores, labels, fraction=0.01)
    n_failed = sum(1 for s in scores if s is None)
    if auc is None:
        return {"usable": False, "auc": None, "ef1": None,
                "n_failed": n_failed,
                "reason": "every dock failed, so nothing could be measured"}

    usable = auc >= MIN_DOCKING_AUC and (ef1 is None or ef1 >= MIN_EF1)
    out = {"usable": usable, "auc": round(auc, 3), "ef1": ef1,
           "n_actives": len(actives), "n_decoys": len(decoys),
           "n_failed": n_failed}
    if not usable:
        out["reason"] = (
            f"docking into this pocket ranks its own known ligands at AUC "
            f"{auc:.2f}"
            + (f" and {ef1}-fold enrichment in the top 1%" if ef1 else "")
            + f", against a bar of {MIN_DOCKING_AUC} AUC and {MIN_EF1}-fold. "
              f"It cannot separate this receptor's binders from random "
              f"approved drugs, so its scores on untested compounds carry "
              f"no information and are not reported.")
    return out


@dataclass
class DockingScreen:
    """A docking campaign against one pocket, with its validation."""
    target: str
    n_docked: int = 0
    n_failed: int = 0
    auc: float | None = None
    ef1: float | None = None
    usable: bool = False
    hits: list = field(default_factory=list)
    error: str | None = None
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def screen_by_docking(receptor: Path, library: list, *, target: str = "?",
                      actives: list | None = None, decoys: list | None = None,
                      autobox: Path | None = None, top_n: int = 25,
                      exhaustiveness: int = 8, seed: int = 0,
                      score_fn=None) -> DockingScreen:
    """Dock an approved library into one pocket, if the pocket validates.

    `library` is `[{"name", "smiles"}]`. The validation runs first and the
    library is not touched unless it passes, because the expensive part is
    worth nothing when the cheap check says the pocket cannot be scored.
    """
    run = DockingScreen(target=target)

    if score_fn is None and _docking_binary() is None:
        run.error = (
            "no docking binary on PATH. Install smina "
            "(https://github.com/mwaskom/smina) or AutoDock Vina and make "
            "it callable as 'smina' or 'vina'. Nothing else in the run is "
            "affected.")
        return run

    val = validate_pocket(receptor, actives or [], decoys or [],
                          autobox=autobox, exhaustiveness=exhaustiveness,
                          seed=seed, score_fn=score_fn)
    run.auc, run.ef1, run.usable = val.get("auc"), val.get("ef1"), \
        bool(val.get("usable"))
    if not run.usable:
        run.error = val.get("reason")
        return run

    scores = _dock_smiles([d["smiles"] for d in library], receptor, autobox,
                          exhaustiveness=exhaustiveness, seed=seed,
                          score_fn=score_fn)
    scored = [{"drug": d.get("name") or "?", "smiles": d["smiles"],
               "score": s}
              for d, s in zip(library, scores) if s is not None]
    scored.sort(key=lambda h: h["score"])
    run.n_docked = len(scored)
    run.n_failed = len(library) - len(scored)
    run.hits = scored[:top_n]

    run.caveats.append(
        f"this pocket validated at AUC {run.auc}"
        + (f" with {run.ef1}-fold enrichment in the top 1%" if run.ef1 else "")
        + ", so docking separates its known binders from random approved "
          "drugs here. That is a statement about enrichment, not about the "
          "order of the hits: a docking function pulls binders toward the "
          "top of a list and ranks them poorly once there")
    run.caveats.append(
        "a score says a molecule can occupy the pocket. It does not say "
        "whether the molecule activates or inhibits the receptor, and that "
        "direction decides whether a hit raises or lowers your gene. "
        "predict_direction above is the tool for that question")
    if run.n_failed:
        run.caveats.append(
            f"{run.n_failed} compound(s) could not be docked, usually "
            f"because no 3D conformer could be generated. They are absent "
            f"from the ranking rather than ranked badly")
    return run


def consensus(docking: DockingScreen, ligand_hits: list, *,
              top_frac: float = 0.2) -> list:
    """Compounds in the top slice of BOTH the dock and the ligand screen.

    The reason to run two methods. A docking score and a ligand-based
    model share no input: one reads the pocket, the other reads chemistry
    measured at the receptor. Agreement between them is worth more than a
    better score from either, and disagreement is informative too -
    ligand-only hits resemble known actives, dock-only hits are the novel
    scaffolds the trained model could not have recognised.
    """
    if not docking or not docking.usable or not ligand_hits:
        return []
    n_d = max(1, int(len(docking.hits) * top_frac)) or 1
    dock_top = {h["drug"].strip().lower(): h for h in docking.hits[:n_d]}
    n_l = max(1, int(len(ligand_hits) * top_frac)) or 1

    out = []
    for h in ligand_hits[:n_l]:
        name = (h.get("drug") or "").strip().lower()
        if name in dock_top:
            out.append({
                "drug": h.get("drug"),
                "docking_score": dock_top[name]["score"],
                "ligand_quadrant": h.get("quadrant"),
                "why": ("picked by both a pocket-shape argument and a "
                        "chemistry argument trained on this receptor's "
                        "measured ligands. The two share no input, so this "
                        "is the strongest computational support available "
                        "here - and still a binding hypothesis, with the "
                        "direction unaddressed by either"),
            })
    return out
