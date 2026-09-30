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
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, asdict, field
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
