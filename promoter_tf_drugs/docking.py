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
import threading
from concurrent.futures import ThreadPoolExecutor
import time
from dataclasses import dataclass, field, asdict, field
from pathlib import Path

log = logging.getLogger(__name__)

# The WSL search runs on every engine lookup; it is reported once.
_WSL_PROBE_REPORTED = False

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


def state_for_product(goal: str = "up") -> str:
    """Which conformation to dock into for the gene's OWN product.

    Simpler than `state_for`, which needs a transcription factor's role
    because a factor can activate or repress the gene it binds. The
    product has no such ambiguity: you are moving the thing itself, so
    wanting more of its activity means the agonist-bound conformation and
    wanting less means the antagonist-bound one.

    Worth its own function rather than calling `state_for("activator")`,
    which gives the same answer for a reason that does not apply.
    """
    return "active" if goal == "up" else "inactive"


def state_for(role: str, goal: str = "up") -> str | None:
    """Which receptor conformation to dock into, from role and goal.

    The direction depends on BOTH, which is the whole point:

        activator + raise  -> you want the factor ON   -> active
        repressor + raise  -> you want the factor OFF  -> inactive
        activator + lower  -> inactive
        repressor + lower  -> active

    `--state` used to be set by hand against help text that only
    described the activator case, so anyone working on a repressor had to
    invert it themselves or silently dock the wrong conformation.

    None when the role is unknown or dual, because this filter DISCARDS
    structures: guess wrong and the only pocket with a co-crystallised
    ligand is thrown away and the dock refused outright. That is stricter
    than `drugs.desired_direction`, which folds dual in with activator -
    there a wrong guess only reorders a table, here it can delete the
    evidence.
    """
    from .drugs import desired_direction
    if (role or "unknown") not in ("activator", "repressor"):
        return None
    want = desired_direction(role, goal)
    if want == "up":
        return "active"
    if want == "down":
        return "inactive"
    return None


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
    out.write_text(text, encoding="utf-8")
    return out


# ------------------------------------------------------------------ docking
def _docking_binary() -> str | None:
    """The best docking binary reachable, Windows PATH or WSL.

    WSL is tried FIRST, and gnina first within it, because the only
    reason to look there is that the better programs are Linux software.
    A Windows smina or vina is the fallback, not the preference.
    """
    try:
        from . import wsldock as _wsl
        name, path = _wsl.find_binary()
        if path:
            return path
    except Exception as exc:                                # noqa: BLE001
        global _WSL_PROBE_REPORTED
        if not _WSL_PROBE_REPORTED:
            _WSL_PROBE_REPORTED = True
            log.warning("the WSL docking-program search failed (%s: %s); "
                        "only Windows programs on PATH are considered",
                        type(exc).__name__, str(exc)[:100])
    for name in ("gnina", "smina", "smina.static", "vina"):
        path = shutil.which(name)
        if path:
            return path
    return None


def dock_one(receptor_pdb: Path, ligand_sdf: Path, *,
             autobox_ligand: Path | None = None,
             exhaustiveness: int = 16, seed: int = 0,
             cpu: int = 1) -> float | None:
    """Run one docking job; return the best score in kcal/mol.

    Goes through `_run_dock` rather than calling the binary directly,
    because the binary may live in WSL. Handing a Windows-side
    subprocess a path like /home/u/bin/gnina raises FileNotFoundError,
    and the caller catches that and records None - so every compound in
    a library scores None and the run reports an empty screen rather
    than a broken one.
    """
    binary = _docking_binary()
    if binary is None:
        raise RuntimeError(
            "no docking binary found. " + SETUP)

    from . import wsldock as _wsl
    conv = _wsl.to_wsl_path if binary.startswith("/") else (lambda q: str(q))

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "docked.sdf"
        args = ["-r", conv(receptor_pdb), "-l", conv(ligand_sdf),
                "-o", conv(out), "--exhaustiveness", str(exhaustiveness),
                "--seed", str(seed), "--cpu", str(cpu)]
        # Box: reuse the crystallographic ligand when we have one, which
        # is the only reliable way to place the box on the orthosteric
        # pocket.
        if autobox_ligand:
            args += ["--autobox_ligand", conv(autobox_ligand),
                     "--autobox_add", "4"]
        if "gnina" in binary:
            args += ["--cnn_scoring", "rescore"]
        stdout, sdf = _run_dock(binary, args, out)
        if stdout is None:
            return None
        poses = _parse_poses(sdf, stdout, engine_name(binary))
        affinities = [q.affinity for q in poses if q.affinity is not None]
        return min(affinities) if affinities else None


@dataclass
class PoseScore:
    """What a docking program said about one ligand's best pose.

    Three numbers because the programs answer different questions and
    disagree often enough that collapsing them would hide it.
    """
    affinity: float | None = None       # Vina-style kcal/mol, lower binds
    cnn_score: float | None = None      # 0-1, how pose-like the geometry is
    cnn_affinity: float | None = None   # predicted pK, larger binds
    engine: str = "?"

    def best(self, metric: str = "affinity") -> float | None:
        """One number for ranking, in the convention the caller asks for.

        `affinity` and `cnn_affinity` point opposite ways - more negative
        is tighter for the first, larger for the other two - so this
        returns everything in the "lower is better" convention the rest
        of the module already uses, by negating the CNN numbers.
        """
        if metric == "cnn_score":
            return None if self.cnn_score is None else -self.cnn_score
        if metric == "cnn_affinity":
            return None if self.cnn_affinity is None else -self.cnn_affinity
        return self.affinity

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def _parse_poses(sdf_text: str, stdout: str, engine: str) -> list:
    """Every pose's scores, in the order the program wrote them.

    Anchored on the SDF tag syntax - a property is `> <NAME>` alone on a
    line with its value on the next - because a loose search for
    `CNNaffinity` also matches `CNNaffinity_variance`, which gnina writes
    for every pose. That read nine poses as eighteen and interleaved each
    affinity with its own variance, so half of them came back near 0.18
    and looked like weak binders.
    """
    def tagged(field):
        pat = re.compile(
            r"^>\s*<" + re.escape(field) + r">[^\n]*\n\s*(-?\d+\.?\d*)",
            re.M)
        return [float(m.group(1)) for m in pat.finditer(sdf_text)]

    aff = tagged("minimizedAffinity")
    cnn = tagged("CNNscore")
    cnn_aff = tagged("CNNaffinity")
    if not aff:
        # Vina and older smina print a table instead of tagging the SDF.
        aff = [float(m.group(1)) for m in
               re.finditer(r"^\s+\d+\s+(-?\d+\.\d+)", stdout, re.M)]
    n = max(len(aff), len(cnn), len(cnn_aff))
    out = []
    for i in range(n):
        out.append(PoseScore(
            affinity=aff[i] if i < len(aff) else None,
            cnn_score=cnn[i] if i < len(cnn) else None,
            cnn_affinity=cnn_aff[i] if i < len(cnn_aff) else None,
            engine=engine))
    return out


def _run_dock(binary: str, args: list, out_sdf: Path, *,
              timeout: float = 1800) -> tuple:
    """Run a docking binary, on Windows or through WSL, and read it back.

    The binary is a WSL one when it is an absolute POSIX path, which is
    what `wsldock.find_binary` returns. Paths in `args` are translated
    for that side; on Windows they are passed through.
    """
    from . import wsldock as _wsl

    is_wsl = binary.startswith("/")
    if is_wsl:
        script = _wsl.dock_script(binary, args)
        try:
            proc = _wsl.run_script(script, timeout=timeout)
        except Exception as exc:                            # noqa: BLE001
            log.warning("WSL docking failed: %s", str(exc)[:160])
            return None, ""
    else:
        try:
            proc = subprocess.run([binary] + [str(a) for a in args],
                                  capture_output=True, text=True,
                                  timeout=timeout, encoding="utf-8", errors="replace")
        except subprocess.TimeoutExpired:
            log.warning("docking timed out")
            return None, ""
    if proc.returncode != 0:
        log.warning("docking failed: %s", (proc.stderr or "").strip()[:220])
        return None, ""
    text = out_sdf.read_text(errors="replace", encoding="utf-8") if out_sdf.exists() else ""
    return proc.stdout or "", text


def _dock_args(binary: str, receptor: Path, ligand: Path, out: Path,
               box: tuple, *, exhaustiveness: int, seed: int, cpu: int,
               num_modes: int = 9) -> list:
    """Command-line arguments, with paths on the side the binary runs."""
    from . import wsldock as _wsl

    conv = _wsl.to_wsl_path if binary.startswith("/") else (lambda p: str(p))
    cx, cy, cz, sx, sy, sz = box
    args = ["-r", conv(receptor), "-l", conv(ligand), "-o", conv(out),
            "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}",
            "--center_z", f"{cz:.3f}",
            "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}",
            "--size_z", f"{sz:.3f}",
            "--exhaustiveness", str(exhaustiveness),
            "--seed", str(seed), "--cpu", str(cpu),
            "--num_modes", str(num_modes)]
    if "gnina" in binary:
        # Rescore every pose, not only the top one, so the CNN can
        # reorder them - which is the entire reason to use gnina.
        args += ["--cnn_scoring", "rescore"]
    return args


def dock_one_scored(receptor_pdb: Path, ligand_sdf: Path, box: tuple, *,
                    exhaustiveness: int = 8, seed: int = 0, cpu: int = 1,
                    binary: str | None = None) -> PoseScore | None:
    """Dock one ligand and return every score the program produced."""
    binary = binary or _docking_binary()
    if binary is None:
        raise RuntimeError(
            "no docking binary found. " + SETUP)
    engine = ("gnina" if "gnina" in binary else
              "smina" if "smina" in binary else "vina")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "docked.sdf"
        args = _dock_args(binary, receptor_pdb, ligand_sdf, out, box,
                          exhaustiveness=exhaustiveness, seed=seed, cpu=cpu)
        stdout, sdf = _run_dock(binary, args, out)
        if stdout is None:
            return None
        poses = _parse_poses(sdf, stdout, engine)
        if not poses:
            return None
        # The program writes poses best-first on its own scoring, so the
        # first is its answer. Re-ranking on the CNN is the caller's
        # decision, not something to bury here.
        return poses[0]


# Organic counter-ions that RDKit's default salt list does not cover, as
# neutral SMILES. Taken from the pieces left over when the approved library is
# split. Anything that can also be the drug itself is left out on purpose:
# piperazine, lysine, caffeine, choline, hippuric acid and salicylic acid
# each appear as an active ingredient, and an entry such as piperazine
# citrate would otherwise lose the wrong half. Those entries stay whole and
# are reported as ambiguous.
_ORGANIC_COUNTER_IONS = (
    "O=C(O)CC(O)(CC(=O)O)C(=O)O",                 # citric acid
    "O=C(O)CCC(=O)O",                             # succinic acid
    "O=C(O)CC(O)C(=O)O",                          # malic acid
    "CC(O)C(=O)O",                                # lactic acid
    "O=C(O)c1ccccc1",                             # benzoic acid
    "O=C(O)CCCCC(=O)O",                           # adipic acid
    "NC(CC(=O)O)C(=O)O",                          # aspartic acid
    "OCC(O)C(O)C(O)C(O)C(=O)O",                   # gluconic acid
    "OCC(O)C(O)C(O)C(O)C(O)C(=O)O",               # glucoheptonic acid
    "O=C(O)C(O)C(O)C(O)C(O)C(=O)O",               # mucic and saccharic acid
    "O=C1CCC(C(=O)O)N1",                          # pidolic acid
    "O=S(=O)(O)c1ccccc1",                         # benzenesulfonic acid
    "O=S(=O)(O)c1ccc2ccccc2c1",                   # naphthalenesulfonic acid
    "O=S(=O)(O)CCO",                              # isethionic acid
    "COS(=O)(=O)O",                               # methyl sulfate
    "O=C(O)c1cc(S(=O)(=O)O)ccc1O",                # sulfosalicylic acid
    "O=[N+]([O-])c1cc([N+](=O)[O-])c(O)c([N+](=O)[O-])c1",   # picric acid
    "CNCC(O)C(O)C(O)C(O)CO",                      # meglumine
    "NC(CO)(CO)CO",                               # tromethamine
    "NCCO", "OCCNCCO", "OCCN1CCCC1",              # ethanolamines
    "CCO", "CC(O)CO",                             # ethanol, propylene glycol
    "O=C(O)c1cc2ccccc2c(Cc2c(O)c(C(=O)O)cc3ccccc23)c1O",      # pamoic acid
    "CCCCCCCCCCCCCCCCCC(=O)O",                    # stearic acid
    "CC1(C)C2CCC1(CS(=O)(=O)O)C(=O)C2",           # camphorsulfonic acid
    "CCS(=O)(=O)O", "O=S(=O)(O)CCS(=O)(=O)O",     # ethane(di)sulfonic acid
    "CCCCCCCCCCCCOS(=O)(=O)O",                    # lauryl sulfate
    "CC(C)(C)N",                                  # tert-butylamine (erbumine)
    "O=C(O)c1ccc(C(=O)O)cc1",                     # terephthalic acid
    "OCC(O)C(O)C(OC1OC(CO)C(O)C(O)C1O)C(O)C(=O)O",            # lactobionic acid
)
_COUNTER_ION_KEYS: set = set()


def _neutral_key(mol) -> str:
    """Charge-neutral, stereo-free SMILES, so an acid and its anion match."""
    from rdkit import Chem, rdBase
    from rdkit.Chem.MolStandardize import rdMolStandardize
    _quiet = rdBase.BlockLogs()     # the Uncharger announces itself each call
    try:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    except Exception as exc:                                # noqa: BLE001
        # Matching then uses the charged form, which can only miss a
        # counter-ion (the entry stays whole), never remove the wrong piece.
        log.debug("could not neutralise a fragment for matching: %s",
                  type(exc).__name__)
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def _is_organic_counter_ion(mol) -> bool:
    from rdkit import Chem
    if not _COUNTER_ION_KEYS:
        _COUNTER_ION_KEYS.update(
            _neutral_key(Chem.MolFromSmiles(x)) for x in _ORGANIC_COUNTER_IONS)
    return _neutral_key(mol) in _COUNTER_ION_KEYS


def split_counter_ions(smiles: str) -> dict:
    """Separate a drug from the counter-ions it is stored with.

    About a third of the approved library is written as a salt: the drug
    plus chloride, calcium, sulfate, maleate and so on, sometimes with two
    copies of the drug. Docked whole, the extra pieces are scored as part
    of the ligand, which is not what binds a receptor (the salt dissociates
    in solution). Returns

        smiles      the drug alone: one copy of each remaining component
        removed     the counter-ions taken out, one entry per distinct kind
        n_fragments how many pieces the input had
        ambiguous   True when more than one component is left

    A piece is a counter-ion when it has no carbon (chloride, a metal,
    sulfate, phosphate), RDKit's default salt list recognises it (acetate,
    maleate, fumarate, tartrate, mesylate and the like), or it is on the
    short list above (citrate, succinate, besylate, meglumine and so on).
    Anything else stays. Two organic components that are both left, such as a drug
    paired with piperazine or two active ingredients sold together, are
    reported as ambiguous and kept, because choosing the "real" one would be
    a guess about pharmacology. Charges are not changed. When every piece is
    a counter-ion (sodium chloride), the largest is kept so that something
    is still returned.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import SaltRemover

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"could not parse SMILES: {smiles}")
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) < 2:
        return {"smiles": smiles, "removed": [], "n_fragments": len(frags),
                "ambiguous": False}

    distinct: dict = {}
    for f in frags:
        distinct.setdefault(Chem.MolToSmiles(f), f)
    remover = SaltRemover.SaltRemover()
    keep, removed = [], []
    for smi, f in distinct.items():
        has_carbon = any(a.GetAtomicNum() == 6 for a in f.GetAtoms())
        stripped = remover.StripMol(f, dontRemoveEverything=False)
        is_counter_ion = ((not has_carbon) or stripped.GetNumAtoms() == 0
                          or _is_organic_counter_ion(f))
        (removed if is_counter_ion else keep).append(smi)
    if not keep:
        biggest = max(distinct, key=lambda k: distinct[k].GetNumHeavyAtoms())
        keep = [biggest]
        removed = [k for k in distinct if k != biggest]
    return {"smiles": ".".join(keep), "removed": removed,
            "n_fragments": len(frags), "ambiguous": len(keep) > 1}


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

    # `n_confs` was accepted and then ignored: EmbedMolecule embeds exactly
    # one whatever was asked for. For a flexible ligand the starting
    # conformer matters, so generate the requested number, minimise them
    # all and write the lowest-energy one. The default stays 1, so no
    # existing result changes unless a caller asks for more.
    if n_confs > 1:
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs,
                                              params=params))
        if not ids:
            raise ValueError("3D embedding failed")
        # MMFFOptimizeMoleculeConfs minimises every conformer in one call
        # and returns (converged, energy) per conformer, in order.
        try:
            ff_res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500,
                                                       numThreads=1)
        except Exception:                                   # noqa: BLE001
            ff_res = None
        energies = ([(e, cid) for (_conv, e), cid in zip(ff_res, ids)]
                    if ff_res and len(ff_res) == len(ids) else [])
        if energies:
            best = min(energies)[1]
            keep = [c for c in ids if c != best]
            for cid in keep:
                mol.RemoveConformer(cid)
        else:
            # No force field for this molecule: keep the first embedding
            # rather than picking arbitrarily on an unscored set.
            for cid in ids[1:]:
                mol.RemoveConformer(cid)
    else:
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
        conf = classify_conformation(path.read_text(encoding="utf-8", errors="replace"), pid)
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
    out.write_text(text, encoding="utf-8")
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
# AUC and enrichment compare the known ligands with the decoys, so the
# decoys that SCORED matter as much as the actives that did: five actives
# ranked against one surviving decoy is a coin flip with a number on it.
# The gate normally docks at least 50 decoys, so twenty scored means
# docking handled well under half of them failing. A floor on sample size,
# not a calibrated threshold.
MIN_VALIDATION_DECOYS = 20

# Above this the gate is sampled rather than run whole. Sixty actives
# against ten times as many decoys puts the standard error on the AUC at
# roughly 0.02, which is far finer than the distinction being made, and
# a target with thousands of measured ligands would otherwise spend more
# docking on the gate than on the library it is gating.
MAX_VALIDATION_ACTIVES = 60


def roc_auc(scores: list, labels: list, *, lower_is_better: bool = True
            ) -> float | None:
    """AUC via Mann-Whitney U, ties at half. None entries are failed docks
    and are dropped rather than scored worst.

    Dropping them is the usual choice and it is NOT neutral: it measures
    recall only over the compounds that docked, so a method that fails on
    most of the known ligands is judged on the few it could handle, and
    looks better than it is. The caller has to bound how many were
    dropped; `validate_pocket` does, and reports the counts."""
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


def enrichment_basis(scores: list, labels: list, *, fraction: float = 0.01,
                     lower_is_better: bool = True) -> dict | None:
    """What the enrichment factor rests on: how many compounds, how many hits.

    With a few hundred compounds the top 1% is a handful, so the factor
    moves in coarse steps (one active more or fewer among seven shifts it
    by about 1.7) while being quoted to two decimals. Two runs of the same
    pocket gave 11.31 and 9.62, which is one compound. Stating the basis
    lets a reader see that, instead of reading the difference as real.
    """
    pairs = [(s, y) for s, y in zip(scores, labels) if s is not None]
    if not pairs:
        return None
    n = len(pairs)
    n_act = sum(1 for _, y in pairs if y)
    if not n_act or n_act == n:
        return None
    k = max(1, int(math.ceil(n * fraction)))
    pairs.sort(key=lambda t: t[0] if lower_is_better else -t[0])
    return {"top_n": k, "actives_in_top_n": sum(1 for _, y in pairs[:k] if y),
            "n_scored": n, "n_actives_scored": n_act}


# The receptor PDBQT is the same file for every compound in a
# screen, so it is converted once and reused.
_REC_CACHE: dict = {}


def _dock_smiles(smiles_list: list, receptor: Path, autobox: Path | None,
                 *, exhaustiveness: int, seed: int, score_fn=None,
                 box: tuple | None = None, cpu: int = 1,
                 jobs: int = 1, detailed: bool = False,
                 tally: dict | None = None,
                 strip_counter_ions: bool = False) -> list:
    """Score each SMILES, None where docking or embedding failed.

    `strip_counter_ions` docks the drug without its salt (see
    `split_counter_ions`). It is off by default, which docks the structure
    exactly as stored, so no earlier result changes.

    `score_fn` exists so the gate and the ranking can be exercised
    without a docking binary; it takes a SMILES and returns a score.

    `detailed` returns a PoseScore per compound rather than one float,
    keeping gnina's CNN numbers. Without it they are computed and thrown
    away, which is the whole reason to run gnina rather than smina.

    `tally`, when given, is filled with the counts of what failed and
    why: a compound that could not be given a 3D structure is a
    different problem from one whose docking crashed, and reporting
    either as the other sends you after the wrong thing.

    `jobs` runs that many compounds at once. Each is a separate
    subprocess with its own temporary directory, so they share nothing
    but the read-only receptor, and the returned list stays in the order
    of `smiles_list` whatever order they finish in. Scores do not depend
    on it: same receptor, same box, same seed.
    """
    n = len(smiles_list)
    tally_lock = threading.Lock()

    def _count(kind):
        if tally is None:
            return
        with tally_lock:
            tally[kind] = tally.get(kind, 0) + 1

    def one(smi):
        """Score a single SMILES. Returns None rather than raising."""
        if strip_counter_ions:
            try:
                smi = split_counter_ions(smi)["smiles"]
            except ValueError:
                pass            # unparseable: fails below, and is counted
        if score_fn is not None:
            try:
                got = score_fn(smi)
            except Exception as exc:                        # noqa: BLE001
                log.debug("scorer raised for %s: %s", smi[:40],
                          type(exc).__name__)
                _count("scorer_failed")
                return None
            # The gate and the ranking are exercised without a binary by
            # passing score_fn, so it has to work in both shapes.
            if detailed and got is not None and not isinstance(got,
                                                               PoseScore):
                return PoseScore(affinity=float(got), engine="score_fn")
            return got
        with tempfile.TemporaryDirectory() as tmp:
            sdf = Path(tmp) / "lig.sdf"
            try:
                smiles_to_sdf(smi, sdf, seed=seed)
            except Exception as exc:                        # noqa: BLE001
                log.debug("embedding failed for %s: %s", smi[:40], exc)
                _count("embedding_failed")
                return None
            try:
                if box is not None and engine_name() == "vina":
                    return dock_one_vina(
                        receptor, smi, box,
                        exhaustiveness=exhaustiveness, seed=seed,
                        cpu=cpu, receptor_pdbqt_cache=_REC_CACHE)
                if box is not None:
                    if detailed:
                        got = dock_one_scored(
                            receptor, sdf, box,
                            exhaustiveness=exhaustiveness, seed=seed,
                            cpu=cpu)
                    else:
                        got = dock_one_boxed(
                            receptor, sdf, box,
                            exhaustiveness=exhaustiveness, seed=seed,
                            cpu=cpu)
                    if got is None:
                        _count("no_pose")
                    return got
                got = dock_one(
                    receptor, sdf, autobox_ligand=autobox,
                    exhaustiveness=exhaustiveness, seed=seed, cpu=cpu)
                if got is None:
                    _count("no_pose")
                return got
            except Exception as exc:                        # noqa: BLE001
                log.debug("docking failed: %s", str(exc)[:80])
                _count("docking_failed")
                return None

    # A screen used to print nothing until it finished, so a working run
    # and a wedged one looked identical for hours. The rate is reported
    # too: it is what tells you whether to wait.
    t0 = time.time()
    step = max(25, n // 20)
    done = 0
    lock = threading.Lock()

    def tick():
        nonlocal done
        with lock:
            done += 1
            d = done
        if d % step == 0 and d < n:
            rate = (time.time() - t0) / d
            log.info("docking %d/%d (%.0f%%), %.1fs each, about %d min "
                     "left", d, n, 100.0 * d / n, rate,
                     int((n - d) * rate // 60))

    def scored(smi):
        try:
            return one(smi)
        finally:
            tick()

    if jobs <= 1:
        return [scored(smi) for smi in smiles_list]

    # map keeps the results in argument order however they complete.
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(scored, smiles_list))


def validate_pocket(receptor: Path, actives: list, decoys: list, *,
                    autobox: Path | None = None, exhaustiveness: int = 8,
                    seed: int = 0, score_fn=None,
                    box: tuple | None = None, cpu: int = 1,
                    jobs: int = 1, strip_counter_ions: bool = False) -> dict:
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
    poses = _dock_smiles(smiles, receptor, autobox,
                         jobs=jobs, detailed=True,
                         exhaustiveness=exhaustiveness, seed=seed,
                         score_fn=score_fn, box=box, cpu=cpu,
                         strip_counter_ions=strip_counter_ions)
    scores = [None if q is None else q.affinity for q in poses]

    # What each metric is worth on THIS pocket, measured rather than
    # assumed. The poses are already scored, so it costs nothing, and
    # known ligands against decoys is the only honest way to find out
    # whether the CNN separates binders better than the empirical
    # function here. Reported, not acted on: the gate still decides on
    # the affinity, so there is a fixed thing to compare against.
    by_metric = {}
    for metric in ("affinity", "cnn_score", "cnn_affinity"):
        vals = [None if q is None else q.best(metric) for q in poses]
        if all(v is None for v in vals):
            continue
        m_auc = roc_auc(vals, labels)
        if m_auc is None:
            continue
        by_metric[metric] = {
            "auc": round(m_auc, 3),
            "ef1": enrichment_factor(vals, labels, fraction=0.01),
        }

    auc = roc_auc(scores, labels)
    ef1 = enrichment_factor(scores, labels, fraction=0.01)
    n_failed = sum(1 for s in scores if s is None)

    # How many of each class actually scored. `n_failed` mixes actives and
    # decoys, and the minimum-actives check above ran on what was SUPPLIED,
    # before docking: a pocket where docking failed for 58 of 60 known
    # ligands passed with an AUC of 1.0 from the two that scored. What the
    # gate measures is recall on known ligands, and a ligand that failed to
    # dock is dropped from that measurement rather than counted as a miss.
    n_act = len(actives)
    a_scored = sum(1 for s in scores[:n_act] if s is not None)
    d_scored = sum(1 for s in scores[n_act:] if s is not None)
    gate = {"n_actives": n_act, "n_actives_scored": a_scored,
            "n_actives_failed": n_act - a_scored,
            "n_decoys": len(decoys), "n_decoys_scored": d_scored,
            "n_decoys_failed": len(decoys) - d_scored,
            "ef1_basis": enrichment_basis(scores, labels, fraction=0.01)}

    if auc is None:
        return {"usable": False, "auc": None, "ef1": None,
                "n_failed": n_failed, "gate": gate,
                "reason": "every dock failed, so nothing could be measured"}
    if a_scored < MIN_VALIDATION_ACTIVES:
        return {"usable": False, "auc": None, "ef1": None,
                "n_failed": n_failed, "gate": gate,
                "reason": (f"docking produced a score for only {a_scored} "
                           f"of {n_act} known ligand(s), below the "
                           f"{MIN_VALIDATION_ACTIVES} needed to say whether "
                           f"docking works in this pocket. The other "
                           f"{n_act - a_scored} failed (3D embedding or "
                           f"pose generation), and an AUC taken from "
                           f"{a_scored} compound(s) would describe those "
                           f"alone. Not reported.")}
    if d_scored < MIN_VALIDATION_DECOYS:
        return {"usable": False, "auc": None, "ef1": None,
                "n_failed": n_failed, "gate": gate,
                "reason": (f"docking produced a score for only {d_scored} "
                           f"of {len(decoys)} decoy(s), below the "
                           f"{MIN_VALIDATION_DECOYS} needed to rank the "
                           f"known ligands against. The other "
                           f"{len(decoys) - d_scored} failed, and an AUC "
                           f"against {d_scored} decoy(s) would be mostly "
                           f"noise. Not reported.")}

    usable = auc >= MIN_DOCKING_AUC and (ef1 is None or ef1 >= MIN_EF1)
    out = {"usable": usable, "auc": round(auc, 3), "ef1": ef1,
           "n_actives": len(actives), "n_decoys": len(decoys),
           "n_failed": n_failed, "gate": gate}
    if n_act - a_scored:
        out["gate_note"] = (
            f"The AUC and enrichment cover the {a_scored} of {n_act} known "
            f"ligands that docked. The {n_act - a_scored} that failed are "
            f"left out, not counted as misses, so the figures describe the "
            f"compounds docking can handle.")
    if len(by_metric) > 1:
        out["auc_by_metric"] = by_metric
        best = max(by_metric, key=lambda k: by_metric[k]["auc"])
        if best != "affinity":
            out["metric_note"] = (
                f"{best} separates this pocket's known ligands better "
                f"than the empirical affinity "
                f"(AUC {by_metric[best]['auc']} against "
                f"{by_metric['affinity']['auc']}). Ranking still uses the "
                f"affinity; this is a measurement on {len(actives)} known "
                f"ligand(s) of one pocket, not a general finding.")
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
    # What each scoring metric was worth on this pocket's own known
    # ligands, when the engine produced more than one.
    auc_by_metric: dict | None = None
    # Why each absent compound is absent, counted rather than assumed.
    failures: dict | None = None
    # Whether agonist/antagonist was predicted for the hits, and on what.
    direction: dict | None = None
    # What the validation gate rested on: how many known ligands and decoys
    # actually scored, and what the top-1% enrichment was counted over.
    gate: dict | None = None
    # How ligands were prepared, when that differs from "as stored": set
    # only when counter-ions were stripped, with how many entries changed.
    ligand_prep: dict | None = None
    hits: list = field(default_factory=list)
    error: str | None = None
    caveats: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def _model_assessment(receptor) -> "PocketAssessment | None":
    """Assess `receptor` when it is a predicted model, else None.

    An experimental structure's B-factor column holds crystallographic B
    factors, not pLDDT, so it must not be read as confidence. AlphaFold
    files are recognised by their own markers rather than by filename:
    every pLDDT sits in 0-100 and the header names the predictor.
    """
    try:
        text = Path(receptor).read_text(errors="replace", encoding="utf-8")
    except OSError:
        return None
    head = text[:4000].upper()
    predicted = ("ALPHAFOLD" in head or "DEEPMIND" in head
                 or "PREDICTED" in head or "ESMFOLD" in head
                 or "ROSETTAFOLD" in head)
    if not predicted:
        return None
    return assess_model(text, uniprot=Path(receptor).stem)


def screen_by_docking(receptor: Path, library: list, *, target: str = "?",
                      actives: list | None = None, decoys: list | None = None,
                      autobox: Path | None = None, top_n: int = 25,
                      exhaustiveness: int = 8, seed: int = 0,
                      score_fn=None, box: tuple | None = None,
                      cpu: int = 1, jobs: int = 1,
                      strip_counter_ions: bool = False) -> DockingScreen:
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

    # A predicted model is checked before anything is docked into it.
    # `assess_model` was written to make this call - "low-pLDDT regions
    # look like open greasy cavities to a scoring function" - and nothing
    # called it, so an AlphaFold model of a transcription factor could be
    # docked into its disordered activation domain and score beautifully.
    model = _model_assessment(receptor)
    if model is not None:
        run.caveats.append(
            f"predicted structure: mean pLDDT {model.mean_plddt}, "
            f"{model.ordered_fraction:.0%} of residues confidently folded")
        run.caveats.extend(model.notes)
        if not model.dockable:
            run.error = (
                f"{receptor.name} is a predicted model with only "
                f"{model.ordered_fraction:.0%} of residues above pLDDT "
                f"{PLDDT_ORDERED:.0f}. That is a disorder prediction, not a "
                f"structure with pockets, and docking into it returns "
                f"confident numbers that mean nothing. Refused. Use an "
                f"experimental structure, or restrict the receptor to its "
                f"folded domain first.")
            return run

    val = validate_pocket(receptor, actives or [], decoys or [],
                          jobs=jobs,
                          autobox=autobox, exhaustiveness=exhaustiveness,
                          seed=seed, score_fn=score_fn, box=box, cpu=cpu,
                          strip_counter_ions=strip_counter_ions)
    run.auc, run.ef1, run.usable = val.get("auc"), val.get("ef1"), \
        bool(val.get("usable"))
    run.gate = val.get("gate")
    if val.get("gate_note"):
        run.caveats.append(val["gate_note"])
    if not run.usable:
        run.error = val.get("reason")
        return run

    tally: dict = {}
    poses = _dock_smiles([d["smiles"] for d in library], receptor, autobox,
                         jobs=jobs, detailed=True, tally=tally,
                         exhaustiveness=exhaustiveness, seed=seed,
                         score_fn=score_fn, box=box, cpu=cpu,
                         strip_counter_ions=strip_counter_ions)
    if strip_counter_ions:
        changed = ambiguous = 0
        for d in library:
            try:
                got = split_counter_ions(d["smiles"])
            except ValueError:
                continue
            changed += bool(got["removed"])
            ambiguous += got["ambiguous"]
        run.ligand_prep = {"strip_counter_ions": True,
                           "n_library": len(library),
                           "n_changed": changed, "n_ambiguous": ambiguous}
    scored = []
    for d, q in zip(library, poses):
        if q is None or q.affinity is None:
            continue
        row = {"drug": d.get("name") or "?", "smiles": d["smiles"],
               "score": q.affinity}
        # Only present when the engine produced them, so a smina run
        # does not carry empty CNN columns that look like a measurement
        # of zero.
        if q.cnn_score is not None:
            row["cnn_score"] = q.cnn_score
        if q.cnn_affinity is not None:
            row["cnn_affinity"] = q.cnn_affinity
        if q.engine and q.engine != "?":
            row["engine"] = q.engine
        scored.append(row)
    scored.sort(key=lambda h: h["score"])
    run.n_docked = len(scored)
    run.n_failed = len(library) - len(scored)
    run.auc_by_metric = val.get("auc_by_metric")
    if val.get("metric_note"):
        run.caveats.append(val["metric_note"])
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
        run.failures = dict(tally)
        # What actually happened, not what usually happens. These are
        # different problems: one is a limit of conformer generation,
        # the other is the docking program failing on input it was
        # given, and they send you to different places.
        parts = []
        if tally.get("embedding_failed"):
            parts.append(f"{tally['embedding_failed']} had no 3D conformer "
                         f"that could be generated, which is typical of "
                         f"macrocycles and peptides")
        if tally.get("docking_failed"):
            parts.append(f"{tally['docking_failed']} raised an error in "
                         f"the docking program")
        if tally.get("scorer_failed"):
            parts.append(f"{tally['scorer_failed']} raised an error in the "
                         f"scoring function")
        if tally.get("no_pose"):
            parts.append(f"{tally['no_pose']} were docked but produced no "
                         f"pose")
        counted = sum(tally.values())
        if counted < run.n_failed:
            parts.append(f"{run.n_failed - counted} for an unrecorded "
                         f"reason")
        run.caveats.append(
            f"{run.n_failed} compound(s) are absent from the ranking "
            f"rather than ranked badly: "
            + ("; ".join(parts) if parts else "cause not recorded") + ".")
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
    # A dict is accepted as well as the dataclass, because the docking
    # run usually happens on another machine and comes back as JSON. The
    # typed object never crosses that boundary, so requiring it meant
    # this function could not be called on a real split run at all.
    if hasattr(docking, "usable"):
        usable, d_hits = docking.usable, docking.hits
    elif isinstance(docking, dict):
        usable, d_hits = docking.get("usable"), docking.get("hits") or []
    else:
        return []
    if not usable or not ligand_hits:
        return []
    n_d = max(1, int(len(d_hits) * top_frac)) or 1
    dock_top = {(h.get("drug") or "").strip().lower(): h
                for h in d_hits[:n_d] if h.get("drug")}
    n_l = max(1, int(len(ligand_hits) * top_frac)) or 1

    out = []
    for h in ligand_hits[:n_l]:
        name = (h.get("drug") or "").strip().lower()
        if name in dock_top:
            out.append({
                "drug": h.get("drug"),
                "docking_score": dock_top[name].get("score"),
                "ligand_quadrant": h.get("quadrant"),
                "why": ("picked by both a pocket-shape argument and a "
                        "chemistry argument trained on this receptor's "
                        "measured ligands. The two share no input, so this "
                        "is the strongest computational support available "
                        "here - and still a binding hypothesis, with the "
                        "direction unaddressed by either"),
            })
    return out


# ===================================================================== #
# Preparing a receptor and finding its pocket.                          #
# ===================================================================== #
# dock_one can place a box from a crystallographic ligand, but nothing
# produced that ligand file, and the receptor was handed to the docking
# binary with the crystal ligand still sitting in the pocket - which means
# docking into an occupied site. These fill both gaps.

# Crystallisation additives, buffers and ions. Present in most structures,
# never the thing that defines the orthosteric pocket.
SOLVENTS = frozenset({
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "1PE",
    "MPD", "ACT", "ACY", "FMT", "TRS", "DMS", "IOD", "NO3", "EPE", "MES",
    "CIT", "TLA", "BME", "IMD", "CAC", "BCT", "AZI", "SCN", "FLC",
    "CL", "BR", "NA", "K", "MG", "ZN", "CA", "MN", "FE", "NI", "CU", "CD",
    "HG", "CO", "SR", "CS", "LI", "RB", "BA",
})

# Below this a HETATM group is an ion or a fragment rather than a ligand
# that marks out a drug-sized pocket.
MIN_LIGAND_ATOMS = 10

# Modified amino acids. These are part of the polypeptide and are written
# as HETATM only because they are non-standard, so a box placed on one
# sits on the protein rather than in a pocket. CCR5 found this the
# expensive way: its N-terminal structures carry sulfated tyrosines, and
# TYS was selected as the co-crystallised ligand.
#
# Hand-maintained and certainly incomplete. A residue missing from here
# is treated as a ligand, which is the failure this exists to catch, so
# it is worth extending when a screen chooses something odd.
MODIFIED_RESIDUES = frozenset({
    "MSE", "SEC", "PYL",                      # selenium, rare encodings
    "TYS", "PTR", "SEP", "TPO", "NEP",        # sulfation, phosphorylation
    "CSO", "CSD", "CSX", "OCS", "CME", "SMC", "SNC", "CAS", "CSS",
    "KCX", "LLP", "MLY", "M3L", "MLZ", "ALY", "DM0",
    "HYP", "PCA", "FME", "CGU", "HIC", "MHO", "TPQ", "DAH", "ABA",
    "AIB", "SAC", "ORN", "DAL", "DLE", "DVA", "DPR", "MVA", "NLE",
    "PHD", "AGM", "2MR", "DDE", "SCY", "YCM", "CZZ", "MGN",
})


# Lipids, detergents and cryoprotectant fragments that come from the
# crystallisation rather than from biology. A membrane protein solved in
# a lipidic cubic phase is full of monoolein, and it lies against the
# outside of the protein, so a box placed on one screens the membrane
# face instead of the pocket.
#
# NOT excluded, only ranked last. For a lipid receptor the fatty acid is
# the real ligand, and a blanket rule would throw away the structure that
# answers the question. Where nothing else is present the choice is still
# made and flagged.
CRYSTALLISATION_ADDITIVES = frozenset({
    "OLC", "OLA", "OLB", "MOH", "D10", "D12", "DD9", "R16", "PLM",
    "STE", "MYR", "DAO", "LDA", "LMT", "LMU", "BNG", "BOG", "OCT",
    "C8E", "9BA", "12P", "PE4", "PEU", "P6G", "2PE", "7PE", "15P",
    "CLR", "CHD", "CHS", "Y01", "HC3",
    "PC1", "PEE", "PGW", "PGV", "POV", "PLC", "LHG", "DGA", "3PE",
    "SQU", "TRD", "UND", "HEX", "HEZ", "DKA", "LI1", "HTG", "SOG",
})


def hetatm_groups(pdb_text: str) -> list:
    """Every non-solvent HETATM group, largest first.

    Each entry is {resname, chain, resi, n_atoms, centre}. The largest is
    almost always the co-crystallised ligand, which is the most reliable
    marker of where the orthosteric pocket is.
    """
    groups: dict = {}
    for line in pdb_text.splitlines():
        # An NMR entry holds many models of the same atoms. Keyed on
        # (resname, chain, resi) they all pool into one group, so a
        # sixteen-atom residue in a ten-model file counts as a
        # hundred and sixty atoms and its centre is averaged across ten
        # different conformations - a point no model occupies. Since the
        # largest group is taken as the ligand, that makes every
        # multi-model structure outrank a real crystal ligand.
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("HETATM"):
            continue
        resn = line[17:20].strip().upper()
        if resn in SOLVENTS or resn in MODIFIED_RESIDUES:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element == "H":
            continue
        key = (resn, line[21], line[22:27].strip())
        try:
            xyz = (float(line[30:38]), float(line[38:46]),
                   float(line[46:54]))
        except ValueError:
            continue
        groups.setdefault(key, []).append(xyz)

    out = []
    for (resn, chain, resi), coords in groups.items():
        n = len(coords)
        out.append({
            "resname": resn, "chain": chain, "resi": resi, "n_atoms": n,
            "centre": tuple(sum(c[i] for c in coords) / n for i in range(3)),
        })
    out.sort(key=lambda g: -g["n_atoms"])
    return out


def best_ligand(pdb_text: str) -> dict | None:
    """The HETATM group most likely to mark the orthosteric pocket.

    A real ligand beats a crystallisation additive of any size. Monoolein
    in a GPCR structure is routinely the largest group in the file and
    lies against the membrane face, so taking the largest outright
    screens the outside of the protein.

    An additive is still returned when it is all there is, carrying
    `additive: True`, because for a lipid receptor it may be the genuine
    ligand. The caller decides what that is worth; it is not silently
    presented as an orthosteric pocket.
    """
    groups = [g for g in hetatm_groups(pdb_text)
              if g["n_atoms"] >= MIN_LIGAND_ATOMS]
    for g in groups:
        if g["resname"] not in CRYSTALLISATION_ADDITIVES:
            g["additive"] = False
            return g
    if groups:
        groups[0]["additive"] = True
        return groups[0]
    return None


def write_ligand_pdb(pdb_text: str, group: dict, path) -> Path:
    """The chosen ligand on its own, for placing the docking box."""
    path = Path(path)
    keep = []
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        if (line[17:20].strip().upper() == group["resname"]
                and line[21] == group["chain"]
                and line[22:27].strip() == group["resi"]):
            keep.append(line)
    path.write_text("\n".join(keep) + "\nEND\n", encoding="utf-8")
    return path


def prepare_receptor(pdb_text: str, path, *, keep_chain: str | None = None
                     ) -> Path:
    """Protein only: HETATMs and waters removed, one chain if asked.

    Docking into a structure with its crystal ligand still in place means
    docking into an occupied pocket, and every score comes back weak for
    reasons that have nothing to do with the compounds.
    """
    path = Path(path)
    keep = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            if keep_chain and line[21] != keep_chain:
                continue
            keep.append(line)
        elif line.startswith("TER"):
            keep.append(line)
    path.write_text("\n".join(keep) + "\nEND\n", encoding="utf-8")
    return path


def box_from_ligand(ligand_pdb_text: str, padding: float = 4.0) -> tuple:
    """(cx, cy, cz, sx, sy, sz) around a ligand, in angstroms.

    Computed here rather than left to the docking binary because
    `--autobox_ligand` is a smina and gnina flag. AutoDock Vina has no
    equivalent and needs the centre and size given explicitly, so a
    command built for smina fails on the binary most people install.
    """
    xs, ys, zs = [], [], []
    for line in ligand_pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            xs.append(float(line[30:38]))
            ys.append(float(line[38:46]))
            zs.append(float(line[46:54]))
        except ValueError:
            continue
    if not xs:
        raise ValueError("no coordinates in the ligand file")
    centre = [(min(v) + max(v)) / 2.0 for v in (xs, ys, zs)]
    size = [max(max(v) - min(v) + 2 * padding, 12.0) for v in (xs, ys, zs)]
    return tuple(centre) + tuple(size)


class _NeedsVinaPath(RuntimeError):
    """Raised when the binary is Vina, so the caller uses the Vina path."""


def engine_name(binary: str | None = None) -> str:
    """Which docking program the binary on PATH is."""
    b = (binary or _docking_binary() or "").lower()
    for n in ("gnina", "smina", "vina"):
        if n in b:
            return n
    return "unknown"


def dock_one_boxed(receptor_pdb: Path, ligand_sdf: Path, box: tuple, *,
                   exhaustiveness: int = 8, seed: int = 0,
                   cpu: int = 1) -> float | None:
    """Dock into an explicit box, so the command works on any of the three.

    smina and gnina accept a PDB receptor and an SDF ligand and do the
    atom typing themselves. AutoDock Vina wants PDBQT for both, which is
    a separate preparation step, so it is detected and refused here with
    an instruction rather than failing inside the binary.

    Returns the empirical affinity, for compatibility with every caller
    that ranks on it. `dock_one_scored` returns the CNN numbers too.
    """
    binary = _docking_binary()
    if binary is None:
        raise RuntimeError(
            "no docking binary found. " + SETUP)
    eng = engine_name(binary)
    if eng == "vina":
        # Vina wants PDBQT for both files. That conversion is available
        # on Windows without WSL, so this is no longer a refusal: the
        # ligand goes through Meeko and the receptor through Meeko or
        # Open Babel.
        raise _NeedsVinaPath(
            "AutoDock Vina needs PDBQT input; use dock_one_vina, which "
            "prepares it.")

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "docked.sdf"
        args = _dock_args(binary, receptor_pdb, ligand_sdf, out, box,
                          exhaustiveness=exhaustiveness, seed=seed,
                          cpu=cpu)
        stdout, sdf = _run_dock(binary, args, out)
        if stdout is None:
            return None
        poses = _parse_poses(sdf, stdout, eng)
        affinities = [q.affinity for q in poses if q.affinity is not None]
        return min(affinities) if affinities else None


# ===================================================================== #
# Docking on Windows without WSL.                                       #
# ===================================================================== #
# smina has no Windows build, which is why WSL was the recommendation.
# AutoDock Vina does ship a Windows executable, and the only thing that
# stopped it being usable here was file preparation: smina takes a PDB
# receptor and an SDF ligand and types the atoms itself, while Vina wants
# PDBQT for both.
#
# Both conversions are available on Windows without WSL:
#
#   ligand    Meeko. Pure Python on top of RDKit, pip-installable, and the
#             tool the AutoDock developers point at for this.
#   receptor  Meeko's receptor preparation where the installed version has
#             it, otherwise Open Babel, which has Windows wheels and does
#             rigid-receptor PDBQT with `-xr`.
#
# So the Windows path is: vina.exe on PATH, plus `pip install meeko` and
# optionally `pip install openbabel-wheel`. Nothing else changes - the
# same box, the same validation gate, the same output.

VINA_SETUP_WINDOWS = """\
Docking on Windows, no WSL needed.

1. AutoDock Vina. Download the Windows executable from
   https://github.com/ccsb-scripps/AutoDock-Vina/releases
   and put vina.exe somewhere on your PATH.

2. File preparation, both pip-installable:

     pip install meeko              (ligands, and receptors on newer meeko)
     pip install openbabel-wheel    (receptor fallback)

Check it with:  python -m promoter_tf_drugs.docking --check

smina is still the better engine where you can run it - it takes PDB and
SDF directly and needs no conversion step - but it has no Windows build.
Vina with Meeko gets the same job done natively.
"""


def _meeko_ligand_pdbqt(smiles: str, seed: int = 0) -> str | None:
    """SMILES to a ligand PDBQT string via Meeko, or None."""
    try:
        from meeko import MoleculePreparation
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem
    except ImportError:
        return None
    try:
        from meeko import PDBQTWriterLegacy as _W
    except ImportError:
        _W = None

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return None
    if _W is not None:
        text, ok, err = _W.write_string(setups[0])
        if not ok:
            log.debug("meeko ligand write failed: %s", str(err)[:80])
            return None
        return text
    # Older meeko exposed the writer on the preparation object.
    try:
        return setups[0].write_pdbqt_string()
    except Exception as exc:                                # noqa: BLE001
        log.debug("meeko could not write the ligand: %s",
                  type(exc).__name__)
        return None


def _receptor_pdbqt(receptor_pdb: Path, out: Path) -> Path | None:
    """Rigid-receptor PDBQT, by whichever converter is installed."""
    # Meeko first: it is already a dependency of the ligand path.
    try:
        import subprocess as _sp
        import sys as _sys
        for mod in ("meeko.cli.mk_prepare_receptor", "mk_prepare_receptor"):
            try:
                proc = _sp.run(
                    [_sys.executable, "-m", mod, "--read_pdb",
                     str(receptor_pdb), "-o", str(out.with_suffix("")),
                     "--skip_gpf"],
                    capture_output=True, text=True, encoding="utf-8",
                    errors="replace", timeout=300)
            except FileNotFoundError:
                continue
            if proc.returncode == 0 and out.exists():
                return out
    except Exception as exc:                                # noqa: BLE001
        # The next converter (Open Babel) is tried; if all fail the caller
        # reports that no receptor file could be made.
        log.info("receptor conversion with the Python tools failed "
                 "(%s: %s); trying Open Babel", type(exc).__name__,
                 str(exc)[:100])

    # Open Babel: -xr writes a rigid receptor.
    obabel = shutil.which("obabel")
    if obabel is None:
        try:
            from openbabel import pybel          # noqa: F401
            obabel = "python-openbabel"
        except ImportError:
            obabel = None
    if obabel == "python-openbabel":
        try:
            from openbabel import pybel
            m = next(pybel.readfile("pdb", str(receptor_pdb)))
            m.OBMol.AddPolarHydrogens()
            m.write("pdbqt", str(out), overwrite=True)
            return out if out.exists() else None
        except Exception as exc:                            # noqa: BLE001
            log.debug("pybel receptor prep failed: %s", str(exc)[:80])
            return None
    if obabel:
        try:
            proc = subprocess.run(
                [obabel, "-ipdb", str(receptor_pdb), "-opdbqt",
                 "-O", str(out), "-xr"],
                capture_output=True, text=True, timeout=300, encoding="utf-8", errors="replace")
            if proc.returncode == 0 and out.exists():
                return out
        except Exception as exc:                            # noqa: BLE001
            log.debug("obabel receptor prep failed: %s", str(exc)[:80])
    return None


def vina_prep_available() -> tuple:
    """(ok, why_not) for the Windows Vina preparation chain."""
    try:
        import meeko                                        # noqa: F401
    except ImportError:
        return False, ("meeko is needed to turn SMILES into Vina ligand "
                       "files: pip install meeko")
    have_receptor = bool(shutil.which("obabel"))
    if not have_receptor:
        try:
            from openbabel import pybel                     # noqa: F401
            have_receptor = True
        except ImportError:
            pass
    if not have_receptor:
        # Newer meeko can do the receptor too; assume it can and let the
        # run report a failure rather than refusing up front.
        return True, ("no Open Babel found; relying on meeko for the "
                      "receptor. If that fails, pip install openbabel-wheel")
    return True, None


def dock_one_vina(receptor_pdb: Path, smiles: str, box: tuple, *,
                  exhaustiveness: int = 8, seed: int = 0, cpu: int = 1,
                  receptor_pdbqt_cache: dict | None = None) -> float | None:
    """One Vina dock from a PDB receptor and a SMILES, on any platform.

    The receptor conversion is cached across calls: it is the same file for
    every compound in a screen, and converting it several thousand times
    would dominate the run.
    """
    binary = shutil.which("vina") or shutil.which("vina.exe")
    if binary is None:
        # A vina inside WSL is still a vina. The PDBQT files are written
        # on the Windows side and WSL reads them under /mnt/c, so the
        # only difference is how the binary is invoked.
        from . import wsldock as _wsl
        try:
            _, binary = _wsl.find_binary(prefer="vina")
        except Exception:                                   # noqa: BLE001
            binary = None
        if binary is None or "vina" not in binary:
            raise RuntimeError("vina not found on PATH. "
                               + VINA_SETUP_WINDOWS)

    cache = receptor_pdbqt_cache if receptor_pdbqt_cache is not None else {}
    key = str(receptor_pdb)
    rec_pdbqt = cache.get(key)
    if rec_pdbqt is None or not Path(rec_pdbqt).exists():
        out = Path(receptor_pdb).with_suffix(".pdbqt")
        got = _receptor_pdbqt(Path(receptor_pdb), out)
        if got is None:
            raise RuntimeError(
                "could not convert the receptor to PDBQT. "
                + VINA_SETUP_WINDOWS)
        cache[key] = str(got)
        rec_pdbqt = str(got)

    lig = _meeko_ligand_pdbqt(smiles, seed=seed)
    if lig is None:
        return None

    cx, cy, cz, sx, sy, sz = box
    with tempfile.TemporaryDirectory() as tmp:
        lp = Path(tmp) / "lig.pdbqt"
        lp.write_text(lig, encoding="utf-8")
        op = Path(tmp) / "out.pdbqt"
        from . import wsldock as _wsl
        conv = (_wsl.to_wsl_path if binary.startswith("/")
                else (lambda q: str(q)))
        args = ["--receptor", conv(rec_pdbqt), "--ligand", conv(lp),
                "--out", conv(op),
                "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}",
                "--center_z", f"{cz:.3f}",
                "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}",
                "--size_z", f"{sz:.3f}",
                "--exhaustiveness", str(exhaustiveness),
                "--seed", str(seed), "--cpu", str(cpu)]
        stdout, _unused = _run_dock(binary, args, op)
        if stdout is None:
            return None

        text = op.read_text(encoding="utf-8", errors="replace") if op.exists() else ""
        scores = [float(m.group(1)) for m in
                  re.finditer(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", text)]
        if not scores:
            scores = [float(m.group(1)) for m in
                      re.finditer(r"^\s+1\s+(-?\d+\.\d+)", stdout, re.M)]
        return min(scores) if scores else None


# ===================================================================== #
# The orchestrator: a gene symbol to a ranked shortlist.                #
# ===================================================================== #
SETUP = """\
To run docking you need two things.

1. A docking binary. There are two tiers, and the difference is real.

   gnina - the better answer
     It docks the way smina does and then RESCORES every pose with a
     convolutional network trained on crystal structures, which is a
     different kind of judgement about whether a pose looks real. You
     get a CNNscore (0-1, how pose-like the geometry is) and a
     CNNaffinity (predicted pK, comparable across ligands in a way a
     Vina score is not) alongside the usual empirical number.

     Linux, or Windows through WSL:
       curl -L -o gnina https://github.com/gnina/gnina/releases/latest/download/gnina
       chmod +x gnina && mkdir -p ~/bin && mv gnina ~/bin/

     The release is static in everything except CUDA, so it still needs
     a few NVIDIA shared libraries. With no root and no toolkit:

       python scripts/install_gnina_libs.py

     which unpacks them from NVIDIA's PyPI wheels into ~/cuda-libs.
     Check what was found with:

       python -m promoter_tf_drugs.wsldock

     On Windows there is no gnina build, so it runs through WSL and this
     package handles the path translation. A GPU changes the speed, not
     the answer; without one the rescoring falls back to CPU.

     Measured on PDB 1DB1: calcitriol redocked from a SMILES-embedded
     conformer, with no crystal coordinate in the input, came back at
     1.04 A from the crystal pose - correct by the usual 2 A convention.
     Affinity and CNNaffinity chose the same pose, and the CNN separated
     it far more sharply, collapsing to 0.63 for a 9.7 A pose the
     empirical score still rated -10.0 kcal/mol.

   smina - the fallback
     A single static file, and the easier install. Empirical scoring
     only, so no CNN rescoring.

     Linux / WSL / macOS
       curl -L -o smina https://sourceforge.net/projects/smina/files/smina.static/download
       chmod +x smina && sudo mv smina /usr/local/bin/

     Windows - no WSL needed
       conda install -c conda-forge smina

       Verified working: smina 2020.12.10 on win-64, docking PDB 1DB1
       with its own vitamin D ligands at AUC 1.0 against twelve
       unrelated approved drugs.

     AutoDock Vina also works and needs PDBQT for both files, which
     `pip install meeko` and `pip install openbabel-wheel` provide.
     There is no reason to prefer that route.

   Speed
     A screen is thousands of independent jobs, so --cpu is the biggest
     lever. Measured: 3.7 s per compound on 8 cores with smina, so a
     600-compound library is about ten minutes on 40 cores and several
     hours on one. gnina on a GPU is comparable per compound.

   None of this changes what a screen means. A better scoring function
   ranks a VALIDATED pocket better; it does not make an unvalidated one
   meaningful, and screen_by_docking still refuses below its AUC and
   enrichment floors. A network confident a molecule binds also says
   nothing about whether binding activates or inhibits, which is what
   decides the direction your gene moves.

2. RDKit, for turning SMILES into 3D ligands.

     pip install rdkit

Check it with:  python -m promoter_tf_drugs.wsldock
           or:  smina --version
"""


@dataclass
class PocketChoice:
    pdb_id: str
    path: str
    ligand_path: str | None
    ligand_resname: str | None
    box: tuple | None
    state: str                      # from classify_conformation
    note: str = ""

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def choose_pocket(client, uniprot: str, out_dir, *, limit: int = 12,
                  want_state: str | None = None) -> PocketChoice | None:
    """Pick a structure whose co-crystallised ligand defines the pocket.

    A holo structure is preferred over an apo one for a simple reason: the
    ligand tells you where the site is and roughly how big it is. Without
    one the box has to be guessed, and a guessed box is the commonest way
    a docking campaign quietly measures nothing.

    `want_state` filters on the agonist/antagonist conformation that
    classify_conformation reads out of the file, which matters because a
    compound that fits the active conformation is the one that might raise
    your gene.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ids = find_structures(client, uniprot, limit=limit)
    if not ids:
        log.info("no PDB structures for %s", uniprot)
        return None

    best = None
    best_atoms = (-1, -1)
    filtered_out: list = []
    for pdb_id in ids:
        path = fetch_pdb(client, pdb_id, out_dir)
        if path is None:
            continue
        text = path.read_text(errors="replace", encoding="utf-8")
        lig = best_ligand(text)
        if lig is None:
            continue
        conf = classify_conformation(text, pdb_id)
        if want_state and conf.state != want_state:
            # Dropped for its conformation, not for lacking a ligand. The
            # caller needs to tell those apart or it reports "the PDB has
            # no usable structure", which is false.
            filtered_out.append(f"{pdb_id} is {conf.state}")
            continue
        lig_path = write_ligand_pdb(
            text, lig, out_dir / f"{pdb_id}_{lig['resname']}.pdb")
        rec_path = prepare_receptor(
            text, out_dir / f"{pdb_id}_receptor.pdb",
            keep_chain=lig["chain"])
        try:
            box = box_from_ligand(lig_path.read_text(encoding="utf-8", errors="replace"))
        except ValueError:
            continue
        choice = PocketChoice(
            pdb_id=pdb_id, path=str(rec_path), ligand_path=str(lig_path),
            ligand_resname=lig["resname"], box=box, state=conf.state,
            note=(f"box placed on {lig['resname']} ({lig['n_atoms']} atoms) "
                  f"in chain {lig['chain']}; conformation {conf.state} "
                  f"({conf.evidence})"
                  + (f". {lig['resname']} is a crystallisation lipid or "
                     f"detergent, not necessarily a ligand: the box may "
                     f"be on the membrane-facing surface rather than in "
                     f"a pocket. Treat a screen against it as unvalidated "
                     f"unless this target genuinely binds lipids."
                     if lig.get("additive") else "")))
        # A structure whose conformation is known beats one that is not.
        if conf.state != "unknown":
            return choice
        # Among the unreadable ones, keep the BEST rather than the first.
        # For a GPCR every structure is unreadable, so "first" meant
        # whichever PDB id happened to sort earliest: CCR5 picked a
        # peptide fragment over 4MBS, which has maraviroc bound. A bigger
        # co-crystallised ligand is the better evidence of a drug-sized
        # pocket, which is the only thing being judged here.
        # A real ligand outranks any additive, regardless of size.
        rank = (0 if lig.get("additive") else 1, lig["n_atoms"])
        if best is None or rank > best_atoms:
            best, best_atoms = choice, rank

    if best is None and filtered_out:
        log.info("%s: %d structure(s) with a ligand were excluded for not "
                 "being in the %r conformation (%s)", uniprot,
                 len(filtered_out), want_state, "; ".join(filtered_out[:6]))
    return best


def direction_labels(client, symbol: str, *, library: list | None = None,
                     target_chembl_id: str | None = None) -> tuple:
    """Ligands of this receptor whose direction is known, and where from.

    A local drug table first, because it needs no network and most runs
    have no network. ChEMBL's mechanism endpoint second.
    """
    from . import agonism

    rows = library
    if rows is None:
        try:
            from . import localdrugs as _ld
            rows = _ld.library(None, approved_only=False)
        except Exception:                                   # noqa: BLE001
            rows = []
    local = agonism.labelled_ligands_local(symbol, rows or [])
    if local:
        return local, "a local drug table's mechanism column"

    if client is not None and target_chembl_id:
        try:
            got = agonism.fetch_labelled_ligands(client, target_chembl_id)
        except Exception as exc:                            # noqa: BLE001
            log.debug("labelled ligand fetch failed: %s", str(exc)[:90])
            got = []
        if got:
            return got, "ChEMBL mechanism annotations"
    return [], "no source"


def annotate_directions(run, symbol: str, labelled: list, *,
                        source: str = "?") -> dict:
    """Predict agonist or antagonist for each hit, in place.

    Returns a summary dict for the report. Every hit keeps its position:
    this annotates, it does not rank or filter. A hit the model will not
    call is left uncalled and counted, rather than defaulted to either
    direction.
    """
    from . import agonism

    hits = getattr(run, "hits", None)
    if hits is None:
        hits = (run or {}).get("hits") or []
    out = {"attempted": True, "source": source,
           "n_labelled": len(labelled), "n_called": 0,
           "n_out_of_domain": 0, "model_auc": None}

    if len(labelled) < 2 * agonism.MIN_PER_CLASS:
        out.update(attempted=False, reason=(
            f"only {len(labelled)} ligand(s) of known direction for "
            f"{symbol}, below what a per-receptor classifier needs. "
            f"Direction is not annotated rather than guessed."))
        return out

    try:
        model = agonism.train_model(labelled)
    except ImportError:
        out.update(attempted=False, reason=(
            "scikit-learn is not installed, so direction cannot be "
            "predicted. pip install 'promoter-tf-drugs[chemistry]'"))
        return out
    except Exception as exc:                                # noqa: BLE001
        out.update(attempted=False,
                   reason=f"direction model failed to train: "
                          f"{str(exc)[:120]}")
        return out

    if model is None:
        n_ag = sum(1 for x in labelled if x.get("label") == "agonist")
        out.update(attempted=False, reason=(
            f"{symbol} has {n_ag} agonist(s) and {len(labelled) - n_ag} "
            f"antagonist(s) of known direction, and a classifier needs at "
            f"least {agonism.MIN_PER_CLASS} of each. Direction is not "
            f"annotated rather than guessed."))
        return out

    out["model_auc"] = model.auc
    for h in hits:
        smi = h.get("smiles")
        if not smi:
            continue
        pred = agonism.predict(model, smi, drug_name=h.get("drug") or "?",
                               target_symbol=symbol)
        h["direction"] = pred.predicted
        h["direction_confidence"] = pred.confidence
        h["direction_in_domain"] = pred.in_domain
        if pred.probability_agonist is not None:
            h["probability_agonist"] = round(pred.probability_agonist, 3)
        if pred.nearest_tanimoto is not None:
            h["nearest_labelled_tanimoto"] = round(pred.nearest_tanimoto, 3)
        if pred.predicted in ("agonist", "antagonist"):
            out["n_called"] += 1
        if not pred.in_domain:
            out["n_out_of_domain"] += 1
    return out


def dock_screen_for_target(client, symbol: str, uniprot: str, library: list,
                           *, actives: list, decoys: list | None = None,
                           work_dir="docking_work", top_n: int = 25,
                           exhaustiveness: int = 8, seed: int = 0,
                           want_state: str | None = None,
                           state_note: str | None = None,
                           score_fn=None, cpu: int = 1,
                           jobs: int = 1,
                           strip_counter_ions: bool = False) -> dict:
    """Gene symbol to a validated, ranked shortlist of approved drugs.

    Returns the DockingScreen as a dict, plus which structure was used.
    Every refusal along the way is reported rather than raised: no
    structure, no co-crystallised ligand, no binary, or a pocket whose
    docking cannot recover its own known ligands.
    """
    import random

    if score_fn is None and _docking_binary() is None:
        return {"target": symbol, "error": "no docking binary on PATH.\n"
                                           + SETUP, "usable": False}

    pocket = choose_pocket(client, uniprot, work_dir, want_state=want_state)
    if pocket is None and want_state:
        # Was it the conformation filter, or is there genuinely no usable
        # structure? Asking again without the filter answers that, and
        # costs nothing: every structure is already cached on disk from
        # the first pass.
        unfiltered = choose_pocket(client, uniprot, work_dir)
        if unfiltered is not None:
            if unfiltered.state == "unknown":
                # Not "the other conformation" - no conformation at all.
                # The classifier reads coactivator and corepressor
                # peptides, which only nuclear receptors carry, so for a
                # GPCR or a kinase it returns unknown for every
                # structure. Reporting that as "none in the inactive
                # conformation" would turn a measurement that was never
                # made into a negative result.
                return {
                    "target": symbol, "usable": False,
                    "excluded_by_conformation": True,
                    "conformation_unreadable": True,
                    "error": (
                        f"{symbol} ({uniprot}) has a structure with a "
                        f"co-crystallised ligand ({unfiltered.pdb_id}), "
                        f"but its conformation could not be read, and "
                        f"the {want_state} conformation was asked for. "
                        f"This is not a finding about the structure: the "
                        f"classifier works by spotting a coactivator or "
                        f"corepressor peptide, which only nuclear "
                        f"receptors carry, so every structure of a GPCR, "
                        f"a kinase or an enzyme comes back unknown. "
                        f"Nothing was measured and nothing is claimed. "
                        f"Dock it anyway with --state any, which reports "
                        f"the conformation as unverified and withholds "
                        f"the direction claim that depends on it."),
                }
            return {
                "target": symbol, "usable": False,
                "excluded_by_conformation": True,
                "error": (
                    f"{symbol} ({uniprot}) has a structure with a "
                    f"co-crystallised ligand ({unfiltered.pdb_id}, "
                    f"conformation {unfiltered.state}), but none in the "
                    f"{want_state} conformation that was asked for. That "
                    f"IS a finding: the structure was classified and it "
                    f"is the other one. Docking it would answer a "
                    f"different question, so it is refused rather than "
                    f"substituted. Use --state any to dock it regardless "
                    f"and have the mismatch recorded."),
            }
    if pocket is None:
        return {
            "target": symbol, "usable": False,
            "error": (f"no usable structure for {symbol} ({uniprot}). Either "
                      f"the PDB has none, or none carries a "
                      f"co-crystallised ligand big enough to place a box on. "
                      f"Docking into a guessed box measures nothing, so this "
                      f"is refused rather than approximated."),
        }

    # Sample the gate when the target has more known ligands than it
    # needs. At random under the seed, never by potency: taking the
    # tightest binders would make the actives easier to recover and
    # raise the AUC, which is tuning the gate to pass rather than
    # measuring whether the pocket works.
    n_actives_available = len(actives)
    if len(actives) > MAX_VALIDATION_ACTIVES:
        actives = random.Random(seed).sample(actives,
                                             MAX_VALIDATION_ACTIVES)
        log.info("validation gate sampled %d of %d known ligand(s)",
                 len(actives), n_actives_available)

    if decoys is None:
        rng = random.Random(seed)
        pool = [d for d in library
                if d["smiles"] not in {a["smiles"] for a in actives}]
        decoys = rng.sample(pool, min(len(pool), max(50, len(actives) * 10)))

    run = screen_by_docking(
        Path(pocket.path), library, target=symbol, actives=actives,
        decoys=decoys, top_n=top_n, exhaustiveness=exhaustiveness,
        cpu=cpu, jobs=jobs,
        seed=seed, score_fn=score_fn,
        box=pocket.box, strip_counter_ions=strip_counter_ions)
    # Which way each hit would push the receptor. An annotation, never a
    # filter: the screen ranks on fit, and direction is a prediction
    # about the protein rather than about the pocket.
    labelled, src = direction_labels(client, symbol)
    run.direction = annotate_directions(run, symbol, labelled, source=src)
    if not run.direction.get("attempted") and run.direction.get("reason"):
        run.caveats.append(run.direction["reason"])

    out = run.to_dict()
    out["pocket"] = pocket.to_dict()
    # Whether the pocket's conformation was READ, not merely accepted.
    # A screen run with --state any against an unreadable structure is
    # still a screen; it just cannot support a claim about direction,
    # and the difference has to survive into the saved result.
    out["n_actives_available"] = n_actives_available
    out["n_actives_used_for_gate"] = len(actives)
    out["conformation_verified"] = pocket.state in ("active", "inactive")
    if not out["conformation_verified"]:
        out["conformation_caveat"] = (
            "The conformation of this structure could not be read, so "
            "nothing here says whether a hit would activate or inhibit. "
            "Binding and direction are different questions, and only the "
            "first was addressed.")

    # What this run was ASKED for, as against what it found.
    #
    # Without it a saved result cannot be audited after the fact, and
    # the case that needs auditing is ordinary: per-factor files are
    # named <FACTOR>.json, so a second workup overwrites or sits beside
    # the first, and a merge globs whatever is in the directory. The
    # requested conformation is derived from the direction and the
    # factor's role - lowering a gene through an activator wants that
    # activator's INACTIVE pocket - so the same factor docked for `up`
    # and for `down` gives different answers under the same filename.
    #
    # `matched` is the question a reader actually has: was the pocket
    # the one asked for. None means nothing specific was requested.
    out["asked_for"] = {
        "state": want_state,
        # Why it is what it is. `state` alone comes back None both for
        # `--state any`, which is an explicit opt-out, and for a factor
        # whose role could not be resolved, where nothing could be
        # derived. Those are different and the note is what separates
        # them. None when the caller supplied no note, rather than a
        # sentence invented on its behalf.
        "state_note": state_note,
        "pocket_state": pocket.state,
        "matched": None if not want_state else pocket.state == want_state,
    }
    try:
        from . import provenance as _prov
        out["produced"] = _prov.run_provenance()
    except Exception as exc:                                # noqa: BLE001
        log.warning("the provenance stamp could not be written to the "
                    "bundle (%s: %s); the bundle has no 'produced' record",
                    type(exc).__name__, str(exc)[:100])
        # Provenance is a record, not a result. A run that produced a
        # validated screen must not be thrown away because the stamp
        # could not be written.
        pass
    return out


# ===================================================================== #
# Prepare here, dock there: a portable bundle.                          #
# ===================================================================== #
# Everything above needs the network - RCSB for the structure, ChEMBL for
# the receptor's known ligands and the approved library. The docking
# itself needs a Linux binary. Those two requirements often do not live on
# the same machine: on Windows the network is here and smina is not, and
# in WSL smina is there and the network may not be.
#
# So the two halves are separable. `write_bundle` does every networked
# step and leaves a folder of plain files. `dock_from_bundle` reads that
# folder and runs with no network at all. Copy the folder, or just point
# WSL at it - Windows drives appear under /mnt/c.

BUNDLE_FILE = "bundle.json"


def write_bundle(client, symbol: str, uniprot: str, out_dir, *,
                 library_limit: int = 600, want_state: str | None = None,
                 state_note: str | None = None,
                 pdb_file=None, library_file=None, actives_file=None
                 ) -> dict:
    """Fetch everything docking needs into one folder. Needs the network.

    Returns the manifest. After this the folder is self-contained and the
    dock can run somewhere else, offline.
    """
    import json as _json

    from . import agonism
    from .drugs import _chembl_targets
    from .screening import approved_library

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- the structure
    if pdb_file:
        text = Path(pdb_file).read_text(errors="replace", encoding="utf-8")
        lig = best_ligand(text)
        if lig is None:
            return {"error": (f"{pdb_file} has no co-crystallised ligand big "
                              f"enough to place a box on. Supply a holo "
                              f"structure, or let --prepare search RCSB.")}
        lig_path = write_ligand_pdb(text, lig, out_dir / "ligand.pdb")
        rec_path = prepare_receptor(text, out_dir / "receptor.pdb",
                                    keep_chain=lig["chain"])
        pocket = PocketChoice(
            pdb_id=Path(pdb_file).stem, path=str(rec_path),
            ligand_path=str(lig_path), ligand_resname=lig["resname"],
            box=box_from_ligand(lig_path.read_text(encoding="utf-8", errors="replace")),
            state=classify_conformation(text, Path(pdb_file).stem).state,
            note=f"supplied structure {pdb_file}")
    else:
        pocket = choose_pocket(client, uniprot, out_dir,
                               want_state=want_state)
        if pocket is None:
            return {"error": (f"no usable structure for {symbol} "
                              f"({uniprot}): none carries a "
                              f"co-crystallised ligand to place a box on.")}
        # Normalise the paths into the bundle so it stays portable.
        for src, name in ((pocket.path, "receptor.pdb"),
                          (pocket.ligand_path, "ligand.pdb")):
            dst = out_dir / name
            if Path(src) != dst:
                dst.write_text(Path(src).read_text(errors="replace", encoding="utf-8"), encoding="utf-8")
        pocket.path = str(out_dir / "receptor.pdb")
        pocket.ligand_path = str(out_dir / "ligand.pdb")

    # --- the validation set and the library.
    # A local table is tried first and the network only if there is none,
    # because on a blocked network the ChEMBL calls are the step that
    # stops the whole analysis.
    from . import localdrugs as ld

    local_rows = ld.library(library_file)
    actives: list = []
    library: list = []
    if local_rows:
        actives = ld.actives_for(symbol, local_rows)
        library = [{"name": d["name"], "smiles": d["smiles"]}
                   for d in local_rows[:library_limit]]
        log.info("using a local drug table: %d library compound(s), "
                 "%d known ligand(s) for %s", len(library), len(actives),
                 symbol)
    if actives_file:
        actives = [{"chembl_id": None, "name": a["name"],
                    "smiles": a["smiles"], "pchembl": None}
                   for a in ld.read_table(actives_file, approved_only=False)]
    if not library:
        library = [{"name": d.get("name"), "smiles": d["smiles"]}
                   for d in approved_library(client)[:library_limit]
                   if d.get("smiles")]
    if not actives:
        tids = _chembl_targets(client, symbol, uniprot)
        for tid in tids[:2]:
            actives.extend(agonism.fetch_actives(client, tid))
        actives = [{"chembl_id": a.get("chembl_id"), "smiles": a["smiles"],
                    "pchembl": a.get("pchembl")}
                   for a in actives if a.get("smiles")]

    manifest = {
        "symbol": symbol, "uniprot": uniprot,
        "pdb_id": pocket.pdb_id, "ligand_resname": pocket.ligand_resname,
        "state": pocket.state, "note": pocket.note,
        # What was ASKED for, beside what was found. A bundle is
        # prepared in one place and docked in another, later, so by the
        # time it is used the command that chose the conformation is
        # gone - which makes this worth more here than in a run that
        # does both halves in one process.
        "asked_for": {
            "state": want_state,
            "state_note": state_note,
            "pocket_state": pocket.state,
            "matched": None if not want_state else pocket.state == want_state,
        },
        "box": list(pocket.box),
        "receptor": "receptor.pdb", "ligand": "ligand.pdb",
        "n_actives": len(actives), "n_library": len(library),
        "actives": actives, "library": library,
    }
    (out_dir / BUNDLE_FILE).write_text(
        _json.dumps(manifest, indent=2), encoding="utf-8")

    manifest["ready"] = len(actives) >= MIN_VALIDATION_ACTIVES
    if not manifest["ready"]:
        manifest["warning"] = (
            f"only {len(actives)} known ligand(s) for {symbol}, below the "
            f"{MIN_VALIDATION_ACTIVES} needed to check that docking into "
            f"this pocket works. The bundle is written, but the dock will "
            f"refuse to report.")
    return manifest


def load_bundle(bundle_dir) -> dict:
    """Read a bundle written by `write_bundle`. No network."""
    import json as _json

    d = Path(bundle_dir)
    f = d / BUNDLE_FILE
    if not f.exists():
        raise FileNotFoundError(
            f"{f} not found. Run the prepare step on a machine with network "
            f"access first.")
    m = _json.loads(f.read_text(encoding="utf-8"))
    m["_dir"] = str(d)
    m["_receptor"] = str(d / m["receptor"])
    m["_box"] = tuple(m["box"])
    return m


def dock_from_bundle(bundle_dir, *, top_n: int = 25,
                     exhaustiveness: int = 8, seed: int = 0,
                     score_fn=None, cpu: int = 1,
                     strip_counter_ions: bool = False) -> dict:
    """Run the validated screen from a bundle, with no network at all."""
    import random

    m = load_bundle(bundle_dir)
    actives = m["actives"]
    library = m["library"]

    rng = random.Random(seed)
    active_smiles = {a["smiles"] for a in actives}
    pool = [d for d in library if d["smiles"] not in active_smiles]
    decoys = rng.sample(pool, min(len(pool), max(50, len(actives) * 10)))

    run = screen_by_docking(
        Path(m["_receptor"]), library, target=m.get("symbol", "?"),
        actives=actives, decoys=decoys, top_n=top_n,
        exhaustiveness=exhaustiveness, seed=seed, score_fn=score_fn,
        box=m["_box"], cpu=cpu, strip_counter_ions=strip_counter_ions)
    out = run.to_dict()
    out["pocket"] = {k: m.get(k) for k in
                     ("pdb_id", "ligand_resname", "state", "note")}
    return out


# ===================================================================== #
# Entry point.  python -m promoter_tf_drugs.docking --check             #
# ===================================================================== #
def check_setup() -> dict:
    """What is installed, what is missing, and what to do about it."""
    import shutil as _sh

    out = {"ok": True, "found": {}, "missing": []}

    binary = _docking_binary()
    out["found"]["docking_binary"] = binary
    out["found"]["engine"] = engine_name(binary) if binary else None
    if binary is None:
        out["ok"] = False
        out["missing"].append(
            "a docking binary (smina preferred) - nothing can be docked")
    elif engine_name(binary) == "vina":
        ok_prep, why = vina_prep_available()
        out["found"]["vina_prep"] = ok_prep
        if not ok_prep:
            out["ok"] = False
            out["missing"].append(str(why))
        elif why:
            out["found"]["vina_prep_note"] = why

    try:
        import rdkit                                        # noqa: F401
        out["found"]["rdkit"] = True
    except ImportError:
        out["ok"] = False
        out["found"]["rdkit"] = False
        out["missing"].append("rdkit - SMILES cannot be turned into 3D "
                              "ligands (pip install rdkit)")

    out["found"]["curl_or_network"] = bool(_sh.which("curl"))
    return out


def _print_check() -> int:
    import json as _json
    c = check_setup()
    print("Docking setup")
    print("-------------")
    print(f"  docking binary : {c['found'].get('docking_binary') or 'NOT FOUND'}"
          + (f"  ({c['found'].get('engine')})"
             if c["found"].get("engine") else ""))
    print(f"  rdkit          : "
          f"{'yes' if c['found'].get('rdkit') else 'NOT INSTALLED'}")
    print()
    if c["ok"]:
        print("Ready. Everything docking needs is present.")
        return 0
    print("Not ready. Missing:")
    for m in c["missing"]:
        print(f"  - {m}")
    print()
    print(SETUP)
    return 1


def load_results(path) -> dict:
    """A docking result written on another machine, ready to fold back in.

    This is the third step of the split: prepare where the network is,
    dock where the binary is, and bring the JSON back here.
    """
    import json as _json

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{p} not found")
    data = _json.loads(p.read_text(encoding="utf-8"))
    runs = data if isinstance(data, list) else [data]
    for r in runs:
        if not isinstance(r, dict) or "usable" not in r:
            raise ValueError(
                f"{p} does not look like a docking result. Write one with "
                f"`--bundle DIR --out results.json`.")
    return data


def _direction_cell(hit: dict) -> str:
    """The predicted direction for one hit, or why there is none.

    A molecule outside the model's applicability domain reads as "-",
    not as a weak call: the model has not seen anything like it and a
    low-confidence guess would look like a measurement.
    """
    d = hit.get("direction")
    if not d or d == "undetermined":
        return "-"
    conf = hit.get("direction_confidence") or ""
    short = {"agonist": "agonist", "antagonist": "antagonist"}.get(d, d)
    return short if conf in ("high", "moderate") else f"{short}?"


def _report(res: dict, top: int = 25) -> None:
    """Print a screen result, from either the one-shot or bundle path."""
    if res.get("error"):
        print(f"Not reported: {res['error']}")
        return
    pk = res.get("pocket") or {}
    print(f"Structure : {pk.get('pdb_id')} "
          f"(ligand {pk.get('ligand_resname')})")

    state = pk.get("state")
    if res.get("conformation_verified"):
        print(f"Conformation: {state}, read from the structure")
    else:
        print(f"Conformation: {state} - NOT read from the structure, so "
              f"the docking addresses fit and not direction. A direction "
              f"column below, if there is one, is a separate ligand-based "
              f"prediction rather than a reading of this structure.")

    print(f"Validation: AUC {res.get('auc')}, "
          f"top-1% enrichment {res.get('ef1')}")
    gate = res.get("gate") or {}
    if gate:
        # What the two numbers rest on. The enrichment is counted over a
        # handful of compounds, so one more or fewer shifts it visibly.
        basis = gate.get("ef1_basis") or {}
        print(f"            {gate.get('n_actives_scored')} of "
              f"{gate.get('n_actives')} known ligand(s) and "
              f"{gate.get('n_decoys_scored')} of {gate.get('n_decoys')} "
              f"decoy(s) scored"
              + (f"; enrichment counts {basis.get('actives_in_top_n')} "
                 f"known ligand(s) among the top {basis.get('top_n')}"
                 if basis else ""))
    prep = res.get("ligand_prep") or {}
    if prep.get("strip_counter_ions"):
        print(f"Ligands   : counter-ions removed from {prep.get('n_changed')} "
              f"of {prep.get('n_library')} library entries"
              + (f"; {prep.get('n_ambiguous')} with a second organic "
                 f"component were kept whole"
                 if prep.get("n_ambiguous") else ""))
    n_av = res.get("n_actives_available")
    n_used = res.get("n_actives_used_for_gate")
    if n_av and n_used and n_av > n_used:
        print(f"            measured on {n_used} of {n_av} known "
              f"ligand(s), sampled at random")

    by_metric = res.get("auc_by_metric") or {}
    if len(by_metric) > 1:
        print()
        print("What each score was worth on this pocket's own ligands:")
        print(f"  {'metric':14s} {'AUC':>6s} {'EF1%':>7s}")
        for name, m in sorted(by_metric.items(),
                              key=lambda kv: -kv[1]["auc"]):
            ef = m.get("ef1")
            print(f"  {name:14s} {m['auc']:6.3f} "
                  f"{(f'{ef:.1f}' if ef is not None else '-'):>7s}")
        print("  Ranking below uses the affinity.")

    # Only show the CNN columns when the engine produced them, so a
    # smina run does not print blanks that look like measurements.
    hits = res.get("hits", [])[:top]
    has_cnn = any("cnn_affinity" in h for h in hits)

    d = res.get("direction") or {}
    if d.get("attempted"):
        auc = d.get("model_auc")
        print()
        print(f"Direction predicted from {d.get('n_labelled')} ligand(s) of "
              f"known direction ({d.get('source')})"
              + (f", cross-validated AUC {auc:.2f}" if auc else "")
              + f". Called for {d.get('n_called')} of {len(hits)} hit(s); "
              + str(d.get("n_out_of_domain"))
              + (" was" if d.get("n_out_of_domain") == 1 else " were")
              + " unlike anything the model was trained on and is left "
                "uncalled.")
    print()
    if has_cnn:
        print(f"{'drug':28s} {'kcal/mol':>9s} {'CNNscore':>9s} "
              f"{'CNNaff':>7s} {'direction':>12s}")
        for h in hits:
            cs = h.get("cnn_score")
            ca = h.get("cnn_affinity")
            print(f"{h['drug'][:28]:28s} {h['score']:9.2f} "
                  f"{(f'{cs:.3f}' if cs is not None else '-'):>9s} "
                  f"{(f'{ca:.2f}' if ca is not None else '-'):>7s} "
                  f"{_direction_cell(h):>12s}")
    else:
        print(f"{'drug':30s} {'kcal/mol':>9s} {'direction':>12s}")
        for h in hits:
            print(f"{h['drug'][:30]:30s} {h['score']:9.2f} "
                  f"{_direction_cell(h):>12s}")

    print()
    if res.get("conformation_caveat"):
        print(f"  - {res['conformation_caveat']}")
    for c in res.get("caveats", []):
        print(f"  - {c}")


def _write_refusal(out_path, target: str, uniprot: str | None,
                   error: str) -> None:
    """Record a refusal that happened before any screen could run.

    A refusal is a result: it says the scores for this target would
    carry no information, and that is different from the target never
    having been tried. Without this the two are indistinguishable in a
    directory of outputs.
    """
    if not out_path:
        return
    import json as _json

    payload = {
        "target": target, "uniprot": uniprot, "usable": False,
        "n_docked": 0, "n_failed": 0, "hits": [],
        "error": error, "refused_before_screening": True,
    }
    try:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        Path(out_path).write_text(_json.dumps(payload, indent=2),
                                  encoding="utf-8")
        print(f"Recorded the refusal in {out_path}")
    except OSError as exc:                                  # noqa: BLE001
        log.warning("could not record the refusal: %s", exc)


def main(argv=None) -> int:
    import argparse
    import json as _json

    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.docking",
        description=(
            "Dock approved drugs into a receptor pocket and report the ones "
            "that fit, after checking that docking into that pocket can "
            "recover its own known ligands."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=SETUP)
    p.add_argument("--check", action="store_true",
                   help="report what is installed and what is missing, then "
                        "exit")
    p.add_argument("--prepare", metavar="DIR",
                   help="fetch the structure, the receptor's known ligands "
                        "and the approved library into DIR, then stop. This "
                        "is the only step that needs the network, so run it "
                        "where you have one. No docking binary required")
    p.add_argument("--bundle", metavar="DIR",
                   help="dock from a folder written by --prepare, with no "
                        "network at all. This is the step that needs smina, "
                        "so run it wherever smina is - WSL sees your "
                        "Windows drives under /mnt/c")
    p.add_argument("--library", metavar="FILE",
                   help="local drug table with SMILES instead of fetching "
                        "the library from ChEMBL. The Broad Repurposing Hub "
                        "export works and also supplies the known ligands "
                        "per target, so one file removes the whole ChEMBL "
                        "dependency. data/repurposing_drugs.txt is picked up "
                        "automatically")
    p.add_argument("--actives", metavar="FILE",
                   help="local known ligands for the validation gate, as a "
                        "name/SMILES table. Defaults to the ones annotated "
                        "against this target inside --library")
    p.add_argument("--pdb", metavar="FILE",
                   help="use this structure instead of searching RCSB. It "
                        "must carry a co-crystallised ligand, which is what "
                        "places the box")
    p.add_argument("--gene", help="receptor gene symbol, e.g. VDR or RARA")
    p.add_argument("--uniprot",
                   help="UniProt accession, if the symbol does not resolve")
    p.add_argument("--state", choices=["active", "inactive", "any"],
                   help="dock into the agonist-bound (active) or "
                        "antagonist-bound (inactive) conformation. Leave it "
                        "out and it is derived from --direction and the "
                        "factor's role, which is what decides it: to raise "
                        "a gene you want an ACTIVATOR active and a "
                        "REPRESSOR inactive. 'any' docks whichever pocket "
                        "exists and records the conformation as unverified "
                        "- needed for GPCRs, kinases and enzymes, whose "
                        "conformation this cannot read, and it withholds "
                        "the direction claim rather than guessing it")
    p.add_argument("--direction", choices=["up", "down"], default="up",
                   help="which way you want YOUR GENE to move (default up). "
                        "Used with the factor's role to pick the "
                        "conformation when --state is not given")
    p.add_argument("--product", action="store_true",
                   help="the target is the gene's own protein, not a "
                        "transcription factor. The conformation then "
                        "follows --direction directly: raising the gene's "
                        "activity wants the agonist-bound structure, "
                        "lowering it wants the antagonist-bound one, with "
                        "no factor role in between")
    p.add_argument("--role", choices=["activator", "repressor", "dual"],
                   help="the factor's role, if it cannot be resolved from "
                        "local annotation. Without it the conformation is "
                        "not filtered rather than guessed")
    p.add_argument("--top", type=int, default=25,
                   help="how many hits to report (default 25)")
    p.add_argument("--jobs", type=int, default=1,
                   help="how many compounds to dock AT ONCE. A screen is "
                        "thousands of independent jobs and a docking "
                        "program spreads one job across cores with heavy "
                        "diminishing returns, so --jobs 8 --cpu 4 beats "
                        "--cpu 32 by about 2.5x. Scores do not change "
                        "(default 1)")
    p.add_argument("--cpu", type=int, default=1,
                   help="cores per docking job. A screen is thousands of "
                        "independent jobs, so raising this is the single "
                        "biggest speed-up available (default 1)")
    p.add_argument("--strip-counter-ions", action="store_true",
                   help="dock each drug without its salt (chloride, "
                        "calcium, maleate and so on) and with one copy "
                        "rather than two. About a third of the library is "
                        "stored as a salt. Off by default, which docks the "
                        "structure as stored; the result records how many "
                        "entries changed, and a drug sold with a second "
                        "organic component is kept whole and counted as "
                        "ambiguous")
    p.add_argument("--exhaustiveness", type=int, default=8,
                   help="docking search effort; 8 is the usual default, 16 "
                        "is slower and steadier (default 8)")
    p.add_argument("--library-limit", type=int, default=600,
                   help="how many approved drugs to dock. Each one is "
                        "seconds of CPU, so the full library is an "
                        "overnight job (default 600)")
    p.add_argument("--work-dir", default="docking_work",
                   help="where fetched structures are cached")
    p.add_argument("--out", help="write the result as JSON here")
    p.add_argument("--cache-dir", help="HTTP cache directory")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s")

    if args.check:
        return _print_check()

    # Offline half: dock from a prepared folder. Needs smina, not network.
    if args.bundle:
        chk = check_setup()
        if not chk["ok"] and args.out is None:
            print("Cannot dock yet:")
            for m in chk["missing"]:
                print(f"  - {m}")
            print()
            print(SETUP)
            return 1
        res = dock_from_bundle(args.bundle, top_n=args.top,
                               exhaustiveness=args.exhaustiveness,
                               cpu=args.cpu,
                               strip_counter_ions=args.strip_counter_ions)
        _report(res, args.top)
        if args.out:
            Path(args.out).write_text(_json.dumps(res, indent=2),
                                      encoding="utf-8")
            print(f"\nWrote {args.out}")
        return 0 if not res.get("error") else 1

    if not args.gene and not args.uniprot:
        p.error("--gene, --uniprot, --bundle or --check is required")

    # Networked half: fetch everything, write a portable folder, stop.
    # The conformation to dock into follows from the factor's role and the
    # direction wanted for the gene. Resolved once here so the two call
    # sites below cannot disagree.
    want_state = args.state
    state_note = "set explicitly with --state"
    if want_state == "any":
        # An explicit opt-out, not an absence of one. Left as None it
        # would fall into the derivation below and come back as exactly
        # the state the user asked to stop filtering on.
        want_state = None
        state_note = ("not filtered, by request: the conformation is "
                      "recorded as unverified and no direction is claimed")
    elif want_state is None and args.product:
        want_state = state_for_product(args.direction)
        state_note = (f"the gene's own product, to move its activity "
                      f"{args.direction}")
    elif want_state is None:
        role = args.role
        if role is None and args.gene:
            try:
                from . import localannot as _la
                _ann = _la.annotate_local([args.gene]).get(args.gene) or {}
                role = _ann.get("role")
            except Exception:                               # noqa: BLE001
                role = None
        want_state = state_for(role or "unknown", args.direction)
        if want_state:
            state_note = (f"derived: a {role} factor, to move the gene "
                          f"{args.direction}")
        else:
            state_note = (
                f"not filtered: the factor's role is "
                f"{role or 'unresolved'}, so there is no defensible "
                f"conformation to prefer. Pass --role or --state to "
                f"choose one")
    print(f"Conformation   : {want_state or 'any'} ({state_note})")

    if args.prepare:
        from .http import Client as _C
        from . import annotate as _ann
        client = _C(cache_dir=args.cache_dir) if args.cache_dir else _C()
        uniprot = args.uniprot
        symbol = args.gene or args.uniprot
        if not uniprot:
            info = (_ann.annotate_tfs(client, [args.gene]) or {}).get(
                args.gene.upper())
            uniprot = (info or {}).get("uniprot")   # a dict KEY, not an attribute
            if not uniprot:
                print(f"Could not resolve {args.gene} to a UniProt "
                      f"accession. Pass --uniprot.")
                return 1
        print(f"Preparing {symbol} ({uniprot}) into {args.prepare} ...")
        man = write_bundle(client, symbol, uniprot, args.prepare,
                           library_limit=args.library_limit,
                           want_state=want_state, state_note=state_note,
                           pdb_file=args.pdb,
                           library_file=args.library,
                           actives_file=args.actives)
        if man.get("error"):
            print(f"Refused: {man['error']}")
            return 1
        print(f"  structure      : {man['pdb_id']} "
              f"(ligand {man['ligand_resname']}, "
              f"conformation {man['state']})")
        print(f"  known ligands  : {man['n_actives']}")
        print(f"  library        : {man['n_library']}")
        if man.get("warning"):
            print(f"\n  WARNING: {man['warning']}")
        print(f"\nBundle written. It needs no network from here. Dock it "
              f"with:\n"
              f"  python -m promoter_tf_drugs.docking --bundle "
              f"{args.prepare}")
        return 0

    chk = check_setup()
    if not chk["ok"]:
        print("Cannot dock yet:")
        for m in chk["missing"]:
            print(f"  - {m}")
        print()
        print(SETUP)
        return 1

    from .http import Client
    from . import annotate as ann_mod, agonism

    client = Client(cache_dir=args.cache_dir) if args.cache_dir else Client()

    uniprot = args.uniprot
    symbol = args.gene or args.uniprot
    if not uniprot:
        info = (ann_mod.annotate_tfs(client, [args.gene]) or {}).get(
            args.gene.upper()) or (
            ann_mod.annotate_tfs(client, [args.gene]) or {}).get(args.gene)
        uniprot = (info or {}).get("uniprot")   # a dict KEY, not an attribute
        if not uniprot:
            msg = (f"could not resolve {args.gene} to a UniProt "
                   f"accession, so there was nothing to dock against")
            print(msg + ". Pass --uniprot explicitly.")
            _write_refusal(args.out, symbol, None, msg)
            return 1

    print(f"Receptor  : {symbol} ({uniprot})")

    # The library first: it is also where this target's known ligands
    # come from, so resolving the actives before loading it read
    # `local_rows` before it was assigned and the whole command died.
    from . import localdrugs as ld
    local_rows = ld.library(args.library)
    if local_rows:
        library = [{"name": d["name"], "smiles": d["smiles"]}
                   for d in local_rows[:args.library_limit]]
        print(f"Library from a local table: {len(library)} compound(s)")
    else:
        from .screening import approved_library
        library = approved_library(client)[:args.library_limit]
    print(f"Approved library: {len(library)} compound(s)")

    # The validation set: this receptor's own measured ligands.
    actives: list = []
    if args.actives:
        actives = [{"smiles": a["smiles"], "pchembl": None}
                   for a in ld.read_table(args.actives, approved_only=False)]
    elif local_rows:
        actives = ld.actives_for(symbol, local_rows)
    if not actives:
        # The local potency tables, before the network. BindingDB holds
        # thousands of measured ligands for targets the drug table
        # annotates with none, and going online for them when they are
        # already on disk makes an offline run refuse a target it could
        # have validated.
        from . import localpotency as _lp
        try:
            rows = _lp.actives_for(uniprot)
        except Exception as exc:                            # noqa: BLE001
            log.warning("local potency lookup failed: %s", str(exc)[:120])
            rows = []
        if rows:
            actives = [{"smiles": r["smiles"], "pchembl": r.get("pchembl")}
                       for r in rows if r.get("smiles")]
            print(f"Known ligands from a local potency table: "
                  f"{len(actives)}")
    if not actives:
        from .drugs import _chembl_targets
        tids = _chembl_targets(client, symbol, uniprot)
        for tid in tids[:2]:
            actives.extend(agonism.fetch_actives(client, tid))
    print(f"Known ligands for validation: {len(actives)}")
    if len(actives) < MIN_VALIDATION_ACTIVES:
        msg = (f"{len(actives)} measured ligand(s) for {symbol}, below the "
               f"{MIN_VALIDATION_ACTIVES} needed to check whether docking "
               f"into this pocket works at all. An unvalidated dock "
               f"produces confident nonsense, so nothing is reported")
        print(f"\nRefused: {msg}.")
        _write_refusal(args.out, symbol, uniprot, msg)
        return 1

    # The gate is sampled, so the estimate has to use the capped
    # count or it reports several times the work to be done.
    n_gate = min(len(actives), MAX_VALIDATION_ACTIVES)
    n_jobs = (len(library) + n_gate
              + min(len(library), max(50, n_gate * 10)))
    print(f"\nThis will run roughly {n_jobs} docking jobs "
          f"({len(library)} library, {n_gate} known ligand(s) and "
          f"their decoys for the validation gate). At a few seconds "
          f"each that is {n_jobs * 4 // 60} minutes or more.\n")

    res = dock_screen_for_target(
        client, symbol, uniprot, library, actives=actives,
        work_dir=args.work_dir, top_n=args.top,
        exhaustiveness=args.exhaustiveness, want_state=want_state,
        state_note=state_note,
        cpu=args.cpu, jobs=args.jobs,
        strip_counter_ions=args.strip_counter_ions)

    _report(res, args.top)

    if args.out:
        Path(args.out).write_text(_json.dumps(res, indent=2),
                                  encoding="utf-8")
        print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
