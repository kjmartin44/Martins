"""Exploratory docking into a pocket that has no co-crystallised ligand.

Why this exists
---------------
The docking stage refuses a target whose pocket cannot first recover its own
known ligands. That is the right default: a score from an unchecked pocket is
a number with no meaning. But it leaves a gap. Some genes have a model in the
AlphaFold database, and a close relative that has been crystallised with a
ligand, and no experimental structure of the target with one. For those the
question "does anything approved fit here at all" is reasonable to ask, as
long as the answer is labelled for what it is.

What this module does
---------------------
1. Downloads the AlphaFold model of the target.
2. Finds the same pocket on that model by transplanting it from a template: a
   related protein whose structure holds a ligand. The template and the
   model are aligned on sequence, superposed on the residues that agree, and
   the template ligand's position in the model's frame gives the docking box.
3. Docks the library into that box.
4. Ranks by how well each compound scores for its size, not by raw score.
5. Checks whether known members of a drug class that binds this kind of
   pocket (approved kinase inhibitors) rise to the top.

What it does not do
-------------------
**It is not validated.** The normal docking gate asks the pocket to recover
known ligands of this protein. Here there are none, so there is no gate. A
result says a compound is shaped to fit a model of a pocket. It does not say
the compound binds the protein, in a cell, at a dose anyone could reach.

**The model is a prediction of an unbound protein.** Side chains are not in
the positions a ligand would induce, and a pocket can be partly closed. The
run reports how many template ligand atoms land inside protein atoms. A high
count means the pocket is narrower in the model than in the template.

**The template is a relative, not the target.** Transplanting a site assumes
the relative's ATP site (or other pocket) is conserved. The aligned motif is
reported so a reader can see where it is not. A variant motif is a reason to
ask whether the domain binds the ligand class at all.

**The class check is not specific to the target.** Approved kinase inhibitors
score well in almost any model of an ATP site. Their enrichment shows the
pocket looks like an ATP site to the scoring function. It does not show that
the target is inhibited by them.

**Nothing here goes into a ranking.** The output is separate from the
evidence ledger and the candidate table, and the main pipeline does not call
this module.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np

from . import docking

log = logging.getLogger(__name__)

BANNER = ("UNVALIDATED EXPLORATORY RESULT. A compound listed here is shaped to "
          "fit a predicted model of a pocket. That is not evidence that it "
          "binds the protein, changes its activity, or changes expression.")

MIN_POCKET_PLDDT = 70.0       # the pocket must be confidently modelled
MAX_CORE_RMSD = 3.0           # superposition quality, in angstroms
MIN_ALIGNED = 100             # residues aligned between model and template
MIN_IDENTITY = 0.25           # fraction of aligned residues that agree
BOX_PADDING = 6.0             # added to the template ligand's extent
CLASH_DISTANCE = 2.0          # a template ligand atom this close is a clash
POCKET_DISTANCE = 4.5
MIN_CLASS_MEMBERS = 5
# Residues below this pLDDT are left out of the receptor: very low confidence
# means predicted disorder, and a disordered tail looks like an open cavity.
MIN_RESIDUE_PLDDT = 50.0
# Compounds are ranked only when their size is within these multiples of the
# template ligand's heavy atoms. A pocket sized for ADP (27 heavy atoms) does
# not hold a 60-atom glycoside, and a 10-atom fragment scores on efficiency.
# A heuristic range, not a measured one.
SIZE_WINDOW = (0.5, 1.6)
# The class check passes when the control class ranks above this AUC after
# the size correction. Chance is 0.5. Below it, the pocket does not look like
# the pocket that class binds, and the ranking is not read as "these fit".
MIN_CONTROL_AUC = 0.6
DFG_LIKE = "D[FLWY]G"         # the kinase motif that starts the activation loop

_THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
          "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
          "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
          "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
_SIMILAR = ("ILVM", "FYW", "KRH", "DE", "NQ", "ST", "AG")


# ---------------------------------------------------------------- structure
def _xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]),
                     float(line[46:54])])


def ca_table(pdb_text: str, chain: str | None = None) -> list:
    """One entry per residue: number, one-letter code, CA position, B-factor."""
    out, seen = [], set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if chain is not None and line[21] != chain:
            continue
        key = (line[21], int(line[22:26]), line[26])
        if key in seen:
            continue
        seen.add(key)
        out.append({"num": int(line[22:26]),
                    "aa": _THREE.get(line[17:20], "X"),
                    "xyz": _xyz(line), "b": float(line[60:66])})
    return out


def heavy_atoms(pdb_text: str) -> list:
    """(residue number, position) for each non-hydrogen protein atom."""
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[76:78].strip() not in ("H", "D"):
            out.append((int(line[22:26]), _xyz(line)))
    return out


def ligand_atoms(pdb_text: str, resname: str, chain: str | None) -> np.ndarray:
    """Positions of one copy of a template ligand's heavy atoms.

    The first residue found on the chain. Hydrogens are skipped (a structure
    refined with them lists about a third more atoms, and they sit close to
    the protein by construction), and so are alternate positions after the
    first.
    """
    keep: list = []
    first = None
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM") or line[17:20].strip() != resname:
            continue
        if chain is not None and line[21] != chain:
            continue
        if line[76:78].strip() in ("H", "D") or line[16] not in (" ", "A"):
            continue
        ident = (line[21], line[22:27])
        first = first or ident
        if ident == first:
            keep.append(_xyz(line))
    return np.array(keep)


# ---------------------------------------------------------------- alignment
def _score(x: str, y: str) -> float:
    if x == y:
        return 4.0
    return 1.0 if any(x in g and y in g for g in _SIMILAR) else -1.0


def align(a: str, b: str, gap: float = 3.0) -> list:
    """Local alignment of two sequences: (index in a, index in b) pairs.

    Smith-Waterman with a linear gap penalty and a small similarity table.
    Two kinase domains at 30 percent identity align well enough for a
    superposition on this, and the result is checked afterwards by the
    superposition's own error, so a poor alignment is refused downstream.
    """
    m, n = len(a), len(b)
    h = np.zeros((m + 1, n + 1))
    move = np.zeros((m + 1, n + 1), dtype=np.int8)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = h[i - 1, j - 1] + _score(a[i - 1], b[j - 1])
            up, left = h[i - 1, j] - gap, h[i, j - 1] - gap
            best = max(0.0, diag, up, left)
            h[i, j] = best
            if best == 0.0:
                move[i, j] = 0
            elif best == diag:
                move[i, j] = 1
            elif best == up:
                move[i, j] = 2
            else:
                move[i, j] = 3
    i, j = (int(v) for v in np.unravel_index(int(h.argmax()), h.shape))
    pairs = []
    while h[i, j] > 0:
        if move[i, j] == 1:
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif move[i, j] == 2:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def _kabsch(p: np.ndarray, q: np.ndarray) -> tuple:
    """Rotation r and shift t that map q onto p: q @ r.T + t."""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((q - qc).T @ (p - pc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, pc - r @ qc


def superpose(p: np.ndarray, q: np.ndarray, rounds: int = 5) -> dict:
    """Fit q onto p using only the residues that agree, dropping the rest.

    Loops, side chains and domain hinges move between related proteins. A
    fit over every aligned residue is pulled off by them, so each round
    keeps the closer part and refits. `rmsd` is over the kept core, and
    `n_core` says how many that was.
    """
    keep = np.arange(len(p))
    r, t = _kabsch(p, q)
    for _ in range(rounds):
        r, t = _kabsch(p[keep], q[keep])
        dist = np.linalg.norm(q @ r.T + t - p, axis=1)
        cut = max(3.0, float(np.percentile(dist, 60)))
        nxt = np.where(dist < cut)[0]
        if len(nxt) < 10 or len(nxt) == len(keep):
            break
        keep = nxt
    r, t = _kabsch(p[keep], q[keep])
    dist = np.linalg.norm(q[keep] @ r.T + t - p[keep], axis=1)
    return {"r": r, "t": t, "n_core": int(len(keep)),
            "rmsd": float(np.sqrt((dist ** 2).mean()))}


# ------------------------------------------------------------------- site
def box_around(points: np.ndarray) -> list:
    """(cx, cy, cz, sx, sy, sz): the points' extent plus padding on each side."""
    centre = points.mean(0)
    size = (points.max(0) - points.min(0)) + 2 * BOX_PADDING
    return [round(float(v), 2) for v in (*centre, *size)]


def transplant_site(model_text: str, template_text: str, *,
                    template_chain: str, template_ligand: str) -> dict:
    """Where the template's ligand sits in the model's frame, and how good.

    Returns the docking `box`, the model residues near the transplanted
    ligand, their mean pLDDT, how many ligand atoms land inside protein
    atoms, and what the model has at the template's DFG motif. An `error`
    entry means the site could not be transplanted and says why.
    """
    model = ca_table(model_text)
    tmpl = ca_table(template_text, template_chain)
    lig = ligand_atoms(template_text, template_ligand, template_chain)
    if not model:
        return {"error": "no residues could be read from the model"}
    if not tmpl:
        return {"error": f"template has no residues on chain "
                         f"{template_chain}"}
    if not len(lig):
        return {"error": f"template has no {template_ligand} on chain "
                         f"{template_chain}"}

    ma = "".join(r["aa"] for r in model)
    ta = "".join(r["aa"] for r in tmpl)
    pairs = align(ma, ta)
    site: dict = {"n_aligned": len(pairs)}
    if len(pairs) < MIN_ALIGNED:
        site["error"] = (f"only {len(pairs)} residues align between the "
                         f"model and the template (need {MIN_ALIGNED}); the "
                         f"two are not close enough to transplant a site")
        return site
    ident = sum(ma[i] == ta[j] for i, j in pairs) / len(pairs)
    site["identity"] = round(ident, 3)
    if ident < MIN_IDENTITY:
        site["error"] = (f"aligned residues are {ident:.0%} identical (need "
                         f"{MIN_IDENTITY:.0%}); too distant to transplant a "
                         f"site")
        return site

    p = np.array([model[i]["xyz"] for i, _ in pairs])
    q = np.array([tmpl[j]["xyz"] for _, j in pairs])
    fit = superpose(p, q)
    site["core_rmsd"] = round(fit["rmsd"], 2)
    site["n_core"] = fit["n_core"]
    if fit["rmsd"] > MAX_CORE_RMSD:
        site["error"] = (f"the superposition core is off by "
                         f"{fit['rmsd']:.1f} A (limit {MAX_CORE_RMSD}); the "
                         f"fold does not match well enough to place a site")
        return site

    moved = lig @ fit["r"].T + fit["t"]
    site["box"] = box_around(moved)

    atoms = heavy_atoms(model_text)
    pos = np.array([a[1] for a in atoms])
    dist = np.linalg.norm(moved[:, None, :] - pos[None], axis=2)
    site["n_ligand_atoms"] = int(len(moved))
    site["n_clashing_atoms"] = int((dist.min(1) < CLASH_DISTANCE).sum())
    near = sorted({atoms[k][0] for k in np.where(dist.min(0)
                                                 < POCKET_DISTANCE)[0]})
    site["pocket_residues"] = near
    plddt = {r["num"]: r["b"] for r in model}
    scored = [plddt[n] for n in near if n in plddt]
    site["pocket_mean_plddt"] = (round(float(np.mean(scored)), 1)
                                 if scored else None)
    if not near:
        site["error"] = (f"no model atom lies within {POCKET_DISTANCE} A of "
                         f"the transplanted ligand, so the site is not in "
                         f"the protein; the alignment placed it wrongly")
        return site
    if not scored or site["pocket_mean_plddt"] < MIN_POCKET_PLDDT:
        site["error"] = (f"the transplanted pocket has mean pLDDT "
                         f"{site['pocket_mean_plddt']} (need "
                         f"{MIN_POCKET_PLDDT:.0f}); docking into a poorly "
                         f"modelled region scores artefacts")
        return site

    # What the model has where the template has its DFG-like motif (DFG in
    # most kinases, DLG in a few, BUB1 among them).
    to_model = {j: i for i, j in pairs}
    for m in re.finditer(DFG_LIKE, ta):
        trio = [to_model.get(m.start() + k) for k in range(3)]
        if None not in trio:
            site["template_dfg"] = m.group()
            site["aligned_dfg"] = "".join(ma[i] for i in trio)
            break
    return site


# ------------------------------------------------------------------ model
def fetch_model(client, uniprot: str, out_dir) -> tuple:
    """The current AlphaFold model of a protein: (path, url), or (None, why).

    The file name carries a model version, and old versions are removed from
    the server, so the address is read from the prediction record rather than
    built from a version number that will go stale.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        url = docking.alphafold_url(client, uniprot)
        path = out_dir / Path(url).name
        if not path.exists():
            path.write_text(client.request(url, parse_json=False),
                            encoding="utf-8")
        return path, url
    except Exception as exc:                                # noqa: BLE001
        return None, (f"could not fetch the AlphaFold model of {uniprot}: "
                      f"{str(exc)[:100]}")


# ----------------------------------------------------------------- ranking
def _heavy_atom_count(smiles: str) -> int | None:
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        return mol.GetNumHeavyAtoms() if mol is not None else None
    except ImportError:
        return None


def rank_by_residual(rows: list, window: tuple | None = None) -> dict:
    """Score each compound against what its size predicts.

    A docking score usually falls (improves) with heavy atom count whatever
    the pocket, so ranking by raw score hands the top of the list to the
    largest molecules. Fit score on heavy atoms and rank by the residual.
    `ligand_efficiency` (score per heavy atom) is added for reading, not for
    ranking.

    `window` is `(low, high)` heavy atoms. Only compounds inside it are
    fitted and ranked, because a line through every size is pulled around by
    molecules that cannot occupy the pocket at all: in a narrow pocket the
    slope turns positive, and the residual then rewards the largest molecule
    that still scored. The others keep their scores and `in_window` False.

    Each row needs `affinity` and `heavy_atoms`; rows missing either are
    left out. Returns the ranked rows (with `residual`, `rank` and
    `ligand_efficiency`), the unranked ones, and the fit used. The fit is
    None when there were too few compounds or too little spread in size to
    fit one, in which case the raw score is the residual.
    """
    scored = [r for r in rows if r.get("affinity") is not None
              and r.get("heavy_atoms")]
    for r in scored:
        r["in_window"] = (window is None
                          or window[0] <= r["heavy_atoms"] <= window[1])
    ok = [r for r in scored if r["in_window"]]
    outside = [r for r in scored if not r["in_window"]]
    fit = None
    if len(ok) >= 10 and len({r["heavy_atoms"] for r in ok}) > 2:
        slope, icpt = np.polyfit([r["heavy_atoms"] for r in ok],
                                 [r["affinity"] for r in ok], 1)
        fit = {"slope": round(float(slope), 4),
               "intercept": round(float(icpt), 3), "n": len(ok)}
        for r in ok:
            r["residual"] = round(float(
                r["affinity"] - (slope * r["heavy_atoms"] + icpt)), 3)
    else:
        for r in ok:
            r["residual"] = r["affinity"]
    for r in ok:
        r["ligand_efficiency"] = round(-r["affinity"] / r["heavy_atoms"], 3)
    ok.sort(key=lambda r: r["residual"])
    for i, r in enumerate(ok, 1):
        r["rank"] = i
    return {"ranked": ok, "unranked": outside, "fit": fit}


def _round(v):
    return None if v is None else round(v, 3)


def class_control(rows: list, class_name: str) -> dict:
    """Do approved members of a pocket-matched drug class rise to the top?

    Reported for both the raw score and the size-corrected residual. This
    tests whether the pocket looks like the pocket that class binds, to the
    scoring function. It cannot say anything about the target specifically.
    """
    from . import systemic
    labels, raw, res = [], [], []
    for r in rows:
        cls = systemic.classify(r["name"])
        labels.append(bool(cls and cls.name == class_name))
        raw.append(r.get("affinity"))
        res.append(r.get("residual"))
    n = sum(labels)
    out = {"class": class_name, "n_members": n, "n_compared": len(rows)}
    if n < MIN_CLASS_MEMBERS:
        out["note"] = (f"only {n} member(s) of the class in the library (need "
                       f"{MIN_CLASS_MEMBERS}), so no enrichment is reported")
        return out
    out["auc_raw_score"] = _round(docking.roc_auc(raw, labels))
    out["auc_size_corrected"] = _round(docking.roc_auc(res, labels))
    out["ef1_size_corrected"] = docking.enrichment_factor(res, labels)
    out["note"] = ("This shows whether the pocket looks like a kinase ATP "
                   "site to the scoring function. It does not show the "
                   "target is inhibited by these drugs.")
    return out


# --------------------------------------------------------------------- run
def run(gene: str, uniprot: str, library: list, *, template_pdb: str,
        template_chain: str, template_ligand: str, client=None,
        model_text: str | None = None, template_text: str | None = None,
        out_dir="exploratory_work", exhaustiveness: int = 16, seed: int = 0,
        jobs: int = 1, cpu: int = 1, strip_counter_ions: bool = False,
        score_fn=None, top_n: int = 25, reference: bool = False,
        control_class: str = "antineoplastic kinase inhibitor") -> dict:
    """Dock a library into a transplanted pocket on an AlphaFold model.

    `model_text` and `template_text` are for tests and offline use; without
    them the structures are downloaded through `client`. `score_fn` stands in
    for the docking binary the way it does in `docking.screen_by_docking`.
    `reference` also docks the same compounds into the template's own pocket,
    which costs as much again and shows whether a failed class check is the
    model's doing. Failures are returned in `error`, never raised.
    """
    res = {"gene": gene, "uniprot": uniprot, "status": "unvalidated",
           "banner": BANNER, "usable": False, "error": None, "caveats": [],
           "site": {}, "model": {}, "ranking": [], "class_control": {}}
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if score_fn is None and docking._docking_binary() is None:
        res["error"] = "no docking binary on PATH.\n" + docking.SETUP
        return res

    if model_text is None:
        path, info = fetch_model(client, uniprot, out_dir)
        if path is None:
            res["error"] = info
            return res
        model_text = path.read_text(encoding="utf-8")
        res["model"]["source"] = info
    if template_text is None:
        tpath = docking.fetch_pdb(client, template_pdb, out_dir)
        if tpath is None:
            res["error"] = f"could not fetch template {template_pdb}"
            return res
        template_text = Path(tpath).read_text(encoding="utf-8")

    assessed = docking.assess_model(model_text, uniprot)
    res["model"].update({"mean_plddt": assessed.mean_plddt,
                         "ordered_fraction": assessed.ordered_fraction})
    if assessed.dockable:
        res["caveats"].extend(assessed.notes)
    else:
        # `docking.screen_by_docking` refuses a model like this outright. The
        # pocket here is judged on its own residues instead, so say so: the
        # standing advice not to dock the model would contradict the run.
        res["caveats"].append(
            f"the model as a whole is mostly disordered: mean pLDDT "
            f"{assessed.mean_plddt}, {assessed.ordered_fraction:.0%} of "
            f"residues confidently folded. The normal docking stage would "
            f"refuse it. This run docks only into the transplanted pocket, "
            f"which is judged on its own residues, and leaves residues "
            f"below pLDDT {MIN_RESIDUE_PLDDT:.0f} out of the receptor.")

    site = transplant_site(model_text, template_text,
                           template_chain=template_chain,
                           template_ligand=template_ligand)
    site.update({"template_pdb": template_pdb,
                 "template_chain": template_chain,
                 "template_ligand": template_ligand})
    res["site"] = site
    if site.get("error"):
        res["error"] = site["error"]
        return res

    n_clash, n_lig = site["n_clashing_atoms"], site["n_ligand_atoms"]
    if n_clash:
        res["caveats"].append(
            f"{n_clash} of {n_lig} template ligand atoms sit within "
            f"{CLASH_DISTANCE} A of model atoms, so the pocket is narrower in "
            f"the unbound model than in the template. Large compounds may "
            f"score poorly for that reason.")
    dfg = site.get("aligned_dfg")
    if dfg and not re.fullmatch(DFG_LIKE, dfg):
        res["caveats"].append(
            f"the model has {dfg} where the template has "
            f"{site['template_dfg']}. A variant motif can mean the domain "
            f"does not bind ATP-site ligands the way the template does; "
            f"whether this one does is not tested here.")

    weak = {n for n, v in docking.residue_plddt(model_text).items()
            if v < MIN_RESIDUE_PLDDT}
    lost = [n for n in site["pocket_residues"] if n in weak]
    if lost:
        res["caveats"].append(
            f"{len(lost)} pocket residue(s) are below pLDDT "
            f"{MIN_RESIDUE_PLDDT:.0f} and are left out of the receptor, so "
            f"the pocket is incomplete: {lost}")
    kept = chr(10).join(
        line for line in model_text.splitlines()
        if not (line.startswith("ATOM") and int(line[22:26]) in weak))
    receptor = docking.prepare_receptor(kept, out_dir / f"{gene}_receptor.pdb")
    res["model"]["residues_left_out"] = len(weak)
    names, smiles = [], []
    for row in library:
        smi = row["smiles"]
        if strip_counter_ions:
            try:
                smi = docking.split_counter_ions(smi)["smiles"]
            except ValueError:
                pass
        names.append(row["name"])
        smiles.append(smi)

    dock = dict(exhaustiveness=exhaustiveness, seed=seed, score_fn=score_fn,
                cpu=cpu, jobs=jobs)
    window = (round(SIZE_WINDOW[0] * n_lig), round(SIZE_WINDOW[1] * n_lig))
    tally: dict = {}
    rows, engine = _dock_rows(names, smiles, receptor, site["box"],
                              tally=tally, **dock)
    got = _summarise(rows, window, control_class)
    passed = _control_passed(got["class_control"])
    res.update({"engine": engine, "n_library": len(library),
                "n_scored": len(rows), "failed": tally,
                "size_window": list(window), "n_ranked": got["n_ranked"],
                "size_fit": got["fit"], "median_affinity": got["median"],
                "ranking": got["ranked"][:top_n], "all_scores": rows,
                "class_control": got["class_control"],
                "control_passed": passed, "usable": bool(got["ranked"])})
    if passed is False:
        res["caveats"].append(
            f"the class check failed (size corrected AUC "
            f"{got['class_control'].get('auc_size_corrected')}, need "
            f"{MIN_CONTROL_AUC}): approved {control_class}s do not rank "
            f"above the rest, so this pocket does not look like the one they "
            f"bind. The ranking below is not a list of compounds that fit.")
    if not got["ranked"]:
        res["error"] = "nothing docked; see `failed` for the counts"
    elif reference:
        res["reference"] = _reference_control(
            template_text, template_pdb, template_chain, template_ligand,
            names, smiles, window, control_class, out_dir, **dock)
    return res


def _dock_rows(names: list, smiles: list, receptor, box: list, *,
               exhaustiveness: int, seed: int, score_fn, cpu: int, jobs: int,
               tally: dict | None = None) -> tuple:
    """Dock every compound into one box: (rows, engine). Failures are left out
    and counted in `tally`."""
    from . import systemic
    scored = docking._dock_smiles(
        smiles, receptor, None, exhaustiveness=exhaustiveness, seed=seed,
        score_fn=score_fn, box=tuple(box), cpu=cpu, jobs=jobs,
        detailed=True, tally=tally)
    rows, engine = [], None
    for name, smi, got in zip(names, smiles, scored):
        if got is None:
            continue
        engine = engine or got.engine
        cls = systemic.classify(name)
        rows.append({"name": name, "smiles": smi,
                     "affinity": got.affinity, "cnn_score": got.cnn_score,
                     "cnn_affinity": got.cnn_affinity,
                     "heavy_atoms": _heavy_atom_count(smi),
                     "class": cls.name if cls else None,
                     "mitotic_disruptor": systemic.is_mitotic_disruptor(name)})
    return rows, engine


def _summarise(rows: list, window: tuple, control_class: str) -> dict:
    """Rank inside the size window and run the class check on that set."""
    ranked = rank_by_residual(rows, window)
    inside = [r["affinity"] for r in ranked["ranked"]]
    return {"ranked": ranked["ranked"], "n_ranked": len(ranked["ranked"]),
            "fit": ranked["fit"],
            "median": round(float(np.median(inside)), 2) if inside else None,
            "class_control": class_control(ranked["ranked"], control_class)}


def _control_passed(cc: dict) -> bool | None:
    """True or False when the class check could be read, None when it could
    not (too few class members), which is not a pass."""
    auc = cc.get("auc_size_corrected")
    return None if auc is None else auc >= MIN_CONTROL_AUC


def _reference_control(template_text: str, template_pdb: str, chain: str,
                       ligand: str, names: list, smiles: list, window: tuple,
                       control_class: str, out_dir, **dock) -> dict:
    """The same compounds in the template's own pocket.

    The template holds the ligand in an open pocket, so this shows what the
    method does when the pocket is known to fit. If the class check passes
    here and fails on the model, the model's pocket is closed to those
    compounds; if it fails here too, the scoring cannot see the class at all
    and nothing about the model can be read from it.
    """
    out = {"template": template_pdb, "error": None}
    lig = ligand_atoms(template_text, ligand, chain)
    if not len(lig):
        out["error"] = f"template has no {ligand} on chain {chain}"
        return out
    receptor = docking.prepare_receptor(
        template_text, Path(out_dir) / f"{template_pdb}_reference.pdb",
        keep_chain=chain)
    rows, _ = _dock_rows(names, smiles, receptor, box_around(lig), **dock)
    got = _summarise(rows, window, control_class)
    out.update({"n_scored": len(rows), "n_ranked": got["n_ranked"],
                "size_fit": got["fit"], "median_affinity": got["median"],
                "class_control": got["class_control"],
                "control_passed": _control_passed(got["class_control"])})
    return out


# ------------------------------------------------------------------ output
def markdown(res: dict) -> str:
    L = [f"# Exploratory docking: {res['gene']}", "",
         f"**{res['banner']}**", ""]
    if res.get("error"):
        L += [f"Not run: {res['error']}", ""]
    s = res.get("site") or {}
    if s:
        L += ["## Pocket", "",
              f"Placed by transplanting {s.get('template_ligand')} from "
              f"{s.get('template_pdb')} chain {s.get('template_chain')}.", ""]
        for key, label in (("n_aligned", "residues aligned"),
                           ("identity", "sequence identity"),
                           ("core_rmsd", "superposition error (A)"),
                           ("pocket_mean_plddt", "pocket mean pLDDT"),
                           ("n_clashing_atoms", "template ligand atoms "
                                                "inside model atoms"),
                           ("template_dfg", "template DFG-like motif"),
                           ("aligned_dfg", "model residues at that motif")):
            if s.get(key) is not None:
                L.append(f"- {label}: {s[key]}")
        L.append("")
    for c in res.get("caveats") or []:
        L.append(f"- {c}")
    if res.get("caveats"):
        L.append("")
    cc = res.get("class_control") or {}
    if cc:
        L += ["## Class check", "", f"Class: {cc['class']} "
              f"({cc['n_members']} of {cc['n_compared']} compounds).", ""]
        if res.get("control_passed") is False:
            L += ["**The class check failed.** Approved members of this class "
                  "do not rank above the rest, so the pocket does not look "
                  "like the one they bind. Read the ranking below as scoring "
                  "noise from a pocket the compounds do not fit, not as a "
                  "list of what fits.", ""]
        elif res.get("control_passed") is None:
            L += ["The class check could not be read, so nothing here says "
                  "the pocket resembles a real one.", ""]
        for key, label in (("auc_raw_score", "AUC, raw score"),
                           ("auc_size_corrected", "AUC, size corrected"),
                           ("ef1_size_corrected", "EF1, size corrected")):
            if cc.get(key) is not None:
                L.append(f"- {label}: {cc[key]}")
        L += ["", cc.get("note", ""), ""]
    ref = res.get("reference")
    if ref:
        L += ["## The same compounds in the template's own pocket", ""]
        if ref.get("error"):
            L += [f"Not run: {ref['error']}", ""]
        else:
            rc = ref["class_control"]
            L += [f"Docked into {ref['template']} with its ligand's box. "
                  f"Median score {ref['median_affinity']} against "
                  f"{res.get('median_affinity')} in the model.", ""]
            for key, label in (("auc_raw_score", "AUC, raw score"),
                               ("auc_size_corrected", "AUC, size corrected")):
                if rc.get(key) is not None:
                    L.append(f"- {label}: {rc[key]}")
            L.append("")
            if ref["control_passed"] and res.get("control_passed") is False:
                L += ["The class is recovered in the open pocket and not in "
                      "the model, so the model's pocket is closed to these "
                      "compounds. That is a property of the model, not a "
                      "finding about the protein.", ""]
            elif ref["control_passed"] is False:
                L += ["The class is not recovered here either, so this "
                      "scoring cannot see the class and nothing about the "
                      "model can be read from it.", ""]
    if res.get("ranking"):
        fit = res.get("size_fit")
        lo, hi = res.get("size_window") or (None, None)
        L += ["## Best fit for size"
              + ("" if res.get("control_passed") is not False
                 else " (the class check failed; not evidence of fit)"), "",
              f"Ranked: {res.get('n_ranked')} of {res.get('n_scored')} "
              f"compounds that docked, those with {lo} to {hi} heavy atoms "
              f"(the template ligand has {res['site']['n_ligand_atoms']}). "
              f"Larger and smaller compounds keep their scores in the JSON "
              f"and are not ranked.", ""]
        if fit:
            L += [f"Score is fitted to heavy atoms over {fit['n']} compounds "
                  f"(slope {fit['slope']} per atom); the rank is by the "
                  f"residual, so size alone does not win.", ""]
        L += ["| Rank | Compound | Score | Heavy atoms | Residual | "
              "Efficiency | Notes |", "|---|---|---|---|---|---|---|"]
        for r in res["ranking"]:
            note = "; ".join(x for x in (
                r.get("class"),
                "mitotic disruptor" if r.get("mitotic_disruptor") else None)
                if x)
            L.append(f"| {r['rank']} | {r['name'].replace('|', '/')} | "
                     f"{r['affinity']} | {r['heavy_atoms']} | "
                     f"{r['residual']} | {r['ligand_efficiency']} | "
                     f"{note} |")
        L.append("")
    return "\n".join(L)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m promoter_tf_drugs.exploratory",
        description=("Dock the approved library into a pocket transplanted "
                     "onto an AlphaFold model. The result is not validated."))
    p.add_argument("--gene", required=True)
    p.add_argument("--uniprot", required=True)
    p.add_argument("--template-pdb", required=True,
                   help="a related protein's structure that holds a ligand")
    p.add_argument("--template-chain", required=True)
    p.add_argument("--template-ligand", required=True,
                   help="residue name of that ligand, for example ADP")
    p.add_argument("--library", default=None,
                   help="drug table; the bundled approved library by default")
    p.add_argument("--library-limit", type=int, default=100000)
    p.add_argument("--exhaustiveness", type=int, default=16)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--cpu", type=int, default=1)
    p.add_argument("--strip-counter-ions", action="store_true")
    p.add_argument("--reference", action="store_true",
                   help="also dock the same compounds into the template's own "
                        "pocket, to calibrate the class check (doubles the "
                        "time)")
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--out", default="results/exploratory")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    if args.library_limit < 1 or args.exhaustiveness < 1 or args.top < 1:
        p.error("--library-limit, --exhaustiveness and --top must be at "
                "least 1")

    from . import localdrugs as ld
    from .http import Client
    client = Client(offline=args.offline)
    rows = ld.library(args.library)
    if rows:
        library = [{"name": d["name"], "smiles": d["smiles"]} for d in rows]
    else:
        from .screening import approved_library
        library = approved_library(client)
    library = docking._limit_library(library, args.library_limit)

    out = Path(args.out) / args.gene
    res = run(args.gene, args.uniprot, library, client=client,
              template_pdb=args.template_pdb,
              template_chain=args.template_chain,
              template_ligand=args.template_ligand, out_dir=out,
              exhaustiveness=args.exhaustiveness, seed=args.seed,
              jobs=args.jobs, cpu=args.cpu,
              strip_counter_ions=args.strip_counter_ions, top_n=args.top,
              reference=args.reference)
    out.mkdir(parents=True, exist_ok=True)
    (out / "exploratory.json").write_text(
        json.dumps(res, indent=2, default=str), encoding="utf-8")
    (out / "exploratory.md").write_text(markdown(res), encoding="utf-8")
    print(markdown(res))
    print(f"Wrote {out / 'exploratory.json'} and {out / 'exploratory.md'}")
    return 1 if res.get("error") else 0


if __name__ == "__main__":
    raise SystemExit(main())
