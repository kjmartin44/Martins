"""Exploratory docking into a pocket transplanted onto a predicted model.

The result is unvalidated by design, so these tests pin what keeps it honest:
the site is only placed when the model and template agree, the ranking is not
a ranking of molecular size, and every output carries the warning.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from promoter_tf_drugs import docking, exploratory as ex

AA = "ACDEFGHIKLMNPQRSTVWY"
_ONE_TO_THREE = {v: k for k, v in ex._THREE.items()}


def _sequence(n=160, seed=1, motif=None):
    rng = np.random.default_rng(seed)
    s = [AA[i] for i in rng.integers(0, 20, n)]
    s = [c if c not in "DG" else "A" for c in s]      # keep motifs unique
    if motif:
        s[100:103] = list(motif)
    return "".join(s)


def _coords(n=160, seed=2):
    rng = np.random.default_rng(seed)
    steps = rng.normal(size=(n, 3))
    steps = steps / np.linalg.norm(steps, axis=1)[:, None] * 3.8
    return np.cumsum(steps, axis=0)


def _rotation(seed=3):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q if np.linalg.det(q) > 0 else -q


def _atom(kind, serial, name, res, chain, num, xyz, b=90.0, element="C"):
    return (f"{kind:<6}{serial:>5} {name:<4} {res:>3} {chain}{num:>4}    "
            f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}{1.0:>6.2f}"
            f"{b:>6.2f}          {element:>2}")


def _pdb(seq, xyz, chain="A", plddt=90.0, ligand=None, first=1):
    lines = [_atom("ATOM", i + 1, "CA", _ONE_TO_THREE[a], chain, first + i,
                   xyz[i], plddt(first + i) if callable(plddt) else plddt)
             for i, a in enumerate(seq)]
    for k, p in enumerate(ligand if ligand is not None else []):
        lines.append(_atom("HETATM", 5000 + k, "C1", "ADP", chain, 900, p, 20.0))
    return "\n".join(lines + ["END"]) + "\n"


def _pair(motif_model="DFG", motif_template="DFG", plddt=90.0):
    """A model and a template of one fold, the template moved and holding a
    ligand beside residue 100. Returns the texts and the ligand's expected
    centre in the model's frame."""
    seq_t = _sequence(motif=motif_template)
    seq_m = _sequence(motif=motif_model)
    xyz = _coords()
    r, t = _rotation(), np.array([40.0, -25.0, 10.0])
    lig = np.array([xyz[99] + [2.6, 0, 0], xyz[100] + [0, 2.6, 0],
                    xyz[101] + [0, 0, 2.6], xyz[100] + [-2.6, 0, 0],
                    xyz[99] + [0, -2.6, 0], xyz[101] + [0, 0, -2.6]])
    template = _pdb(seq_t, xyz @ r.T + t, ligand=lig @ r.T + t)
    model = _pdb(seq_m, xyz, plddt=plddt)
    return model, template, lig.mean(0)


# -------------------------------------------------------------- alignment
def test_identical_sequences_align_one_to_one():
    s = _sequence()
    assert ex.align(s, s) == [(i, i) for i in range(len(s))]


def test_an_insertion_is_bridged_not_fatal():
    s = _sequence()
    longer = s[:70] + "WWWWW" + s[70:]
    pairs = dict(ex.align(longer, s))
    assert pairs[0] == 0 and pairs[len(longer) - 1] == len(s) - 1
    assert pairs[75] == 70, "residues after the insertion map back"


# ---------------------------------------------------------- superposition
def test_a_rigid_move_is_recovered_exactly():
    p = _coords()
    r, t = _rotation(), np.array([5.0, 6.0, -7.0])
    q = (p - t) @ r      # so that q @ r.T + t == p
    fit = ex.superpose(p, q)
    assert fit["rmsd"] < 1e-6 and fit["n_core"] == len(p)
    assert np.allclose(q @ fit["r"].T + fit["t"], p, atol=1e-6)


def test_residues_that_moved_do_not_drag_the_fit():
    p = _coords()
    q = p.copy()
    q[::5] += 25.0                     # a fifth of the chain has moved away
    fit = ex.superpose(p, q)
    assert fit["rmsd"] < 0.5
    assert fit["n_core"] < len(p), "the moved residues are left out of the core"


# ------------------------------------------------------------------- site
def test_template_ligand_hydrogens_and_alternate_positions_are_not_atoms():
    def het(serial, name, element, alt=" ", x=0.0):
        line = _atom("HETATM", serial, name, "ADP", "A", 900, (x, 0, 0),
                     element=element)
        return line[:16] + alt + line[17:]
    text = chr(10).join([het(1, "C1", "C"), het(2, "N1", "N", x=1.5),
                         het(3, "H1", "H", x=2.0),
                         het(4, "C2", "C", alt="B", x=3.0),
                         het(5, "C3", "C", alt="A", x=4.0)]) + chr(10)
    got = ex.ligand_atoms(text, "ADP", "A")
    assert got.shape == (3, 3), "two heavy atoms and the A position; no H, no B"
    assert sorted(got[:, 0]) == [0.0, 1.5, 4.0]



def test_the_box_lands_on_the_template_ligand_in_the_model_frame():
    model, template, centre = _pair()
    site = ex.transplant_site(model, template, template_chain="A",
                              template_ligand="ADP")
    assert "error" not in site
    assert site["n_aligned"] == 160 and site["identity"] == 1.0
    assert site["core_rmsd"] < 0.1
    assert np.allclose(site["box"][:3], centre, atol=0.1)
    assert site["box"][3] >= 2 * ex.BOX_PADDING, "padding is added to the extent"
    assert 100 + 0 in site["pocket_residues"] or 101 in site["pocket_residues"]
    assert site["pocket_mean_plddt"] == 90.0


def test_a_poorly_modelled_pocket_is_refused():
    model, template, _ = _pair(plddt=50.0)
    site = ex.transplant_site(model, template, template_chain="A",
                              template_ligand="ADP")
    assert "pLDDT" in site["error"] and "box" in site, \
        "the box is still reported so the refusal can be inspected"


def test_a_template_without_the_ligand_is_refused():
    model, template, _ = _pair()
    site = ex.transplant_site(model, template, template_chain="A",
                              template_ligand="ATP")
    assert "no ATP" in site["error"]
    site = ex.transplant_site(model, template, template_chain="B",
                              template_ligand="ADP")
    assert "no residues on chain B" in site["error"]
    lone = _pdb(_sequence(), _coords())          # residues but no ligand
    assert "no ADP on chain A" in ex.transplant_site(
        model, lone, template_chain="A", template_ligand="ADP")["error"]


def test_an_unrelated_template_is_not_used():
    model, template, _ = _pair()
    other = _pdb(_sequence(seed=99), _coords(seed=7),
                 ligand=np.array([[1.0, 2, 3]] * 6))
    site = ex.transplant_site(model, other, template_chain="A",
                              template_ligand="ADP")
    assert "error" in site and "box" not in site


def test_a_template_of_the_same_fold_but_distant_sequence_is_refused():
    model, template, _ = _pair()
    swap = str.maketrans("ILVMFYWKRHNQSTAGE", "LVMIYWFRHKQNTSGAD")
    lines = []
    for line in template.splitlines():
        if line.startswith("ATOM"):
            aa = ex._THREE[line[17:20]].translate(swap)
            line = line[:17] + _ONE_TO_THREE[aa] + line[20:]
        lines.append(line)
    site = ex.transplant_site(model, chr(10).join(lines) + chr(10),
                              template_chain="A", template_ligand="ADP")
    assert site["n_aligned"] >= ex.MIN_ALIGNED
    assert "identical" in site["error"] and "box" not in site


def test_a_fold_that_does_not_match_is_refused_even_when_sequences_do():
    model, template, _ = _pair()
    rng = np.random.default_rng(5)
    scrambled = _pdb(_sequence(), rng.normal(size=(160, 3)) * 30,
                     ligand=np.zeros((6, 3)))
    site = ex.transplant_site(model, scrambled, template_chain="A",
                              template_ligand="ADP")
    assert "superposition" in site["error"] and "box" not in site


# ---------------------------------------------------------------- ranking
def test_rank_is_by_score_for_size_not_by_size():
    rows = [{"name": f"c{n}", "affinity": -0.3 * n + (0.01 * (n % 3)),
             "heavy_atoms": n} for n in range(10, 40)]
    rows.append({"name": "tight", "affinity": -0.3 * 20 - 2.0,
                 "heavy_atoms": 20})
    out = ex.rank_by_residual(rows)
    assert out["ranked"][0]["name"] == "tight"
    biggest = max(rows, key=lambda r: r["heavy_atoms"])
    raw_best = min(rows, key=lambda r: r["affinity"])
    assert raw_best is biggest, "by raw score the largest molecule leads"
    assert out["ranked"][0]["rank"] == 1
    assert out["fit"]["slope"] < 0 and out["fit"]["n"] == len(rows)


def test_compounds_outside_the_size_window_are_scored_but_not_ranked():
    """In a narrow pocket big molecules score badly. A fit through every size
    then has a positive slope, and the residual rewards the one huge molecule
    that still scored well. Fitting only inside the window keeps it out."""
    rows = [{"name": f"c{n}", "affinity": -3.0 + 0.1 * (n % 5),
             "heavy_atoms": n} for n in range(15, 40)]
    rows += [{"name": f"clash{n}", "affinity": 10.0, "heavy_atoms": n}
             for n in (70, 75, 80, 85)]
    rows.append({"name": "bighit", "affinity": -6.0, "heavy_atoms": 90})
    rows.append({"name": "tiny", "affinity": -3.0, "heavy_atoms": 6})
    wide = ex.rank_by_residual([dict(r) for r in rows])
    assert wide["fit"]["slope"] > 0
    assert wide["ranked"][0]["name"] == "bighit", "the trap the window avoids"
    out = ex.rank_by_residual(rows, window=(14, 43))
    assert "bighit" not in [r["name"] for r in out["ranked"]]
    assert sorted(r["name"] for r in out["unranked"]) == [
        "bighit", "clash70", "clash75", "clash80", "clash85", "tiny"]
    assert all(r["in_window"] is False for r in out["unranked"])
    assert out["fit"]["n"] == 25, "the line is fitted to the window only"


def test_rows_without_a_score_or_size_are_left_out_not_ranked_last():
    rows = [{"name": "a", "affinity": None, "heavy_atoms": 20},
            {"name": "b", "affinity": -7.0, "heavy_atoms": None},
            {"name": "c", "affinity": -7.0, "heavy_atoms": 20}]
    out = ex.rank_by_residual(rows)
    assert [r["name"] for r in out["ranked"]] == ["c"]
    assert out["fit"] is None, "three compounds are too few to fit a size trend"


KINASE = ["imatinib", "gefitinib", "erlotinib", "dasatinib", "nilotinib",
          "sunitinib", "sorafenib", "lapatinib"]


def _library_rows(kinase_residual):
    rows = [{"name": n, "affinity": -9.0, "residual": kinase_residual}
            for n in KINASE]
    rows += [{"name": f"other{i}", "affinity": -7.0, "residual": 0.0 + i * .01}
             for i in range(60)]
    return rows


def test_the_class_check_reports_enrichment_when_the_class_rises():
    cc = ex.class_control(_library_rows(-1.0), "antineoplastic kinase inhibitor")
    assert cc["n_members"] == len(KINASE)
    assert cc["auc_size_corrected"] == 1.0 and cc["auc_raw_score"] == 1.0
    assert "does not show the target is inhibited" in cc["note"]


def test_the_class_check_reports_no_enrichment_when_it_does_not():
    cc = ex.class_control(_library_rows(+1.0), "antineoplastic kinase inhibitor")
    assert cc["auc_size_corrected"] == 0.0


def test_the_class_check_is_withheld_for_too_few_members():
    rows = _library_rows(-1.0)[:3] + _library_rows(-1.0)[8:]
    cc = ex.class_control(rows, "antineoplastic kinase inhibitor")
    assert cc["n_members"] == 3 and "auc_raw_score" not in cc
    assert "need 5" in cc["note"]


# ------------------------------------------------------------------ model
class _Client:
    def __init__(self, record):
        self.record, self.urls = record, []

    def get(self, url, params=None, **kw):
        self.urls.append(url)
        if isinstance(self.record, Exception):
            raise self.record
        return self.record

    def request(self, url, **kw):
        self.urls.append(url)
        return "ATOM\n"


def test_the_model_address_is_read_from_the_record_not_built_from_a_version(
        tmp_path):
    c = _Client([{"pdbUrl": "https://example.org/files/AF-P1-F1-model_v9.pdb"}])
    path, url = ex.fetch_model(c, "P1", tmp_path)
    assert path.name == "AF-P1-F1-model_v9.pdb" and url.endswith("v9.pdb")
    assert not any("v4" in u for u in c.urls)
    c.urls.clear()
    ex.fetch_model(c, "P1", tmp_path)
    assert c.urls == [docking.ALPHAFOLD_API.format(uniprot="P1")], \
        "the file is cached, only the record is asked for again"


def test_the_shared_fetch_also_follows_the_record(tmp_path):
    c = _Client({"pdbUrl": "https://example.org/files/AF-P2-F1-model_v9.pdb"})
    got = docking.fetch_alphafold(c, "P2", tmp_path)
    assert got.name == "AF-P2.pdb" and got.exists()
    assert not any("v4" in u for u in c.urls)
    assert docking.fetch_alphafold(_Client(RuntimeError("503")), "P3",
                                   tmp_path) is None


def test_a_failed_download_comes_back_as_text_not_an_exception(tmp_path):
    path, why = ex.fetch_model(_Client(RuntimeError("503")), "P1", tmp_path)
    assert path is None and "could not fetch" in why and "503" in why


# -------------------------------------------------------------------- run
def _library(n=24):
    return [{"name": f"drug{i}", "smiles": "C" * (4 + i)} for i in range(n)]


def _score(smiles):
    return -0.3 * len(smiles)


def test_a_run_labels_itself_unvalidated_and_serialises(tmp_path):
    model, template, _ = _pair()
    res = ex.run("G1", "P1", _library(), template_pdb="TMPL",
                 template_chain="A", template_ligand="ADP",
                 model_text=model, template_text=template, out_dir=tmp_path,
                 score_fn=_score)
    assert res["usable"] and res["status"] == "unvalidated"
    assert res["banner"].startswith("UNVALIDATED")
    assert res["n_scored"] == 24 and res["ranking"][0]["rank"] == 1
    assert json.loads(json.dumps(res, default=str))["banner"] == res["banner"]
    md = ex.markdown(res)
    assert md.count("UNVALIDATED") == 1 and "## Best fit for size" in md


def test_the_box_handed_to_the_docker_is_the_transplanted_one(tmp_path,
                                                              monkeypatch):
    model, template, centre = _pair()
    seen = {}
    real = docking._dock_smiles

    def spy(smiles, receptor, autobox, **kw):
        seen.update(box=kw["box"], detailed=kw["detailed"])
        return real(smiles, receptor, autobox, **kw)
    monkeypatch.setattr(docking, "_dock_smiles", spy)
    ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
           template_ligand="ADP", model_text=model, template_text=template,
           out_dir=tmp_path, score_fn=_score)
    assert np.allclose(seen["box"][:3], centre, atol=0.1)
    assert seen["detailed"] is True


def test_a_mostly_disordered_model_is_docked_only_at_the_pocket(tmp_path):
    model, template, _ = _pair()
    ordered = lambda n: 90.0 if 90 <= n <= 125 else 40.0
    text = _pdb(_sequence(), _coords(), plddt=ordered)
    res = ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
                 template_ligand="ADP", model_text=text,
                 template_text=template, out_dir=tmp_path, score_fn=_score)
    assert res["usable"], "the pocket is what is judged, not the whole model"
    said = " ".join(res["caveats"])
    assert "mostly disordered" in said and "would refuse it" in said
    assert "Do not dock this model" not in said
    kept = {int(line[22:26]) for line in
            (tmp_path / "G1_receptor.pdb").read_text().splitlines()
            if line.startswith("ATOM")}
    assert kept == set(range(90, 126)), "disordered residues are not docked into"
    assert res["model"]["residues_left_out"] == 160 - 36
    lost = [n for n in res["site"]["pocket_residues"] if not 90 <= n <= 125]
    assert bool(lost) == any("left out of the receptor, so" in c
                             for c in res["caveats"])


def test_a_variant_motif_is_reported_as_a_caveat(tmp_path):
    model, template, _ = _pair(motif_model="DFS", motif_template="DLG")
    res = ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
                 template_ligand="ADP", model_text=model,
                 template_text=template, out_dir=tmp_path, score_fn=_score)
    assert res["site"]["template_dfg"] == "DLG"
    assert res["site"]["aligned_dfg"] == "DFS"
    assert any("DFS where the template has DLG" in c for c in res["caveats"])


def test_a_conserved_motif_gives_no_motif_caveat(tmp_path):
    model, template, _ = _pair(motif_model="DFG", motif_template="DLG")
    res = ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
                 template_ligand="ADP", model_text=model,
                 template_text=template, out_dir=tmp_path, score_fn=_score)
    assert not any("motif" in c for c in res["caveats"])


def test_template_ligand_atoms_inside_the_model_are_counted_and_reported(
        tmp_path):
    model, template, _ = _pair()
    res = ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
                 template_ligand="ADP", model_text=model,
                 template_text=template, out_dir=tmp_path, score_fn=_score)
    n = res["site"]["n_clashing_atoms"]
    assert n > 0, "the synthetic ligand sits tight against the chain"
    assert any(f"{n} of 6 template ligand atoms" in c and "narrower" in c
               for c in res["caveats"])


def test_a_site_placed_outside_the_protein_is_refused():
    model, template, _ = _pair()
    far = chr(10).join(
        line[:30] + f"{float(line[30:38]) + 500:>8.3f}" + line[38:]
        if line.startswith("HETATM") else line
        for line in template.splitlines()) + chr(10)
    site = ex.transplant_site(model, far, template_chain="A",
                              template_ligand="ADP")
    assert "not in the protein" in site["error"]


def test_a_refused_site_docks_nothing_and_says_why(tmp_path, monkeypatch):
    model, template, _ = _pair(plddt=40.0)
    monkeypatch.setattr(docking, "_dock_smiles",
                        lambda *a, **k: pytest.fail("docked into a refused site"))
    res = ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
                 template_ligand="ADP", model_text=model,
                 template_text=template, out_dir=tmp_path, score_fn=_score)
    assert not res["usable"] and res["error"]
    assert "Not run:" in ex.markdown(res)
    assert res["banner"].startswith("UNVALIDATED")


def test_without_a_docking_binary_it_stops_before_any_download(tmp_path,
                                                               monkeypatch):
    monkeypatch.setattr(docking, "_docking_binary", lambda: None)
    res = ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
                 template_ligand="ADP", client=_Client(RuntimeError("no net")),
                 out_dir=tmp_path)
    assert "no docking binary" in res["error"] and not res["usable"]


def test_compounds_that_fail_to_dock_are_counted_not_ranked(tmp_path):
    model, template, _ = _pair()
    lib = _library()

    def score(smi):
        return None if len(smi) == 6 else _score(smi)
    res = ex.run("G1", "P1", lib, template_pdb="T", template_chain="A",
                 template_ligand="ADP", model_text=model,
                 template_text=template, out_dir=tmp_path, score_fn=score)
    assert res["n_scored"] == 23 and res["n_library"] == 24
    assert all(r["name"] != "drug2" for r in res["all_scores"])


# -------------------------------------------------------------------- cli
def test_the_command_line_rejects_nonsense_before_any_work(capsys):
    with pytest.raises(SystemExit):
        ex.main(["--gene", "G", "--uniprot", "P", "--template-pdb", "T",
                 "--template-chain", "A", "--template-ligand", "ADP",
                 "--top", "0"])
    assert "at least 1" in capsys.readouterr().err


# ------------------------------------------- the class check as a gate
def test_the_gate_reads_the_size_corrected_auc_and_not_a_missing_one_as_a_pass():
    assert ex._control_passed({"auc_size_corrected": 0.6}) is True
    assert ex._control_passed({"auc_size_corrected": 0.59}) is False
    assert ex._control_passed({"auc_size_corrected": 0.0}) is False
    assert ex._control_passed({"note": "only 3 members"}) is None


def _class_library():
    rows = [{"name": n, "smiles": "N" + "C" * (3 + i % 6)}
            for i, n in enumerate(KINASE)]
    rows += [{"name": f"other{i}", "smiles": "C" * (4 + i % 6)}
             for i in range(40)]
    return rows


def _likes_n(smi):
    return -6.0 if "N" in smi else -2.0


def _dislikes_n(smi):
    return -2.0 if "N" in smi else -6.0


def _run_class(tmp_path, score, **kw):
    model, template, _ = _pair()
    return ex.run("G1", "P1", _class_library(), template_pdb="T",
                  template_chain="A", template_ligand="ADP",
                  model_text=model, template_text=template, out_dir=tmp_path,
                  score_fn=score, **kw)


def test_a_failed_class_check_is_said_loudly_in_the_data_and_the_report(
        tmp_path):
    res = _run_class(tmp_path, _dislikes_n)
    assert res["class_control"]["n_members"] == len(KINASE)
    assert res["control_passed"] is False
    assert any("class check failed" in c and "not a list of compounds that fit"
               in c for c in res["caveats"])
    md = ex.markdown(res)
    assert "**The class check failed.**" in md
    assert "(the class check failed; not evidence of fit)" in md


def test_a_passed_class_check_carries_no_failure_text(tmp_path):
    res = _run_class(tmp_path, _likes_n)
    assert res["control_passed"] is True
    md = ex.markdown(res)
    assert "class check failed" not in md and "## Best fit for size\n" in md


def test_too_few_class_members_is_not_a_pass(tmp_path):
    model, template, _ = _pair()
    res = ex.run("G1", "P1", _library(), template_pdb="T", template_chain="A",
                 template_ligand="ADP", model_text=model,
                 template_text=template, out_dir=tmp_path, score_fn=_score)
    assert res["control_passed"] is None
    assert "could not be read" in ex.markdown(res)


# --------------------------------------------------- the reference pocket
def test_the_reference_docks_the_same_compounds_into_the_templates_own_box(
        tmp_path, monkeypatch):
    model, template, _ = _pair()
    calls = []
    real = docking._dock_smiles

    def spy(smiles, receptor, autobox, **kw):
        calls.append((len(smiles), str(receptor), kw["box"]))
        return real(smiles, receptor, autobox, **kw)
    monkeypatch.setattr(docking, "_dock_smiles", spy)
    res = _run_class(tmp_path, _likes_n, reference=True)
    assert len(calls) == 2 and calls[0][0] == calls[1][0]
    assert calls[1][1].endswith("T_reference.pdb")
    lig = ex.ligand_atoms(template, "ADP", "A")
    assert np.allclose(calls[1][2], ex.box_around(lig))
    assert not np.allclose(calls[0][2][:3], calls[1][2][:3]), \
        "the reference uses the template's frame, not the transplanted one"
    assert res["reference"]["control_passed"] is True
    text = (tmp_path / "T_reference.pdb").read_text()
    assert "HETATM" not in text, "the template's ligand is not left in"


def test_the_reference_receptor_is_the_template_chain_alone(tmp_path):
    model, template, _ = _pair()
    other = chr(10).join(line[:21] + "B" + line[22:] for line in
                         template.splitlines() if line.startswith("ATOM"))
    ex.run("G1", "P1", _class_library(), template_pdb="T", template_chain="A",
           template_ligand="ADP", model_text=model,
           template_text=template + other + chr(10), out_dir=tmp_path,
           score_fn=_likes_n, reference=True)
    chains = {line[21] for line in
              (tmp_path / "T_reference.pdb").read_text().splitlines()
              if line.startswith("ATOM")}
    assert chains == {"A"}


def test_the_reference_is_off_unless_asked_for(tmp_path):
    assert "reference" not in _run_class(tmp_path, _likes_n)


def test_the_report_says_whose_doing_a_failed_check_is():
    base = {"gene": "G", "banner": ex.BANNER, "caveats": [], "ranking": [],
            "class_control": {"class": "k", "n_members": 8, "n_compared": 40,
                              "auc_size_corrected": 0.2},
            "control_passed": False, "median_affinity": -1.0}
    ok = {"template": "T", "error": None, "median_affinity": -6.0,
          "control_passed": True,
          "class_control": {"auc_size_corrected": 0.9}}
    md = ex.markdown({**base, "reference": ok})
    assert "recovered in the open pocket and not in the model" in md
    assert "a property of the model, not a finding about the protein" in md
    bad = {**ok, "control_passed": False,
           "class_control": {"auc_size_corrected": 0.4}}
    assert "cannot see the class" in ex.markdown({**base, "reference": bad})
    gone = {"template": "T", "error": "template has no ADP on chain A"}
    assert "Not run: template has no ADP" in ex.markdown(
        {**base, "reference": gone})


def test_the_command_line_passes_the_reference_flag_through(tmp_path,
                                                            monkeypatch):
    seen = {}

    def fake(gene, uniprot, library, **kw):
        seen.update(kw)
        return {"gene": gene, "banner": ex.BANNER, "error": None,
                "caveats": []}
    monkeypatch.setattr(ex, "run", fake)
    code = ex.main(["--gene", "G", "--uniprot", "P", "--template-pdb", "T",
                    "--template-chain", "A", "--template-ligand", "ADP",
                    "--reference", "--library-limit", "5", "--offline",
                    "--out", str(tmp_path)])
    assert code == 0 and seen["reference"] is True
    assert (tmp_path / "G" / "exploratory.json").exists()
