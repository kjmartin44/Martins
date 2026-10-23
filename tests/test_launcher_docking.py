"""Direct tests for the docking and WSL functions the launcher calls.

The suite did not enter these: the docking binary, WSL and the PDBQT
converters are external, so each test replaces the one process or file
boundary with a fake and checks what the function does on its own side of
it. Nothing here starts a program or touches the network.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from generoutes import docking as dk
from generoutes import wsldock

NL = chr(10)

# Two poses as gnina writes them: tags alone on a line, value on the next.
_GNINA_SDF = NL.join([
    "pose1", "", "", "  0  0  0     0  0  0  0  0  0999 V2000", "M  END",
    "> <minimizedAffinity>", "-7.5", "",
    "> <CNNscore>", "0.81", "",
    "> <CNNaffinity>", "6.2", "",
    "$$$$",
    "pose2", "", "", "  0  0  0     0  0  0  0  0  0999 V2000", "M  END",
    "> <minimizedAffinity>", "-6.0", "",
    "> <CNNscore>", "0.40", "",
    "> <CNNaffinity>", "5.1", "",
    "$$$$", ""])


def _lf_flag(cmd, name):
    return cmd[cmd.index(name) + 1]


class _LfRun:
    """A stand-in for subprocess.run that records the command and writes
    whatever the real program would have written to its output file."""

    def __init__(self, *, out_text="", out_flag="-o", stdout="",
                 returncode=0, stderr=""):
        self.calls = []
        self.out_text, self.out_flag = out_text, out_flag
        self.stdout, self.returncode, self.stderr = stdout, returncode, stderr

    def __call__(self, cmd, **kw):
        self.calls.append(list(cmd))
        if self.returncode == 0 and self.out_flag in cmd:
            Path(_lf_flag(cmd, self.out_flag)).write_text(
                self.out_text, encoding="utf-8")
        return subprocess.CompletedProcess(
            cmd, self.returncode, self.stdout, self.stderr)


@pytest.fixture
def fake_gnina(monkeypatch):
    monkeypatch.setattr(dk, "_docking_binary", lambda: "C:/fake/gnina.exe")
    run = _LfRun(out_text=_GNINA_SDF)
    monkeypatch.setattr(dk.subprocess, "run", run)
    return run


# --------------------------------------------------------------- _run_dock
def test_run_dock_returns_the_programs_stdout_and_its_output_file(
        monkeypatch, tmp_path):
    run = _LfRun(out_text="SDF TEXT", stdout="table")
    monkeypatch.setattr(dk.subprocess, "run", run)
    out = tmp_path / "o.sdf"
    got = dk._run_dock("C:/fake/smina.exe", ["-o", str(out)], out)
    assert got == ("table", "SDF TEXT")
    assert run.calls[0][0] == "C:/fake/smina.exe"


def test_run_dock_reports_failure_as_none_not_as_an_empty_score(
        monkeypatch, tmp_path, caplog):
    import logging
    monkeypatch.setattr(dk.subprocess, "run",
                        _LfRun(returncode=1, stderr="cuda init failed"))
    with caplog.at_level(logging.WARNING):
        got = dk._run_dock("C:/fake/gnina.exe", [], tmp_path / "o.sdf")
    assert got == (None, "")
    assert any("cuda init failed" in r.message for r in caplog.records)


def test_run_dock_treats_a_timeout_as_a_failure(monkeypatch, tmp_path):
    def slow(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, 1)
    monkeypatch.setattr(dk.subprocess, "run", slow)
    assert dk._run_dock("C:/fake/gnina.exe", [], tmp_path / "o.sdf") \
        == (None, "")


def test_run_dock_sends_a_posix_binary_through_wsl(monkeypatch, tmp_path):
    """An absolute POSIX path is a WSL program; handing it to a Windows
    subprocess raises FileNotFoundError and scores every compound None."""
    seen = {}
    out = tmp_path / "o.sdf"
    out.write_text("FROM WSL", encoding="utf-8")
    monkeypatch.setattr(wsldock, "dock_script",
                        lambda b, a: seen.setdefault("script", f"{b} {a}"))

    def fake_script(script, timeout=0):
        seen["ran"] = script
        return subprocess.CompletedProcess([], 0, "stdout", "")
    monkeypatch.setattr(wsldock, "run_script", fake_script)
    got = dk._run_dock("/home/u/bin/gnina", ["-r", "x"], out)
    assert got == ("stdout", "FROM WSL")
    assert "/home/u/bin/gnina" in seen["ran"]


def test_run_dock_survives_wsl_itself_failing(monkeypatch, tmp_path):
    monkeypatch.setattr(wsldock, "dock_script", lambda b, a: "s")

    def boom(script, timeout=0):
        raise OSError("wsl.exe vanished")
    monkeypatch.setattr(wsldock, "run_script", boom)
    assert dk._run_dock("/home/u/bin/gnina", [], tmp_path / "o") == (None, "")


# ------------------------------------------------------------ dock_one*
def test_dock_one_returns_the_best_affinity_among_the_poses(fake_gnina,
                                                           tmp_path):
    got = dk.dock_one(tmp_path / "r.pdb", tmp_path / "l.sdf",
                      autobox_ligand=tmp_path / "ref.sdf", exhaustiveness=4,
                      seed=7, cpu=2)
    assert got == -7.5
    cmd = fake_gnina.calls[0]
    assert _lf_flag(cmd, "--exhaustiveness") == "4"
    assert _lf_flag(cmd, "--seed") == "7"
    assert _lf_flag(cmd, "--cpu") == "2"
    assert "--autobox_ligand" in cmd, "the box must sit on the known ligand"
    assert _lf_flag(cmd, "--cnn_scoring") == "rescore", "gnina rescoring"


def test_dock_one_without_a_reference_ligand_passes_no_box(fake_gnina,
                                                           tmp_path):
    dk.dock_one(tmp_path / "r.pdb", tmp_path / "l.sdf")
    assert "--autobox_ligand" not in fake_gnina.calls[0]


def test_dock_one_says_how_to_set_up_when_there_is_no_binary(monkeypatch,
                                                            tmp_path):
    monkeypatch.setattr(dk, "_docking_binary", lambda: None)
    with pytest.raises(RuntimeError, match="no docking binary"):
        dk.dock_one(tmp_path / "r.pdb", tmp_path / "l.sdf")


def test_dock_one_is_none_when_the_program_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(dk, "_docking_binary", lambda: "C:/fake/gnina.exe")
    monkeypatch.setattr(dk.subprocess, "run", _LfRun(returncode=2))
    assert dk.dock_one(tmp_path / "r.pdb", tmp_path / "l.sdf") is None


def test_dock_one_scored_returns_all_three_numbers_of_the_first_pose(
        fake_gnina, tmp_path):
    got = dk.dock_one_scored(tmp_path / "r.pdb", tmp_path / "l.sdf",
                             (1, 2, 3, 20, 20, 20))
    assert (got.affinity, got.cnn_score, got.cnn_affinity) == (-7.5, 0.81, 6.2)
    assert got.engine == "gnina"
    cmd = fake_gnina.calls[0]
    assert _lf_flag(cmd, "--center_x") == "1.000"
    assert _lf_flag(cmd, "--size_z") == "20.000"


def test_dock_one_scored_is_none_without_poses_or_on_failure(monkeypatch,
                                                            tmp_path):
    box = (0, 0, 0, 10, 10, 10)
    monkeypatch.setattr(dk.subprocess, "run", _LfRun(out_text="no tags"))
    assert dk.dock_one_scored(tmp_path / "r", tmp_path / "l", box,
                              binary="C:/fake/gnina.exe") is None
    monkeypatch.setattr(dk.subprocess, "run", _LfRun(returncode=1))
    assert dk.dock_one_scored(tmp_path / "r", tmp_path / "l", box,
                              binary="C:/fake/gnina.exe") is None


def test_dock_one_scored_names_the_engine_from_the_binary(monkeypatch,
                                                         tmp_path):
    monkeypatch.setattr(
        dk.subprocess, "run",
        _LfRun(out_text=NL.join(["> <minimizedAffinity>", "-5.5", ""])))
    got = dk.dock_one_scored(tmp_path / "r", tmp_path / "l",
                             (0, 0, 0, 10, 10, 10),
                             binary="C:/fake/smina.exe")
    assert got.engine == "smina" and got.affinity == -5.5
    assert got.cnn_score is None, "a missing score is None, not 0"


def test_dock_one_scored_needs_a_binary(monkeypatch, tmp_path):
    monkeypatch.setattr(dk, "_docking_binary", lambda: None)
    with pytest.raises(RuntimeError, match="no docking binary"):
        dk.dock_one_scored(tmp_path / "r", tmp_path / "l", (0,) * 6)


# ------------------------------------------------------------- dock_one_vina
@pytest.fixture
def vina_world(monkeypatch, tmp_path):
    """vina on PATH, a receptor converter and a ligand preparer that work."""
    monkeypatch.setattr(
        dk.shutil, "which",
        lambda n, *a, **k: "C:/fake/vina.exe" if n == "vina" else None)
    state = {"receptor_calls": 0, "ligand": "LIGAND PDBQT"}

    def receptor(pdb, out):
        state["receptor_calls"] += 1
        Path(out).write_text("REC", encoding="utf-8")
        return Path(out)
    monkeypatch.setattr(dk, "_receptor_pdbqt", receptor)
    monkeypatch.setattr(dk, "_meeko_ligand_pdbqt",
                        lambda smi, seed=0: state["ligand"])
    run = _LfRun(out_flag="--out", stdout="",
                 out_text=NL.join([
                     "REMARK VINA RESULT:   -8.2      0.000      0.000",
                     "REMARK VINA RESULT:   -7.1      1.200      2.000", ""]))
    monkeypatch.setattr(dk.subprocess, "run", run)
    state["run"] = run
    state["pdb"] = tmp_path / "rec.pdb"
    state["pdb"].write_text("ATOM" + NL, encoding="utf-8")
    return state


def test_dock_one_vina_returns_the_best_result_line(vina_world):
    got = dk.dock_one_vina(vina_world["pdb"], "CCO", (0, 0, 0, 10, 10, 10),
                           exhaustiveness=3, seed=5, cpu=2)
    assert got == -8.2
    cmd = vina_world["run"].calls[0]
    assert _lf_flag(cmd, "--exhaustiveness") == "3"
    assert _lf_flag(cmd, "--size_x") == "10.000"


def test_dock_one_vina_converts_the_receptor_once_per_screen(vina_world):
    cache = {}
    box = (0, 0, 0, 10, 10, 10)
    dk.dock_one_vina(vina_world["pdb"], "CCO", box,
                     receptor_pdbqt_cache=cache)
    dk.dock_one_vina(vina_world["pdb"], "CCN", box,
                     receptor_pdbqt_cache=cache)
    assert vina_world["receptor_calls"] == 1


def test_the_bundle_path_passes_jobs_through():
    """`--jobs` is the documented speed lever - its own help says
    "--jobs 8 --cpu 4 beats --cpu 32 by about 2.5x".

    `dock_from_bundle` had no `jobs` parameter at all, so the whole
    "prepare here, dock there" route - the one that needs WSL and takes
    longest - silently ran one compound at a time against a flag the CLI
    accepted without a word.
    """
    import inspect

    assert "jobs" in inspect.signature(dk.dock_from_bundle).parameters, (
        "the bundle path has to accept it before it can honour it")
    body = inspect.getsource(dk.dock_from_bundle)
    assert "jobs=jobs" in body, "and hand it to the screen"

    cli = inspect.getsource(dk.main)
    i = cli.index("dock_from_bundle(")
    assert "jobs=args.jobs" in cli[i:i + 400], "and the CLI has to pass it"


def test_parallel_compounds_never_convert_the_receptor_at_the_same_time(
        vina_world, monkeypatch):
    """`_dock_smiles` runs compounds through a ThreadPoolExecutor.

    The first --jobs threads miss the cache together and used to all
    launch a converter onto one fixed path. The only guard was that the
    output file existed, which is true the moment any writer creates it,
    so a thread could hand Vina a receptor another process was still
    writing - and Vina exits 0 on a short receptor, so that reaches the
    gate and the ranking as a wrong number rather than as an error.
    """
    import threading as _th
    import time as _time
    from concurrent.futures import ThreadPoolExecutor

    jobs = 4
    gate = _th.Barrier(jobs, timeout=10)
    tally = _th.Lock()
    state = {"inside": 0, "most": 0, "calls": 0}

    def slow_receptor(pdb, out):
        with tally:
            state["inside"] += 1
            state["calls"] += 1
            state["most"] = max(state["most"], state["inside"])
        _time.sleep(0.05)          # the window a real converter leaves open
        Path(out).write_text("REC", encoding="utf-8")
        with tally:
            state["inside"] -= 1
        return Path(out)

    monkeypatch.setattr(dk, "_receptor_pdbqt", slow_receptor)
    cache = {}
    box = (0, 0, 0, 10, 10, 10)

    def one(smi):
        gate.wait()                # all of them miss the cache together
        return dk.dock_one_vina(vina_world["pdb"], smi, box,
                                receptor_pdbqt_cache=cache)

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        got = list(pool.map(one, ["CCO", "CCN", "CCC", "CCF"]))

    assert got == [-8.2] * jobs, "every compound still gets its score"
    assert state["most"] == 1, (
        "two converters wrote the same receptor file at once, so a thread "
        "could dock against a half-written one")
    assert state["calls"] == 1, "and the cache still saves the other three"


def test_dock_one_vina_fails_loudly_when_the_receptor_cannot_be_converted(
        vina_world, monkeypatch):
    monkeypatch.setattr(dk, "_receptor_pdbqt", lambda pdb, out: None)
    with pytest.raises(RuntimeError, match="could not convert the receptor"):
        dk.dock_one_vina(vina_world["pdb"], "CCO", (0, 0, 0, 1, 1, 1))


def test_dock_one_vina_is_none_for_a_ligand_that_cannot_be_prepared(
        vina_world):
    vina_world["ligand"] = None
    assert dk.dock_one_vina(vina_world["pdb"], "CCO",
                            (0, 0, 0, 1, 1, 1)) is None


def test_dock_one_vina_without_vina_anywhere_says_so(monkeypatch, tmp_path):
    monkeypatch.setattr(dk.shutil, "which", lambda *a, **k: None)
    monkeypatch.setattr(wsldock, "find_binary", lambda prefer=None: (None, None))
    with pytest.raises(RuntimeError, match="vina not found"):
        dk.dock_one_vina(tmp_path / "r.pdb", "CCO", (0, 0, 0, 1, 1, 1))


# --------------------------------------------------- ligand and receptor prep
def test_meeko_ligand_is_none_when_meeko_is_not_installed(monkeypatch):
    monkeypatch.setitem(sys.modules, "meeko", None)
    assert dk._meeko_ligand_pdbqt("CCO") is None


def test_meeko_ligand_makes_a_pdbqt_and_refuses_bad_smiles():
    pytest.importorskip("meeko")
    pytest.importorskip("rdkit")
    text = dk._meeko_ligand_pdbqt("CCO", seed=1)
    assert text and "ROOT" in text and "ATOM" in text
    assert dk._meeko_ligand_pdbqt("not a smiles ((") is None


def test_receptor_pdbqt_uses_the_python_tools_first(monkeypatch, tmp_path):
    out = tmp_path / "rec.pdbqt"

    def fake(cmd, **kw):
        assert cmd[0] == sys.executable
        out.write_text("REC", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(subprocess, "run", fake)
    assert dk._receptor_pdbqt(tmp_path / "rec.pdb", out) == out


def test_receptor_pdbqt_falls_back_to_open_babel(monkeypatch, tmp_path):
    out = tmp_path / "rec.pdbqt"
    calls = []

    def fake(cmd, **kw):
        calls.append(cmd[0])
        if cmd[0] == "/usr/bin/obabel":
            out.write_text("REC", encoding="utf-8")
            return subprocess.CompletedProcess(cmd, 0, "", "")
        return subprocess.CompletedProcess(cmd, 1, "", "no module")
    monkeypatch.setattr(subprocess, "run", fake)
    monkeypatch.setattr(dk.shutil, "which",
                        lambda n, *a, **k: "/usr/bin/obabel"
                        if n == "obabel" else None)
    assert dk._receptor_pdbqt(tmp_path / "rec.pdb", out) == out
    assert calls[-1] == "/usr/bin/obabel"


def test_receptor_pdbqt_is_none_when_no_converter_works(monkeypatch,
                                                        tmp_path):
    monkeypatch.setattr(
        subprocess, "run",
        lambda cmd, **kw: subprocess.CompletedProcess(cmd, 1, "", "no"))
    monkeypatch.setattr(dk.shutil, "which", lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "openbabel", None)
    assert dk._receptor_pdbqt(tmp_path / "rec.pdb",
                              tmp_path / "rec.pdbqt") is None


# ----------------------------------------------------- printing and records
def test_direction_cell_marks_a_weak_call_with_a_question_mark():
    cell = dk._direction_cell
    assert cell({}) == "-"
    assert cell({"direction": "undetermined"}) == "-"
    assert cell({"direction": "agonist", "direction_confidence": "high"}) \
        == "agonist"
    assert cell({"direction": "antagonist",
                 "direction_confidence": "moderate"}) == "antagonist"
    assert cell({"direction": "agonist", "direction_confidence": "low"}) \
        == "agonist?"
    assert cell({"direction": "agonist"}) == "agonist?"


def _lf_screen_result(with_cnn=True):
    hits = [
        {"drug": "Drug A", "score": -9.1, "direction": "agonist",
         "direction_confidence": "high"},
        {"drug": "Drug B", "score": -8.0, "direction": "antagonist"},
    ]
    if with_cnn:
        hits[0].update(cnn_score=0.81, cnn_affinity=6.2)
        hits[1].update(cnn_score=None, cnn_affinity=None)
    return {
        "pocket": {"pdb_id": "1ABC", "ligand_resname": "XYZ",
                   "state": "active"},
        "conformation_verified": True, "auc": 0.91, "ef1": 8.5,
        "gate": {"n_actives": 60, "n_actives_scored": 58, "n_decoys": 600,
                 "n_decoys_scored": 590,
                 "ef1_basis": {"actives_in_top_n": 5, "top_n": 7}},
        "n_actives_available": 100, "n_actives_used_for_gate": 60,
        "auc_by_metric": {"affinity": {"auc": 0.9, "ef1": 8.0},
                          "cnn_score": {"auc": 0.95, "ef1": None}},
        "direction": {"attempted": True, "model_auc": 0.8, "n_labelled": 30,
                      "source": "ChEMBL", "n_called": 2,
                      "n_out_of_domain": 1},
        "hits": hits, "caveats": ["read the gate first"],
        "conformation_caveat": "not a functional assay",
    }


def test_report_prints_the_gate_basis_the_hits_and_the_caveats(capsys):
    dk._report(_lf_screen_result())
    out = capsys.readouterr().out
    assert "Structure : 1ABC (ligand XYZ)" in out
    assert "Conformation: active, read from the structure" in out
    assert "AUC 0.91" in out and "top-1% enrichment 8.5" in out
    assert "58 of 60 known ligand(s) and 590 of 600 decoy(s) scored" in out
    assert "enrichment counts 5 known ligand(s) among the top 7" in out
    assert "measured on 60 of 100 known ligand(s), sampled at random" in out
    assert "What each score was worth" in out
    assert "CNNscore" in out and "Drug A" in out
    assert "agonist" in out and "antagonist?" in out
    assert "read the gate first" in out and "not a functional assay" in out


def test_report_without_cnn_columns_prints_no_blank_measurements(capsys):
    res = _lf_screen_result(with_cnn=False)
    res["conformation_verified"] = False
    dk._report(res)
    out = capsys.readouterr().out
    assert "CNNscore" not in out
    assert "NOT read from the structure" in out
    assert "Drug B" in out


def test_report_of_a_refusal_prints_the_reason_and_nothing_else(capsys):
    dk._report({"error": "the gate refused"})
    assert capsys.readouterr().out.strip() == "Not reported: the gate refused"


def test_print_check_exit_code_follows_readiness(monkeypatch, capsys):
    monkeypatch.setattr(dk, "check_setup", lambda: {
        "ok": True, "missing": [],
        "found": {"docking_binary": "/x/gnina", "engine": "gnina",
                  "rdkit": True}})
    assert dk._print_check() == 0
    assert "Ready" in capsys.readouterr().out
    monkeypatch.setattr(dk, "check_setup", lambda: {
        "ok": False, "missing": ["a docking binary"],
        "found": {"rdkit": False}})
    assert dk._print_check() == 1
    out = capsys.readouterr().out
    assert "NOT FOUND" in out and "a docking binary" in out
    assert "NOT INSTALLED" in out


def test_write_refusal_records_a_refusal_as_a_result(tmp_path, capsys):
    path = tmp_path / "sub" / "r.json"
    dk._write_refusal(path, "ESR1", "P03372", "no pocket")
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["usable"] is False and data["target"] == "ESR1"
    assert data["refused_before_screening"] is True
    assert data["error"] == "no pocket" and data["hits"] == []
    assert "Recorded the refusal" in capsys.readouterr().out
    # and a run that refused is loadable like any other result
    assert dk.load_results(path)["usable"] is False


def test_write_refusal_with_no_path_does_nothing(tmp_path):
    dk._write_refusal(None, "ESR1", None, "x")
    assert list(tmp_path.iterdir()) == []


def test_write_refusal_survives_an_unwritable_location(tmp_path, caplog):
    import logging
    blocker = tmp_path / "file"
    blocker.write_text("x", encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        dk._write_refusal(blocker / "r.json", "ESR1", None, "x")
    assert any("could not record the refusal" in r.message
               for r in caplog.records)


# ------------------------------------------------------------------- wsldock
@pytest.fixture
def fresh_wsl(monkeypatch):
    monkeypatch.setattr(wsldock, "_CACHE", {})
    return wsldock


def _lf_cp(rc=0, out=""):
    return subprocess.CompletedProcess([], rc, out, "")


def test_wsl_available_needs_wsl_exe_and_a_script_that_runs(fresh_wsl,
                                                           monkeypatch):
    monkeypatch.setattr(wsldock.shutil, "which", lambda n: None)
    assert wsldock.available() is False

    monkeypatch.setattr(wsldock, "_CACHE", {})
    monkeypatch.setattr(wsldock.shutil, "which", lambda n: "C:/wsl.exe")
    monkeypatch.setattr(wsldock, "run_script", lambda s, timeout=0: _lf_cp(0))
    assert wsldock.available() is True


def test_wsl_available_is_false_when_the_probe_raises_and_is_cached(
        fresh_wsl, monkeypatch):
    calls = []

    def boom(script, timeout=0):
        calls.append(1)
        raise OSError("no distro")
    monkeypatch.setattr(wsldock.shutil, "which", lambda n: "C:/wsl.exe")
    monkeypatch.setattr(wsldock, "run_script", boom)
    assert wsldock.available() is False
    assert wsldock.available() is False
    assert len(calls) == 1, "the answer is remembered"


def test_find_binary_returns_the_first_program_found_and_remembers_it(
        fresh_wsl, monkeypatch):
    scripts = []

    def fake(script, timeout=0):
        scripts.append(script)
        return _lf_cp(0, "/home/u/bin/gnina" + NL) \
            if "gnina" in script else _lf_cp(1)
    monkeypatch.setattr(wsldock, "available", lambda: True)
    monkeypatch.setattr(wsldock, "run_script", fake)
    assert wsldock.find_binary() == ("gnina", "/home/u/bin/gnina")
    assert wsldock.find_binary() == ("gnina", "/home/u/bin/gnina")
    assert len(scripts) == 1


def test_find_binary_honours_a_preference_and_falls_through(fresh_wsl,
                                                           monkeypatch):
    asked = []

    def fake(script, timeout=0):
        asked.append(script.split()[2])          # "command -v NAME ..."
        return _lf_cp(0, "/usr/bin/smina" + NL) if "smina" in script \
            else _lf_cp(1)
    monkeypatch.setattr(wsldock, "available", lambda: True)
    monkeypatch.setattr(wsldock, "run_script", fake)
    assert wsldock.find_binary(prefer="vina") == ("smina", "/usr/bin/smina")
    assert asked[0] == "vina", "the preferred program is asked for first"


def test_find_binary_is_none_without_wsl_or_without_programs(fresh_wsl,
                                                            monkeypatch):
    monkeypatch.setattr(wsldock, "available", lambda: False)
    assert wsldock.find_binary() == (None, None)
    monkeypatch.setattr(wsldock, "_CACHE", {})
    monkeypatch.setattr(wsldock, "available", lambda: True)
    monkeypatch.setattr(wsldock, "run_script", lambda s, timeout=0: _lf_cp(1))
    assert wsldock.find_binary() == (None, None)


def test_library_path_finds_the_cuda_directory_once(fresh_wsl, monkeypatch):
    seen = []
    monkeypatch.setattr(wsldock, "available", lambda: True)
    monkeypatch.setattr(wsldock, "_first_existing",
                        lambda c, t: seen.append(t) or "/home/u/cuda-libs")
    assert wsldock.library_path() == "/home/u/cuda-libs"
    assert wsldock.library_path() == "/home/u/cuda-libs"
    assert len(seen) == 1 and "libcudnn" in seen[0]


def test_library_path_is_none_without_wsl(fresh_wsl, monkeypatch):
    monkeypatch.setattr(wsldock, "available", lambda: False)
    assert wsldock.library_path() is None


def test_dock_script_exports_the_cuda_libraries_only_when_there_are_some(
        monkeypatch):
    monkeypatch.setattr(wsldock, "library_path", lambda: "/home/u/cuda-libs")
    script = wsldock.dock_script("/home/u/bin/gnina", ["-r", "a b", 3])
    assert "/home/u/cuda-libs:" in script
    assert '"$LD_LIBRARY_PATH"' in script, "the variable must still expand"
    monkeypatch.setattr(wsldock, "library_path", lambda: None)
    assert wsldock.dock_script("/x/vina", ["-r"]) == "/x/vina -r" + NL


def test_dock_script_keeps_a_path_with_a_space_as_one_argument(monkeypatch):
    """A Windows username with a space split one argument into two.

    gnina then exited non-zero, `_run_dock` recorded the compound as
    unscored, and a whole screen came back empty rather than broken - for
    a username. The script must parse back into exactly what was asked.
    """
    import shlex

    monkeypatch.setattr(wsldock, "library_path", lambda: None)
    args = ["-r", "/mnt/c/Users/Ada Lovelace/rec.pdb", "-l", "lig.sdf",
            "--seed", 0]
    script = wsldock.dock_script("/home/u/bin/gnina", args)
    assert shlex.split(script) == ["/home/u/bin/gnina"] + [
        str(a) for a in args]


# ------------------------------------------------ the library limit says so
def test_a_library_limit_that_drops_drugs_says_it_took_the_start_of_the_list(
        caplog):
    """A limit of 600 over an alphabetical table docked A to C and never
    reached maraviroc, the approved CCR5 drug and the pocket's own ligand."""
    import logging
    rows = [f"drug{i:04d}" for i in range(1000)]
    with caplog.at_level(logging.WARNING):
        got = dk._limit_library(rows, 600)
    assert got == rows[:600]
    msg = " ".join(r.getMessage() for r in caplog.records)
    assert "first 600 of 1000" in msg
    assert "not a random sample" in msg and "--library-limit" in msg


def test_a_limit_that_drops_nothing_is_silent_and_a_copy(caplog):
    import logging
    rows = ["a", "b", "c"]
    with caplog.at_level(logging.WARNING):
        got = dk._limit_library(rows, 100000)
        none = dk._limit_library(rows, None)
    assert got == rows and none == rows
    assert got is not rows, "a copy, not the caller's list"
    assert not caplog.records
