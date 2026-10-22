"""Docking a drug as stored (the full structure, salt included) and without
its counter-ions.

About a third of the approved library is written as a salt. By default the
structure is docked exactly as stored, so no earlier result changes. With
`strip_counter_ions` each drug is docked without its salt. These tests run
both and check that the two differ only where they should.
"""

from __future__ import annotations

import csv
import json
import random
from pathlib import Path

import pytest

from generoutes import docking as dk

NL = chr(10)
TAB = chr(9)

ACAMPROSATE_CA = ("CC(=O)NCCCS(=O)(=O)[O-].CC(=O)NCCCS(=O)(=O)[O-]."
                  "[Ca+2]")
ACAMPROSATE = "CC(=O)NCCCS(=O)(=O)[O-]"
ABACAVIR = "Nc1nc(NC2CC2)c2ncn([C@H]3C=C[C@@H](CO)C3)c2n1"
ABACAVIR_SULFATE = ".".join([ABACAVIR, ABACAVIR, "O=S(=O)(O)O"])
ACEFYLLINE = "Cn1c(=O)c2c(ncn2CC(=O)O)n(C)c1=O"
ACEFYLLINE_PIPERAZINE = ".".join(["C1CNCCN1", ACEFYLLINE, ACEFYLLINE])


# ------------------------------------------------------ split_counter_ions
@pytest.mark.parametrize("stored,parent,removed,ambiguous", [
    ("CCO", "CCO", [], False),
    ("CCN.Cl", "CCN", ["Cl"], False),
    (ACAMPROSATE_CA, ACAMPROSATE, ["[Ca+2]"], False),
    (ABACAVIR_SULFATE, ABACAVIR, ["O=S(=O)(O)O"], False),
    ("CC(=O)O.CCN", "CCN", ["CC(=O)O"], False),                # acetate
    ("OCC(O)C(O)C(O)C(O)C(=O)O.CCN", "CCN",                    # gluconate
     ["O=C(O)C(O)C(O)C(O)C(O)CO"], False),     # RDKit's canonical form
    ("CCN.O", "CCN", ["O"], False),                            # a hydrate
])
def test_a_drug_is_separated_from_its_salt(stored, parent, removed,
                                           ambiguous):
    got = dk.split_counter_ions(stored)
    assert got["smiles"] == parent
    assert got["removed"] == removed
    assert got["ambiguous"] is ambiguous


def test_two_copies_of_the_drug_become_one():
    got = dk.split_counter_ions(ABACAVIR_SULFATE)
    assert got["n_fragments"] == 3
    assert got["smiles"].count(".") == 0, "a dimer is not what binds"


def test_a_second_organic_component_is_kept_and_flagged_not_guessed():
    got = dk.split_counter_ions(ACEFYLLINE_PIPERAZINE)
    assert got["ambiguous"] is True
    assert got["removed"] == []
    parts = set(got["smiles"].split("."))
    assert parts == {"C1CNCCN1", ACEFYLLINE}, "both kept, one copy each"
    pair = dk.split_counter_ions("CC(=O)N[C@@H](CS)C(=O)O."
                                 "NCCCC[C@H](N)C(=O)O")
    assert pair["ambiguous"] is True and len(pair["smiles"].split(".")) == 2


def test_a_counter_ion_that_is_also_a_drug_does_not_lose_the_wrong_half():
    """Piperazine citrate is piperazine; sodium citrate is citrate. Neither
    may come back empty, and piperazine must not be thrown away as if it
    were the salt."""
    got = dk.split_counter_ions("C1CNCCN1.O=C(O)CC(O)(CC(=O)O)C(=O)O")
    assert got["smiles"] == "C1CNCCN1"
    assert got["removed"] == ["O=C(O)CC(O)(CC(=O)O)C(=O)O"]
    sodium = dk.split_counter_ions(
        "[Na+].[Na+].[Na+].O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]")
    assert sodium["smiles"] == "O=C([O-])CC(O)(CC(=O)[O-])C(=O)[O-]"
    assert sodium["removed"] == ["[Na+]"]


def test_when_everything_is_a_counter_ion_something_is_still_returned():
    got = dk.split_counter_ions("[Na+].[Cl-]")
    assert got["smiles"] in ("[Na+]", "[Cl-]")
    assert len(got["removed"]) == 1


def test_charges_are_left_as_stored():
    got = dk.split_counter_ions(ACAMPROSATE_CA)
    assert "[O-]" in got["smiles"], "protonation is not changed"


def test_an_unparseable_structure_raises():
    with pytest.raises(ValueError, match="could not parse"):
        dk.split_counter_ions("not a smiles ((")


def test_the_real_library_is_mostly_unchanged_and_few_entries_are_ambiguous():
    path = Path(__file__).resolve().parents[1] / "data" / "approved_library.tsv"
    if not path.exists():
        pytest.skip("data/approved_library.tsv is not in this checkout")
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter=TAB))
    changed = ambiguous = 0
    for r in rows:
        got = dk.split_counter_ions(r["smiles"])
        assert got["smiles"], r["name"]
        changed += bool(got["removed"])
        ambiguous += got["ambiguous"]
    assert changed > 0.25 * len(rows), "about a third are stored as salts"
    assert ambiguous < 0.03 * len(rows), "and almost all can be separated"


# ------------------------------------------------ the full structure vs stripped
def _lf_record():
    seen = []

    def scorer(smiles):
        seen.append(smiles)
        # A salt docks worse than the drug alone, which is the artefact the
        # option exists to remove.
        return -4.0 if "." in smiles else -9.0
    return seen, scorer


def test_by_default_the_structure_is_docked_exactly_as_stored(tmp_path):
    seen, scorer = _lf_record()
    out = dk._dock_smiles([ACAMPROSATE_CA, "CCO"], tmp_path / "r.pdb", None,
                          exhaustiveness=1, seed=0, score_fn=scorer)
    assert seen == [ACAMPROSATE_CA, "CCO"]
    assert out == [-4.0, -9.0]


def test_with_the_option_the_drug_is_docked_without_its_salt(tmp_path):
    seen, scorer = _lf_record()
    out = dk._dock_smiles([ACAMPROSATE_CA, "CCO"], tmp_path / "r.pdb", None,
                          exhaustiveness=1, seed=0, score_fn=scorer,
                          strip_counter_ions=True)
    assert seen == [ACAMPROSATE, "CCO"]
    assert out == [-9.0, -9.0]


def test_an_unparseable_entry_is_passed_through_and_counted_as_a_failure(
        tmp_path):
    seen, scorer = _lf_record()
    dk._dock_smiles(["not a smiles (("], tmp_path / "r.pdb", None,
                    exhaustiveness=1, seed=0, score_fn=scorer,
                    strip_counter_ions=True)
    assert seen == ["not a smiles (("], "left as it was; the dock decides"


def test_the_3d_structure_of_a_salt_has_the_extra_pieces_unless_stripped(
        tmp_path):
    """What the docking program is actually handed."""
    pytest.importorskip("rdkit")
    from rdkit import Chem

    def fragments(smiles, name):
        sdf = tmp_path / name
        dk.smiles_to_sdf(smiles, sdf, seed=1)
        mol = Chem.MolFromMolFile(str(sdf), removeHs=False)
        return len(Chem.GetMolFrags(mol))

    assert fragments(ACAMPROSATE_CA, "full.sdf") == 3
    parent = dk.split_counter_ions(ACAMPROSATE_CA)["smiles"]
    assert fragments(parent, "stripped.sdf") == 1


def test_the_validation_gate_docks_known_ligands_both_ways(tmp_path):
    rec = tmp_path / "r.pdb"
    rec.write_text("ATOM" + NL, encoding="utf-8")
    actives = [{"smiles": f"CCN{'C' * i}.Cl"} for i in range(6)]
    decoys = [{"smiles": f"CCO{'C' * i}"} for i in range(30)]
    got = {}
    for flag in (False, True):
        seen = []

        def scorer(smiles, seen=seen):
            seen.append(smiles)
            return -9.0 if smiles.startswith("CCN") else -5.0
        out = dk.validate_pocket(rec, actives, decoys, score_fn=scorer,
                                 box=(0, 0, 0, 20, 20, 20),
                                 strip_counter_ions=flag)
        got[flag] = (out, seen)
    assert got[False][0]["usable"] and got[True][0]["usable"]
    assert any(s.endswith(".Cl") for s in got[False][1])
    assert not any("." in s for s in got[True][1]), "no salt reached the dock"


def _lf_library():
    return [{"name": "Plain", "smiles": "CCO"},
            {"name": "Acamprosate calcium", "smiles": ACAMPROSATE_CA},
            {"name": "Acefylline piperazine", "smiles": ACEFYLLINE_PIPERAZINE}]


def _lf_screen(tmp_path, strip):
    rec = tmp_path / "r.pdb"
    rec.write_text("ATOM" + NL, encoding="utf-8")
    rng = random.Random(2)
    actives = [{"smiles": f"CCN{'C' * i}"} for i in range(8)]
    decoys = [{"smiles": f"CCO{'C' * i}"} for i in range(40)]

    def scorer(smiles):
        if smiles.startswith("CCN"):
            return -10.0 - rng.random()
        return -4.0 if "." in smiles else -8.0
    return dk.screen_by_docking(
        rec, _lf_library(), target="T", actives=actives, decoys=decoys,
        score_fn=scorer, box=(0, 0, 0, 20, 20, 20),
        strip_counter_ions=strip)


def test_a_full_structure_screen_records_nothing_about_ligand_preparation(
        tmp_path):
    run = _lf_screen(tmp_path, strip=False)
    assert run.usable and run.ligand_prep is None
    by = {h["drug"]: h for h in run.hits}
    assert by["Acamprosate calcium"]["score"] == -4.0
    assert by["Acefylline piperazine"]["score"] == -4.0


def test_a_stripped_screen_scores_the_drug_and_says_what_it_changed(tmp_path):
    run = _lf_screen(tmp_path, strip=True)
    assert run.usable
    assert run.ligand_prep == {"strip_counter_ions": True, "n_library": 3,
                               "n_changed": 1, "n_ambiguous": 1}
    by = {h["drug"]: h for h in run.hits}
    assert by["Acamprosate calcium"]["score"] == -8.0, "docked as the drug"
    assert by["Acefylline piperazine"]["score"] == -4.0, (
        "a mixed entry is kept whole, not guessed at")
    assert by["Acamprosate calcium"]["smiles"] == ACAMPROSATE_CA, (
        "the result still names the structure that was stored")
    assert run.to_dict()["ligand_prep"]["n_changed"] == 1


def _lf_bundle(tmp_path):
    (tmp_path / "rec.pdb").write_text("ATOM" + NL, encoding="utf-8")
    manifest = {
        "receptor": "rec.pdb", "box": [0, 0, 0, 20, 20, 20],
        "symbol": "T", "pdb_id": "1ABC",
        "actives": [{"smiles": f"CCN{'C' * i}"} for i in range(8)],
        "library": _lf_library() + [{"name": f"Filler{i}",
                                     "smiles": f"CCO{'C' * i}"}
                                    for i in range(60)],
    }
    (tmp_path / dk.BUNDLE_FILE).write_text(json.dumps(manifest),
                                           encoding="utf-8")
    return tmp_path


def test_a_bundle_run_passes_the_option_through_to_the_screen(tmp_path):
    bundle = _lf_bundle(tmp_path)
    rng = random.Random(3)

    def scorer(smiles):
        if smiles.startswith("CCN"):
            return -10.0 - rng.random()
        return -4.0 if "." in smiles else -8.0
    full = dk.dock_from_bundle(bundle, score_fn=scorer, top_n=100)
    stripped = dk.dock_from_bundle(bundle, score_fn=scorer, top_n=100,
                                   strip_counter_ions=True)
    assert full["usable"] and stripped["usable"]
    assert full.get("ligand_prep") is None
    assert stripped["ligand_prep"]["n_changed"] == 1
    score = {h["drug"]: h["score"] for h in full["hits"]}
    score2 = {h["drug"]: h["score"] for h in stripped["hits"]}
    assert score["Acamprosate calcium"] == -4.0
    assert score2["Acamprosate calcium"] == -8.0
    assert score["Plain"] == score2["Plain"], "a single-piece drug is unmoved"


def test_the_command_line_flag_reaches_the_bundle_run(monkeypatch, capsys):
    seen = {}
    monkeypatch.setattr(dk, "check_setup",
                        lambda: {"ok": True, "missing": [], "found": {}})

    def fake(bundle, **kw):
        seen.update(kw)
        return {"error": None, "hits": []}
    monkeypatch.setattr(dk, "dock_from_bundle", fake)
    monkeypatch.setattr(dk, "_report", lambda res, top: None)
    assert dk.main(["--bundle", "somewhere"]) == 0
    assert seen["strip_counter_ions"] is False, "off unless asked for"
    assert dk.main(["--bundle", "somewhere", "--strip-counter-ions"]) == 0
    assert seen["strip_counter_ions"] is True


def test_the_report_says_how_ligands_were_prepared(capsys):
    res = {"pocket": {"pdb_id": "1ABC", "ligand_resname": "X"},
           "auc": 0.9, "ef1": 5.0, "hits": [],
           "ligand_prep": {"strip_counter_ions": True, "n_library": 600,
                           "n_changed": 190, "n_ambiguous": 7}}
    dk._report(res)
    out = capsys.readouterr().out
    assert "counter-ions removed from 190 of 600 library entries" in out
    assert "7 with a second organic component were kept whole" in out
    res["ligand_prep"] = None
    dk._report(res)
    assert "counter-ions" not in capsys.readouterr().out
