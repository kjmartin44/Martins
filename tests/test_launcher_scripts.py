"""Direct tests for the helper scripts the launcher runs between stages and
for the result records every stage writes into results.json.
"""

from __future__ import annotations

import dataclasses
import importlib
import importlib.util
import json
import typing
from pathlib import Path

import pytest

NL = chr(10)


def _lf_script(name):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"_lf_script_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------ merge_docking
def _lf_run(target, usable=True):
    return {"target": target, "usable": usable, "hits": []}


def test_merge_docking_expands_a_glob_the_shell_did_not(tmp_path):
    md = _lf_script("merge_docking")
    for n in ("b.json", "a.json", "c.txt"):
        (tmp_path / n).write_text("{}", encoding="utf-8")
    got = md.expand([str(tmp_path / "*.json"), str(tmp_path / "c.txt")])
    assert [p.name for p in got] == ["a.json", "b.json", "c.txt"]
    plain = md.expand([str(tmp_path / "no_such.json")])
    assert [p.name for p in plain] == ["no_such.json"], "kept, skipped later"


def test_merge_docking_flattens_lists_and_keeps_refusals(tmp_path, capsys):
    md = _lf_script("merge_docking")
    (tmp_path / "one.json").write_text(
        json.dumps(_lf_run("ESR1")), encoding="utf-8")
    (tmp_path / "two.json").write_text(
        json.dumps([_lf_run("AR", usable=False), {"not": "a run"}]),
        encoding="utf-8")
    (tmp_path / "bad.json").write_text("{not json", encoding="utf-8")
    runs = md.merge([tmp_path / "one.json", tmp_path / "two.json",
                     tmp_path / "bad.json", tmp_path / "absent.json"])
    assert [r["target"] for r in runs] == ["ESR1", "AR"]
    assert runs[1]["usable"] is False, "a refusal is a result and is kept"
    err = capsys.readouterr().err
    assert "skipping bad.json" in err
    assert "skipping a record in two.json" in err
    assert "absent.json" not in err, "a path that is not a file is skipped quietly"


def test_merge_docking_main_writes_the_combined_file(tmp_path, capsys):
    md = _lf_script("merge_docking")
    (tmp_path / "a.json").write_text(json.dumps(_lf_run("ESR1")),
                                     encoding="utf-8")
    (tmp_path / "b.json").write_text(json.dumps(_lf_run("AR", False)),
                                     encoding="utf-8")
    out = tmp_path / "all" / "merged.out"
    out.parent.mkdir()
    rc = md.main([str(out), str(tmp_path / "a.json"), str(tmp_path / "b.json")])
    assert rc == 0
    assert [r["target"] for r in json.loads(out.read_text(encoding="utf-8"))] \
        == ["ESR1", "AR"]
    text = capsys.readouterr().out
    assert "2 run(s)" in text and "1 validated, 1 refused" in text


def test_merge_docking_main_fails_when_there_is_nothing_to_merge(tmp_path,
                                                                 capsys):
    md = _lf_script("merge_docking")
    out = tmp_path / "o.json"
    assert md.main([str(out), str(tmp_path / "none.json")]) == 1
    assert not out.exists()
    assert "no docking runs to merge" in capsys.readouterr().err


# -------------------------------------------------------------- top_factors
def _lf_results(tmp_path, tfs, name="results.json"):
    p = tmp_path / name
    p.write_text(json.dumps({"results": [{"tfs": tfs}]}), encoding="utf-8")
    return p


_TFS = [
    {"symbol": "NR3C1", "role": "activator", "chip_peaks": 4, "score": 9.0},
    {"symbol": "SP1", "role": "unknown", "chip_peaks": 9, "score": 12.0},
    {"symbol": "REST", "role": "repressor", "chip_peaks": 2, "score": 7.0},
    {"symbol": "GENE1", "role": "activator", "chip_peaks": 5, "score": 99.0},
    {"symbol": "MYC", "role": "activator", "chip_peaks": 0, "score": 50.0},
]


def test_top_factors_main_prints_symbol_and_role_best_score_first(tmp_path,
                                                                  capsys):
    tf = _lf_script("top_factors")
    path = _lf_results(tmp_path, _TFS)
    assert tf.main([str(path), "--top", "5"]) == 0
    # SP1 has no resolved role and MYC has no ChIP peak: both left out.
    assert capsys.readouterr().out.splitlines() == [
        "GENE1 activator", "NR3C1 activator", "REST repressor"]


def test_top_factors_never_names_the_gene_being_studied(tmp_path, capsys):
    tf = _lf_script("top_factors")
    p = tmp_path / "r.json"
    p.write_text(json.dumps({"query": {"gene": "gene1"},
                             "results": [{"tfs": _TFS}]}), encoding="utf-8")
    assert tf.main([str(p), "--symbols-only", "--top", "5"]) == 0
    assert "GENE1" not in capsys.readouterr().out


def test_top_factors_main_options(tmp_path, capsys):
    tf = _lf_script("top_factors")
    path = _lf_results(tmp_path, _TFS)
    assert tf.main([str(path), "--symbols-only", "--top", "1"]) == 0
    assert capsys.readouterr().out.splitlines() == ["GENE1"]

    assert tf.main([str(path), "--any-role", "--min-peaks", "0",
                    "--top", "2", "--symbols-only"]) == 0
    assert capsys.readouterr().out.splitlines() == ["GENE1", "MYC"]

    assert tf.main([str(path), "--verbose", "--top", "1"]) == 0
    captured = capsys.readouterr()
    assert "role=activator" in captured.err and "GENE1 activator" in captured.out


def test_top_factors_main_writes_a_plain_file_with_no_byte_order_mark(
        tmp_path, capsys):
    tf = _lf_script("top_factors")
    path = _lf_results(tmp_path, _TFS)
    out = tmp_path / "factors.txt"
    assert tf.main([str(path), "--out", str(out), "--top", "2"]) == 0
    raw = out.read_bytes()
    assert not raw.startswith(bytes([0xEF, 0xBB, 0xBF])), "no BOM"
    assert raw.endswith(b"\n"), "newline-terminated"
    assert raw.decode("utf-8").splitlines() == [
        "GENE1 activator", "NR3C1 activator"]
    assert capsys.readouterr().out == "", "nothing goes to stdout"


def test_top_factors_main_fails_clearly(tmp_path, capsys):
    tf = _lf_script("top_factors")
    assert tf.main([str(tmp_path / "missing.json")]) == 1
    assert "no such file" in capsys.readouterr().err

    bad = tmp_path / "bad.json"
    bad.write_text("{oops", encoding="utf-8")
    assert tf.main([str(bad)]) == 1
    assert "could not read" in capsys.readouterr().err

    none = _lf_results(tmp_path, [{"symbol": "SP1", "role": "unknown",
                                   "chip_peaks": 3, "score": 1}], "n.json")
    assert tf.main([str(none)]) == 1
    assert "no factor met the bar" in capsys.readouterr().err


# ------------------------------------------------- result records: to_dict
def _lf_value(hint, name):
    """A plausible value for a required dataclass field, from its type."""
    text = str(hint)
    if "bool" in text:
        return True
    if "int" in text:
        return 1
    if "float" in text:
        return 0.5
    if "dict" in text:
        return {}
    if "list" in text or "tuple" in text:
        return []
    return name


def _lf_build(cls):
    hints = typing.get_type_hints(cls, include_extras=False) \
        if False else {f.name: f.type for f in dataclasses.fields(cls)}
    kwargs = {}
    for f in dataclasses.fields(cls):
        if f.default is dataclasses.MISSING \
                and f.default_factory is dataclasses.MISSING:
            kwargs[f.name] = _lf_value(hints[f.name], f.name)
    return cls(**kwargs)


_RECORDS = [
    ("docking", "DockingResult"), ("docking", "PocketAssessment"),
    ("docking", "PoseScore"), ("eqtl", "TransHit"), ("eqtl", "TransSupport"),
    ("agonism", "AgonismPrediction"), ("analogs", "Analog"),
    ("signatures", "SignatureHit"), ("abundance", "AbundanceHit"),
    ("cellcycle", "PhaseProfile"), ("delivery", "DeliveryAssessment"),
    ("liabilities", "Liability"), ("regulon", "RegulonEvidence"),
]
# PanelGene and Prediction write a deliberate summary rather than every
# field; each has its own test below.


@pytest.mark.parametrize("module,name", _RECORDS)
def test_every_result_record_serialises_to_json_with_all_its_fields(module,
                                                                    name):
    cls = getattr(importlib.import_module(f"promoter_tf_drugs.{module}"), name)
    record = _lf_build(cls)
    d = record.to_dict()
    assert isinstance(d, dict)
    json.dumps(d, default=str)
    field_names = {f.name for f in dataclasses.fields(cls)}
    missing = field_names - set(d)
    assert not missing, f"{name}.to_dict drops {sorted(missing)}"


def test_panel_gene_to_dict_carries_the_verdict_and_the_reason_for_exclusion():
    from promoter_tf_drugs.panel import PanelGene
    g = PanelGene("CCR5", result={"evidence_verdict": {
        "level": "degraded", "reasons": ["no ChIP"]}},
        excluded="not interpretable")
    assert g.to_dict() == {"gene": "CCR5", "level": "degraded",
                           "error": None, "excluded": "not interpretable",
                           "reasons": ["no ChIP"]}
    bare = PanelGene("X", error="Ensembl died").to_dict()
    assert bare["level"] == "unknown" and bare["reasons"] == []
    assert bare["error"] == "Ensembl died"


def test_prediction_to_dict_writes_the_fields_the_report_reads():
    """It is a summary: `mechanism` and `p_range` are on the record and are
    not written (the outcome block carries p_low and p_high). That is a
    choice about saved output, so it is pinned here, not changed silently."""
    from promoter_tf_drugs.predictions import Prediction
    p = Prediction("Dexamethasone", target="NR3C1", category="agonist",
                   mechanism="m", tier=2, grade="A", score=61.2,
                   p_move=0.4, p_range=(0.2, 0.6), limiting_step="delivery")
    d = p.to_dict()
    assert d["drug"] == "Dexamethasone" and d["tier"] == 2
    assert d["p_move"] == 0.4 and d["limiting_step"] == "delivery"
    assert set(d) == {"drug", "target", "bbb", "what_it_is", "treats",
                      "category", "tier", "grade", "score", "p_move",
                      "limiting_step", "why", "caveats", "falsifier"}
    json.dumps(d)


def test_to_dict_returns_a_copy_that_does_not_alias_the_record():
    from promoter_tf_drugs.docking import DockingResult
    r = DockingResult("Drug", "ESR1", notes=["a"])
    d = r.to_dict()
    d["notes"].append("b")
    assert r.notes == ["a"], "editing the dict must not edit the record"


# ------------------------------------------------------ launcher: Orange Book
def test_the_launcher_offers_the_orange_book_only_when_it_is_missing():
    """The prompt follows the gnomAD one: skipped when either file exists,
    default yes, and a download that is not the Orange Book is deleted."""
    bs = chr(92)
    src = (Path(__file__).resolve().parents[1] / "scripts"
           / "run_everything.bat").read_text(encoding="utf-8")
    zip_path, txt_path = f"data{bs}EOBZIP.zip", f"data{bs}products.txt"
    assert f'if exist {zip_path} set "HAVEOB=1"' in src
    assert f'if exist {txt_path} set "HAVEOB=1"' in src
    assert "if not defined HAVEOB (" in src
    assert 'set /p GETOB="Download it? Y/n: "' in src
    assert 'if /i not "!GETOB!"=="n" (' in src, "yes unless told otherwise"
    assert "https://www.fda.gov/media/76860/download?attachment" in src
    # the check that what arrived is the Orange Book, and the clean-up
    assert "endswith('products.txt')" in src
    assert f"del {zip_path}" in src
    # it sits with the other data downloads, before the inputs are asked
    assert src.index("GETOB") < src.index('set "GENE=%~1"')
    assert src.index("GETGNOMAD") < src.index("GETOB")


def test_the_orange_book_zip_the_launcher_downloads_is_one_the_package_reads(
        tmp_path, monkeypatch):
    """Same layout as the FDA file: three tables in one zip."""
    import zipfile

    from promoter_tf_drugs import orangebook
    # The module caches what it loaded; start empty and put it back after.
    monkeypatch.setattr(orangebook, "_PRODUCTS", {})
    monkeypatch.setattr(orangebook, "_LOADED", None)
    header = ("Ingredient~DF;Route~Trade_Name~Applicant~Strength~Appl_Type~"
              "Appl_No~Product_No~TE_Code~Approval_Date~RLD~RS~Type~"
              "Applicant_Full_Name")
    row = ("TETRABENAZINE~TABLET;ORAL~XENAZINE~LUNDBECK~12.5MG~N~021894~001~"
           "AB~Aug 15, 2008~Yes~No~RX~Lundbeck")
    z = tmp_path / "EOBZIP.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("products.txt", NL.join([header, row]) + NL)
        zf.writestr("patent.txt", "x" + NL)
    got =orangebook.lookup("tetrabenazine", path=str(z))
    assert got is not None and got["systemic_routes"] == ["ORAL"]
