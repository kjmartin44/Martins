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
                      "p_engages", "p_activates", "knockout",
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


def test_the_launcher_runs_its_expensive_processes_by_default():
    """The whole approved library, the steadier search setting, more factors,
    and the activity test with an empirical null. Each costs time and none
    of them is optional in the sense of being skipped unless asked for."""
    bs = chr(92)
    src = (Path(__file__).resolve().parents[1] / "scripts"
           / "run_everything.bat").read_text(encoding="utf-8")
    assert 'set "LIBLIMIT=100000"' in src and 'set "EXH=16"' in src
    # both docking commands and the exploratory dock carry them
    assert src.count("--library-limit !LIBLIMIT! --exhaustiveness !EXH!") == 3
    assert 'set "NTF=5"' in src and 'if "!NTF!"=="" set "NTF=5"' in src
    assert "--activity-permutations 2000" in src
    # the permutation flag only exists together with the activity test
    line = [ln for ln in src.splitlines()
            if "--activity-permutations" in ln and "set " in ln][0]
    assert "--tf-activity" in line and "--activity-tissue" in line
    # what the launcher says it will do matches what it does
    assert "the whole" in src and "exhaustiveness !EXH!" in src
    assert "600 compounds" not in src


def test_the_launcher_ends_with_a_combined_list_that_uses_only_real_flags():
    """Stage 5 fuses the docking results with the panel's survivors. It must
    run after the panel, take the per-target docking files rather than the
    merged one (which holds the same hits twice), and pass flags the command
    accepts."""
    import subprocess
    import sys
    bs = chr(92)
    src = (Path(__file__).resolve().parents[1] / "scripts"
           / "run_everything.bat").read_text(encoding="utf-8")
    start = src.index("-m promoter_tf_drugs.combined")
    optional = src.index("REM Optional: the dose check")
    block = src[start:src.rindex("REM ---", 0, optional)]
    assert src.index("promoter_tf_drugs.panel") < start
    assert f'for %%F in ("!DOCKDIR!{bs}*.json")' in src
    assert "docking_all.json" not in block
    assert f'--results "!OUT!{bs}results.json"' in block
    help_text = subprocess.run(
        [sys.executable, "-m", "promoter_tf_drugs.combined", "--help"],
        capture_output=True, text=True, cwd=str(
            Path(__file__).resolve().parents[1])).stdout
    flags = {t for t in block.replace("^", " ").split() if t.startswith("--")}
    flags |= {"--docking", "--panel"}
    assert flags, "no flags found"
    for f in flags:
        assert f in help_text, f"{f} is not a flag of the combined command"
    assert f"{bs}panel.json" in src and "!OUT!_combined" in src


# ------------------------------------------------- stage count and extra steps
def _launcher():
    return (Path(__file__).resolve().parents[1] / "scripts"
            / "run_everything.bat").read_text(encoding="utf-8")


def _command_flags(src, module, marker):
    """Flags the launcher passes to `module`, read from the lines that follow
    its `-m` call until the next blank or non-continued line."""
    import re
    start = src.index(f"-m promoter_tf_drugs.{module}")
    lines = src[start:].splitlines()
    block = []
    for ln in lines:
        block.append(ln)
        if not ln.rstrip().endswith("^"):
            break
    return set(re.findall(r"--[a-z][a-z0-9-]*", " ".join(block)))


def test_the_launcher_counts_five_stages_and_labels_each():
    src = _launcher()
    for n in range(1, 6):
        assert f"=== {n} of 5:" in src, n
    assert "of 4:" not in src


def test_the_launcher_says_which_data_gated_routes_it_will_skip():
    """A route whose data is not in data/ used to vanish without a word. The
    measured phenocopy route needs both matrices, so one alone must not turn
    the flag on."""
    src = _launcher()
    assert "if defined HAVEPRISM if defined HAVEDEPMAP set \"SENS=" in src
    for needle in ("no PRISM matrix", "no DepMap CRISPRGeneEffect.csv",
                   "no GWAS Catalog table", "no tissue given for factor activity"):
        assert needle in src, needle
    # every file name the launcher tests for is one the package also reads
    from promoter_tf_drugs import depmap, localannot, sensitivity
    for name in ("Repurposing_Public_24Q2_LFC_COLLAPSED.csv", "prism_lfc.csv"):
        assert f"data/{name}" in sensitivity.DEFAULT_PRISM_PATHS
        assert name in src
    assert "data/CRISPRGeneEffect.csv" in depmap.DEFAULT_EFFECT_PATHS
    assert "CRISPRGeneEffect.csv" in src
    assert "data/gwas_catalog_associations.tsv" in localannot.DEFAULT_GWAS_PATHS
    assert "gwas_catalog_associations.tsv" in src


def test_the_optional_steps_run_after_the_combined_list_with_real_flags():
    import subprocess
    import sys
    root = Path(__file__).resolve().parents[1]
    src = _launcher()
    assert (src.index("-m promoter_tf_drugs.combined")
            < src.index("-m promoter_tf_drugs.exposure")
            < src.index("echo === Done ===")
            < src.rindex("-m promoter_tf_drugs.exploratory")
            or src.index("-m promoter_tf_drugs.exploratory")
            < src.index("echo === Done ==="))
    for module in ("exposure", "exploratory"):
        help_text = subprocess.run(
            [sys.executable, "-m", f"promoter_tf_drugs.{module}", "--help"],
            capture_output=True, text=True, cwd=str(root)).stdout
        flags = _command_flags(src, module, module)
        assert flags, module
        for f in flags:
            assert f in help_text, f"{f} is not a flag of {module}"
    # nothing is asked for: the exploratory dock finds its own template
    block = src[src.index("-m promoter_tf_drugs.exploratory"):]
    block = block[:block.index("echo === Done ===")]
    for flag in ("--template-pdb", "--template-chain", "--template-ligand"):
        assert flag not in block
    assert "Template PDB id" not in src and "set /p XPDB" not in src
    # a declined dose check or a failed optional step leaves the analysis alone
    assert 'if /i "!DOSECHK_ANS!"=="n" set "DOSECHK="' in src
    assert "The dose check failed. The analysis above is unaffected." in src
    assert "The exploratory dock failed. The analysis above is unaffected." in src


def test_the_exploratory_mode_runs_only_for_a_refused_product_pocket():
    src = _launcher()
    assert 'if /i "!PSTATE!"=="refused" if not "!XUNI!"=="-"' in src
    assert "scripts" + chr(92) + "product_status.py" in src


def test_the_exploratory_dock_is_off_unless_the_user_says_yes():
    """It is unvalidated and takes hours, so it must never start by default."""
    src = _launcher()
    # the question exists, defaults to no, and comes before any stage runs
    line = next(ln for ln in src.splitlines() if "set /p XDOCK_ANS=" in ln)
    assert "y/N" in line
    assert src.index("set /p XDOCK_ANS=") < src.index("=== 1 of 5")
    # EXPLORE is set in exactly one place, inside the branch for an explicit y
    assert src.count('set "EXPLORE=1"') == 1
    yes = src.index('if /i "!XDOCK_ANS!"=="y" (')
    assert yes < src.index('set "EXPLORE=1"') < src.index(") else (", yes)
    # a refused pocket with no yes prints how to run it later
    assert "is available and is not run" in src
    assert "python -m promoter_tf_drugs.exploratory --gene !GENE!" in src
    # nothing else starts it: one command, and it is gated by EXPLORE
    assert src.count("-m promoter_tf_drugs.exploratory") == 2
    assert "if defined EXPLORE (" in src


def test_product_status_reads_usable_refused_and_missing(tmp_path, capsys):
    mod = _lf_script("product_status")
    ok = tmp_path / "ok.json"
    ok.write_text(json.dumps({"target": "G", "uniprot": "P12345",
                              "usable": True}), encoding="utf-8")
    refused = tmp_path / "refused.json"
    refused.write_text(json.dumps({"target": "G", "uniprot": "P78363",
                                   "usable": False}), encoding="utf-8")
    nouni = tmp_path / "nouni.json"
    nouni.write_text(json.dumps({"usable": False}), encoding="utf-8")
    as_list = tmp_path / "list.json"
    as_list.write_text(json.dumps([{"uniprot": "Q1", "usable": True}]),
                       encoding="utf-8")
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    assert mod.status(ok) == ("usable", "P12345")
    assert mod.status(refused) == ("refused", "P78363")
    assert mod.status(nouni) == ("refused", "-")
    assert mod.status(as_list) == ("usable", "Q1")
    assert mod.status(bad) == ("refused", "-")
    assert mod.status(tmp_path / "missing.json") == ("refused", "-")
    assert mod.main([str(refused)]) == 0
    assert capsys.readouterr().out.strip() == "refused P78363"
    assert mod.main([]) == 2


# ------------------------------------------- tissues, GWAS unpack, notices
def _table(tmp_path, rows, header="Gene\tGene name\tTissue\tNX"):
    p = tmp_path / "expr.tsv"
    p.write_text(header + NL + NL.join(rows) + NL, encoding="utf-8")
    return p


def test_tissues_lists_distinct_names_sorted_and_numbered(tmp_path, capsys):
    mod = _lf_script("tissues")
    p = _table(tmp_path, ["G1\tA\tliver\t1", "G1\tA\tB-cells\t2",
                          "G2\tB\tliver\t3", "G2\tB\tlymph node\t4"])
    assert mod.read_tissues(p) == ["B-cells", "liver", "lymph node"]
    assert mod.main([str(p)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert [ln.split(None, 1) for ln in out] == [
        ["1", "B-cells"], ["2", "liver"], ["3", "lymph node"]]


def test_tissues_pick_by_number_name_or_unique_part(tmp_path, capsys):
    mod = _lf_script("tissues")
    p = _table(tmp_path, ["G\tA\tliver\t1", "G\tA\tlymph node\t2",
                          "G\tA\tthyroid gland\t3", "G\tA\tsalivary gland\t4"])
    names = mod.read_tissues(p)
    assert mod.pick(names, "1") == "liver"
    assert mod.pick(names, "LYMPH NODE") == "lymph node"
    assert mod.pick(names, "thyroid") == "thyroid gland"
    for bad, why in (("gland", "2 tissues"), ("retina", "no tissue"),
                     ("9", "between 1 and 4"), ("", "nothing entered")):
        with pytest.raises(ValueError, match=why):
            mod.pick(names, bad)
    # the launcher reads stdout, so a failed pick prints nothing there
    assert mod.main([str(p), "--pick", "gland"]) == 1
    captured = capsys.readouterr()
    assert captured.out == "" and "matches 2 tissues" in captured.err
    assert mod.main([str(p), "--pick", "2"]) == 0
    assert capsys.readouterr().out.strip() == "lymph node"


def test_tissues_refuses_a_table_without_a_tissue_column(tmp_path, capsys):
    mod = _lf_script("tissues")
    p = _table(tmp_path, ["G\tA\t1"], header="Gene\tGene name\tNX")
    assert mod.main([str(p)]) == 2
    assert "no Tissue column" in capsys.readouterr().err
    assert mod.main([str(tmp_path / "missing.tsv")]) == 2


def test_unpack_gwas_writes_the_table_the_loader_reads(tmp_path):
    import zipfile
    mod = _lf_script("unpack_gwas")
    header = "DISEASE/TRAIT\tMAPPED_GENE\tP-VALUE\tPUBMEDID"
    row = "Stargardt disease\tABCA4\t3e-12\t123"
    zpath = tmp_path / "g.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("readme.txt", "not a table")
        zf.writestr("assoc.tsv", header + NL + row + NL)
    dest = tmp_path / "out.tsv"
    assert mod.unpack(zpath, dest) == "assoc.tsv"
    assert dest.read_text(encoding="utf-8").splitlines()[1] == row
    from promoter_tf_drugs import localannot
    assert mod.usable_header(header)
    saved = (dict(localannot._GWAS), dict(localannot._LOADED))
    try:
        localannot._GWAS.clear()
        assert localannot.load_gwas(dest) == 1
    finally:
        localannot._GWAS.clear()
        localannot._GWAS.update(saved[0])
        localannot._LOADED.update(saved[1])


def test_unpack_gwas_leaves_nothing_when_the_zip_has_no_usable_table(tmp_path,
                                                                    capsys):
    import zipfile
    mod = _lf_script("unpack_gwas")
    zpath = tmp_path / "bad.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("page.txt", "<html>captive portal</html>")
        zf.writestr("other.tsv", "a\tb" + NL + "1\t2" + NL)
    dest = tmp_path / "out.tsv"
    assert mod.main([str(zpath), str(dest)]) == 1
    assert not dest.exists() and not list(tmp_path.glob("*.part"))
    assert "no table with a gene column" in capsys.readouterr().err
    notzip = tmp_path / "notzip.zip"
    notzip.write_text("<html>", encoding="utf-8")
    assert mod.main([str(notzip), str(dest)]) == 1
    assert not dest.exists()


def test_the_launcher_offers_the_gwas_download_and_the_tissue_list():
    src = _launcher()
    bs = chr(92)
    # the offer runs only when the table is missing, and only after a yes
    assert "if not defined HAVEGWAS (" in src
    assert 'if /i not "!GETGWAS!"=="n" (' in src
    assert "ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/" in src
    assert f"scripts{bs}unpack_gwas.py" in src
    # one definition of HAVEGWAS, set before the notices that read it
    assert src.count('set "HAVEGWAS="') == 1
    assert src.index('set "HAVEGWAS="') < src.index("no GWAS Catalog table")
    # the tissue is resolved against the table, and what reaches the command
    # line is the resolved name
    assert f"scripts{bs}tissues.py" in src
    assert '--activity-tissue "!APICK!"' in src
    assert "--activity-tissue" + ' "!ATISSUE!"' not in src
    # a declined or failed offer still says where to get the files by hand
    assert "figshare.com/articles/dataset/DepMap_Public_24Q2/25880521" in src
    assert "figshare.com/articles/dataset/Repurposing_Public_24Q2/25917643" in src


def test_check_download_accepts_a_matching_checksum_and_refuses_the_rest(
        tmp_path, capsys):
    import hashlib
    mod = _lf_script("check_download")
    f = tmp_path / "t.csv.part"
    f.write_bytes(b"a,b" + NL.encode() + b"1,2" + NL.encode())
    good = hashlib.md5(f.read_bytes()).hexdigest()
    assert mod.main([str(f), good]) == 0
    assert mod.main([str(f), good.upper()]) == 0
    assert mod.main([str(f), "0" * 32]) == 1
    assert "expected" in capsys.readouterr().err
    assert mod.main([str(tmp_path / "missing"), good]) == 1
    assert mod.main([str(f)]) == 2


def test_the_launcher_fetches_depmap_and_prism_and_keeps_only_verified_files():
    """The files are downloaded under a .part name and renamed only after the
    published checksum matches. The checksums below are the ones figshare
    publishes for DepMap 24Q2 and PRISM 24Q2."""
    src = _launcher()
    bs = chr(92)
    for name, md5 in (
            ("CRISPRGeneEffect.csv", "d10d73692672380b3792231549e97d5a"),
            ("Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
             "042b15d3ebcc138a80a571b9286451ae")):
        part = f"data{bs}{name}.part"
        assert f"-o {part}" in src, name
        assert f"scripts{bs}check_download.py {part} {md5}" in src, name
        assert f"move /y {part} data{bs}{name}" in src, name
        assert f"del {part}" in src, name
    # the offers sit before the notices that read what is on disk
    assert (src.index("DepMap's CRISPR knockout scores")
            < src.index("no DepMap CRISPRGeneEffect.csv"))
    assert src.count("ndownloader.figshare.com/files/") == 4
    # the companion tables the package reads beside them
    assert "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv" in src
    assert f"data{bs}Model.csv" in src


def test_check_gmt_accepts_a_gene_set_library_and_refuses_a_page(tmp_path,
                                                                  capsys):
    mod = _lf_script("check_gmt")
    good = tmp_path / "ok.gmt"
    good.write_text(NL.join(
        f"T{i} Up" + chr(9) * 2 + chr(9).join(["G1", "G2"])
        for i in range(30)) + NL, encoding="utf-8")
    assert mod.main([str(good), "20"]) == 0
    assert mod.main([str(good), "31"]) == 1, "too few lines"
    page = tmp_path / "page.gmt"
    page.write_text("<html><body>503</body></html>" + NL * 40,
                    encoding="utf-8")
    assert mod.main([str(page), "1"]) == 1
    assert "gene sets" in capsys.readouterr().err
    assert mod.main([str(tmp_path / "missing.gmt"), "1"]) == 1
    assert mod.main([str(good)]) == 2 and mod.main([str(good), "x"]) == 2


def test_check_gmt_tolerates_a_few_malformed_lines(tmp_path):
    mod = _lf_script("check_gmt")
    lines = [f"T{i} Up" + chr(9) * 2 + "G1" for i in range(99)] + ["junk"]
    p = tmp_path / "m.gmt"
    p.write_text(NL.join(lines) + NL, encoding="utf-8")
    assert mod.main([str(p), "50"]) == 0
    lines = [f"T{i} Up" + chr(9) * 2 + "G1" for i in range(50)] + ["junk"] * 50
    p.write_text(NL.join(lines) + NL, encoding="utf-8")
    assert mod.main([str(p), "50"]) == 1


def test_the_launcher_offers_the_knockout_libraries_and_keeps_only_good_files():
    src = _launcher()
    bs = chr(92)
    assert "if not defined HAVEKO (" in src
    assert 'if /i not "!GETKO!"=="n" (' in src
    for lib, name, floor in (
            ("LINCS_L1000_CRISPR_KO_Consensus_Sigs",
             "knockout_lincs_crispr_ko.gmt", "1000"),
            ("TF_Perturbations_Followed_by_Expression",
             "knockout_tf_perturbations.gmt", "500")):
        part = f"data{bs}{name}.part"
        assert f"libraryName={lib}" in src
        assert f"scripts{bs}check_gmt.py {part} {floor}" in src
        assert f"move /y {part} data{bs}{name}" in src
        assert f"del {part}" in src
    # the file names the launcher saves are the ones the package looks for
    from promoter_tf_drugs import knockout
    for name in ("knockout_lincs_crispr_ko.gmt", "knockout_tf_perturbations.gmt"):
        assert f"data/{name}" in knockout.DEFAULT_KNOCKOUT_PATHS
