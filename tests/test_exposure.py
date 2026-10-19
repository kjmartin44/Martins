"""Exposure against need: can a person taking the drug reach the concentration
at which the evidence was seen?

The comparison is only worth making if it never turns missing data into an
answer. These tests pin that: a drug with no exposure value, or no known
effective concentration, is "not assessable", never "far below".
"""

from __future__ import annotations

import json
import logging

import pytest

from promoter_tf_drugs import exposure as ex

_REAL_WEIGHT = ex.molecular_weight      # the autouse fixture below replaces it


# ------------------------------------------------------------------ units
def test_molar_and_mass_units_convert_to_nm():
    assert ex.to_nm(5, "nM") == 5
    assert ex.to_nm(2, "uM") == 2000
    assert ex.to_nm(2, "µM") == 2000
    assert ex.to_nm(1, "mM") == 1e6
    # 100 ng/mL of a 500 g/mol drug is 0.2 uM
    assert ex.to_nm(100, "ng/mL", 500) == pytest.approx(200)
    assert ex.to_nm(1, "ug.mL-1", 300) == pytest.approx(1e6 / 300)
    assert ex.to_nm(1, "mcg/mL", 300) == pytest.approx(1e6 / 300)
    assert ex.to_nm(1, "mg/L", 300) == pytest.approx(1e6 / 300)


def test_a_mass_unit_without_a_molecular_weight_is_left_out_not_guessed():
    assert ex.to_nm(100, "ng/mL") is None
    assert ex.to_nm(100, "ng/mL", 0) is None
    assert ex.to_nm(100, "furlongs", 300) is None


def test_zero_negative_and_text_values_are_not_concentrations():
    for bad in (0, -3, "", None, "n/a"):
        assert ex.to_nm(bad, "nM") is None


def test_a_structure_that_cannot_be_read_gives_no_weight_and_says_so(caplog):
    with caplog.at_level(logging.DEBUG):
        assert _REAL_WEIGHT("not a smiles ((") is None
        assert _REAL_WEIGHT("CCO") == pytest.approx(46.07, abs=0.05)
    assert "no molecular weight" in caplog.text


# ------------------------------------------------------------- the table
def _table(tmp_path, text, name="cmax.tsv"):
    p = tmp_path / name
    p.write_text(text.replace("|", chr(9)), encoding="utf-8")
    return p


def test_a_table_gives_exposure_with_its_source(tmp_path):
    p = _table(tmp_path,
               "name|cmax|unit|mw|fu|dose|source" + chr(10)
               + "Imatinib|2.6|ug/mL|493.6|0.05|400 mg|label X" + chr(10)
               + "Aspirin|5000|nM|||325 mg|" + chr(10))
    got = ex.read_cmax_table(p)
    assert got["imatinib"]["cmax_nm"] == pytest.approx(2.6e6 / 493.6)
    assert got["imatinib"]["fu"] == 0.05 and got["imatinib"]["dose"] == "400 mg"
    assert got["imatinib"]["source"] == "label X"
    assert got["aspirin"]["source"] == "table cmax.tsv", "a default, not blank"


def test_the_table_accepts_a_cmax_nm_column_and_commas(tmp_path):
    p = _table(tmp_path, "Name,Cmax_nM" + chr(10) + "Tretinoin,300" + chr(10),
               "c.csv")
    p.write_text(p.read_text().replace(chr(9), ","), encoding="utf-8")
    assert ex.read_cmax_table(p)["tretinoin"]["cmax_nm"] == 300


def test_unreadable_rows_are_skipped_loudly_and_never_read_as_zero(tmp_path,
                                                                    caplog):
    p = _table(tmp_path, "name|cmax|unit|mw" + chr(10)
               + "Good|10|nM|" + chr(10)
               + "NoUnit|10|||" + chr(10)
               + "MassNoMw|10|ng/mL|" + chr(10)
               + "Blank|||" + chr(10))
    with caplog.at_level(logging.WARNING):
        got = ex.read_cmax_table(p)
    assert list(got) == ["good"]
    assert "NoUnit" in caplog.text and "MassNoMw" in caplog.text


def test_a_drug_listed_twice_keeps_the_larger_value(tmp_path):
    p = _table(tmp_path, "name|cmax_nM" + chr(10) + "Drug|100" + chr(10)
               + "Drug|900" + chr(10) + "Drug|300" + chr(10))
    assert ex.read_cmax_table(p)["drug"]["cmax_nm"] == 900


def test_an_unbound_fraction_outside_zero_to_one_is_ignored(tmp_path, caplog):
    p = _table(tmp_path, "name|cmax_nM|fu" + chr(10) + "D|100|5" + chr(10))
    with caplog.at_level(logging.WARNING):
        assert ex.read_cmax_table(p)["d"]["fu"] is None
    assert "not in (0, 1]" in caplog.text


# --------------------------------------------------------------- ChEMBL
class _Chembl:
    def __init__(self, activities, parent="CHEMBLP", fail=False):
        self.activities, self.parent, self.fail = activities, parent, fail
        self.calls = []

    def get(self, url, params=None, **kw):
        self.calls.append((url, params))
        if self.fail:
            raise OSError("503")
        if "/molecule/" in url:
            return {"molecule_hierarchy": {"parent_chembl_id": self.parent}}
        return {"activities": self.activities}


def _act(desc, value, units="nM"):
    return {"assay_description": desc, "standard_value": value,
            "standard_units": units}


def test_only_human_plasma_records_are_kept():
    c = _Chembl([
        _act("Cmax in human at 10 mg, po", 100),
        _act("Cmax in chronic myeloid leukemia patient", 300),
        _act("Cmax in Sprague-Dawley rat at 10 mg/kg", 5000),
        _act("Cmax in mouse plasma", 9000),
        _act("Cmax in human CSF at 500 mg", 7),
        _act("Cmax in human brain tissue", 8),
        _act("Cmax in Beagle dog", 11),
        _act("Cmax in plasma", 12),          # does not say human
    ])
    got = ex.chembl_human_cmax(c, "CHEMBL1", 400.0)
    assert got["n_records"] == 2 and got["cmax_nm"] == 200
    assert got["range_nm"] == [100, 300]
    assert got["dose"] == "10 mg"
    assert c.calls[1][1]["parent_molecule_chembl_id"] == "CHEMBLP", \
        "records are asked for under the parent, not a salt"


def test_a_record_that_says_human_and_names_an_animal_is_not_human_exposure():
    c = _Chembl([_act("Cmax in humanized mouse at 10 mg/kg", 100),
                 _act("Cmax in human xenograft rat model", 200),
                 _act("Cmax in human volunteers at 5 mg", 300)])
    got = ex.chembl_human_cmax(c, "C", 400.0)
    assert got["n_records"] == 1 and got["cmax_nm"] == 300


def test_interaction_and_special_population_records_are_not_the_approved_use():
    """Calcitriol records in ChEMBL that came from ketoconazole and ritonavir
    interaction studies read as 4.5 and 18 uM, thousands of times the real
    level, and made a 'reaches' that was not true."""
    c = _Chembl([
        _act("Cmax in healthy human co-treated with 200 mg ketoconazile", 4500),
        _act("Cmax in healthy human co-treated with 600 mg ritonavir", 18000),
        _act("Cmax in human with hepatic impairment at 10 mg", 900),
        _act("Cmax in human after rifampin pretreatment", 700),
        _act("Cmax in human volunteers at 5 mg", 20)])
    got = ex.chembl_human_cmax(c, "C", 400.0)
    assert got["n_records"] == 1 and got["cmax_nm"] == 20


def test_records_that_disagree_widely_are_not_averaged_into_a_value():
    c = _Chembl([_act("Cmax in human at 5 mg", 10),
                 _act("Cmax in human at 50 mg", 5000)])
    got = ex.chembl_human_cmax(c, "C", 400.0)
    assert got["cmax_nm"] is None and got["disagree"]
    assert got["range_nm"] == [10, 5000]


def test_disagreeing_records_leave_the_drug_not_assessable_and_say_why():
    block = _block(_routes(), hits=[
        {"drug": "Tretinoin", "dose": "10", "cell_line": "A"}])
    c = _Chembl([_act("Cmax in human at 5 mg", 10),
                 _act("Cmax in human at 50 mg", 5000)])
    d = ex.assess(block, client=c, library=LIB)["drugs"]["Tretinoin"]
    assert d["verdict"] == "not assessable"
    assert "more than 30-fold, so none is used" in d["reason"]
    assert "no human exposure value found" not in d["reason"]


def test_mass_units_use_the_weight_and_are_counted_when_they_cannot():
    rows = [_act("Cmax in human at 10 mg, po", 0.00281, "ug.mL-1")]
    got = ex.chembl_human_cmax(_Chembl(rows), "C", 558.6)
    assert got["cmax_nm"] == pytest.approx(0.00281 * 1e6 / 558.6)
    none = ex.chembl_human_cmax(_Chembl(rows), "C", None)
    assert none["cmax_nm"] is None and none["n_skipped"] == 1


def test_no_usable_record_is_none_and_a_failed_call_raises():
    assert ex.chembl_human_cmax(_Chembl([_act("Cmax in mouse", 5)]),
                                "C", 300.0) is None
    assert ex.chembl_human_cmax(_Chembl([]), None, 300.0) is None
    with pytest.raises(RuntimeError, match="ChEMBL lookup failed"):
        ex.chembl_human_cmax(_Chembl([], fail=True), "C", 300.0)


# ---------------------------------------------------------------- needed
def _block(routes=(), hits=(), grades=()):
    return {"all_routes": {"routes": [dict(r) for r in routes]},
            "lincs": {"hits": [dict(h) for h in hits]},
            "drug_grades": [dict(g) for g in grades]}


def test_the_lowest_lincs_dose_is_the_needed_concentration_in_nm():
    b = _block(hits=[{"drug": "Vorinostat", "dose": "10.0", "cell_line": "A375"},
                     {"drug": "vorinostat", "dose": "1.11", "cell_line": "HA1E"},
                     {"drug": "Other", "dose": "0", "cell_line": "X"},
                     {"drug": "NoDose", "dose": None}])
    got = ex.needed_concentrations(b)
    assert got["vorinostat"]["needed_nm"] == pytest.approx(1110)
    assert "HA1E" in got["vorinostat"]["needed_basis"]
    assert "other" not in got and "nodose" not in got


def test_potency_is_the_fallback_and_never_overrides_a_lincs_dose():
    b = _block(hits=[{"drug": "A", "dose": "10", "cell_line": "X"}],
               grades=[{"drug_name": "A", "tf_symbol": "T",
                        "evidence": {"pchembl_median": 9.0}},
                       {"drug_name": "B", "tf_symbol": "RXRA",
                        "evidence": {"pchembl_median": 7.0}},
                       {"drug_name": "C", "tf_symbol": "T",
                        "evidence": {"pchembl_median": None}}])
    got = ex.needed_concentrations(b)
    assert got["a"]["needed_nm"] == 10000, "the observed dose on the gene wins"
    assert got["b"]["needed_nm"] == pytest.approx(100)
    assert "the factor, not the gene" in got["b"]["needed_basis"]
    assert "c" not in got


# --------------------------------------------------------------- verdicts
def test_the_verdict_boundaries():
    assert ex.verdict(1000, 1000) == ("reaches", 1.0)
    assert ex.verdict(999, 1000)[0] == "within tenfold"
    assert ex.verdict(100, 1000)[0] == "within tenfold"
    assert ex.verdict(99, 1000)[0] == "far below"


def test_a_missing_side_is_not_assessable_never_far_below():
    assert ex.verdict(None, 1000) == ("not assessable", None)
    assert ex.verdict(1000, None) == ("not assessable", None)
    assert ex.verdict(0, 1000) == ("not assessable", None)


# ----------------------------------------------------------------- assess
LIB = [{"name": "IMATINIB MESYLATE", "smiles": "CC", "chembl_id": "CHEMBL1"},
       {"name": "VORINOSTAT", "smiles": "CC", "chembl_id": "CHEMBL2"},
       {"name": "TRETINOIN", "smiles": "CC", "chembl_id": "CHEMBL3"}]


@pytest.fixture(autouse=True)
def _weight(monkeypatch):
    monkeypatch.setattr(ex, "molecular_weight", lambda smi: 400.0)


def _routes():
    return [{"drug": "Imatinib Mesylate", "tier": 2, "score": 50},
            {"drug": "Imatinib", "tier": 3, "score": 90},
            {"drug": "Vorinostat", "tier": 1, "score": 10},
            {"drug": "Tretinoin", "tier": 2, "score": 70}]


def test_a_table_value_against_a_lincs_dose_gives_a_verdict(tmp_path):
    block = _block(_routes(), hits=[
        {"drug": "Vorinostat", "dose": "10", "cell_line": "A375"}])
    table = {"vorinostat": {"cmax_nm": 1200.0, "fu": None, "n_records": 1,
                            "range_nm": [1200, 1200], "dose": "400 mg",
                            "source": "label X"}}
    got = ex.assess(block, table=table, library=LIB)
    d = got["drugs"]["Vorinostat"]
    assert d["ratio"] == pytest.approx(0.12, abs=1e-3)
    assert d["verdict"] == "within tenfold"
    assert d["level"] == "total" and d["source"] == "label X"


def test_the_unbound_fraction_lowers_the_level_and_says_free():
    block = _block(_routes(), hits=[
        {"drug": "Vorinostat", "dose": "10", "cell_line": "A"}])
    table = {"vorinostat": {"cmax_nm": 2000.0, "fu": 0.1, "n_records": 1,
                            "range_nm": [2000, 2000], "dose": None,
                            "source": "s"}}
    d = ex.assess(block, table=table, library=LIB)["drugs"]["Vorinostat"]
    assert d["reached_nm"] == pytest.approx(200) and d["level"] == "free"
    assert d["verdict"] == "far below"


def test_no_exposure_or_no_need_is_not_assessable_with_the_reason():
    block = _block(_routes(), hits=[
        {"drug": "Vorinostat", "dose": "10", "cell_line": "A"}])
    got = ex.assess(block, table={}, library=LIB)["drugs"]
    v = got["Vorinostat"]
    assert v["verdict"] == "not assessable"
    assert "no human exposure value found" in v["reason"]
    assert "no effective concentration" not in v["reason"]
    t = got["Tretinoin"]
    assert "no effective concentration is known" in t["reason"]


def test_chembl_is_asked_only_for_drugs_the_table_lacks():
    block = _block(_routes(), hits=[
        {"drug": "Vorinostat", "dose": "10", "cell_line": "A"},
        {"drug": "Tretinoin", "dose": "10", "cell_line": "A"}])
    table = {"vorinostat": {"cmax_nm": 20000.0, "fu": None, "n_records": 1,
                            "range_nm": [2e4, 2e4], "dose": None,
                            "source": "t"}}
    c = _Chembl([_act("Cmax in human at 5 mg", 50)])
    got = ex.assess(block, client=c, table=table, library=LIB)
    asked = {p["parent_molecule_chembl_id"] for _, p in c.calls if p}
    assert asked == {"CHEMBLP"} and len(c.calls) == 4, \
        "imatinib and tretinoin, two calls each; vorinostat came from the table"
    assert got["drugs"]["Vorinostat"]["verdict"] == "reaches"
    assert got["drugs"]["Tretinoin"]["verdict"] == "far below"


def test_a_failed_lookup_is_reported_and_is_not_a_negative():
    block = _block(_routes(), hits=[
        {"drug": "Tretinoin", "dose": "10", "cell_line": "A"}])
    got = ex.assess(block, client=_Chembl([], fail=True), library=LIB)
    d = got["drugs"]["Tretinoin"]
    assert d["verdict"] == "not assessable" and "lookup failed" in d["reason"]
    assert got["errors"] and got["counts"]["far below"] == 0


def test_salt_forms_are_one_drug_and_top_limits_the_count():
    got = ex.assess(_block(_routes()), table={}, library=LIB, top=3)
    assert list(got["drugs"]) == ["Vorinostat", "Tretinoin",
                                  "Imatinib Mesylate"], \
        "best tier first, then score; 'Imatinib' and its mesylate are one drug"
    assert got["n_assessed"] == 3 and sum(got["counts"].values()) == 3
    assert list(ex.assess(_block(_routes()), table={}, library=LIB,
                          top=2)["drugs"]) == ["Vorinostat", "Tretinoin"]


# ----------------------------------------------------------------- output
def test_the_report_orders_by_verdict_and_keeps_the_caveats():
    block = _block(_routes(), hits=[
        {"drug": "Vorinostat", "dose": "1", "cell_line": "A"},
        {"drug": "Tretinoin", "dose": "10", "cell_line": "A"}])
    table = {"vorinostat": {"cmax_nm": 5000.0, "fu": None, "n_records": 1,
                            "range_nm": [5000, 5000], "dose": "400 mg",
                            "source": "a|b"},
             "tretinoin": {"cmax_nm": 10.0, "fu": None, "n_records": 1,
                           "range_nm": [10, 10], "dose": None, "source": "s"}}
    md = ex.markdown(ex.assess(block, table=table, library=LIB), "G1")
    assert "Total plasma drug is not free drug" in md
    assert "not the same as unreachable" in md
    assert "changes no rank" in md
    rows = [ln for ln in md.splitlines() if ln.startswith("| ")][1:]
    assert [r.split("|")[2].strip() for r in rows] == [
        "reaches", "far below", "not assessable"]
    assert "a/b" in md, "a pipe in a source name cannot break the table"


# -------------------------------------------------------------------- cli
def test_the_command_reads_a_saved_result_and_writes_both_files(tmp_path):
    res = {"query": {"gene": "G1"},
           "results": [_block(_routes(), hits=[
               {"drug": "Vorinostat", "dose": "1", "cell_line": "A"}])]}
    rp = tmp_path / "results.json"
    rp.write_text(json.dumps(res), encoding="utf-8")
    tp = _table(tmp_path, "name|cmax_nM" + chr(10) + "Vorinostat|3000" + chr(10))
    out = tmp_path / "out"
    lib = tmp_path / "lib.tsv"
    lib.write_text(chr(9).join(["name", "smiles", "chembl_id"]) + chr(10)
                   + chr(9).join(["VORINOSTAT", "CC", "CHEMBL2"]) + chr(10),
                   encoding="utf-8")
    code = ex.main([str(rp), "--cmax-table", str(tp), "--library", str(lib),
                    "--no-chembl", "--out", str(out)])
    assert code == 0
    saved = json.loads((out / "exposure.json").read_text())
    assert saved["drugs"]["Vorinostat"]["verdict"] == "reaches"
    assert "Exposure against need: G1" in (out / "exposure.md").read_text()


def test_promoters_are_merged_by_each_drugs_best_verdict(tmp_path):
    """The exposure does not depend on which promoter a route came from."""
    one = _block(_routes(), hits=[
        {"drug": "Vorinostat", "dose": "1", "cell_line": "A"}])
    two = _block(_routes())                       # no dose known here
    res = {"query": {"gene": "G1"}, "results": [two, one]}
    rp = tmp_path / "results.json"
    rp.write_text(json.dumps(res), encoding="utf-8")
    tp = _table(tmp_path, "name|cmax_nM" + chr(10) + "Vorinostat|3000" + chr(10))
    lib = tmp_path / "lib.tsv"
    lib.write_text(chr(9).join(["name", "smiles", "chembl_id"]) + chr(10)
                   + chr(9).join(["VORINOSTAT", "CC", "CHEMBL2"]) + chr(10),
                   encoding="utf-8")
    ex.main([str(rp), "--cmax-table", str(tp), "--library", str(lib),
             "--no-chembl", "--out", str(tmp_path)])
    saved = json.loads((tmp_path / "exposure.json").read_text())
    assert saved["drugs"]["Vorinostat"]["verdict"] == "reaches"
    assert saved["counts"]["reaches"] == 1


def test_the_command_rejects_nonsense_before_any_work(tmp_path, capsys):
    with pytest.raises(SystemExit):
        ex.main([str(tmp_path / "missing.json")])
    assert "not found" in capsys.readouterr().err
    rp = tmp_path / "r.json"
    rp.write_text("{}", encoding="utf-8")
    with pytest.raises(SystemExit):
        ex.main([str(rp), "--top", "0"])
    assert "at least 1" in capsys.readouterr().err


# ------------------------------------------- why a drug is not assessable
def test_a_patient_population_is_not_read_as_a_tissue():
    """Bexarotene's only human record read "patient with solid tumor", and
    "tumor" on its own was taken for the matrix."""
    c = _Chembl([
        _act("Cmax in patient with solid tumor at 300 mg/m2, po", 3390),
        _act("Cmax in patients with lung cancer at 100 mg", 100),
        _act("Cmax in human tumor tissue at 10 mg", 5),
        _act("Cmax in human liver at 10 mg", 6),
        _act("Cmax in human brain", 7),
        _act("Cmax in human brain concentration at 10 mg", 8),
    ])
    got = ex.chembl_human_cmax(c, "CHEMBL1", 400.0)
    assert got["n_records"] == 2
    assert got["range_nm"] == [100, 3390]


def test_the_lookup_reports_how_many_records_it_saw_and_kept():
    diag = {}
    c = _Chembl([_act("Cmax in human at 10 mg, po", 100),
                 _act("Cmax in Sprague-Dawley rat", 5000),
                 _act("Cmax in mouse", 9000)])
    ex.chembl_human_cmax(c, "CHEMBL1", 400.0, diag)
    assert diag == {"n_total": 3, "n_human": 1}
    diag = {}
    ex.chembl_human_cmax(_Chembl([]), "CHEMBL1", 400.0, diag)
    assert diag["n_total"] == 0


def test_the_molecular_weight_comes_from_chembl_when_the_table_has_none():
    class WithWeight(_Chembl):
        def get(self, url, params=None, **kw):
            if "/molecule/" in url:
                return {"molecule_hierarchy": {"parent_chembl_id": "P"},
                        "molecule_properties": {"full_mwt": "500.0"}}
            return super().get(url, params, **kw)
    c = WithWeight([_act("Cmax in human at 10 mg", 100, "ng/mL")])
    got = ex.chembl_human_cmax(c, "CHEMBL1", None)
    assert got["cmax_nm"] == pytest.approx(200)      # 100 ng/mL of 500 g/mol


def test_a_screens_measured_potency_is_a_needed_concentration():
    b = _block()
    b["screens"] = {
        "RXRA": {"hits": [{"drug": "Tretinoin", "measured_pchembl": 9.4},
                          {"drug": "Alitretinoin", "measured_pchembl": None}]},
        "RARB": {"hits": [{"drug": "Tretinoin", "measured_pchembl": 8.0},
                          {"drug": "Vorinostat", "measured_pchembl": 6.0}]}}
    got = ex.needed_concentrations(b)
    assert got["tretinoin"]["needed_nm"] == pytest.approx(10 ** -0.4)
    assert "RXRA" in got["tretinoin"]["needed_basis"], "the most potent wins"
    assert "alitretinoin" not in got
    assert got["vorinostat"]["needed_nm"] == pytest.approx(1000)


def test_a_lincs_dose_and_a_graded_potency_still_come_before_a_screen():
    b = _block(hits=[{"drug": "Vorinostat", "dose": "10", "cell_line": "A"}],
               grades=[{"drug_name": "Tretinoin", "tf_symbol": "RARA",
                        "evidence": {"pchembl_median": 7.0}}])
    b["screens"] = {"RXRA": {"hits": [
        {"drug": "Vorinostat", "measured_pchembl": 9.0},
        {"drug": "Tretinoin", "measured_pchembl": 9.0}]}}
    got = ex.needed_concentrations(b)
    assert "LINCS" in got["vorinostat"]["needed_basis"]
    assert "RARA" in got["tretinoin"]["needed_basis"]


def test_approved_drugs_are_assessed_before_research_compounds():
    routes = [{"drug": "Tool-1", "tier": 1, "score": 99, "approved": False},
              {"drug": "Tool-2", "tier": 1, "score": 98, "approved": False},
              {"drug": "Vorinostat", "tier": 4, "score": 5, "approved": True}]
    got = ex.candidate_drugs(_block(routes), 2)
    assert [n for _, n in got] == ["Vorinostat", "Tool-1"]


def test_a_drug_missing_from_the_table_uses_the_id_the_screen_found():
    block = _block([{"drug": "Obeticholic Acid", "tier": 2, "score": 5,
                     "approved": True}])
    block["screens"] = {"NR1H4": {"hits": [
        {"drug": "Obeticholic Acid", "chembl_id": "CHEMBL566315",
         "measured_pchembl": 6.8}]}}
    c = _Chembl([_act("Cmax in human at 25 mg", 700)], parent="CHEMBL566315")
    got = ex.assess(block, client=c, library=[])["drugs"]["Obeticholic Acid"]
    assert got["chembl_id"] == "CHEMBL566315" and got["cmax_nm"] == 700
    assert got["verdict"] == "reaches"
    assert not any("search" in url for url, _ in c.calls), \
        "the id was in hand, so no name search was made"


def test_each_kind_of_gap_is_named_and_counted():
    routes = [{"drug": f"Zzqx-{i}", "tier": 2, "score": 10 - i,
               "approved": True} for i in range(3)]
    block = _block(routes)
    lib = [{"name": "ZZQX-1", "smiles": "CC", "chembl_id": "CHEMBL1"},
           {"name": "ZZQX-2", "smiles": "CC", "chembl_id": "CHEMBL2"}]

    class PerDrug(_Chembl):
        def get(self, url, params=None, **kw):
            if "/molecule/CHEMBL1" in url or "/molecule/CHEMBL2" in url:
                return {"molecule_hierarchy": {
                    "parent_chembl_id": url.rsplit("/", 1)[1][:-5]}}
            if "/molecule/" in url:
                return {"molecules": []}
            parent = params["parent_molecule_chembl_id"]
            return {"activities": [] if parent == "CHEMBL1"
                    else [_act("Cmax in Sprague-Dawley rat", 5)]}
    got = ex.assess(block, client=PerDrug([]), library=lib)
    kinds = {n: d["exposure_gap"] for n, d in got["drugs"].items()}
    assert kinds == {"Zzqx-0": "no_id", "Zzqx-1": "no_records",
                     "Zzqx-2": "no_human_records"}
    assert got["gaps"]["exposure"] == {"no_records": 1, "no_human_records": 1,
                                       "no_id": 1}
    assert got["gaps"]["need"] == 3
    md = ex.markdown(got, "G")
    assert "Why the not assessable drugs could not be compared" in md
    assert "holds no Cmax record" in md and "none from people in plasma" in md
    assert "3 drug(s) have no concentration" in md
