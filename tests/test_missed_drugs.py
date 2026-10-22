"""Fixes for drugs that should have been found and were not.

Found by asking, for rare disease genes whose drugs are known, where each one was
lost: a ChEMBL target type the index dropped (hydroxyurea at ribonucleotide
reductase), an approved HDAC inhibitor missing from the curated chromatin list
(givinostat), a factor two steps from the gene (regulators of regulators), and a
drug with a measured activity at a receptor but no curated mechanism (rifampin at
PXR).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

from generoutes import catalogue as cat
from generoutes import drugs as drugs_mod
from generoutes import regulators as reg
from generoutes.http import HttpError
from generoutes.network import Edge, build_graph

from test_upstream_search import _LfChembl, _lf_drug, _lf_library, _lf_target

ROOT = Path(__file__).resolve().parents[1]


def _script(name):
    spec = importlib.util.spec_from_file_location(
        f"_{name}", ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------ the target types
def test_group_fusion_and_selectivity_targets_are_indexed_and_marked():
    mechs = [
        {"molecule_chembl_id": "CHEMBL0", "target_chembl_id": "TG",
         "action_type": "INHIBITOR"},
        {"molecule_chembl_id": "CHEMBL1", "target_chembl_id": "TC",
         "action_type": "INHIBITOR"},
        {"molecule_chembl_id": "CHEMBL3", "target_chembl_id": "TS",
         "action_type": "INHIBITOR"},
        {"molecule_chembl_id": "CHEMBL4", "target_chembl_id": "TX",
         "action_type": "INHIBITOR"},
    ]
    targets = [
        _lf_target("TG", "RRM1", kind="PROTEIN COMPLEX GROUP"),
        _lf_target("TC", "ABL1", kind="CHIMERIC PROTEIN"),
        _lf_target("TS", "ACHE", kind="SELECTIVITY GROUP"),
        _lf_target("TX", "NUCL", kind="NUCLEIC-ACID"),
    ]
    idx = drugs_mod.approved_target_index(_LfChembl(mechs, targets),
                                          _lf_library(5))
    assert set(idx) == {"RRM1", "ABL1", "ACHE"}, "a nucleic acid is not a protein"
    assert idx["RRM1"][0].target_type == "PROTEIN COMPLEX GROUP"
    assert idx["ABL1"][0].target_type == "CHIMERIC PROTEIN"
    assert idx["RRM1"][0].source.endswith("group target)")


def test_a_single_protein_is_still_unmarked_and_a_family_keeps_its_label():
    mechs = [{"molecule_chembl_id": "CHEMBL0", "target_chembl_id": "T1",
              "action_type": "INHIBITOR"},
             {"molecule_chembl_id": "CHEMBL1", "target_chembl_id": "T2",
              "action_type": "INHIBITOR"}]
    idx = drugs_mod.approved_target_index(
        _LfChembl(mechs, [_lf_target("T1", "ESR1"),
                          _lf_target("T2", "PTGS2", kind="PROTEIN FAMILY")]),
        _lf_library(3))
    assert idx["ESR1"][0].target_type is None
    assert idx["PTGS2"][0].source.endswith("family target)")


# ------------------------------------------------------ chromatin agents
def _hdac_index():
    return {
        "HDAC1": [_lf_drug("Givinostat Hydrochloride", "HDAC1", "down"),
                  _lf_drug("Vorinostat", "HDAC1", "down"),
                  _lf_drug("Activator", "HDAC1", "up")],
        "HDAC2": [_lf_drug("Givinostat Hydrochloride", "HDAC2", "down")],
        "DNMT1": [_lf_drug("Decitabine", "DNMT1", "down"),
                  _lf_drug("Newstat", "DNMT1", "down")],
        "EZH2": [_lf_drug("Tazemetostat Hydrobromide", "EZH2", "down")],
        "ESR1": [_lf_drug("Notchromatin", "ESR1", "down")],
    }


def test_approved_inhibitors_of_chromatin_enzymes_are_found_from_the_index():
    extra = cat.chromatin_agents_from_index(_hdac_index())
    names = [n for n, *_ in extra]
    assert "Givinostat Hydrochloride" in names and "Newstat" in names
    assert names.count("Givinostat Hydrochloride") == 1, "listed once"
    assert "Vorinostat" not in names and "Decitabine" not in names, (
        "already on the curated list")
    assert "Tazemetostat Hydrobromide" not in names, "salt form of a curated one"
    assert "Activator" not in names and "Notchromatin" not in names
    assert dict((n, c) for n, c, *_ in extra)["Newstat"] == "DNMT inhibitor"


def test_the_extra_agents_become_chromatin_routes_that_say_where_they_came_from():
    extra = cat.chromatin_agents_from_index(_hdac_index())
    routes = cat.chromatin_routes("up", extra)
    giv = [r for r in routes if r.drug.startswith("Givinostat")]
    assert len(giv) == 1 and giv[0].category == "chromatin-level"
    assert any("ChEMBL's mechanism table" in n for n in giv[0].notes)
    curated = [r for r in routes if r.drug == "Vorinostat"]
    assert not any("ChEMBL's mechanism" in n for n in curated[0].notes)
    assert len(routes) == len(cat.CHROMATIN_DRUGS) + len(extra)
    assert cat.chromatin_routes("down", extra) == []
    out = cat.build({"gene_symbol": "G", "chromatin_extra": extra}, goal="up")
    assert any(r["drug"].startswith("Givinostat") for r in out["routes"])


def test_no_index_means_the_curated_list_alone():
    assert cat.chromatin_agents_from_index(None) == []
    assert len(cat.chromatin_routes("up")) == len(cat.CHROMATIN_DRUGS)


# ----------------------------------------------- regulators of regulators
SIGNED = {
    "AAA": {"GENE": 1},             # AAA raises the gene
    "BBB": {"GENE": -1},            # BBB lowers it
    "UP1": {"AAA": 1},              # raises AAA, so raises the gene
    "DN1": {"AAA": -1},             # lowers AAA, so lowers the gene
    "UP2": {"BBB": -1},             # lowers BBB, so raises the gene
    "BOTH": {"AAA": 1, "BBB": 1},   # raises the gene via AAA, lowers it via BBB
    "NOSIGN": {"AAA": 0},
    "OTHER": {"ZZZ": 1},
}


def test_a_factor_one_step_removed_gets_the_product_of_the_two_signs():
    got = {r.symbol: r for r in reg.candidates("GENE", SIGNED, steps=2,
                                               use_knockout=False)}
    assert got["AAA"].sign == 1 and got["AAA"].via is None
    assert got["UP1"].sign == 1 and got["UP1"].via == "AAA"
    assert got["DN1"].sign == -1
    assert got["UP2"].sign == 1 and got["UP2"].via == "BBB"
    assert "BOTH" not in got, "opposite signs by two routes: no single direction"
    assert "NOSIGN" not in got and "OTHER" not in got
    assert got["UP1"].sources == ["CollecTRI, one step removed"]


def test_one_step_is_the_default_and_adds_nothing():
    got = {r.symbol for r in reg.candidates("GENE", SIGNED,
                                            use_knockout=False)}
    assert got == {"AAA", "BBB"}


def test_a_skipped_regulator_is_still_a_way_in_but_is_not_repeated():
    got = {r.symbol: r for r in reg.candidates(
        "GENE", SIGNED, steps=2, use_knockout=False, skip=["AAA"])}
    assert "AAA" not in got, "it has a route already"
    assert got["UP1"].via == "AAA", "its regulators are still found"


def test_a_two_step_route_says_so_and_sits_one_step_further_out():
    index = {"UP1": [_lf_drug("Agonist", "UP1", "up")]}
    regs = reg.candidates("GENE", SIGNED, steps=2, use_knockout=False)
    found = reg.find_drugs(regs, "up", index)
    assert found["UP1"]["via"] == "AAA"
    (r,) = cat.from_regulators(found)
    assert r.directness == 2 and r.axes.strength == 0.4
    assert "regulates AAA" in r.mechanism
    assert any("by way of AAA" in n for n in r.notes)


# ------------------------------------------------------ measured binders
class _Activities:
    def __init__(self, pages, target_type="SINGLE PROTEIN"):
        self.pages, self.calls, self.target_type = pages, [], target_type

    def get(self, url, params=None, **kw):
        self.calls.append((url.rsplit("/", 1)[-1], dict(params or {})))
        if url.endswith("/target/search") or url.endswith("/target"):
            return {"targets": [{
                "target_chembl_id": "TPXR", "organism": "Homo sapiens",
                "target_type": self.target_type,
                "target_components": [{"accession": "O75469",
                                       "target_component_synonyms": [
                                           {"syn_type": "GENE_SYMBOL",
                                            "component_synonym": "NR1I2"}]}]}]}
        if url.endswith("/activity"):
            i = (params or {}).get("offset", 0) // 1000
            acts = self.pages[i]
            nxt = "more" if i + 1 < len(self.pages) else None
            return {"activities": acts, "page_meta": {"next": nxt}}
        raise AssertionError(url)


LIB = [{"chembl_id": f"CHEMBL{i}", "name": f"DRUG {i}"} for i in range(12)]


def _act(i, p):
    return {"molecule_chembl_id": f"CHEMBL{i}", "pchembl_value": p}


def test_measured_binders_keep_library_drugs_best_first_with_unclear_direction():
    client = _Activities([[_act(1, "5.5"), _act(2, "7.0"), _act(99, "9.0"),
                           {"molecule_chembl_id": "CHEMBL3",
                            "pchembl_value": None}],
                          [_act(2, "6.0"), _act(4, "6.5")]])
    got = drugs_mod.measured_binders(client, "NR1I2", LIB, top=2)
    assert [d.name for d in got] == ["Drug 2", "Drug 4"]
    assert got[0].mechanism == "measured activity at NR1I2, pChEMBL 7"
    assert all(d.direction == "unclear" and d.action_type == "MEASURED ACTIVITY"
               for d in got)
    asked = [p for n, p in client.calls if n == "activity"]
    assert asked[0]["pchembl_value__gte"] == 5.0 and len(asked) == 2


def test_measured_binders_need_a_single_protein_and_report_a_failed_request():
    assert drugs_mod.measured_binders(
        _Activities([[]], target_type="PROTEIN FAMILY"), "NR1I2", LIB) == []

    class Down(_Activities):
        def get(self, url, params=None, **kw):
            if url.endswith("/activity"):
                raise HttpError("500")
            return super().get(url, params, **kw)
    with pytest.raises(HttpError):
        drugs_mod.measured_binders(Down([[]]), "NR1I2", LIB)


def test_binders_are_asked_only_when_the_mechanism_index_has_nothing():
    asked = []

    def binders(sym):
        asked.append(sym)
        return [drugs_mod.Drug(
            name="Rifampin", chembl_id="C", target_symbol=sym,
            action_type=drugs_mod.MEASURED, mechanism="measured activity",
            approved=True, max_phase=4, first_approval=None, withdrawn=False,
            oral=None, direction="unclear", source="ChEMBL (measured activity)")]
    index = {"HAS": [_lf_drug("Agonist", "HAS", "up")]}
    regs = [reg.Regulator("HAS", 1, ["CollecTRI"]),
            reg.Regulator("NONE", 1, ["CollecTRI"])]
    found = reg.find_drugs(regs, "up", index, binders=binders)
    assert asked == ["NONE"]
    assert [d["name"] for d in found["HAS"]["drugs"]] == ["Agonist"]
    assert [d["name"] for d in found["NONE"]["drugs"]] == ["Rifampin"]
    assert reg.find_drugs(regs, "up", index) .keys() == {"HAS"}
    (r,) = [x for x in cat.from_regulators({"NONE": found["NONE"]})]
    assert any("only a measured activity at NONE" in n for n in r.notes)


# ---------------------------------------------------- the run and the CLI
def test_the_run_and_the_command_line_take_the_new_options():
    import inspect
    from generoutes import cli, pipeline
    sig = inspect.signature(pipeline.run).parameters
    assert sig["regulator_steps"].default == 1
    assert sig["regulator_measured"].default is False
    kw = cli.run_kwargs_from_args(cli.build_parser().parse_args(
        ["--gene", "X", "--regulator-route", "--regulator-steps", "2",
         "--regulator-measured"]))
    assert kw["regulator_steps"] == 2 and kw["regulator_measured"] is True
    with pytest.raises(Exception):
        pipeline.run("DEMOGENE", demo=True, regulator_steps=3)


def test_the_launcher_asks_for_two_steps_and_not_for_measured_activities():
    launcher = (ROOT / "scripts" / "run_everything.bat").read_text(
        encoding="utf-8")
    assert launcher.count("--regulator-route --regulator-steps 2") == 2
    assert "--regulator-measured" not in launcher


# ----------------------------------------------------------- drug_trace
def _trace(drug, *, library=None, index=None, regs=None, edges=(), ctd=None,
           routes=None, goal="up"):
    dt = _script("drug_trace")
    incoming, _, _ = build_graph([Edge(a, b, s, 1, "") for a, b, s in edges])
    return dt.trace("G", goal, drug,
                    library=library if library is not None else [
                        {"name": "DRUGX", "chembl_id": "C1"}],
                    index=index or {}, regs=regs or {}, incoming=incoming,
                    ctd=ctd, run_routes=routes)


def _steps(lines):
    return {s: ok for s, ok, _ in lines}


def test_the_trace_stops_at_the_first_missing_link(monkeypatch):
    from generoutes import fda
    monkeypatch.setattr(fda, "is_approved", lambda name: True)
    assert _steps(_trace("nothere", library=[{"name": "A", "chembl_id": "x"}])
                  ) == {1: False}
    assert _steps(_trace("drugx"))[3] is False, "no mechanism recorded"


def test_the_trace_names_the_regulator_and_checks_the_direction(monkeypatch):
    from generoutes import fda
    monkeypatch.setattr(fda, "is_approved", lambda name: True)
    index = {"TFA": [_lf_drug("Drugx", "TFA", "down")]}
    regs = {"TFA": {"sources": {"CollecTRI": "raises"}, "sign": 1}}
    lines = _trace("drugx", index=index, regs=regs)
    assert _steps(lines)[4] is True and _steps(lines)[5] is False
    assert any("the wrong way" in t for _, _, t in lines)
    right = _trace("drugx", index={"TFA": [_lf_drug("Drugx", "TFA", "up")]},
                   regs=regs)
    assert _steps(right)[5] is True
    unsigned = _trace("drugx", index=index, regs={
        "TFA": {"sources": {"ChIP occupancy": "near"}, "sign": 0}})
    assert any("no table gives the sign" in t for _, _, t in unsigned)


def test_the_trace_says_when_a_drug_is_within_reach_of_a_regulator(monkeypatch):
    from generoutes import fda
    monkeypatch.setattr(fda, "is_approved", lambda name: True)
    index = {"UPSTREAM": [_lf_drug("Drugx", "UPSTREAM", "up")]}
    regs = {"TFA": {"sources": {"CollecTRI": "raises"}, "sign": 1}}
    near = _trace("drugx", index=index, regs=regs,
                  edges=[("UPSTREAM", "TFA", 1)])
    assert _steps(near)[4] is False and _steps(near)[6] is True
    far = _trace("drugx", index=index, regs=regs, edges=[("OTHER", "TFA", 1)])
    assert _steps(far)[6] is False
