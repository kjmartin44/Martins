"""The funnel that says where known drugs are lost.

The classification is only worth reading if each stage means one thing, a drug
with several targets takes the stage nearest to being found, and a lookup that
failed is never read as a drug with no target.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from promoter_tf_drugs.network import Edge, build_graph


def _script():
    path = Path(__file__).resolve().parents[1] / "scripts" / "positive_fate.py"
    spec = importlib.util.spec_from_file_location("_positive_fate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pf = _script()


def graph(*edges):
    return build_graph([Edge(a, b, 1, 1, "") for a, b in edges])


def test_reach_is_counted_upstream_and_stops_at_three_steps():
    incoming, _, _ = graph(("A", "B"), ("B", "C"), ("C", "D"), ("D", "E"))
    reach = pf.upstream_reach(incoming, {"E"})
    assert reach == {"E": 0, "D": 1, "C": 2, "B": 3}, "A is four steps away"
    assert pf.upstream_reach(incoming, {"E"}, max_hops=1) == {"E": 0, "D": 1}


def test_reach_walks_against_the_edge_direction_only():
    incoming, _, _ = graph(("X", "F"), ("F", "Y"))
    assert "Y" not in pf.upstream_reach(incoming, {"F"}), \
        "Y is downstream of F, not upstream"
    assert "X" in pf.upstream_reach(incoming, {"F"})


@pytest.mark.parametrize("targets,expected", [
    (set(), pf.NO_TARGET),
    ({"RANKED"}, pf.FOUND_FACTOR),
    ({"NEAR"}, pf.IN_REACH),
    ({"FARTF"}, pf.FACTOR_MISSED),
    ({"FAR"}, pf.TOO_FAR),
    ({"ISLAND"}, pf.OFF_NETWORK),
])
def test_each_target_lands_at_one_stage(targets, expected):
    got = pf.classify(targets, ranked={"RANKED"}, tfs={"FARTF", "RANKED"},
                      reach={"NEAR": 2, "RANKED": 0},
                      in_network={"FAR", "NEAR", "RANKED", "FARTF"})
    assert got == expected


def test_a_drug_with_several_targets_takes_the_nearest_stage():
    got = pf.classify({"ISLAND", "NEAR", "FAR"}, ranked=set(), tfs=set(),
                      reach={"NEAR": 1}, in_network={"NEAR", "FAR"})
    assert got == pf.IN_REACH
    got = pf.classify({"ISLAND", "FAR"}, ranked=set(), tfs=set(), reach={},
                      in_network={"FAR"})
    assert got == pf.TOO_FAR


def test_a_drug_acting_on_a_ranked_factor_beats_every_other_target():
    got = pf.classify({"RANKED", "NEAR"}, ranked={"RANKED"}, tfs=set(),
                      reach={"NEAR": 1}, in_network={"NEAR", "RANKED"})
    assert got == pf.FOUND_FACTOR


class _Client:
    """Serves one mechanism and one target; can fail on demand."""

    def __init__(self, fail_mechanism=False):
        self.fail_mechanism = fail_mechanism
        self.calls = 0

    def get(self, url, params=None, **kw):
        self.calls += 1
        if "mechanism" in url:
            if self.fail_mechanism:
                raise pf.HttpError("503")
            mid = params["molecule_chembl_id"]
            return {"mechanisms": [{"target_chembl_id": "T1"}]
                    if mid == "CHEMBL_SALT" else []}
        return {"target_components": [{"target_component_synonyms": [
            {"syn_type": "GENE_SYMBOL", "component_synonym": "esr1"},
            {"syn_type": "GENE_SYMBOL_OTHER", "component_synonym": "ESR"}]}]}


def test_a_mechanism_on_the_salt_form_is_found_from_the_bare_name(monkeypatch):
    monkeypatch.setattr(pf.time, "sleep", lambda s: None)
    fails: list = []
    got = pf.target_symbols(_Client(), {"CHEMBL_BASE", "CHEMBL_SALT"}, {},
                            fails)
    assert got == {"ESR1"} and not fails


def test_a_failed_lookup_is_not_cached_as_no_target(monkeypatch):
    monkeypatch.setattr(pf.time, "sleep", lambda s: None)
    cache: dict = {}
    fails: list = []
    got = pf.target_symbols(_Client(fail_mechanism=True), {"CHEMBL_SALT"},
                            cache, fails)
    assert got == set() and len(fails) == 1 and not cache
    again = pf.target_symbols(_Client(), {"CHEMBL_SALT"}, cache, fails)
    assert again == {"ESR1"}, "the next attempt asks again"


def test_the_summary_counts_stages_and_names_the_targets():
    runs = [{"gene": "G", "goal": "up", "n_positives": 3, "n_ranked_factors": 5,
             "drugs": [
                 {"drug": "a", "fate": pf.FOUND, "targets": []},
                 {"drug": "b", "fate": pf.IN_REACH, "targets": ["PTGS2"]},
                 {"drug": "c", "fate": pf.IN_REACH, "targets": ["PTGS2",
                                                                "PTGS1"]}]}]
    s = pf.summarise(runs)
    assert s["n"] == 3 and s["total"] == {pf.FOUND: 1, pf.IN_REACH: 2}
    assert s["top_targets"][pf.IN_REACH][0] == ("PTGS2", 2)
    md = pf.markdown(runs, s)
    assert "| in reach upstream | 2 | 66.7% |" in md
    assert "G (up)" in md
