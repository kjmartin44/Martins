"""The combined ranking: docking and the other routes in one list.

A fused list is only honest if it cannot count one thing twice, cannot let an
unvalidated docking run speak, and does not reward a tie with a lucky place.
These tests pin those, and that every docking row carries its caution.
"""

from __future__ import annotations

import json

import pytest

from promoter_tf_drugs import combined as cb


def _route(drug, tier=2, score=50.0, category="TF: orthosteric ligand",
           source="screening", gene=None, **kw):
    r = {"drug": drug, "tier": tier, "score": score, "category": category,
         "acts_on": "T", "source": source}
    if gene:
        r["gene"] = gene
    r.update(kw)
    return r


def _doc(drugs, usable=True, **kw):
    d = {"target": "G1", "usable": usable, "auc": 0.9, "ef1": 8.0,
         "n_docked": 100, "conformation_verified": False,
         "pocket": {"pdb_id": "4ABC"},
         "hits": [{"drug": n, "score": -10.0 - i, "cnn_score": 0.5}
                  for i, n in enumerate(drugs)]}
    d.update(kw)
    return d


# ---------------------------------------------------------- the routes list
def test_routes_are_ordered_by_tier_then_score_and_ties_share_a_rank():
    src = cb.routes_source([_route("A", 3, 90), _route("B", 2, 40),
                            _route("C", 2, 40), _route("D", 2, 10)])
    assert src["ranks"] == {"b": 1, "c": 1, "d": 3, "a": 4}


def test_a_drugs_best_route_stands_for_it_and_salts_are_one_drug():
    src = cb.routes_source([_route("Imatinib Mesylate", 3, 80),
                            _route("Imatinib", 2, 30)])
    assert list(src["ranks"]) == ["imatinib"]
    assert src["best"]["imatinib"]["tier"] == 2


def test_docking_made_routes_are_left_out_so_docking_is_not_counted_twice():
    src = cb.routes_source([_route("A"), _route("Dockhit", source="docking")])
    assert set(src["ranks"]) == {"a"}


def test_other_genes_routes_in_a_panel_are_left_out():
    rs = [_route("A", gene="G1"), _route("B", gene="ACTB"),
          _route("C", gene="g1"), _route("D")]
    assert set(cb.routes_source(rs, "G1")["ranks"]) == {"a", "c", "d"}


# --------------------------------------------------------- the docking list
def test_a_docking_run_that_passed_its_gate_ranks_in_its_own_order():
    s = cb.docking_source(_doc(["X", "Y", "Z"]), "G1")
    assert s["ranks"] == {"x": 1, "y": 2, "z": 3}
    assert s["auc"] == 0.9 and s["pocket"] == "4ABC"


def test_a_docking_run_that_failed_its_gate_is_not_used_and_says_why():
    s = cb.docking_source(_doc(["X"], usable=False, error="no pocket"), "G1")
    assert "ranks" not in s and s["unused"] == "no pocket"
    s2 = cb.docking_source(_doc(["X"], usable=False), "G1")
    assert "validation gate did not pass" in s2["unused"]


def test_a_run_with_no_hits_is_not_used():
    assert "no hits" in cb.docking_source(_doc([]), "G1")["unused"]


# ------------------------------------------------------------------ fusion
def _fuse(route_list, docs):
    routes = cb.routes_source(route_list)
    routes["label"] = "x"
    return cb.fuse(routes, [cb.docking_source(d, "G1") for d in docs])


def test_a_drug_in_both_lists_beats_one_in_a_single_list_at_the_same_rank():
    rows = _fuse([_route("Both", 2, 60), _route("RouteOnly", 2, 50)],
                 [_doc(["Both", "DockOnly"])])
    order = [r["key"] for r in rows]
    assert order[0] == "both"
    assert rows[0]["n_sources"] == 2 and rows[0]["rank"] == 1


def test_the_score_is_reciprocal_rank_with_the_conventional_constant():
    rows = _fuse([_route("A", 2, 60)], [_doc(["A"])])
    assert rows[0]["score"] == pytest.approx(2 / (cb.RRF_K + 1), abs=1e-6)
    assert cb.RRF_K == 60


def test_equal_scores_share_a_rank_and_the_next_rank_skips():
    rows = _fuse([_route("A", 2, 60), _route("B", 2, 50)],
                 [_doc(["C", "D"])])
    # A and C are each first in their list; B and D each second
    assert [(r["key"], r["rank"]) for r in rows] == [
        ("a", 1), ("c", 1), ("b", 3), ("d", 3)]


def test_ties_are_broken_by_name_not_by_input_order():
    rows = _fuse([_route("B", 2, 50), _route("A", 2, 50)], [])
    assert [r["key"] for r in rows] == ["a", "b"]
    assert [r["rank"] for r in rows] == [1, 1]


def test_a_failed_docking_run_changes_nothing():
    base = _fuse([_route("A", 2, 60), _route("B", 2, 50)], [])
    withbad = _fuse([_route("A", 2, 60), _route("B", 2, 50)],
                    [_doc(["Z"], usable=False)])
    assert [(r["key"], r["score"]) for r in base] == [
        (r["key"], r["score"]) for r in withbad]


def test_several_docking_files_each_add_a_list():
    rows = _fuse([], [_doc(["A", "B"]), _doc(["A", "C"])])
    assert rows[0]["key"] == "a" and rows[0]["n_sources"] == 2
    assert set(rows[0]["source_ranks"]) == {"docking: G1", "docking: G1 (2)"},         "two files for one target are two lists, not one overwriting the other"


# ----------------------------------------------------------------- cautions
def test_a_docking_row_says_binding_is_not_direction():
    rows = _fuse([], [_doc(["X"])])
    assert "docking: binding, not direction" in cb.flags(rows[0])
    verified = _fuse([], [_doc(["X"], conformation_verified=True)])
    assert "docking: binding, not direction" not in cb.flags(verified[0])


def test_route_cautions_reach_the_row():
    r = _route("A", breadth_penalty=0.5, mitotic_penalty=0.5, burden="high")
    rows = _fuse([r], [])
    f = cb.flags(rows[0])
    assert "moves many genes" in f and "halts mitosis" in f
    assert "high clinical burden" in f and "one source only" in f


def test_a_systemic_class_is_flagged_even_for_a_drug_found_only_by_docking():
    rows = _fuse([], [_doc(["Plicamycin", "Maraviroc"])])
    by = {r["key"]: cb.flags(r) for r in rows}
    assert any("acts across the body" in f for f in by["plicamycin"])
    assert not any("acts across the body" in f for f in by["maraviroc"])


# ------------------------------------------------------------------- report
def test_the_report_names_each_source_and_the_unused_one():
    res = cb.build(
        panel={"query": {"genes": ["G1", "ACTB"]},
               "surviving": [_route("A", gene="G1"), _route("B", gene="ACTB")]},
        dockings=[("good", _doc(["A", "X"])),
                  ("bad", _doc(["Y"], usable=False, error="no ligand"))])
    md = cb.markdown(res)
    assert "the panel's surviving routes for G1" in md
    assert "Docking `good`: used. AUC 0.9" in md
    assert "Docking `bad`: **not used**, no ligand" in md
    assert "binder and not known to activate or block" in md
    assert res["gene"] == "G1" and res["n_ranked"] == 2   # A and X, not B
    assert "reciprocal rank" in md and "not a probability" in md


def test_without_a_panel_the_runs_own_routes_are_used():
    results = {"query": {"gene": "G1"},
               "results": [{"all_routes": {"routes": [
                   _route("A"), _route("Dock", source="docking")]}}]}
    res = cb.build(results=results, dockings=[])
    assert "every route in the run for G1" in cb.markdown(res)
    assert [t["drug"] for t in res["table"]] == ["A"]


# ---------------------------------------------------------------------- cli
def test_the_command_writes_all_three_files(tmp_path):
    rp = tmp_path / "results.json"
    rp.write_text(json.dumps({"query": {"gene": "G1"}, "results": [
        {"all_routes": {"routes": [_route("A"), _route("B", 3, 20)]}}]}))
    dp = tmp_path / "product.json"
    dp.write_text(json.dumps(_doc(["B", "C"])))
    out = tmp_path / "out"
    code = cb.main(["--results", str(rp), "--docking", str(dp),
                    "--out", str(out), "--top", "5"])
    assert code == 0
    saved = json.loads((out / "combined.json").read_text())
    assert [t["drug"] for t in saved["table"]][0] == "B"
    tsv = (out / "combined.tsv").read_text().splitlines()
    assert tsv[0].split(chr(9))[:3] == ["rank", "drug", "score"]
    assert len(tsv) == 4 and "# Combined ranking: G1" in (
        out / "combined.md").read_text()


def test_the_command_rejects_nothing_to_read_and_bad_options(tmp_path, capsys):
    with pytest.raises(SystemExit):
        cb.main([])
    assert "at least one of" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        cb.main(["--results", str(tmp_path / "missing.json")])
    assert "not found" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        cb.main(["--results", "x", "--top", "0"])
    assert "at least 1" in capsys.readouterr().err


def test_only_a_class_that_really_acts_across_the_body_is_flagged(monkeypatch):
    from promoter_tf_drugs import systemic

    def fake(burden):
        return lambda name: systemic.SystemicClass(
            "a class", burden, "why", frozenset())
    rows = _fuse([], [_doc(["Anything"])])
    monkeypatch.setattr(systemic, "classify", fake(0.49))
    assert not any("acts across the body" in f for f in cb.flags(rows[0]))
    monkeypatch.setattr(systemic, "classify", fake(0.5))
    assert any("acts across the body" in f for f in cb.flags(rows[0]))
