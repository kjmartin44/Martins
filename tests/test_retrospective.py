"""The retrospective check: do predicted routes recover drugs known to move a
gene?

A check like this is only worth running if it cannot flatter the pipeline. So
these tests pin the ways it could: grading routes against their own input,
counting a drug that was never tested as a negative, reading a mixed
direction as a clean one, and counting reruns of one gene as many genes.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

from generoutes import retrospective as rt


# --------------------------------------------------------------- metrics
def _arr(*v):
    return np.array(v, dtype=float)


def test_auc_is_one_for_perfect_zero_for_reversed_and_half_for_ties():
    y = np.array([1, 1, 0, 0, 0])
    assert rt.auc(_arr(9, 8, 3, 2, 1), y) == 1.0
    assert rt.auc(_arr(1, 2, 3, 8, 9), y) == 0.0
    assert rt.auc(_arr(5, 5, 5, 5, 5), y) == 0.5
    assert rt.auc(_arr(1, 1, 1), np.array([1, 1, 1])) is None
    assert rt.auc(_arr(1, 1, 1), np.array([0, 0, 0])) is None


def test_tied_ranks_are_shared():
    assert list(rt._ranks(_arr(3, 1, 3, 2))) == [3.5, 1.0, 3.5, 2.0]


def test_the_permutation_p_is_small_for_a_real_signal_and_large_for_none():
    y = np.array([1] * 6 + [0] * 24)
    perfect = np.array([9.0] * 6 + [1.0] * 24)
    a = rt.auc(perfect, y)
    assert rt.permutation_p(perfect, y, a, 500, 0) == 1 / 501,         "never zero: the observed labelling counts as one of the shuffles"
    flat = np.zeros(30)
    assert rt.permutation_p(flat, y, rt.auc(flat, y), 500, 0) == 1.0
    assert rt.permutation_p(perfect, y, a, 0, 0) is None
    assert (rt.permutation_p(perfect, y, a, 200, 3)
            == rt.permutation_p(perfect, y, a, 200, 3)), "seeded"


# ----------------------------------------------------------- route groups
def _route(drug, category, tier=2, score=50.0, acts_on="T", mech="m"):
    return {"drug": drug, "category": category, "tier": tier, "score": score,
            "acts_on": acts_on, "mechanism": mech}


ROUTES = [
    _route("Dex", "TF: orthosteric ligand", tier=2, score=40),
    _route("Dex Sodium Phosphate", "upstream signalling", tier=3, score=90),
    _route("Statin", "TF: abundance change", tier=2, score=70),
    _route("Statin", "TF: inferred activity change", tier=2, score=65),
    _route("Vorino", "empirical signature", tier=1, score=99),
    _route("Binder", "direct: predicted product binder", tier=1, score=99),
    _route("Hit", "virtual screen hypothesis", tier=1, score=99),
    _route("Chrom", "chromatin-level", tier=4, score=10),
    _route("Dna", "promoter DNA binder", tier=3, score=20),
]


def test_each_group_takes_only_its_own_categories():
    pred = rt.route_scores(ROUTES, "predicted")
    assert set(pred) == {"dex", "statin", "chrom", "dna"}
    no_l = rt.route_scores(ROUTES, "predicted_no_lincs")
    assert "statin" not in no_l, "built from LINCS, so not graded against it"
    assert set(no_l) == {"dex", "chrom", "dna"}
    assert set(rt.route_scores(ROUTES, "measured")) == {"vorino"}


def test_the_gene_product_and_virtual_screen_routes_are_in_no_group():
    for g in rt.GROUPS:
        assert not {"binder", "hit"} & set(rt.route_scores(ROUTES, g))
    anyr = {d["key"] for d in rt.drug_ranking(ROUTES)}
    assert not {"binder", "hit"} & anyr


def test_a_better_tier_beats_a_better_score_and_salts_are_one_drug():
    sc = rt.route_scores(ROUTES, "predicted")
    assert sc["dex"] == (10 - 2) * 1000 + 40, "tier 2 beats tier 3 at 90"
    assert "dexsodiumphosphate" not in sc


def test_nothing_graded_against_lincs_is_built_from_lincs():
    for g in rt.KEY_GROUPS["lincs"]:
        for cat in rt.LINCS_DERIVED + rt.MEASURED_CATEGORIES:
            assert not rt._in_group(cat, g)


def test_the_route_count_baseline_counts_distinct_routes():
    r = [_route("A", "c1", acts_on="x"), _route("A", "c1", acts_on="x"),
         _route("A", "c2", acts_on="x"), _route("B", "c1", acts_on="x")]
    assert rt.route_counts(r) == {"a": 2.0, "b": 1.0}


def test_drug_ranking_orders_by_tier_then_score_and_names_the_route():
    order = rt.drug_ranking(ROUTES, "predicted")
    assert [d["key"] for d in order] == ["statin", "dex", "dna", "chrom"]
    assert order[0]["category"] == "TF: abundance change"
    assert order[1]["category"] == "TF: orthosteric ligand",         "the better of a salt form's routes names the drug"


# ---------------------------------------------------------- the two tests
def _ks(name="ctd", wanted=(), opposite=(), unclear=(), universe=None):
    return rt.KeySet(name, wanted=set(wanted), opposite=set(opposite),
                     unclear=set(unclear), n_records=10, universe=universe)


def _run(routes, keys):
    return rt.Run(gene="G", goal="up", path="x", routes=routes,
                  keys={k.name: k for k in keys}, skipped={})


def _lib(n=60):
    return {f"d{i}" for i in range(n)}


def _routes_for(drugs, category="upstream signalling", base=50):
    return [_route(d.upper(), category, tier=2, score=base + i)
            for i, d in enumerate(drugs)]


def test_known_drugs_with_the_best_routes_score_one():
    pos = [f"d{i}" for i in range(6)]
    run = _run(_routes_for(pos), [_ks(wanted=pos)])
    got = rt.evaluate_key(run, run.keys["ctd"], _lib(), permutations=200,
                          seed=0)
    t = got["groups"]["predicted"]["reported_or_not"]
    assert t["auc"] == 1.0 and t["n_pos"] == 6 and t["n_neg"] == 54
    assert t["n_pos_with_route"] == 6 and t["p"] < 0.01
    assert got["groups"]["predicted"]["top"]["wanted"] == 6


def test_a_mixed_direction_drug_is_neither_a_positive_nor_a_negative():
    pos = [f"d{i}" for i in range(6)]
    run = _run([], [_ks(wanted=pos + ["d7"], opposite=["d7", "d8"],
                        unclear=["d9"])])
    got = rt.evaluate_key(run, run.keys["ctd"], _lib(), permutations=0,
                          seed=0)
    t = got["groups"]["predicted"]["reported_or_not"]
    assert t["n_pos"] == 6, "d7 is in both lists and is not counted"
    assert t["n_neg"] == 60 - 6 - 3, "d7, d8 and d9 are recorded, so not negatives"
    d = got["groups"]["predicted"]["right_direction"]
    assert d["n_neg"] == 1, "d8 alone moved it the other way; d7 moved both"


def test_the_direction_test_uses_only_drugs_seen_moving_the_gene():
    pos = [f"d{i}" for i in range(6)]
    opp = [f"d{i}" for i in range(10, 17)]
    routes = _routes_for(pos) + _routes_for(opp, base=1)
    run = _run(routes, [_ks(wanted=pos, opposite=opp)])
    t = rt.evaluate_key(run, run.keys["ctd"], _lib(), permutations=0,
                        seed=0)["groups"]["predicted"]["right_direction"]
    assert t["auc"] == 1.0 and t["n_pos"] == 6 and t["n_neg"] == 7
    # the same routes with the labels swapped: the wrong way round scores 0
    run2 = _run(routes, [_ks(wanted=opp, opposite=pos)])
    t2 = rt.evaluate_key(run2, run2.keys["ctd"], _lib(), permutations=0,
                         seed=0)["groups"]["predicted"]["right_direction"]
    assert t2["auc"] == 0.0


def test_a_class_with_too_few_members_is_not_scored_and_says_so():
    run = _run(_routes_for(["d0", "d1"]), [_ks(wanted=["d0", "d1"])])
    t = rt.evaluate_key(run, run.keys["ctd"], _lib(), permutations=0,
                        seed=0)["groups"]["predicted"]["reported_or_not"]
    assert t["auc"] is None and "needs 5 of each class" in t["reason"]


def test_a_drug_never_tested_is_not_a_negative():
    pos = [f"d{i}" for i in range(6)]
    tested = set(pos) | {f"d{i}" for i in range(20, 30)}
    run = _run(_routes_for(pos), [_ks("lincs", wanted=pos, universe=tested)])
    t = rt.evaluate_key(run, run.keys["lincs"], _lib(), permutations=0,
                        seed=0)["groups"]["predicted_no_lincs"][
        "reported_or_not"]
    assert t["n_neg"] == 10, "only drugs the signature files contain count"


def test_positives_reached_only_by_lincs_built_routes_have_no_route_for_lincs():
    pos = [f"d{i}" for i in range(6)]
    routes = _routes_for(pos, category="TF: abundance change")
    run = _run(routes, [_ks("lincs", wanted=pos, universe=_lib())])
    t = rt.evaluate_key(run, run.keys["lincs"], _lib(), permutations=0,
                        seed=0)["groups"]["predicted_no_lincs"][
        "reported_or_not"]
    assert t["n_pos_with_route"] == 0 and t["auc"] == 0.5


def test_a_baseline_that_knows_nothing_of_the_gene_is_reported():
    pos = [f"d{i}" for i in range(6)]
    run = _run(_routes_for(pos), [_ks(wanted=pos)])
    t = rt.evaluate_key(run, run.keys["ctd"], _lib(), permutations=0,
                        seed=0)["groups"]["predicted"]["reported_or_not"]
    assert t["baseline_auc"] == 1.0, "drugs with routes are the positives here"
    assert "baseline_auc" not in rt.evaluate_key(
        run, run.keys["ctd"], _lib(), permutations=0, seed=0)["groups"][
        "predicted"]["right_direction"]


# -------------------------------------------------------------- the keys
def _blocks(**imp):
    base = {"reported_links": {
        "wanted_direction": [{"chemical": "Dexamethasone"},
                             {"chemical": "Aspirin Sodium"}],
        "opposite_direction": [{"chemical": "Arsenic"}],
        "unclear_direction": []}, "ctd_links_found": 5}
    base.update(imp)
    return [{"improvements": base}]


def test_ctd_links_become_normalised_drug_names():
    ks, why = rt.ctd_keyset(_blocks())
    assert why is None and ks.wanted == {"dexamethasone", "aspirin"}
    assert ks.opposite == {"arsenic"} and ks.n_records == 5


def test_an_older_result_without_the_flag_counts_as_read():
    assert rt.ctd_keyset(_blocks())[0] is not None


def test_ctd_that_was_not_read_or_holds_nothing_is_not_built():
    assert rt.ctd_keyset(_blocks(ctd_assessed=False))[1] \
        == "CTD was not read for this run"
    assert "no expression record" in rt.ctd_keyset(
        _blocks(ctd_links_found=0))[1]
    assert "no CTD record" in rt.ctd_keyset([{"improvements": {}}])[1]


def _gmt(path, rows):
    path.write_text(chr(10).join(
        chr(9).join([name, ""] + genes) for name, genes in rows) + chr(10),
        encoding="utf-8")
    return path


@pytest.fixture
def gmts(tmp_path):
    up = _gmt(tmp_path / "up.gmt", [
        ("LJP001 A375 24H-RAISER-10", ["G1", "X"]),
        ("LJP002 HA1E 24H-RAISER-10", ["G1"]),
        ("LJP003 HA1E 24H-RAISER-10", ["G1"]),
        ("LJP004 HA1E 24H-MIXED-10", ["G1"]),
        ("LJP005 HA1E 24H-MIXED-10", ["G1"]),
        ("LJP006 HA1E 24H-ONCE-10", ["G1"]),
        ("LJP007 HA1E 24H-LOWERER-10", ["X"]),
        ("LJP008 HA1E 24H-NEVERSEEN-10", ["X"])])
    down = _gmt(tmp_path / "down.gmt", [
        ("LJP001 A375 24H-LOWERER-10", ["G1"]),
        ("LJP002 HA1E 24H-LOWERER-10", ["G1"]),
        ("LJP003 HA1E 24H-MIXED-10", ["G1"]),
        ("LJP004 HA1E 24H-MIXED-10", ["G1"]),
        ("LJP005 HA1E 24H-RAISER-10", ["X"]),
        ("LJP008 HA1E 24H-NEVERSEEN-10", ["X"])])
    return str(up), str(down)


def test_lincs_calls_a_direction_only_when_it_is_repeated_and_lopsided(gmts):
    ks, why = rt.lincs_keyset("G1", "up", *gmts)
    assert why is None
    assert ks.wanted == {"raiser"}
    assert ks.opposite == {"lowerer"}
    assert ks.unclear == {"mixed", "once"}, "mixed signs and a single sighting"
    assert "neverseen" in ks.universe and "raiser" in ks.universe


def test_the_goal_decides_which_side_is_wanted(gmts):
    ks, _ = rt.lincs_keyset("G1", "down", *gmts)
    assert ks.wanted == {"lowerer"} and ks.opposite == {"raiser"}


def test_lincs_is_not_built_without_files_or_for_a_gene_it_lacks(gmts,
                                                                 tmp_path):
    missing = str(tmp_path / "nope.gmt")
    assert rt.lincs_keyset("G1", "up", missing, missing)[1] \
        == "no local LINCS .gmt files"
    assert "no LINCS signature contains NOPE" in rt.lincs_keyset(
        "NOPE", "up", *gmts)[1]
    assert "not up or down" in rt.lincs_keyset("G1", "sideways", *gmts)[1]


# -------------------------------------------------------------- the runs
def _write(tmp_path, name, gene="G1", goal="up", routes=None, **imp):
    d = tmp_path / name
    d.mkdir()
    res = {"query": {"gene": gene, "goal": goal},
           "results": [{"all_routes": {"routes": routes or ROUTES},
                        **{"improvements": _blocks(**imp)[0]["improvements"]}}]}
    p = d / "results.json"
    p.write_text(json.dumps(res), encoding="utf-8")
    return p


def test_a_saved_result_loads_with_both_keys_or_says_why_not(tmp_path, gmts):
    p = _write(tmp_path, "r1")
    both = rt.load_run(p.parent, gmt_up=gmts[0], gmt_down=gmts[1])
    assert set(both.keys) == {"ctd", "lincs"} and both.gene == "G1"
    nogmt = rt.load_run(p, gmt_up=str(tmp_path / "x"),
                        gmt_down=str(tmp_path / "y"))
    assert set(nogmt.keys) == {"ctd"}
    assert nogmt.skipped["lincs"] == "no local LINCS .gmt files"


def test_an_unreadable_or_empty_result_is_skipped_not_fatal(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    assert rt.load_run(bad) is None
    empty = tmp_path / "empty.json"
    empty.write_text(json.dumps({"results": []}), encoding="utf-8")
    assert rt.load_run(empty) is None


def test_reruns_of_one_gene_are_one_gene_and_the_newest_is_kept(tmp_path,
                                                                gmts):
    old = _write(tmp_path, "g1_old")
    new = _write(tmp_path, "g1_new")
    other = _write(tmp_path, "g1_down", goal="down")
    os.utime(old, (1_000_000, 1_000_000))
    os.utime(new, (2_000_000, 2_000_000))
    runs = [rt.load_run(p, gmt_up=gmts[0], gmt_down=gmts[1])
            for p in (old, new, other)]
    kept, dropped = rt.latest_per_gene_goal(runs)
    assert sorted(r.path for r in kept) == sorted([str(new), str(other)])
    assert dropped == [str(old)]


# ------------------------------------------------------ positive controls
def _control_run(gene, goal, routes):
    return rt.Run(gene=gene, goal=goal, path="p", routes=routes, keys={},
                  skipped={})


def test_a_control_reports_where_the_known_drug_landed_and_by_which_route():
    routes = [_route("Other", "TF: orthosteric ligand", tier=1, score=90),
              _route("Rifampin", "TF: orthosteric ligand", tier=2, score=60,
                     acts_on="NR1I2"),
              _route("Rifampicin SV", "empirical signature", tier=1, score=99)]
    c = rt.control_result(_control_run("CYP3A4", "up", routes))
    any_r, pred = c["views"]["any route"], c["views"]["predicted routes"]
    assert any_r["best_rank"] == 3 and any_r["n_drugs"] == 3,         "a tier 1 signature route and a tier 1 route rank above tier 2"
    assert any_r["best"]["acts_on"] == "NR1I2"
    assert pred["best_rank"] == 2 and pred["n_drugs"] == 2, \
        "the signature route is not a prediction"


def test_a_control_where_the_drug_is_absent_says_absent_not_zero():
    c = rt.control_result(_control_run(
        "FKBP5", "up", [_route("Other", "upstream signalling")]))
    v = c["views"]["any route"]
    assert v["found"] is False and v["best_rank"] is None and v["best"] is None


def test_only_the_right_gene_and_goal_is_a_control():
    r = [_route("Rifampin", "upstream signalling")]
    assert rt.control_result(_control_run("CYP3A4", "down", r)) is None
    assert rt.control_result(_control_run("ACTB", "up", r)) is None
    assert rt.control_result(_control_run("cyp3a4", "up", r)) is not None


def test_every_control_names_its_goal_drugs_and_mechanism():
    for gene, (goal, names, why) in rt.POSITIVE_CONTROLS.items():
        assert gene == gene.upper() and goal in ("up", "down")
        assert names and why


# ------------------------------------------------------------- the report
def test_the_summary_counts_runs_above_chance_and_positives_reached(tmp_path,
                                                                    gmts):
    good = [f"d{i}" for i in range(6)]
    runs = []
    for name, wanted, routes in (("a", good, _routes_for(good)),
                                 ("b", good, [])):
        run = _run(routes, [_ks(wanted=wanted)])
        run.gene = name
        runs.append(rt.evaluate_run(run, _lib(), permutations=0))
    s = rt.summarise(runs)["ctd"]["predicted"]["reported_or_not"]
    assert s["n_runs"] == 2 and s["n_above_chance"] == 1
    assert s["n_positives"] == 12 and s["n_positives_with_route"] == 6
    assert s["min_auc"] == 0.5 and s["max_auc"] == 1.0


def test_the_report_states_its_limits_and_coverage():
    pos = [f"d{i}" for i in range(6)]
    run = _run(_routes_for(pos[:2]), [_ks(wanted=pos)])
    run.skipped = {"lincs": "no local LINCS .gmt files"}
    rep = {"runs": [rt.evaluate_run(run, _lib(), permutations=0)],
           "dropped": ["old/results.json"], "controls": [
               rt.control_result(_control_run(
                   "CYP3A4", "up", [_route("Rifampin", "upstream signalling",
                                           acts_on="NR1I2")]))]}
    rep["summary"] = rt.summarise(rep["runs"])
    md = rt.markdown(rep)
    assert "2 of 6 positives have a route" in md
    assert "not built, no local LINCS .gmt files" in md
    assert "A low AUC with few positives that have a route is a coverage" in md
    assert "unstudied, not inert" in md
    assert "the unrecorded candidates are right" in md
    assert "## Positive controls" in md and "rank 1 of 1" in md
    assert "old/results.json" in md


# -------------------------------------------------------------------- cli
def test_the_command_scores_saved_results_and_writes_both_files(tmp_path,
                                                                gmts):
    p = _write(tmp_path, "run1", gene="CYP3A4", routes=[
        _route("Rifampin", "TF: orthosteric ligand", acts_on="NR1I2")])
    lib = tmp_path / "lib.tsv"
    lib.write_text(chr(10).join(
        [chr(9).join(["name", "smiles", "chembl_id"])]
        + [chr(9).join([f"DRUG{i}", "CC", f"CHEMBL{i}"]) for i in range(8)])
        + chr(10), encoding="utf-8")
    out = tmp_path / "out"
    code = rt.main([str(p.parent), "--library", str(lib), "--gmt-up", gmts[0],
                    "--gmt-down", gmts[1], "--permutations", "10",
                    "--out", str(out)])
    assert code == 0
    saved = json.loads((out / "retrospective.json").read_text())
    assert saved["controls"][0]["gene"] == "CYP3A4"
    assert saved["controls"][0]["views"]["any route"]["best_rank"] == 1
    assert "# Retrospective check" in (out / "retrospective.md").read_text()


def test_the_command_rejects_unreadable_input_and_bad_options(tmp_path,
                                                              capsys):
    lib = tmp_path / "lib.tsv"
    lib.write_text(chr(10).join([chr(9).join(["name", "smiles"]),
                                 chr(9).join(["A", "CC"])]) + chr(10),
                   encoding="utf-8")
    with pytest.raises(SystemExit):
        rt.main([str(tmp_path / "missing"), "--library", str(lib)])
    assert "none of the given paths held a readable result" in \
        capsys.readouterr().err
    with pytest.raises(SystemExit):
        rt.main([str(tmp_path), "--library", str(lib), "--permutations",
                 "-1"])
    assert "cannot be negative" in capsys.readouterr().err


# ------------------------------------------------------------- drug fame
def _ctd_run(gene, wanted, opposite=(), routes=()):
    return rt.Run(gene=gene, goal="up", path="p", routes=list(routes),
                  keys={"ctd": _ks(wanted=wanted, opposite=opposite)},
                  skipped={})


def test_fame_counts_how_many_genes_record_a_drug_and_waits_for_enough_genes():
    runs = [_ctd_run(f"G{i}", ["a", "b"] if i < 3 else ["a"]) for i in range(5)]
    got = rt.fame_counts(runs)
    assert got == {"a": 5, "b": 3}
    assert rt.fame_counts(runs[:4]) == {}, "four genes say nothing about fame"
    assert rt.fame_counts([]) == {}


def test_a_genes_own_record_is_taken_out_of_its_fame():
    common = [f"d{i}" for i in range(6)]
    runs = [_ctd_run(f"G{i}", common) for i in range(5)]
    runs[0].keys["ctd"].wanted.add("d7")
    fame = rt.fame_counts(runs)
    assert fame["d7"] == 1 and fame["d0"] == 5
    got = rt.evaluate_run(runs[0], _lib(), permutations=0, fame=fame)
    t = got["keys"]["ctd"]["groups"]["predicted"]["reported_or_not"]
    # d7 is recorded by this gene alone, so once its own record is taken out
    # its fame is 0 and it is not "known"; d0 to d5 keep four other genes.
    # With fame left in, d7 would look known and the matched set would hold
    # seven positives. Taken out, six positives are known.
    assert t["fame_auc"] is not None
    assert t.get("matched_n_pos") in (None, 6) and t["n_pos"] == 7


def test_a_score_that_only_prefers_famous_drugs_is_exposed_by_matching_on_fame():
    keys = [f"d{i}" for i in range(80)]
    famous = set(keys[:40])
    fame = {k: (10 if k in famous else 0) for k in keys}
    # positives are mostly the famous drugs, and so is the score
    labels = [1 if (k in famous and i % 2 == 0) or i in (60, 61, 62, 63, 64)
              else 0 for i, k in enumerate(keys)]
    scores = {k: 100.0 for k in famous}
    t = rt._score_test(keys, labels, scores, permutations=0, seed=0,
                       fame=fame)
    assert t["auc"] > 0.7 and t["fame_auc"] > 0.7
    # among the famous drugs alone the score is a constant: no skill
    assert t["matched_auc"] == 0.5
    assert t["matched_n_pos"] == 20 and t["matched_n_neg"] == 20


def test_matching_on_fame_is_withheld_when_too_few_drugs_are_known():
    keys = [f"d{i}" for i in range(40)]
    fame = {k: 0 for k in keys}
    fame["d0"] = fame["d1"] = 9
    labels = [1 if i < 8 else 0 for i in range(40)]
    t = rt._score_test(keys, labels, {}, permutations=0, seed=0, fame=fame)
    assert t["matched_auc"] is None


def test_the_report_adds_a_fame_section_only_when_enough_genes_were_scored():
    pos = [f"d{i}" for i in range(8)]
    runs = [_ctd_run(f"G{i}", pos, routes=_routes_for(pos)) for i in range(5)]
    library = {f"d{i}" for i in range(60)}
    fame = rt.fame_counts(runs)
    per = [rt.evaluate_run(r, library, permutations=0, fame=fame)
           for r in runs]
    rep = {"runs": per, "summary": rt.summarise(per), "controls": [],
           "dropped": [], "fame_used": True}
    md = rt.markdown(rep)
    assert "### Drug fame" in md
    assert "Too few drugs recorded for at least 3 other genes" in md,         "every drug here is well known, so nothing is left to compare"
    s = rep["summary"]["ctd"]["predicted"]["reported_or_not"]
    assert s["mean_fame_auc"] is not None
    few = [rt.evaluate_run(r, library, permutations=0) for r in runs[:2]]
    rep2 = {"runs": few, "summary": rt.summarise(few), "controls": [],
            "dropped": [], "fame_used": False}
    assert "### Drug fame" not in rt.markdown(rep2)


def test_the_report_tells_the_reader_which_figure_to_trust_when_it_has_one():
    rep = {"runs": [], "controls": [], "dropped": [], "fame_used": True,
           "summary": {k: {g: {t: {"n_runs": 4, "mean_auc": 0.62,
                                   "min_auc": 0.5, "max_auc": 0.7,
                                   "n_above_chance": 3, "n_positives": 40,
                                   "n_positives_with_route": 20,
                                   "mean_baseline_auc": 0.6,
                                   "mean_fame_auc": 0.9,
                                   "mean_matched_auc": 0.47,
                                   "n_matched_runs": 8}
                               for t in ("reported_or_not",
                                         "right_direction")}
                           for g in groups}
                       for k, groups in rt.KEY_GROUPS.items()}}
    md = rt.markdown(rep)
    assert "averages 0.9 on the reported-or-not test, against 0.62" in md
    assert "average 0.47 over 8 genes" in md
    assert "That is the figure to trust" in md


def test_a_genes_own_record_cannot_make_its_own_positives_look_famous():
    pos = [f"d{i}" for i in range(6)]
    run = _ctd_run("G0", pos)
    # counts as the run set would have them: only this gene records them
    got = rt.evaluate_run(run, _lib(), permutations=0,
                          fame={k: 1 for k in pos})
    t = got["keys"]["ctd"]["groups"]["predicted"]["reported_or_not"]
    assert t["fame_auc"] == 0.5, "with its own record out, nobody knows them"


def test_the_fame_baseline_is_scored_on_fame_and_not_on_the_routes_score():
    keys = [f"d{i}" for i in range(40)]
    labels = [1 if i < 8 else 0 for i in range(40)]
    fame = {k: (5 if i < 8 else 0) for i, k in enumerate(keys)}
    t = rt._score_test(keys, labels, {}, permutations=0, seed=0, fame=fame)
    assert t["auc"] == 0.5 and t["fame_auc"] == 1.0


def test_a_drug_recorded_for_exactly_the_floor_number_of_genes_is_known():
    keys = [f"d{i}" for i in range(30)]
    labels = [1 if i < 10 else 0 for i in range(30)]
    at = {k: rt.FAME_FLOOR for k in keys[:5] + keys[10:15]}
    t = rt._score_test(keys, labels, {}, permutations=0, seed=0, fame=at)
    assert t["matched_n_pos"] == 5 and t["matched_n_neg"] == 5
    below = {k: rt.FAME_FLOOR - 1 for k in keys[:5] + keys[10:15]}
    assert rt._score_test(keys, labels, {}, permutations=0, seed=0,
                          fame=below)["matched_auc"] is None


def test_fame_counts_a_drug_recorded_in_any_direction_or_as_unclear():
    runs = []
    for i in range(5):
        r = _ctd_run(f"G{i}", [])
        r.keys["ctd"].opposite = {"opp"} if i < 2 else set()
        r.keys["ctd"].unclear = {"unc"} if i < 4 else set()
        r.keys["ctd"].wanted = {"want"} if i < 3 else set()
        runs.append(r)
    assert rt.fame_counts(runs) == {"opp": 2, "unc": 4, "want": 3}


def test_the_summary_averages_the_fame_matched_scores_it_was_given():
    def rec(matched):
        t = {"auc": 0.6, "n_pos": 8, "n_neg": 40, "n_pos_with_route": 4,
             "baseline_auc": 0.55, "fame_auc": 0.9, "matched_auc": matched}
        t2 = {"auc": None, "n_pos": 0, "n_neg": 0, "n_pos_with_route": 0}
        g = {"reported_or_not": t, "right_direction": t2}
        return {"gene": "G", "goal": "up", "path": "p", "skipped": {},
                "keys": {"ctd": {"groups": {"predicted": g,
                                            "predicted_no_lincs": g,
                                            "measured": g}}}}
    s = rt.summarise([rec(0.4), rec(0.6), rec(None)])
    e = s["ctd"]["predicted"]["reported_or_not"]
    assert e["mean_matched_auc"] == 0.5 and e["n_matched_runs"] == 2
    assert e["mean_fame_auc"] == 0.9
