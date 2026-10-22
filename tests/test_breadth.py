"""How broadly a compound acts, from the LINCS signature files.

A .gmt signature lists the genes one compound moved in one condition. The
lists differ in length, so the length says how many genes were moved. That
is a caution about specificity, never a score: nothing here changes a rank.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from generoutes import pipeline, report
from generoutes import signatures as sg

NL = chr(10)
TAB = chr(9)


def _gmt(path, sets):
    """`sets` is {signature name: [genes]}; one line each, as Enrichr writes."""
    path.write_text(NL.join(TAB.join([name, ""] + genes)
                            for name, genes in sets.items()) + NL,
                    encoding="utf-8")
    return path


def _files(tmp_path):
    """Twelve signatures. `beta` moves many genes, `alpha` few, `gamma` has a
    single signature."""
    def genes(prefix, n):
        return [f"{prefix}{i}" for i in range(n)]
    up, down = {}, {}
    for i in range(4):
        n = f"P{i} A549 24H-alpha-10.0"
        up[n], down[n] = genes("UA", 5), genes("DA", 5) + ["TESTG"]
    for i in range(3):
        n = f"P{i} MCF7 24H-beta hydrochloride-10.0"
        up[n], down[n] = genes("UB", 80), genes("DB", 80) + ["TESTG"]
    n = "P9 PC3 6H-gamma-1.11"
    up[n], down[n] = genes("UG", 40), genes("DG", 40) + ["TESTG"]
    for i in range(4):
        n = f"P{i} HA1E 24H-filler-{i}.0"
        up[n], down[n] = genes("UF", 10), genes("DF", 10)
    return _gmt(tmp_path / "up.gmt", up), _gmt(tmp_path / "down.gmt", down)


@pytest.fixture(autouse=True)
def _fresh_cache(monkeypatch):
    monkeypatch.setattr(sg, "_BREADTH_CACHE", {})


# ----------------------------------------------------------- the name key
def test_spellings_of_one_compound_meet_on_one_key():
    key = sg.breadth_key
    assert key("Trichostatin A") == key("trichostatin-a") == "trichostatina"
    assert key("Emetine Dihydrochloride Hydrate") == key("emetine")
    assert key("NVP-AUY922") == key("nvp auy922")
    assert key("") == "" and key(None) == ""


# ------------------------------------------------------- the calculation
def test_each_compound_gets_its_median_genes_moved_and_a_percentile(tmp_path):
    up, down = _files(tmp_path)
    b = sg.signature_breadth(up, down)
    assert b["n_signatures"] == 12
    c = b["compounds"]
    assert c["alpha"]["median_genes"] == 11 and c["alpha"]["n_signatures"] == 4
    assert c["beta"]["median_genes"] == 161, "up and down together"
    assert c["gamma"]["median_genes"] == 81 and c["gamma"]["n_signatures"] == 1
    assert c["beta"]["percentile"] > c["gamma"]["percentile"] \
        > c["alpha"]["percentile"]
    assert c["beta"]["percentile"] == pytest.approx(100 * 9 / 12, abs=0.1), (
        "the share of all signatures that moved fewer genes")


def test_a_missing_file_gives_none_not_an_empty_answer(tmp_path):
    up, _ = _files(tmp_path)
    assert sg.signature_breadth(up, tmp_path / "nope.gmt") is None
    assert sg.breadth_for(["alpha"], None) == {}


def test_the_files_are_read_once_and_again_when_they_change(tmp_path):
    up, down = _files(tmp_path)
    first = sg.signature_breadth(up, down)
    assert sg.signature_breadth(up, down) is first
    _gmt(down, {"P0 A549 24H-alpha-10.0": ["X"] * 3})
    second = sg.signature_breadth(up, down)
    assert second is not first


def test_flags_say_broad_thin_or_neither(tmp_path):
    up, down = _files(tmp_path)
    b = sg.signature_breadth(up, down)
    got = sg.breadth_for(["Beta Hydrochloride", "alpha", "gamma", "unknown"], b)
    assert got["beta"]["broad"] is True and got["beta"]["thin"] is False
    assert got["alpha"]["broad"] is False and got["alpha"]["thin"] is False
    assert got["gamma"]["thin"] is True and got["gamma"]["broad"] is False, (
        "one signature is too few to call a compound broad")
    assert "unknown" not in got, "absent from LINCS is not the same as narrow"


def test_the_real_files_flag_the_broad_agents_when_they_are_present():
    root = Path(__file__).resolve().parents[1] / "data"
    up, down = root / "lincs_chem_pert_up.gmt", root / "lincs_chem_pert_down.gmt"
    if not (up.exists() and down.exists()):
        pytest.skip("the LINCS .gmt files are not in this checkout")
    b = sg.signature_breadth(up, down)
    got = sg.breadth_for(["vorinostat", "trichostatin A", "emetine"], b)
    assert len(got) == 3 and all(r["broad"] for r in got.values())
    assert 3000 < len(b["compounds"]) < 6000


# ------------------------------------------------------ the report cell
def test_the_cell_reads_plainly_in_each_case():
    cell = report._breadth_cell
    assert cell(None) == "—"
    assert cell({"thin": True}) == "few signatures"
    assert cell({"thin": False, "broad": False, "percentile": 41.0,
                 "median_genes": 98}) == "41th pct, 98 genes"
    assert cell({"thin": False, "broad": True, "percentile": 78.0,
                 "median_genes": 193}).startswith("⚠ broad (78th pct")


def _res_with_breadth(drug_stem_value):
    res = pipeline.run("SMN2", demo=True, offline=True)
    rows = report.corroboration(res)
    drug = rows[0]["drug"]
    res["compound_breadth"] = {
        "broad_percentile": 75, "min_signatures": 3, "n_signatures": 33132,
        "median_genes_all": 119,
        "compounds": {sg.breadth_key(drug): drug_stem_value}}
    return res, drug


def test_the_candidate_table_has_a_breadth_column_and_explains_it():
    res, drug = _res_with_breadth({
        "median_genes": 193, "n_signatures": 105, "percentile": 77.7,
        "thin": False, "broad": True})
    md = report.summary_report(res)
    assert "| Breadth |" in md
    assert "⚠ broad (78th pct, 193 genes)" in md
    assert "**Breadth** is how many genes the compound moved" in md
    assert "not a verdict" in md and "33,132" in md


def test_without_breadth_data_the_column_is_dashes_and_nothing_is_explained():
    res = pipeline.run("SMN2", demo=True, offline=True)
    md = report.summary_report(res)
    assert "| Breadth |" in md
    assert "**Breadth** is how many genes" not in md


def test_the_summary_marks_a_broad_compound_in_the_top_list():
    res, drug = _res_with_breadth({
        "median_genes": 300, "n_signatures": 50, "percentile": 95.0,
        "thin": False, "broad": True})
    routes = [r for b in res["results"]
              for r in (b["all_routes"] or {}).get("routes", [])]
    first = sorted(routes, key=lambda r: (r["tier"], -(r.get("score") or 0)))[0]
    res["compound_breadth"]["compounds"] = {
        sg.breadth_key(first["drug"]): res["compound_breadth"]["compounds"][
            sg.breadth_key(drug)]}
    text = report.summary_report(res)
    assert f"**{first['drug']}** ⚠broad" in text
    assert "A caution about specificity, not a score" in text


def test_a_breadth_flag_changes_no_rank_and_no_score():
    plain = pipeline.run("SMN2", demo=True, offline=True)
    flagged, drug = _res_with_breadth({
        "median_genes": 400, "n_signatures": 90, "percentile": 99.0,
        "thin": False, "broad": True})
    a = [(r["drug"], r["tier"], r["score"]) for b in plain["results"]
         for r in (b["all_routes"] or {}).get("routes", [])]
    b = [(r["drug"], r["tier"], r["score"]) for b in flagged["results"]
         for r in (b["all_routes"] or {}).get("routes", [])]
    assert a == b
    assert [r["drug"] for r in report.corroboration(plain)] \
        == [r["drug"] for r in report.corroboration(flagged)]


# --------------------------------------------- through a whole offline run
def test_a_run_with_lincs_files_records_breadth_for_its_candidates(
        tmp_path, monkeypatch):
    from test_upstream_search import _lf_offline_run, _lf_tf_edges

    up, down = _files(tmp_path)
    monkeypatch.setattr(sg, "DEFAULT_GMT_PATHS", dict(sg.DEFAULT_GMT_PATHS))
    res, rows = _lf_offline_run(
        tmp_path, monkeypatch, edges=_lf_tf_edges(),
        library=[{"chembl_id": "C1", "name": "x"}], index={},
        use_lincs=True, lincs_gmt=(str(up), str(down)), full=True)
    routes = {r["drug"].lower() for b in res["results"]
              for r in (b["all_routes"] or {}).get("routes", [])}
    assert {"alpha", "beta hydrochloride", "gamma"} <= routes, routes
    info = res["compound_breadth"]
    assert info["n_signatures"] == 12 and info["broad_percentile"] == 75
    comp = info["compounds"]
    assert comp["beta"]["broad"] is True and comp["beta"]["median_genes"] == 161
    assert comp["alpha"]["broad"] is False
    assert comp["gamma"]["thin"] is True
    md = report.summary_report(res)
    assert "⚠ broad" in md and "few signatures" in md


def test_a_run_without_lincs_files_has_an_empty_breadth_block(
        tmp_path, monkeypatch):
    from test_upstream_search import _lf_offline_run, _lf_tf_edges

    monkeypatch.setattr(sg, "DEFAULT_GMT_PATHS",
                        {"up": (str(tmp_path / "no_up.gmt"),),
                         "down": (str(tmp_path / "no_down.gmt"),)})
    res, rows = _lf_offline_run(
        tmp_path, monkeypatch, edges=_lf_tf_edges(),
        library=[{"chembl_id": "C1", "name": "x"}], index={}, full=True)
    assert res["compound_breadth"] == {}


# ------------------------------------------- very broad compounds lose score
def _extreme_files(tmp_path):
    """Thirty-two narrow filler signatures put `beta` (161 genes) in the top
    decile; `alpha` is narrow; `gamma` has a single, very large signature
    that would be very broad if one signature were enough."""
    def genes(prefix, n):
        return [f"{prefix}{i}" for i in range(n)]
    up, down = {}, {}
    for i in range(4):
        n = f"P{i} A549 24H-alpha-10.0"
        up[n], down[n] = genes("UA", 5), genes("DA", 5) + ["TESTG"]
    for i in range(3):
        n = f"P{i} MCF7 24H-beta hydrochloride-10.0"
        up[n], down[n] = genes("UB", 80), genes("DB", 80) + ["TESTG"]
    n = "P9 PC3 6H-gamma-1.11"
    up[n], down[n] = genes("UG", 200), genes("DG", 200) + ["TESTG"]
    for i in range(32):
        n = f"P{i} HA1E 24H-filler{i}-1.0"
        up[n], down[n] = genes("UF", 10), genes("DF", 10)
    return _gmt(tmp_path / "up.gmt", up), _gmt(tmp_path / "down.gmt", down)


def test_only_the_top_decile_with_enough_signatures_is_very_broad(tmp_path):
    up, down = _extreme_files(tmp_path)
    b = sg.signature_breadth(up, down)
    got = sg.breadth_for(["beta", "alpha", "gamma"], b)
    assert got["beta"]["percentile"] >= sg.VERY_BROAD_PERCENTILE == 90
    assert got["beta"]["broad"] and got["beta"]["very_broad"]
    assert not got["alpha"]["very_broad"]
    assert got["gamma"]["thin"] and not got["gamma"]["very_broad"], (
        "one signature is never enough to cost a route score")


def test_the_quartile_that_is_only_broad_keeps_its_score(tmp_path):
    """The threshold that flags and the one that scores down are different."""
    up, down = _files(tmp_path)               # beta sits at the 75th percentile
    got = sg.breadth_for(["beta"], sg.signature_breadth(up, down))["beta"]
    assert got["broad"] is True and got["very_broad"] is False


def _block():
    res = pipeline.run("SMN2", demo=True, offline=True)
    return res["results"][0]


def test_a_very_broad_compound_has_its_score_halved_inside_its_tier():
    from generoutes import catalogue as cat
    block = _block()
    base = cat.build(block, goal="up")
    target = next(r for r in base["routes"] if r["score"] > 0)
    rec = {"median_genes": 300, "n_signatures": 80, "percentile": 96.0,
           "thin": False, "broad": True, "very_broad": True}
    pen = cat.build(block, goal="up",
                    breadth={sg.breadth_key(target["drug"]): rec})
    after = [r for r in pen["routes"] if r["drug"] == target["drug"]]
    before = [r for r in base["routes"] if r["drug"] == target["drug"]]
    assert len(after) == len(before) >= 1
    for a, b in zip(sorted(after, key=lambda r: r["acts_on"] or ""),
                    sorted(before, key=lambda r: r["acts_on"] or "")):
        assert a["score"] == round(b["score"] * cat.BROAD_SCORE_FACTOR, 1)
        assert a["breadth_penalty"] == cat.BROAD_SCORE_FACTOR == 0.5
        assert a["tier"] == b["tier"], "the penalty never moves a route's tier"
        assert any("score multiplied by 0.5" in n for n in a["notes"])
    others_a = {(r["drug"], r["acts_on"], r["category"]): r["score"]
                for r in pen["routes"] if r["drug"] != target["drug"]}
    others_b = {(r["drug"], r["acts_on"], r["category"]): r["score"]
                for r in base["routes"] if r["drug"] != target["drug"]}
    assert others_a == others_b, "no other route's score changes"
    tiers = [r["tier"] for r in pen["routes"]]
    assert tiers == sorted(tiers), "still ordered tier first"


def test_without_breadth_nothing_is_penalised_and_scores_are_unchanged():
    from generoutes import catalogue as cat
    block = _block()
    plain = cat.build(block, goal="up")
    again = cat.build(block, goal="up", breadth={})
    assert [r["score"] for r in plain["routes"]] \
        == [r["score"] for r in again["routes"]]
    assert all(r["breadth_penalty"] is None for r in plain["routes"])


def test_a_whole_run_penalises_the_very_broad_compound_and_only_that_one(
        tmp_path, monkeypatch):
    from test_upstream_search import _lf_offline_run, _lf_tf_edges

    up, down = _extreme_files(tmp_path)
    monkeypatch.setattr(sg, "DEFAULT_GMT_PATHS", dict(sg.DEFAULT_GMT_PATHS))
    kw = dict(edges=_lf_tf_edges(), library=[{"chembl_id": "C1", "name": "x"}],
              index={}, use_lincs=True, lincs_gmt=(str(up), str(down)),
              full=True)
    pen, _ = _lf_offline_run(tmp_path, monkeypatch, **kw)
    # the same run with nothing counted as very broad
    monkeypatch.setattr(sg, "VERY_BROAD_PERCENTILE", 101)
    sg._BREADTH_CACHE.clear()
    base, _ = _lf_offline_run(tmp_path, monkeypatch, **kw)

    def scores(res):
        out = {}
        for b in res["results"]:
            for r in (b["all_routes"] or {}).get("routes", []):
                out[(r["drug"].lower(), r["category"])] = r
        return out
    sp, sb = scores(pen), scores(base)
    beta = [k for k in sp if k[0].startswith("beta")]
    assert beta, sp.keys()
    for k in beta:
        assert sp[k]["breadth_penalty"] == 0.5
        assert sp[k]["score"] == round(sb[k]["score"] * 0.5, 1)
        assert sb[k]["breadth_penalty"] is None
    for k in sp:
        if k not in beta:
            assert sp[k]["score"] == sb[k]["score"] and not sp[k]["breadth_penalty"]
    assert pen["compound_breadth"]["compounds"]["beta"]["very_broad"] is True


def test_a_compound_that_is_only_broad_is_not_penalised_by_the_catalogue():
    from generoutes import catalogue as cat
    block = _block()
    base = cat.build(block, goal="up")
    target = next(r for r in base["routes"] if r["score"] > 0)
    rec = {"median_genes": 190, "n_signatures": 100, "percentile": 78.0,
           "thin": False, "broad": True, "very_broad": False}
    out = cat.build(block, goal="up",
                    breadth={sg.breadth_key(target["drug"]): rec})
    assert [r["score"] for r in out["routes"]]         == [r["score"] for r in base["routes"]]
    assert all(r["breadth_penalty"] is None for r in out["routes"])
