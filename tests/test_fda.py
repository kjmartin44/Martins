"""Reports list only FDA approved drugs, and results.json keeps the rest.

The rule is narrow and has a stated cost, so these tests pin both: what counts
as approved, that a missing Orange Book never empties a report, and that the
full data is untouched.
"""

from __future__ import annotations

import copy
import json

import pytest

from promoter_tf_drugs import cli, combined, exposure, fda, orangebook, pipeline
from promoter_tf_drugs import report

APPROVED = {"calcitriol", "tretinoin", "bexarotene"}
_RECORD = {"ingredient": "x", "reach": "systemic", "routes": ["ORAL"],
           "reach_note": "approved systemically (oral)",
           "generic_available": True, "forms": ["TABLET"], "n_products": 1,
           "n_applicants": 1, "trade_names": [], "systemic_routes": ["ORAL"],
           "local_routes": [], "unclassified_routes": [],
           "combination_only": False, "fda_approved": True}


@pytest.fixture
def table(monkeypatch):
    """An Orange Book that knows three drugs, and a filter that is switched on."""
    fda.reset()
    monkeypatch.setattr(fda, "available", lambda: True)
    monkeypatch.setattr(orangebook, "ingredients",
                        lambda: ["tretinoin", "bexarotene", "calcitriol",
                                 "lapatinib ditosylate"])
    monkeypatch.setattr(
        orangebook, "lookup",
        lambda name, path=None: (_RECORD if name.lower() in APPROVED
                                 else None))
    yield
    fda.reset()


def _demo_applied():
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["fda_filter"] = {"requested": True, "applied": True, "note": None}
    return res


# ------------------------------------------------------------------- rule
def test_a_drug_is_approved_when_the_orange_book_lists_it(table):
    assert fda.is_approved("Tretinoin") and fda.is_approved("  BEXAROTENE ")
    assert not fda.is_approved("Pemafibrate"), "approved only outside the US"
    assert not fda.is_approved("roscovitine") and not fda.is_approved("")


def test_a_bare_name_finds_the_salt_the_orange_book_spells_out(table):
    assert fda.is_approved("Lapatinib") and fda.is_approved("lapatinib")
    assert fda.base_name("Rezafungin Acetate") == "rezafungin"
    assert fda.base_name("Sodium Fluoride") == "sodium fluoride",         "a name made only of ions is the drug itself"
    assert not fda.is_approved("Lapatinibx")


def test_a_curated_gene_therapy_counts_though_the_orange_book_has_no_biologics(
        table):
    assert fda.is_approved("Voretigene neparvovec")
    assert fda.is_approved("Onasemnogene abeparvovec")
    assert not fda.is_approved("Ivacaftor"), \
        "a small molecule has to be in the Orange Book"


def test_the_status_says_why_the_filter_is_or_is_not_applied(monkeypatch):
    monkeypatch.setattr(fda, "available", lambda: True)
    assert fda.status(demo=False)["applied"] is True
    demo = fda.status(demo=True)
    assert not demo["applied"] and "synthetic" in demo["note"]
    off = fda.status(demo=False, requested=False)
    assert not off["applied"] and not off["requested"]
    monkeypatch.setattr(fda, "available", lambda: False)
    gone = fda.status(demo=False)
    assert not gone["applied"] and "no Orange Book" in gone["note"]


def test_a_missing_orange_book_never_empties_a_report(monkeypatch):
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    res["fda_filter"] = fda.status(demo=False)       # the autouse fixture: none
    assert not res["fda_filter"]["applied"]
    assert fda.apply(res) is res
    md = report.markdown_report(res)
    assert "FDA filter not applied" in md and "Tretinoin" in md


# ------------------------------------------------------------------- view
def test_the_view_drops_unapproved_routes_and_rebuilds_the_counts(table):
    res = _demo_applied()
    before = copy.deepcopy(res)
    view = fda.apply(res)
    names = {r["drug"] for r in view["results"][0]["all_routes"]["routes"]}
    assert names and all(n.lower() in APPROVED for n in names)
    ar = view["results"][0]["all_routes"]
    assert ar["n_routes"] == len(ar["routes"])
    assert sum(ar["by_category"].values()) == ar["n_routes"]
    assert set(ar["best_per_category"]) == set(ar["by_category"])
    assert res == before, "the original result is not changed"
    h = view["fda_filter"]["hidden"]
    assert h["entries"] > 0 and "Pioglitazone" in h["drugs"]


def test_the_view_reaches_every_list_of_drugs(table):
    res = _demo_applied()
    blk = res["results"][0]
    blk["lincs"] = {"hits": [{"drug": "Tretinoin"}, {"drug": "BX-795"}]}
    blk["screens"] = {"RXRA": {"hits": [{"drug": "Bexarotene"},
                                        {"drug": "Pemafibrate"}]}}
    blk["improvements"] = {"reported_links": {"wanted_direction": [
        {"chemical": "Calcitriol"}, {"chemical": "Bisphenol A"}]}}
    blk["upstream"] = {"CREB1": [{"drugs": [{"name": "Tretinoin"},
                                            {"name": "Forskolin"}]}]}
    res["direct_therapies"] = [{"drug": "Voretigene neparvovec"},
                               {"drug": "Ivacaftor"}]
    view = fda.apply(res)
    v = view["results"][0]
    assert [h["drug"] for h in v["lincs"]["hits"]] == ["Tretinoin"]
    assert [h["drug"] for h in v["screens"]["RXRA"]["hits"]] == ["Bexarotene"]
    assert [h["chemical"] for h in
            v["improvements"]["reported_links"]["wanted_direction"]] == [
        "Calcitriol"]
    assert [d["name"] for d in v["upstream"]["CREB1"][0]["drugs"]] == [
        "Tretinoin"]
    assert [t["drug"] for t in view["direct_therapies"]] == [
        "Voretigene neparvovec"]
    for t in v["tfs"]:
        assert all(d["name"].lower() in APPROVED for d in t["drugs"])
    # a transcription factor's own name is never read as a drug
    assert [t["name"] for t in v["tfs"]] == [t["name"] for t in blk["tfs"]]


def test_a_result_that_did_not_apply_the_filter_is_returned_as_it_is(table):
    res = pipeline.run("DEMOGENE", demo=True, offline=True)
    assert fda.apply(res) is res
    assert fda.apply({"results": []}) == {"results": []}


# ---------------------------------------------------------------- reports
def test_the_report_and_summary_omit_unapproved_drugs_and_say_so(table):
    res = _demo_applied()
    md, summ = report.markdown_report(res), report.summary_report(res)
    for text in (md, summ):
        assert "Only FDA approved drugs are shown" in text
        assert "Pioglitazone" not in text and "Rosiglitazone" not in text
        assert "Tretinoin" in text
    assert "Pioglitazone" in json.dumps(res), "results.json keeps them"


def test_write_outputs_keeps_everything_in_results_json(table, tmp_path):
    res = _demo_applied()
    pipeline.write_outputs(res, tmp_path)
    saved = json.loads((tmp_path / "results.json").read_text(encoding="utf-8"))
    assert "Pioglitazone" in json.dumps(saved)
    assert "Pioglitazone" not in (tmp_path / "report.md").read_text(
        encoding="utf-8")
    assert "Pioglitazone" not in (tmp_path / "tf_summary.tsv").read_text(
        encoding="utf-8")
    assert saved["fda_filter"]["applied"] is True


def test_run_records_the_filter_and_the_flag_reaches_it(monkeypatch):
    monkeypatch.setattr(fda, "available", lambda: True)
    got = pipeline.run("DEMOGENE", demo=True, offline=True, fda_only=False)
    assert got["query"]["fda_only"] is False
    assert got["fda_filter"]["requested"] is False


def test_the_command_line_flag_turns_the_filter_off():
    parser = cli.build_parser()
    on = cli.run_kwargs_from_args(parser.parse_args(["--gene", "G"]))
    off = cli.run_kwargs_from_args(
        parser.parse_args(["--gene", "G", "--include-non-fda"]))
    assert on["fda_only"] is True and off["fda_only"] is False


# --------------------------------------------------- the other readers
def test_the_dose_check_can_keep_only_approved_drugs(table):
    block = {"all_routes": {"routes": [
        {"drug": "Tretinoin", "tier": 2, "score": 5, "approved": True},
        {"drug": "Pemafibrate", "tier": 1, "score": 9, "approved": True},
        {"drug": "roscovitine", "tier": 1, "score": 8, "approved": False}]}}
    assert [n for _, n in exposure.candidate_drugs(block, 10, True)] == [
        "Tretinoin"]
    assert len(exposure.candidate_drugs(block, 10)) == 3, "default unchanged"
    assert list(exposure.assess(block, fda_only=True)["drugs"]) == ["Tretinoin"]


def test_docking_hits_are_dropped_before_they_take_a_rank(table):
    doc = {"usable": True, "hits": [{"drug": "Pemafibrate"},
                                    {"drug": "Tretinoin"},
                                    {"drug": "Bexarotene"}]}
    src = combined.docking_source(doc, "RXRA", keep=fda.is_approved)
    assert src["ranks"] == {"tretinoin": 1, "bexarotene": 2}
    assert len(combined.docking_source(doc, "RXRA")["ranks"]) == 3


def test_the_combined_list_uses_the_filter_the_run_recorded(table):
    res = _demo_applied()
    got = combined.build(results=res, dockings=[("RXRA", {
        "usable": True, "hits": [{"drug": "Pemafibrate"},
                                 {"drug": "Bexarotene"}]})], gene="DEMOGENE")
    text = json.dumps(got["table"])
    assert "Bexarotene" in text and "Tretinoin" in text
    for gone in ("Pemafibrate", "Pioglitazone", "Rosiglitazone", "Prednisolone"):
        assert gone not in text, gone
