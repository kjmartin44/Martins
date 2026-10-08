"""Approved treatments for the condition a gene causes, shown beside the
ranking and never in it.
"""

from __future__ import annotations

import json

import pytest

from promoter_tf_drugs import pipeline, report
from promoter_tf_drugs import standard_of_care as soc


def test_huntington_disease_lists_the_three_approved_chorea_drugs():
    got = soc.for_gene("HTT")
    assert got["condition"] == "Huntington disease" and got["gene"] == "HTT"
    by = {t["generic"]: t for t in got["treatments"]}
    assert set(by) == {"tetrabenazine", "deutetrabenazine", "valbenazine"}
    assert by["tetrabenazine"]["brand"] == "Xenazine"
    assert by["tetrabenazine"]["approved"] == "2008-08"
    assert by["deutetrabenazine"]["brand"] == "Austedo"
    assert by["deutetrabenazine"]["approved"] == "2017-04"
    assert by["valbenazine"]["brand"] == "Ingrezza"
    assert by["valbenazine"]["approved"] == "2023-08"
    json.dumps(got)


def test_none_of_them_is_recorded_as_acting_on_the_gene():
    """They block VMAT2 (SLC18A2) and treat the movement symptom. Listing
    them as drugs for HTT would be wrong, and the table says so."""
    got = soc.for_gene("htt")
    assert all(t["target_gene"] == "SLC18A2" for t in got["treatments"])
    assert not any(t["acts_on_gene"] for t in got["treatments"])
    assert "None of these acts on HTT" in got["note"]
    assert "SLC18A2" in got["note"]


def test_a_gene_with_no_entry_gets_none_not_an_empty_claim():
    assert soc.for_gene("CCR5") is None
    assert soc.for_gene("") is None and soc.for_gene(None) is None


def test_a_treatment_that_does_act_on_the_gene_is_marked_so(monkeypatch):
    cond = soc.Condition(
        name="Example", genes=("GENEX",), source="test",
        treatments=(soc.Treatment("drugone", "Brand", "2020-01", "use",
                                  "GENEX", "inhibitor"),))
    monkeypatch.setattr(soc, "CONDITIONS", (cond,))
    got = soc.for_gene("genex")
    assert got["treatments"][0]["acts_on_gene"] is True
    assert got["note"] == "drugone acts on GENEX itself."


def test_the_report_lists_them_apart_from_the_ranking():
    res = pipeline.run("HTT", demo=True, offline=True, goal="down")
    assert res["standard_of_care"]["condition"] == "Huntington disease"
    md = report.markdown_report(res)
    assert "## Approved treatments for Huntington disease" in md
    for name in ("Xenazine", "Austedo", "Ingrezza", "2023-08", "SLC18A2"):
        assert name in md
    assert "They are context, not candidates" in md
    assert "none of them lowers HTT" in md
    head, section = md.split("## Approved treatments for Huntington disease")
    assert "Ingrezza" not in head, "not mixed into anything before it"


def test_the_summary_says_it_in_one_line():
    res = pipeline.run("HTT", demo=True, offline=True, goal="down")
    text = report.summary_report(res)
    assert "**Already approved for Huntington disease:** tetrabenazine, " \
        "deutetrabenazine, valbenazine." in text


def test_a_gene_without_an_entry_gets_no_section():
    res = pipeline.run("SMN2", demo=True, offline=True)
    assert res["standard_of_care"] is None
    assert "Approved treatments for" not in report.markdown_report(res)
    assert "Already approved for" not in report.summary_report(res)


def test_they_are_not_in_any_ranked_route():
    res = pipeline.run("HTT", demo=True, offline=True, goal="down")
    names = {r["drug"].lower() for b in res["results"]
             for r in (b.get("all_routes") or {}).get("routes") or []}
    assert not names & {"tetrabenazine", "deutetrabenazine", "valbenazine"}
