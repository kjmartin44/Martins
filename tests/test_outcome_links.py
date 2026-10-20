"""The outcome chain is reported as the links that come from evidence.

Direction defaults to the database label and responsiveness is a placeholder,
so each is the same for every row. A constant cannot separate candidates, and
multiplying it in only made a small number that read like a measurement. These
tests pin that the two evidence links are shown, the two constants are said once,
and a link that some row did measure comes back as a column.
"""

from __future__ import annotations

import copy

from promoter_tf_drugs import outcome, pipeline, report


def est(drug="D", tf="T", engage=0.8, direction=0.75, activate=0.6,
        responsive=0.4, direction_measured=False, responsive_measured=False):
    notes = []
    if not direction_measured:
        notes.append("direction from database annotation only, not "
                     "independently predicted")
    if not responsive_measured:
        notes.append("gene responsiveness in your cell type is unknown and "
                     "assumed 0.40; this is the link the package cannot see")
    return {"drug": drug, "tf_symbol": tf, "p_engages": engage,
            "p_direction": direction, "p_activates": activate,
            "p_responsive": responsive, "p_point": 0.1, "p_low": 0.05,
            "p_high": 0.2, "limiting_step": "gene responsiveness in your "
            "cells", "tissue_dependent": False, "verdict": "speculative",
            "next_experiment": "x", "notes": notes}


def test_the_two_constant_links_are_assumed_and_the_other_two_are_shown():
    ests = [est("A", engage=0.9), est("B", engage=0.5)]
    assert outcome.assumed_links(ests) == ["p_direction", "p_responsive"]
    assert outcome.shown_links(ests) == ["p_engages", "p_activates"]


def test_one_row_is_enough_to_know_which_links_were_assumed():
    """Constancy is read from how the value was made, not from how many rows
    there are, so a single estimate is not called constant by accident."""
    assert outcome.shown_links([est()]) == ["p_engages", "p_activates"]
    assert outcome.assumed_links([]) == []


def test_a_measured_direction_becomes_a_column():
    ests = [est("A", direction=0.9, direction_measured=True),
            est("B", direction=0.2, direction_measured=True)]
    assert outcome.assumed_links(ests) == ["p_responsive"]
    assert outcome.shown_links(ests) == ["p_engages", "p_direction",
                                         "p_activates"]


def test_a_link_measured_for_only_some_rows_is_shown_not_assumed():
    ests = [est("A"), est("B", direction=0.9, direction_measured=True)]
    assert "p_direction" not in outcome.assumed_links(ests)
    assert "p_direction" in outcome.shown_links(ests)


def test_the_weaker_link_is_the_lowest_of_the_links_shown():
    e = est(engage=0.9, activate=0.5)
    keys = outcome.shown_links([e])
    assert outcome.weaker_link(e, keys) == "promoter activation"
    e2 = est(engage=0.4, activate=0.5)
    assert outcome.weaker_link(e2, keys) == "drug engagement"
    # an assumed constant is never named as the weaker link
    assert outcome.weaker_link(e2, keys) != "responsiveness"


def test_nothing_is_said_about_links_that_were_not_measured():
    """A constant is left out of the tables and not explained at length."""
    assert not hasattr(outcome, "assumed_caveat")
    assert outcome.knockout_caveat([est()]) == []


def test_experiments_exist_for_every_link():
    for key in outcome.STEP_OF_LINK:
        assert outcome.experiment_for(key)


def test_the_report_shows_a_measured_direction_as_a_column():
    res = pipeline.run("SMN2", demo=True, offline=True)
    blk = res["results"][0]
    for i, e in enumerate(blk["outcome_estimates"]):
        e["notes"] = [n for n in e["notes"] if "annotation only" not in n]
        e["p_direction"] = 0.3 + 0.1 * (i % 5)
    md = report.markdown_report(res)
    start = md.index("Strength of evidence for each factor route")
    section = md[start:md.index("\n### ", start + 10)]
    header = [l for l in section.splitlines() if "Drug engagement" in l][0]
    assert "Direction" in header and "Responsiveness" not in header
    assert "Responsiveness" not in section


def test_saved_fields_are_unchanged():
    """results.json keeps the whole chain for anything that reads it."""
    res = pipeline.run("SMN2", demo=True, offline=True)
    e = res["results"][0]["outcome_estimates"][0]
    for key in ("p_engages", "p_direction", "p_activates", "p_responsive",
                "p_point", "p_low", "p_high", "limiting_step", "verdict",
                "next_experiment", "notes"):
        assert key in e


def test_no_report_prints_a_chance_that_the_gene_moves():
    res = pipeline.run("SMN2", demo=True, offline=True)
    before = copy.deepcopy(res)
    for text in (report.markdown_report(res), report.summary_report(res)):
        for phrase in ("P(gene moves)", "P(moves)", "P(responsive)",
                       "magnitude of uncertainty"):
            assert phrase not in text, phrase
    assert res == before, "rendering does not change the result"


# ------------------------------------------- a model's call is shown, if any
class _Pred:
    predicted = "agonist"
    probability_agonist = 0.9
    confidence = "high"
    protrusion: dict = {}


def test_a_model_call_replaces_the_database_label_for_that_row():
    given = outcome.estimate("D", "T", "G", prediction=_Pred())
    plain = outcome.estimate("D", "T", "G")
    assert plain.p_direction == 0.75
    assert any("annotation only" in n for n in plain.notes)
    assert given.p_direction > 0.75
    assert not any("annotation only" in n for n in given.notes)


def test_an_undetermined_call_leaves_the_label_in_place():
    class Undetermined(_Pred):
        predicted = "undetermined"
    got = outcome.estimate("D", "T", "G", prediction=Undetermined())
    assert got.p_direction == 0.75
    assert any("annotation only" in n for n in got.notes)


def test_the_summary_table_shows_direction_only_when_a_model_made_a_call():
    res = pipeline.run("SMN2", demo=True, offline=True)
    lines = pipeline.corroboration_table(res)
    header = [l for l in lines if l.startswith("| Drug | Methods")][0]
    assert "Direction" not in header
    for e in res["results"][0]["outcome_estimates"]:
        e["notes"] = [n for n in e["notes"] if "annotation only" not in n]
        e["p_direction"] = 0.9
    lines = pipeline.corroboration_table(res)
    header = [l for l in lines if l.startswith("| Drug | Methods")][0]
    assert "| Direction |" in header
    sep = lines[lines.index(header) + 1]
    body = [l for l in lines[lines.index(header) + 2:] if l.startswith("| **")]
    assert body and sep.count("---") == header.count("|") - 1 == 13
    assert all(b.count("|") == header.count("|") for b in body)
    assert any("| 0.90 |" in b for b in body)
