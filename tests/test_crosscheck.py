"""CTD against LINCS and against the routes, and the mitotic caution.

A route predicts that a drug should move the gene the wanted way. CTD records
which way someone measured it, and LINCS which way a set of cell lines did. These
tests pin that the three are compared and shown beside the drug, that nothing is
demoted for the comparison, and that a spindle poison on a mitotic gene is said
to be a different case from a switch at the promoter.
"""

from __future__ import annotations

import pytest

from promoter_tf_drugs import catalogue as cat
from promoter_tf_drugs import crosscheck as xc_mod
from promoter_tf_drugs import report, systemic
from promoter_tf_drugs import retrospective as rt


def _links(wanted=(), opposite=()):
    def side(names, d):
        return [{"chemical": n, "direction": d, "kind": "chemical",
                 "actions": "", "in_candidates": False} for n in names]
    return {"reported_links": {"wanted_direction": side(wanted, "up"),
                               "opposite_direction": side(opposite, "down")}}


@pytest.fixture
def lincs(monkeypatch):
    def install(wanted=(), opposite=(), reason=None):
        def fake(gene, goal, up_path=None, down_path=None):
            if reason:
                return None, reason
            return rt.KeySet("lincs", wanted=set(wanted),
                             opposite=set(opposite)), None
        monkeypatch.setattr(rt, "lincs_keyset", fake)
    return install


def test_no_ctd_record_means_no_crosscheck():
    assert xc_mod.build("G1", "up", None) is None
    assert xc_mod.build("G1", "up", {}) is None
    assert xc_mod.build("G1", "sideways", _links(["a"])) is None


def test_a_drug_ctd_records_the_other_way_is_listed(lincs):
    lincs()
    xc = xc_mod.build("G1", "up", _links(["Aspirin"], ["Tretinoin"]))
    assert xc["ctd_opposite"] == ["Tretinoin"]
    assert xc["ctd_wanted"] == ["Aspirin"]


def test_exposures_and_headings_are_not_drugs(lincs):
    lincs()
    imp = _links()
    imp["reported_links"]["opposite_direction"] = [
        {"chemical": "Formaldehyde", "kind": "exposure"},
        {"chemical": "Anti-HIV Agents", "kind": "category"}]
    assert xc_mod.build("G1", "up", imp)["ctd_opposite"] == []


def test_ctd_and_lincs_that_differ_are_a_conflict_and_the_rest_agree(lincs):
    lincs(wanted={"aspirin"}, opposite={"tretinoin", "caffeine"})
    xc = xc_mod.build("G1", "up", _links(["Aspirin", "Tretinoin"],
                                         ["Caffeine"]))
    assert xc["n_compared"] == 3 and xc["n_agree"] == 2
    assert xc["conflicts"] == [{"drug": "Tretinoin", "ctd": "up",
                                "lincs": "down"}]


def test_a_down_goal_reports_directions_as_up_and_down(lincs):
    lincs(wanted={"x"})
    xc = xc_mod.build("G1", "down", _links([], ["X"]))
    assert xc["conflicts"] == [{"drug": "X", "ctd": "up", "lincs": "down"}]


def test_a_drug_lincs_does_not_call_is_not_compared(lincs):
    lincs(wanted={"other"})
    xc = xc_mod.build("G1", "up", _links(["Aspirin"]))
    assert xc["n_compared"] == 0 and xc["conflicts"] == []


def test_missing_lincs_is_said_and_is_not_agreement(lincs):
    lincs(reason="no local LINCS .gmt files")
    xc = xc_mod.build("G1", "up", _links(["a"], ["b"]))
    assert xc["lincs"]["checked"] is False
    assert xc["lincs"]["reason"] == "no local LINCS .gmt files"
    assert xc["ctd_opposite"] == ["b"], "the CTD half is still returned"
    text = "\n".join(xc_mod.markdown(xc))
    assert "LINCS was not compared with CTD: no local LINCS .gmt files" in text
    off = xc_mod.build("G1", "up", _links(["a"]), use_lincs=False)
    assert off["lincs"]["reason"] == "LINCS was not requested for this run"


def test_a_route_for_a_drug_ctd_records_the_other_way_carries_a_caution(lincs):
    lincs()
    xc = xc_mod.build("G1", "up", _links([], ["Tretinoin"]))
    known = {"drug": "Tretinoin Sodium", "category": "TF: orthosteric ligand",
             "notes": []}
    assert xc_mod.cautions_for(known, xc) == [
        "CTD records the opposite direction"]
    unknown = dict(known, notes=["which way X moves RXRA is not recorded, so "
                                 "it may push the gene the wrong way"])
    assert "own direction is unknown" in xc_mod.cautions_for(unknown, xc)[0]
    assert xc_mod.cautions_for(dict(known, drug="Aspirin"), xc) == []
    literature = dict(known, category="literature: reported effect")
    assert xc_mod.cautions_for(literature, xc) == [], "that is CTD itself"
    assert xc_mod.cautions_for(known, None) == []


def test_the_catalogue_notes_the_caution_and_changes_no_tier_or_score(lincs):
    lincs()
    xc = xc_mod.build("G1", "up", _links([], ["Tretinoin"]))
    block = {"lincs": {"hits": [{"drug": "Tretinoin", "direction": "up",
                                 "score": 2.0, "is_landmark": True,
                                 "feature_space": "landmark",
                                 "cell_line": "HA1E"}]}}
    plain = cat.build(block, goal="up", include_chromatin=False)["routes"][0]
    marked = cat.build(dict(block, crosscheck=xc), goal="up",
                       include_chromatin=False)["routes"][0]
    assert any("CTD records the opposite direction" in n
               for n in marked["notes"])
    assert (marked["tier"], marked["score"]) == (plain["tier"], plain["score"])


def _summary_result(routes, xc):
    return {"query": {"gene": "G1", "goal": "up"},
            "results": [{"all_routes": {"routes": routes}, "crosscheck": xc}]}


def _route(drug, tier=2, notes=()):
    return {"drug": drug, "tier": tier, "score": 50.0 - tier, "source": "t",
            "category": "TF: orthosteric ligand", "acts_on": "RXRA",
            "mechanism": "m", "notes": list(notes),
            "axes": {"evidence": "measured binding"}}


def test_the_summary_marks_the_drug_and_names_it_under_the_table(lincs):
    lincs(wanted={"aspirin"}, opposite={"tretinoin"})
    xc = xc_mod.build("G1", "up", _links(["Tretinoin"], ["Aspirin"]))
    text = report.summary_report(_summary_result(
        [_route("Aspirin"), _route("Caffeine", 3)], xc))
    rows = [ln for ln in text.splitlines()
            if ln.startswith("| ") and "**" in ln]
    assert any("Aspirin" in r and "⚠CTD opposite" in r for r in rows)
    assert not any("Caffeine" in r and "⚠CTD opposite" in r for r in rows)
    assert "CTD records Aspirin moving G1 the other way" in text
    assert "CTD and LINCS disagree on 2 drugs" in text
    assert "Tretinoin" in text.split("CTD and LINCS disagree")[1]


def test_a_summary_without_a_crosscheck_is_unchanged():
    text = report.summary_report(_summary_result([_route("Aspirin")], None))
    assert "⚠CTD opposite" not in text and "disagree" not in text


def test_the_full_report_section_lists_flagged_drugs_and_conflicts(lincs):
    lincs(wanted={"aspirin"}, opposite={"tretinoin"})
    xc = xc_mod.build("G1", "up", _links(["Tretinoin"], ["Aspirin"]))
    lines = xc_mod.markdown(xc, [_route("Aspirin", 2, [
        "which way A moves B is not recorded"])])
    text = "\n".join(lines)
    assert "## Where the sources disagree" in text
    assert ("| Aspirin | 2 | RXRA | CTD records the opposite direction and "
            "this route's own direction is unknown |") in text
    assert "| Tretinoin | up | down |" in text
    assert "Nothing is demoted" in text


def test_no_disagreement_and_nothing_to_flag_writes_no_section(lincs):
    lincs()
    xc = xc_mod.build("G1", "up", _links(["a"]))
    assert xc_mod.markdown(xc, [_route("a")]) == []


# ------------------------------------------------------------- the mitotic case
def test_spindle_poisons_and_the_ctd_protein_compounds_count_as_mitotic():
    for name in ("Paclitaxel", "Benomyl", "Rotenone",
                 "2'-Fluoro-6,7-Methylenedioxy-2-Phenyl-4-Quinolone"):
        assert systemic.is_mitotic_disruptor(name), name
    assert not systemic.is_mitotic_disruptor("Estradiol")


def _ctd_block(name):
    return {"literature": {"gene": "BUB1B", "papers": {name.lower(): 1},
                           "reported_links": {"wanted_direction": [{
                               "chemical": name, "kind": "chemical",
                               "direction": "up", "actions": ""}]}}}


def test_a_ctd_report_for_a_spindle_poison_on_a_mitotic_gene_says_so():
    got = cat.build(_ctd_block("Paclitaxel"), goal="up", mitotic_gene=True,
                    include_chromatin=False)["routes"][0]
    assert got["mitotic_penalty"] == cat.MITOTIC_SCORE_FACTOR
    assert any("cells held in mitosis" in n for n in got["notes"])


def test_the_mitotic_caution_is_only_for_mitotic_genes_and_disruptors():
    other_gene = cat.build(_ctd_block("Paclitaxel"), goal="up",
                           include_chromatin=False)["routes"][0]
    assert not any("held in mitosis" in n for n in other_gene["notes"])
    other_drug = cat.build(_ctd_block("Estradiol"), goal="up",
                           mitotic_gene=True,
                           include_chromatin=False)["routes"][0]
    assert not any("held in mitosis" in n for n in other_drug["notes"])
