"""Responsiveness read from knockout experiments instead of a fixed 0.40.

The data are lists of the genes that changed most when a factor was lost or
raised. What matters is how a list is read: which way a gene moved, how a gain
differs from a loss, what a mutation or a missing factor does NOT say, and that
absence from a list is weak evidence.
"""

from __future__ import annotations

import pytest

from generoutes import knockout as ko
from generoutes import outcome, pipeline, report


def gmt(tmp_path, name, lines):
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(p)


def row(term, *genes):
    return term + "\t\t" + "\t".join(genes)


# ----------------------------------------------------------------- parsing
@pytest.mark.parametrize("term,expected", [
    ("FOXA1 Up", ("FOXA1", "FOXA1", "loss", "up", "LINCS CRISPR KO")),
    ("FOXA1 Down", ("FOXA1", "FOXA1", "loss", "down", "LINCS CRISPR KO")),
    ("NR3C1 KO MOUSE GSE1 CREEDSID GENE 5 UP",
     ("NR3C1", "NR3C1 KO MOUSE GSE1 CREEDSID GENE 5", "loss", "up",
      "GEO TF perturbation")),
    ("MYC OE HUMAN GSE2 CREEDSID GENE 7 DOWN",
     ("MYC", "MYC OE HUMAN GSE2 CREEDSID GENE 7", "gain", "down",
      "GEO TF perturbation")),
    ("STAT3 SIRNA HELA GSE3 RNASEQ UP",
     ("STAT3", "STAT3 SIRNA HELA GSE3 RNASEQ", "loss", "up",
      "GEO TF perturbation")),
])
def test_a_library_term_is_read_as_a_loss_or_a_gain(term, expected):
    assert ko._parse(term) == expected


@pytest.mark.parametrize("term", [
    "TP53 MUT HUMAN GSE1 CREEDSID GENE 1 UP",       # may raise or lower
    "ESR1 TRUNCATED MOUSE GSE1 CREEDSID GENE 1 UP",
    "EWS EWSFUSION HUMAN GSE1 CREEDSID GENE 1 UP",
    "no direction word",
    "",
])
def test_a_mutation_a_fusion_or_a_term_with_no_direction_is_not_read(term):
    assert ko._parse(term) is None


def test_loss_down_and_gain_up_both_mean_the_factor_activates():
    assert ko._vote("loss", "down") == "activates"
    assert ko._vote("gain", "up") == "activates"
    assert ko._vote("loss", "up") == "represses"
    assert ko._vote("gain", "down") == "represses"


# --------------------------------------------------------------- assessing
@pytest.fixture
def files(tmp_path):
    crispr = gmt(tmp_path, "crispr.gmt", [
        row("AAA Up", "G1", "X1"), row("AAA Down", "X2"),
        row("BBB Up", "X1"), row("BBB Down", "X2"),
        row("CCC Up", "G1"), row("CCC Down", "G1")])
    tfp = gmt(tmp_path, "tfp.gmt", [
        row("AAA KO MOUSE GSE1 CREEDSID GENE 1 UP", "X1"),
        row("AAA KO MOUSE GSE1 CREEDSID GENE 1 DOWN", "G1"),
        row("AAA OE HUMAN GSE2 CREEDSID GENE 2 UP", "G1"),
        row("AAA OE HUMAN GSE2 CREEDSID GENE 2 DOWN", "X1"),
        row("DDD MUT HUMAN GSE4 CREEDSID GENE 4 UP", "G1")])
    return [crispr, tfp]


def test_a_gene_that_falls_on_loss_and_rises_on_gain_is_activated(files):
    got = ko.assess_many("G1", ["AAA"], files)["AAA"]
    # CRISPR: up on loss (represses). GEO KO: down on loss (activates). GEO
    # OE: up on gain (activates). Three experiments, two of them agree.
    assert got["class"] == ko.MEASURED
    assert got["n_experiments"] == 3 and got["n_responded"] == 3
    assert (got["activates"], got["represses"]) == (2, 1)
    assert ko.direction_word(got) == "disagrees"
    assert got["species"] == ["human", "human cell lines", "mouse"]


def test_a_factor_perturbed_with_no_response_is_weak_evidence_and_said_so(files):
    got = ko.assess_many("G1", ["BBB"], files)["BBB"]
    assert got["class"] == ko.NO_RESPONSE and got["n_experiments"] == 1
    assert "not among top changes" in ko.cell(got)


def test_a_factor_nobody_perturbed_is_not_tested_and_not_a_no(files):
    got = ko.assess_many("G1", ["ZZZ", "DDD"], files)
    for f in ("ZZZ", "DDD"):     # DDD has only a mutation, which is not read
        assert got[f]["class"] == ko.UNTESTED and got[f]["n_experiments"] == 0
    assert ko.cell(got["ZZZ"]) == "—"


def test_a_gene_that_moves_both_ways_in_one_experiment_is_not_called_one_way(files):
    got = ko.assess_many("G1", ["CCC"], files)["CCC"]
    assert got["n_experiments"] == 1 and got["n_responded"] == 1
    assert ko.direction_word(got) == "disagrees"
    assert "both ways" in ko.cell(got)


def test_the_gene_and_factor_names_are_matched_in_any_case(files):
    got = ko.assess_many("g1", ["aaa"], files)
    assert got["AAA"]["class"] == ko.MEASURED


def test_the_cell_says_which_way_and_how_many(files):
    one = ko.assess_many("X2", ["AAA"], files)["AAA"]
    assert ko.direction_word(one) == "activates" and "falls when" in ko.cell(one)
    rise = ko.assess_many("X1", ["BBB"], files)["BBB"]
    assert "rises when" in ko.cell(rise)
    assert ko.cell(None) == "not checked"


def test_files_are_streamed_and_missing_ones_are_skipped(tmp_path, files):
    got = ko.assess_many("G1", ["AAA"], files + [str(tmp_path / "gone.gmt")])
    assert got["AAA"]["n_experiments"] == 3
    assert ko.resolve_paths(files + [str(tmp_path / "gone.gmt")]) == files
    assert ko.available(files) and not ko.available([str(tmp_path / "x")])
    assert ko.assess_many("G1", ["AAA"], [])["AAA"]["class"] == ko.UNTESTED


def test_a_measured_direction_that_contradicts_the_role_is_flagged():
    act = {"class": ko.MEASURED, "activates": 2, "represses": 0}
    rep = {"class": ko.MEASURED, "activates": 0, "represses": 2}
    assert "annotated as a repressor" in ko.role_conflict(act, "repressor")
    assert "annotated as an activator" in ko.role_conflict(rep, "activator")
    assert ko.role_conflict(act, "activator") is None
    assert ko.role_conflict(rep, "repressor") is None
    mixed = {"class": ko.MEASURED, "activates": 1, "represses": 1}
    assert ko.role_conflict(mixed, "activator") is None
    assert ko.role_conflict(act, None) is None
    assert ko.role_conflict(None, "activator") is None


def test_the_real_libraries_parse_when_present():
    paths = ko.resolve_paths()
    if not paths:
        pytest.skip("no knockout files in data/")
    got = ko.assess_many("PPARG", ["PPARA", "RXRA"], paths)
    assert got["PPARA"]["n_experiments"] > 0
    assert got["RXRA"]["class"] in (ko.MEASURED, ko.NO_RESPONSE)


# ------------------------------------------------------- the outcome chain
def est(tf, klass=None, **extra):
    e = {"drug": "D" + tf, "tf_symbol": tf, "p_engages": 0.8,
         "p_direction": 0.75, "p_activates": 0.6, "p_responsive": 0.4,
         "tissue_dependent": False,
         "notes": ["direction from database annotation only",
                   "gene responsiveness is unknown and assumed 0.40"]}
    if klass is not None:
        e["knockout"] = {"class": klass, "n_experiments": 4,
                         "n_responded": 1 if klass == ko.MEASURED else 0,
                         "activates": 1 if klass == ko.MEASURED else 0,
                         "represses": 0, "experiments": [
                             {"experiment": "T KO MOUSE GSE9", "reads":
                              "activates"}] if klass == ko.MEASURED else []}
    e.update(extra)
    return e


def test_responsiveness_stays_assumed_until_some_factor_has_data():
    plain = [est("A"), est("B")]
    assert "p_responsive" in outcome.assumed_links(plain)
    assert not outcome.knockout_evidence(plain)
    untested = [est("A", ko.UNTESTED)]
    assert not outcome.knockout_evidence(untested)
    assert "p_responsive" in outcome.assumed_links(untested)
    measured = [est("A", ko.MEASURED), est("B", ko.UNTESTED)]
    assert outcome.knockout_evidence(measured)
    assert "p_responsive" not in outcome.assumed_links(measured)


def test_responsiveness_is_never_a_numeric_column():
    for rows in ([est("A")], [est("A", ko.MEASURED)]):
        assert "p_responsive" not in outcome.shown_links(rows)


def test_the_knockout_caveat_is_one_paragraph_and_only_with_evidence():
    lines = outcome.knockout_caveat([est("A", ko.MEASURED),
                                     est("B", ko.UNTESTED),
                                     est("C", ko.NO_RESPONSE)])
    assert len(lines) == 1
    text = lines[0]
    assert "not that it will in your cells" in text
    assert "Not among top changes" in text and "weak" in text
    assert "No knockout data" not in text, "a gap is just a dash in the table"
    assert outcome.knockout_caveat([est("A")]) == []
    assert outcome.knockout_caveat([est("A", ko.UNTESTED)]) == []


# ------------------------------------------------------------- the reports
def _demo_with_knockout():
    res = pipeline.run("SMN2", demo=True, offline=True)
    blk = res["results"][0]
    tfs = sorted({e["tf_symbol"] for e in blk["outcome_estimates"]})
    for i, e in enumerate(blk["outcome_estimates"]):
        k = [ko.MEASURED, ko.NO_RESPONSE, ko.UNTESTED][tfs.index(
            e["tf_symbol"]) % 3]
        e["knockout"] = est(e["tf_symbol"], k)["knockout"]
    return res, tfs


def test_the_factor_table_gets_a_knockout_column_and_no_placeholder():
    res, _ = _demo_with_knockout()
    md = report.markdown_report(res)
    start = md.index("Strength of evidence for each factor route")
    section = md[start:md.index("\n### ", start + 10)]
    header = [l for l in section.splitlines() if "Drug engagement" in l][0]
    assert header.rstrip().endswith("Knockout evidence |")
    assert "Responsiveness" not in header
    assert "falls when the factor is lost" in section
    assert "| — |" in section, "a factor with no data is a dash"
    assert "**Knockout evidence** says whether your gene" in section
    assert "set to 0.4" not in section


def test_without_knockout_data_the_report_simply_has_no_such_column():
    res = pipeline.run("SMN2", demo=True, offline=True)
    md = report.markdown_report(res)
    assert "Knockout evidence" not in md
    for phrase in ("set to 0.4", "not measured by this package",
                   "same for every row"):
        assert phrase not in md, phrase


def test_a_role_that_the_knockout_data_contradict_is_shown_under_the_table():
    res, tfs = _demo_with_knockout()
    blk = res["results"][0]
    target = tfs[0]                        # index 0 -> MEASURED, activates
    for t in blk["tfs"]:
        if t["symbol"] == target:
            t["role"] = "repressor"
    md = report.markdown_report(res)
    assert f"{target}: annotated as a repressor" in md
    assert "one of them is wrong for this gene" in md


def test_the_summary_table_has_a_knockout_column_of_the_right_width():
    res, _ = _demo_with_knockout()
    lines = pipeline.corroboration_table(res)
    header = [l for l in lines if l.startswith("| Drug | Methods")][0]
    assert "Knockout evidence" in header
    sep = lines[lines.index(header) + 1]
    body = [l for l in lines[lines.index(header) + 2:] if l.startswith("| **")]
    assert body and sep.count("---") == header.count("|") - 1 == 13
    assert all(b.count("|") == header.count("|") for b in body)


def test_the_saved_estimates_keep_their_numeric_fields_and_gain_the_knockout():
    res, _ = _demo_with_knockout()
    e = res["results"][0]["outcome_estimates"][0]
    for key in ("p_responsive", "p_point", "p_engages", "knockout"):
        assert key in e
    assert e["p_responsive"] == 0.4, "the placeholder is still saved, unused"
