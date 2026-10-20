"""Reading CTD's expression sentences instead of its action list.

CTD's `InteractionActions` lists every action in a sentence, so a chemical
that raises another gene, or that blocks an inducer of this one, carries
`increases^expression` too. The sentence says which gene and whose effect.
"""

from __future__ import annotations

import pytest

from promoter_tf_drugs import discovery as d


def eff(chem, text, gene="CYP3A4"):
    return d.expression_effects(chem, text, gene)


def test_a_simple_clause_gives_direction_and_form():
    assert eff("Rifampin", "Rifampin results in increased expression of "
               "CYP3A4 mRNA") == [("up", "mRNA", "direct")]
    assert eff("X", "X results in decreased expression of CYP3A4 protein") == [
        ("down", "protein", "direct")]
    assert eff("X", "X affects the expression of CYP3A4 mRNA") == [
        ("affects", "mRNA", "direct")]
    assert eff("X", "X results in increased expression of CYP3A4") == [
        ("up", None, "direct")]


def test_another_genes_expression_is_not_this_genes():
    text = ("[2,2',4,4'-tetrabromodiphenyl ether results in increased "
            "activity of CYP3A4 protein] which results in increased "
            "expression of CYP2B6 mRNA")
    assert eff("2,2',4,4'-tetrabromodiphenyl ether", text) == []


def test_a_shared_object_after_a_coordinate_clause_is_found():
    text = ("Rotenone results in increased expression of and affects the "
            "localization of BUB1B protein")
    assert eff("Rotenone", text, "BUB1B") == [("up", "protein", "direct")]


def test_the_chemical_that_blocks_an_inducer_is_not_read_as_raising_the_gene():
    text = ("Ketoconazole inhibits the reaction [Rifampin results in "
            "increased expression of CYP3A4 mRNA]")
    assert eff("Ketoconazole", text) == [("affects", "mRNA", "modulated")]
    assert eff("Rifampin", text) == [("up", "mRNA", "direct")]
    promoted = "W promotes the reaction [V results in decreased expression of CYP3A4 mRNA]"
    assert eff("W", promoted) == [("affects", "mRNA", "modulated")]
    assert eff("V", promoted) == [("down", "mRNA", "direct")]


def test_a_protein_that_mediates_leaves_the_inducer_direct():
    text = ("NR1I2 protein affects the reaction [Z results in increased "
            "expression of CYP3A4 mRNA]")
    assert eff("Z", text) == [("up", "mRNA", "direct")]


def test_an_effect_through_another_gene_is_marked_indirect():
    text = ("[A results in decreased expression of NR1I2 mRNA] which results "
            "in increased expression of CYP3A4 mRNA")
    assert eff("A", text) == [("up", "mRNA", "indirect")]


def test_an_unknown_chemical_or_empty_sentence_gives_nothing():
    assert eff("", "X results in increased expression of CYP3A4 mRNA") == []
    assert eff("X", "") == []
    assert eff("Absent", "Y results in increased expression of CYP3A4 mRNA"
               ) == [("up", "mRNA", "unparsed")]


ROWS = [
    {"ChemicalName": "Rifampin", "GeneSymbol": "CYP3A4",
     "InteractionActions": "increases^expression",
     "Interaction": "Rifampin results in increased expression of CYP3A4 mRNA"},
    {"ChemicalName": "Ketoconazole", "GeneSymbol": "CYP3A4",
     "InteractionActions": "decreases^reaction|increases^expression",
     "Interaction": "Ketoconazole inhibits the reaction [Rifampin results in "
                    "increased expression of CYP3A4 mRNA]"},
    {"ChemicalName": "Flame retardant", "GeneSymbol": "CYP3A4",
     "InteractionActions": "increases^activity|increases^expression",
     "Interaction": "[Flame retardant results in increased activity of "
                    "CYP3A4 protein] which results in increased expression "
                    "of CYP2B6 mRNA"},
    {"ChemicalName": "Protein only", "GeneSymbol": "CYP3A4",
     "InteractionActions": "decreases^expression",
     "Interaction": "Protein only results in decreased expression of CYP3A4 "
                    "protein"},
    {"ChemicalName": "Lowerer", "GeneSymbol": "CYP3A4",
     "InteractionActions": "decreases^expression",
     "Interaction": "Lowerer results in decreased expression of CYP3A4 mRNA"},
    {"ChemicalName": "No text", "GeneSymbol": "CYP3A4",
     "InteractionActions": "increases^expression", "Interaction": ""},
    {"ChemicalName": "Not expression", "GeneSymbol": "CYP3A4",
     "InteractionActions": "increases^activity",
     "Interaction": "Not expression results in increased activity of CYP3A4"},
]


def test_the_links_keep_only_the_chemicals_that_change_the_transcript():
    links = d.strict_expression_links(ROWS)
    assert links["rifampin"] == ["increases^expression"]
    assert links["lowerer"] == ["decreases^expression"]
    assert links["ketoconazole"] == ["affects^expression"]
    assert "flame retardant" not in links
    assert "protein only" not in links
    assert "not expression" not in links
    assert links["no text"] == ["increases^expression"], (
        "a row with no sentence reads as it always did")
    assert d.link_direction(links["ketoconazole"]) == "unclear"


def test_protein_rows_can_be_asked_for():
    links = d.strict_expression_links(ROWS, forms=("mRNA", "protein", None))
    assert links["protein only"] == ["decreases^expression"]


class _Client:
    def get(self, url, params=None, **kw):
        return ROWS


def test_the_fetch_reads_sentences_unless_asked_for_the_old_reading():
    strict = d.ctd_expression_links(_Client(), "CYP3A4")
    legacy = d.ctd_expression_links(_Client(), "CYP3A4", strict=False)
    assert "flame retardant" not in strict and "flame retardant" in legacy
    assert d.link_direction(strict["ketoconazole"]) == "unclear"
    assert d.link_direction(legacy["ketoconazole"]) == "up"


def test_a_failed_lookup_is_still_none_not_empty():
    class Boom:
        def get(self, *a, **k):
            raise d.HttpError("down")
    assert d.ctd_expression_links(Boom(), "CYP3A4") is None


# ------------------------------------------------- the literature route
from promoter_tf_drugs import catalogue as cat  # noqa: E402


def _literature(entries, papers=None, gene="CYP3A4"):
    return {"gene": gene, "papers": papers or {},
            "reported_links": {"wanted_direction": entries}}


def _entry(name, kind="chemical"):
    return {"chemical": name, "direction": "up", "kind": kind, "actions": ""}


def test_papers_are_counted_for_the_rows_the_links_count():
    rows = [dict(r, PubMedIDs=p) for r, p in zip(ROWS, [
        "1|2", "3", "4", "5", "6|7|8", "9", "10"])]
    got = d.expression_papers(rows)
    assert got["rifampin"] == {"1", "2"}
    assert got["ketoconazole"] == {"3"}, "a modulator is still a known link"
    assert "flame retardant" not in got and "protein only" not in got
    assert got["no text"] == {"9"}, "a row with no sentence is counted as before"


def test_the_papers_lookup_reads_counts_and_survives_a_failed_request():
    rows = [dict(r, PubMedIDs="1|2") for r in ROWS[:1]]

    class Client_:
        def get(self, url, params=None, **kw):
            return rows
    assert d.ctd_expression_papers(Client_(), "CYP3A4") == {"rifampin": 2}

    class Boom:
        def get(self, *a, **k):
            raise d.HttpError("down")
    assert d.ctd_expression_papers(Boom(), "CYP3A4") == {}


def test_each_chemical_ctd_records_the_wanted_way_becomes_a_tier_four_route():
    lit = _literature([_entry("Rifampin"), _entry("Formaldehyde", "exposure"),
                       _entry("Anti-HIV Agents", "category"),
                       _entry("Single")], {"rifampin": 3, "single": 1})
    routes = cat.from_ctd(lit, "up")
    assert [r.drug for r in routes] == ["Rifampin", "Single"], (
        "exposures and headings are not drugs")
    r = routes[0]
    assert r.category == "literature: reported effect" and r.tier == 4
    assert r.directness == 2 and r.source == "ctd" and r.acts_on == "not stated"
    assert r.mechanism == "reported to raise CYP3A4 mRNA"
    assert r.axes.strength == 0.5 and routes[1].axes.strength == 0.35
    assert any("in 3 papers" in n for n in r.notes)
    assert any("in 1 paper" in n and "papers" not in n for n in routes[1].notes)
    assert any("not a prediction" in n for n in r.notes)
    assert cat.from_ctd(lit, "down")[0].mechanism == "reported to lower CYP3A4 mRNA"


def test_no_literature_means_no_routes_and_the_category_is_documented():
    assert cat.from_ctd(None) == [] and cat.from_ctd({}) == []
    assert "literature: reported effect" in cat.CATEGORIES
    guide = cat.CATEGORY_GUIDE["literature: reported effect"]
    assert "3%" in guide["cost"]


def test_the_catalogue_builds_them_from_the_result_block():
    block = {"gene_symbol": "CYP3A4", "literature": _literature(
        [_entry("Rifampin")], {"rifampin": 2})}
    out = cat.build(block, goal="up")
    assert [r["drug"] for r in out["routes"]
            if r["category"] == "literature: reported effect"] == ["Rifampin"]


def test_the_retrospective_never_scores_the_answer_key_against_itself():
    from promoter_tf_drugs import retrospective as rt
    routes = [{"drug": "Rifampin", "category": "literature: reported effect",
               "tier": 4, "score": 50},
              {"drug": "Other", "category": "TF: known regulator", "tier": 4,
               "score": 40}]
    assert [x["drug"] for x in rt.drug_ranking(routes)] == ["Other"]
    for group in rt.GROUPS:
        assert rt.breadth_key("Rifampin") not in rt.route_scores(routes, group)


def test_the_run_and_the_command_line_can_turn_the_route_off():
    import inspect
    from promoter_tf_drugs import cli, pipeline
    assert inspect.signature(pipeline.run).parameters["ctd_route"].default is True
    kw = cli.run_kwargs_from_args(cli.build_parser().parse_args(
        ["--gene", "X", "--no-ctd-route"]))
    assert kw["ctd_route"] is False


# ------------------------------------------------------------ protein level
def test_the_protein_form_is_read_apart_from_the_transcript():
    rows = [
        {"ChemicalName": "A", "GeneSymbol": "G1", "InteractionActions":
         "increases^expression", "PubMedIDs": "1|2",
         "Interaction": "A results in increased expression of G1 protein"},
        {"ChemicalName": "B", "GeneSymbol": "G1", "InteractionActions":
         "increases^expression", "PubMedIDs": "3",
         "Interaction": "B results in increased expression of G1 mRNA"}]

    class Client_:
        def get(self, url, params=None, **kw):
            return rows
    got = d.ctd_protein_effects(Client_(), "G1")
    assert got["links"] == {"a": ["increases^expression"]}
    assert got["papers"] == {"a": 2}
    assert "a" not in d.strict_expression_links(rows), (
        "a protein report is not a transcript report")


def test_a_failed_protein_lookup_is_none_and_an_empty_one_is_not():
    class Boom:
        def get(self, *a, **k):
            raise d.HttpError("down")

    class Empty:
        def get(self, *a, **k):
            return []
    assert d.ctd_protein_effects(Boom(), "G1") is None
    assert d.ctd_protein_effects(Empty(), "G1") == {"links": {}, "papers": {}}


def _plit(mrna, protein, opposite_m=(), opposite_p=(), pp=None, mp=None):
    def side(want, opp):
        return {"wanted_direction": [_entry(n) for n in want],
                "opposite_direction": [_entry(n) for n in opp]}
    return {"gene": "G1", "reported_links": side(mrna, opposite_m),
            "protein_links": side(protein, opposite_p),
            "papers": mp or {}, "protein_papers": pp or {}}


def test_a_protein_only_report_is_a_route_that_says_so():
    routes = cat.from_ctd(_plit([], ["Prot"], pp={"prot": 2}), "up")
    assert [r.drug for r in routes] == ["Prot"]
    r = routes[0]
    assert r.mechanism == "reported to raise G1 protein"
    assert r.tier == 4 and r.axes.strength == 0.5
    assert any("without any change at the promoter" in n for n in r.notes)
    assert any("G1 protein the way the goal needs, in 2 papers" in n
               for n in r.notes)


def test_a_drug_reported_at_both_levels_is_one_route():
    routes = cat.from_ctd(_plit(["Both", "Rna"], ["Both"], mp={"both": 1},
                                pp={"both": 1}), "up")
    assert [r.drug for r in routes] == ["Both", "Rna"]
    both = routes[0]
    assert both.mechanism == "reported to raise G1 mRNA and protein"
    assert both.axes.strength == 0.5, "papers at both levels are added"
    assert not any("without any change" in n for n in both.notes)
    assert routes[1].mechanism == "reported to raise G1 mRNA"


def test_a_report_the_other_way_at_the_other_level_is_noted():
    routes = cat.from_ctd(_plit(["X"], [], opposite_p=["X"]), "up")
    assert any("protein the other way" in n for n in routes[0].notes)
    routes = cat.from_ctd(_plit([], ["Y"], opposite_m=["Y"]), "up")
    assert any("mRNA the other way" in n for n in routes[0].notes)


def test_a_literature_block_from_before_protein_reports_still_builds():
    lit = {"gene": "G1", "reported_links": {"wanted_direction": [
        _entry("Old")]}, "papers": {"old": 1}}
    assert [r.mechanism for r in cat.from_ctd(lit, "up")] == [
        "reported to raise G1 mRNA"]
