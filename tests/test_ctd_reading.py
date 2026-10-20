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
