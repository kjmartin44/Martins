"""Mitotic genes, mitosis-disrupting compounds, and brake exons.

For a gene whose job is mitosis, a compound that disturbs mitosis is the wrong
kind of candidate: it moves the gene by changing the cell's cycle. And for a
gene with one working copy, an exon that only decay-prone isoforms use is a
place an antisense oligonucleotide might raise output. Neither is a design or
a verdict; these tests pin what each one reports.
"""

from __future__ import annotations

import pytest

from promoter_tf_drugs import catalogue as cat
from promoter_tf_drugs import cellcycle, pipeline, poisonexon, report
from promoter_tf_drugs import signatures as sg
from promoter_tf_drugs import systemic


# ----------------------------------------------------- the compound tag
def test_the_tag_matches_spellings_and_salts():
    tag = systemic.is_mitotic_disruptor
    assert tag("BI-2536") and tag("BI 2536") and tag("bi2536")
    assert tag("Paclitaxel") and tag("Docetaxel Anhydrous")
    assert tag("Alisertib") and tag("MLN8237") and tag("Tozasertib")
    assert not tag("Vorinostat") and not tag("Aspirin")
    assert not tag("") and not tag(None)


def test_the_tag_is_not_a_class_so_the_classes_stay_disjoint():
    """Paclitaxel is cytotoxic and also a mitotic disruptor. A class would
    have taken it out of one; the tag leaves its class alone."""
    assert systemic.classify("paclitaxel").name == "cytotoxic / antiproliferative"
    assert systemic.is_mitotic_disruptor("paclitaxel")
    assert systemic.classify("bi-2536").name == "antineoplastic kinase inhibitor"
    assert systemic.is_mitotic_disruptor("bi-2536")


# ------------------------------------------------------ the gene set
def test_mitotic_genes_are_recognised_whatever_the_case():
    assert cellcycle.is_mitotic_gene("BUB1B")
    assert cellcycle.is_mitotic_gene("bub1b") and cellcycle.is_mitotic_gene(" PLK1 ")
    assert not cellcycle.is_mitotic_gene("HTT")
    assert not cellcycle.is_mitotic_gene("CCR5")
    assert not cellcycle.is_mitotic_gene("") and not cellcycle.is_mitotic_gene(None)


# ------------------------------------------------ the score demotion
def _block():
    return pipeline.run("SMN2", demo=True, offline=True)["results"][0]


def test_a_mitotic_disruptor_loses_half_its_score_for_a_mitotic_gene(
        monkeypatch):
    block = _block()
    base = cat.build(block, goal="up")
    target = next(r for r in base["routes"] if r["score"] > 0)
    monkeypatch.setattr(systemic, "is_mitotic_disruptor",
                        lambda n: (n or "").lower() == target["drug"].lower())
    hit = cat.build(block, goal="up", mitotic_gene=True)
    miss = cat.build(block, goal="up", mitotic_gene=False)
    assert [r["score"] for r in miss["routes"]] \
        == [r["score"] for r in base["routes"]], "off unless the gene is mitotic"
    for a in (r for r in hit["routes"] if r["drug"] == target["drug"]):
        b = next(r for r in base["routes"] if r["drug"] == a["drug"]
                 and r["acts_on"] == a["acts_on"]
                 and r["category"] == a["category"])
        assert a["score"] == round(b["score"] * cat.MITOTIC_SCORE_FACTOR, 1)
        assert a["mitotic_penalty"] == 0.5 and a["tier"] == b["tier"]
        assert any("suspected to hold cells in mitosis" in n
                   for n in a["notes"])
    others = [(r["drug"], r["score"]) for r in hit["routes"]
              if r["drug"] != target["drug"]]
    assert others == [(r["drug"], r["score"]) for r in base["routes"]
                      if r["drug"] != target["drug"]]


def test_a_compound_that_is_both_broad_and_mitotic_is_halved_not_quartered(
        monkeypatch):
    block = _block()
    base = cat.build(block, goal="up")
    target = next(r for r in base["routes"] if r["score"] > 0)
    monkeypatch.setattr(systemic, "is_mitotic_disruptor",
                        lambda n: (n or "").lower() == target["drug"].lower())
    rec = {"median_genes": 400, "n_signatures": 90, "percentile": 99.0,
           "thin": False, "broad": True, "very_broad": True}
    both = cat.build(block, goal="up", mitotic_gene=True,
                     breadth={sg.breadth_key(target["drug"]): rec})
    a = next(r for r in both["routes"] if r["drug"] == target["drug"])
    b = next(r for r in base["routes"] if r["drug"] == target["drug"]
             and r["acts_on"] == a["acts_on"] and r["category"] == a["category"])
    assert a["score"] == round(b["score"] * 0.5, 1)
    assert a["mitotic_penalty"] == 0.5 and a["breadth_penalty"] == 0.5


# -------------------------------------------- the pipeline and the report
def test_a_mitotic_gene_is_flagged_in_the_result_and_the_report():
    res = pipeline.run("BUB1B", demo=True, offline=True)
    assert res["mitotic_gene"] is True
    md = report.markdown_report(res)
    assert "## A mitotic gene" in md
    assert "BUB1B is a mitotic gene" in md
    assert "synchronised or phase-resolved" in md


def test_a_gene_that_is_not_mitotic_gets_no_caution():
    res = pipeline.run("SMN2", demo=True, offline=True)
    assert res["mitotic_gene"] is False
    assert "## A mitotic gene" not in report.markdown_report(res)


def test_the_candidate_table_and_summary_mark_a_mitotic_disruptor(monkeypatch):
    res = pipeline.run("BUB1B", demo=True, offline=True)
    routes = [r for b in res["results"]
              for r in (b["all_routes"] or {}).get("routes", [])]
    first = sorted(routes, key=lambda r: (r["tier"], -(r.get("score") or 0)))[0]
    table_first = report.corroboration(res)[0]["drug"]
    names = {first["drug"].lower(), table_first.lower()}
    monkeypatch.setattr(systemic, "is_mitotic_disruptor",
                        lambda n: (n or "").lower() in names)
    text = report.summary_report(res)
    assert "+ halts mitosis" in text
    assert f"**{first['drug']}** ⚠halts mitosis" in text
    assert "⚠halts mitosis: the compound is suspected to hold cells" in text


# ------------------------------------------------------- brake exons
def _tx(tid, biotype, *spans):
    return {"id": tid, "biotype": biotype,
            "Exon": [{"start": s, "end": e} for s, e in spans]}


def _gene():
    return {"display_name": "G1", "Transcript": [
        _tx("C1", "protein_coding", (100, 200), (300, 400), (500, 600)),
        _tx("C2", "protein_coding", (100, 200), (500, 600)),
        # novel exon 330-360 (inside the intron between 300-400? no: 450-480)
        _tx("N1", "nonsense_mediated_decay",
            (100, 200), (450, 480), (500, 600)),
        _tx("N2", "nonsense_mediated_decay",
            (100, 200), (450, 480), (500, 600)),
        # a coding exon that starts earlier: extended, not novel. Three
        # isoforms use it, more than use the novel exon, so a sort on support
        # alone would put it first.
        _tx("N3", "nonsense_mediated_decay", (100, 200), (290, 400)),
        _tx("N5", "nonsense_mediated_decay", (100, 200), (290, 400)),
        _tx("N6", "nonsense_mediated_decay", (100, 200), (290, 400)),
        # decay by splicing pattern only: every exon is one coding forms use
        _tx("N4", "nonsense_mediated_decay", (100, 200), (500, 600)),
        _tx("R1", "retained_intron", (100, 700)),
    ]}


def test_a_novel_exon_in_decay_isoforms_is_found_and_counted_once():
    got = poisonexon.find(_gene())
    assert got["n_coding"] == 2 and got["n_nmd"] == 6
    novel = [c for c in got["candidates"] if c["kind"] == "novel exon"]
    assert [(c["start"], c["end"], c["n_transcripts"]) for c in novel] \
        == [(450, 480, 2)], "two isoforms share one exon: one candidate"
    assert novel[0]["length"] == 31
    assert sorted(novel[0]["transcripts"]) == ["N1", "N2"]


def test_an_extended_coding_exon_is_listed_apart_from_a_novel_one():
    got = poisonexon.find(_gene())
    kinds = [c["kind"] for c in got["candidates"]]
    assert kinds[0] == "novel exon", "novel exons come first"
    ext = [c for c in got["candidates"] if c["kind"] != "novel exon"]
    assert [(c["start"], c["end"]) for c in ext] == [(290, 400)]
    assert ext[0]["n_transcripts"] == 3 > got["candidates"][0]["n_transcripts"]
    assert ext[0]["kind"] == "extended or shifted coding exon"


def test_a_decay_isoform_with_no_candidate_exon_is_counted_not_listed():
    got = poisonexon.find(_gene())
    assert got["n_nmd_without_candidate"] == 1          # N4
    assert "3 of them" in got["note"] or "of them with an exon" in got["note"]


def test_a_gene_with_no_decay_isoform_says_what_that_does_and_does_not_mean():
    g = {"display_name": "G2", "Transcript": [
        _tx("C1", "protein_coding", (1, 10), (20, 30))]}
    got = poisonexon.find(g)
    assert got["candidates"] == [] and got["error"] is None
    assert "no nonsense-mediated-decay isoform" in got["note"]
    assert "'basic' GTF omits most of them" in got["note"]


def test_a_record_without_exons_is_an_error_not_an_empty_answer():
    assert poisonexon.find({"Transcript": []})["error"]
    g = {"display_name": "G3", "Transcript": [
        {"id": "T", "biotype": "protein_coding"}]}
    got = poisonexon.find(g)
    assert got["error"] and "no exon coordinates" in got["error"]


# --------------------------------- the pipeline step and the report section
def test_brake_exons_are_only_sought_when_the_goal_is_to_raise_the_gene(
        monkeypatch):
    from promoter_tf_drugs import genes
    monkeypatch.setattr(genes, "_lookup_symbol",
                        lambda client, symbol, species: _gene())
    up = pipeline._brake_exons(None, "G1", "homo_sapiens", "up")
    assert up["candidates"] and up["gene"] == "G1"
    assert pipeline._brake_exons(None, "G1", "homo_sapiens", "down") == {}


def test_a_failed_isoform_lookup_is_reported_not_read_as_no_brake_exon(
        monkeypatch):
    from promoter_tf_drugs import genes
    from promoter_tf_drugs.http import HttpError

    def down(client, symbol, species):
        raise HttpError("Ensembl 500")
    monkeypatch.setattr(genes, "_lookup_symbol", down)
    got = pipeline._brake_exons(None, "G1", "homo_sapiens", "up")
    assert got["candidates"] == [] and "could not be read" in got["error"]


def test_a_loaded_gtf_does_not_replace_the_ensembl_isoform_list(monkeypatch):
    """The GTF records have no exons, so reading them reported an error for
    every gene in a launcher run, which loads one."""
    from promoter_tf_drugs import genes
    monkeypatch.setitem(genes._TX_BY_GENE, "G1", [
        {"id": "T1", "biotype": "protein_coding", "start": 1, "end": 9,
         "strand": 1, "chrom": "1"}])
    monkeypatch.setattr(genes, "_lookup_symbol",
                        lambda client, symbol, species: _gene())
    got = pipeline._brake_exons(None, "G1", "homo_sapiens", "up")
    assert got["error"] is None and got["candidates"]


def _res_with(brake):
    res = pipeline.run("SMN2", demo=True, offline=True)
    res["brake_exons"] = brake
    return res


def test_the_report_lists_novel_exons_and_says_it_is_not_a_design():
    md = report.markdown_report(_res_with(poisonexon.find(_gene())))
    assert "## Exons that may brake SMN2 output" in md
    assert "| 450–480 | 31 bp | novel exon | 2 |" in md
    assert "extended or shifted coding exon" in md
    assert "a list of places to look, not a design" in md
    assert "Nothing here picks an oligonucleotide sequence" in md


def test_the_report_says_when_the_isoforms_could_not_be_read():
    err = {"gene": "SMN2", "candidates": [], "n_nmd": 0, "note": "",
           "error": "the isoform list could not be read (Ensembl 500)"}
    md = report.markdown_report(_res_with(err))
    assert "Not assessed: the isoform list could not be read" in md
    assert "not a finding that SMN2 has no such exon" in md


def test_a_down_goal_or_a_demo_run_has_no_brake_section():
    assert "may brake" not in report.markdown_report(_res_with({}))
    assert pipeline.run("SMN2", demo=True, offline=True)["brake_exons"] == {}
