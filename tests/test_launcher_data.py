"""Direct tests for the lookup, annotation and table-reading functions the
launcher calls and the suite never entered.

Each test hands the function a small real file or a fake client and checks
what it returns, including the cases where it must say it could not tell
rather than return an empty answer.
"""

from __future__ import annotations

import gzip
import sqlite3

import pytest

from generoutes import (agonism, analogs, annotate, direct, drugs,
                               eqtl, genes, localpotency, modulators,
                               network, screening, sensitivity, signatures,
                               stabilizers)
from generoutes.drugs import Drug
from generoutes.http import HttpError

NL = chr(10)
TAB = chr(9)


class _LfClient:
    """Answers `get` and `post` from a list of responses, or raises them."""

    def __init__(self, *responses):
        self.responses = list(responses)
        self.calls = []

    def _next(self, kind, url, params):
        self.calls.append((kind, url, params))
        r = self.responses[min(len(self.calls), len(self.responses)) - 1]
        if isinstance(r, Exception):
            raise r
        return r(url, params) if callable(r) else r

    def get(self, url, params=None, **kw):
        return self._next("get", url, params)

    def post(self, url, body=None, **kw):
        return self._next("post", url, body)


def _lf_drug(name, action, *, approved=True, fa=2000, withdrawn=False,
             direction="up"):
    return Drug(name=name, chembl_id=None, target_symbol="X",
                action_type=action, mechanism=None, approved=approved,
                max_phase=4, first_approval=fa, withdrawn=withdrawn,
                oral=None, direction=direction, source="test")


# ====================================================================== genes
def test_lookup_symbol_asks_for_mane_and_does_not_share_its_parameters():
    client = _LfClient({"id": "ENSG1"})
    assert genes._lookup_symbol(client, "CCR5", "human") == {"id": "ENSG1"}
    kind, url, params = client.calls[0]
    assert url.endswith("/lookup/symbol/human/CCR5")
    assert params["mane"] == 1 and params["expand"] == 1
    params["mane"] = 0
    assert genes._LOOKUP_PARAMS["mane"] == 1, "each call gets its own copy"


def test_lookup_symbol_names_the_symbol_when_ensembl_cannot_resolve_it():
    with pytest.raises(HttpError, match="could not resolve 'NOPE'"):
        genes._lookup_symbol(_LfClient(HttpError("404")), "NOPE", "human")


def test_lookup_id_uses_the_id_endpoint():
    client = _LfClient({"id": "ENSG00000205571"})
    genes._lookup_id(client, "ENSG00000205571")
    assert client.calls[0][1].endswith("/lookup/id/ENSG00000205571")
    assert client.calls[0][2]["mane"] == 1


def test_fetch_cdna_returns_upper_case_sequence_or_none(caplog):
    assert genes.fetch_cdna(_LfClient({"seq": "acgt"}), "ENST1") == "ACGT"
    assert genes.fetch_cdna(_LfClient({"seq": ""}), "ENST1") is None
    assert genes.fetch_cdna(_LfClient({}), "ENST1") is None
    assert genes.fetch_cdna(_LfClient(HttpError("down")), "ENST1") is None


_GTF = NL.join([
    "#!genome-build GRCh38",
    TAB.join(["chr1", "H", "gene", "100", "900", ".", "+", ".",
              'gene_id "G1"; gene_name "TP53";']),
    TAB.join(["chr1", "H", "transcript", "100", "500", ".", "+", ".",
              'gene_id "G1"; gene_name "TP53"; transcript_id "ENST1.2"; '
              'transcript_type "protein_coding"; tag "MANE_Select"; '
              'tag "Ensembl_canonical";']),
    TAB.join(["chr1", "H", "transcript", "300", "900", ".", "-", ".",
              'gene_id "G1"; gene_name "TP53"; transcript_id "ENST2.1";']),
    TAB.join(["chr7", "H", "transcript", "50", "60", ".", "+", ".",
              'gene_id "G2"; gene_name "EGFR"; transcript_id "ENST3.1";']),
    TAB.join(["chr7", "H", "transcript", "50", "60", ".", "+", ".",
              'gene_id "G3"; transcript_id "ENST4.1";']),   # no gene_name
    "short" + TAB + "line", ""])


@pytest.fixture
def clean_genes(monkeypatch):
    monkeypatch.setattr(genes, "_TX_BY_GENE", {})
    monkeypatch.setattr(genes, "_TX_SOURCE", {"path": None})
    monkeypatch.setattr(genes, "_GENOME", {"path": None, "index": {}})
    return genes


def test_load_transcripts_gtf_keeps_the_flags_transcript_selection_needs(
        clean_genes, tmp_path):
    gtf = tmp_path / "g.gtf"
    gtf.write_text(_GTF, encoding="utf-8")
    assert genes.load_transcripts_gtf(gtf) == 2, "TP53 and EGFR only"
    tp53 = genes.local_gene("tp53")
    assert tp53["id"] == "LOCAL:TP53" and tp53["display_name"] == "TP53"
    assert tp53["seq_region_name"] == "1" and tp53["strand"] == 1
    first, second = tp53["Transcript"]
    assert first["id"] == "ENST1", "the version suffix is dropped"
    assert first["is_mane_select"] is True and first["is_canonical"] == 1
    assert first["biotype"] == "protein_coding" and first["strand"] == 1
    assert second["strand"] == -1 and second["is_mane_select"] is False
    assert second["biotype"] == "unknown"
    assert genes.local_gene("MISSING") is None
    assert genes.local_gene("") is None


def test_load_transcripts_gtf_reads_gzip_and_is_idempotent(clean_genes,
                                                           tmp_path):
    gz = tmp_path / "g.gtf.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as fh:
        fh.write(_GTF)
    assert genes.load_transcripts_gtf(gz) == 2
    gz.unlink()                                  # a second call must not reread
    assert genes.load_transcripts_gtf(tmp_path / "g.gtf.gz", force=False) \
        if (tmp_path / "g.gtf.gz").exists() else True
    genes._TX_SOURCE["path"] = str(tmp_path / "g.gtf.gz")
    assert len(genes._TX_BY_GENE) == 2


def test_load_transcripts_gtf_says_when_the_file_is_missing(clean_genes,
                                                           tmp_path):
    with pytest.raises(FileNotFoundError):
        genes.load_transcripts_gtf(tmp_path / "nope.gtf")


_FASTA = NL.join([">chr1 test", "ACGTACGTAC", "GTACGTACGT", "AC",
                  ">chr2", "TTTTTGGGGG", "CC", ""])


def test_index_genome_matches_the_samtools_fai_layout(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_bytes(_FASTA.encode("ascii"))
    idx = genes.index_genome(fa)
    # name -> (length, offset of first base, bases per line, bytes per line)
    first = len(">chr1 test" + NL)
    assert idx["1"] == (22, first, 10, 11)
    assert idx["2"][0] == 12 and idx["2"][2:] == (10, 11)
    fai = tmp_path / "g.fa.fai"
    assert fai.exists()
    rows = [r.split(TAB) for r in fai.read_text(encoding="utf-8").split(NL)
            if r]
    assert [r[0] for r in rows] == ["chr1", "chr2"]
    assert genes.index_genome(fa) == idx, "the index is read back unchanged"


def test_load_genome_makes_windows_readable_on_both_strands(clean_genes,
                                                           tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_bytes(_FASTA.encode("ascii"))
    assert genes.load_genome(fa) == 2
    assert genes.load_genome(fa) == 2, "idempotent"
    # bases 9-13 span a line break: positions 9,10 | 11,12,13
    assert genes.local_sequence_from_genome("chr1", 9, 13) == "ACGTA"
    assert genes.local_sequence_from_genome("1", 1, 4) == "ACGT"
    assert genes.local_sequence_from_genome("1", 1, 4, strand=-1) == "ACGT"
    assert genes.local_sequence_from_genome("1", 1, 3, strand=-1) == "ACG"[::-1] \
        .translate(str.maketrans("ACGT", "TGCA"))
    assert genes.local_sequence_from_genome("1", 20, 23) is None, "past the end"
    assert genes.local_sequence_from_genome("9", 1, 4) is None
    with pytest.raises(FileNotFoundError):
        genes.load_genome(tmp_path / "missing.fa")


_MINUS_GTF = NL.join([
    "#!genome-build GRCh38",
    TAB.join(["chr1", "H", "transcript", "11", "32", ".", "-", ".",
              'gene_id "G9"; gene_name "MINUSG"; transcript_id "ENSTM.1"; '
              'transcript_type "protein_coding"; tag "MANE_Select"; '
              'tag "Ensembl_canonical";']),
    ""])

# Runs of five, so the plus strand and its reverse complement cannot be
# mistaken for each other the way a palindrome like ACGT can.
_ASYMMETRIC_FASTA = NL.join([">chr1", "AAAAACCCCC", "GGGGGTTTTT",
                             "AAAAACCCCC", "GGGGGTTTTT", ""])


def test_a_minus_strand_promoter_from_the_local_genome_is_gene_oriented(
        clean_genes, tmp_path):
    """The window has to come back 5'->3' along the gene, as Ensembl's does.

    `local_sequence_from_genome` defaults to the plus strand, so a caller
    that leaves the argument off gets a silently wrong window for every
    minus-strand gene - half the genome - with the TSS no longer at index
    `upstream` and no error anywhere.
    """
    gtf = tmp_path / "m.gtf"
    gtf.write_text(_MINUS_GTF, encoding="utf-8")
    fa = tmp_path / "m.fa"
    fa.write_bytes(_ASYMMETRIC_FASTA.encode("ascii"))
    genes.load_transcripts_gtf(gtf)
    genes.load_genome(fa)

    up, down = 5, 3
    # `None` for the client: with the GTF and the FASTA loaded nothing may
    # reach Ensembl, and an attempt would raise rather than pass quietly.
    prom, = genes.get_promoters(None, "MINUSG", upstream=up, downstream=down,
                                mask="soft")
    assert prom.strand == -1 and prom.tss == 32
    assert (prom.start, prom.end) == (30, 37)

    plus = genes.local_sequence_from_genome("1", prom.start, prom.end)
    assert plus == "CGGGGGTT", "the window as the FASTA stores it"
    assert prom.sequence == "AACCCCCG", "the same window, along the gene"

    # What every offset in the run depends on: the TSS sits at index
    # `upstream`, and the base there is the gene's, not the genome's.
    assert prom.offset_to_tss(up) == 1
    tss_base = genes.local_sequence_from_genome("1", prom.tss, prom.tss)
    assert tss_base == "G"
    assert prom.sequence[up] == tss_base.translate(
        str.maketrans("ACGT", "TGCA"))


def test_autoload_offline_uses_the_files_it_is_given(clean_genes, tmp_path):
    gtf = tmp_path / "g.gtf"
    gtf.write_text(_GTF, encoding="utf-8")
    fa = tmp_path / "g.fa"
    fa.write_bytes(_FASTA.encode("ascii"))
    assert genes.autoload_offline(gtf=gtf, genome=fa) == {
        "gtf": True, "genome": True}


def test_autoload_offline_reports_what_it_could_not_read(clean_genes,
                                                        tmp_path, caplog,
                                                        monkeypatch):
    import logging
    monkeypatch.chdir(tmp_path)
    assert genes.autoload_offline() == {"gtf": False, "genome": False}
    broken = tmp_path / "bad.gtf.gz"
    broken.write_bytes(b"this is not gzip")
    with caplog.at_level(logging.WARNING):
        got = genes.autoload_offline(gtf=broken)
    assert got["gtf"] is False
    assert any("could not read" in r.message for r in caplog.records)


# ====================================================================== eqtl
_EQTL_HEADER = TAB.join(["SNP", "SNPChr", "SNPPos", "GeneSymbol", "Zscore",
                         "Pvalue", "FDR", "NrSamples"])


def _lf_eqtl_text():
    rows = [
        ["rs1", "3", "46000000", "ccr5", "5.2", "1e-9", "1e-5", "30000"],
        ["rs2", "3", "46050000", "CCR5", "-4.1", "1e-6", "0.01", "30000"],
        ["rs3", "5", "100", "ACTB", "3.0", "1e-4", "0.04", "30000"],
        ["rsX", "5", "notanumber", "ACTB", "3.0", "1e-4", "0.04", "30000"],
        ["short", "5"],
    ]
    return NL.join([_EQTL_HEADER] + [TAB.join(r) for r in rows]) + NL


@pytest.fixture
def clean_eqtl(monkeypatch):
    monkeypatch.setattr(eqtl, "_CACHE", {"rows": None, "targets": None})
    return eqtl


def test_load_trans_parses_plain_and_gzipped_tables_and_skips_bad_rows(
        clean_eqtl, tmp_path, monkeypatch):
    plain = tmp_path / "t.txt"
    plain.write_text(_lf_eqtl_text(), encoding="utf-8")
    rows = eqtl.load_trans(path=str(plain))
    assert [r.snp for r in rows] == ["rs1", "rs2", "rs3"]
    assert rows[0].gene == "CCR5", "gene symbols are upper-cased"
    assert rows[0].snp_pos == 46000000 and rows[0].n == 30000
    assert eqtl.load_trans() is rows, "the second call is served from memory"

    monkeypatch.setattr(eqtl, "_CACHE", {"rows": None, "targets": None})
    gz = tmp_path / "t.txt.gz"
    gz.write_bytes(gzip.compress(_lf_eqtl_text().encode("utf-8")))
    assert len(eqtl.load_trans(path=str(gz))) == 3


def test_load_trans_reads_a_table_saved_with_windows_line_endings(
        clean_eqtl, tmp_path):
    """A stray CR on the last header name made every row fail the column
    lookup, so the table loaded empty and every gene read as having no
    trans-eQTL at all."""
    p = tmp_path / "crlf.txt"
    p.write_bytes(_lf_eqtl_text().replace(NL, chr(13) + NL).encode("utf-8"))
    rows = eqtl.load_trans(path=str(p))
    assert [r.snp for r in rows] == ["rs1", "rs2", "rs3"]
    assert rows[0].n == 30000


def test_a_renamed_column_is_refused_rather_than_read_as_an_empty_table(
        clean_eqtl, tmp_path):
    """Every row reads its columns by name inside a `continue` handler, so
    one renamed column drops all 59,800 of them.

    The empty table that leaves is what `trans_support` reports as "this
    gene has NO significant trans-eQTL (n=31,684 whole blood)" - a parse
    failure quoting a cohort size, which reads as a measurement.
    """
    renamed = _EQTL_HEADER.replace("GeneSymbol", "Gene Symbol")
    body = TAB.join(["rs1", "3", "46000000", "CCR5", "5.2", "1e-9",
                     "1e-5", "30000"])
    p = tmp_path / "renamed.txt"
    p.write_text(NL.join([renamed, body]) + NL, encoding="utf-8")

    with pytest.raises(ValueError) as exc:
        eqtl.load_trans(path=str(p))
    assert "GeneSymbol" in str(exc.value), "it names the column that is gone"
    assert eqtl._CACHE["rows"] is None, "and nothing empty is cached"


def test_a_refused_table_reports_the_question_as_unanswered(clean_eqtl,
                                                            tmp_path):
    """The whole point: `checked` stays False, so the ledger and the report
    say the question was not answered instead of answering it wrongly."""
    renamed = _EQTL_HEADER.replace("Zscore", "Z-score")
    body = TAB.join(["rs1", "3", "46000000", "CCR5", "5.2", "1e-9",
                     "1e-5", "30000"])
    p = tmp_path / "renamed.txt"
    p.write_text(NL.join([renamed, body]) + NL, encoding="utf-8")

    got = eqtl.trans_support("CCR2", "CCR5", path=str(p))
    assert got.checked is False and got.supported is False
    assert "NOT answered" in got.note
    assert "NO significant trans-eQTL" not in got.note


def test_rows_that_do_not_match_a_good_header_are_refused(clean_eqtl,
                                                          tmp_path):
    """The header can name every column and the rows still be unreadable -
    a truncated copy, or one a spreadsheet re-delimited. An empty table
    must not reach a caller that will read it as a result."""
    rows = [",".join(["rs1", "3", "46000000", "CCR5", "5.2", "1e-9",
                      "1e-5", "30000"])] * 3
    p = tmp_path / "commas.txt"
    p.write_text(NL.join([_EQTL_HEADER] + rows) + NL, encoding="utf-8")

    with pytest.raises(ValueError) as exc:
        eqtl.load_trans(path=str(p))
    assert "none of them parsed" in str(exc.value)


def test_a_header_only_table_is_an_answer_not_a_failure(clean_eqtl, tmp_path):
    """No data lines is not the same as data lines that would not parse."""
    p = tmp_path / "empty.txt"
    p.write_text(_EQTL_HEADER + NL, encoding="utf-8")
    assert eqtl.load_trans(path=str(p)) == []


def test_load_trans_needs_a_source(clean_eqtl):
    with pytest.raises(ValueError, match="client or a local path"):
        eqtl.load_trans()


def test_gene_has_any_trans_eqtl(clean_eqtl, tmp_path):
    p = tmp_path / "t.txt"
    p.write_text(_lf_eqtl_text(), encoding="utf-8")
    assert eqtl.gene_has_any_trans_eqtl("ccr5", path=str(p)) is True
    assert eqtl.gene_has_any_trans_eqtl("TP53") is False


def test_trans_support_finds_variants_near_the_factor(clean_eqtl, tmp_path):
    p = tmp_path / "t.txt"
    p.write_text(_lf_eqtl_text(), encoding="utf-8")
    got = eqtl.trans_support("ccr2", "ccr5",
                             tf_locus=("chr3", 46_100_000, 46_200_000),
                             path=str(p))
    assert got.checked and got.supported
    assert got.tf == "CCR2" and got.gene == "CCR5"
    assert [h.snp for h in got.hits] == ["rs1", "rs2"], "best p first"
    assert "Genetic support" in got.note
    assert got.to_dict()["hits"][0]["snp"] == "rs1"


def test_trans_support_says_which_of_four_things_it_could_not_conclude(
        clean_eqtl, tmp_path):
    p = tmp_path / "t.txt"
    p.write_text(_lf_eqtl_text(), encoding="utf-8")
    none = eqtl.trans_support("SP1", "TP53", path=str(p))
    assert none.checked and not none.supported
    assert "NO significant trans-eQTL" in none.note

    no_locus = eqtl.trans_support("SP1", "CCR5", path=str(p))
    assert not no_locus.supported and "was not supplied" in no_locus.note

    far = eqtl.trans_support("SP1", "CCR5", tf_locus=("9", 1, 100),
                             path=str(p))
    assert far.checked and not far.supported and not far.hits
    assert "none within 1 Mb" in far.note


def test_trans_support_is_not_checked_when_the_table_cannot_be_read(
        clean_eqtl, tmp_path):
    got = eqtl.trans_support("SP1", "CCR5", path=str(tmp_path / "gone.txt"))
    assert got.checked is False and got.supported is False
    assert "NOT answered" in got.note


def _lf_support(checked=True, hits=(), supported=False, note=""):
    return eqtl.TransSupport(tf="T", gene="G", checked=checked,
                             supported=supported, hits=list(hits), note=note)


def test_summarise_keeps_four_different_conclusions_apart():
    assert eqtl.summarise("g", []) == ""
    assert "not retrieved" in eqtl.summarise("g", [_lf_support(checked=False)])
    assert "No factor can be supported" in eqtl.summarise(
        "g", [_lf_support(note="NO significant trans-eQTL here")])
    assert "have human genetic support" in eqtl.summarise(
        "g", [_lf_support(supported=True, note="x")])
    assert "none near any of the factors" in eqtl.summarise(
        "g", [_lf_support(note="has some, none near")])


# ================================================================== annotate
def test_annotation_accessors_read_role_and_receptor_flags():
    ann = annotate._blank("ESR1")
    assert ann.role == "unknown"
    assert ann.is_nuclear_receptor is True and ann["source"] == "curated-NR-list"
    other = annotate._blank("ZZZ3")
    assert other.is_nuclear_receptor is False and other["source"] == "none"
    other["role"] = "dual"
    assert other.is_activator is True
    other["role"] = "repressor"
    assert other.is_activator is False


def test_go_ids_accepts_lists_dicts_and_nothing():
    assert annotate._go_ids({}) == []
    assert annotate._go_ids({"go": {"MF": {"id": "a"}, "BP": [{"id": "b"}]}}) \
        == [{"id": "a"}, {"id": "b"}]
    assert annotate._go_ids({"go": {"MF": [{"id": "a"}]}}) == [{"id": "a"}]


def test_from_hit_derives_role_from_go_terms():
    act = sorted(annotate.ACTIVATOR_GO)[0]
    rep = sorted(annotate.REPRESSOR_GO)[0]

    def hit(*ids):
        return {"name": "N", "entrezgene": "7", "go": {"MF": [
            {"id": i} for i in ids]},
            "uniprot": {"Swiss-Prot": "P12345"}}
    assert annotate._from_hit("X", hit(act))["role"] == "activator"
    assert annotate._from_hit("X", hit(rep))["role"] == "repressor"
    assert annotate._from_hit("X", hit(act, rep))["role"] == "dual"
    plain = annotate._from_hit("X", hit("GO:0000000"))
    assert plain["role"] == "unknown" and plain["uniprot"] == "P12345"
    assert plain["name"] == "N" and plain["source"] == "mygene/GO"
    assert annotate._from_hit("ESR1", hit())["ligand_activated"] is True


@pytest.fixture
def no_local_annotation(monkeypatch):
    from generoutes import localannot
    monkeypatch.setattr(localannot, "available",
                        lambda: {"hgnc": False, "goa": False})


def test_annotate_tfs_uses_mygene_and_keeps_every_requested_symbol(
        no_local_annotation):
    act = sorted(annotate.ACTIVATOR_GO)[0]
    client = _LfClient([
        {"query": "SP1", "go": {"MF": [{"id": act}]}, "name": "Sp1"},
        {"query": "GHOST", "notfound": True},
        {"query": "UNREQUESTED", "go": {}},
    ])
    out = annotate.annotate_tfs(client, ["sp1", "ghost", "", "ESR1"])
    assert sorted(out) == ["ESR1", "GHOST", "SP1"]
    assert out["SP1"]["role"] == "activator"
    assert out["GHOST"]["role"] == "unknown"
    assert out["ESR1"]["nuclear_receptor"] is True
    assert client.calls[0][2]["q"] == ["ESR1", "GHOST", "SP1"]


def test_annotate_tfs_degrades_to_the_curated_list_when_mygene_is_down(
        no_local_annotation):
    out = annotate.annotate_tfs(_LfClient(HttpError("blocked")),
                                ["ESR1", "SP1"])
    assert out["ESR1"]["source"] == "curated-NR-list"
    assert out["SP1"]["role"] == "unknown"


def test_annotate_tfs_prefers_the_informative_duplicate(no_local_annotation):
    act = sorted(annotate.ACTIVATOR_GO)[0]
    client = _LfClient([
        {"query": "SP1", "go": {}},
        {"query": "SP1", "go": {"MF": [{"id": act}]}},
    ])
    assert annotate.annotate_tfs(client, ["SP1"])["SP1"]["role"] == "activator"


# ================================================================== agonism
def _lf_activity(mol, pch, smiles="CCO"):
    return {"molecule_chembl_id": mol, "pchembl_value": pch,
            "canonical_smiles": smiles}


def test_fetch_actives_keeps_the_strongest_measurement_per_molecule():
    client = _LfClient({"activities": [
        _lf_activity("CHEMBL1", "6.5"),
        _lf_activity("CHEMBL1", "8.0"),
        _lf_activity("CHEMBL2", "not a number"),
        {"canonical_smiles": "C", "pchembl_value": "7"},      # no molecule id
        _lf_activity("CHEMBL3", "7.1", smiles="CCN"),
    ], "page_meta": {"next": None}})
    out = agonism.fetch_actives(client, "CHEMBL206")
    by_id = {a["chembl_id"]: a for a in out}
    assert sorted(by_id) == ["CHEMBL1", "CHEMBL3"]
    assert by_id["CHEMBL1"]["pchembl"] == 8.0
    assert client.calls[0][2]["pchembl_value__gte"] == agonism.PCHEMBL_CUTOFF


def test_fetch_actives_follows_pages_and_fills_missing_structures():
    page1 = {"activities": [_lf_activity("CHEMBL1", "7", smiles=None)],
             "page_meta": {"next": "more"}}
    page2 = {"activities": [_lf_activity("CHEMBL2", "7")], "page_meta": {}}
    fill = {"molecules": [{"molecule_chembl_id": "CHEMBL1",
                           "molecule_structures": {
                               "canonical_smiles": "CCCC"}}]}
    client = _LfClient(page1, page2, fill)
    out = {a["chembl_id"]: a["smiles"] for a in
           agonism.fetch_actives(client, "CHEMBL206")}
    assert out == {"CHEMBL1": "CCCC", "CHEMBL2": "CCO"}


def test_fetch_actives_drops_a_molecule_it_cannot_draw_and_keeps_going():
    page = {"activities": [_lf_activity("CHEMBL1", "7", smiles=None)],
            "page_meta": {}}
    assert agonism.fetch_actives(_LfClient(page, HttpError("x")),
                                 "CHEMBL206") == []


def test_fetch_actives_returns_what_it_has_when_a_page_fails():
    page = {"activities": [_lf_activity("CHEMBL1", "7")],
            "page_meta": {"next": "more"}}
    out = agonism.fetch_actives(_LfClient(page, HttpError("500")), "T")
    assert [a["chembl_id"] for a in out] == ["CHEMBL1"]


def test_fetch_labelled_ligands_labels_by_action_and_joins_structures():
    mech = {"mechanisms": [
        {"molecule_chembl_id": "C1", "action_type": "agonist"},
        {"molecule_chembl_id": "C2", "action_type": "ANTAGONIST"},
        {"molecule_chembl_id": "C2", "action_type": "AGONIST"},   # first wins
        {"molecule_chembl_id": "C3", "action_type": "BINDING AGENT"},
        {"action_type": "AGONIST"},
    ]}
    mols = {"molecules": [
        {"molecule_chembl_id": "C1", "pref_name": "ONE",
         "molecule_structures": {"canonical_smiles": "CCO"},
         "molecule_properties": {"full_mwt": "46.07"}},
        {"molecule_chembl_id": "C2", "molecule_structures": None},
    ]}
    out = agonism.fetch_labelled_ligands(_LfClient(mech, mols), "T")
    assert [(o["chembl_id"], o["label"], o["name"], o["mw"]) for o in out] \
        == [("C1", "agonist", "ONE", "46.07")]


def test_fetch_labelled_ligands_is_empty_when_chembl_fails():
    assert agonism.fetch_labelled_ligands(_LfClient(HttpError("x")), "T") == []


def _lf_pred(name="D"):
    return agonism.AgonismPrediction(name, "ESR1", "undetermined")


def test_predict_for_drug_trains_one_model_per_target(monkeypatch):
    trained = []
    monkeypatch.setattr(agonism, "fetch_labelled_ligands",
                        lambda client, tid: [{"x": tid}])
    monkeypatch.setattr(agonism, "train_model",
                        lambda ligands: trained.append(ligands) or "MODEL")
    seen = []
    monkeypatch.setattr(
        agonism, "predict",
        lambda model, smi, drug_name, target_symbol:
        seen.append((model, smi)) or _lf_pred(drug_name))
    drug = _lf_drug("Estradiol", "AGONIST")
    drug.target_chembl_id = "CHEMBL206"
    drug.properties = {"canonical_smiles": "CCO"}
    cache = {}
    agonism.predict_for_drug(None, drug, _model_cache=cache)
    agonism.predict_for_drug(None, drug, _model_cache=cache)
    assert seen == [("MODEL", "CCO"), ("MODEL", "CCO")]
    assert len(trained) == 1, "the model is reused for the same target"


def test_predict_for_drug_without_a_target_or_structure_says_why(monkeypatch):
    drug = _lf_drug("Mystery", "AGONIST")
    out = agonism.predict_for_drug(None, drug)
    assert out.predicted == "undetermined"
    assert "no ChEMBL target id" in out.notes[0]

    monkeypatch.setattr(agonism, "fetch_labelled_ligands", lambda c, t: [])
    monkeypatch.setattr(agonism, "train_model", lambda l: None)
    drug.target_chembl_id = "CHEMBL1"
    drug.chembl_id = "CHEMBL9"
    out = agonism.predict_for_drug(_LfClient(HttpError("down")), drug)
    assert out.notes == ["no structure available"]


def test_predict_for_drug_fetches_a_missing_structure(monkeypatch):
    monkeypatch.setattr(agonism, "fetch_labelled_ligands", lambda c, t: [])
    monkeypatch.setattr(agonism, "train_model", lambda l: "M")
    got = {}
    monkeypatch.setattr(
        agonism, "predict",
        lambda model, smi, drug_name, target_symbol:
        got.setdefault("smi", smi) and _lf_pred())
    drug = _lf_drug("Drug", "AGONIST")
    drug.target_chembl_id, drug.chembl_id = "CHEMBL1", "CHEMBL9"
    client = _LfClient({"molecule_structures": {"canonical_smiles": "CCN"}})
    agonism.predict_for_drug(client, drug)
    assert got["smi"] == "CCN"


# ================================================================= screening
def _lf_molecule(cid, name, smiles="CCO"):
    return {"molecule_chembl_id": cid, "pref_name": name, "max_phase": 4,
            "molecule_structures": {"canonical_smiles": smiles}
            if smiles else None}


def test_approved_library_pages_through_and_skips_unstructured_molecules():
    page1 = {"molecules": [_lf_molecule("C1", "aspirin"),
                           _lf_molecule("C2", None),
                           _lf_molecule("C3", "nostructure", smiles=None)],
             "page_meta": {"next": "x"}}
    page2 = {"molecules": [_lf_molecule("C4", "warfarin")],
             "page_meta": {"next": None}}
    out = screening.approved_library(_LfClient(page1, page2), page=3)
    assert [d["chembl_id"] for d in out] == ["C1", "C2", "C4"]
    assert out[0]["name"] == "Aspirin", "names are title-cased"
    assert out[1]["name"] == "C2", "an unnamed molecule is shown by its id"


def test_approved_library_stops_at_the_limit_on_error_and_on_empty():
    many = {"molecules": [_lf_molecule(f"C{i}", "x") for i in range(5)],
            "page_meta": {"next": "x"}}
    assert len(screening.approved_library(_LfClient(many), limit=5)) == 5
    assert screening.approved_library(_LfClient(HttpError("x"))) == []
    assert screening.approved_library(_LfClient({"molecules": []})) == []


def test_local_library_reads_a_table_and_is_empty_without_one(tmp_path):
    p = tmp_path / "lib.tsv"
    p.write_text(TAB.join(["name", "smiles", "chembl_id"]) + NL
                 + TAB.join(["Aspirin", "CC(=O)Oc1ccccc1C(=O)O", "CHEMBL25"])
                 + NL, encoding="utf-8")
    got = screening.local_library(p)
    assert got and got[0]["name"].lower() == "aspirin"
    assert screening.local_library(tmp_path / "missing.tsv") == []


class _LfHit:
    def __init__(self, cid, drug):
        self.chembl_id, self.drug, self.notes = cid, drug, []

    def to_dict(self):
        return {"drug": self.drug, "chembl_id": self.chembl_id,
                "notes": list(self.notes)}


def _lf_lib(n):
    return [{"chembl_id": f"L{i}", "name": f"D{i}", "smiles": "CCO"}
            for i in range(n)]


@pytest.fixture
def screen_world(monkeypatch):
    """screen_receptor with the learning and the web replaced."""
    actives = [{"chembl_id": f"A{i}", "smiles": "CCO", "pchembl": 7.0}
               for i in range(screening.MIN_ACTIVES)]
    monkeypatch.setattr(agonism, "fetch_actives", lambda c, t: actives)
    monkeypatch.setattr(agonism, "fetch_labelled_ligands", lambda c, t: [])
    monkeypatch.setattr(agonism, "train_model", lambda l: None)
    monkeypatch.setattr(screening, "build_binding_model",
                        lambda a, d, seed=0: ("MODEL", 0.9))
    monkeypatch.setattr(screening, "screen",
                        lambda *a, **k: [_LfHit("L1", "D1")])
    monkeypatch.setattr(screening, "enrichment_check", lambda hits: "ENR")
    return {"actives": actives}


def test_screen_receptor_reports_a_complete_screen(screen_world):
    out = screening.screen_receptor(None, "CHEMBL206", library=_lf_lib(100))
    assert out["task"] == "binding" and out["model_auc"] == 0.9
    assert out["n_actives"] == screening.MIN_ACTIVES
    assert out["n_decoys"] == screening.MIN_ACTIVES * 3
    assert out["library_size"] == 100
    assert [h["drug"] for h in out["hits"]] == ["D1"]
    assert out["enrichment"] == "ENR"
    assert "Direction not annotated" in out["direction_note"]


def test_screen_receptor_refuses_without_a_library():
    got = screening.screen_receptor(None, "CHEMBL206", library=[])
    assert got["hits"] == [] and "no screening library" in got["error"]


def test_screen_receptor_refuses_when_there_are_too_few_actives(monkeypatch):
    monkeypatch.setattr(agonism, "fetch_actives",
                        lambda c, t: [{"chembl_id": "A1"}])
    got = screening.screen_receptor(None, "T", library=_lf_lib(10))
    assert got["n_actives"] == 1 and "data gap" in got["error"]


def test_screen_receptor_refuses_a_model_that_cannot_discriminate(
        screen_world, monkeypatch):
    monkeypatch.setattr(screening, "build_binding_model",
                        lambda a, d, seed=0: ("MODEL", 0.55))
    got = screening.screen_receptor(None, "T", library=_lf_lib(100))
    assert got["hits"] == [] and got["model_auc"] == 0.55
    assert "Screen not reported" in got["error"]
    monkeypatch.setattr(screening, "build_binding_model",
                        lambda a, d, seed=0: (None, None))
    assert "could not be trained" in screening.screen_receptor(
        None, "T", library=_lf_lib(100))["error"]


def test_screen_receptor_flags_a_direction_opposite_to_the_goal(
        screen_world, monkeypatch):
    monkeypatch.setattr(agonism, "train_model", lambda l: "DIRMODEL")
    pred = agonism.AgonismPrediction("D1", "T", "antagonist",
                                     probability_agonist=0.1,
                                     confidence="high")
    monkeypatch.setattr(agonism, "predict", lambda *a, **k: pred)
    out = screening.screen_receptor(None, "T", goal="up",
                                    library=_lf_lib(100))
    notes = out["hits"][0]["notes"]
    assert any("predicted antagonist" in n for n in notes)
    assert any("OPPOSITE to your goal" in n for n in notes)
    assert out["direction_note"] is None
    pred2 = agonism.AgonismPrediction("D1", "T", "antagonist",
                                      probability_agonist=0.1,
                                      confidence="high")
    monkeypatch.setattr(agonism, "predict", lambda *a, **k: pred2)
    out = screening.screen_receptor(None, "T", goal="down",
                                    library=_lf_lib(100))
    assert not any("OPPOSITE" in n for n in out["hits"][0]["notes"])


# ================================================================ sensitivity
def test_load_compound_names_maps_ids_to_names(tmp_path):
    p = tmp_path / "compounds.csv"
    p.write_text(NL.join(["IDvar,Broad_ID,Drug_Name", "x,BRD-1,Aspirin",
                          "x,BRD-1,Duplicate", "x,BRD-2,",
                          "x,,Nobody", ""]), encoding="utf-8")
    assert sensitivity.load_compound_names(p) == {"BRD-1": "Aspirin"}


def test_load_compound_names_without_the_columns_or_file(tmp_path, caplog):
    import logging
    p = tmp_path / "c.csv"
    p.write_text("a,b" + NL + "1,2" + NL, encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        assert sensitivity.load_compound_names(p) == {}
    assert any("no id/name column pair" in r.message for r in caplog.records)
    assert sensitivity.load_compound_names(tmp_path / "gone.csv") == {}


# ================================================================ signatures
def test_autoload_gmt_available_checks_the_direction_asked_for(
        monkeypatch, tmp_path):
    up = tmp_path / "up.gmt"
    up.write_text("x", encoding="utf-8")
    monkeypatch.setattr(signatures, "DEFAULT_GMT_PATHS",
                        {"up": (str(up),), "down": (str(tmp_path / "d"),)})
    assert signatures.autoload_gmt_available("up") is True
    assert signatures.autoload_gmt_available("down") is False
    assert signatures.autoload_gmt_available("sideways") is False


# ===================================================== direct / drugs / others
def test_protein_level_drugs_keeps_the_wanted_direction_once_per_name(
        monkeypatch):
    monkeypatch.setattr(direct, "chembl_drugs", lambda c, g, u, a: [
        _lf_drug("Agon", "AGONIST", fa=1999),
        _lf_drug("Blocker", "INHIBITOR"),
        _lf_drug("Stab", "STABILISER"),
    ])
    monkeypatch.setattr(direct, "dgidb_drugs", lambda c, syms, approved_only: [
        _lf_drug("AGON", "AGONIST"),                 # same drug, other source
    ])
    up = direct.protein_level_drugs(None, "gene1", goal="up")
    assert [t.drug for t in up] == ["Agon", "Stab"]
    agon, stab = up
    assert agon.gene == "GENE1" and agon.direction == "up"
    assert agon.raises_protein is False, "an agonist raises activity only"
    assert stab.raises_protein is True
    assert agon.approval_year == 1999 and agon.modality == "small molecule"
    down = direct.protein_level_drugs(None, "gene1", goal="down",
                                      use_dgidb=False)
    assert [t.drug for t in down] == ["Blocker"]
    assert down[0].direction == "down"


def test_chembl_molecules_batches_by_fifty_and_survives_a_failed_batch():
    ids = [f"C{i}" for i in range(120)]

    def batch(url, params):
        chunk = params["molecule_chembl_id__in"].split(",")
        return {"molecules": [{"molecule_chembl_id": c} for c in chunk]}
    client = _LfClient(batch, HttpError("500"), batch)
    out = drugs._chembl_molecules(client, ids)
    assert len(client.calls) == 3
    assert set(out) == set(ids[:50]) | set(ids[100:])
    assert drugs._chembl_molecules(_LfClient(batch), []) == {}


def _lf_chembl_db(path):
    con = sqlite3.connect(path)
    con.executescript("""
        CREATE TABLE component_sequences (component_id INT, accession TEXT);
        CREATE TABLE target_components (component_id INT, tid INT);
        CREATE TABLE assays (tid INT, assay_id INT);
        CREATE TABLE activities (assay_id INT, molregno INT,
            pchembl_value REAL, standard_type TEXT);
        CREATE TABLE molecule_dictionary (molregno INT, chembl_id TEXT);
        CREATE TABLE compound_structures (molregno INT,
            canonical_smiles TEXT);
        INSERT INTO component_sequences VALUES (1, 'P12345');
        INSERT INTO target_components VALUES (1, 10);
        INSERT INTO assays VALUES (10, 100);
        INSERT INTO molecule_dictionary VALUES (1, 'CHEMBL1'), (2, 'CHEMBL2'),
            (3, 'CHEMBL3');
        INSERT INTO compound_structures VALUES (1, 'CCO'), (2, 'CCN'),
            (3, 'CCC');
        INSERT INTO activities VALUES
            (100, 1, 7.0, 'Ki'), (100, 1, 9.0, 'IC50'), (100, 1, 8.0, 'Kd'),
            (100, 2, 6.0, 'EC50'), (100, 3, 4.0, 'Ki'),
            (100, 2, NULL, 'Ki'), (100, 2, 7.5, 'Potency');
    """)
    con.commit()
    con.close()


def test_load_chembl_sqlite_takes_the_median_and_ranks_by_potency(tmp_path):
    db = tmp_path / "chembl.db"
    _lf_chembl_db(db)
    got = localpotency.load_chembl_sqlite(["p12345"], path=str(db))
    rows = got["P12345"]
    assert [r["chembl_id"] for r in rows] == ["CHEMBL1", "CHEMBL2"], \
        "CHEMBL3 is below the 5.0 floor and Potency is not a counted type"
    assert rows[0]["pchembl"] == 8.0 and rows[0]["n_measurements"] == 3
    assert rows[1]["pchembl"] == 6.0
    strict = localpotency.load_chembl_sqlite(["P12345"], path=str(db),
                                             min_pchembl=7.0)
    assert [r["chembl_id"] for r in strict["P12345"]] == ["CHEMBL1"]


def test_load_chembl_sqlite_is_empty_for_no_ids_no_file_or_a_bad_file(
        tmp_path, caplog):
    import logging
    assert localpotency.load_chembl_sqlite([], path="x") == {}
    assert localpotency.load_chembl_sqlite(
        ["P12345"], path=str(tmp_path / "none.db")) == {}
    bad = tmp_path / "bad.db"
    bad.write_bytes(b"this is not sqlite")
    with caplog.at_level(logging.WARNING):
        assert localpotency.load_chembl_sqlite(["P12345"],
                                               path=str(bad)) == {}
    assert any("ChEMBL SQLite query failed" in r.message
               for r in caplog.records)


def test_find_modulators_groups_by_mechanism_and_puts_allosteric_first(
        monkeypatch):
    monkeypatch.setattr(modulators, "chembl_drugs", lambda c, s, u, a: [
        _lf_drug("Plain", "AGONIST", fa=2001),
        _lf_drug("Allo", "POSITIVE ALLOSTERIC MODULATOR", fa=1995),
        _lf_drug("Blocker", "ANTAGONIST"),
        _lf_drug("Old", "AGONIST", fa=1980, withdrawn=True),
    ])
    monkeypatch.setattr(modulators, "dgidb_drugs",
                        lambda c, syms, approved_only: [
                            _lf_drug("PLAIN", "AGONIST")])
    hits = modulators.find_modulators(None, "ESR1", goal="up")
    assert [h.drug for h in hits] == ["Allo", "Plain", "Old"], \
        "allosteric first, then the live drug before the withdrawn one"
    assert hits[0].is_allosteric and hits[0].mechanism == "PAM"
    assert all(h.target == "ESR1" for h in hits)
    down = modulators.find_modulators(None, "ESR1", goal="down",
                                      use_dgidb=False)
    assert [h.drug for h in down] == ["Blocker"]


def _lf_row(src, tgt, stim=1, inhib=0, refs="1;2"):
    return {"source_genesymbol": src, "target_genesymbol": tgt,
            "is_stimulation": stim, "is_inhibition": inhib,
            "references": refs, "sources": "SIGNOR"}


def test_fetch_network_keeps_only_signed_directed_referenced_edges():
    rows = [
        _lf_row("tp53", "cdkn1a"),
        _lf_row("a", "b", stim=0, inhib=1),
        _lf_row("c", "c"),                          # a self loop
        _lf_row("d", "e", stim=0, inhib=0),         # unsigned
        _lf_row("f", "g", stim=1, inhib=1),         # contradictory
        _lf_row("h", "i", refs=""),                 # no reference
        "not a row",
    ]
    edges = network.fetch_network(_LfClient(rows), min_references=1)
    assert [(e.source, e.target, e.sign, e.n_references) for e in edges] == [
        ("TP53", "CDKN1A", 1, 2), ("A", "B", -1, 2)]
    assert edges[0].sources == "SIGNOR"


def test_fetch_network_accepts_the_wrapped_form_and_a_reference_floor():
    wrapped = {"data": [_lf_row("a", "b", refs="1"),
                        _lf_row("c", "d", refs="1;2;3")]}
    edges = network.fetch_network(_LfClient(wrapped), min_references=2)
    assert [e.source for e in edges] == ["C"]


def test_fetch_network_is_empty_when_omnipath_is_unreachable():
    assert network.fetch_network(_LfClient(HttpError("down"))) == []


def test_chembl_stabilisers_keeps_only_stabilising_mechanisms(monkeypatch):
    monkeypatch.setattr(drugs, "chembl_drugs", lambda c, s, u, a: [
        _lf_drug("Tafamidis", "STABILISER"),
        _lf_drug("Chaperone", "pharmacological CHAPERONE", approved=False),
        _lf_drug("Other", "INHIBITOR"),
    ])
    out = stabilizers.chembl_stabilisers(None, "ttr", goal="up")
    assert [(s.drug, s.status, s.tf_symbol) for s in out] == [
        ("Tafamidis", "approved", "TTR"),
        ("Chaperone", "clinical", "TTR")]
    assert out[0].mechanism == "conformational stabiliser"


def test_analog_search_seeds_from_the_best_drug_per_target(monkeypatch):
    seeds = []

    def fake(client, seed, target, already_known, threshold):
        seeds.append((seed, target, sorted(already_known)))
        if target == "T1":
            return [analogs.Analog("Cousin", "C9", 0.8, seed, target)]
        return []
    monkeypatch.setattr(analogs, "find_analogs", fake)
    graded = [
        {"tf_symbol": "T1", "drug_name": "Lead"},
        {"tf_symbol": "T1", "drug_name": "Second"},
        {"tf_symbol": "T2", "drug_name": "Other"},
        {"tf_symbol": None, "drug_name": "Orphan"},
        {"tf_symbol": "T3", "drug_name": ""},
    ]
    out = analogs.search(None, graded)
    assert list(out) == ["T1"]
    assert out["T1"]["seed"] == "Lead"
    assert out["T1"]["analogs"][0]["drug"] == "Cousin"
    assert seeds[0] == ("Lead", "T1", ["Lead", "Second"])
    assert analogs.search(None, graded, limit_seeds=0) == {}


# ============================================ the eQTL step inside pipeline.run
def test_a_whole_offline_run_performs_the_trans_eqtl_lookup(tmp_path,
                                                            monkeypatch):
    """`run.<locals>._trans_eqtl` is a closure inside the pipeline, so it can
    only be reached by running the pipeline. This does that with no network:
    a local GTF, genome and motif file, and a local eQTLGen table placed
    where the client looks for its download."""
    import random

    from generoutes import neighbours, pipeline

    # Module state that run() fills from the files below, restored after.
    monkeypatch.setattr(genes, "_TX_BY_GENE", {})
    monkeypatch.setattr(genes, "_TX_SOURCE", {"path": None})
    monkeypatch.setattr(genes, "_GENOME", {"path": None, "index": {}})
    monkeypatch.setattr(neighbours, "_GTF_GENES", [])
    monkeypatch.setattr(neighbours, "_GTF_SOURCE", None)
    monkeypatch.setattr(eqtl, "_CACHE", {"rows": None, "targets": None})

    rng = random.Random(1)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    # A twelve-base site inside the scanned window. A four-base matrix cannot
    # reach the scan's p-value, and the scan now leaves such matrices out.
    seq = seq[:1800] + "ACGTACGTACGT" + seq[1812:]
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1" + NL + NL.join(seq[i:i + 60]
                                         for i in range(0, 6000, 60)) + NL,
                  encoding="utf-8")
    gtf = tmp_path / "g.gtf"
    gtf.write_text(NL.join([
        TAB.join(["chr1", "H", "gene", "2001", "4000", ".", "+", ".",
                  'gene_id "G1"; gene_name "TESTG"; '
                  'gene_type "protein_coding";']),
        TAB.join(["chr1", "H", "transcript", "2001", "4000", ".", "+", ".",
                  'gene_id "G1"; gene_name "TESTG"; '
                  'transcript_id "ENST1.1"; '
                  'transcript_type "protein_coding"; tag "MANE_Select"; '
                  'tag "Ensembl_canonical";'])]) + NL, encoding="utf-8")
    pfm = tmp_path / "m.txt"
    pfm.write_text(NL.join([
        ">MA0001.1" + TAB + "TEST",
        "A [ 20 0 0 0 20 0 0 0 20 0 0 0 ]",
        "C [ 0 20 0 0 0 20 0 0 0 20 0 0 ]",
        "G [ 0 0 20 0 0 0 20 0 0 0 20 0 ]",
        "T [ 0 0 0 20 0 0 0 20 0 0 0 20 ]"]) + NL, encoding="utf-8")

    # A table over the 100 kB the client requires before it trusts a file.
    rows = [_EQTL_HEADER, TAB.join(["rs1", "3", "46000000", "TESTG", "5.2",
                                    "1e-9", "1e-5", "30000"])]
    for i in range(7000):
        rows.append(TAB.join([f"rs{i + 10}", str(1 + i % 22),
                              str(rng.randrange(1, 10 ** 8)), f"G{i}",
                              f"{rng.uniform(-9, 9):.4f}",
                              f"{rng.random():.6e}", f"{rng.random():.6e}",
                              str(30000 + i)]))
    cache = tmp_path / "cache"
    (cache / "eqtlgen").mkdir(parents=True)
    table = cache / "eqtlgen" / "trans-eqtl-fdr0.05.txt.gz"
    table.write_bytes(gzip.compress((NL.join(rows) + NL).encode("utf-8")))
    assert table.stat().st_size > 100_000, "the fixture must clear the floor"

    res = pipeline.run("TESTG", gtf=str(gtf), genome=str(fa),
                       pfm_file=str(pfm), offline=True, use_genetics=True,
                       include_direct=False, upstream=500, downstream=100,
                       cache_dir=str(cache), grade_drugs=False,
                       estimate_outcome=False, find_improvements=False,
                       compare_mechanisms=False, include_stabilizers=False)
    genetic = res["results"][0]["genetics"]
    assert genetic["trans_eqtl"], "the closure ran and returned its rows"
    first = genetic["trans_eqtl"][0]
    assert first["checked"] is True and first["gene"] == "TESTG"
    assert first["supported"] is False
    assert "locus of TEST was not supplied" in first["note"], (
        "offline, Ensembl cannot place the factor, and the note says that "
        "rather than reporting the factor as unsupported")
    assert "trans-eQTL" in genetic["trans_eqtl_summary"]
    ledger = {row["source"]: row for row in res["evidence_sources"]["sources"]
              } if isinstance(res["evidence_sources"]["sources"], list) \
        else res["evidence_sources"]["sources"]
    assert "trans_eqtl" in ledger
