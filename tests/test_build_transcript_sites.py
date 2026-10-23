"""The builder for data/transcript_sites.tsv, on one row of each input."""

import importlib.util
from pathlib import Path

import pytest

from generoutes import neighbours, utr

SCRIPT = (Path(__file__).resolve().parents[1] / "scripts"
          / "build_transcript_sites.py")
spec = importlib.util.spec_from_file_location("build_transcript_sites", SCRIPT)
build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)

TARGETSCAN = (
    "Gene ID\tGene Symbol\tTranscript ID\tGene Tax ID\tmiR Family\tSite Type\n"
    "ENSG1\tBUB1B\tENST1\t9606\tmiR-15-5p\t7mer-m8\n"
    "ENSG1\tBub1b\tENSMUST1\t10090\tmiR-15-5p\t7mer-m8\n"
    "ENSG2\tTP53\tENST2\t9606\tmiR-125b-5p\t8mer\n"
)
MIRTARBASE = (
    "miRTarBase ID\tmiRNA\tSpecies (miRNA)\tTarget Gene\t"
    "Species (Target Gene)\tExperiments\tSupport Type\n"
    "MIRT1\thsa-miR-21-5p\tHomo sapiens\tBUB1B\tHomo sapiens\t"
    "Luciferase reporter assay//Western blot\tFunctional MTI\n"
    "MIRT2\tmmu-miR-21-5p\tMus musculus\tBub1b\tMus musculus\tCLIP-seq\t"
    "Non-Functional MTI\n"
    # The csv release writes the species as "hsa" rather than spelling it out.
    "MIRT3\thsa-miR-17-5p\thsa\tBUB1B\thsa\tCLIP-seq\tFunctional MTI\n"
    # A screen hit, not a direct assay: 7.3M of miRTarBase's 7.34M rows.
    "MIRT4\thsa-miR-99-5p\thsa\tBUB1B\thsa\tCLIP-seq\t"
    "Functional MTI (Weak)\n"
    # Measured NOT to act, so it is never a named regulator.
    "MIRT5\thsa-miR-7-5p\thsa\tBUB1B\thsa\tLuciferase reporter assay\t"
    "Non-Functional MTI\n"
)
GTF = (
    "##description: test\n"
    'chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_name "BUB1B"; '
    'gene_type "protein_coding";\n'
    'chr1\tHAVANA\tgene\t5000\t6000\t.\t-\t.\tgene_name "OTHER"; '
    'gene_type "protein_coding";\n'
)
ECLIP = (
    "chr1\t1200\t1300\tRBFOX2_K562_rep01\t1000\t+\n"
    "chr1\t5100\t5200\tRBFOX2_K562_rep01\t1000\t+\n"   # wrong strand for OTHER
    "chr1\t9000\t9100\tRBFOX2_K562_rep01\t1000\t+\n"   # no gene there
)
M6A = "chr1\t1500\t1510\t.\t0\t+\n"


@pytest.fixture(autouse=True)
def _forget_gtf():
    yield
    neighbours._GTF_GENES.clear()
    neighbours._GTF_SOURCE = None


def _inputs(tmp_path):
    paths = {}
    for name, text in (("ts.txt", TARGETSCAN), ("mtb.tsv", MIRTARBASE),
                       ("ann.gtf", GTF), ("RBFOX2_K562_ENCFF1.bed", ECLIP),
                       ("m6a.bed", M6A)):
        p = tmp_path / name
        p.write_text(text, encoding="utf-8")
        paths[name] = str(p)
    return paths


def test_builds_a_table_the_reader_can_use(tmp_path):
    p = _inputs(tmp_path)
    out = tmp_path / "transcript_sites.tsv"
    assert build.main(["--targetscan", p["ts.txt"],
                       "--mirtarbase", p["mtb.tsv"],
                       "--eclip", p["RBFOX2_K562_ENCFF1.bed"],
                       "--m6a", p["m6a.bed"],
                       "--annotation", p["ann.gtf"],
                       "--out", str(out)]) == 0

    sites = utr.annotated_sites("BUB1B", path=str(out))
    assert sites["n"] == 5
    kinds = {k: sorted({r["factor"] for r in rows})
             for k, rows in sites["by_kind"].items()}
    assert kinds == {"microRNA": ["hsa-miR-17-5p", "hsa-miR-21-5p",
                                  "miR-15-5p"],
                     "RBP": ["RBFOX2"], "m6A": ["m6A"]}

    # The prediction stays a prediction, and the assay keeps its own name.
    evidence = {r["factor"]: r["evidence"]
                for rows in sites["by_kind"].values() for r in rows}
    assert evidence["miR-15-5p"] == "predicted"
    assert evidence["hsa-miR-21-5p"] == ("Luciferase reporter assay, "
                                         "Western blot")
    summary = utr.summarise_sites(sites, "BUB1B")
    assert "2 of 3 measured" in summary


def test_drops_other_species_other_strands_and_genes_not_asked_for(tmp_path):
    p = _inputs(tmp_path)
    out = tmp_path / "sites.tsv"
    genes = tmp_path / "genes.txt"
    genes.write_text("# analysed genes\nBUB1B\n", encoding="utf-8")
    assert build.main(["--targetscan", p["ts.txt"],
                       "--mirtarbase", p["mtb.tsv"],
                       "--eclip", p["RBFOX2_K562_ENCFF1.bed"],
                       "--annotation", p["ann.gtf"],
                       "--genes", str(genes), "--out", str(out)]) == 0
    rows = out.read_text(encoding="utf-8").splitlines()[1:]
    assert len(rows) == 4                      # TP53 and OTHER are gone
    assert all(r.split("\t")[0] == "BUB1B" for r in rows)


def test_screen_hits_are_left_out_unless_asked_for(tmp_path):
    """A weak MTI is a screen hit; a non-functional one says it does NOT act."""
    p = _inputs(tmp_path)
    out = tmp_path / "weak.tsv"
    assert build.main(["--mirtarbase", p["mtb.tsv"], "--mirtarbase-weak",
                       "--out", str(out)]) == 0
    rows = utr.annotated_sites("BUB1B", path=str(out))["by_kind"]["microRNA"]
    by_factor = {r["factor"]: r["evidence"] for r in rows}
    assert "hsa-miR-99-5p" in by_factor
    assert by_factor["hsa-miR-99-5p"].endswith("(weak)")
    assert "hsa-miR-7-5p" not in by_factor, "non-functional is never a site"


def test_a_folder_stands_for_the_bed_files_in_it(tmp_path):
    """252 eCLIP peak files is longer than a command line may be."""
    p = _inputs(tmp_path)
    out = tmp_path / "folder.tsv"
    folder = tmp_path / "eclip"
    folder.mkdir()
    (folder / "RBFOX2_K562_ENCFF1.bed").write_text(ECLIP, encoding="utf-8")
    assert build.main(["--eclip", str(folder),
                       "--annotation", p["ann.gtf"], "--out", str(out)]) == 0
    sites = utr.annotated_sites("BUB1B", path=str(out))
    assert sorted({r["factor"] for r in sites["by_kind"]["RBP"]}) == ["RBFOX2"]


def test_bed_inputs_need_an_annotation(tmp_path):
    p = _inputs(tmp_path)
    with pytest.raises(SystemExit):
        build.main(["--eclip", p["RBFOX2_K562_ENCFF1.bed"],
                    "--out", str(tmp_path / "x.tsv")])
