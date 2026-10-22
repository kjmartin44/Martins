"""Every table the analysis can read: where it comes from, and whether you have it.

The problem this solves
-----------------------
Twenty-odd tables feed this package and each one is optional. That is
the right design — a missing file costs you a section and the report
says which — but it makes one question hard to answer: *have I got
everything I meant to get?* A run that quietly skipped half its sources
looks much like one that used them, until you read the evidence table
carefully enough.

So this is the inventory. One entry per table, naming what it buys,
where to get it, every filename the loader will accept, and how many
records came out when it was read.

Why the paths are imported rather than listed
---------------------------------------------
A second copy of `DEFAULT_PRISM_PATHS` would be wrong within a release.
Every entry here imports the constant the loader actually searches, so a
filename added in `sensitivity.py` appears here without anyone
remembering to add it, and a filename removed cannot linger. A test
asserts that every `DEFAULT_*_PATHS` constant in the package is claimed
by an entry, which is what stops a new table being added and silently
missed.

Counting is deliberately cheap
------------------------------
`survey()` reports presence and size for everything, and record counts
only when asked, because counting means parsing — and some of these
files are gigabytes. `python -m generoutes.datasets` lists what
is on disk; `--count` opens them.

What this does not do
---------------------
It does not download anything. The URLs are here so you can see them
beside what you have, and several of these sources want a click-through
or a licence acknowledgement that a script has no business automating.
[The data guide](../docs/data.md) has the full instructions.
"""

from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)


# Supply routes, which change what "missing" means. A table passed by
# flag is not missing when it is absent — it was not asked for.
AUTO = "found by name under data/"
FLAG = "passed with a flag"
BUILT = "built from other files"


@dataclass(frozen=True)
class Dataset:
    """One table, and everything needed to go and get it."""
    key: str
    title: str
    buys: str                  # what having it adds, in one line
    supply: str                # AUTO / FLAG / BUILT
    paths: tuple               # every filename the loader accepts
    url: str | None = None
    size: str | None = None    # approximate, as downloaded
    flag: str | None = None    # the option that passes it, if any
    build: str | None = None   # the command that makes it, if built
    without: str | None = None # what the run loses without it
    counter: object = None     # (path) -> int, or None
    # What the counter counts. "line" means the file was read end to
    # end but not parsed, which catches a truncated download and does
    # not prove the content is what it claims; anything else means the
    # package's own loader read it and produced that many records.
    unit: str = "line"
    notes: tuple = field(default_factory=tuple)

    def found(self) -> Path | None:
        """The first accepted filename that exists, or None."""
        for p in self.paths:
            if Path(p).is_file():
                return Path(p)
        return None


def _paths(module: str, const: str) -> tuple:
    """The loader's own search list, imported rather than restated.

    Flattened, because the constants are not all the same shape:
    `signatures.DEFAULT_GMT_PATHS` is a dict keyed by direction, with a
    list of accepted names under each.
    """
    import importlib
    mod = importlib.import_module(f".{module}", __package__)

    def flat(v):
        if isinstance(v, str):
            yield v
        elif isinstance(v, dict):
            for vv in v.values():
                yield from flat(vv)
        elif isinstance(v, (list, tuple)):
            for vv in v:
                yield from flat(vv)
        else:
            yield str(v)

    return tuple(dict.fromkeys(flat(getattr(mod, const))))


# --------------------------------------------------------------- counters
# Each returns a record count, or raises. Kept tiny: the point is to
# prove the file parses, not to validate it.
def _count_csv_rows(path):
    """Data rows in a CSV, header excluded."""
    import csv
    with open(path, "r", encoding="utf-8", errors="replace",
              newline="") as fh:
        r = csv.reader(fh)
        next(r, None)
        return sum(1 for _ in r)


def _count_gaf(path):
    """Annotation lines in a GAF, which is '!' commented."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        for line in fh:
            if line and not line.startswith(b"!"):
                n += 1
    return n


def _count_bed(path):
    """Peaks in a BED, plain or gzipped, skipping track lines."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        for line in fh:
            if line and not line.startswith((b"#", b"track", b"browser")):
                n += 1
    return n


def _count_sqlite_tables(path):
    """Tables in a SQLite file. Opening it is the check."""
    import sqlite3
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        return con.execute(
            "SELECT count(*) FROM sqlite_master WHERE type='table'"
        ).fetchone()[0]
    finally:
        con.close()


def _count_fasta(path):
    """Sequences in a FASTA."""
    n = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
    return n


def _count_orangebook(path):
    """Lines of `products.txt`, which is the file the loader reads.

    The download is a zip of three tables, so counting its members would
    report "3" for a healthy file and for one that lacks the products table.
    """
    if not str(path).lower().endswith(".zip"):
        return _count_lines(path)
    import zipfile
    with zipfile.ZipFile(path) as z:
        name = next((n for n in z.namelist()
                     if n.lower().endswith("products.txt")), None)
        if name is None:
            raise ValueError("the zip has no products.txt")
        with z.open(name) as fh:
            return sum(1 for _ in fh)


def _count_pfm(path):
    from . import motifs
    return len(motifs.parse_jaspar_pfm_text(
        Path(path).read_text(encoding="utf-8", errors="replace")))


def _count_gnomad(path):
    from . import constraint
    constraint.reset_cache()
    return len(constraint.load_gnomad(str(path)))


def _count_clingen(path):
    from . import constraint
    constraint.reset_cache()
    return len(constraint.load_clingen(str(path)))


def _count_library(path):
    from . import localdrugs
    # approved_only=False, because the count should say how big the
    # table is, not how much of it one run's filter would keep.
    return len(localdrugs.read_table(path, approved_only=False))


def _count_fame(path):
    """Chemicals in the fame table, header and metadata line excluded."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return sum(1 for ln in fh
                   if ln.strip() and not ln.startswith(("#", "chemical	")))


def _count_regulon(path):
    rows = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                rows += 1
    return rows


def _count_lines(path):
    n = 0
    with open(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n


def _count_gzip_lines(path):
    import gzip
    n = 0
    with gzip.open(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n


def registry() -> tuple:
    """Every table, in the order the analysis would miss them."""
    return (
        Dataset(
            key="motifs",
            title="JASPAR CORE vertebrate motifs",
            buys="879 position matrices instead of the small bundled set",
            supply=FLAG, flag="--pfm-file",
            paths=("data/jaspar2024.txt",),
            url="https://jaspar.elixir.no/download/data/2024/CORE/"
                "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
            size="284 KB",
            without="the scan runs on a handful of bundled matrices, so "
                    "most real factors are never tested at all",
            counter=_count_pfm, unit="motif",
        ),
        Dataset(
            key="chip",
            title="ReMap 2022 ChIP-seq peaks",
            buys="measured occupancy, the single largest gain of any file",
            supply=FLAG, flag="--chip-bed",
            paths=("data/remap2022.bed.gz",),
            url="https://remap.univ-amu.fr/",
            size="1.5 GB",
            without="factor ranking rests on sequence alone, and most of a "
                    "motif hit list is chance",
            notes=("The file wanted is remap2022_nr_macs2_hg38_v1_0.bed.gz. "
                   "ChIP-Atlas works too.",),
            counter=_count_bed, unit="peak",
        ),
        Dataset(
            key="gtf",
            title="GENCODE basic annotation",
            buys="transcript lookup and the neighbour check with no network",
            supply=AUTO,
            paths=_paths("neighbours", "DEFAULT_GTF_PATHS"),
            url="https://www.gencodegenes.org/human/",
            size="34 MB",
            flag="--gtf",
            without="TSS lookup needs Ensembl, and the window is never "
                    "checked for other transcription units",
            counter=_count_gzip_lines, unit="line",
        ),
        Dataset(
            key="library",
            title="Broad Drug Repurposing Hub",
            buys="the screening library, clinical phase, and targets per drug",
            supply=AUTO,
            paths=_paths("localdrugs", "DEFAULT_LIBRARY_PATHS"),
            url="https://repurposing.broadinstitute.org/download",
            size="5 MB",
            flag="--library",
            without="the docking screen and the drug library both need "
                    "ChEMBL, so a blocked network stops them",
            counter=_count_library, unit="drug",
            notes=("One file covers the library, the clinical filter and "
                   "the known ligands the docking gate needs.",),
        ),
        Dataset(
            key="regulon_curated",
            title="TRRUST or CollecTRI regulons",
            buys="how many genes each factor controls, from curation",
            supply=AUTO,
            paths=_paths("tfactivity", "DEFAULT_REGULON_PATHS"),
            url="https://www.grnpedia.org/trrust/",
            size="291 KB",
            flag="--regulon-file",
            without="regulon breadth falls back to bundled tiers, so most "
                    "factors are 'not counted' rather than narrow or broad",
            counter=_count_regulon, unit="edge",
            notes=("CollecTRI signs every edge and covers about three times "
                   "as many factors; --tf-activity wants it.",),
        ),
        Dataset(
            key="drug_protein_abundance",
            title="Drug to protein abundance",
            buys="compounds measured to change the gene's own PROTEIN, "
                 "where every other abundance route reads the transcript",
            supply=AUTO,
            paths=_paths("proteomics", "DEFAULT_PATHS"),
            url="https://decryptE.kusterlab.org",
            size="varies",
            flag="on by default once the file exists",
            without="every abundance claim in the run is about mRNA, so a "
                    "stabiliser or a degradation blocker is invisible",
            counter=_count_csv_rows, unit="row",
            notes=("decryptE and ProTargetMiner answer this. LINCS P100 "
                   "(~90 phosphopeptides) and CPTAC (no drug arm) do not.",
                   "Columns: drug, protein, log2fc [, cell_line, dose, "
                   "time, source].",),
        ),
        Dataset(
            key="transcript_expression",
            title="Isoform-level expression",
            buys="how much of the gene's message each isoform carries, "
                 "which decides whether a brake exon is worth an "
                 "oligonucleotide or sits on a transcript nobody makes",
            supply=AUTO,
            paths=_paths("isoforms", "DEFAULT_PATHS"),
            url="https://gtexportal.org/api/v2",
            size="varies",
            flag="fetched from GTEx when absent; the file wins if present",
            without="GTEx is queried live instead, which needs the network "
                    "and gives one tissue of one dataset",
            counter=_count_csv_rows, unit="row",
            notes=("Columns: transcript, tpm [, tissue, gene, source].",
                   "recount3 covers far more samples than GTEx but ships "
                   "as R objects, so it belongs here as a table rather "
                   "than as a live lookup.",),
        ),
        Dataset(
            key="transcript_sites",
            title="Named regulators of the transcript",
            buys="which microRNA or RNA-binding protein acts on the "
                 "transcript, where the UTR scans can only say an element "
                 "is in excess",
            supply=AUTO,
            paths=_paths("utr", "DEFAULT_SITES_PATHS"),
            url="https://www.targetscan.org/",
            size="varies",
            flag="on by default once the file exists",
            without="the UTR scans report motifs with no factor attached, "
                    "so there is nothing to knock down or antagonise",
            counter=_count_csv_rows, unit="row",
            notes=("Columns: gene, factor, kind [, region, evidence, "
                   "source].",
                   "TargetScan (predicted miRNA sites), miRTarBase "
                   "(measured ones), ENCODE eCLIP (RBP footprints) and "
                   "an m6A atlas all reshape into this one table.",),
        ),
        Dataset(
            key="ctd_fame",
            title="Drug fame, from your own saved runs",
            buys="how many other genes the literature records each "
                 "candidate as moving",
            supply=BUILT,
            paths=_paths("fame", "DEFAULT_PATHS"),
            build="python scripts/build_fame_table.py",
            size="30 KB built",
            flag="on by default once the file exists",
            without="no candidate is checked for being a compound the "
                    "literature records for every gene, which is the "
                    "control the retrospective found mattered most",
            counter=_count_fame, unit="chemical",
            notes=("Built from results/, so its resolution is however "
                   "many genes you have run; five is the minimum before "
                   "a count means anything.",
                   "A compound recorded for most genes is telling you "
                   "about the literature, not about yours.",),
        ),
        Dataset(
            key="regulon_chip",
            title="ChIP-derived regulon",
            buys="a target count for the factors curation never reached",
            supply=BUILT,
            paths=("data/chip_regulon.tsv",),
            build="python -m generoutes.chipregulon "
                  "--bed data/remap2022.bed.gz "
                  "--gtf data/gencode.basic.gtf.gz "
                  "--out data/chip_regulon.tsv",
            size="27 MB built",
            flag="--regulon-file",
            without="nothing, if a curated table is present — this covers "
                    "the factors that one misses",
            counter=_count_regulon, unit="edge",
            notes=("Occupancy, not curation. A peak near a gene is binding, "
                   "not regulation, and the report says so.",),
        ),
        Dataset(
            key="expression",
            title="Human Protein Atlas consensus tissue RNA",
            buys="whether a factor is expressed where your gene is",
            supply=FLAG, flag="--expression-file",
            paths=("data/rna_tissue_consensus.tsv",),
            url="https://www.proteinatlas.org/about/download",
            size="46 MB",
            without="--cell-type-check cannot run, and the TF tissue "
                    "breadth axis has no denominator",
            counter=_count_lines,
            notes=("The tissue you pass to --activity-tissue must be "
                   "spelled as this file spells it: 'lymph node', "
                   "'skeletal muscle'.",),
        ),
        Dataset(
            key="potency",
            title="BindingDB measured affinities",
            buys="measured potency, and known ligands for the docking gate",
            supply=AUTO,
            paths=_paths("localpotency", "DEFAULT_BINDINGDB_PATHS"),
            url="https://www.bindingdb.org/bind/chemsearch/marvin/"
                "Download.jsp",
            size="9 GB unpacked, 268 MB distilled",
            build="python -m generoutes.localpotency "
                  "--distil data/BindingDB_All.tsv",
            without="the potency axis in drug scoring has nothing to read "
                    "and the established-ligand ranking loses its order",
            notes=("The full release is 640 columns and seven are read, so "
                   "distil it once: 268 MB afterwards, and read in six "
                   "seconds rather than three minutes. Check it landed "
                   "with `python -m generoutes.localpotency "
                   "--uniprot P51681`.",
                   "Leave `--min-pchembl` alone unless you know why. A "
                   "row-level threshold drops weak Ki measurements and "
                   "lets a stronger IC50 win the assay preference, which "
                   "changes the answers rather than only the file size.",
                   "A censored value ('>10000 nM') is dropped rather than "
                   "converted, because recording it as a number would "
                   "turn a failed experiment into a weak measurement.",),
            counter=_count_lines,
        ),
        Dataset(
            key="chembl",
            title="ChEMBL SQLite release",
            buys="the same potency as BindingDB, authoritative and larger",
            supply=AUTO,
            paths=_paths("localpotency", "DEFAULT_CHEMBL_SQLITE_PATHS"),
            url="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
            size="tens of GB unpacked",
            without="nothing, if BindingDB is present",
            notes=("Queried with the standard library's sqlite3, so it "
                   "adds no dependency.",),
            counter=_count_sqlite_tables, unit="table",
        ),
        Dataset(
            key="lincs_gmt",
            title="LINCS L1000 chemical perturbation signatures",
            buys="the only empirical drug-to-expression link in the package",
            supply=AUTO,
            paths=_paths("signatures", "DEFAULT_GMT_PATHS"),
            url="https://maayanlab.cloud/sigcom-lincs",
            size="16 MB each, up and down",
            flag="--lincs-gmt",
            without="every route rests on annotation, and nothing in the "
                    "run has a measured effect size to compare against",
            counter=_count_lines,
        ),
        Dataset(
            key="lincs_geneinfo",
            title="LINCS gene info",
            buys="maps L1000 probe identifiers to gene symbols",
            supply=AUTO,
            paths=_paths("signatures", "DEFAULT_GENEINFO_PATHS"),
            url="https://clue.io",
            size="1 MB",
            flag="--lincs-geneinfo",
            without="landmark identifiers in the GMT stay unresolved",
            counter=_count_lines,
        ),
        Dataset(
            key="depmap_effect",
            title="DepMap CRISPR gene effect",
            buys="what a cell cannot do without the factor, measured",
            supply=AUTO,
            paths=_paths("depmap", "DEFAULT_EFFECT_PATHS"),
            url="https://depmap.org/portal/download/",
            size="500 MB",
            without="collateral risk rests on a literature-curated regulon "
                    "count, which measures study effort as much as biology",
            notes=("Chronos scores. Zero is no effect, -1 is the median "
                   "common-essential gene.",),
            counter=_count_lines,
        ),
        Dataset(
            key="depmap_model",
            title="DepMap model metadata",
            buys="which tissue each cell line came from",
            supply=AUTO,
            paths=_paths("depmap", "DEFAULT_MODEL_PATHS"),
            url="https://depmap.org/portal/download/",
            size="2 MB",
            without="dependency cannot be read per tissue, only overall",
            counter=_count_lines,
        ),
        Dataset(
            key="depmap_essential",
            title="DepMap common essentials",
            buys="the reference set a selective dependency is judged against",
            supply=AUTO,
            paths=_paths("depmap", "DEFAULT_ESSENTIAL_PATHS"),
            url="https://depmap.org/portal/download/",
            size="100 KB",
            without="'essential everywhere' has nothing to be compared to",
            counter=_count_lines,
        ),
        Dataset(
            key="prism",
            title="PRISM repurposing viability",
            buys="drugs that behave as though they remove your factor, "
                 "with no annotation involved",
            supply=AUTO,
            paths=_paths("sensitivity", "DEFAULT_PRISM_PATHS"),
            url="https://depmap.org/portal/download/",
            size="1 GB",
            flag="--prism-file",
            without="--sensitivity-scan cannot run, so the one route that "
                    "needs no prior annotation is unavailable",
            counter=_count_lines,
        ),
        Dataset(
            key="prism_compounds",
            title="PRISM compound list",
            buys="names and mechanisms for the PRISM columns",
            supply=AUTO,
            paths=_paths("sensitivity", "DEFAULT_COMPOUND_PATHS"),
            url="https://depmap.org/portal/download/",
            size="5 MB",
            without="PRISM hits are identifiers with no drug names",
            counter=_count_lines,
        ),
        Dataset(
            key="ccle_expression",
            title="DepMap/CCLE cell line expression",
            buys="whether sensitivity tracks how much of the gene a line has",
            supply=AUTO,
            paths=_paths("sensitivity", "DEFAULT_EXPRESSION_PATHS"),
            url="https://depmap.org/portal/download/",
            size="1.5 GB",
            flag="--depmap-file",
            without="the expression-dependence half of the scan is skipped",
            counter=_count_lines,
        ),
        Dataset(
            key="gnomad",
            title="gnomAD constraint metrics",
            buys="whether this gene's dose is under selection",
            supply=AUTO,
            paths=_paths("constraint", "DEFAULT_GNOMAD_PATHS"),
            url="https://storage.googleapis.com/gcp-public-data--gnomad/"
                "release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv",
            size="91 MB",
            without="the report says dosage was not assessed, which is not "
                    "the same as the gene tolerating any dose",
            counter=_count_gnomad, unit="gene",
            notes=("One row per transcript. The gene's value comes from "
                   "its MANE Select transcript, falling back to canonical.",),
        ),
        Dataset(
            key="clingen",
            title="ClinGen dosage sensitivity curation",
            buys="a curated panel's verdict on haplo- and triplosensitivity",
            supply=AUTO,
            paths=_paths("constraint", "DEFAULT_CLINGEN_PATHS"),
            url="https://ftp.clinicalgenome.org/"
                "ClinGen_gene_curation_list_GRCh38.tsv",
            size="250 KB",
            without="dosage rests on population constraint alone, which "
                    "cannot see a recessive condition at all",
            counter=_count_clingen, unit="gene",
            notes=("Score 30 means recessive and 40 means dosage "
                   "sensitivity is UNLIKELY — neither is a point on the "
                   "evidence scale, and 0 means nobody has looked.",),
        ),
        Dataset(
            key="hgnc",
            title="HGNC complete set",
            buys="symbol, alias and UniProt resolution with no network",
            supply=AUTO,
            paths=_paths("localannot", "DEFAULT_HGNC_PATHS"),
            url="https://www.genenames.org/download/archive/",
            size="15 MB",
            without="symbol resolution needs Ensembl, and the potency "
                    "tables have no UniProt accession to key on",
            counter=_count_lines,
        ),
        Dataset(
            key="goa",
            title="GO annotations for human",
            buys="whether a factor is annotated as an activator or repressor",
            supply=AUTO,
            paths=_paths("localannot", "DEFAULT_GOA_PATHS"),
            url="https://current.geneontology.org/annotations/goa_human.gaf.gz",
            size="15 MB",
            without="a factor whose role is unresolved cannot be docked "
                    "against a chosen conformation",
            counter=_count_gaf, unit="annotation",
        ),
        Dataset(
            key="gwas",
            title="GWAS Catalog associations",
            buys="human genetic support for moving the gene at all",
            supply=AUTO,
            paths=_paths("localannot", "DEFAULT_GWAS_PATHS"),
            url="https://www.ebi.ac.uk/gwas/docs/file-downloads",
            size="700 MB",
            without="--genetics falls back to Open Targets, or to nothing",
            notes=("A local catalogue reports the gene nearest a signal, "
                   "which is weaker than fine-mapped support. The report "
                   "keeps them apart rather than adding the counts.",),
            counter=_count_lines,
        ),
        Dataset(
            key="orangebook",
            title="FDA Orange Book products",
            buys="route of administration, which no structure source has",
            supply=AUTO,
            paths=_paths("orangebook", "DEFAULT_PATHS"),
            url="https://www.fda.gov/drugs/drug-approvals-and-databases/"
                "approved-drug-products-therapeutic-equivalence-"
                "evaluations-orange-book",
            size="5 MB zipped",
            without="route is known for about seventy hand-curated drugs "
                    "and unknown for the rest",
            notes=("Take the EOBZIP download and keep `products.txt` "
                   "from inside it.",
                   "20,000 products collapse to roughly 1,500 distinct "
                   "ingredients, so it cannot enlarge a screening "
                   "library — route is the reason to load it.",),
            counter=_count_orangebook,
        ),
        Dataset(
            key="knockout",
            title="Factor knockout signatures (Enrichr)",
            buys="whether your gene changes when a factor is lost, which "
                 "replaces a fixed placeholder for responsiveness",
            supply=AUTO,
            paths=_paths("knockout", "DEFAULT_KNOCKOUT_PATHS"),
            flag="--knockout-gmt",
            url="https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text"
                "&libraryName=LINCS_L1000_CRISPR_KO_Consensus_Sigs",
            size="16 MB + 3 MB",
            without="the report says responsiveness is not measured for "
                    "every factor",
            notes=("Two libraries: LINCS_L1000_CRISPR_KO_Consensus_Sigs "
                   "(data/knockout_lincs_crispr_ko.gmt) and "
                   "TF_Perturbations_Followed_by_Expression "
                   "(data/knockout_tf_perturbations.gmt). Each list holds "
                   "only the few hundred genes that changed most.",),
            counter=_count_lines,
        ),
        Dataset(
            key="fda_label_cmax",
            title="Human peak plasma levels from FDA labels",
            buys="a human exposure to set against the concentration the "
                 "evidence was seen at, which ChEMBL has for few drugs",
            supply=BUILT,
            paths=_paths("exposure", "DEFAULT_CMAX_PATHS"),
            build="python scripts/build_fda_cmax.py",
            size="under 100 KB",
            without="the dose check uses ChEMBL's human Cmax records, which "
                    "exist for few approved drugs, and a table you name",
            notes=("Built from openFDA label text. Most labels give no "
                   "usable number, so most drugs get no row; those are "
                   "listed in fda_label_cmax_review.tsv.",),
            counter=_count_lines,
        ),
        Dataset(
            key="genome",
            title="hg38 genome FASTA",
            buys="promoter sequence without calling Ensembl",
            supply=AUTO,
            paths=_paths("genes", "DEFAULT_GENOME_PATHS"),
            url="https://www.gencodegenes.org/human/",
            size="3 GB",
            without="sequence comes from Ensembl, which needs a network",
            counter=_count_fasta, unit="sequence",
        ),
    )


def survey(*, count: bool = False) -> list:
    """One row per table: what it is, whether it is here, how big.

    `count` opens each present file and reports records, which means
    parsing — several of these are gigabytes, so it is off by default.
    """
    rows = []
    for ds in registry():
        found = ds.found()
        row = {
            "key": ds.key,
            "title": ds.title,
            "supply": ds.supply,
            "present": found is not None,
            "path": str(found) if found else None,
            "bytes": found.stat().st_size if found else None,
            "records": None,
            "unit": ds.unit,
            "error": None,
        }
        if count and found is not None and ds.counter is not None:
            try:
                row["records"] = ds.counter(found)
            except Exception as exc:                        # noqa: BLE001
                # A table that is present and unreadable is worse news
                # than one that is absent, so it is reported rather than
                # left looking fine.
                row["error"] = f"{type(exc).__name__}: {str(exc)[:120]}"
        rows.append(row)
    return rows


def _human(n) -> str:
    if n is None:
        return ""
    for unit, scale in (("GB", 1e9), ("MB", 1e6), ("KB", 1e3)):
        if n >= scale:
            return f"{n / scale:.1f} {unit}"
    return f"{n} B"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="python -m generoutes.datasets",
        description="What every analysis table is, where to get it, and "
                    "whether this checkout has it.")
    p.add_argument("--count", action="store_true",
                   help="also open each present file and report how many "
                        "records it holds. This parses them, and some are "
                        "gigabytes")
    p.add_argument("--missing", action="store_true",
                   help="list only what is absent, with where to get it")
    p.add_argument("--key", help="report one table by key")
    args = p.parse_args(argv)

    by_key = {d.key: d for d in registry()}
    if args.key and args.key not in by_key:
        print(f"no such table: {args.key}. Known: "
              f"{', '.join(sorted(by_key))}", file=sys.stderr)
        return 2

    everything = survey(count=args.count)
    # Counted against the whole registry, not against the filtered view.
    # "0 of 12 present" under --missing is true of the filter and false
    # of the checkout, which is the number a reader is actually after.
    have, total = sum(1 for r in everything if r["present"]), len(everything)

    rows = everything
    if args.key:
        rows = [r for r in rows if r["key"] == args.key]
    if args.missing:
        rows = [r for r in rows if not r["present"]]

    if not rows:
        print(f"All {total} tables are present."
              if args.missing else "Nothing to show.")
        return 0

    width = max(len(r["key"]) for r in rows)
    # Named, always. Every loader resolves "data/..." against the
    # working directory, so from the wrong one this prints 1 of 24 and
    # reads as "your downloads are missing" rather than "you are
    # somewhere else" - and the first reading sends someone off to
    # re-download ninety gigabytes.
    print(f"{have} of {total} table(s) present"
          + (f"; {len(rows)} shown" if len(rows) != total else ""))
    print(f"searched {Path('data').resolve()}")
    if have == 0:
        print("  Nothing found. Paths are relative, so run this from the "
              "checkout root.")
    print()
    for r in rows:
        ds = by_key[r["key"]]
        mark = "yes" if r["present"] else "no "
        size = _human(r["bytes"]) or (ds.size or "")
        extra = ""
        if r["error"]:
            extra = f"  !! {r['error']}"
        elif r["records"] is not None:
            # The unit, always. "879 motifs" came from the motif parser
            # and proves the file is a PFM file; "2,225,583 lines" means
            # the stream opened and ran to the end, which catches a
            # truncated download and says nothing about the content.
            extra = f"  {r['records']:,} {ds.unit}(s)"
            if ds.unit == "line":
                extra += ", not parsed"
        print(f"  [{mark}] {r['key']:<{width}}  {ds.title}")
        print(f"        {ds.buys}")
        print(f"        {ds.supply}"
              + (f", {ds.flag}" if ds.flag else "")
              + (f"  ({size})" if size else "") + extra)
        if not r["present"]:
            if ds.url:
                print(f"        get: {ds.url}")
            if ds.build:
                print(f"        build: {ds.build}")
            if ds.without:
                print(f"        without it: {ds.without}")
        print()
    return 0


if __name__ == "__main__":        # pragma: no cover
    logging.basicConfig(level=logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    raise SystemExit(main())
