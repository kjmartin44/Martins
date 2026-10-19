# Setting up the data, from an empty folder

[Documentation index](index.md)

This page takes you from a fresh clone to a machine that can reproduce a full
run. It lists every public source the pipeline reads, what each one is used
for, the exact download command, any step that turns a download into the file
the code expects, and a check that the file arrived intact.

The repository ships none of these tables. You download them once into `data/`.
Every table is optional, and a missing one costs a section of the report that
the report names. [Getting the analysis tables](datasets.md) is the complete
reference for accepted filenames. This page is the walkthrough.

## What you need

- Python 3.10 or later, and about 15 GB of free disk for the recommended set
  (about 30 GB if you add the large optional files).
- `curl`. On Windows 10 and later it is installed as `curl.exe`. In PowerShell
  type `curl.exe`, because plain `curl` is an alias for a different command.
  The commands below use `curl -fL`, so substitute `curl.exe` on Windows.
- A network connection for the downloads. After that, a run can be offline.

Run every command from the repository root (the folder that contains
`pyproject.toml`), because the loaders look for `data/` relative to the working
directory.

```sh
python -m venv .venv
.venv/bin/python -m pip install -e ".[dev,chemistry]"     # Windows: .\.venv\Scripts\python.exe
mkdir -p data
```

On Windows PowerShell use `.\.venv\Scripts\python.exe` in place of `python` in
every command below, and `New-Item -ItemType Directory -Force data` for the
folder.

## Try it with no downloads first

The demo uses synthetic fixtures and never touches the network. If this works,
the install is good and any later problem is in the data.

```sh
python -m promoter_tf_drugs --demo --gene DEMOGENE --out results/demo
```

## The three groups

| Group | Contents | Disk | What you get |
| --- | --- | --- | --- |
| A. Core | motifs, ChIP peaks, gene models, regulons, expression, drug library | about 1.7 GB | the main promoter to factor to drug route |
| B. Evidence | LINCS, dosage, genes and GO, Orange Book, CollecTRI | about 0.3 GB | measured drug effects, dosage limits, route of administration, signed regulons |
| C. Large optional | BindingDB, DepMap, PRISM, GWAS Catalog, ChEMBL database | 12 GB or more | measured potency, dependency, drug sensitivity, human genetics |

Take A and B for a complete run. Add C items only for the sections they unlock.
The run records which sources ran, so skipping one is visible and never silent.

## Group A: core files

### 1. JASPAR 2024 motifs (284 KB)

Used to scan the promoter for transcription factor binding sites.

```sh
curl -fL -o data/jaspar2024.txt \
  "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
head -2 data/jaspar2024.txt        # first line should start with >MA
```

### 2. ReMap 2022 ChIP-seq peaks (1.4 GB)

Measured evidence that a factor binds near the gene. This is the single largest
improvement to the factor ranking. Without it, most motif hits are chance.

```sh
curl -fL --retry 3 -o data/remap2022.bed.gz \
  "https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz"
python -c "import gzip; print(gzip.open('data/remap2022.bed.gz','rt').readline())"
```

The check prints one BED line. On Windows you can run `scripts\download_remap.bat`
instead, which downloads to a temporary name and confirms the result is gzip
before keeping it. ChIP-Atlas peaks also work.

### 3. GENCODE basic annotation (34 MB)

Gene models, transcript start sites, and the check for neighbouring genes, with
no network needed.

```sh
curl -fL -o data/gencode.basic.gtf.gz \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz"
```

The pipeline was developed against release 47. A newer release also works. Keep
the assembly at hg38 so it matches ReMap.

### 4. TRRUST regulons (298 KB)

How many genes each factor controls, from literature curation.

```sh
curl -fL -o data/trrust.tsv \
  "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
```

### 5. ChIP-derived regulon (built, 27 MB)

TRRUST covers few factors. This table counts each factor's targets from the
ReMap peaks, so most factors get a breadth number. It is built, not downloaded,
and needs steps 2 and 3 first. It takes a few minutes.

```sh
python -m promoter_tf_drugs.chipregulon \
  --bed data/remap2022.bed.gz --gtf data/gencode.basic.gtf.gz \
  --out data/chip_regulon.tsv
```

It records occupancy. A peak near a gene is binding, not proof of regulation,
and the report says so.

### 6. Human Protein Atlas tissue expression (5 MB zipped)

Where the gene and each factor are expressed. This is the table the
`--cell-type-check` and factor activity tests read.

```sh
curl -fL -o data/rna_tissue_consensus.tsv.zip \
  "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip"
python -m zipfile -e data/rna_tissue_consensus.tsv.zip data/
head -2 data/rna_tissue_consensus.tsv   # Gene, Gene name, Tissue, nTPM
```

The tracked file `data/hpa.tsv` is a small tissue to organ lookup. It is not an
expression matrix, so do not use it in place of this table.

### 7. Broad Drug Repurposing Hub (about 0.5 MB)

The screening library of approved and clinical compounds, with SMILES structures,
clinical phase, mechanism and target genes. One file replaces most ChEMBL
requests and the known ligands the docking check needs.

The Hub publishes drug annotations and structures as two files, and both begin
with comment lines. Download both, then join them into one table.

```sh
curl -fL -o data/repurposing_drugs.txt \
  "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_drugs_20200324.txt"
curl -fL -o data/repurposing_samples.txt \
  "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/repurposing_samples_20200324.txt"
python scripts/build_hub_library.py data/repurposing_drugs.txt \
  data/repurposing_samples.txt --out data/repurposing_hub.tsv
mkdir -p data/hub_raw && mv data/repurposing_drugs.txt data/repurposing_samples.txt data/hub_raw/
```

The script prints how many drugs it wrote with structures (6,798 for the
2020 release, of which 2,427 are marked Launched). The loaders cannot read the
two raw files, because of the comment lines and because the structures are in
the samples file. `data/repurposing_hub.tsv` is a name the loaders look for.

The last line matters. The loaders search for `data/repurposing_drugs.txt`
before `data/repurposing_hub.tsv`, and the launcher passes the raw file to
`--library` if it is there. Left in place, the raw file is picked up and the
library fails to load with "has no SMILES column". On Windows PowerShell use
`New-Item -ItemType Directory -Force data/hub_raw` and `Move-Item` for that line.

The file is a 2020 snapshot. A drug approved or withdrawn since is not reflected.
Absence from it is not evidence that nothing exists. The full approved library
used in development, about 3,300 compounds, came from ChEMBL (`max_phase 4`),
so a run on the Hub table screens a smaller and older set.

## Group B: evidence files

### 8. LINCS L1000 drug signatures (about 30 MB)

The only source in the pipeline where a drug was observed changing gene
expression. Two Enrichr libraries list the genes each compound moved up and
down. The gene info file says which genes L1000 measured directly and which it
inferred.

```sh
curl -fL -o data/lincs_chem_pert_up.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LINCS_L1000_Chem_Pert_up"
curl -fL -o data/lincs_chem_pert_down.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LINCS_L1000_Chem_Pert_down"
curl -fL -o data/lincs_geneinfo.txt \
  "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt"
```

### 9. CollecTRI signed regulons (2.4 MB)

A regulon where every edge is labelled activating or repressing. The factor
activity test needs it to say which way a factor moved. It covers about three
times as many factors as TRRUST.

```sh
python -m pip install decoupler
python -c "import decoupler as dc; dc.op.collectri(organism='human').to_csv('data/collectri.tsv', sep='\t', index=False)"
```

If your decoupler is version 1.x, use `dc.get_collectri(organism="human")`
instead of `dc.op.collectri(...)`. The file should have about 42,990 edges.

### 10. gnomAD constraint and ClinGen dosage curation (91 MB and 250 KB)

Whether the gene's dose is under selection, and whether a curated panel judges
losing or gaining a copy to be harmful. Both are read from disk only. Without
them the report says dosage was not assessed, which is not the same as the gene
tolerating any dose.

```sh
curl -fL -o data/gnomad.v4.1.constraint_metrics.tsv \
  "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
curl -fL -o data/ClinGen_gene_curation_list_GRCh38.tsv \
  "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
```

### 11. FDA Orange Book (1 MB)

Route of administration for every approved ingredient. It decides whether a
candidate can physically reach the target tissue, and no structure database
records it. Keep the zip as it is, because the loader reads it directly.

```sh
curl -fL -A "Mozilla/5.0" -o data/EOBZIP.zip \
  "https://www.fda.gov/media/76860/download?attachment"
python -m zipfile -l data/EOBZIP.zip     # should list products.txt
```

### 12. HGNC and Gene Ontology annotations (17 MB and 15 MB)

HGNC resolves gene symbols, aliases and UniProt IDs without the network. Gene
Ontology annotations label each factor as an activator, a repressor or unknown.
If the role is unknown the factor cannot be matched to a drug direction, and
that whole route goes quiet when the web services are blocked.

```sh
curl -fL -o data/hgnc_complete_set.txt \
  "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
curl -fL -o data/goa_human.gaf.gz \
  "https://current.geneontology.org/annotations/goa_human.gaf.gz"
```

## Group C: large optional files

Each of these unlocks one specific section. None is needed for a first run.

### 13. BindingDB measured potency (9 GB unpacked, 268 MB distilled)

Measured binding strengths and the known ligands the docking check uses. The
full release has 640 columns and the code reads seven, so distil it once.

Download the current `BindingDB_All` TSV zip from
<https://www.bindingdb.org/rwd/bind/downloads/> (the file name carries the release
date, for example `BindingDB_All_202610_tsv.zip`, 594 MB).

```sh
python -m zipfile -e data/BindingDB_All_202610_tsv.zip data/
python -m promoter_tf_drugs.localpotency --distil data/BindingDB_All.tsv
python -m promoter_tf_drugs.localpotency --uniprot P51681     # CCR5, a quick check
```

After distilling you can delete the zip and the 9 GB TSV. Keep
`data/bindingdb_distilled.tsv`.

### 14. DepMap and PRISM (about 3 GB)

What a cancer cell line loses when a factor is knocked out, and how those lines
respond to about 6,000 drugs. The portal generates release-specific links, so
download these by hand from <https://depmap.org/portal/data_page/> (current
release) and save them under these names:

| Portal file | Save as |
| --- | --- |
| CRISPRGeneEffect.csv | `data/CRISPRGeneEffect.csv` |
| Model.csv | `data/Model.csv` |
| CRISPRInferredCommonEssentials.csv | `data/CRISPRInferredCommonEssentials.csv` |
| PRISM Repurposing Public data matrix | `data/Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` |
| PRISM Repurposing compound list | `data/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` |
| OmicsExpressionProteinCodingGenesTPMLogp1.csv | `data/OmicsExpressionProteinCodingGenesTPMLogp1.csv` |

These measure cancer cell lines growing in a dish. They bound the cost of
acting on a factor and say nothing about a non-dividing tissue.

### 15. GWAS Catalog associations (700 MB)

Human genetic evidence linking variation near a gene or factor to traits. Get
the associations download from <https://www.ebi.ac.uk/gwas/docs/file-downloads>
and save it as `data/gwas_catalog_associations.tsv`. Without it the genetics
stage falls back to Open Targets over the network.

### 16. ChEMBL database (tens of GB)

Only for a fully offline run that needs ChEMBL potency. BindingDB covers the
same need at a fraction of the size, so skip this unless you know you want it.
The SQLite release is at <https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/>
(`chembl_37_sqlite.tar.gz` at the time of writing). Unpack it and save the
database as `data/chembl.db`.

## Sources read over the network at run time

These are not stored in `data/`. Each is requested during a run and cached for
30 days under `~/.cache/promoter_tf_drugs`. A run with `--offline` uses only the
cache and your files.

| Source | Used for |
| --- | --- |
| Ensembl REST | transcripts, start sites, promoter sequence, isoforms for brake exons |
| ChEMBL web API | drug mechanisms, approval phase, potency, warnings, known ligands |
| DGIdb | broader drug and gene interactions |
| Open Targets | drug and target links, tissue expression, genetics evidence |
| OmniPath | signed upstream signalling network |
| Comparative Toxicogenomics Database (CTD) | published chemical effects on expression, to mark what is already known |
| RCSB PDB | experimental structures for docking |
| AlphaFold DB | predicted structures, exploratory mode only |

The first run for a new gene needs the network. A later run for the same gene can
use `--offline`.

## Docking software

Docking is separate from the data tables. It needs a docking program and
RDKit. See [Docking](docking.md) for gnina and smina. To see what your machine
can run:

```sh
python -m promoter_tf_drugs.docking --check
```

## Check everything

Count what is on disk and open each file with the same loader the analysis uses.
This catches a truncated download and a download that is really an error page.

```sh
python -m promoter_tf_drugs.datasets --count
python -m promoter_tf_drugs.datasets --missing      # a shopping list of what is absent
```

Then run a real gene and read the evidence table at the top of `report.md`.
Each source is listed as `ok`, `ok, nothing found`, `FAILED` or `not requested`.
A failed source is never recorded as a zero.

Promoter sequence comes from Ensembl unless you have the hg38 genome FASTA or a
`promoter.fa` from an earlier run, so a first run for a gene needs the network.
This command runs the promoter scan with only the files above:

```sh
python -m promoter_tf_drugs --gene BUB1B --direction up --all-tss \
  --pfm-file data/jaspar2024.txt \
  --chip-bed data/remap2022.bed.gz \
  --regulon-file data/chip_regulon.tsv \
  --expression-file data/rna_tissue_consensus.tsv \
  --gtf data/gencode.basic.gtf.gz \
  --upstream 1000 --downstream 200 --pvalue 1e-5 --tf-weight 0.7 \
  --no-grade --no-outcome --no-improvements --no-mechanisms --no-direct --no-stabilizers \
  --out results/check
```

To run it with no network at all, add `--offline --promoter-fasta
results/check/promoter.fa` after the first run has written that file.

## A tested reference result

This guide was run on an empty folder. For BUB1B (goal up, window 1,000 bases
before to 200 after the start site, p-value 1e-5) a correct setup gives 879
motifs, 68,655,741 ReMap peaks, 37 motif hits, 34 candidate factors, and a
verdict of `interpretable`. If your counts differ, a file is wrong or from another
release.

On Windows, `scripts\run_everything.bat` fetches the small required files for you
(JASPAR, TRRUST, ClinGen, and on request gnomAD and the Orange Book), then runs
docking, the analysis, the funnel and the combined list.

## Final layout

After groups A and B, `data/` should hold these files.

```text
data/
  jaspar2024.txt
  remap2022.bed.gz
  gencode.basic.gtf.gz
  trrust.tsv
  chip_regulon.tsv
  rna_tissue_consensus.tsv
  repurposing_drugs.txt
  repurposing_samples.txt
  repurposing_hub.tsv
  lincs_chem_pert_up.gmt
  lincs_chem_pert_down.gmt
  lincs_geneinfo.txt
  collectri.tsv
  gnomad.v4.1.constraint_metrics.tsv
  ClinGen_gene_curation_list_GRCh38.tsv
  EOBZIP.zip
  hgnc_complete_set.txt
  goa_human.gaf.gz
```

## When a download fails

- **A file is tiny or starts with `<html`.** A captive portal, expired link or
  rate limit returned a web page with a success code. Delete the file and
  download again. Always check the content, because file size alone does not
  catch this.
- **`datasets` shows almost everything missing.** You ran it from outside the
  repository root. Change to the folder that contains `pyproject.toml`.
- **A source list says `FAILED`.** The service was down or blocked. Everything
  fetched is cached, so rerun and the work resumes.
- **A link has moved.** Releases change. Use the source page in
  [Getting the analysis tables](datasets.md), and keep the file name the loaders
  expect.

## Terms of use

Several sources set their own terms. CTD, ReMap, TRRUST and CLUE restrict
commercial use, and others require attribution. Read the terms at each source
before you share a `data/` folder or a results archive. A derived table such as
`chip_regulon.tsv` or `bindingdb_distilled.tsv` carries the terms of its inputs,
and `results.json` embeds values taken from these tables.

Record the release, download date and hash of each file you use. The run itself
records the size, modification time and SHA-256 of every file you pass by flag.
