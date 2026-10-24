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
python -m generoutes --demo --gene DEMOGENE --out results/demo
```

## The four groups

| Group | Contents | Disk | What you get |
| --- | --- | --- | --- |
| A. Core | motifs, ChIP peaks, gene models, regulons, expression, drug library | about 1.7 GB | the promoter motif route, and the tables most routes read |
| B. Evidence | LINCS, dosage, genes and GO, Orange Book, CollecTRI | about 0.3 GB | measured drug effects, dosage limits, route of administration, signed regulons |
| C. Large optional | BindingDB, DepMap, PRISM, GWAS Catalog, ChEMBL database, transcript site catalogues | 13 GB or more | measured potency, dependency, drug sensitivity, human genetics, the microRNAs and RBPs named on a transcript |
| D. Built and curated | tables built from your own runs or from openFDA, and the curated tables inside the package | under 1 MB | the fame check, human peak plasma levels, protein abundance, and the checked curated tables |

Take A and B for a complete run. Add C items only for the sections they unlock.
Group D needs no large download; read it after your first runs.
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
python -m generoutes.chipregulon \
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
python -m generoutes.localpotency --distil data/BindingDB_All.tsv
python -m generoutes.localpotency --uniprot P51681     # CCR5, a quick check
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

### 17. Named regulators of the transcript (about 1 GB downloaded)

The UTR scans can say an AU-rich element or a DRACH site is in excess, but not
which microRNA or RNA-binding protein uses it, and the name is what you can act
on. This table supplies the names. It is built from three downloads, each under
its own terms, and `data/transcript_sites.tsv` keeps the measured rows and the
predicted ones apart.

```sh
curl -fL -o Predicted_Targets_Info.default_predictions.txt.zip \
  "https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip"
unzip Predicted_Targets_Info.default_predictions.txt.zip
curl -fL -C - --retry 6 --retry-delay 10 --retry-all-errors -o hsa_MTI.csv \
  "https://awi.cuhk.edu.cn/miRTarBase/downloads/files/11.0/hsa_MTI.csv"
```

The miRTarBase file is about 974 MB, and the server closed the connection
partway through (at about 444 MB on 2026-10-23). `-C -` continues from what is
already on disk, and the retry flags do it again after a drop. At 100 to 250 KB/s
the file takes one to two hours. The server answers 403 for a while after
several quick requests, so do not restart the command in a loop.

The ENCODE eCLIP peaks are 252 files, one per RNA-binding protein and cell
line. Ask the portal which they are, then save each under a name that starts
with its protein, because ENCODE leaves the BED name column as `.`:

```sh
curl -fL -H "Accept: application/json" -o eclip.json \
  "https://www.encodeproject.org/search/?type=File&assay_title=eCLIP&file_type=bed+narrowPeak&assembly=GRCh38&status=released&limit=all&format=json&field=href&field=accession&field=preferred_default&field=target.label&field=biosample_ontology.term_name"
python - <<'PY'
import json, re, urllib.request
for f in json.load(open("eclip.json"))["@graph"]:
    if not f.get("preferred_default"):
        continue
    name = "_".join(re.sub(r"[^A-Za-z0-9]", "", s) for s in
                    (f["target"]["label"],
                     f["biosample_ontology"]["term_name"], f["accession"]))
    urllib.request.urlretrieve(
        "https://www.encodeproject.org" + f["href"], f"eclip/{name}.bed.gz")
PY
```

Then build the table. The builder fetches nothing; it reads the files you pass.

```sh
python scripts/build_transcript_sites.py \
  --targetscan Predicted_Targets_Info.default_predictions.txt \
  --mirtarbase hsa_MTI.csv \
  --eclip eclip/ \
  --annotation data/gencode.basic.gtf.gz \
  --out data/transcript_sites.tsv
```

`--genes a-file-of-symbols.txt` keeps the table to the genes you actually run,
which matters because TargetScan alone is millions of rows. `--m6a` takes a
methylation atlas in the same BED shape, which is what replaces the DRACH
motif scan with measured sites. `setup_data.py --only transcript_sites` does
all of this for you.

miRTarBase's `Functional MTI (Weak)` class — screen hits rather than direct
assays, and 7.3 million of its 7.34 million human rows — is left out, which
is what keeps the table near 115 MB instead of 880 MB. `--mirtarbase-weak`
puts it back, with every row labelled weak. `Non-Functional MTI` is always
left out, because it records that the pair does not act.

## Group D: built tables and the curated tables

Two kinds of table are not downloads.

A **curated table** is typed into the Python source, so it comes with the
install and needs no setup. They include the drug classes, the nuclear
receptors and the docking residue codes. [Curated tables](curated-tables.md)
lists each one, says what a complete version would be, and gives the steps to
rebuild or check it. After installing, check them:

```sh
python -m generoutes.audit
```

That needs no network. Add `--online` to also check every drug name against
ChEMBL, every gene symbol against HGNC, the nuclear receptor list against the
HGNC receptor group, and the docking residue codes against the PDB. Run it from
the repository root, because two of its checks read `data/`.

A **built table** is made from files you already have, from a web service, or
from your own earlier runs. Git ignores them, so each checkout builds its own.
Every one is optional, and a run names the section it loses without one.

| File | Built from | Command | Without it |
| --- | --- | --- | --- |
| `data/chip_regulon.tsv` | ReMap peaks and GENCODE | step 5 above | factors outside TRRUST have no target count |
| `data/repurposing_hub.tsv` | the Hub's two downloads | step 7 above | the library comes from `approved_library.tsv`, or ChEMBL one request at a time |
| `data/approved_library.tsv` | the ChEMBL web API | `python scripts/setup_data.py --only library` | the docking screen and the upstream drug index query ChEMBL one compound at a time |
| `data/bindingdb_distilled.tsv` | the BindingDB release | step 13 above | the potency axis has nothing to read |
| `data/transcript_sites.tsv` | TargetScan, miRTarBase and ENCODE eCLIP | step 17 above | UTR scans show motifs with no named factor |
| `data/ctd_fame.tsv` | your saved runs in `results/` | `python scripts/build_fame_table.py` | no candidate is checked for being recorded for most genes |
| `data/fda_label_cmax.tsv` and `data/fda_label_cmax_review.tsv` | openFDA drug labels | `python scripts/build_fda_cmax.py` | the dose check uses ChEMBL's human peak levels, which exist for few drugs |
| `data/drug_protein_abundance.tsv` | a decryptE or ProTargetMiner export | `python scripts/build_protein_abundance.py EXPORT.csv` | the protein abundance route does not run |

Notes on the table:

- `approved_library.tsv` holds the approved small molecules in ChEMBL that have
  a structure, about 3,300 of them, with a withdrawn flag and an oral flag.
  `repurposing_hub.tsv` holds about 6,800 and is read first when both exist,
  because the loaders try `repurposing_drugs.txt`, then `repurposing_hub.tsv`,
  then `approved_library.tsv`. Move the Hub files out of `data/` to screen the
  ChEMBL library instead.
- `ctd_fame.tsv` counts, for each chemical, how many distinct genes in your
  saved runs have CTD recording it. It needs runs for at least five genes, and
  its resolution is however many genes you have run. It makes no network
  request. Rebuild it after new runs.
- `fda_label_cmax.tsv` reads the peak concentration sentences in the
  pharmacokinetics section of each label. It looks up the drugs in
  `results/*/results.json`, or the ones you name with `--drugs`, and keeps rows
  from earlier runs. Most labels give no usable number, and those drugs go to
  the review file with the reason, which is a result and not a failure. openFDA
  allows 1,000 requests a day without a key. Put a key in `OPENFDA_API_KEY` for
  more.
- `drug_protein_abundance.tsv` has no public download in the shape the package
  reads. Download the decryptE or ProTargetMiner export by hand, then convert it.
  The script renames columns and does not average replicates.
- A table another table needs is built after it: `chip_regulon.tsv` and
  `transcript_sites.tsv` need `gencode.basic.gtf.gz` first. `setup_data.py`
  orders them for you.

## Sources read over the network at run time

These are not stored in `data/`. Each is requested during a run and cached for
30 days under `~/.cache/generoutes`. A run with `--offline` uses only the
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
python -m generoutes.docking --check
```

## Check everything

Count what is on disk and open each file with the same loader the analysis uses.
This catches a truncated download and a download that is really an error page.

```sh
python -m generoutes.datasets --count
python -m generoutes.datasets --missing      # a shopping list of what is absent
```

Then run a real gene and read the evidence table at the top of `report.md`.
Each source is listed as `ok`, `ok, nothing found`, `FAILED` or `not requested`.
A failed source is never recorded as a zero.

Promoter sequence comes from Ensembl unless you have the hg38 genome FASTA or a
`promoter.fa` from an earlier run, so a first run for a gene needs the network.
This command runs the promoter scan with only the files above:

```sh
python -m generoutes --gene BUB1B --direction up --all-tss \
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

This guide was run on an empty folder on 2026-10-23. For BUB1B (goal up, window
1,000 bases before to 200 after the start site, p-value 1e-5) a correct setup
gives 879 motifs loaded, 68,655,741 ReMap peaks (step 5 prints the count when it
builds the regulon), 17 motif hits, 17 candidate factors, and a verdict of
`interpretable`. GENCODE release 47 and release 50 gave the same result, so
`setup_data.py`, which takes the newest release, and the step 3 command agree.
If your counts differ, a file is wrong or from another release. Earlier versions
of this page gave 37 hits and 34 factors, which the current code no longer
reproduces.

On Windows, one command fetches every reference file that git does not carry:

```bat
scripts\setup_data.bat
```

It lists what is missing, asks once, and fetches the rest. Files already in
`data/` are left alone, so it is safe to rerun after a failure. Options:

| Option | Effect |
| --- | --- |
| `--yes` | Fetch without asking. |
| `--list` | Show what is present and missing, and stop. |
| `--small` | Skip the files over about 500 MB (ReMap, BindingDB, DepMap expression, the transcript site catalogues). |
| `--only remap,collectri` | Fetch only the named items. |
| `--skip bindingdb` | Leave an item out. |
| `--with-optional` | Also fetch the CREEDS and DrugMatrix signatures. |

The items are JASPAR, TRRUST, CollecTRI, ClinGen, gnomAD, the Orange Book, GENCODE,
HPA tissue RNA, HGNC, GO annotations, LINCS, the knockout libraries, the approved
drug library from ChEMBL, DepMap, PRISM, the GWAS Catalog, ReMap, the regulon built
from ReMap, BindingDB (about 9 GB while it is unpacked, 270 MB after it is
distilled), and the named transcript regulators built from TargetScan,
miRTarBase and the ENCODE eCLIP peaks (about 1.4 GB downloaded). The whole set
is roughly 4.5 GB. Each download goes to a `.part` file and
is kept only if it passes a check on its content or its published checksum.
The same holds for the tables that are expanded from a zip or built by a
second pass, such as the HPA expression table, the ChIP regulon and the
transcript sites: they are written beside their final name and moved into
place only once the step finishes. That matters because a file is counted
as present when it merely exists and is not empty, so a half-written one
would be accepted for ever. `scripts/run_everything.bat` fetches its own
small reference files the same way, reading curl's exit code as well as
checking the content — a truncated JASPAR file still contains its first
matrix header, so the content check alone cannot see a dropped connection.

Four files are not fetched. The ChEMBL SQLite release and the hg38 genome are
optional and very large. `fda_label_cmax.tsv` is built from a finished run by
`scripts/build_fda_cmax.py`, which `run_everything.bat` does at the end. And
`drug_protein_abundance.tsv` has no public download in the shape this package
reads: decryptE and ProTargetMiner are manual downloads, and
`scripts/build_protein_abundance.py` converts either one. Without it the
protein abundance route does not run, and the report says so rather than
reporting that nothing moves the protein.

CollecTRI comes from `decoupler` when it is installed, and otherwise from OmniPath's
CollecTRI dataset. The OmniPath table has more signed edges (about 58,000 against
43,000) and agrees with decoupler's table on the sign of 99.9% of the edges they
share. The script says which source it used. Install `decoupler` first if you need
the same table as an earlier run.

`scripts\run_everything.bat` offers to run the setup script when a key file is
missing, then runs docking, the analysis, the funnel and the combined list.

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
- **`HTTP Error 502` or 503 from one host.** The server is down, and nothing on
  your side will fix it. `setup_data.py` carries on with the other items and
  lists the failures at the end. ReMap answered 502 for over an hour on
  2026-10-23. Rerun later with `--only remap,chip_regulon`, or put a copy of
  the file in `data/` yourself.
- **`SSL: SSLV3_ALERT_HANDSHAKE_FAILURE`.** Some servers, miRTarBase among them,
  refuse the TLS handshake from Python and answer `curl`. `setup_data.py`
  retries that download with `curl` when it is on `PATH`, and says so if that
  fails too.
- **`403` after a few requests to one host.** miRTarBase answers 403 for a while
  after several requests in a short time. Wait ten minutes and rerun. The
  download starts again from the beginning of the file.
- **A link has moved.** Releases change. Use the source page in
  [Getting the analysis tables](datasets.md), and keep the file name the loaders
  expect.

## Terms of use

Several sources set their own terms. CTD, ReMap, TRRUST and CLUE restrict
commercial use, TargetScan, miRTarBase and ENCODE each publish their own
conditions on reuse and attribution, and others require attribution. Read the terms at each source
before you share a `data/` folder or a results archive. A derived table such as
`chip_regulon.tsv` or `bindingdb_distilled.tsv` carries the terms of its inputs,
and `results.json` embeds values taken from these tables.

This repository is licensed CC BY 4.0, and it tracks only one file under
`data/`, the small tissue lookup `hpa.tsv`. Six sources state their terms on a
page that was read on 2026-10-23. The others were not read for this page.

| Source | What its page says | Page |
| --- | --- | --- |
| JASPAR | licensed under CC BY 4.0 | jaspar.elixir.no/about |
| openFDA | public domain with a CC0 1.0 dedication, and a warning that some third party works may not be | open.fda.gov/terms |
| ClinGen | curated content free of restriction under CC0 1.0; attribution and the access date are requested | clinicalgenome.org/docs/terms-of-use |
| ChEMBL | CC BY-SA 3.0, so a derived table must carry the same terms | ebi.ac.uk/chembl |
| Human Protein Atlas | CC BY 4.0 for the copyrightable parts; data from other sources may carry its own terms | proteinatlas.org/about/licence |
| TRRUST | freely available for non-commercial research; the site is CC BY-SA 4.0 | grnpedia.org/trrust |

Cite the source when you share anything built from these tables. The run
records the hash of every file passed by flag, which lets a reader find the
release it used.

Record the release, download date and hash of each file you use. The run itself
records the size, modification time and SHA-256 of every file you pass by flag.
