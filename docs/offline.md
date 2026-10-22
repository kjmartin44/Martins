# Running with no network

[Documentation index](index.md)

Everything this package does from a web service it can also do from a
file. This is the list of files, where each one comes from, where to put
it, and what to do to it first.

Most of them need no conversion at all: download, drop under `data/`,
done. The three that do are marked and their commands are given.

## Two different ways to be offline

**A warm cache** replays what a previous run fetched. `http.Client`
stores every response under `~/.cache/generoutes` (or `$GENEROUTES_CACHE`),
so a gene that has been run once will run again with `--offline`. It is
the easiest route and the most fragile: it covers exactly the genes and
factors you happened to ask for.

**Local files** are the real answer. The package looks for them by name
under `data/` and uses them in preference to the network, so a run is
reproducible, fast, and does not depend on what a service returns today.

The two compose. Nothing stops you using files for the bulk sources and
a warm cache for the rest.

## What to get, in order of what it buys

You do not need all of this. The first three change the answer most; the
rest fill in specific sections and the report says plainly when one is
missing.

| Rank | File | Buys |
| --- | --- | --- |
| 1 | ChIP peaks | measured occupancy — the largest single gain |
| 2 | JASPAR motifs | 879 motifs instead of the small bundled set |
| 3 | a drug table | the screening library, and targets per drug |
| 4 | a regulon table | how many genes each factor controls |
| 5 | HPA expression | whether a factor is expressed where your gene is |
| 6 | GENCODE GTF | TSS lookup and the neighbour check, with no network |
| 7 | BindingDB | measured potency, and known ligands for docking |
| 8 | the rest | LINCS, DepMap, PRISM, Orange Book, GWAS, GOA, HGNC |

## The files

Put each one under `data/` with the name in the second column. The names
are what the loaders search for; alternatives are listed where several
are accepted.

This section is the short form, ordered by what each file buys. [Getting
the analysis tables](datasets.md) is the complete per-table reference —
every accepted filename, the exact source, the size, the build step
where there is one, and what the run loses without it. It is generated
from the loaders themselves, so it cannot drift from what the code
searches for. To find out what this checkout is missing right now:

```sh
python -m generoutes.datasets --missing
```

### Transcription factor binding

**JASPAR motifs** — `data/jaspar2024.txt`

```sh
curl -L -o data/jaspar2024.txt \
  "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
```

The bulk PFM file in JASPAR format. Passed with `--pfm-file`. No
conversion.

**ChIP peaks** — `data/remap2022.bed.gz`

ReMap 2022, all human TF ChIP-seq, non-redundant peaks:
<https://remap.univ-amu.fr/> → `remap2022_nr_macs2_hg38_v1_0.bed.gz`.
ChIP-Atlas (<https://chip-atlas.org/>) also works. Passed with
`--chip-bed`. Gzip is read directly; no conversion.

This is the single most valuable file. Without it, factor ranking rests
on sequence alone and most of a motif hit list is chance.

### Which genes each factor controls

**Curated regulon** — `data/trrust.tsv` or `data/collectri.tsv`

TRRUST: <https://www.grnpedia.org/trrust/> → `trrust_rawdata.human.tsv`.
Four columns: TF, target, mode, PubMed id. CollecTRI is the better table
— it signs every edge and covers about three times as many factors —
and is what `--tf-activity` wants.

**ChIP-derived regulon** — `data/chip_regulon.tsv` *(needs building)*

For the factors curation never reached. Built from the ReMap BED and a
GTF:

```sh
python -m generoutes.chipregulon \
    --bed data/remap2022.bed.gz --gtf data/gencode.basic.gtf.gz \
    --out data/chip_regulon.tsv
```

Takes a few minutes and about 0.05 GB of memory. The result is occupancy,
not curation, and the report says so: a peak near a gene is binding, not
regulation.

### Gene coordinates

**GENCODE GTF** — `data/gencode.basic.gtf.gz`

<https://www.gencodegenes.org/human/> → the "basic" annotation GTF. Gives
transcript lookup and the neighbour check with no network at all.

**Genome FASTA** — `data/hg38.fa` *(optional)*

Only needed to pull promoter sequence without Ensembl. Any hg38 primary
assembly FASTA.

### Expression

**Human Protein Atlas** — `data/rna_tissue_consensus.tsv`

<https://www.proteinatlas.org/about/download> → "RNA consensus tissue
gene data". Long format: Gene, Gene name, Tissue, nTPM. Passed with
`--expression-file`. A GTEx median-TPM matrix also works.

The tissue name you pass to `--activity-tissue` must be spelled as this
file spells it — "lymph node", "bone marrow", "skeletal muscle".

### Drugs

**Drug table** — `data/repurposing_drugs.txt`

The Broad Repurposing Hub: <https://repurposing.broadinstitute.org/download>.
This one file removes the whole ChEMBL dependency: SMILES, clinical
phase, mechanism of action and target gene symbols. `data/repurposing_hub.tsv`
and `data/approved_library.tsv` are also accepted.

The mechanism column is what lets a docking screen predict whether a hit
is an agonist or an antagonist offline, so prefer a source that has one.
DrugCentral (<https://drugcentral.org/download>) and TTD
(<https://db.idrblab.net/ttd/full-data-download>) are alternatives.

**Measured potency** — `data/bindingdb_distilled.tsv` *(needs building)*

<https://www.bindingdb.org/bind/chemsearch/marvin/Download.jsp> → the
full `BindingDB_All.tsv` release, about 9 GB unpacked across 640 columns.
Seven of them are read, so project it down once:

```sh
python -m generoutes.localpotency --distil data/BindingDB_All.tsv
```

268 MB afterwards, and read in six seconds rather than three minutes.
Check it landed:

```sh
python -m generoutes.localpotency --uniprot P51681
```

Leave `--min-pchembl` alone unless you know why you are changing it. A
row-level threshold drops weak Ki measurements and lets a stronger IC50
win the assay preference, which changes the answers rather than only the
file size.

**FDA Orange Book** — `data/products.txt`

<https://www.fda.gov/drugs/drug-approvals-and-databases/approved-drug-products-therapeutic-equivalence-evaluations-orange-book>
→ the EOBZIP download, `products.txt` inside it. Gives route of
administration, which nothing else here does. Put the zip itself in `data/`
(`data/EOBZIP.zip`) or unzip `products.txt` there; both are read. About 1 MB.
`scripts/run_everything.bat` offers to download it when neither file is present.

### Measured perturbation

**LINCS L1000** — `data/lincs_chem_pert_up.gmt`, `data/lincs_chem_pert_down.gmt`,
`data/lincs_geneinfo.txt`

The two GMT libraries come from Enrichr:

```sh
curl -L -o data/lincs_chem_pert_up.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LINCS_L1000_Chem_Pert_up"
curl -L -o data/lincs_chem_pert_down.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=LINCS_L1000_Chem_Pert_down"
curl -L -o data/lincs_geneinfo.txt \
  "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt"
```

`geneinfo_beta.txt` matters more than it looks: it says which genes L1000
measured directly and which it inferred, and an inferred gene is much
weaker evidence.

**Other corpora** — optional, same GMT shape:

```sh
curl -L -o data/creeds_drug_up.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Drug_Perturbations_from_GEO_up"
curl -L -o data/drugmatrix_raw.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=DrugMatrix"
```

### Dependency and sensitivity

All from <https://depmap.org/portal/data_page/>, current release.

| File | Keep as |
| --- | --- |
| CRISPR gene effect | `data/CRISPRGeneEffect.csv` |
| cell line metadata | `data/Model.csv` |
| common essentials | `data/CRISPRInferredCommonEssentials.csv` |
| PRISM repurposing data matrix | `data/Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` |
| PRISM compound list | `data/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` |
| expression matrix | `data/OmicsExpressionProteinCodingGenesTPMLogp1.csv` |

No conversion. The names above are what the loaders look for first;
`data/crispr_gene_effect.csv` and similar lowercase forms also work.

### Annotation and human genetics

| Source | Keep as | From |
| --- | --- | --- |
| HGNC | `data/hgnc_complete_set.txt` | <https://www.genenames.org/download/archive/> |
| GO annotations | `data/goa_human.gaf.gz` | <https://current.geneontology.org/annotations/goa_human.gaf.gz> |
| GWAS Catalog | `data/gwas_catalog_associations.tsv` | <https://www.ebi.ac.uk/gwas/docs/file-downloads> |

A local GWAS catalogue reports the gene nearest a signal, which is weaker
than the fine-mapped support Open Targets gives. The report distinguishes
them rather than adding the counts together.

### Dosage tolerance

Whether there is room to move the gene, and whether moving it is safe.
Both are read from disk only — there is no API fallback, so without these
files the report says dosage was not assessed.

| Source | Keep as | From |
| --- | --- | --- |
| gnomAD v4.1 constraint | `data/gnomad.v4.1.constraint_metrics.tsv` | <https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv> |
| ClinGen dosage curation | `data/ClinGen_gene_curation_list_GRCh38.tsv` | <https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv> |

`scripts/setup_data.bat` fetches both, as does `scripts/run_everything.bat`,
asking first about gnomAD
because it is 91 MB. ClinGen is 250 KB and re-dated weekly upstream.

Two things about reading them that the files do not announce:

gnomAD is one row per transcript — 211,523 rows for 18,204 genes — and
the gene's value comes from its MANE Select transcript, falling back to
canonical. Taking the lowest LOEUF across transcripts instead looks
neutral and is not: a minimum over a dozen draws lands below the typical
one, so every multi-transcript gene drifts toward "constrained". It moves
9,207 of 18,204 genes, 922 of them by more than 0.2 LOEUF, which is wider
than the gap between this package's own sensitive and tolerant lines.

ClinGen's score is not a severity scale. Two of its values are category
codes sitting above the evidence grade: `30` means the gene has a
recessive phenotype, which says nothing about dose, and `40` means
dosage sensitivity is **unlikely** — the opposite of the top of the
scale. `0` means "no evidence available", which is 1,306 of the genes in
the current release and is read as not assessed rather than as a
negative.

### Docking

A binary rather than a data file. [The docking guide](docking.md) covers
installation; the short version is gnina on Linux or through WSL, smina
everywhere:

```sh
python -m generoutes.docking --check
```

Structures are fetched from RCSB and cached, so the first run for a
target needs a network and later ones do not. `--prepare` and `--bundle`
split the networked and compute halves across machines.

## Checking what you have

Every loader reports what it found. The fastest check is a real run with
`--offline`, which forbids the network entirely, and then reading the
evidence table at the top of the report:

```sh
python -m generoutes --gene CCR5 --direction down --offline \
    --pfm-file data/jaspar2024.txt \
    --chip-bed data/remap2022.bed.gz \
    --regulon-file data/chip_regulon.tsv \
    --expression-file data/rna_tissue_consensus.tsv \
    --gtf data/gencode.basic.gtf.gz \
    --out results/offline_check
```

The table names every source as `ok`, `ok, nothing found`, `FAILED`, or
`not requested`, with what each failure cost. **A source that failed is
recorded as failed, not as zero** — which is the whole point of reading
that table before the ranking underneath it. The second state is the one
that is easy to skim past: the source was reached, it answered, and the
answer was nothing. That is a measurement, and it is not a find.

Individual checks, when you want one:

```sh
python -m generoutes.localpotency              # potency tables
python -m generoutes.docking --check           # binary and rdkit
python -m generoutes.wsldock                   # docking through WSL
python -m generoutes.audit                     # curated tables
```

## What offline cannot give you

**Structures for a target nobody has solved.** RCSB has what it has, and
nothing local substitutes for it.

**Up-to-date approvals.** A drug table is a snapshot. A drug withdrawn
after the file was published is still marked available in it, and that is
the dangerous direction of the error. `python -m generoutes.audit
--online` cross-checks the curated tables against ChEMBL and is worth
running when you next have a network.

**Anything the file does not cover.** A regulon table with 590 factors
cannot say anything about the 591st, and a factor missing from it is
untested rather than clean. The reports distinguish those two
consistently, and it is the distinction most worth keeping in mind when
reading a run assembled from files.
