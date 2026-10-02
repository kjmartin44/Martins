# Inputs and outputs

[Documentation index](index.md)

## Reference inputs

| Input | Format and behavior |
| --- | --- |
| `--pfm-file` | JASPAR text matrices: a `>matrix_id name` header followed by A, C, G, and T count rows. Without a local file, motifs come from JASPAR. |
| `--matrix-id` | Repeat to select matrix IDs. |
| `--chip-bed` | Plain or gzip BED with target information in the name column. See `evidence.py` for supported names. Coordinates must match the promoter assembly. |
| `--regulon-file` | TRRUST TSV, a two-column TF/target TSV, or a CollecTRI or DoRothEA table with source/target/weight columns, comma or tab separated. Used for target counts, and for signed regulons when `--tf-activity` is on. |
| `--lincs-geneinfo` | CLUE geneinfo table. Records whether each gene was measured on the L1000 array or imputed from the 978 landmarks, which sets the evidence tier of a LINCS hit. |
| `--expression-file` | HPA long table with `Gene name`, `Tissue`, and `nTPM`, or a wide table whose first column contains gene symbols and remaining columns contain numeric tissue values. Gzip and ZIP are detected from file content. |
| `--lincs-file` | Comma- or tab-separated table with `drug`, `gene`, and signed `score`; optional `cell_line`, `dose`, `time`, and `is_landmark`. Gzip files use a `.gz` suffix. |
| `--signature-file` | Two columns: gene and score. Currently used only to count reference genes in the LINCS section. |

Raw GTEx files with Ensembl IDs must be converted to gene symbols before using
the wide expression format. Missing landmark status in a signature table is
treated conservatively as inferred. The current signature loader keeps effects
in the requested direction with absolute score at least 1.5.

`data/jaspar2024.txt` and `data/trrust.tsv` are reference files in this checkout.
They are not included in the installed wheel; pass paths to your own files or
use the network defaults. **`data/hpa.tsv` contains `Tissue` and `Organ` columns.
It is a tissue lookup and cannot be used as an expression matrix.**

Record the source release, download date, assembly, and file hash when obtaining
new datasets. The repository does not currently provide a complete provenance
manifest or redistribution-license inventory for these files.

## Regulon tables for TF activity

`--tf-activity` infers which factors are active from whether their target
genes moved. That needs a signed regulon, meaning each edge is labelled as
activation or repression. Two tables work, and they are not equivalent.

TRRUST is in this checkout at `data/trrust.tsv`. It covers 329 factors and
labels 47 percent of its human edges; the rest are recorded as `Unknown`.
An unlabelled edge still supports the two-tailed test that asks whether the
regulon moved, but it cannot say which way the factor went.

CollecTRI is the better input. It signs every edge and covers roughly three
times as many factors, which is the main limit on how many can be tested at
all. Fetch it once:

```python
import decoupler as dc
dc.get_collectri(organism="human").to_csv(
    "data/collectri.tsv", sep="\t", index=False)
```

Anything at `data/collectri.tsv` or `data/collectri.csv` is picked up
automatically, ahead of TRRUST. Pass `--regulon-file` to override. The
separator is read from the first line rather than the file name, so either
suffix works.

## LINCS signature files for the drug scan

The drug scan asks which compounds move a factor's regulon, so it needs
both Enrichr LINCS gene set libraries present:

```
data/lincs_chem_pert_up.gmt     16 MB
data/lincs_chem_pert_down.gmt   14 MB
```

Both are picked up automatically when they exist, and the scan is skipped
without comment when they do not. Each factor costs about a second, so the
scan covers the top eight ranked factors rather than all of them.

## Drug descriptions

Report tables name a compound and its mechanism against one receptor, which
is not enough to judge whether a row is applicable. `indications.py` adds
what the molecule is and what it is approved to treat, from three sources,
and every entry records which one it used.

| source | what it is | limits |
| --- | --- | --- |
| `curated` | A short hand-written entry for the compounds these reports return. | Exact, and incomplete by construction. |
| `ChEMBL` | `/drug_indication`, taken at the highest phase reached. Needs the network; cached like every request. | The list includes uses studied but never approved. Approved entries are preferred and studied-only ones are labelled. |
| `class` | The systemic class from `systemic.py`. | Gives the kind of drug, not its indication. |

None of this is a label. Check the prescribing information before acting on
a row. Every route in a report is off-label by construction, so the approved
indication is the baseline the proposal departs from rather than a
recommendation.

## DepMap dependency data

`--depmap-file` reads a DepMap CRISPR knockout release and reports what a
cell loses without each candidate factor. The regulon count answers the
same collateral question from literature curation; this answers it by
measurement, and it covers the factors no curator has written up.

Download `CRISPRGeneEffect.csv` from the current release at
https://depmap.org/portal/data_page/ and put it in `data/`, where it is
picked up automatically. `Model.csv` from the same download adds the
lineage column and is optional; it is looked for beside the effect file
first, then in `data/`.

The effect file is cell lines by genes and around 500 MB. It is streamed
once and only the requested columns are kept, so the size costs a pass over
the file rather than memory. Column names carry the Entrez id, as in
`TP53 (7157)`, and the loader strips it.

Scores are on DepMap's Chronos scale. Zero means knockout changed nothing
and -1 is the median of known common-essential genes. A cell line counts as
dependent below -0.5, and a gene is pan-essential when at least 90 percent
of lines are.

Three limits decide whether any of it transfers. The panel is cancer cell
lines dividing in vitro, so a factor whose role is in a post-mitotic or
differentiated cell looks dispensable for reasons unrelated to its
importance. Knockout is not inhibition: the protein is removed entirely,
including scaffolding roles no ligand would touch, and essential genes are
drugged successfully at partial occupancy. And nothing here concerns your
gene, so it bounds the cost of acting on a factor and says nothing about
the benefit.

## Docking

Docking is opt-in and separate from the main pipeline, because it needs an
external binary and fetched structures. Check whether your machine can run
it:

```
python -m promoter_tf_drugs.docking --check
```

It reports what is present and what to install. Two things are needed: a
docking binary on PATH, and RDKit.

smina is the binary to get. It is a single static file and takes a PDB
receptor with an SDF ligand directly, which is why it is preferred here.
On Linux, WSL or macOS:

```
curl -L -o smina https://sourceforge.net/projects/smina/files/smina.static/download
chmod +x smina && sudo mv smina /usr/local/bin/
```

On Windows, conda-forge ships smina for win-64, so no WSL is needed:

```
conda install -c conda-forge smina
```

An earlier version of this page said smina had no Windows build. That was
wrong.

AutoDock Vina is also supported for cases where you have it instead. It
needs PDBQT for both receptor and ligand, which `pip install meeko` and
`pip install openbabel-wheel` provide; the receptor is converted once per
screen and reused. There is no reason to prefer that route when smina
installs in one command.

A screen is thousands of independent docking jobs, so `--cpu` is the
largest speed lever available. Measured at 3.7 seconds per compound on
eight cores, which puts a 600-compound library at roughly ten minutes on
forty cores and several hours on one.

A GPU does not help much here. The GPU-accelerated engines are
Linux-only (unidock on conda-forge is linux-64; gnina and Vina-GPU are not
packaged for conda at all), and with enough CPU cores the library finishes
in minutes anyway. The case for a GPU would be DiffDock, which is a
different method rather than a faster one, and whose confidence score
rates the predicted pose rather than binding strength.

Then run a screen:

```
python -m promoter_tf_drugs.docking --gene VDR --state active --top 25
```

### Running the networked and Linux halves separately

The network and the docking binary are often not on the same machine: on
Windows the network is there and smina is not, in WSL the reverse. The two
halves are separable, so nothing has to be online and Linux at once.

Step one, where the network is. Fetches the structure, the receptor's
known ligands and the approved library into one folder. No docking binary
needed.

```
python -m promoter_tf_drugs.docking --gene VDR --state active --prepare dock_vdr
```

Step two, where smina is, with no network at all. WSL sees Windows drives
under `/mnt/c`, so the folder does not need copying.

```
python -m promoter_tf_drugs.docking --bundle /mnt/c/Users/you/dock_vdr --out dock_vdr/results.json
```

Step three, back where the report is built. The hits join the route tables
as tier 4, `TF: predicted pocket binder`.

```
promoter-tf-drugs --gene YOURGENE --direction up --docking-results dock_vdr/results.json
```

A run whose pocket failed validation contributes nothing, so a failed dock
cannot quietly add rows. `--pdb FILE` supplies your own structure instead
of searching RCSB; it must carry a co-crystallised ligand, which is what
places the box.

What happens, in order. The gene symbol resolves to a UniProt accession,
RCSB is searched for structures, and one with a co-crystallised ligand is
chosen; a structure without one gives no reliable box and is skipped, and
if none has one the run is refused rather than docking into a guessed
site. The chosen structure is stripped of its ligand, waters and ions,
because docking into an occupied pocket scores every compound weakly for
reasons unconnected to the compounds. The box is placed on the removed
ligand.

Then the gate. The receptor's own measured ligands are docked alongside
decoys, and the ranking is scored by ROC AUC and by enrichment in the top
1 percent. Both must clear a bar before any score on an untested compound
is reported. Below it, nothing is reported and the reason is given.

`--state active` restricts to agonist-bound conformations, read from the
file by detecting a bound LXXLL coactivator peptide. That matters when you
want a compound that raises your gene: an agonist has to fit the active
conformation.

Costs roughly four seconds per compound per CPU core, so `--library-limit`
defaults to 600 rather than the whole approved set. Structures are cached
under `--work-dir`.

Two limits worth knowing before investing the time. Docking needs a pocket
and a structure, so it does not reach transcription factors that have
neither, which is most of them; the receptors it works on are the
nuclear-receptor ligand-binding domains already prominent in the report.
And a docking score is not an affinity: these functions enrich for binders
and rank them poorly, so the output is a shortlist to assay rather than a
predicted order.

## Perturbation corpora

The drug scan asks whether a compound's own expression signature moves a
factor's regulon. Run against LINCS alone that is one platform, one set of
mostly cancer cell lines, one dose regime and one pipeline. Every
systematic error of that kind is shared by all 33,000 LINCS signatures, so
no amount of internal replication detects it. A second corpus built
differently is what does.

Check what is on disk:

```
python -m promoter_tf_drugs.corpora
```

| corpus | species | system | the main limit |
| --- | --- | --- | --- |
| LINCS L1000 | human | cancer cell lines, hours at micromolar | 978 genes measured, the rest inferred |
| CREEDS | human and mouse | published GEO experiments | heterogeneous quality, few replicates per signature |
| DrugMatrix | rat | tissue dosed in vivo, mostly liver | rat symbols, liver-dominated, toxicology doses |
| Open TG-GATEs | rat and human | rat in vivo, rat and human hepatocytes | liver only, about 170 compounds |

Each is a pair of Enrichr gene-set libraries in the same `.gmt` format the
LINCS files already use, so nothing new has to be parsed. Fetch what you
want into `data/`:

```
curl -L -o data/creeds_drug_up.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Drug_Perturbations_from_GEO_up"
curl -L -o data/creeds_drug_down.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Drug_Perturbations_from_GEO_down"
```

DrugMatrix ships as a single library with the direction in the term name,
so it needs splitting first:

```
curl -L -o data/drugmatrix_raw.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=DrugMatrix"
python -m promoter_tf_drugs.corpora --split data/drugmatrix_raw.gmt --as drugmatrix
```

A term carrying no direction marker cannot be assigned and is dropped
rather than guessed onto one side; the count of dropped terms is reported.

TG-GATEs is not on Enrichr. Open TG-GATEs publishes the raw arrays at
https://dbarchive.biosciencedbc.jp/en/open-tggates/ and a derived up/down
pair placed at `data/tggates_{up,down}.gmt` is picked up automatically.

Nothing is required. A corpus that is absent is skipped and reported as
absent, which is a different statement from a corpus that ran and found
nothing.

### What agreement across corpora is worth

A compound that shifts a factor's regulon in LINCS and in DrugMatrix has
done so in a human cell line and in rat liver, on different platforms, in
different laboratories. They share few systematic errors, so agreement is
worth more than a stronger result in either alone. That is the headline
the report leads with.

Disagreement usually means nothing. These corpora differ in species,
tissue, dose, duration and platform all at once, so a compound active in
one and not another may simply not have been tested comparably. Conflicting
directions are shown rather than hidden, because suppressing them would
make the agreeing rows look more uniform than the data is.

ARCHS4 is deliberately not offered here. It is uniformly reprocessed
RNA-seq rather than perturbation data, and its useful output is
co-expression — the genes that move with a factor across thousands of
samples. That is a regulon, not a drug effect, so listing it beside LINCS
would invite reading one as the other.

None of this fixes attribution. Every corpus measures the regulon, so all
of them share the same confound: a factor that shares targets with another
cannot be separated from it, and a compound that arrests the cell cycle
moves every proliferation-linked regulon. More corpora make the observation
more robust, not more specific.

## Running without a usable network

Several sources are web APIs, and on a corporate network they are often
blocked or intermittent. Two things make that survivable.

The HTTP client gives each host one full retry budget and then treats it as
down for the rest of the run, so a dead endpoint costs one request's worth
of waiting rather than one per call. ChEMBL is called hundreds of times in a
run, so without that a 500 would dominate the wall clock. Everything not
depending on that host is unaffected.

That makes a broken endpoint cheap. It does not make the analysis possible,
and for that the data has to be on disk. The docking and screening paths
need three things from ChEMBL - the gene's target id, the measured ligands
for the validation gate, and the approved-drug library - and one annotated
file supplies all three.

```
# Broad Drug Repurposing Hub: https://repurposing.broadinstitute.org/download
# save the drugs table as data/repurposing_drugs.txt
python -m promoter_tf_drugs.docking --gene VDR --pdb 1DB1.pdb --library data/repurposing_drugs.txt
```

Around 6,800 compounds with SMILES, clinical phase, mechanism and target
gene symbols, in a few megabytes. `data/repurposing_drugs.txt` is picked up
automatically; any table with a name column and a SMILES column works, and
the delimiter is read from the header rather than the file name because the
Hub ships tab-separated `.txt` that people save as `.csv`.

The local path is tried first and the network only as a fallback, because on
a blocked network the API calls are the step that stops everything.
`--actives FILE` overrides the validation set separately when you have a
better one than the table's own annotation.

Verified end to end with every network entry point replaced by a function
that raises: a local table, a local structure, smina, AUC 1.0 on the
vitamin D ligands, no network touched.

What a local table cannot supply is potency. The Hub records that a compound
acts on a target, not how hard it binds, so an `actives` set from it carries
no pChEMBL and anything ranking on potency degrades to unranked rather than
scoring zero. The docking gate needs only the identity of the known binders,
so docking is unaffected. It is also a snapshot: a drug approved after
publication is absent, and absence is not evidence that nothing exists.

## Output files

Every completed analysis writes six files:

| File | Contents |
| --- | --- |
| `summary.md` | One page: the shortlist, what the evidence is worth, and what would change the answer. Start here. |
| `report.md` | The full analysis. Every number qualified, every stage accounted for. |
| `results.json` | Query metadata, scoring weights, direct therapies, and per-promoter analysis blocks. |
| `tf_summary.tsv` | TF rows with motif, role, ChIP, score, and drug fields. |
| `sites.bed` | Predicted motif sites in genomic coordinates with a UCSC track header. |
| `promoter.fa` | Promoter sequences oriented 5′ to 3′ along the gene. |

`summary.md` carries nothing that `report.md` does not also say. It exists
because `report.md` is long on purpose: every number in it is qualified,
since an unqualified number from this pipeline reads as stronger than it
is. That makes it the wrong document to open when you want to know what to
do next. The summary leads with how much the ranking is worth, because in
a degraded run the order means little and burying that under the table
would be the one mistake the page could make. It never drops a caveat to
make a row fit.

`results.json` has top-level `query`, `results`, `direct_therapies`, `weights`,
and `drug_weights`. Each entry in `results` includes `promoter`, `sequence`,
`hits`, `tfs`, `drug_grades`, `all_routes`, and optional stage outputs. Empty
sections can mean no result, a disabled stage, or a service failure; consult the
run log. There is no versioned JSON schema yet.

The TF TSV combines promoter blocks without transcript identifiers. Use the
JSON when distinguishing alternative promoters. The BED contains a `track`
header; remove that header before passing it to tools that accept data rows only.

## Coordinates

Promoter `start`, `end`, and `tss` use 1-based genomic coordinates with inclusive
ends. Hit `start` and `end` index the gene-oriented sequence from zero with an
exclusive end. A hit's JSON strand is relative to that sequence.

BED uses zero-based genomic start and exclusive end, as specified by
[UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). Its strand is genomic.
For a negative-strand promoter, the exporter reverses the motif strand and maps
the interval back to the forward genomic coordinate system. For example, on a
promoter spanning bases 101–120, sequence interval `[2, 7)` maps to BED
`[102, 107)` on a positive-strand gene or `[113, 118)` on a negative-strand gene.

## HTTP cache

Responses are gzip JSON files under `~/.cache/promoter_tf_drugs` by default.
The online TTL is 30 days. Offline mode accepts older entries. Keys include URL,
method, body, parsing mode, and supplied headers. This prevents a text response
from being reused as parsed JSON or a GET from being reused as a POST.

**Compatibility:** cache keys changed during the reliability update. Old cache
files remain on disk but are not reused; populate the new cache online before
expecting an offline run to work. The package does not delete old entries.

Concurrent writers use separate temporary files and atomic replacement. A cache
write failure is logged and does not discard a successful HTTP response.
Service requests retry transient failures, honor valid `Retry-After` seconds or
dates, and bound manually followed POST redirects. No authenticated service
configuration is provided by the CLI.

## Running without the web services

Ensembl, ChEMBL and the maayanlab hosts all return HTTP 500 while merely
busy, so a 500 is treated as a backoff signal rather than a verdict: six
attempts with capped, jittered backoff, tunable with `PTD_RETRIES` and
`PTD_RETRY_CAP`. `Retry-After` is honoured in both its seconds and HTTP-date
forms.

Every stage has a local-file route. Ensembl sequence was the last one
without, so a 500 there used to take the whole run down — but the package
already *writes* every window it fetched to `promoter.fa`, which makes a
previous run its own offline source:

```bash
python -m promoter_tf_drugs --gene YOURGENE     --promoter-fasta results\YOURGENE_up\promoter.fa ...
```

Matching is on exact coordinates, so a window is reused only when it is
genuinely the same window; a near-miss is rejected rather than silently
analysing the wrong sequence. Any FASTA whose headers carry a `chrN:start-end`
span works.

With `--all-tss`, one unreachable window no longer aborts the run: it is
skipped with a warning and the rest continue, with the lost windows reported.
Only a run where *every* window failed raises, because an empty promoter list
would read downstream as "this gene has no promoter".
