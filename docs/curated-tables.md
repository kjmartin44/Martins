# Curated tables

[Documentation index](index.md)

Some of what the pipeline knows was typed into the source and not read from a
database. These are the curated tables. This page lists each one, says what a
complete version would be, and gives the steps to rebuild or check it.

## Read this first

The repository does not record who entered each entry or from what source.
Treat an entry as unreferenced until someone has checked it against the source
named here. The checks below show that names exist, that gene symbols are
current, and that a few tables agree with a public list. They do not show that
an entry is right, and they cannot find what is missing from a table that no
public list defines.

"Complete" means one of three things, and each table below says which:

- **Defined by a public list.** The table should equal a list that someone
  else maintains. The audit compares the two and names every difference.
- **Defined by a rule.** The table holds the answer to a question, such as
  "approved drugs that bind tubulin". A public database can answer the
  question, but nothing runs it for you. The steps below show how.
- **Defined by judgment.** A person chose the entries for a purpose, and no
  public list says what the right set is. The table is complete only in the
  sense that its stated rule was applied. Add an entry only with a reference.

A missing entry is never evidence of absence. The tables' own notes say the
same: a drug or gene a table does not name is unassessed, not cleared.

## Check the tables

Run these from the checkout root. The checks that read a downloaded file look
for it under `data/` relative to the working directory.

```sh
python -m generoutes.audit
python -m generoutes.audit --online
```

The first needs no network. It checks that no drug sits in two classes, that
names survive salt stripping, that text fields are not too thin, and that
`KNOWN_LANDMARKS` and the Orange Book route tables agree with the CLUE gene
table and the Orange Book when those files are on disk. A check whose file is
absent says so and is not counted as a pass.

The second adds four checks against public services. The first run makes
about 600 ChEMBL requests, 9 HGNC searches, 4 PDB requests and one HGNC
download. The HTTP cache keeps them for 30 days, so a repeat run is quick.
`--limit N` checks only N names per table for a quick look.

| Check | Tables | Service |
| --- | --- | --- |
| Every drug name resolves, and ChEMBL's withdrawn flag | 13 tables, 604 names | ChEMBL |
| Every gene symbol is an approved symbol | 19 lists, 333 symbols | HGNC search |
| The nuclear receptor list equals the HGNC group | `NUCLEAR_RECEPTORS` | HGNC group download |
| Residue codes exist and have the right type | 4 tables, 177 codes | PDB chemical component dictionary |

The exit status is 1 when a check finds a problem and 0 otherwise. A check
that could not run prints `UNKNOWN, not confirmed` and does not change the
status.

A name that resolves in ChEMBL has matched something. The search returns its
best match, so resolution shows the name exists and not that it is the
molecule the table means.

### Result of the last run, 2026-10-23

| Check | Result |
| --- | --- |
| Internal consistency | no problems |
| Drug names | 604 of 604 resolved |
| Withdrawn in ChEMBL | clofibrate, etretinate, rosiglitazone (burden class), etretinate and rosiglitazone (indications), berberine (promoter DNA agents, already marked not approved), rosiglitazone (positive controls) |
| Our `withdrawn` flag against ChEMBL | no entry we present as available is withdrawn; we flag mipomersen and ChEMBL does not |
| Nuclear receptors | 48 of 48 match HGNC |
| Gene symbols | 2 stale symbols found and fixed (see below); the rest are approved |
| Residue codes | all known and correctly typed |
| `KNOWN_LANDMARKS` | 11 of 32 genes were not landmarks; removed |
| Orange Book routes | 28 routes in neither table, left unclassified on purpose |

## Tables defined by a public list

### Nuclear receptors

`annotate.NUCLEAR_RECEPTORS`, 48 genes. `localannot.py` imports it.

The reference is the HGNC gene group "Nuclear hormone receptors" (group 71).
Rebuild the set from the protein coding members, which leaves out the
pseudogene NR1H5P:

```python
import csv, io, urllib.request
url = "https://www.genenames.org/cgi-bin/genegroup/download?id=71&type=branch"
text = urllib.request.urlopen(url, timeout=60).read().decode("utf-8")
rows = csv.DictReader(io.StringIO(text), delimiter="\t")
print(sorted(r["Approved symbol"] for r in rows
             if r["Locus type"] == "gene with protein product"))
```

Paste the printed list into `annotate.py`. `python -m generoutes.audit
--online` then reports any difference in either direction.

### Residue codes for docking

`docking.MODIFIED_RESIDUES` (50), `docking.SOLVENTS` (48),
`docking.CRYSTALLISATION_ADDITIVES` (50) and `exploratory.SUGARS` (29).

The reference is the PDB chemical component dictionary, served by
`https://data.rcsb.org/graphql`. A modified residue is a peptide linking
component that a structure writes as HETATM. A solvent or additive is a
component that is not peptide linking. A sugar is a saccharide. The audit
checks every code against that type.

The dictionary holds far more modified residues than the table, so this table
is not complete in the sense of covering every component. It lists the ones
met in practice. Add a code when a screen picks something odd as its pocket
ligand, as happened with sulfotyrosine in CCR5 structures. `WAT` and five
disaccharide codes (`CBK`, `LAK`, `MAL`, `SUC`, `TRE`) are retired from the
dictionary but still appear in older structures, so the audit allows them.

### Landmark fallback

`signatures.KNOWN_LANDMARKS`, 21 genes.

The reference is the CLUE gene table, `data/lincs_geneinfo.txt`, where the
column `feature_space` is `landmark` for the 978 measured genes. Every name
must be a landmark there. The audit checks this when the file is present.
The list is a fallback for runs without that file, so keep it short.

### Orange Book routes

`orangebook.SYSTEMIC_ROUTES` (13) and `orangebook.LOCAL_ONLY_ROUTES` (17).

The reference is the route field of the Orange Book `products.txt`. The audit
prints the routes the file uses that neither table classifies, and the table
entries the file no longer uses. Today 28 routes are in neither table, for
example `INTRAVESICAL`, `SPINAL` and `INTRADERMAL`. A drug whose routes are
all unclassified is reported as "unclear", which is the intended answer when
nobody has decided.

Two entries, `PARENTERAL` and `RESPIRATORY (INHALATION)`, are not used by
the current file. They come from older releases and stay.

One placement is open for the project owner to decide. `BUCCAL`, `SUBLINGUAL`,
`RECTAL` and `INHALATION` are listed as local only, and the default report
hides a drug approved only by those routes. Some drugs reach the bloodstream
this way, for example sublingual buprenorphine and rectal diazepam. This page
did not change them, because that would change which drugs the reports show.

## Tables defined by a rule

### Drug names

Each table below holds drug names chosen by a pharmacological rule. The audit
resolves every name in ChEMBL, which shows the name still exists and whether
ChEMBL calls it withdrawn. To rebuild a table from a public source, ask
ChEMBL the rule and compare:

1. Find the target with `/target/search` by gene symbol, or take the drug
   list for a class from `/atc_class` or `/mechanism`.
2. List approved molecules (`max_phase` 4) for that target or class with
   `/mechanism?target_chembl_id=...`.
3. Check each result against the table, and add what is missing with the
   reference in a comment.

| Table | Entries | Rule | Where to ask |
| --- | --- | --- | --- |
| `systemic.CLASSES` | 11 classes, 257 names | Approved drugs whose class acts on the whole body, such as glucocorticoids, retinoids and sex steroids. Each class has a burden score between 0 and 1. | ChEMBL `/mechanism` for the class target, and the ATC class |
| `systemic.MITOTIC_DISRUPTORS` | 45 | Compounds that bind tubulin or inhibit PLK1, Aurora kinases, Eg5, CDK1 or MPS1. A tag, so a drug may also sit in a class above. | ChEMBL `/mechanism` for those targets |
| `systemic.GENOTOXIC` | 46 | Approved drugs whose mechanism is DNA damage and whose label carries a mutagenicity or secondary malignancy warning. | Drug labels through openFDA (`api.fda.gov/drug/label.json`) |
| `delivery.PGP_SUBSTRATES` | 49 | P-glycoprotein substrates, lower case. | The FDA table of transporter substrates; this page did not compare the two |
| `indications.CURATED` | 74 | What the molecule is and what it treats, for compounds the reports return. | ChEMBL `/drug_indication`, highest phase first |
| `standard_of_care.CONDITIONS` | 1 condition, 3 drugs | Approved treatments for the condition a gene causes. | Drugs@FDA through openFDA |
| `catalogue.CHROMATIN_DRUGS` | 8 | Approved HDAC, DNMT and EZH2 inhibitors, with approval year. | Drugs@FDA; the years were not checked |
| `retrospective.POSITIVE_CONTROLS` | 6 genes, 32 names | A drug that acts on a factor known to drive the gene. | The primary literature |

`indications.CURATED` is deliberately short. Of 403 distinct drug names in the
other curated drug tables, 63 have an entry. The rest are described from
ChEMBL when the network is on, and from the class text when it is not. A
report states which source each description used. To see the gap today:

```python
from generoutes import audit, indications
from generoutes.congeners import salt_free_name
names = {salt_free_name(n).lower()
         for k, v in audit.all_curated_names().items()
         if k != "indications.CURATED" for n in v}
print(len(names - set(indications.CURATED)), "names have no curated entry")
```

The approval dates in `standard_of_care.py` were supplied by the project
owner. They were checked on 2026-10-23 against Drugs@FDA through openFDA, with
one request for each brand:

```sh
curl -s "https://api.fda.gov/drug/drugsfda.json?search=products.brand_name:%22AUSTEDO%22&limit=5"
```

Xenazine (NDA 021894) was approved on 2008-08-15, Austedo (NDA 208082) on
2017-04-03, and Ingrezza (NDA 209241) first on 2017-04-11 for tardive
dyskinesia, with an efficacy supplement on 2023-08-18. The record does not
name the indication of that supplement, so the match to Huntington's chorea
rests on the date. Only Huntington disease has an entry. A gene with no entry
has no standard of care recorded, which does not mean it has none. Adding one
means taking the indication from the label's "Indications and Usage" section
and the date from Drugs@FDA, and writing both down.

## Gene symbol lists

The audit checks that every symbol in these lists is an HGNC approved symbol.
It found two stale symbols, now fixed:

| Symbol found | Problem | Now |
| --- | --- | --- |
| `TAZ` in `systemic.ONCOGENES` | The old symbol of TAFAZZIN, and an alias of WWTR1. The factor meant is WWTR1. | `WWTR1` |
| `SLC22A18` in `reactivation.IMPRINTED` | Renamed by HGNC. | `SLC67A1` |

Two entries are not gene symbols and stay: `SNORD116` names a cluster of
snoRNA genes, and `NESP` is a transcript of GNAS, which is listed beside it.
They are named in `audit._NOT_SYMBOLS`.

| List | Entries | Rule | Reference |
| --- | --- | --- | --- |
| `systemic.ONCOGENES`, `TUMOUR_SUPPRESSORS` | 63, 37 | Genes whose role in cancer is not in dispute, so that moving them the wrong way is flagged. Short on purpose. | The COSMIC Cancer Gene Census, which needs registration, or a published driver gene list |
| `cellcycle.MITOTIC_GENES` | 63 | Canonical G2/M and checkpoint genes. A data driven version was tried and could not work. | GO terms for mitosis give support but a far longer list |
| `proliferation.MARKERS` | 26 | The well measured core of the G2/M programme, mostly L1000 landmarks. | The same |
| `reactivation.IMPRINTED` | 26 | Human imprinted loci worth an allele specific question. | A published imprinting catalogue |
| `regulon` tiers | 29, 29 and 23 | Coarse breadth of well known factors, used only when no regulon file is given. | Target counts per factor in `data/collectri.tsv` or `data/chip_regulon.tsv` |
| `cellcycle.PROFILES` | 31 | Cell cycle phase dependence of factors. | The primary literature |

Changing a membership list changes which genes get a caution, so it changes a
report. Record the reference in a comment on the same line.

## Tables defined by judgment

These hold a person's choices. The audit resolves their drug names and checks
their gene symbols, but nothing public says what the right set is. Each entry
should name its source in a comment. None does today.

| Table | Entries | What an entry claims |
| --- | --- | --- |
| `stabilizers.CURATED`, `EXEMPLARS`, `DEGRON`, `DUBS` | 11 factors, 4, 7, 3 | A drug stabilises, or a ligase and a deubiquitinase act on, a named protein. A wrong pair sends someone after the wrong target, which is why `DUBS` is short. |
| `ligands.ROUTES` | 8 receptors, 27 drugs | A drug changes how much of a receptor's own ligand is made, broken down or carried. |
| `direct.CURATED` | 17 genes, 29 entries | A drug or oligonucleotide acts on the gene's own product or transcript. The `withdrawn` flag is checked against ChEMBL. |
| `paralogs.FAMILIES` | 21 families | Receptors that share a binding pocket. HGNC gene groups give a partial check. |
| `modulators` | 8 mechanism profiles | What each mechanism class does and costs. |
| `promoterdna.AGENTS` | 12 agents | A compound has a documented sequence preference. Membership is not a recommendation. |
| `congeners.SALTS`, `ESTERS` | 78, 29 | Words that mark a salt or ester of the same active molecule. |
| `modalities.DELIVERY` | 9 tissues | How an oligonucleotide or gene therapy can be delivered to each tissue. |

Weights, labels and score vocabularies (`pipeline.WEIGHTS`, `drug_scoring`,
`provenance`) are design choices. [Methods and limitations](methods.md) lists
them. `demo_data.py` holds synthetic fixtures for `--demo` and is not evidence.

## Adding or changing an entry

1. Find the source: a database record, a label, or a paper. Write its name and
   identifier in a comment on the same line.
2. Check the name or symbol with `python -m generoutes.audit --online`.
3. Run `python -m pytest -q`. Tests that pin a table will fail if an entry
   moves a result they expect.
4. If a report would change, add a compatibility note to
   [Inputs and outputs](data.md). Say which routes or cautions move.

## What changed on 2026-10-23

These edits fix errors the checks found. They change a few results.

- `systemic.ONCOGENES` lists `WWTR1` in place of `TAZ`. TAZ is an alias of
  WWTR1 but also the old symbol of TAFAZZIN, so the entry matched only a
  source that used the alias. A route that raises WWTR1 now carries the
  proto-oncogene caution.
- `reactivation.IMPRINTED` lists `SLC67A1` in place of `SLC22A18`. A run for
  that gene now gets the imprinting question. Run it under the new symbol.
- `signatures.KNOWN_LANDMARKS` lost `DDIT3`, `ATF3`, `SOD1`, `NQO1`, `CEBPB`,
  `SREBF1`, `VDR`, `ESR1`, `AR`, `E2F1` and `MKI67`, which the CLUE table
  calls "best inferred". Without a CLUE gene table, a LINCS hit at one of those
  genes is now weighted as inferred where it was weighted as measured. Runs
  with `data/lincs_geneinfo.txt` read the file and are unchanged.
- `orangebook` reads `ORAL-20`, `ORAL-21` and `ORAL-28` as `ORAL`, and
  `IM-IV` as systemic. A drug sold only as an oral contraceptive pack, or only
  as `IM-IV`, was "unclear" and now reaches the body. The `routes` list in a
  lookup shows `ORAL` where it showed `ORAL-28`.
- `audit.py` checks more tables and exits with status 1 for the new problems.
