# Validation record

[Documentation index](index.md)

## Current record, 20 October 2026

Run on Windows with Python 3.12, against a working tree that is not committed.

| Check | Result |
| --- | --- |
| Offline test suite | 1,821 passed |
| The same suite in a copy holding only what a commit would contain (no `data/`, empty cache) | 1,811 passed, 10 skipped; the skipped tests need data files |
| Wheel built, installed with NumPy only into a new environment, then the demo, the panel demo and `--help` for the four installed commands, all run outside the checkout | all ran |
| Python 3.10 grammar | every file in the package, scripts and tests parses |
| Linux, and Python 3.10 and 3.14 | not run here; the CI matrix covers them |
| `approved_library.tsv` rebuilt from ChEMBL by `setup_data.py` | identical to the existing file for 3,311 molecules; ChEMBL no longer lists the other 4 |
| Reference files fetched by `setup_data.py` into an empty folder (JASPAR, TRRUST, CollecTRI, ClinGen, HPA, HGNC, GO, LINCS, knockout libraries) | all fetched and passed their checks. ReMap, BindingDB, DepMap, PRISM and the GWAS Catalog were checked for their addresses and checksums, not downloaded |
| Live run, CYP24A1 raised, current code, full options without docking | completed; sample output in `examples/CYP24A1_up` |
| Live runs, ABCD2, LPL and FXN raised, with two-step regulators | completed; the expected drugs were found for ABCD2 and LPL, and for FXN all but omaveloxolone. See [checks](checks.md#rare-disease-genes-whose-drugs-are-known) |

The CYP24A1 run is a positive control. The report's verdict is **degraded**: the
motif scan found 13 sites against about 11 expected from shuffled windows. The vitamin
D analogues, the textbook answer, are found through the known regulators route and rank
10 to 15 of 65 approved drugs, behind nine tier-3 routes. That ordering is what the
retrospective check in [checks](checks.md) predicts: the routes reach about a tenth of
the approved drugs CTD records as raising a gene, and rank them no better than chance
once a drug's fame is held equal.

What this record does not show. No drug in any list has been shown to change the
gene. The 36 saved runs behind the figures in `checks.md` predate the motif scan
correction. The BUB1B and CYP3A4 runs made earlier in the same week also predate it, so
they are not part of this record. See [improvements](improvements.md#open-at-the-time-of-writing)
for the open items.

## Earlier record, 30 September 2026

This predates the motif scan correction, so its site counts describe the earlier scan.

Validation on Windows with Python 3.12.12, 2026-09-30. Live checks below used
the current public services, the checkout's 879 JASPAR matrices, its TRRUST file,
the default 1,000/200-base promoter window, and `pvalue=1e-4`. These were real
analyses with `demo=False`. Optional ChIP, chemistry screening, LINCS, and
genetics stages were not enabled.

| Gene | Requested direction | Transcript | Motif sites | TFs | Drug grades | Catalogue routes | Scan verdict |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| ABCA4 | down | ENST00000370225 | 246 | 199 | 68 | 68 | weak |
| CCR5 | down | ENST00000292303 | 241 | 200 | 96 | 97 | weak |
| SMN1 | up | ENST00000380707 | 93 | 87 | 14 | 23 | no signal |

All three commands exited successfully and produced the five documented output
files. All resolved to GRCh38 and selected an Ensembl canonical transcript.
ABCA4 uses the negative gene strand; CCR5 and SMN1 use the positive strand.
The pipeline's null expectations were 181.1, 187.1, and 100.6 sites respectively.
The route counts do not indicate validated drug effects.

**Incomplete evidence:** CTD lookups returned HTTP 302 for every gene. The
pipeline logged the failures and continued, so these runs do not establish the
absence of prior drug/gene links. There was no ChIP evidence supplied. The SMN1
scan produced no enrichment above the package's rough chance expectation.

The local outputs and logs are under `results/validation/`, which Git ignores.
Each run used the following command with the appropriate gene and direction:

```sh
python -m promoter_tf_drugs --gene ABCA4 --direction down --pfm-file data/jaspar2024.txt --regulon-file data/trrust.tsv --cache-dir results/validation/cache --out results/validation/ABCA4_down -v
```

The full offline suite passed: **323 tests**, including 37 new regression cases.
Each live result was also checked for output completeness, promoter sequence
length, every BED interval and genomic strand, and TF/drug direction consistency
(68 ABCA4 pairs, 96 CCR5 pairs, and 14 SMN1 pairs).

Software regressions cover numeric validation, missing input files, report
limits, positive and negative BED coordinates and strands, zero-valued shortlist
options, local matrix filtering, demo network isolation, cache identity and
concurrent writes, malformed cache entries, failed cache persistence, query
parameters, retry delays, and redirect cycles. Markdown links and fenced blocks,
Python 3.10 syntax compatibility, and `git diff --check` also passed. Execution
on other Python versions and operating systems is left to the added CI matrix;
that remote workflow has not been run as part of this local task.
