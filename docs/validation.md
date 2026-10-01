# Validation record

[Documentation index](index.md)

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
