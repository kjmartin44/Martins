# Inputs and outputs

[Documentation index](index.md)

## Reference inputs

| Input | Format and behavior |
| --- | --- |
| `--pfm-file` | JASPAR text matrices: a `>matrix_id name` header followed by A, C, G, and T count rows. Without a local file, motifs come from JASPAR. |
| `--matrix-id` | Repeat to select matrix IDs. |
| `--chip-bed` | Plain or gzip BED with target information in the name column. See `evidence.py` for supported names. Coordinates must match the promoter assembly. |
| `--regulon-file` | TRRUST TSV, a two-column TF/target TSV, or a DoRothEA CSV with TF/source and target/gene headers. Used for target counts. |
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

## Output files

Every completed analysis writes five files:

| File | Contents |
| --- | --- |
| `report.md` | Human-readable rankings, explanations, and caveats. |
| `results.json` | Query metadata, scoring weights, direct therapies, and per-promoter analysis blocks. |
| `tf_summary.tsv` | TF rows with motif, role, ChIP, score, and drug fields. |
| `sites.bed` | Predicted motif sites in genomic coordinates with a UCSC track header. |
| `promoter.fa` | Promoter sequences oriented 5′ to 3′ along the gene. |

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
