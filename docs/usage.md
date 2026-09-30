# Usage

[Documentation index](index.md)

Commands below assume an installed package and an activated environment.
They use one line so they work in PowerShell and POSIX shells.

## Start with the demo

```sh
python -m promoter_tf_drugs --demo --gene DEMOGENE --out results/demo
python -m promoter_tf_drugs --demo --gene DEMOGENE --direction down --out results/demo_down
```

Demo mode uses synthetic sequence, motifs, annotations, and pharmacology. The
gene argument labels the run; a real gene name does not turn the demo into a
real analysis. The report marks the mode and the assembly is `DEMO`. Demo mode
forces the HTTP client offline, including when an optional stage requests data.

## Analyze a gene

```sh
python -m promoter_tf_drugs --gene ABCA4 --direction down --pfm-file data/jaspar2024.txt --regulon-file data/trrust.tsv --out results/ABCA4_down -v
python -m promoter_tf_drugs --gene CCR5 --direction down --pfm-file data/jaspar2024.txt --out results/CCR5_down -v
python -m promoter_tf_drugs --gene SMN1 --direction up --pfm-file data/jaspar2024.txt --out results/SMN1_up -v
```

These specify a desired expression change. They do not assert that the genes
are dysregulated in a particular sample. `up` is the default. It seeks activation
of activating TFs and inhibition of repressors; `down` reverses that logic.

The default promoter window is 1,000 bases upstream and 200 downstream.
`--all-tss` selects additional TSS clusters; nearby starts can be grouped.
`--upstream 2000` changes the window. `--upstream-network` enables causal-network
search: these are separate options.

Local motif and regulon files do not make a real run fully offline. Gene
resolution, sequence, annotation, and drug lookups still need live or cached
responses. The package defaults to human data. Its `--species` option does not
make all drug, annotation, and curated evidence stages species independent.

## Add evidence selectively

| Option | Behavior |
| --- | --- |
| `--chip-bed peaks.bed.gz` | Count local ChIP peaks over the promoter. Match the genome assembly. |
| `--chip-api` | Request ReMap peaks; a local BED takes precedence. |
| `--regulon-file data/trrust.tsv` | Use target counts from a local TF-target table. |
| `--cell-type-check --expression-file expression.tsv` | Use local tissue expression for the overlap check. |
| `--tissue cns` | Apply the package's CNS delivery heuristics. |
| `--proliferation post-mitotic` | Apply the curated cell-cycle compatibility rules. |
| `--lincs-file signatures.tsv` | Read local drug/gene expression effects. |
| `--lincs` | Query expression signatures from a service. |
| `--genetics` | Add genetics context and lookup guidance. |
| `--upstream-network --upstream-depth 3` | Search signed upstream paths. |
| `--predict-direction` | Fit optional receptor direction models when data permits. |
| `--screen` | Run optional virtual screening; can require many requests. |
| `--design-screen` / `--design-modality` | Add experimental or sequence-design suggestions. |

`--signature-file` currently only records the number of genes in a reference
signature when the LINCS stage is enabled. It does not run a connectivity search.
See [methods and limitations](methods.md) before interpreting optional sections.

## Reuse cached data

```sh
python -m promoter_tf_drugs --gene SMN1 --pfm-file data/jaspar2024.txt --cache-dir results/cache --out results/SMN1
python -m promoter_tf_drugs --gene SMN1 --pfm-file data/jaspar2024.txt --cache-dir results/cache --offline --out results/SMN1_cached
```

The second command forbids network requests. It can fail if a required response
is absent. Optional stages may log a warning and return partial or empty results.
A successful offline run therefore does not prove that every source was available.
Use `--cache-dir` or set `PTD_CACHE` before starting Python to choose the cache.

## Read and narrow results

```sh
python -m promoter_tf_drugs.shortlist results/SMN1/results.json --preset strict --out results/SMN1/shortlist.md
python -m promoter_tf_drugs.shortlist results/SMN1/results.json --max-steps 0
python -m promoter_tf_drugs.explain results/SMN1/results.json "compound name"
```

The installed aliases are `promoter-tf-drugs-shortlist` and
`promoter-tf-drugs-explain`. Both read saved JSON without network calls.
Shortlist presets are `strict`, `mechanism`, `empirical`, `specific`, and `novel`.
Individual supplied filters override preset values. A zero step limit selects
direct routes. An empty shortlist is a valid result.

`--top 10` on the analysis command limits report tables in both stdout and
`report.md`. It does not truncate the JSON, TSV, BED, or FASTA exports. Some
optional stages have their own internal limits, so JSON is not a complete dump
of every source record. Output filenames are fixed and rerunning in the same
directory overwrites them; use a new directory to preserve a previous run.
