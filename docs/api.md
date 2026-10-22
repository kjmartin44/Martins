# Python API

[Documentation index](index.md)

The package exposes three functions:

```python
from generoutes import run, markdown_report, write_outputs

result = run("DEMOGENE", demo=True, cache_dir="results/cache", goal="down")
text = markdown_report(result, top_n=10)
paths = write_outputs(result, "results/python_demo", top_n=10)
```

`run(gene, **options)` returns a JSON-serializable dictionary. The Python keyword
for expression direction is `goal`; the CLI calls it `--direction`.
`markdown_report(result, top_n=25)` returns a string.
`write_outputs(result, outdir, *, top_n=25)` creates the directory and returns
the five `pathlib.Path` objects written there. Existing output files are replaced.

## Frequently used run options

| Keyword | Default | Meaning |
| --- | --- | --- |
| `upstream`, `downstream` | `1000`, `200` | Promoter window in bases. |
| `species` | `"homo_sapiens"` | Ensembl gene lookup species; later stages remain human-oriented. |
| `pvalue`, `min_ic` | `1e-4`, `6.0` | Motif threshold and minimum information content. |
| `pfm_file`, `chip_bed`, `regulon_file` | `None` | Local reference files. |
| `all_tss`, `chip_api` | `False` | Additional TSS clusters and ReMap lookup. |
| `mask_repeats` | `True` | Skip lowercase repeat sequence outside demo mode. |
| `goal` | `"up"` | Desired expression direction: `"up"` or `"down"`. |
| `grade_drugs`, `include_direct` | `True` | Drug grading and direct-route lookup. |
| `tf_weight` | `0.5` | TF confidence share of the combined grade. |
| `approved_only` | `True` | Use the package's source-specific approval filters. |
| `demo`, `offline` | `False` | Synthetic fixtures or cached service responses only. |
| `cache_dir` | `None` | Override the default cache directory. |

For the complete signature, use `help(run)` or inspect
[`pipeline.py`](../generoutes/pipeline.py). Boolean CLI switches are
translated explicitly in [`cli.py`](../generoutes/cli.py); for example,
`--no-grade` means `grade_drugs=False`.

## Errors and partial results

Invalid numeric ranges, empty gene names, unsupported direction, and unsupported
tissue states raise `ValueError` before the HTTP client is created. So does a
scenario that contradicts the direction, such as `scenario="haploinsufficiency"`
with `goal="down"`: the requirement and the ranking would point opposite ways.
The scenarios `null` and `modifier` accept either direction. Missing
explicit input files raise `FileNotFoundError`. Failed required requests can
raise `generoutes.http.HttpError`. Optional stages often log warnings and
return defaults instead of aborting. Retain logs with saved results when using
the API for batch analyses.

The package is synchronous. Local expression data currently resides in module
state, so independent analyses that use different expression tables should run
in separate processes until that state is isolated per run.
