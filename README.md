# promoter-tf-drugs

A Python research tool that maps a gene's promoter to predicted transcription
factor binding sites and candidate drugs that may raise or lower gene expression.
It combines motif scans, annotations, optional experimental evidence, and heuristic
rankings. Results are hypotheses for further study, not validated treatment advice.

## Quick start

Requires Python 3.10 or later. From this checkout:

```sh
python -m venv .venv
```

Activate the environment (`.venv\Scripts\Activate.ps1` in PowerShell, or
`source .venv/bin/activate` on Linux/macOS), then run:

```sh
python -m pip install -e .
python -m promoter_tf_drugs --demo --gene DEMOGENE --out results/demo
```

The demo runs offline on synthetic data. Open `results/demo/report.md`.
A completed run also writes `results.json`, `tf_summary.tsv`, `sites.bed`,
and `promoter.fa`.

To analyze a real gene with live services and the repository's motif file:

```sh
python -m promoter_tf_drugs --gene CCR5 --direction down --pfm-file data/jaspar2024.txt --out results/CCR5_down
```

See the [documentation](docs/index.md) for installation details, input formats,
methods, limitations, the Python API, and development instructions.
