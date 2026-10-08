# Installation

[Documentation index](index.md)

Use Python 3.10 or later. The core package requires NumPy. Run these commands
from the checkout directory containing `pyproject.toml`.

## Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m promoter_tf_drugs --demo --gene DEMOGENE --out results/demo
```

Calling the environment's executable directly avoids activation policy issues
and makes it clear which Python is being used. `python -m pip` and a bare `pip`
can refer to different installations; use the former.

## Linux and macOS

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -m promoter_tf_drugs --demo --gene DEMOGENE --out results/demo
```

Once the environment is activated, `promoter-tf-drugs` is equivalent to
`python -m promoter_tf_drugs`.

## Optional dependencies

```sh
python -m pip install -e ".[chemistry]"
python -m pip install -e ".[dev,chemistry]"
```

`chemistry` installs RDKit and scikit-learn for structure features, direction
classification, and virtual screening. `dev` adds pytest and the build frontend.
The full test suite currently exercises chemistry features and needs both.
Docking helpers additionally require separately installed tools such as Vina
and receptor preparation software; there is no docking command in the main CLI.

## Windows launchers

The batch files under `scripts/` change to the repository root before running.
They use `python` from `PATH`, so run them from a terminal with the intended
environment activated. `run_demo.bat` uses synthetic fixtures.
`run_analysis.bat` and `run_everything.bat` use live services. The latter may
install missing chemistry dependencies and download reference files.

`run_everything.bat` runs every stage at full effort, and a run takes a long
time. It docks the whole approved library (about 3,300 compounds, not the first
600) at exhaustiveness 16, docks five transcription factors unless you say
otherwise, and runs the factor activity test with 2000 permutations when you
name a tissue. A pocket that cannot be validated is refused before the library
is touched, so a gene whose structure or ligands are missing finishes the
docking stages in minutes. A pocket that validates costs about two hours per
target. The stages that need a large optional file (PRISM and DepMap for the
sensitivity scan) run only when the file is in `data/`.

The launchers are convenience wrappers. Prefer the documented Python commands
for reproducible runs and when diagnosing failures.

## Troubleshooting

- **No module named promoter_tf_drugs:** install the checkout with the same
  Python used to run the command. An installed package can run outside the repo.
- **No module named pip:** select a standard Python installation with pip, then
  create the virtual environment from that interpreter.
- **A chemistry import fails:** install the `chemistry` extra in that environment.
- **A live service fails:** inspect stderr or repeat with `-v` for a traceback.
  Cached responses can reduce repeated work, but they are not stage checkpoints.
- **Expression input is rejected:** see [accepted formats](data.md). The included
  `data/hpa.tsv` is not a gene-expression table.

See the [Python Packaging User Guide](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/)
for the metadata and optional dependency conventions used in `pyproject.toml`.
