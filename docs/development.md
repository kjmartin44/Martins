# Development

[Documentation index](index.md)

## Setup and checks

```sh
python -m pip install -e ".[dev,chemistry]"
python -m pytest -q
python -m build
git diff --check
```

Use the environment's Python explicitly when it is not activated. The test suite
uses synthetic data and mocked services. Live checks should be separate commands
with outputs and logs under `results/`; do not add network calls to unit tests.

After packaging changes, install the wheel into another environment and run
`python -m promoter_tf_drugs --demo --gene DEMOGENE --out demo` outside the
checkout. This checks installed imports and confirms that repository files are
not accidentally required. Also exercise each installed command's `--help`.

## Structure

`pipeline.py` owns orchestration, report formatting, and export. The main CLI
translates user options to its `run` function. Source adapters and scientific
helpers live in individual modules. `http.py` centralizes external requests.
`shortlist.py` and `explain.py` operate entirely on saved result dictionaries.

The existing `test_pipeline.py` covers many domains in one file. New operational
regressions are in `test_reliability.py`. Keep tests focused on observable behavior
and use fixed seeds for models. For request behavior, mock `urlopen` and sleep;
do not rely on service availability.

## Documentation maintenance

Keep the README short. Update the appropriate guide when changing behavior.
The command reference contains parser help; regenerate it after adding or
changing options with `python scripts/generate_cli_docs.py`. The generated file
includes the analysis, shortlist, and explain interfaces.

The CI workflow tests supported Python versions on Linux and Windows, installs
the chemistry extra, builds distributions, and smoke-tests a wheel outside the
checkout. A local passing test run does not mean CI has run or passed remotely.
