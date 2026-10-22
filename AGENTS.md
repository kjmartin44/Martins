# Working in this repository

This package ranks research hypotheses about changing gene expression. Keep
measured evidence, database annotations, synthetic fixtures, and predictions
distinct in code and documentation. A successful run does not validate a drug
or establish that it changes expression in a particular tissue.

## Before making changes

Read `pyproject.toml`, the relevant modules and tests, and `docs/index.md`.
Check `git status --short` and preserve work already present. Keep changes
focused on the request. Do not commit, push, publish, or change external data
unless the user asks.

## Package layout

- `generoutes/pipeline.py` orchestrates the analysis and writes reports.
- `cli.py` parses options; `shortlist.py` and `explain.py` read saved results.
- `http.py` owns caching, offline behavior, throttling, and request retries.
- The other package modules implement data adapters and analysis stages.
- `tests/` contains offline tests. `data/` holds repository reference files.
- `examples/` contains saved analyses. Put new run outputs under `results/`.

## Implementation rules

Support Python 3.10 and later. Keep NumPy as the only required runtime
dependency. Import optional chemistry libraries only where they are needed.
Use `Client` for external requests and preserve `--offline`. Demo mode must
never send a network request. Validate user input before expensive work.

Keep the public functions `run`, `markdown_report`, and `write_outputs`
compatible unless the task explicitly calls for an API change. Preserve the
meaning of saved result fields. Document any output or cache compatibility
change. Do not tune scientific scores just to make an example look stronger.

Promoter coordinates are 1-based and inclusive. Motif indices are 0-based and
relative to the sequence oriented along the gene. BED coordinates are 0-based
with an exclusive end, and its strand is genomic. Test both gene strands when
changing coordinate conversion.

## Validation

Create a virtual environment and install `python -m pip install -e
".[dev,chemistry]"`. The current full suite requires the chemistry extra.
Run `python -m pytest -q` and add focused regression tests for changed behavior.
Use temporary directories and mocked services in tests. Do not turn unit tests
into network checks. On Windows, use the virtual environment's Python if the
default `python` resolves to another installation.

For packaging changes, build the source distribution and wheel with
`python -m build`, then install the wheel in a separate environment and run the
demo from outside the checkout. Check `git diff --check` before handing off.
Report checks that could not be run and distinguish offline tests from live
service checks.

## Documentation

Keep the README limited to what the project is and how to quick start. Put
installation details, usage, inputs, outputs, methods, and maintenance guidance
under `docs/`. Keep command examples consistent with the actual parsers.

Write in plain, natural prose. Avoid compound-adjective stacks and unnecessary
hyphenation. Prefer verbs and short sentences to labels. Document limitations
as concrete behavior, and do not present heuristic scores as calibrated
probabilities or regulatory approval checks.
