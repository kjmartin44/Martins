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

[How the package is laid out](architecture.md) is the map: the chain from a
gene symbol to a ranked table, which of the sixty-one modules owns each step,
and a table of every module with its purpose.

The short version. `pipeline.py` owns orchestration, report formatting and
export, and at about 4,600 lines it is the largest file by a wide margin —
splitting it is the first item under maintainability in the
[improvement review](improvements.md). The main CLI translates user options
into its `run` function, and both CLIs derive those options from `run`'s own
signature through `cli.run_kwargs_from_args`, so a new flag reaches both
without anyone maintaining a list. `http.py` centralises external requests.
`shortlist.py` and `explain.py` work entirely on saved result dictionaries.

The module table in `architecture.md` is generated, not hand-written:

```sh
python scripts/generate_module_table.py           # rewrite it
python scripts/generate_module_table.py --check   # fail if stale
```

The `--check` form runs as a test, so a new module cannot be added without
the map noticing.

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

### Why a green local run can be red in CI

The package declares Python 3.10 and is usually written on 3.12. Python 3.12
accepts f-string forms that 3.10 rejects: a newline inside the braces of a
single-quoted f-string, the outer quote character reused inside them, and a
comment inside them. A file using one imports without complaint on 3.12 and
fails to import on 3.10, so every test in that job errors before it runs.
This happened on 2026-10-01 and the 3.10 jobs failed from then on.

Two checks now catch it. The first step of each CI job compiles every file on
that job's Python and names the file and line. On 3.12 and later,
`test_no_f_string_needs_python_3_12_to_parse` finds the same constructs with the
tokenizer before you push. The fix is always to lift the expression out of the
f-string and name it first.

The test suite cannot import `tomllib`, which exists only from 3.11. Read
`pyproject.toml` with a regular expression instead.

### Reproducing CI on your own machine

Test a clean clone, not your working copy. A runner has no `data/` directory
beyond the one tracked file, an empty cache, and no files you downloaded, so a
test that quietly reads one of them passes on your machine and fails there.

```bash
git clone . /tmp/ci-check
cd /tmp/ci-check
PTD_CACHE="$(mktemp -d)" python -m pytest -q
```

To make any test that reaches for the network fail loudly instead of being
answered from your warm cache, point the proxy variables at a closed port:

```bash
HTTP_PROXY=http://127.0.0.1:9 HTTPS_PROXY=http://127.0.0.1:9 \
PTD_CACHE="$(mktemp -d)" python -m pytest -q
```

On Windows keep the clone and the temporary directory under a short path such
as `C:\Users\you\AppData\Local\Temp\ci`. Several tests write cache files
three or four directories deep, and a long parent path pushes them past the
260 character limit. That produces `WinError 3` failures that CI, which uses
short paths, never sees.

To check a different Python, compile the source with that interpreter. It
needs no dependencies and finds syntax the running interpreter accepts:

```bash
python3.10 -m compileall -q promoter_tf_drugs scripts tests
```

`python scripts/generate_cli_docs.py` pins the terminal width and writes LF line
endings, so the command reference is the same whatever your shell exports and
whichever platform you run it on.

## Auditing the curated tables

Several tables are hand-entered: the systemic drug classes, the receptor
families, the curated stabilisers, the P-gp substrates. Each says in its
docstring that it is incomplete, which is honest and useless on its own,
because nothing says *when* it has drifted. A table that is quietly wrong is
worse than one known to be incomplete — the report presents both with the
same confidence.

```bash
python -m promoter_tf_drugs.audit            # consistency, no network
python -m promoter_tf_drugs.audit --online   # also resolve names in ChEMBL
```

The offline pass runs as a test. It catches changes that would silently alter
behaviour: a drug in two systemic classes makes the burden it reports depend
on table declaration order, and a name stored with a counterion can never be
found, because lookup is on the salt-free form.

The online pass resolves every curated name against ChEMBL and cross-checks
our own `withdrawn` flags against theirs. The dangerous direction is
we-say-available / ChEMBL-says-withdrawn, which would put a withdrawn drug at
the top of a report unflagged. Names ChEMBL cannot answer for are counted as
**lookup errors, not as missing**.

Neither pass can tell you what is *absent* from a table. Only reading the
literature does that. Re-run `--online` every few months.
