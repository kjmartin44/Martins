# What this package needs

[Documentation index](index.md)

Everything the package imports or runs, in one place. The reference data
is a separate question with its own page: see
[getting the analysis tables](datasets.md).

NumPy is the only thing a run requires. Everything below it is optional,
and the report says which routes it lost when something is missing rather
than reporting that those routes found nothing.

## Python packages

| Package | Install with | What stops without it |
| --- | --- | --- |
| numpy | required | nothing runs |
| rdkit | `.[chemistry]` | structure features, clustering, counter-ion handling, the direction classifier, every virtual screen |
| scikit-learn | `.[chemistry]` | the agonist/antagonist classifier behind `--predict-direction` |
| pytest | `.[dev]` | the test suite |
| build | `.[dev]` | `python -m build` |
| meeko | `.[docking]` | the Vina path: it writes the ligand PDBQT, so every compound in a screen returns nothing |
| decoupler | `.[data]` | `setup_data.py` falls back to OmniPath for CollecTRI, which is a different table |
| openpyxl | `.[data]` | `build_transcript_sites.py` cannot read the TargetScan workbook |
| openbabel | `.[openbabel]` | nothing, while meeko is installed: it is the fallback receptor converter |

```sh
python -m pip install -e ".[chemistry]"            # the usual one
python -m pip install -e ".[dev,chemistry]"        # to run the tests
python -m pip install -e ".[docking,data]"         # docking and data prep
```

Open Babel has an extra of its own because it often has no wheel for the
running platform and builds from source. Inside `docking` a failure there
would take meeko down with it, and meeko is the one that matters. On
conda, `conda install -c conda-forge openbabel` is usually easier than
pip.

The full test suite needs `chemistry`. Roughly ten tests skip when an
optional file is absent from `data/`, which on a clean checkout is
always, so CI runs `pytest -rs` to name them.

## External programs

None of these is a Python package, so none can be declared in
`pyproject.toml` and none arrives with `pip install`.

| Program | Needed for | Notes |
| --- | --- | --- |
| AutoDock Vina | docking on Windows | `vina` on PATH. The package builds the box itself, because `--autobox_ligand` is a smina and gnina flag that Vina does not have |
| gnina | docking with a CNN rescore | Linux only; run through WSL from Windows |
| smina | docking | Linux only, same WSL route |
| `obabel` | receptor preparation | the Open Babel command line, tried after meeko |
| WSL | running the Linux binaries from Windows | `generoutes.wsldock` quotes every argument it passes |

Docking is the one part that is not self-contained, which is why it
splits in two: `python -m generoutes.docking --prepare` does every
networked step and leaves a folder, and `--bundle` docks that folder with
no network at all, on whichever machine has the binary. See
[docking](docking.md) for the validation gate a pocket has to pass before
any score is reported.

`generoutes-dock` is installed by the wheel and is the same thing as
`python -m generoutes.docking`.

## Checking what you have

```sh
python -m generoutes.docking --check
python scripts/setup_data.py --list
```

The first says which docking pieces are present and which are missing.
The second says the same for the reference tables. Neither installs
anything.
