# Tables that ship in this folder

Everything else under `data/` is downloaded or built by
`python scripts/setup_data.py` and is ignored by git. See
[the data setup guide](../docs/data-setup.md). The files below are tracked on
purpose, so a fresh clone can screen a drug library and run the dose check
before any download. Each one carries the terms of its source, and the rest of
the repository is CC BY 4.0.

| File | What it holds | Source and terms | Rebuild |
| --- | --- | --- | --- |
| `approved_library.tsv` | 3,315 approved small molecules that have a structure: name, SMILES, ChEMBL id, highest clinical phase, a withdrawn flag and an oral flag. File dated 2026-10-05. | The ChEMBL web API, <https://www.ebi.ac.uk/chembl/>. ChEMBL is licensed CC BY-SA 3.0, so this table is too: credit ChEMBL, and distribute it or anything built from it under CC BY-SA 3.0 or a compatible licence. Cite Zdrazil et al., "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods". The ChEMBL release is not recorded in the file. | `python scripts/setup_data.py --only library` |
| `fda_label_cmax.tsv` | Human peak blood levels read from the pharmacokinetics section of FDA drug labels, for 42 drugs. Every row names its label, the date and the openFDA set id, and quotes the sentence it used. | The openFDA drug label endpoint, <https://open.fda.gov/terms/>. openFDA content is public domain under a CC0 1.0 dedication, with a warning that some third party works on it may not be. Attribution is not required and is welcome. | `python scripts/build_fda_cmax.py` |
| `fda_label_cmax_review.tsv` | The 150 drugs the label reader set aside, with the reason and an example sentence. | The same. | `python scripts/build_fda_cmax.py` |
| `hpa.tsv` | A small tissue to organ lookup. It is not an expression table. | Not recorded. | none |

The two label tables grow as `build_fda_cmax.py` looks up more drugs, so a run
can leave them modified. A label can change after a row is read, and the row
keeps the label date so a reader can tell.

A table built from your own runs, `ctd_fame.tsv`, is not here. It holds counts
taken from CTD, which allows research and educational use only, so it stays
out of a repository that anyone may reuse.
