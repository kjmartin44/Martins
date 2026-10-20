# Example outputs

Two saved runs, so the output can be read without installing anything or
downloading the reference files. `results.json` is left out (about 1.3 MB, and
git ignores it); the other files of each run are here.

| Folder | What it is |
|---|---|
| `demo/` | The offline demo on synthetic data: `python -m promoter_tf_drugs --demo --gene DEMOGENE --out results/demo`. Nothing in it is a real finding. |
| `CYP24A1_up/` | A real run, raising CYP24A1, with live services and the full reference set. Read its README first. |

The demo was made on 20 October 2026 and the CYP24A1 run on 21 October 2026, from an
uncommitted working tree of version 0.1.0. A real run uses live services, so rerunning it later can give a different
list.
