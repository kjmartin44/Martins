# promoter-tf-drugs

A Python research tool that maps a gene's promoter to predicted transcription
factor binding sites and candidate drugs that may raise or lower gene expression.
It combines motif scans, annotations, optional experimental evidence, and heuristic
rankings. Results are hypotheses for further study, not validated treatment advice.

## The pipeline

One run, from a gene and a direction to ranked drug hypotheses. The [compact flowchart](docs/pipeline_detail.mmd) adds the branches and checks, and [how it works](docs/how-it-works.md) explains each step.

```mermaid
timeline
    title promoter-tf-drugs: one run, from a gene and a direction to ranked drug hypotheses
    section Before
        Setup : setup_data.bat fetches reference files : Input is a gene, raise or lower, a dosage scenario
    section Factor route
        1 Promoter : Ensembl or GENCODE : One window per start
        2 Scan : JASPAR at p 1e-5 : Chance measured on shuffled windows : Rescans wider if no signal
        3-4 Factors : GO role per factor : ReMap ChIP peaks : Neighbouring genes
        5-7 Drugs : Rank the factors : Approved drugs by direction : Grade each pair
    section Beside it
        8 Twelve routes : Own product : Upstream network : LINCS : Docking : Known regulators : Genetics and CTD as evidence
    section Merge
        9 Catalogue : Tier 1 to 5, three axes : Tier orders first
        10-11 Verdict and files : Ledger marks each source : FDA approved only shown : report, summary, results.json
    section After
        Checks : Retrospective against CTD : Where drugs are lost : Dose check
    section What each route returns
        Rows in the list : a Own product, tier 1 drug on the gene or RNA : b Upstream, tier 4 drug with its path : c LINCS, tier 2 to 4 compounds that moved the gene : d Factor activity, tier 3 compounds : e Screens, tier 2 ligands and tier 5 hypotheses : f Docking, tier 4 binders if the pocket validates : i DepMap and PRISM, tier 4 phenocopy drugs : l Known regulators, tier 4 drug at a regulator
        Evidence only : g Genetics, support for moving the gene : h CTD, known or novel label : k Dosage, how far and which way
        Own section : j Modalities, a designed ASO or CRISPRa sketch
```

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

To analyze a real gene with live services, first fetch the reference files. Git
does not carry them, and the full set is about 3 GB (`--small` skips the largest;
`--list` shows what is missing):

```sh
python scripts/setup_data.py
python -m promoter_tf_drugs --gene CCR5 --direction down --pfm-file data/jaspar2024.txt --out results/CCR5_down
```

On Windows, `scripts\run_everything.bat` does both and runs every stage.

To compare several candidate genes and narrow to the one with the least
collateral effect:

```sh
python -m promoter_tf_drugs.panel --genes GENE1,GENE2,GENE3 --control-gene ACTB --demo --out results/panel
```

Two things worth reading before acting on any output. Every report opens with
an **Evidence available in this run** table, because a run where the lookups
silently failed is otherwise structurally identical to one where they worked.
And the ranking is narrowed on four independent axes of collateral — how many
genes the factor controls, how body-wide it is, how many body systems the
drug moves, and whether a subtype-selective ligand can exist at all.

Saved example outputs, including a real run, are in [examples](examples/README.md).

See the [documentation](docs/index.md) for installation details, input formats,
methods, limitations, the Python API, and development instructions —
in particular [specificity and evidence](docs/specificity.md) and
[comparing several genes](docs/panel.md).

## License

This repository is licensed under [Creative Commons Attribution 4.0 International](LICENSE) (CC BY 4.0). The public data sources it reads have their own terms; see [Setting up the data](docs/data-setup.md).
