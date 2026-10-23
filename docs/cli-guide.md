# Running it from the command line

[Documentation index](index.md)

This page says what you have to type, what the program fills in for you, and what
each switch costs. The full option list is in the [command reference](cli.md), which
is generated from the parser.

## The least you need

```sh
python -m generoutes --gene CCR5
```

Only `--gene` is required. Everything else has a default. That one command asks for a
list of ways to **raise** CCR5 (the default direction), fetches the motif file from the
web, and writes `results/report.md` and its companions. It works, but it leaves out most
of the evidence, so a real run adds a few options. Three levels:

| Level | Command | What it gives you |
|---|---|---|
| Try it | `python -m generoutes --demo --gene DEMOGENE --out results/demo` | Synthetic data, no network, nothing real |
| A basic real run | `python -m generoutes --gene CCR5 --direction down --pfm-file data/jaspar2024.txt --out results/CCR5_down -v` | The promoter scan, factors and their drugs |
| Everything | `scripts\run_everything.bat` (Windows), or the long command below | Every route, the funnel and the combined list |

The long command is what the launcher runs for the analysis, without docking. Fill in the
gene, the direction and a window, and every file named must exist (`python scripts/setup_data.py`
fetches them):

```sh
python -m generoutes --gene GENE --direction up \
  --pfm-file data/jaspar2024.txt --regulon-file data/chip_regulon.tsv \
  --chip-bed data/remap2022.bed.gz --gtf data/gencode.basic.gtf.gz \
  --expression-file data/rna_tissue_consensus.tsv --cell-type-check \
  --upstream 1000 --downstream 200 --pvalue 1e-5 --all-tss --tf-weight 0.7 \
  --predict-direction --upstream-network --upstream-depth 3 \
  --regulator-route --regulator-steps 2 \
  --lincs --lincs-gmt data/lincs_chem_pert_up.gmt data/lincs_chem_pert_down.gmt \
  --lincs-geneinfo data/lincs_geneinfo.txt --genetics --screen --design-screen \
  --design-modality --promoter-dna --screen-product --find-analogs \
  --sensitivity-scan --out results/GENE_up -v
```

## What you must fill in

| Setting | When | Why |
|---|---|---|
| `--gene` | always | The symbol to analyse. The only required option. |
| `--activity-tissue` | if you use `--tf-activity` | The run refuses without it: the test infers factor activity from one tissue's expression and has nothing to infer from otherwise. Spell the tissue as the expression table does ("liver", "lymph node"). |
| `--direction` | if you use `--scenario` | The two must agree (a haploinsufficiency plan needs `up`). The run refuses a contradiction before it starts. |

## Files you point it at

None are required, but a feature that needs a file does nothing without one. Files in `data/`
under their usual names are found without a flag (marked auto).

| Feature | Option | What happens without it |
|---|---|---|
| The motif set | `--pfm-file data/jaspar2024.txt` | Motifs are fetched from JASPAR, which is slower and needs the network. Recommended. |
| Measured binding | `--chip-bed data/remap2022.bed.gz` | The factor ranking rests on sequence alone. `--chip-api` asks the web instead. |
| How many genes a factor controls | `--regulon-file` | The bundled tiers are used. |
| Tissue overlap | `--expression-file` with `--cell-type-check` | The Open Targets web service is used instead. |
| Offline gene lookup | `--gtf`, `--genome` | Ensembl is asked. |
| Measured drug effects | `--lincs` with `--lincs-gmt UP DOWN` and `--lincs-geneinfo` | The route is skipped, or a web service is asked. |
| Does the gene move when a factor is lost | `--knockout-gmt` | Found in `data/` (auto). Without them the report says responsiveness is not measured. |
| Known regulators | CollecTRI at `data/collectri.tsv` (auto) | `--regulator-route` records that no signed table was found. |
| Cell dependency, drug sensitivity | `--depmap-file`, `--prism-file` | Found in `data/` (auto). `--sensitivity-scan` needs both. |
| Docking results | `--docking-results FILE` | No docking enters the run. |
| How much of the order is drug fame | `data/ctd_fame.tsv` (auto), built by `python scripts/build_fame_table.py` from your own saved runs | No candidate is checked for being a compound the literature records for every gene. That is the control the retrospective found mattered most. |
| Which miRNA or RBP acts on the transcript | `data/transcript_sites.tsv` (auto), built by `python scripts/setup_data.py --only transcript_sites` from TargetScan, miRTarBase and ENCODE eCLIP | The UTR scans report motifs with no factor attached, so there is nothing to knock down or antagonise. |
| Drugs measured on the PROTEIN | `data/drug_protein_abundance.tsv` (auto) | Every abundance claim in the run is about mRNA, so a stabiliser or a degradation blocker is invisible. decryptE and ProTargetMiner fit; LINCS P100 and CPTAC do not. |
| What is wrong with the gene | `--scenario`, with `--broken-allele-makes`, `--target-fraction`, `--condition` | The ranking states no requirement. `haploinsufficiency` assumes the surviving copy is wild type; for a recessive condition where it is not, use `biallelic_hypomorph`. |

## Defaults you can change

<!-- defaults:start -->
| Option | Default | Meaning |
|---|---|---|
| `--direction` | `up` | raise (`up`) or lower (`down`) the gene |
| `--upstream` | `1000` | bases before the start site that are scanned |
| `--downstream` | `200` | bases after it |
| `--pvalue` | `0.0001` | per-position threshold for calling a motif site (the launcher uses 1e-5) |
| `--min-ic` | `6.0` | drops motifs below this many bits of information |
| `--tf-weight` | `0.5` | share of a factor's grade that rests on the promoter call (the launcher uses 0.7 with ChIP data) |
| `--proliferation` | `unknown` | cycling, post-mitotic, mixed or unknown |
| `--null-shuffles` | `30` | shuffled windows used to estimate chance sites; 0 uses a formula |
| `--upstream-depth` | `3` | longest signalling path searched with `--upstream-network` |
| `--upstream-factors` | `8` | factors the upstream search covers |
| `--upstream-drugs-per-node` | `5` | drugs listed per upstream node |
| `--regulator-steps` | `1` | with `--regulator-route`: 1 for a gene's own regulators, 2 to add theirs (the launcher uses 2) |
| `--top` | `25` | rows in the printed tables |
| `--species` | `homo_sapiens` | human data throughout |
| `--out` | `results` | folder for the report and its files |
<!-- defaults:end -->

## On by default (turn off with)

| Switched on | Turn it off with |
|---|---|
| Skip repeat-masked sequence | `--no-mask-repeats` |
| Grade each factor and drug pair | `--no-grade` |
| Drugs on the gene's own product | `--no-direct` |
| Estimate each route's chain of evidence | `--no-outcome` |
| Ask CTD and list its records | `--no-improvements` |
| Compare orthosteric and allosteric options | `--no-mechanisms` |
| Stabiliser and degradation routes | `--no-stabilizers` |
| Drug-gene interactions from DGIdb | `--no-dgidb` |
| Drug links from Open Targets | `--no-opentargets` |
| Approved drugs only | `--include-investigational` |
| Reports list only FDA approved drugs | `--include-non-fda` |
| Retry a no-signal scan at -5000/+1000 | `--no-widen` |
| CTD chemicals as a route | `--no-ctd-route` |

## Off by default (turn on with)

| Switch | What it adds, and what it costs or needs |
|---|---|
| `--all-tss` | One window per distinct start site. More work per gene. |
| `--cell-type-check` | Whether the gene and each factor are expressed in your tissue. Needs an expression table, or asks the web. |
| `--predict-direction` | Fits a model of whether a drug turns a receptor on or off. Needs `pip install -e ".[chemistry]"`. |
| `--upstream-network` | Drugs at proteins upstream of a factor. About 80 requests the first time. |
| `--regulator-route` | Drugs at the factors recorded as regulating the gene. Needs CollecTRI. |
| `--regulator-measured` | With the route: drugs with only a measured activity at a regulator. About 250 lookups per gene on a first run, and no direction. |
| `--lincs` | Compounds measured to move the gene. |
| `--genetics` | Human genetic support for moving the gene. |
| `--tf-activity` | Infers which factors are active. Needs `--activity-tissue`. |
| `--screen` | Virtual screens against the receptors. Many requests. |
| `--design-screen` | The experiment this analysis prioritises. A promoter reporter by default; a screen of the approved drug library against the gene's PROTEIN in patient cells when `--scenario` says the protein is the problem, or with `--has-ptc` or `--broken-allele-makes truncated`. |
| `--design-modality` | siRNA, ASO and CRISPRa candidates designed from the transcript. |
| `--promoter-dna` | DNA-directed agents that prefer sequence features of this promoter. |
| `--screen-product` | Screens the gene's own protein for binders. |
| `--find-analogs` | Close analogues of the candidates. |
| `--sensitivity-scan` | Drugs that behave as if they removed a factor. Needs the DepMap and PRISM tables. |
| `--chip-api` | Ask the web for ChIP peaks when you have no BED file. |
| `--include-unknown-direction` | Keep drug-gene links that record no direction. |
| `--include-local-only` | Keep drugs approved only for local use (topical, eye, inhaled). |
| `--has-ptc` | Say the gene carries a premature stop codon in the alleles you care about. |
| `--avoid-genotoxic` | Say the patient already has chromosomal instability or is cancer-predisposed. Compounds that work by damaging DNA or the spindle drop to the bottom tier and say why. They are not deleted. |
| `--offline` | Never use the network; cached answers only. |
| `--demo` | Synthetic data, never the network. |
| `--quiet`, `-v` | Fewer or more progress and log messages. |

## What the program refuses

It stops before the run starts, rather than at the end, when you ask for:
`--tf-activity` without `--activity-tissue`; `--activity-permutations` without
`--tf-activity`; a `--scenario` whose direction contradicts `--direction`; or a value out
of range (for example `--null-shuffles -1` or `--upstream-factors 0`).

## How the launcher differs

`scripts\run_everything.bat` asks for the gene, the direction, the tissue state and the
dosage scenario, then runs the long command above with a window of -1000/+200, a p-value
of 1e-5, `--all-tss`, `--tf-weight 0.7` when ChIP data is present, and `--regulator-steps 2`.

Three of its questions change what the run produces rather than how it is tuned. For any
scenario where a copy is impaired it asks what the broken copy makes, and whether the
transcript carries a premature stop; either answer turns `--design-screen` from a promoter
reporter into a screen of the approved library against the gene's protein in patient cells.
It then asks whether DNA-damaging compounds should be demoted, which passes
`--avoid-genotoxic`. That one defaults to yes, because much of what the literature records
as raising a gene is a carcinogen, a tumour promoter or an irritant, and for a patient with
a cancer predisposition those are not candidates at any score.
It adds docking before the analysis, a funnel against a control gene after it, and a dose
check. It uses a local file for each option only if that file exists in `data/`. It also
asks how many factors to search for upstream drugs. Blank keeps the default of 8, and a
number is passed as `--upstream-factors` to the analysis and to the funnel. A larger number
adds drugs that reach the gene only through the network, all of them tier 4, and takes a few
more minutes.

## What runs without being asked

Some checks have no switch. They cost nothing, they only ever lower a
route or add a caution, and the report says when a missing file stopped
one from running.

| Check | What it does |
|---|---|
| Coverage audit | Marks each of the six steps between a gene and a protein as searched, listed only, or not covered. A full table of transcriptional routes is not evidence that the other five were considered. |
| Scan quality | Above a 25% false discovery proportion, routes that rest on the motif scan are capped at tier 3. |
| Proliferation | For a gene whose transcript tracks the cell cycle, a compound that moves the rest of the mitotic programme with it is halved. Needs the LINCS `.gmt` files. |
| Drug fame | A compound the literature records for most genes loses 15% of its score. Needs `data/ctd_fame.tsv`. |
| UTR scan | Upstream ORFs in the 5'UTR and AU-rich elements in the 3'UTR, each against a shuffle of the same sequence. |
| Degradation | Names the ligase that destroys the protein and the enzyme that reverses it, and asks whether any approved drug acts on either. |

## Reading the result

Open `report.md`, and read its **Evidence available in this run** table first. A source
that failed is recorded as failed and never as zero. Then read **How much of the problem
this run examined**, which says which of the six steps were searched and which were only
named — a full ranking is not evidence that the rest were considered. `summary.md` is the
one-page version, `results.json` keeps everything, and `tf_summary.tsv` is the factor table.

The order within a tier carries less than it looks. Tested against CTD, the ranking orders
drugs no better than chance once drug fame is held equal, so read the list as a stratified
pool rather than a priority order, and use the tiers, the method count and the exclusions
rather than the score.
