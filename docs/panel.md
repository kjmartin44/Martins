# Comparing several candidate genes

[Documentation index](index.md)

```bash
python -m promoter_tf_drugs.panel --genes CCR5,PCSK9,LDLR --direction down \
    --pfm-file data\jaspar2024.txt --chip-bed data\remap2022.bed.gz \
    --cell-type-check --expression-file data\rna_tissue_consensus.tsv \
    --control-gene ACTB --max-regulon 600 --max-burden moderate \
    --out results\panel
```

The rest of the package asks "which drugs move THIS gene". If you have
several candidate targets and must pick one, that is the wrong unit: you need
the gene+drug pair with the least collateral effect on everything else, and
that comparison can only be made across genes.

**The genes are treated as alternatives, not a panel to move together.** That
inverts how a recurring hit reads: a drug appearing against four of five
candidates is not four times validated, it is a drug that moves a lot of
things and one of them happens to be whichever gene you pick. `1 gene` is the
good case.

## A funnel, not a score

Collapsing specificity into one number hides which criterion did the work.
Each stage applies one criterion and reports what it removed.

| # | stage | axis |
|---|---|---|
| 0 | admissible genes | provenance — removes *genes*, not rows |
| 1 | evidence floor | tier, grade, direction conflicts |
| 2 | not genome-wide | chromatin agents: specificity ~0.05, no grade |
| 3 | narrow regulon | [axis 1](specificity.md) |
| 4 | narrow TF tissue footprint | [axis 2](specificity.md) |
| 5 | selective ligand achievable | [axis 4](specificity.md) |
| 6 | tolerable clinical burden | [axis 3](specificity.md) |
| 7 | annotate recurrence | diagnostic — removes nothing |

Three properties keep the accounting honest.

**Missing data passes, and is counted.** Dropping a candidate for lacking a
field its route never had is backwards. But a stage reporting "removed 0"
when it had nothing to work with reads like a clean bill of health, so
`not assessable` is printed beside it — and a stage nobody set a threshold
for says *that* rather than pretending to have run.

**Correlated stages are exposed.** Global regulators tend to have large
regulons *and* broad expression *and* systemic drug classes — three views of
one property, which displayed as three independent tests overstates the
evidence. Each stage reports **unique kills**: rows it removed that no later
stage would have removed anyway. Large removed with small unique means the
stage is nearly redundant with the next.

**Genes, not rows.** Counts are of distinct (active moiety, target) pairs per
gene. With `--all-tss` a gene contributes several near-identical windows —
CCR5's two overlapped by 71% — and counting rows would let the
most-annotated gene win on annotation volume. `--per-gene-top` caps the table
for the same reason.

The deliverable is the **gene × stage grid**: the stage where a candidate
runs out of routes is the answer to which gene is cleanest.

## The negative control

Without one, **an empty shortlist and a misconfigured filter look
identical.** `--control-gene` runs an extra gene — a housekeeping gene is the
natural choice — through exactly the same stages, and the report compares its
survival against the candidate median. If the control survives as well as
your candidates, the funnel is passing routes on properties every gene
shares, and the surviving list is not evidence about your genes specifically.
The report says so in those words.

The control is excluded from the ranking and from recurrence, so it cannot
contaminate the comparison it exists to validate.

## Re-narrowing is free; re-running is not

A fresh run costs tens of minutes per gene; tuning a cutoff should cost a
second.

```bash
python -m promoter_tf_drugs.panel --renarrow results\panel\panel.json \
    --max-regulon 300 --out results\panel-tight
```

No network at all. The control gene is preserved.

## Known biases, corrected or declared

**The winner is a maximum.** Selecting the best of N genes × M drugs inflates
the top row, because you are looking at it *because* its predicted score came
out high. The report states the number of comparisons and, when the margin
over the runner-up is within noise, says the two are not separable.

**Genes are not run identically unless you make them so.** A gene run without
`--chip-bed` loses every chip-dependent stage for reasons unrelated to its
biology. `parity()` reports which options differ; mixed `--direction` is a
hard error rather than a warning.

**Tractability is not therapeutic value.** The funnel narrows on how cleanly
a gene can be drugged. Human genetic support is printed *above* it and gates
nothing, so the easiest gene to drug cannot quietly win the argument.

## Presets

| preset | question |
|---|---|
| `narrow` | what survives the three specificity axes? |
| `narrow+clean` | what would I defend in a meeting? |
| `exploratory` | what is there at all? |

Individual flags override a preset. See the [command
reference](cli.md#python--m-promoter_tf_drugspanel) for every option.
