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

## One list from docking and the rest

The run ranks routes by tier and then score. Docking hits against the gene's
own protein enter that ranking as tier 4, below every curated factor route
however well the pocket validated (for CCR5 the pocket recovered its own known
ligands at AUC 0.93). The two kinds of evidence answer different questions, so
`combined` keeps them as separate rankings and fuses them:

```bash
python -m promoter_tf_drugs.combined --results results/CCR5_down/results.json \
    --panel results/CCR5_down_panel/panel.json \
    --docking results/CCR5_down/docking/product.json --top 40
```

`--docking` can be repeated, and the launcher passes the product dock and each
factor dock. The routes list is the panel's surviving routes when a panel file
is given, since the funnel has already applied its off-target criteria, and
otherwise every route in the run. Routes made by docking are left out of it so
docking is not counted twice. A docking result whose validation gate did not
pass is not used and is listed with its reason.

A drug's score is the sum of `1 / (60 + rank)` over the lists it is in
(reciprocal rank fusion). It uses only ranks, so the different scales never
meet. The constant 60 is the conventional one and is not tuned here, no list is
weighted, and tied drugs share the best rank of their group. A drug in two lists
beats one in a single list at the same rank.

Close analogues are folded into the best-ranked one so a family does not fill
the top. Two drugs are folded when their structures, without counter-ions, have
a Tanimoto similarity of at least 0.75 on Morgan fingerprints (radius 2, 2,048
bits), or when one is the same drug under an ester name. The leader is the
better-ranked member, a drug joins only the first leader it is close to, and a
chain of analogues does not merge into one row. The folded drugs are named
beside the leader with their own ranks, and the full list stays in
`combined.json`. `--similarity` changes the cut and `--no-fold` lists every
drug. Similar structure is not equivalent behaviour: analogues can differ in
potency, selectivity, tissue reach and safety, so a folded drug is a neighbour
to check, not a duplicate to drop. A compound with no structure in the approved
drug table stays on its own.

Each row has one sentence on how that drug would act on the gene, written from
the route's own fields (what it acts on, whether its binding is measured, the
factor's role at the promoter). It states the mechanism the route assumes. For
most routes the step from the drug's target to the gene is not measured, and
the sentence says so. The paragraph above the table gives the dosage
requirement when the run declared a scenario. Nothing in it says that moving
the gene treats the condition.

The score orders drugs and is not a probability. A drug high in both lists is
high in two rankings, each as weak as it was alone, and lists that share an
input are not independent. Every docking row says binding is not direction
unless the run recovered the conformation, and a drug found only by docking
carries its systemic class if it has one. The command writes `combined.md`,
`combined.json` and `combined.tsv`. The launcher runs it last and writes
`results/GENE_DIR_combined`.
