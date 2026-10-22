# A full workup, step by step

[Documentation index](index.md)

Five steps from a gene symbol to a shortlist you can test. Each one says
what it produces and what would make you stop.

```
1. Set up            ──▶  a machine that runs      --check lists the gaps
        │
2. Analyse the gene  ──▶  report.md                every route, with its tier
        │
3. Dock the pocket   ──▶  dock.json                refused unless it validates
   (optional)
        │
4. Narrow it         ──▶  panel.md                 the funnel, stage by stage
        │
5. Read it, then test──▶  a shortlist              and the experiment to run
```

---

## 1. Set up

```sh
python -m pip install -e ".[dev,chemistry]"
python -m generoutes --demo --gene DEMOGENE --out results/demo
```

The demo runs offline on synthetic data and never touches the network.
If it writes a report, the installation is sound.

For a real gene you want local data files, because they make the run
reproducible and let it work with no network at all. The ones that change
the answer most, in order:

| File | What it buys |
| --- | --- |
| `jaspar2024.txt` | 879 motifs instead of the small bundled set |
| `remap2022.bed.gz` | measured ChIP occupancy — the largest single gain |
| a regulon table | how many genes each factor controls |
| `rna_tissue_consensus.tsv` | whether a factor is expressed where your gene is |
| an approved-drug table | the screening library, and known ligands per target |

[Inputs and outputs](data.md) says where each comes from. Docking needs a
binary as well:

```sh
python -m generoutes.docking --check
python -m generoutes.wsldock
```

**Stop if:** the demo does not produce a report. Everything below assumes
it does.

## 2. Analyse the gene

One command does the whole chain: promoter, motif scan, which factors are
plausible, which have measured occupancy, and which approved drugs act on
them.

```sh
python -m generoutes --gene CCR5 --direction down \
    --pfm-file data/jaspar2024.txt \
    --chip-bed data/remap2022.bed.gz \
    --regulon-file data/chip_regulon.tsv \
    --expression-file data/rna_tissue_consensus.tsv \
    --cell-type-check --genetics --upstream-network \
    --tf-activity --activity-tissue "lymph node" \
    --out results/ccr5
```

`--direction` is the one argument you cannot get wrong by accident: it
decides whether an activator is wanted on or off, and reversing it
inverts every recommendation.

Writes `report.md`, `summary.md`, `results.json`, `tf_summary.tsv`,
`sites.bed` and `promoter.fa`.

While it works it prints a checkpoint per stage to stderr, with the
elapsed time and what the stage found:

```
[  0:00] CCR5: lowering
[  0:04]  2/10  motifs loaded - 879 at p <= 0.0001; scanning
[  0:31]  4/10  measured occupancy - peaks for 103 of 200 factor(s)
[  3:18]  6/10  drug/factor pairs graded - 412 pair(s)
```

They go to stderr, so piping the report to a file leaves them on the
terminal and the markdown untouched. `--quiet` turns them off; `-v`
leaves them on and adds the full lookup log underneath, which is a
debugging tool rather than progress.

**Stop if:** the evidence table at the top of the report says a source you
asked for failed. A failed lookup is recorded as failed, never as zero,
and the verdict degrades to say so — but it is still on you to decide
whether to proceed without it.

## 3. Dock the pocket

Optional, and the only step that looks at a structure rather than an
annotation. Worth running when a factor near the top has a structure with
a ligand in it, or when you want to screen the gene's own protein.

```sh
python -m generoutes.docking --gene VDR --state active \
    --jobs 8 --cpu 4 --library data/approved_library.tsv \
    --out results/ccr5/dock.json
```

Then fold it back into the analysis:

```sh
python -m generoutes --gene CCR5 --direction down \
    --docking-results results/ccr5/dock.json --out results/ccr5
```

**Stop if:** the pocket does not validate. Docking is tested first on
ligands whose answer is known, and when it cannot separate those from
decoys the run reports the refusal instead of a ranking, because the
scores would carry no information. [The docking guide](docking.md) covers
what the score is and is not worth.

## 4. Narrow it

Two ways, depending on whether you are choosing between genes.

**One gene.** Re-narrow the saved result — no re-run, no network:

```sh
python -m generoutes.shortlist results/ccr5/results.json \
    --preset specific
```

**Several candidate genes.** The panel runs each one and puts them
through a staged funnel, with a housekeeping gene as a control:

```sh
python -m generoutes.panel --genes GENE1,GENE2,GENE3 \
    --control-gene ACTB --direction down \
    --max-regulon 600 --max-burden moderate \
    --out results/panel
```

The genes are treated as **alternatives you are choosing between**, which
inverts how a recurring drug reads: a compound hitting several of them is
unselective, not corroborated.

Tuning a cutoff never re-runs anything:

```sh
python -m generoutes.panel --renarrow results/panel/panel.json \
    --max-burden low --out results/panel_strict
```

**Stop if:** the control retains as many routes as your candidates. That
means the funnel is not discriminating at those settings, and a clean
shortlist would be telling you nothing. Tighten the thresholds or accept
that the axes have no purchase here.

## 5. Read it, then test

Read the report in this order. It is not the order it is printed in by
accident.

1. **The evidence table.** Which sources produced evidence and which
   failed. A run where the lookups silently failed is otherwise
   structurally identical to one where they worked and found nothing.
2. **The funnel**, if you ran a panel. Where each gene lost its routes,
   and whether each stage removed anything a later stage would not have
   removed anyway.
3. **The routes.** Ordered by evidence tier first, then score — a tier-1
   route with a mediocre score beats a tier-5 route with a great one.
   Scores are comparable within a tier and not across one.
4. **The caveats under each table.** They are specific, not boilerplate.

Then the part the package cannot do for you:

```sh
python -m generoutes.explain results/ccr5/results.json maraviroc
```

`explain` puts everything known about one compound in one place, and the
report's experiment section describes the measurement that would settle
what the rest of this only predicts.

**Nothing here establishes that a drug changes your gene in your tissue.**
Every route is a hypothesis with its reasoning attached, which is what
makes it worth testing rather than worth believing.
