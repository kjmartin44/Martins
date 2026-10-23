# How the package is laid out

[Documentation index](index.md)

Eighty-seven modules, about forty-five thousand lines. This is the map: what the
run does in order, which module owns each step, and which parts you can
read on their own.

If you are looking for one thing in particular, the table at the bottom
lists every module with its own one-line description.

## The routes, and what they share

The package answers one question — which approved drugs plausibly move
this gene — by asking it several independent ways and putting the answers
on one scale. The routes do not feed each other. They read different data,
fail for unrelated reasons, and meet only in the catalogue.

```
gene symbol
  │
  ├─ the promoter motif route ──────────────────────────────────┐
  │    genes.py          resolve to a transcript and a TSS      │
  │    motifs.py, pwm.py scan the promoter with matrices        │
  │      → motif hits    most of these are chance; see below    │
  │    evidence.py       is the site real? ChIP over the window │
  │    annotate.py       is the factor an activator?            │
  │    neighbours.py     what else is transcribed here?         │
  │    drugs.py          approved drugs annotated at the factor │
  │                                                             │
  ├─ direct.py        drugs on the gene's own product ──────────┤
  ├─ regulators.py    drugs at its recorded regulators ─────────┤
  ├─ signatures.py    compounds measured to move the gene ──────┤
  ├─ network.py       further upstream: TF ← kinase ← receptor ─┤
  ├─ sensitivity.py   compounds that phenocopy losing a factor ─┤
  ├─ screening.py     ligand screens at the factor ─────────────┤
  ├─ docking.py       poses, once the pocket passes a control ──┤
  ├─ abundance.py     how much of the factor itself exists ─────┤
  ├─ stabilizers.py   hold the factor in its working shape ─────┤
  ├─ ligands.py       move its endogenous ligand instead ───────┤
  ├─ promoterdna.py   agents that bind the promoter itself ─────┤
  └─ discovery.py     what the literature already records ──────┤
                                                                ▼
                                                             routes
  │  catalogue.py      every route on one scale, with an evidence tier
  │  drug_scoring.py   grade each route on factor and drug confidence
  ▼
ranked routes
  │  shortlist.py      narrow one gene's result
  │  panel.py          narrow across several candidate genes
  ▼
a table you can act on
```

The honest part of that diagram is the note next to the motif hits. A
promoter scan at p ≤ 1e-4 over a kilobase returns far more sites than
exist, and the report says so with a number: for CCR5, 241 hits against
187 expected by chance. Everything on that branch inherits it, which is
why so much of the package is evidence that a site is real rather than
predicted, and why a noisy scan now caps that branch's tiers.

That weakness belongs to the motif route and not to the others. The routes that rest on a measurement of the gene — the literature
record, the LINCS signature, the sensitivity phenocopy — carry better
evidence that the gene moves, and worse evidence about why.

## The groups

### Entry points

| Module | What it is |
| --- | --- |
| `cli.py` | the single-gene command, `python -m generoutes` |
| `panel.py` | several candidate genes, narrowed through a staged funnel |
| `pipeline.py` | `run()`, the orchestration everything else feeds |
| `shortlist.py` | re-narrow a saved `results.json`, no re-run |
| `explain.py` | everything known about one compound, in one place |
| `audit.py` | are the hand-curated tables still true? |
| `validation.py` | what `run()` accepts, as a table rather than conditions |

`pipeline.run` is the one function to read first. It is also the largest
file in the package by a wide margin, which is a known problem.

### Finding the factors

| Module | What it is |
| --- | --- |
| `genes.py` | gene symbol to transcription start site and promoter sequence |
| `motifs.py`, `pwm.py` | load binding motifs, scan a sequence, score a match |
| `annotate.py` | which matched factors are plausible activators |
| `evidence.py` | orthogonal evidence that a predicted site is real |
| `neighbours.py` | what else is transcribed in the window you scanned |
| `chipregulon.py` | regulons built from ChIP peaks, for factors curation missed |
| `regulon.py` | how many genes a factor controls — the specificity that counts |

### Finding the drugs

| Module | What it is |
| --- | --- |
| `drugs.py` | factor to approved drug, filtered to the right direction |
| `localdrugs.py` | the same from a file instead of an API |
| `direct.py` | drugs acting on the gene itself, skipping factors entirely |
| `network.py` | multi-hop: promoter ← TF ← kinase ← receptor ← drug |
| `ligands.py` | change the natural ligand rather than the receptor |
| `abundance.py` | drugs that change how much factor exists, not what it does |
| `stabilizers.py` | locking a factor active without a ligand pocket |
| `modalities.py` | design a direct modality rather than looking one up |
| `analogs.py` | approved drugs resembling a known modulator |
| `posttranscriptional.py`, `reactivation.py`, `compensation.py` | routes that are not transcriptional at all |

### Measured evidence

These are the modules that carry something someone actually observed, as
opposed to something inferred from annotation.

| Module | What it is |
| --- | --- |
| `signatures.py` | LINCS L1000 — a compound's measured effect on this gene |
| `depmap.py` | what happens to a cell when the factor is removed |
| `sensitivity.py` | PRISM: drugs behaving as though they act on the factor |
| `tfactivity.py` | which factors are demonstrably active, from their targets |
| `eqtl.py`, `genetics.py` | does human variation in the factor move the gene |
| `corpora.py` | other perturbation corpora, and what agreement is worth |
| `localpotency.py` | measured binding potency from a file (BindingDB, ChEMBL) |
| `localannot.py` | annotation and genetic evidence from files, not APIs |

### Structure

| Module | What it is |
| --- | --- |
| `docking.py` | dock a library into a pocket, after checking the pocket works |
| `wsldock.py` | run a Linux docking binary from Windows, through WSL |
| `screening.py` | screen the approved library against a receptor |
| `agonism.py` | predict whether a compound activates or inhibits |
| `promoterdna.py` | molecules aimed at the promoter DNA rather than a protein |

See [the docking guide](docking.md) for how these fit together and what a
docking score is and is not worth.

### Collateral: what else the drug moves

Four axes, each a separate module, described in full in [specificity and
evidence](specificity.md).

| Module | What it is |
| --- | --- |
| `delivery.py` | does the drug reach the cells that express your gene |
| `systemic.py` | drugs and receptors whose action is systemic by design |
| `paralogs.py` | how many close relatives share the pocket |
| `cellcycle.py` | phase dependence, a separate axis from the others |
| `liabilities.py` | what makes a predicted hit not worth believing |
| `indications.py`, `orangebook.py` | what a drug is and what it treats |
| `congeners.py` | collapse salt forms so the table shows distinct chemistry |

### Scoring and reporting

| Module | What it is |
| --- | --- |
| `catalogue.py` | every route in one table, on one scale, with a tier |
| `drug_scoring.py` | grade a factor/drug pair on both kinds of confidence |
| `predictions.py` | the headline answer and its reasoning |
| `discovery.py` | find routes that plausibly beat what already exists |
| `outcome.py` | from "the drug binds" to "my gene moves": the chain, explicitly |
| `modulators.py` | orthosteric against allosteric, a comparison not a ranking |
| `screens.py` | design the experiment that would settle it |

### Infrastructure

| Module | What it is |
| --- | --- |
| `http.py` | the cached client. Stdlib only, so NumPy stays the only dependency |
| `provenance.py` | did each evidence source actually produce evidence |
| `demo_data.py` | fixtures so the demo runs with no network at all |

## Two ideas worth knowing before you read the code

**A failed lookup is not a measured negative.** This is the single
constraint that shapes the most code. When a source cannot be reached,
the result records that it was requested and failed — it does not record
zero. `provenance.py` holds that ledger, every report opens with it, and
the evidence verdict degrades when a pillar is missing. A run where the
lookups silently failed would otherwise be structurally identical to one
where they worked and found nothing.

**Evidence tiers, not scores, decide the order.** A tier-1 route with a
mediocre score beats a tier-5 route with a great one, because the tiers
describe what kind of claim is being made: clinical, measured, curated,
predicted, hypothesis. Scores are comparable within a tier and not
across one. `catalogue.TIERS` is the definition.

## Everything, alphabetically

Generated from each module's own first docstring line.

| Module | Purpose |
| --- | --- |
| `abundance.py` | Drugs that change how much of the factor exists, not what it does. |
| `agonism.py` | Predict the *functional direction* of a drug at a receptor from its structure. |
| `analogs.py` | Approved drugs that resemble a known modulator but are not annotated at it. |
| `annotate.py` | Decide which of the motif-matched TFs are plausible *activators*. |
| `audit.py` | Are the hand-curated tables still true? |
| `catalogue.py` | Every route that acts on your gene, in one table, on one scale. |
| `cellcycle.py` | Cell-cycle phase dependence of TF activity — a separate axis from everything else. |
| `chemblid.py` | Find the ChEMBL id of a drug that arrived without one. |
| `chipregulon.py` | Regulons built from ChIP peaks, for the factors curation never reached. |
| `cli.py` | Command-line interface. |
| `clusters.py` | Fold near-identical drugs into one row of a ranking. |
| `combined.py` | One ranked list from docking and the rest of the analysis. |
| `compensation.py` | You may not need to raise your gene at all. |
| `congeners.py` | Collapse salt forms so the tables show distinct chemistry. |
| `constraint.py` | How much dosage this gene tolerates, in either direction. |
| `corpora.py` | Other perturbation corpora, and what agreement between them is worth. |
| `coverage.py` | Which of the six steps to a protein this run actually looked at. |
| `crosscheck.py` | Where CTD and LINCS speak about the same drug, and where either contradicts a route. |
| `datasets.py` | Every table the analysis can read: where it comes from, and whether you have it. |
| `delivery.py` | Does the drug reach the cells that express your gene? |
| `demo_data.py` | Self-contained fixtures so the pipeline runs with zero network access. |
| `depmap.py` | What happens to a cell when a factor is removed, measured not curated. |
| `direct.py` | Drugs that act on the target gene directly, bypassing transcription factors. |
| `discovery.py` | Find interventions that plausibly BEAT what already exists for this gene. |
| `docking.py` | Structure-based prediction of functional direction, done the way that can |
| `drug_scoring.py` | Grade each TF/drug pair by combining TF-level confidence with drug quality. |
| `drugs.py` | TF -> approved drug lookup, filtered to pharmacology that moves the TF in |
| `eqtl.py` | Trans-eQTL evidence: does variation in the factor actually move your gene? |
| `evidence.py` | Orthogonal evidence that a predicted site is real. |
| `explain.py` | Everything the pipeline knows about one compound, in one place. |
| `exploratory.py` | Exploratory docking into a pocket that has no co-crystallised ligand. |
| `exposure.py` | Can the drug reach the concentration at which the evidence was seen? |
| `fame.py` | How many other genes the literature records a compound as moving. |
| `fda.py` | Show only drugs the FDA has approved. |
| `fdalabel.py` | Human peak plasma concentrations read from FDA drug labels. |
| `genes.py` | Resolve a gene symbol to a transcription start site and pull the promoter. |
| `genetics.py` | Human genetic evidence, via eQTLs — and what GWAS/PheWAS can and cannot tell you. |
| `http.py` | Tiny cached HTTP client. Stdlib only, so the package needs only numpy. |
| `indications.py` | What each drug actually is, and what it is approved to treat. |
| `isoforms.py` | How much message each isoform actually carries. |
| `knockout.py` | Does your gene change when a factor is lost? Read from perturbation data. |
| `liabilities.py` | Structural liabilities: what makes a predicted hit not worth believing. |
| `ligands.py` | Change the natural ligand, not the receptor. |
| `localannot.py` | Gene annotation and genetic evidence from files, not web APIs. |
| `localdrugs.py` | Approved drugs and their known targets, from a file instead of an API. |
| `localpotency.py` | Measured binding potency from a file, not from the ChEMBL API. |
| `magnitude.py` | How much the gene has to move, and which routes have ever been measured. |
| `modalities.py` | Design a direct modality for YOUR gene, rather than looking up existing ones. |
| `modulators.py` | Orthosteric agonists vs allosteric modulators: a comparison, not a ranking. |
| `motifs.py` | Load TF binding motifs and scan a promoter with them. |
| `neighbours.py` | What else lives in the promoter window you just scanned? |
| `network.py` | Multi-hop upstream search: promoter -> TF -> kinase -> receptor -> drug. |
| `orangebook.py` | The FDA's 20,000 approved products, used for what they actually add. |
| `outcome.py` | From "the drug binds the receptor" to "my gene goes up": the chain, explicitly. |
| `panel.py` | Compare several candidate genes and narrow to the cleanest route. |
| `paralogs.py` | How many close relatives share the pocket you are trying to hit? |
| `pipeline.py` | Orchestration: run every route for one gene and rank what they return. |
| `poisonexon.py` | Exons that only decay-prone isoforms use: candidate brakes on output. |
| `posttranscriptional.py` | Raising a gene's output after transcription has already happened. |
| `predictions.py` | The headline answer: which medications the model predicts, and why. |
| `progress.py` | Telling the user what the run is doing while it does it. |
| `proliferation.py` | Did the compound move the gene, or just move the cell cycle? |
| `promoterdna.py` | Small molecules aimed at the promoter DNA rather than at a protein. |
| `proteomics.py` | Drug to PROTEIN abundance, where the transcript side has LINCS. |
| `provenance.py` | Did each evidence source actually produce evidence? |
| `pwm.py` | Position-weight-matrix mathematics. |
| `rationale.py` | A plain sentence on how a candidate would act on the condition. |
| `reactivation.py` | Reactivating a silenced allele, which is a different problem from a low gene. |
| `regulators.py` | Start from the factors known to regulate your gene, then look for drugs. |
| `regulon.py` | How many genes does this TF control? The specificity that actually matters. |
| `report.py` | Turning a finished run into something a person reads. |
| `retrospective.py` | Does the pipeline recover drugs already known to move a gene? |
| `screening.py` | Screen the whole approved-drug library against a receptor. |
| `screenlist.py` | Every compound worth putting on a plate, with nothing ranked. |
| `screens.py` | Design the experiment that settles what this whole package only predicts. |
| `sensitivity.py` | Which approved drugs behave as though they act on your factor. |
| `services.py` | Do the web services still return the shape the parsers expect? |
| `shortlist.py` | Narrow an existing results.json to a shortlist, without re-running anything. |
| `signatures.py` | LINCS L1000 transcriptomic signatures — the only empirical link in the package. |
| `stabilizers.py` | Locking a transcription factor into an active state — without a ligand pocket. |
| `standard_of_care.py` | Approved treatments for the condition a gene causes. |
| `stats.py` | Statistics shared by more than one analysis stage. |
| `strata.py` | What to test first, grouped by what was actually measured. |
| `systemic.py` | Drugs and receptors whose action is systemic by design, not by accident. |
| `tfactivity.py` | Which factors are demonstrably active, inferred from their targets. |
| `utr.py` | Translation and mRNA stability, read off the transcript's own UTRs. |
| `validation.py` | Declarative argument checking for the entry points. |
| `wsldock.py` | Run a Linux docking binary from Windows, through WSL. |

The descriptions come from the modules themselves.
`scripts/generate_module_table.py` regenerates this table, and a test
fails when it is stale.
