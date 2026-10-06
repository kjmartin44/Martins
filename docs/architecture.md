# How the package is laid out

[Documentation index](index.md)

Sixty-one modules, about thirty thousand lines. This is the map: what the
run does in order, which module owns each step, and which parts you can
read on their own.

If you are looking for one thing in particular, the table at the bottom
lists every module with its own one-line description.

## The chain, in order

The package answers one question — which approved drugs plausibly move
this gene — by walking a chain. Every stage narrows or annotates what the
stage before produced.

```
gene symbol
  │  genes.py          resolve to a transcript and a TSS, pull the promoter
  ▼
promoter sequence
  │  motifs.py, pwm.py scan it with position-weight matrices
  ▼
motif hits                      ← most of these are chance; see below
  │  evidence.py       is the site real? ChIP peaks over the same window
  │  annotate.py       is the factor an activator or a repressor?
  │  neighbours.py     what else is transcribed in this window?
  ▼
candidate factors
  │  drugs.py          approved drugs annotated at each factor
  │  network.py        and drugs further upstream: TF ← kinase ← receptor
  │  direct.py         and drugs acting on the gene's own product
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
187 expected by chance. Everything downstream inherits that, which is why
so much of the package is evidence that a site is real rather than
predicted.

## The groups

### Entry points

| Module | What it is |
| --- | --- |
| `cli.py` | the single-gene command, `python -m promoter_tf_drugs` |
| `panel.py` | several candidate genes, narrowed through a staged funnel |
| `pipeline.py` | `run()`, the orchestration everything else feeds |
| `shortlist.py` | re-narrow a saved `results.json`, no re-run |
| `explain.py` | everything known about one compound, in one place |
| `audit.py` | are the hand-curated tables still true? |
| `validation.py` | what `run()` accepts, as a table rather than conditions |

`pipeline.run` is the one function to read first. It is also the largest
file in the package by a wide margin, which is a known problem and the
first item under maintainability in the [improvement
review](improvements.md).

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

| Module | Lines | Purpose |
| --- | --- | --- |
| `abundance.py` | 325 | Drugs that change how much of the factor exists, not what it does. |
| `agonism.py` | 416 | Predict the *functional direction* of a drug at a receptor from its structure. |
| `analogs.py` | 210 | Approved drugs that resemble a known modulator but are not annotated at it. |
| `annotate.py` | 199 | Decide which of the motif-matched TFs are plausible *activators*. |
| `audit.py` | 343 | Are the hand-curated tables still true? |
| `catalogue.py` | 1263 | Every route that acts on your gene, in one table, on one scale. |
| `cellcycle.py` | 247 | Cell-cycle phase dependence of TF activity — a separate axis from everything else. |
| `chipregulon.py` | 419 | Regulons built from ChIP peaks, for the factors curation never reached. |
| `cli.py` | 411 | Command-line interface. |
| `compensation.py` | 229 | You may not need to raise your gene at all. |
| `congeners.py` | 159 | Collapse salt forms so the tables show distinct chemistry. |
| `corpora.py` | 415 | Other perturbation corpora, and what agreement between them is worth. |
| `delivery.py` | 791 | Does the drug reach the cells that express your gene? |
| `demo_data.py` | 228 | Self-contained fixtures so the pipeline runs with zero network access. |
| `depmap.py` | 433 | What happens to a cell when a factor is removed, measured not curated. |
| `direct.py` | 230 | Drugs that act on the target gene directly, bypassing transcription factors. |
| `discovery.py` | 516 | Find interventions that plausibly BEAT what already exists for this gene. |
| `docking.py` | 2577 | Structure-based prediction of functional direction, done the way that can |
| `drug_scoring.py` | 454 | Grade each TF/drug pair by combining TF-level confidence with drug quality. |
| `drugs.py` | 492 | TF -> approved drug lookup, filtered to pharmacology that moves the TF in |
| `eqtl.py` | 266 | Trans-eQTL evidence: does variation in the factor actually move your gene? |
| `evidence.py` | 341 | Orthogonal evidence that a predicted site is real. |
| `explain.py` | 289 | Everything the pipeline knows about one compound, in one place. |
| `genes.py` | 564 | Resolve a gene symbol to a transcription start site and pull the promoter. |
| `genetics.py` | 212 | Human genetic evidence, via eQTLs — and what GWAS/PheWAS can and cannot tell you. |
| `http.py` | 255 | Tiny cached HTTP client. Stdlib only, so the package needs only numpy. |
| `indications.py` | 397 | What each drug actually is, and what it is approved to treat. |
| `liabilities.py` | 229 | Structural liabilities: what makes a predicted hit not worth believing. |
| `ligands.py` | 306 | Change the natural ligand, not the receptor. |
| `localannot.py` | 378 | Gene annotation and genetic evidence from files, not web APIs. |
| `localdrugs.py` | 971 | Approved drugs and their known targets, from a file instead of an API. |
| `localpotency.py` | 563 | Measured binding potency from a file, not from the ChEMBL API. |
| `modalities.py` | 330 | Design a direct modality for YOUR gene, rather than looking up existing ones. |
| `modulators.py` | 301 | Orthosteric agonists vs allosteric modulators: a comparison, not a ranking. |
| `motifs.py` | 296 | Load TF binding motifs and scan a promoter with them. |
| `neighbours.py` | 286 | What else lives in the promoter window you just scanned? |
| `network.py` | 366 | Multi-hop upstream search: promoter -> TF -> kinase -> receptor -> drug. |
| `orangebook.py` | 324 | The FDA's 20,000 approved products, used for what they actually add. |
| `outcome.py` | 181 | From "the drug binds the receptor" to "my gene goes up": the chain, explicitly. |
| `panel.py` | 1201 | Compare several candidate genes and narrow to the cleanest route. |
| `paralogs.py` | 229 | How many close relatives share the pocket you are trying to hit? |
| `pipeline.py` | 4649 | Orchestration: gene -> promoter -> motifs -> activators -> approved drugs. |
| `posttranscriptional.py` | 229 | Raising a gene's output after transcription has already happened. |
| `predictions.py` | 456 | The headline answer: which medications the model predicts, and why. |
| `promoterdna.py` | 342 | Small molecules aimed at the promoter DNA rather than at a protein. |
| `provenance.py` | 233 | Did each evidence source actually produce evidence? |
| `pwm.py` | 161 | Position-weight-matrix mathematics. |
| `reactivation.py` | 247 | Reactivating a silenced allele, which is a different problem from a low gene. |
| `regulon.py` | 252 | How many genes does this TF control? The specificity that actually matters. |
| `screening.py` | 1035 | Screen the whole approved-drug library against a receptor. |
| `screens.py` | 294 | Design the experiment that settles what this whole package only predicts. |
| `sensitivity.py` | 631 | Which approved drugs behave as though they act on your factor. |
| `shortlist.py` | 575 | Narrow an existing results.json to a shortlist, without re-running anything. |
| `signatures.py` | 575 | LINCS L1000 transcriptomic signatures — the only empirical link in the package. |
| `stabilizers.py` | 448 | Locking a transcription factor into an active state — without a ligand pocket. |
| `systemic.py` | 360 | Drugs and receptors whose action is systemic by design, not by accident. |
| `tfactivity.py` | 969 | Which factors are demonstrably active, inferred from their targets. |
| `validation.py` | 173 | Declarative argument checking for the entry points. |
| `wsldock.py` | 324 | Run a Linux docking binary from Windows, through WSL. |

Line counts are indicative and go stale; the descriptions come from the
modules themselves and `docs/development.md` explains how to regenerate
this table.
