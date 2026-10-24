# BUB1B: how to raise it - summary

Objective: **raise BUB1B**. Full detail, with every qualification, in `report.md`.

> **Only FDA approved drugs are shown.** 16 entries are left out. 11 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 4 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

**Evidence: interpretable.**
- About **46% of the 17 motif sites** (range 24% to 65%) would be expected on shuffled sequence of the same window. A site is a sequence match, not evidence of binding; ChIP peaks are what separate the real ones.

## What to test first

Grouped by what was actually measured, in 1 group(s), strongest kind of evidence first. **Within a group nothing is ranked** - the rows are ordered by how many independent methods found the drug, and beyond that alphabetically. That is deliberate: tested against CTD, this package's scoring orders drugs no better than chance once drug fame is held equal, so a numbered list would imply a precision it does not have. Use the groups, the method count and the exclusions; do not read position within a group as priority.

Drugs that reach the gene through the same factor by the same kind of route, and close analogues acting through a shared target, are folded into one row, so a family is one idea to test rather than several.

The aim is to raise BUB1B. Each sentence below says how that drug would move the gene. That is the mechanism the route assumes, and for most routes the step from the drug's target to the gene has not been measured. Nothing here shows that moving the gene treats the condition.

### Predicted throughout (9)

Neither link is measured here. Worth testing only when something above has run out, or when a group above is empty for a reason the evidence table explains.

| Drug | Approved to treat | Acts on | Tier | How it would work |
|---|---|---|---|---|
| **Azacitidine** | myelodysplastic syndrome, acute myeloid leukaemia | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Belinostat** ⚠broad | not established here | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Decitabine** | myelodysplastic syndrome | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Givinostat Hydrochloride** | not established here | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Panobinostat** | not established here | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Romidepsin** | cutaneous and peripheral T-cell lymphoma | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Tazemetostat** | not established here | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Valproic acid** | not established here | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Vorinostat** ⚠broad | cutaneous T-cell lymphoma | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |

⚠broad: the compound moves a large share of the transcriptome (top quartile of LINCS signatures), so lowering this gene is one of many changes it makes. A caution about specificity, not a score; the full table in `report.md` says how it is measured.

**Every row is off-label.** The third column is what the drug was approved for; using it to move your gene is a different purpose, and the distance between the two is part of what you are judging.

## All 9 candidates, and which methods found each

One row per compound rather than per route, so a drug reached several ways is visible as one candidate with several supports.

| Drug | Methods | Tier | Acts on | What it is | Systems it acts on | Breadth | Reach | BBB? | BUB1B change | Engagement / promoter | Measured effect |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Vorinostat** | 1 | 3 | chromatin (genome-wide) | histone deacetylase inhibitor | chromatin-level agent | ⚠ broad (89th pct, 246 genes) | systemic | — | increase (predicted) | not estimated | — |
| **Romidepsin** | 1 | 3 | chromatin (genome-wide) | histone deacetylase inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | — |
| **Belinostat** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | chromatin-level agent | ⚠ broad (76th pct, 185 genes) | systemic | — | increase (predicted) | not estimated | — |
| **Panobinostat** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | — |
| **Azacitidine** | 1 | 3 | chromatin (genome-wide) | DNA methyltransferase inhibitor | chromatin-level agent | 41th pct, 101 genes | systemic | — | increase (predicted) | not estimated | — |
| **Decitabine** | 1 | 3 | chromatin (genome-wide) | DNA methyltransferase inhibitor | chromatin-level agent | 21th pct, 67 genes | systemic | — | increase (predicted) | not estimated | — |
| **Tazemetostat** | 1 | 3 | chromatin (genome-wide) | EZH2 inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | — |
| **Givinostat Hydrochloride** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | no systemic class flagged | — | systemic | — | increase (predicted) | not estimated | — |
| **Valproic acid** | 1 | 3 | chromatin (genome-wide) | weak HDAC inhibitor | chromatin-level agent | few signatures | systemic | — | increase (predicted) | not estimated | — |

**Breadth** is how many genes the compound moved across its LINCS signatures (up and down together), against all 33,132 signatures. ⚠ broad means the top quartile. It is a caution, not a verdict: it rises with dose and depends on the cell line, and the files hold set membership without scores, so they cannot say whether BUB1B stands out among the genes a compound changed. A compound moving hundreds of genes lowers BUB1B among many others. Compounds with fewer than 3 signatures are not judged.

**There is no predicted fold-change here, and that is not an omission.** Deriving one would need the share of BUB1B's transcription driven by that factor, the factor's occupancy change at the dose achieved, and the dose reached in your tissue. None of the three is available to this package, so any number in a fold-change column would be fabricated. The two numbers below are the ones it can actually compute.

**No compound was found by more than one method.** With 1 method(s) contributing, that is expected rather than informative: corroboration needs two methods to have run and both to reach the same compound. `report.md` lists the flags that turn the others on.

**The method count is not a vote.** These methods differ enormously in strength, so two weak agreements do not equal one measured result. Tier still decides the order and the count breaks ties within a tier.

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

**— is not a no.** 9 of 9 row(s) returned no physicochemical properties — common for DGIdb- and Open Targets-sourced rows, which carry a name and a mechanism but no structure. Nothing was measured and nothing is being claimed; look the compound up before excluding it on this column.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

**Reach** is the FDA route of administration. Every candidate here has a systemic route, so none is ruled out for being topical-only.

**The BUB1B change column is a direction.** `measured` means a compound moved BUB1B in a cell-based assay; `established` means it is an approved therapy acting on the gene or its product; `predicted` means it follows from the factor's role and the drug's action, with no measurement of BUB1B behind it.

## What the answer rests on

- **Ran, and found something:** chip, depmap, fame, knockout, regulon, utr
- **Ran, and found nothing:** factor_abundance, proliferation, protein_abundance. That is a measurement, not a gap: these were asked and answered no.

## Read these before the table

- The drug tables are dominated by nuclear receptors because 2 of 17 factors with a site here have any ligand at all. That is attrition, not a ranking.
- 15 of the best-scoring factors have no approved ligand (OSR2, THAP11, KLF3, SP3) and appear in no drug table. The gene-specific signal is partly in them.
- 6 factor(s) here have no measured occupancy because the ChIP table never covers them, not because they were looked for and not found - including NR1H4, which has a drug. They are scored on sequence alone and sit below factors with a peak for that reason, which is a fact about what has been assayed rather than about this promoter.

## Next step

Every row above is a hypothesis. The experiment that settles one is the same in each case: a promoter-reporter assay for the element, a dose-response with RT-qPCR or protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it. Without the third, a change in BUB1B does not establish the mechanism this report proposes.

*Not medical advice. Hypotheses for laboratory testing.*
