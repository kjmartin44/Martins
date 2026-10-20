# CYP24A1: how to raise it - summary

Objective: **raise CYP24A1**. Full detail, with every qualification, in `report.md`.

> **Only FDA approved drugs are shown.** 1129 entries are left out. 232 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 22 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

**Evidence: degraded.** Read the ranking as a reading order, not as evidence. A pillar of the analysis is missing.
- the motif scan is 'weak' — about 73% of the hit list is expected to be chance (range 55% to 92%)

## Top 10 to look at first

Ordered by evidence tier, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. Close analogues that act through the same target are folded into the best-ranked one, so a family of near identical drugs takes one row.

The aim is to raise CYP24A1. Each sentence below says how that drug would move the gene. That is the mechanism the route assumes, and for most routes the step from the drug's target to the gene has not been measured. Nothing here shows that moving the gene treats the condition.

| # | Drug | Approved to treat | Acts on | Tier | How it would work |
|---|---|---|---|---|---|
| 1 | **teniposide** | not established here | STAT3 | 3 | Changes how much STAT3 the cell makes, which should raise CYP24A1 through the sites STAT3 has in its promoter. |
| 2 | **selumetinib** | neurofibromatosis type 1 plexiform neurofibroma | STAT3 | 3 | Changes how much STAT3 the cell makes, which should raise CYP24A1 through the sites STAT3 has in its promoter. |
| 3 | **Cholecalciferol** (+1 similar) | vitamin D deficiency and supplementation | VDR, not stated | 3 | Acts on VDR, which is recorded as regulating CYP24A1, so changing VDR activity should raise CYP24A1. The record is from a database or other cells, not from yours. |
| 4 | **Paricalcitol** | secondary hyperparathyroidism in chronic kidney disease | CYP24A1, VDR, not stated | 3 | Acts on VDR, which is recorded as regulating CYP24A1, so changing VDR activity should raise CYP24A1. The record is from a database or other cells, not from yours. |
| 5 | **Doxercalciferol** (+1 similar) | secondary hyperparathyroidism in dialysis | VDR | 3 | Acts on VDR, which is recorded as regulating CYP24A1, so changing VDR activity should raise CYP24A1. The record is from a database or other cells, not from yours. |
| 6 | **Calcitriol** | hypocalcaemia and secondary hyperparathyroidism in renal… | VDR, not stated | 3 | Acts on VDR, which is recorded as regulating CYP24A1, so changing VDR activity should raise CYP24A1. The record is from a database or other cells, not from yours. |
| 7 | **dasatinib** | not established here | STAT3 | 3 | Changes how much STAT3 the cell makes, which should raise CYP24A1 through the sites STAT3 has in its promoter. |
| 8 | **vorinostat** ⚠broad | cutaneous T-cell lymphoma | FOSL2, HNF4A, JUND, STAT3, chromatin (genome-wide) | 3 | Changes how much FOSL2 the cell makes, which should raise CYP24A1 through the sites FOSL2 has in its promoter. |
| 9 | **mitoxantrone** | acute myeloid leukaemia, prostate cancer, multiple… | FOS, FOSL2, JUNB, JUND, KLF5, STAT3, promoter DNA | 3 | Changes how much FOS the cell makes, which should raise CYP24A1 through the sites FOS has in its promoter. |
| 10 | **vemurafenib** | not established here | SOX2 | 3 | Changes how much SOX2 the cell makes, which should raise CYP24A1 through the sites SOX2 has in its promoter. |

**Folded into the rows above** (structure within Tanimoto 0.75 of the leader, same target). They are not duplicates: analogues can differ in potency, selectivity, tissue reach and safety.

- 3. Cholecalciferol: Calcifediol
- 5. Doxercalciferol: Ergocalciferol

⚠broad: the compound moves a large share of the transcriptome (top quartile of LINCS signatures), so lowering this gene is one of many changes it makes. A caution about specificity, not a score; the full table in `report.md` says how it is measured.

**Every row is off-label.** The third column is what the drug was approved for; using it to move your gene is a different purpose, and the distance between the two is part of what you are judging.

## All 20 candidates, and which methods found each

One row per compound rather than per route, so a drug reached several ways is visible as one candidate with several supports.

| Drug | Methods | Tier | Acts on | What it is | Systems it acts on | Breadth | Reach | BBB? | CYP24A1 change | Engagement / promoter | Measured effect |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Cholecalciferol** | 2 | 3 | VDR, not stated | vitamin D3, a precursor requiring… | vitamin D receptor agonist | — | systemic | — | increase (predicted) | not estimated | — |
| **Paricalcitol** | 2 | 3 | CYP24A1, VDR, not stated | selective VDR activator | vitamin D receptor agonist | — | systemic | — | increase (predicted) | not estimated | — |
| **Calcifediol** | 2 | 3 | VDR, not stated | 25-hydroxyvitamin D | vitamin D receptor agonist | — | systemic | — | increase (predicted) | not estimated | — |
| **Calcitriol** | 2 | 3 | VDR, not stated | active vitamin D, VDR agonist | vitamin D receptor agonist | 20th pct, 66 genes | systemic | — | increase (predicted) | not estimated | — |
| **vorinostat** | 2 | 3 | FOSL2, HNF4A, JUND | histone deacetylase inhibitor | chromatin-level agent | ⚠ broad (89th pct, 246 genes) | systemic | — | increase (predicted) | not estimated | — |
| **teniposide** | 1 | 3 | STAT3 | raises STAT3 transcript in 10… | no systemic class flagged | 72th pct, 174 genes | systemic | — | increase (predicted) | not estimated | — |
| **selumetinib** | 1 | 3 | STAT3 | MEK1/2 inhibitor | antineoplastic kinase… | 60th pct, 142 genes | systemic | — | increase (predicted) | not estimated | — |
| **dasatinib** | 1 | 3 | STAT3 | raises STAT3 transcript in 7… | antineoplastic kinase… | 50th pct, 119 genes | systemic | — | increase (predicted) | not estimated | — |
| **mitoxantrone** | 1 | 3 | FOS, FOSL2, JUNB | topoisomerase II inhibitor, DNA… | cytotoxic /… | 74th pct, 180 genes | systemic | — | increase (predicted) | not estimated | — |
| **vemurafenib** | 1 | 3 | SOX2 | raises SOX2 transcript in 5… | antineoplastic kinase… | 48th pct, 114 genes | systemic | — | increase (predicted) | not estimated | — |
| **Valproic acid** | 2 | 3 | chromatin (genome-wide), not stated | weak HDAC inhibitor | chromatin-level agent | few signatures | systemic | — | increase (predicted) | not estimated | — |
| **Doxercalciferol** | 1 | 3 | VDR | vitamin D prodrug | vitamin D receptor agonist | — | systemic | — | increase (predicted) | not estimated | — |
| **Ergocalciferol** | 1 | 3 | VDR | vitamin D receptor agonist class | vitamin D receptor agonist | — | systemic | — | increase (predicted) | not estimated | — |
| **Romidepsin** | 1 | 3 | chromatin (genome-wide) | histone deacetylase inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | — |
| **Belinostat** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | chromatin-level agent | ⚠ broad (76th pct, 185 genes) | systemic | — | increase (predicted) | not estimated | — |
| **Panobinostat** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | — |
| **Azacitidine** | 1 | 3 | chromatin (genome-wide) | DNA methyltransferase inhibitor | chromatin-level agent | 41th pct, 101 genes | systemic | — | increase (predicted) | not estimated | — |
| **Decitabine** | 1 | 3 | chromatin (genome-wide) | DNA methyltransferase inhibitor | chromatin-level agent | 21th pct, 67 genes | systemic | — | increase (predicted) | not estimated | — |
| **Tazemetostat** | 1 | 3 | chromatin (genome-wide) | EZH2 inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | — |
| **Givinostat Hydrochloride** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | no systemic class flagged | — | systemic | — | increase (predicted) | not estimated | — |

**Breadth** is how many genes the compound moved across its LINCS signatures (up and down together), against all 33,132 signatures. ⚠ broad means the top quartile. It is a caution, not a verdict: it rises with dose and depends on the cell line, and the files hold set membership without scores, so they cannot say whether CYP24A1 stands out among the genes a compound changed. A compound moving hundreds of genes lowers CYP24A1 among many others. Compounds with fewer than 3 signatures are not judged.

**There is no predicted fold-change here, and that is not an omission.** Deriving one would need the share of CYP24A1's transcription driven by that factor, the factor's occupancy change at the dose achieved, and the dose reached in your tissue. None of the three is available to this package, so any number in a fold-change column would be fabricated. The two numbers below are the ones it can actually compute.

**6 of 20 were found by more than one method.** These methods share little - a curated ligand annotation, a measured expression signature, a regulon statistic, a pocket shape - so they fail for unrelated reasons, and agreement between them is worth more than a better score from any one.

- **Cholecalciferol**: factor recorded as regulating the gene; literature record of an effect (CTD)
- **Paricalcitol**: factor recorded as regulating the gene; literature record of an effect (CTD)
- **Calcifediol**: factor recorded as regulating the gene; literature record of an effect (CTD)
- **Calcitriol**: factor recorded as regulating the gene; literature record of an effect (CTD)
- **vorinostat**: measured change in the factor's own transcript; chromatin-level agent
- **Valproic acid**: chromatin-level agent; literature record of an effect (CTD)

**The method count is not a vote.** These methods differ enormously in strength, so two weak agreements do not equal one measured result. Tier still decides the order and the count breaks ties within a tier.

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

**— is not a no.** 20 of 20 row(s) returned no physicochemical properties — common for DGIdb- and Open Targets-sourced rows, which carry a name and a mechanism but no structure. Nothing was measured and nothing is being claimed; look the compound up before excluding it on this column.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

**Reach** is the FDA route of administration. Every candidate here has a systemic route, so none is ruled out for being topical-only.

**The CYP24A1 change column is a direction.** `measured` means a compound moved CYP24A1 in a cell-based assay; `established` means it is an approved therapy acting on the gene or its product; `predicted` means it follows from the factor's role and the drug's action, with no measurement of CYP24A1 behind it.

## What the answer rests on

- **Ran:** chip, ctd, depmap, expression, factor_abundance, genetics, knockout, lincs, product_screen, promoter_dna, regulators, regulon, sensitivity, upstream
- **Did not run:** direction_model, screen, trans_eqtl. A source that failed is not a negative result; `report.md` names the flag that would enable each.

## Read these before the table

- 25 of the best-scoring factors have no approved ligand (SOX2, POU2F1, FOS, JUNB) and appear in no drug table. The gene-specific signal is partly in them.
- The promoter windows overlap each other, so counts across them are not independent.

## Next step

Every row above is a hypothesis. The experiment that settles one is the same in each case: a promoter-reporter assay for the element, a dose-response with RT-qPCR or protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it. Without the third, a change in CYP24A1 does not establish the mechanism this report proposes.

*Not medical advice. Hypotheses for laboratory testing.*
