# DEMOGENE: how to raise it - summary

Objective: **raise DEMOGENE**. Demo data - not a real analysis. Full detail, with every qualification, in `report.md`.

> **FDA filter not applied:** not applied to demo data, whose drugs are synthetic

**Evidence: degraded.** Read the ranking as a reading order, not as evidence. A pillar of the analysis is missing.
- no source of measured evidence is available (ChIP, LINCS, genetics all absent or failed)
- About **14% of the 14 motif sites** (range 0% to 29%) would be expected on shuffled sequence of the same window. A site is a sequence match, not evidence of binding; ChIP peaks are what separate the real ones.

## Top 10 to look at first

Ordered by evidence tier, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. Close analogues that act through the same target are folded into the best-ranked one, so a family of near identical drugs takes one row.

The aim is to raise DEMOGENE. Each sentence below says how that drug would move the gene. That is the mechanism the route assumes, and for most routes the step from the drug's target to the gene has not been measured. Nothing here shows that moving the gene treats the condition.

| # | Drug | Approved to treat | Acts on | Tier | How it would work |
|---|---|---|---|---|---|
| 1 | **Calcitriol** | hypocalcaemia and secondary hyperparathyroidism in renal… | VDR | 2 | Binds VDR, a factor predicted to switch on DEMOGENE; activating VDR should raise DEMOGENE. Its binding is measured; the step from VDR to DEMOGENE is not. |
| 2 | **Paricalcitol** | secondary hyperparathyroidism in chronic kidney disease | VDR | 2 | Binds VDR, a factor predicted to switch on DEMOGENE; activating VDR should raise DEMOGENE. Its binding is measured; the step from VDR to DEMOGENE is not. |
| 3 | **Dexamethasone** | inflammation, cerebral oedema, myeloma, and as an… | NR3C1 | 2 | Binds NR3C1, a factor predicted to act on DEMOGENE; activating NR3C1 should raise DEMOGENE. Its binding is measured; the step from NR3C1 to DEMOGENE is not. |
| 4 | **Rosiglitazone** | type 2 diabetes; restricted for cardiovascular risk | PPARG | 2 | Binds PPARG, a factor predicted to switch on DEMOGENE; activating PPARG should raise DEMOGENE. Its binding is measured; the step from PPARG to DEMOGENE is not. |
| 5 | **Tretinoin** | acute promyelocytic leukaemia; acne topically | RARA, — | 2 | Binds RARA, a factor predicted to switch on DEMOGENE; activating RARA should raise DEMOGENE. Its binding is measured; the step from RARA to DEMOGENE is not. |
| 6 | **Bexarotene** | cutaneous T-cell lymphoma | RXRA | 2 | Binds RXRA, a factor predicted to switch on DEMOGENE; activating RXRA should raise DEMOGENE. Its binding is measured; the step from RXRA to DEMOGENE is not. |
| 7 | **Prednisolone** | asthma, autoimmune and inflammatory disease | NR3C1 | 2 | Binds NR3C1, a factor predicted to act on DEMOGENE; activating NR3C1 should raise DEMOGENE. Its binding is measured; the step from NR3C1 to DEMOGENE is not. |
| 8 | **Alitretinoin** | Kaposi sarcoma skin lesions; chronic hand eczema | RARA | 2 | Binds RARA, a factor predicted to switch on DEMOGENE; activating RARA should raise DEMOGENE. Its binding is measured; the step from RARA to DEMOGENE is not. |
| 9 | **Hydrocortisone** | adrenal insufficiency; topical dermatitis | NR3C1, — | 2 | Binds NR3C1, a factor predicted to act on DEMOGENE; activating NR3C1 should raise DEMOGENE. Its binding is measured; the step from NR3C1 to DEMOGENE is not. |
| 10 | **Pioglitazone** | type 2 diabetes | PPARG | 2 | Binds PPARG, a factor predicted to switch on DEMOGENE; activating PPARG should raise DEMOGENE. Its binding is measured; the step from PPARG to DEMOGENE is not. |

**Every row is off-label.** The third column is what the drug was approved for; using it to move your gene is a different purpose, and the distance between the two is part of what you are judging.

## All 20 candidates, and which methods found each

One row per compound rather than per route, so a drug reached several ways is visible as one candidate with several supports.

| Drug | Methods | Tier | Acts on | What it is | Systems it acts on | Breadth | Reach | BBB? | DEMOGENE change | Engagement / promoter | Measured effect |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Tretinoin** | 2 | 2 | RARA, — | RAR agonist, all-trans retinoic… | retinoid | — | — | — | increase (predicted) | 0.92 / 0.68 | — |
| **Hydrocortisone** | 2 | 2 | NR3C1, — | glucocorticoid and… | glucocorticoid | — | — | — | increase (predicted) | 0.79 / 0.57 | — |
| **Calcitriol** | 1 | 2 | VDR | active vitamin D, VDR agonist | vitamin D receptor agonist | — | — | no | increase (predicted) | 0.95 / 0.59 | — |
| **Paricalcitol** | 1 | 2 | VDR | selective VDR activator | vitamin D receptor agonist | — | — | — | increase (predicted) | 0.95 / 0.59 | — |
| **Dexamethasone** | 1 | 2 | NR3C1 | potent glucocorticoid receptor… | glucocorticoid | — | — | no | increase (predicted) | 0.89 / 0.57 | — |
| **Rosiglitazone** | 1 | 2 | PPARG | PPARgamma agonist… | PPAR agonist | — | — | — | increase (predicted) | 0.86 / 0.65 | — |
| **Bexarotene** | 1 | 2 | RXRA | RXR-selective agonist (rexinoid) | retinoid | — | — | — | increase (predicted) | 0.85 / 0.68 | — |
| **Prednisolone** | 1 | 2 | NR3C1 | glucocorticoid | glucocorticoid | — | — | — | increase (predicted) | 0.82 / 0.57 | — |
| **Alitretinoin** | 1 | 2 | RARA | RAR/RXR pan-agonist | retinoid | — | — | — | increase (predicted) | 0.88 / 0.68 | — |
| **Pioglitazone** | 1 | 2 | PPARG | PPARgamma agonist | PPAR agonist | — | — | — | increase (predicted) | 0.74 / 0.65 | — |
| **Valproic Acid** | 2 | 3 | REST, chromatin (genome-wide) | ligand at the transcription factor | chromatin-level agent | — | — | yes | increase (predicted) | 0.55 / 0.57 | — |
| **Demo-Pam-2** | 1 | 3 | VDR | ligand at the transcription factor | no systemic class flagged | — | — | — | increase (predicted) | 0.78 / 0.59 | — |
| **Demo-Pam-1** | 1 | 3 | PPARG | PAM | no systemic class flagged | — | — | — | increase (predicted) | 0.78 / 0.65 | — |
| **Carbenoxolone** | 1 | 3 | HSD11B2 | inhibit catabolism (cortisol) via… | no systemic class flagged | — | — | — | increase (predicted) | not estimated | — |
| **Calcifediol** | 1 | 3 | — | 25-hydroxyvitamin D | vitamin D receptor agonist | — | — | — | increase (predicted) | not estimated | — |
| **Ketoconazole** | 1 | 3 | CYP24A1 | inhibit catabolism (calcitriol)… | no systemic class flagged | — | — | — | increase (predicted) | not estimated | — |
| **Alfacalcidol** | 1 | 3 | — | vitamin D prodrug | vitamin D receptor agonist | — | — | — | increase (predicted) | not estimated | — |
| **Cholecalciferol** | 1 | 3 | — | vitamin D3, a precursor requiring… | vitamin D receptor agonist | — | — | — | increase (predicted) | not estimated | — |
| **Vorinostat** | 1 | 3 | chromatin (genome-wide) | histone deacetylase inhibitor | chromatin-level agent | — | — | — | increase (predicted) | not estimated | — |
| **Romidepsin** | 1 | 3 | chromatin (genome-wide) | histone deacetylase inhibitor | chromatin-level agent | — | — | — | increase (predicted) | not estimated | — |

**There is no predicted fold-change here, and that is not an omission.** Deriving one would need the share of DEMOGENE's transcription driven by that factor, the factor's occupancy change at the dose achieved, and the dose reached in your tissue. None of the three is available to this package, so any number in a fold-change column would be fabricated. The two numbers below are the ones it can actually compute.

**Engagement / promoter** are the two links that come from evidence: the drug engaging its target at a dose that can be reached, and the factor acting on this promoter. Each is a 0 to 1 heuristic and they differ from row to row. **Direction** and **Knockout evidence** appear only when a model predicted the direction or an experiment perturbed the factor. Order on tier and method count.

**3 of 20 were found by more than one method.** These methods share little - a curated ligand annotation, a measured expression signature, a regulon statistic, a pocket shape - so they fail for unrelated reasons, and agreement between them is worth more than a better score from any one.

- **Tretinoin**: annotated ligand (ChEMBL/DGIdb/Open Targets); endogenous ligand availability
- **Hydrocortisone**: annotated ligand (ChEMBL/DGIdb/Open Targets); endogenous ligand availability
- **Valproic Acid**: annotated ligand (ChEMBL/DGIdb/Open Targets); chromatin-level agent

**The method count is not a vote.** These methods differ enormously in strength, so two weak agreements do not equal one measured result. Tier still decides the order and the count breaks ties within a tier.

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

A **yes** predicts *passive diffusion only*, and is the weaker of the two answers. Active efflux by P-glycoprotein at the endothelium pumps an otherwise suitable molecule straight back out; known substrates are checked here and a hit overrides a pass outright.

**— is not a no.** 17 of 20 row(s) returned no physicochemical properties — common for DGIdb- and Open Targets-sourced rows, which carry a name and a mechanism but no structure. Nothing was measured and nothing is being claimed; look the compound up before excluding it on this column.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

**Reach** is blank for every row because no Orange Book is on disk. That is unknown route, not systemic - a topical-only drug looks identical. See `orangebook.SOURCE` for the download.

**The DEMOGENE change column is a direction.** `measured` means a compound moved DEMOGENE in a cell-based assay; `established` means it is an approved therapy acting on the gene or its product; `predicted` means it follows from the factor's role and the drug's action, with no measurement of DEMOGENE behind it.

## Read these before the table

- The drug tables are dominated by nuclear receptors because 6 of 10 factors with a site here have any ligand at all. That is attrition, not a ranking.
- 4 of the best-scoring factors have no approved ligand (SP1, CREB1, NFKB1, HNF4A) and appear in no drug table. The gene-specific signal is partly in them.

## Next step

Every row above is a hypothesis. The experiment that settles one is the same in each case: a promoter-reporter assay for the element, a dose-response with RT-qPCR or protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it. Without the third, a change in DEMOGENE does not establish the mechanism this report proposes.

*Not medical advice. Hypotheses for laboratory testing.*
