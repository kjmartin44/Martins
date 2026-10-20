# CYP24A1: how to raise it - summary

Objective: **raise CYP24A1**. Full detail, with every qualification, in `report.md`.

> **Only FDA approved drugs are shown.** 798 entries are left out. 262 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 42 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

**Evidence: degraded.** Read the ranking as a reading order, not as evidence. A pillar of the analysis is missing.
- the motif scan is 'no signal' — about 86% of the hit list is expected to be chance (range 46% to 100%)

## Top 10 to look at first

Ordered by evidence tier, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. Close analogues that act through the same target are folded into the best-ranked one, so a family of near identical drugs takes one row.

The aim is to raise CYP24A1. Each sentence below says how that drug would move the gene. That is the mechanism the route assumes, and for most routes the step from the drug's target to the gene has not been measured. Nothing here shows that moving the gene treats the condition.

| # | Drug | Approved to treat | Acts on | Tier | How it would work |
|---|---|---|---|---|---|
| 1 | **Dimethyl fumarate** (+1 similar) | not established here | KEAP1 | 3 | Changes how long KEAP1 lasts in the cell (degradation blocker), which should raise CYP24A1. |
| 2 | **Omaveloxolone** | not established here | KEAP1, NFE2L2 | 3 | Changes how long KEAP1 lasts in the cell (degradation blocker), which should raise CYP24A1. |
| 3 | **vorinostat** ⚠broad | cutaneous T-cell lymphoma | ELK1, NFE2L2, NFIB, chromatin (genome-wide) | 3 | Changes how much NFIB the cell makes, which should raise CYP24A1 through the sites NFIB has in its promoter. |
| 4 | **dasatinib** | not established here | CEBPD, TBX2 | 3 | Changes how much CEBPD the cell makes, which should raise CYP24A1 through the sites CEBPD has in its promoter. |
| 5 | **mitoxantrone** | acute myeloid leukaemia, prostate cancer, multiple… | CEBPA, promoter DNA | 3 | Changes how much CEBPA the cell makes, which should raise CYP24A1 through the sites CEBPA has in its promoter. |
| 6 | **vemurafenib** | not established here | SOX2 | 3 | Changes how much SOX2 the cell makes, which should raise CYP24A1 through the sites SOX2 has in its promoter. |
| 7 | **sirolimus** | transplant rejection prophylaxis, lymphangioleiomyomatosis | TBX2 | 3 | Changes how much TBX2 the cell makes, which should raise CYP24A1 through the sites TBX2 has in its promoter. |
| 8 | **selumetinib** | neurofibromatosis type 1 plexiform neurofibroma | TBX2 | 3 | Changes how much TBX2 the cell makes, which should raise CYP24A1 through the sites TBX2 has in its promoter. |
| 9 | **Romidepsin** | cutaneous and peripheral T-cell lymphoma | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move CYP24A1 together with thousands of other genes. |
| 10 | **Belinostat** ⚠broad | not established here | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move CYP24A1 together with thousands of other genes. |

**Folded into the rows above** (structure within Tanimoto 0.75 of the leader, same target). They are not duplicates: analogues can differ in potency, selectivity, tissue reach and safety.

- 1. Dimethyl fumarate: Monomethyl fumarate

⚠broad: the compound moves a large share of the transcriptome (top quartile of LINCS signatures), so lowering this gene is one of many changes it makes. A caution about specificity, not a score; the full table in `report.md` says how it is measured.

**Every row is off-label.** The third column is what the drug was approved for; using it to move your gene is a different purpose, and the distance between the two is part of what you are judging.

## All 20 candidates, and which methods found each

One row per compound rather than per route, so a drug reached several ways is visible as one candidate with several supports.

| Drug | Methods | Tier | Acts on | What it is | Systems it acts on | Breadth | Reach | BBB? | CYP24A1 change | Engagement / promoter | Knockout evidence | Measured effect |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Omaveloxolone** | 2 | 3 | KEAP1, NFE2L2 | degradation blocker | no systemic class flagged | — | systemic | no | increase (predicted) | 0.78 / 0.56 | gene rises when the factor is lost, 1 of 27 | — |
| **vorinostat** | 2 | 3 | ELK1, NFE2L2, NFIB | histone deacetylase inhibitor | chromatin-level agent | ⚠ broad (89th pct, 246 genes) | systemic | — | increase (predicted) | not estimated | not checked | — |
| **dasatinib** | 1 | 3 | CEBPD, TBX2 | raises CEBPD transcript in 51… | antineoplastic kinase… | 50th pct, 119 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **mitoxantrone** | 1 | 3 | CEBPA, promoter DNA | topoisomerase II inhibitor, DNA… | cytotoxic /… | 74th pct, 180 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **vemurafenib** | 1 | 3 | SOX2 | raises SOX2 transcript in 5… | antineoplastic kinase… | 48th pct, 114 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **sirolimus** | 1 | 3 | TBX2 | mTOR inhibitor | broad immunosuppressant | 43th pct, 105 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **selumetinib** | 1 | 3 | TBX2 | MEK1/2 inhibitor | antineoplastic kinase… | 60th pct, 142 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Dimethyl fumarate** | 1 | 3 | KEAP1 | degradation blocker | no systemic class flagged | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Monomethyl fumarate** | 1 | 3 | KEAP1 | degradation blocker | no systemic class flagged | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Romidepsin** | 1 | 3 | chromatin (genome-wide) | histone deacetylase inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Belinostat** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | chromatin-level agent | ⚠ broad (76th pct, 185 genes) | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Panobinostat** | 1 | 3 | chromatin (genome-wide) | HDAC inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Valproic acid** | 1 | 3 | chromatin (genome-wide) | weak HDAC inhibitor | chromatin-level agent | few signatures | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Azacitidine** | 1 | 3 | chromatin (genome-wide) | DNA methyltransferase inhibitor | chromatin-level agent | 41th pct, 101 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Decitabine** | 1 | 3 | chromatin (genome-wide) | DNA methyltransferase inhibitor | chromatin-level agent | 21th pct, 67 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Tazemetostat** | 1 | 3 | chromatin (genome-wide) | EZH2 inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Progesterone** | 1 | 4 | PGR | sex steroid class | sex steroid | 21th pct, 67 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Norethindrone** | 1 | 4 | PGR | sex steroid class | sex steroid | 46th pct, 111 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Ethynodiol Diacetate** | 1 | 4 | PGR | activate PGR → FOXO1 (1 hop) | no systemic class flagged | — | unclear | — | increase (predicted) | not estimated | not checked | — |
| **Hydroxyprogesterone Caproate** | 1 | 4 | PGR | activate PGR → FOXO1 (1 hop) | no systemic class flagged | — | systemic | — | increase (predicted) | not estimated | not checked | — |

**Breadth** is how many genes the compound moved across its LINCS signatures (up and down together), against all 33,132 signatures. ⚠ broad means the top quartile. It is a caution, not a verdict: it rises with dose and depends on the cell line, and the files hold set membership without scores, so they cannot say whether CYP24A1 stands out among the genes a compound changed. A compound moving hundreds of genes lowers CYP24A1 among many others. Compounds with fewer than 3 signatures are not judged.

**There is no predicted fold-change here, and that is not an omission.** Deriving one would need the share of CYP24A1's transcription driven by that factor, the factor's occupancy change at the dose achieved, and the dose reached in your tissue. None of the three is available to this package, so any number in a fold-change column would be fabricated. The two numbers below are the ones it can actually compute.

**Engagement / promoter** are the two links that come from evidence: the drug engaging its target at a dose that can be reached, and the factor acting on this promoter. Each is a 0 to 1 heuristic and they differ from row to row. **Direction** and **Knockout evidence** appear only when a model predicted the direction or an experiment perturbed the factor. Order on tier and method count.

**2 of 20 were found by more than one method.** These methods share little - a curated ligand annotation, a measured expression signature, a regulon statistic, a pocket shape - so they fail for unrelated reasons, and agreement between them is worth more than a better score from any one.

- **Omaveloxolone**: protein stability / degradation; annotated ligand (ChEMBL/DGIdb/Open Targets)
- **vorinostat**: measured change in the factor's own transcript; chromatin-level agent

**The method count is not a vote.** These methods differ enormously in strength, so two weak agreements do not equal one measured result. Tier still decides the order and the count breaks ties within a tier.

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

**— is not a no.** 19 of 20 row(s) returned no physicochemical properties — common for DGIdb- and Open Targets-sourced rows, which carry a name and a mechanism but no structure. Nothing was measured and nothing is being claimed; look the compound up before excluding it on this column.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

**Reach** is the FDA route of administration. Every candidate here has a systemic route, so none is ruled out for being topical-only.

**The CYP24A1 change column is a direction.** `measured` means a compound moved CYP24A1 in a cell-based assay; `established` means it is an approved therapy acting on the gene or its product; `predicted` means it follows from the factor's role and the drug's action, with no measurement of CYP24A1 behind it.

## What the answer rests on

- **Ran:** chip, ctd, depmap, expression, factor_abundance, genetics, knockout, lincs, product_screen, promoter_dna, regulators, regulon, screen, sensitivity, upstream
- **Did not run:** direction_model, trans_eqtl. A source that failed is not a negative result; `report.md` names the flag that would enable each.

## Read these before the table

- 14 of the best-scoring factors have no approved ligand (SOX2, POU5F1, POU2F1, ZNF530) and appear in no drug table. The gene-specific signal is partly in them.

## Next step

Every row above is a hypothesis. The experiment that settles one is the same in each case: a promoter-reporter assay for the element, a dose-response with RT-qPCR or protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it. Without the third, a change in CYP24A1 does not establish the mechanism this report proposes.

*Not medical advice. Hypotheses for laboratory testing.*
