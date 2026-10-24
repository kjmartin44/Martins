# BUB1B: how to raise it - summary

Objective: **raise BUB1B**. Full detail, with every qualification, in `report.md`.

> **Only FDA approved drugs are shown.** 490 entries are left out. 223 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 14 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

**Evidence: degraded.** Read the ranking as a reading order, not as evidence. A pillar of the analysis is missing.
- these tests could not be completed: product_screen — which is not the same as their having returned a negative result
- About **46% of the 17 motif sites** (range 24% to 65%) would be expected on shuffled sequence of the same window. A site is a sequence match, not evidence of binding; ChIP peaks are what separate the real ones.

## What to test first

Grouped by what was actually measured, in 2 group(s), strongest kind of evidence first. **Within a group nothing is ranked** - the rows are ordered by how many independent methods found the drug, and beyond that alphabetically. That is deliberate: tested against CTD, this package's scoring orders drugs no better than chance once drug fame is held equal, so a numbered list would imply a precision it does not have. Use the groups, the method count and the exclusions; do not read position within a group as priority.

Drugs that reach the gene through the same factor by the same kind of route, and close analogues acting through a shared target, are folded into one row, so a family is one idea to test rather than several.

For BUB1B: one weak allele, one dead. The target is about 2x from the weak allele. Each sentence below says how that drug would move the gene. That is the mechanism the route assumes, and for most routes the step from the drug's target to the gene has not been measured. Nothing here shows that moving the gene treats the condition.

### Measured on the drug's target (6)

The compound's effect on its target is observed; the step from that target to your gene is inferred. Strong about mechanism, unverified about outcome.

| Drug | Approved to treat | Acts on | Tier | How it would work |
|---|---|---|---|---|
| **Bexarotene** (+5 similar) | cutaneous T-cell lymphoma | RXRA, RXRB, RXRG | 3 | Binds RXRA, a factor with a site at the BUB1B promoter whose role there is not known, so whether to activate or block it to raise BUB1B is not known. Its binding is measured; the step from RXRA to BUB1B is not. |
| **vorinostat** ⚠broad (+4 similar) | cutaneous T-cell lymphoma | HDAC1, HDAC3, NRL, chromatin (genome-wide) | 3 | Changes how much NRL the cell makes, which should raise BUB1B through the sites NRL has in its promoter. |
| **Obeticholic Acid** (+5 similar) | not established here | NR1H4 | 3 | Binds NR1H4, a factor with a site at the BUB1B promoter whose role there is not known, so whether to activate or block it to raise BUB1B is not known. Its binding is measured; the step from NR1H4 to BUB1B is not. |
| **dasatinib** ⚠databases disagree | not established here | ETS2, FOXA1, SOX17, TFAP2C | 3 | Changes how much TFAP2C the cell makes, which should raise BUB1B through the sites TFAP2C has in its promoter. |
| **Isotretinoin** | severe nodular acne | RARA | 4 | Acts on RARA, which is recorded as regulating BUB1B, so changing RARA activity should raise BUB1B. The record is from a database or other cells, not from yours. |
| **Raloxifene Hydrochloride** (+5 similar) | osteoporosis; breast cancer risk reduction | ESR2 | 4 | Acts on ESR2, which is recorded as regulating BUB1B, so changing ESR2 activity should raise BUB1B. The record is from a database or other cells, not from yours. |

### Predicted throughout (4)

Neither link is measured here. Worth testing only when something above has run out, or when a group above is empty for a reason the evidence table explains.

| Drug | Approved to treat | Acts on | Tier | How it would work |
|---|---|---|---|---|
| **Balsalazide Disodium** | not established here | PPARG | 4 | Acts on PPARG, upstream of factors that regulate BUB1B; the effect reaches BUB1B only by inference through the network. |
| **Olsalazine Sodium** | not established here | PPARG | 4 | Acts on PPARG, upstream of factors that regulate BUB1B; the effect reaches BUB1B only by inference through the network. |
| **Tazemetostat** | not established here | EZH2, chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |
| **Valproic acid** ⚠databases disagree | not established here | chromatin (genome-wide) | 3 | Changes chromatin across the genome, so it can move BUB1B together with thousands of other genes. |

### Found and excluded (6)

These were reached by the search and then rejected. They are shown because a reader who expected one of them should learn it was considered, not be left wondering whether the search missed it.

| Compound | Why it is not a candidate |
|---|---|
| Azacitidine | works by damaging DNA or the spindle, which is the defect being treated |
| dasatinib | ETS2 is a proto-oncogene and this route raises it |
| Decitabine | works by damaging DNA or the spindle, which is the defect being treated |
| mitoxantrone | works by damaging DNA or the spindle, which is the defect being treated |
| Paclitaxel | works by damaging DNA or the spindle, which is the defect being treated |
| thioridazine hydrochloride | FOXM1 is a proto-oncogene and this route raises it |

**Folded into the rows above** (the same factor reached the same way, or structure within Tanimoto 0.75 of the leader with a shared target). They are not duplicates: members of a family can differ in potency, selectivity, tissue reach and safety, and each has its own row in `report.md`.

- Obeticholic Acid: Odevixibat, Chenodiol, Ursodiol, Simvastatin, Fulvestrant
- Bexarotene: Acitretin, Etretinate, Levothyroxine, Sulindac, Tretinoin
- vorinostat: belinostat, Romidepsin, Panobinostat, Givinostat Hydrochloride
- Raloxifene Hydrochloride: Estramustine Phosphate Sodium, Toremifene Citrate, Chlorotrianisene, Ospemifene, Bazedoxifene Acetate

⚠databases disagree: CTD records Fulvestrant, Tretinoin, dasatinib, Valproic acid moving BUB1B the other way. A route with no recorded direction at its factor is lowered to tier 4; the others are only flagged, since CTD mixes cell types and doses. `report.md` lists them all.

⚠broad: the compound moves a large share of the transcriptome (top quartile of LINCS signatures), so lowering this gene is one of many changes it makes. A caution about specificity, not a score; the full table in `report.md` says how it is measured.

**Every row is off-label.** The third column is what the drug was approved for; using it to move your gene is a different purpose, and the distance between the two is part of what you are judging.

## The strongest 20 of 64 candidates, and which methods found each

One row per compound rather than per route, so a drug reached several ways is visible as one candidate with several supports.

| Drug | Methods | Tier | Acts on | What it is | Systems it acts on | Breadth | Reach | BBB? | BUB1B change | Engagement / promoter | Knockout evidence | Measured effect |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Bexarotene** | 4 | 3 | RXRA, RXRB, RXRG | RXR-selective agonist (rexinoid) | retinoid | — | systemic | yes? | increase (predicted) | 0.70 / 0.79 | not among top changes, 0 of 2 | — |
| **Acitretin** | 3 | 3 | RARA, RXRA, RXRB | retinoid | retinoid | few signatures | systemic | yes? | increase (predicted) | 0.70 / 0.79 | not among top changes, 0 of 2 | — |
| **Etretinate** | 3 | 3 | RARA, RXRA, RXRB | retinoid | retinoid | — | systemic | yes? | increase (predicted) | 0.70 / 0.79 | not among top changes, 0 of 2 | — |
| **Obeticholic Acid** | 2 | 3 | NR1H4 | binds NR1H4 (measured), and NR1H4… | no systemic class flagged | — | systemic | no | increase (predicted) | 0.74 / 0.51 | — | — |
| **vorinostat** | 3 | 3 | HDAC1, HDAC3, NRL | histone deacetylase inhibitor | chromatin-level agent | ⚠ broad (89th pct, 246 genes) | systemic | — | increase (predicted) | not estimated | not checked | — |
| **belinostat** | 3 | 3 | HDAC1, HDAC3, HDAC9 | lowers SOX17 transcript in 2… | chromatin-level agent | ⚠ broad (76th pct, 185 genes) | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Levothyroxine** | 2 | 3 | RXRA, THRB | synthetic T4 | thyroid hormone | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Chenodiol** | 1 | 3 | NR1H4 | ligand at the transcription factor | no systemic class flagged | — | systemic | no | increase (predicted) | 0.74 / 0.51 | — | — |
| **Ursodiol** | 1 | 3 | NR1H4 | ligand at the transcription factor | no systemic class flagged | — | systemic | no | increase (predicted) | 0.74 / 0.51 | — | — |
| **Sulindac** | 1 | 3 | RXRA | ligand at the transcription factor | no systemic class flagged | few signatures | systemic | — | increase (predicted) | 0.70 / 0.79 | not among top changes, 0 of 2 | — |
| **Romidepsin** | 2 | 3 | HDAC1, HDAC3, HDAC9 | histone deacetylase inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Panobinostat** | 2 | 3 | HDAC1, HDAC3, HDAC9 | HDAC inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Tazemetostat** | 2 | 3 | EZH2, chromatin (genome-wide) | EZH2 inhibitor | chromatin-level agent | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Givinostat Hydrochloride** | 2 | 3 | HDAC1, HDAC3, HDAC9 | HDAC inhibitor | no systemic class flagged | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **dasatinib** | 1 | 3 | ETS2, FOXA1, SOX17 | raises TFAP2C transcript in 7… | antineoplastic kinase… | 50th pct, 119 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Odevixibat** | 1 | 3 | NR1H4 | binds NR1H4 (measured), and NR1H4… | no systemic class flagged | — | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Simvastatin** | 1 | 3 | NR1H4 | HMG-CoA reductase inhibitor | no systemic class flagged | 38th pct, 96 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Valproic acid** | 1 | 3 | chromatin (genome-wide) | weak HDAC inhibitor | chromatin-level agent | few signatures | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Tretinoin** | 2 | 4 | RARA, RXRA | RAR agonist, all-trans retinoic… | retinoid | 31th pct, 84 genes | systemic | — | increase (predicted) | not estimated | not checked | — |
| **Fulvestrant** | 1 | 4 | NR1H4 | sex steroid class | sex steroid | 26th pct, 75 genes | systemic | — | increase (predicted) | not estimated | not checked | — |

**Breadth** is how many genes the compound moved across its LINCS signatures (up and down together), against all 33,132 signatures. ⚠ broad means the top quartile. It is a caution, not a verdict: it rises with dose and depends on the cell line, and the files hold set membership without scores, so they cannot say whether BUB1B stands out among the genes a compound changed. A compound moving hundreds of genes lowers BUB1B among many others. Compounds with fewer than 3 signatures are not judged.

**There is no predicted fold-change here, and that is not an omission.** Deriving one would need the share of BUB1B's transcription driven by that factor, the factor's occupancy change at the dose achieved, and the dose reached in your tissue. None of the three is available to this package, so any number in a fold-change column would be fabricated. The two numbers below are the ones it can actually compute.

**Engagement / promoter** are the two links that come from evidence: the drug engaging its target at a dose that can be reached, and the factor acting on this promoter. Each is a 0 to 1 heuristic and they differ from row to row. **Direction** and **Knockout evidence** appear only when a model predicted the direction or an experiment perturbed the factor. Order on tier and method count.

**12 of the 20 listed here were found by more than one method.** These methods share little - a curated ligand annotation, a measured expression signature, a regulon statistic, a pocket shape - so they fail for unrelated reasons, and agreement between them is worth more than a better score from any one.

- **Bexarotene**: ligand-based model (trained on this receptor); annotated ligand (ChEMBL/DGIdb/Open Targets); factor recorded as regulating the gene; upstream signalling network
- **Acitretin**: annotated ligand (ChEMBL/DGIdb/Open Targets); factor recorded as regulating the gene; upstream signalling network
- **Etretinate**: annotated ligand (ChEMBL/DGIdb/Open Targets); factor recorded as regulating the gene; upstream signalling network
- **Obeticholic Acid**: ligand-based model (trained on this receptor); annotated ligand (ChEMBL/DGIdb/Open Targets)
- **vorinostat**: measured change in the factor's own transcript; chromatin-level agent; factor recorded as regulating the gene
- **belinostat**: measured change in the factor's own transcript; chromatin-level agent; factor recorded as regulating the gene
- **Levothyroxine**: ligand-based model (trained on this receptor); upstream signalling network
- **Romidepsin**: chromatin-level agent; factor recorded as regulating the gene

**The method count is not a vote.** These methods differ enormously in strength, so two weak agreements do not equal one measured result. Tier still decides the order and the count breaks ties within a tier.

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

A **yes** predicts *passive diffusion only*, and is the weaker of the two answers. Active efflux by P-glycoprotein at the endothelium pumps an otherwise suitable molecule straight back out; known substrates are checked here and a hit overrides a pass outright.

**A question mark means undecided** — the score landed within noise of the threshold, and 3 of 20 row(s) here did. This is not pedantry. **Atenolol clears a complete five-parameter MPO by 0.04 and is a textbook poor penetrant**, so a margin that thin carries no information. Read a question-marked cell as "look this one up", never as the answer it leans toward.

**— is not a no.** 14 of 20 row(s) returned no physicochemical properties — common for DGIdb- and Open Targets-sourced rows, which carry a name and a mechanism but no structure. Nothing was measured and nothing is being claimed; look the compound up before excluding it on this column.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

**Reach** is the FDA route of administration. Every candidate here has a systemic route, so none is ruled out for being topical-only.

**The BUB1B change column is a direction.** `measured` means a compound moved BUB1B in a cell-based assay; `established` means it is an approved therapy acting on the gene or its product; `predicted` means it follows from the factor's role and the drug's action, with no measurement of BUB1B behind it.

## What the answer rests on

- **Ran, and found something:** chip, ctd, depmap, expression, factor_abundance, fame, genetics, knockout, lincs, proliferation, regulators, regulon, screen, sensitivity, upstream, utr
- **Ran, and found nothing:** direction_model, docking, promoter_dna, protein_abundance, trans_eqtl. That is a measurement, not a gap: these were asked and answered no.
- **Could not run:** product_screen. A source that broke is not a negative result, so nothing below should be read as evidence that it would have found nothing.

## Read these before the table

- The drug tables are dominated by nuclear receptors because 2 of 17 factors with a site here have any ligand at all. That is attrition, not a ranking.
- 15 of the best-scoring factors have no approved ligand (OSR2, THAP11, KLF3, SP3) and appear in no drug table. The gene-specific signal is partly in them.
- 6 factor(s) here have no measured occupancy because the ChIP table never covers them, not because they were looked for and not found - including NR1H4, which has a drug. They are scored on sequence alone and sit below factors with a peak for that reason, which is a fact about what has been assayed rather than about this promoter.

## Next step

Every row above is a hypothesis. The experiment that settles one is the same in each case: a promoter-reporter assay for the element, a dose-response with RT-qPCR or protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it. Without the third, a change in BUB1B does not establish the mechanism this report proposes.

*Not medical advice. Hypotheses for laboratory testing.*
