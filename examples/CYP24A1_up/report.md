# Upregulation strategy for CYP24A1 — promoter / TF / drug analysis

Objective: **increase** transcription of CYP24A1.

> **Only FDA approved drugs are shown.** 1129 entries are left out. 232 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 22 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

## What would have to happen

**CYP24A1 dosage: tolerant.** LOEUF 1.51: the population carries about as many loss-of-function variants as chance predicts, so losing ONE copy is tolerated. That is not evidence the gene is unimportant. A recessive condition is invisible to this measure entirely - carriers are healthy, so selection never saw the phenotype. Mosaic variegated aneuploidy works this way: biallelic BUB1B, with unaffected heterozygous parents. One copy of this gene is dispensable in the population, so if your condition is recessive the relevant question is how much function the variant alleles retain - not what the population tolerates.

*Read from gnomAD constraint v4.1; ClinGen dosage curation 07 Oct,2026.*

## Evidence available in this run

**⚠ DEGRADED** — Read the ranking as a reading order, not as evidence. At least one pillar of the analysis is missing. Measured evidence that DID arrive: chip, lincs, ctd, genetics.

- the motif scan is 'weak' — about 73% of the hit list is expected to be chance (range 55% to 92%)

| Source | State | What it returned |
|---|---|---|
| chip | ok | local BED: peaks for 421 factor(s) over the window |
| ctd | ok | 180 chemical(s) with curated expression effects on this gene; 85 move it the way you want |
| depmap | ok | 80 of 81 factor(s) scored in CRISPRGeneEffect.csv |
| direction_model | ok, nothing found | 0 of 3 grade(s) got a model call; too few labelled agonists and antagonists to train |
| expression | ok | 16 tissue(s) for CYP24A1 |
| factor_abundance | ok | 48 compound(s) move a factor's own transcript the needed way; 1 factor(s) skipped for an unrecorded role |
| genetics | ok | 13 associated phenotype(s) (Open Targets); 206 genome-wide significant GWAS Catalog association(s), mapped or nearest gene, which is weaker than fine-mapped support |
| knockout | ok | 41 of 81 factor(s) were perturbed in the knockout data; the gene moved for 6 |
| lincs | ok | 3 signature hit(s) for this gene |
| product_screen | ok | 1 unannotated candidate(s) at the CYP24A1 protein, cross-validated AUC 1.0 |
| promoter_dna | ok | 4 sequence feature(s) enriched over 200 shuffles; 12 agent(s) prefer them |
| regulators | ok | 170 recorded regulator(s) of the gene had no drug route yet; 7 have an approved drug that moves them the needed way (48 route(s)) |
| regulon | ok | 590 TFs with targets called from ChIP occupancy. A peak near a gene is binding, not regulation |
| screen | ok, nothing found | no factor has an approved drug that names its target, so there was no receptor to screen |
| sensitivity | ok | 7 correlated compound(s) across 6 of 6 factor(s) scanned |
| trans_eqtl | ok, nothing found | 0 of 6 factor(s) have trans-eQTL support |
| upstream | ok | druggable routes for 7 factor(s); 3668 of 7138 reachable upstream node(s) checked across 8 factor(s) |

## Drugs targeting CYP24A1 directly

None found in the curated table or in ChEMBL/DGIdb for this direction. The curated list is hand-maintained and goes out of date; check FDA/EMA and clinicaltrials.gov before concluding nothing exists.

## Exons that may brake CYP24A1 output

Some genes make part of their transcript with an extra exon that puts a stop codon early, so the message is destroyed. An antisense oligonucleotide that makes the cell skip that exon turns those transcripts into working ones and raises the protein, from the working copy only. Ensembl lists no nonsense-mediated-decay isoform of CYP24A1, so there is no annotated brake exon to target. That depends on the annotation: a local 'basic' GTF omits most of them.

**This is a list of places to look, not a design.** An annotated exon is one Ensembl has seen in some tissue. Whether it is included in your cells, how often, and whether skipping it raises protein is a question for RNA sequencing and a splicing assay. Nothing here picks an oligonucleotide sequence or checks it against the rest of the transcriptome.

## Model-predicted medications for CYP24A1

Compounds this analysis predicts would **raise** CYP24A1, with the reasoning for each. Everything here is a **hypothesis generated from annotation**, not a recommendation — the tier on each row tells you how far the claim is reaching.

> **⚠ Read as a reading order, not as evidence.** The run's evidence verdict is *degraded*: at least one pillar of the analysis is missing. See the evidence table above.

| # | Medication | What it is | Approved to treat | Acts on | BBB? | Tier | Grade | Score |
|---|---|---|---|---|---|---|---|---|
| 1 | **Cholecalciferol** | vitamin D3, a precursor requiring activation | vitamin D deficiency and supplementation | VDR | — | 3 | — | 45.2 |
| 2 | **Paricalcitol** | selective VDR activator | secondary hyperparathyroidism in chronic kidney disease | VDR | — | 3 | — | 45.2 |
| 3 | **Doxercalciferol** | vitamin D prodrug | secondary hyperparathyroidism in dialysis | VDR | — | 3 | — | 45.2 |
| 4 | **Ergocalciferol** | vitamin D receptor agonist class | secondary hyperparathyroidism in renal disease, osteoporosis, and plaque psoriasis topically | VDR | — | 3 | — | 45.2 |
| 5 | **Calcifediol** | 25-hydroxyvitamin D | vitamin D deficiency | VDR | — | 3 | — | 45.2 |
| 6 | **Calcitriol** | active vitamin D, VDR agonist | hypocalcaemia and secondary hyperparathyroidism in renal failure | VDR | — | 3 | — | 45.2 |
| 7 | **vorinostat** | histone deacetylase inhibitor | cutaneous T-cell lymphoma | FOSL2 | — | 3 | — | 44.8 |
| 8 | **mitoxantrone** | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis | FOS | — | 3 | — | 44.8 |
| 9 | **mitoxantrone** | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis | JUNB | — | 3 | — | 44.8 |
| 10 | **mitoxantrone** | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis | FOSL2 | — | 3 | — | 44.8 |

**Every row is off-label.** The drug was approved for what the fourth column says, and is listed here for something else. That column is the baseline this proposal departs from, and the distance between it and your gene is part of what you are judging - a topical acne gel and an infused cytotoxic can sit on adjacent rows with similar scores.

*Ordered by evidence tier first, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. **BBB?** is not part of the ordering — it is an applicability filter you apply yourself, since whether brain exposure is wanted or unwanted depends on where your gene needs to move. See the note under the graded table for what a yes does and does not establish.*

### Why each one

Every route below is **evidence tier 3** — a human curator annotated this mechanism.

**The direction rule**, which applies to every group below: you want to raise the gene, so you agonise an activator and block a repressor.

#### Through **VDR** — 6 candidates

- The factor has no activator/repressor annotation, so the direction of its effect on your gene is assumed, not known.
- Family: VDR shares its structurally diverse lipophiles pocket with NR1I2, NR1I3, so a ligand here will tend to engage them too — though subtype selectivity has been achieved (calcitriol (VDR, does not meaningfully engage PXR/CAR)).

- ⚠ 1 knockout experiment(s) agree and 0 disagree
- ⚠ VDR raises the gene according to CollecTRI and knockout data

1. **Cholecalciferol**
2. **Paricalcitol**
3. **Doxercalciferol**
4. **Ergocalciferol**
5. **Calcifediol**
6. **Calcitriol**

#### Through **FOSL2** — 2 candidates

- Role: **activator**
- Binding: 4 predicted site(s), the best scoring 1.0 -4812 bp from the TSS; **4 ChIP peak(s) confirm it actually binds here** — the strongest link in this chain.

- ⚠ measured against FOSL2's own transcript, not against your gene. That FOSL2 controls your gene is the inference this route rests on

7. **vorinostat** — raises FOSL2 transcript in 88 profile(s)
   - ⚠ moves 4 of the 7 factors asked about, so it is not a regulator of any one of them. Compounds that shut down transcription or proteostasis broadly - HSP90 and CDK inhibitors above all - produce this pattern
8. **mitoxantrone** — raises FOSL2 transcript in 21 profile(s)
   - ⚠ moves 5 of the 7 factors asked about, so it is not a regulator of any one of them. Compounds that shut down transcription or proteostasis broadly - HSP90 and CDK inhibitors above all - produce this pattern

#### Through **FOS**

- Role: **activator**
- Binding: 3 predicted site(s), the best scoring 1.0 -4811 bp from the TSS; **1 ChIP peak(s) confirm it actually binds here** — the strongest link in this chain.

9. **mitoxantrone**
   - ⚠ moves 5 of the 7 factors asked about, so it is not a regulator of any one of them. Compounds that shut down transcription or proteostasis broadly - HSP90 and CDK inhibitors above all - produce this pattern

#### Through **JUNB**

- Role: **activator**
- Binding: 5 predicted site(s), the best scoring 1.0 -4811 bp from the TSS; **5 ChIP peak(s) confirm it actually binds here** — the strongest link in this chain.

10. **mitoxantrone**
   - ⚠ moves 5 of the 7 factors asked about, so it is not a regulator of any one of them. Compounds that shut down transcription or proteostasis broadly - HSP90 and CDK inhibitors above all - produce this pattern

*Not medical advice. "Approved" means approved for its own indication — using one of these to move an unrelated gene is off-label repurposing with the full systemic pharmacology attached. The experimental sections at the end of this report are how you would actually test any of it.*

> **⚠ These promoter windows overlap each other.** these two TSS windows share 3966 bp (66%) — largely one analysis counted twice, which inflates this gene's route count Counts and route lists below are therefore not independent across the windows.

## Promoter of CYP24A1 (ENST00000395954)

- Transcript chosen: **distinct TSS**
- Assembly: from GTF  |  Strand: -
- TSS: chr20:54171953
- Window: chr20:54170954-54176953 (-5000 / +1000, 6000 bp)
- Motifs scanned: 879 at p ≤ 1e-05 per position
- **The window was widened.** The requested window (-1000/+200) showed no signal above chance (13 sites), so it was rescanned at -5000/+1000, which gave 89 sites against about 65.0 expected by chance. Sites further out say less about this gene than core promoter sites do. `--no-widen` keeps the requested window.
- 322 of those motifs were not scanned: in this window their best possible match is more common than p = 1e-05, so any hit would be a chance match at a higher rate than stated. A looser `--pvalue` scans them, with more chance hits.
- Non-redundant sites: 89 across 81 TFs

> **⚠ Weak scan.** 89 sites against 65 expected by chance over 4600 scannable bases. Most of this list is noise; treat only the highest-scoring sites with ChIP support as real.

*about 73% of these 89 sites (range 55% to 92%) would be expected on shuffled sequence of this window, so roughly 24 are above chance; among the best-scoring sites the share is lower: 100% of the top 5, 83% of the top 10, 100% of the top 20. Estimated from 30 shuffled copies that keep the repeat mask and every dinucleotide count. The formula that assumes each motif matches at exactly its stated p-value would expect 51 and is not used. Nothing here says which sites are real; only ChIP evidence separates them.*

#### What else is in this window

Every other specificity check here asks whether a drug or factor affects other genes. This asks the prior question: is the sequence that was scanned uniquely this gene's?

The window overlaps 1 other transcription unit(s), including ENSG00000286587 on the opposite strand. Sites called here are not unambiguously your gene's, and for the antisense overlap the direction logic does not carry.

| Gene | Biotype | Location | Strand | Relation |
|---|---|---|---|---|
| ENSG00000286587 | lncRNA | chr20:54059120-54190210 | + | antisense overlap |

- **ENSG00000286587** — window lies inside this gene on the OPPOSITE strand. A site called here may drive that transcript instead, and an antisense RNA can stabilise or repress your gene — so the pipeline's direction logic does not reliably carry

*An antisense overlap matters more than it looks. The direction logic in this report assumes a factor acts on YOUR promoter; where the window lies inside an antisense transcript, a site may drive that instead, and antisense RNA can stabilise or repress the sense gene — so the sign does not carry.*

*This checks annotation, not sequence uniqueness. BLAT the window if you need to know it appears nowhere else in the genome.*

#### Why the drug tables below look the way they do

**81 factor(s) had predicted sites here. 0 of them have an approved drug (0 of those are nuclear receptors).**

Roughly 48 of the ~1,600 human transcription factors have a ligand-binding pocket, and the steroid, retinoid, thyroid and vitamin-D families are most of them. So the drug tables are not a ranking that hormone drugs won — they are what is left after every factor without a ligand has been dropped.

**25 of the top 25 factors by score have no drug and appear in no drug table:** SOX2, POU2F1, FOS, JUNB, JUND, FOSL2, JUN, KLF5, STAT3, NFYB, FOSB, CEBPA ….

That is where the gene-specific signal is. A half-site like AGGTCA occurs every few kilobases by chance, so nuclear receptors surface for almost any promoter; these factors surfaced for *this* one. The routes that can reach them are the upstream-network, stabilisation and protein-level sections, not the ligand tables.

### Top-ranked transcription factors

| Rank | TF | Role | NR | Other genes it controls | Cells need it? | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SOX2 | activator |  | ~1200 | a few lines | 2 | 0.913 | +134 | 2 | 0 | 79.8 |
| 2 | POU2F1 | activator |  | not counted | some lines | 1 | 0.883 | +133 | 2 | 0 | 75.6 |
| 3 | FOS | activator |  | ~5000 | a few lines | 3 | 1.0 | -4811 | 1 | 0 | 75.0 |
| 4 | JUNB | activator |  | ~5000 | some lines | 5 | 1.0 | -4811 | 5 | 0 | 75.0 |
| 5 | JUND | activator |  | not counted | a few lines | 3 | 1.0 | -4811 | 6 | 0 | 75.0 |
| 6 | FOSL2 | activator |  | not counted | some lines | 4 | 1.0 | -4812 | 4 | 0 | 75.0 |
| 7 | JUN | dual |  | ~5000 | some lines | 3 | 1.0 | -4811 | 8 | 0 | 69.0 |
| 8 | KLF5 | activator |  | ~5000 | some lines | 1 | 1.0 | -1675 | 1 | 0 | 68.5 |
| 9 | STAT3 | activator |  | ~1200 | a few lines | 1 | 1.0 | -1725 | 5 | 0 | 68.5 |
| 10 | NFYB | activator |  | ~5000 | yes, most | 1 | 1.0 | -2093 | 2 | 0 | 68.4 |
| 11 | FOSB | activator |  | 189 | a few lines | 1 | 1.0 | -4811 | 3 | 0 | 68.3 |
| 12 | CEBPA | activator |  | not counted | a few lines | 1 | 0.947 | -2405 | 2 | 0 | 66.8 |
| 13 | E2F4 | activator |  | ~5000 | a few lines | 1 | 0.915 | -1359 | 2 | 0 | 66.3 |
| 14 | POU5F1 | repressor |  | ~1200 | a few lines | 1 | 0.913 | +134 | 6 | 0 | 65.9 |
| 15 | HOXB13 | repressor |  | not counted | a few lines | 1 | 0.96 | +193 | 1 | 0 | 65.9 |
| 16 | NFIB | activator |  | not counted | a few lines | 1 | 0.84 | -935 | 1 | 0 | 65.0 |
| 17 | ETV2 | activator |  | 1979 | a few lines | 3 | 0.96 | +193 | 0 | 0 | 63.1 |
| 18 | BACH1 | dual |  | not counted | a few lines | 1 | 1.0 | -4811 | 1 | 0 | 62.3 |
| 19 | SOX6 | repressor |  | not counted | a few lines | 1 | 0.976 | +456 | 1 | 0 | 61.9 |
| 20 | FOSL1 | unknown |  | not counted | some lines | 3 | 1.0 | -4811 | 1 | 0 | 60.0 |
| 21 | GLIS3 | dual |  | not counted | no | 1 | 0.864 | -2293 | 1 | 0 | 58.3 |
| 22 | CTCF | repressor |  | ~5000 | yes, nearly all | 2 | 0.803 | -1506 | 4 | 0 | 55.6 |
| 23 | KLF9 | unknown |  | not counted | a few lines | 1 | 0.938 | -716 | 1 | 0 | 54.0 |
| 24 | ERF | repressor |  | not counted | no | 1 | 0.852 | -2465 | 1 | 0 | 53.4 |
| 25 | STAT1 | unknown |  | ~1200 | no | 1 | 1.0 | -4754 | 4 | 0 | 53.3 |

**Other genes it controls** is the collateral question. A drug acting through a factor moves everything that factor regulates, so this column says how wide the effect will be. Read it against the score: a strong route through a factor that controls two thousand genes is a worse starting point than a weaker route through one that controls thirty.

A count written `~N` is a bundled estimate for a factor the regulon file does not list (11 here), not a count of curated edges. Treat it as an order of magnitude.

`not counted` means the regulon table has no entry and there is no bundled estimate either (46 here). It is not a small number. Most of the roughly 1,600 human transcription factors have no curated regulon at all.

Counts come from literature curation, so they measure what has been published as much as what the factor does. SP1 has thousands of recorded targets partly because it regulates many genes and partly because thousands of people assayed it.

### What a cell loses without each factor (DepMap)

The regulon column counts how many genes a factor controls, from literature curation, which measures how much anyone has published as much as what the factor does. This is the same collateral question answered by measurement: CRISPR knockout of each factor across roughly 1,100 cancer cell lines.

It also reaches factors the curated tables cannot. A factor nobody has written a regulon for still has a knockout phenotype here.

80 of 81 factor(s) were scored in DepMap. 1 (CTCF) is needed by nearly every cell line tested, so strong inhibition is likely cytotoxic. 62 (BACH1, BACH2, CEBPA, DRGX, E2F4) show a dependency confined to a minority of lines, which is the informative case: the requirement has a context, and whether your tissue shares it is the question to ask next. 14 (ERF, FIGLA, GLIS1, GLIS3, NHLH1) were dispensable throughout. That bounds the proliferative cost of removing them and says nothing about a differentiated tissue.

| Factor | Cells need it? | Mean effect | Lines dependent | Most dependent lineage |
|---|---|---|---|---|
| CTCF | essential in nearly every line | -1.68 | 1147/1150 (100%) | Biliary Tract (-1.80) |
| NFYB | essential in most lines | -0.70 | 961/1150 (84%) | Thyroid (-0.84) |
| KLF5 | essential in a substantial minority | -0.41 | 307/1150 (27%) | Bowel (-0.99) |
| FOSL1 | essential in a substantial minority | -0.41 | 368/1150 (32%) | Pancreas (-0.64) |
| TFAP4 | essential in a minority of lines | -0.35 | 263/1150 (23%) | Myeloid (-0.54) |
| POU2F1 | essential in a minority of lines | -0.34 | 146/1150 (13%) | Kidney (-0.44) |
| GLI1 | essential in a minority of lines | -0.30 | 85/1150 (7%) | Soft Tissue (-0.43) |
| JUN | essential in a minority of lines | -0.29 | 192/1150 (17%) | Pleura (-0.54) |
| JUNB | essential in a minority of lines | -0.27 | 163/1150 (14%) | Pancreas (-0.43) |
| TCF7L1 | essential in a handful of lines | -0.26 | 36/1150 (3%) | Biliary Tract (-0.32) |
| STAT5A | essential in a handful of lines | -0.25 | 50/1150 (4%) | Peripheral Nervous System (-0.31) |
| MGA | essential in a handful of lines | -0.21 | 44/1150 (4%) | Head and Neck (-0.29) |
| MTF1 | essential in a minority of lines | -0.21 | 110/1150 (10%) | Myeloid (-0.43) |
| STAT3 | essential in a handful of lines | -0.19 | 39/1150 (3%) | Cervix (-0.31) |
| STAT5B | essential in a handful of lines | -0.19 | 27/1150 (2%) | Myeloid (-0.49) |
| FOSL2 | essential in a minority of lines | -0.17 | 60/1150 (5%) | Ovary/Fallopian Tube (-0.29) |
| TBX3 | essential in a handful of lines | -0.16 | 50/1150 (4%) | Prostate (-0.53) |
| JUND | essential in a handful of lines | -0.16 | 28/1150 (2%) | Thyroid (-0.25) |
| PATZ1 | essential in a handful of lines | -0.15 | 9/1150 (1%) | Peripheral Nervous System (-0.29) |
| E2F4 | essential in a handful of lines | -0.13 | 14/1150 (1%) | Lymphoid (-0.20) |
| STAT6 | essential in a handful of lines | -0.12 | 3/1150 (0%) | Eye (-0.20) |
| POU5F1 | essential in a handful of lines | -0.11 | 17/1150 (2%) | Biliary Tract (-0.17) |
| ZNF281 | essential in a handful of lines | -0.11 | 29/1150 (2%) | Myeloid (-0.26) |
| FIGLA | dispensable across the panel | -0.10 | 0/1150 (0%) | — |
| ZNF684 | essential in a handful of lines | -0.10 | 8/1150 (1%) | Head and Neck (-0.19) |

**Mean effect** is on DepMap's own scale: 0 means knockout changed nothing, and -1 is the median of known common-essential genes. A line counts as dependent below -0.5.

#### What this does not mean

- **These are cancer cell lines dividing in a dish.** Essentiality here is a statement about proliferating transformed cells. A factor whose job is in a post-mitotic or differentiated cell looks dispensable for reasons that have nothing to do with its importance in your tissue.
- **Knockout is not inhibition.** The protein is removed entirely, including scaffolding roles no ligand would touch. Essential genes are drugged successfully at partial occupancy all the time, so a dependency is a reason for care rather than a veto.
- **Nothing here is about your gene.** This bounds the cost of acting on a factor. Whether acting on it moves your gene is the rest of the report.
- **A dispensable result is not permission.** It bounds the proliferative cost and nothing else.

**The selective ones are the interesting ones.** BACH1, BACH2, CEBPA, DRGX, E2F4, EBF2 matter in a minority of lines rather than all of them, so the requirement has a context. The lineage column is where to look for what that context is, and whether your tissue shares it is the question this panel cannot answer.

### Drugs that behave like losing the factor

Every other route in this report starts from something somebody wrote down - a curated ligand, a motif, a structure, a signed edge. None of them can surface a drug nobody has connected to the factor. This one reads a measured phenotype instead: if a compound hurts the same cell lines that knocking the factor out hurts, it behaves as though it removes the factor.

40740 compound(s) were correlated across 6 factor(s) and none cleared |r| >= 0.2 at q <= 0.1. That is the common outcome. It bounds what this method contributes here rather than saying no such drug exists: the panel is proliferating cancer lines, and a factor that matters elsewhere leaves no trace in it.

### Every route that acts on this gene

235 routes across 7 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

**Why established binders sit at tier 2, just under the direct agents.** A drug already measured to bind one of your factors has half its chain settled: the drug-to-receptor step has assay data behind it, so the pharmacology is de-risked and only the receptor-to-gene step is a hypothesis. A curated tier-3 ligand is annotation on *both* halves. But they stay below tier 1, because a direct agent acts on your gene itself and skips the receptor-to-gene inference entirely — a stronger claim however well-characterised the binder is. The cost of the indirect route is in the row: everything else that factor regulates moves with your gene.

- 23 × curated — mechanism annotated by a human curator
- 207 × predicted — a model says so
- 5 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 3 | TF: known regulator | Cholecalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Paricalcitol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Doxercalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Ergocalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Calcifediol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Calcitriol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: abundance change | vorinostat | FOSL2 | 1 | raises FOSL2 transcript in 88 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | FOS | 1 | raises FOS transcript in 29 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | JUNB | 1 | raises JUNB transcript in 28 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | FOSL2 | 1 | raises FOSL2 transcript in 21 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | JUND | 1 | raises JUND transcript in 14 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | vorinostat | JUND | 1 | raises JUND transcript in 13 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | KLF5 | 1 | raises KLF5 transcript in 9 profile(s) | 0.9 | 0.1 | 0.5 | 41.6 |
| 3 | TF: abundance change | vemurafenib | SOX2 | 1 | raises SOX2 transcript in 5 profile(s) | 0.5 | 0.45 | 0.5 | 38.6 |
| 3 | chromatin-level | Vorinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Romidepsin | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Belinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Panobinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Valproic acid | chromatin (genome-wide) | 3 | weak HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Azacitidine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Decitabine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Tazemetostat | chromatin (genome-wide) | 3 | EZH2 inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Givinostat Hydrochloride | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 4 | upstream signalling | Testosterone Propionate | AR | 2 | activate AR → JUN (1 hop) | 0.75 | 0.513 | 0.6 | 42.0 |
| 4 | upstream signalling | Testosterone Enanthate | AR | 2 | activate AR → JUN (1 hop) | 0.75 | 0.513 | 0.6 | 42.0 |
| 4 | upstream signalling | Oxandrolone | AR | 2 | activate AR → JUN (1 hop) | 0.75 | 0.513 | 0.6 | 42.0 |
| 4 | upstream signalling | Nandrolone Phenpropionate | AR | 2 | activate AR → JUN (1 hop) | 0.75 | 0.513 | 0.6 | 42.0 |
| 4 | upstream signalling | Oxymetholone | AR | 2 | activate AR → JUN (1 hop) | 0.75 | 0.513 | 0.6 | 42.0 |
| 4 | TF: known regulator | Bexarotene | RXRA | 1 | activate RXRA, a recorded regulator of the ge | 0.5 | 0.5 | 0.6 | 42.0 |
| 4 | TF: known regulator | Acitretin | RXRA | 1 | activate RXRA, a recorded regulator of the ge | 0.5 | 0.5 | 0.6 | 42.0 |

Best in each category:

- **TF: known regulator** — Cholecalciferol (tier 3, score 45.2)
- **TF: abundance change** — vorinostat (tier 3, score 44.8)
- **chromatin-level** — Vorinostat (tier 3, score 23.4)
- **upstream signalling** — Testosterone Propionate (tier 4, score 42.0)
- **literature: reported effect** — Acetaminophen (tier 4, score 33.3)
- **direct: predicted product binder** — PARICALCITOL (tier 5, score 62.3)
- **promoter DNA binder** — Plicamycin (tier 5, score 31.2)

*No route found in: direct: gene or transcript, direct: gene product, TF: orthosteric ligand, TF: inferred activity change, TF: predicted pocket binder, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, TF: degradation blocker, off-pathway pocket binder, measured phenocopy, empirical signature, virtual screen hypothesis. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

| Category | What it changes | Your drugs in it | Assay that sees it |
|---|---|---|---|
| TF: known regulator | **abundance** | Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol … +30 | RT-qPCR for your gene after treatment; check the factor's own activity (a reporter or a known target) moved first |
| TF: abundance change | **abundance** | vorinostat, mitoxantrone, vemurafenib | RT-qPCR and a western for the FACTOR first, then RT-qPCR for your gene |
| chromatin-level | **abundance, genome-wide** | Vorinostat, Romidepsin, Belinostat, Panobinostat … +5 | RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened |
| upstream signalling | **abundance** | Testosterone Propionate, Testosterone Enanthate, Oxandrolone, Nandrolone Phenpropionate … +59 | RT-qPCR then western/ELISA; confirm the pathway node moved too, or you cannot attribute the effect |
| literature: reported effect | **abundance (reported)** | Acetaminophen, Calcifediol, Calcitriol, Cholecalciferol … +14 | RT-qPCR in your cells; read the paper for the dose, cell type and species first |
| direct: predicted product binder | **activity** | PARICALCITOL | a binding assay against the protein first, then a functional assay for direction - RT-qPCR will not see this route at all, and a western blot reads it as no effect |
| promoter DNA binder | **abundance** | Plicamycin, Mitoxantrone, Doxorubicin, Pentamidine | an EMSA or footprint on this promoter first, then a reporter, then RT-qPCR - and a counter-screen on an unrelated promoter carrying the same feature, because selectivity is the whole question here |

*Every category present here works through transcription, so all of them move your gene's product **abundance** and an RT-qPCR is the right first readout. The one exception in the taxonomy — a drug on the gene product itself, which changes activity only — did not appear in this run.*

**TF: known regulator** (36 found)

- *Yours:* Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol, Calcifediol, Calcitriol, Bexarotene, Acitretin … +26 more
- Drugs that act on a factor CollecTRI or perturbation data records as raising or lowering your gene, whether or not the motif scan ranked that factor.
- *Why:* Reaches factors the motif scan misses, such as those that bind through a partner or a distant enhancer.
- *Cost:* The link comes from a database or a knockout in other cells, so it may not hold in yours. The drug's effect on the factor is taken from its recorded action.
- *Differs:* Starts from the literature on your gene, not from its promoter sequence. CollecTRI and CTD share sources, so agreement between them is not independent.

**TF: abundance change** (8 found)

- *Yours:* vorinostat, mitoxantrone, vemurafenib
- A compound measured to change how much of the factor's own transcript exists, rather than what the factor does.
- *Why:* The only route here that needs the factor to be expressed rather than druggable, so it reaches factors with no ligand, no structure and no curated upstream kinase - which is most of them.
- *Cost:* More factor is not more factor ACTIVITY. For a ligand-activated receptor the two are nearly unrelated, and transcript is not protein. Broad agents that shut down transcription move every factor and are flagged.
- *Differs:* Unlike 'empirical signature', which is a compound measured against YOUR GENE, this was measured against the FACTOR - one further inferential step from what you care about.

**chromatin-level** (9 found)

- *Yours:* Vorinostat, Romidepsin, Belinostat, Panobinostat, Valproic acid, Azacitidine, Decitabine, Tazemetostat … +1 more
- HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.
- *Why:* Reliable positive control that your reporter *can* be induced.
- *Cost:* Specificity ~0.05. These move thousands of promoters. Almost never an acceptable therapeutic route on that basis alone.
- *Differs:* Not targeted at your gene in any sense; it is in the table so you can see what a non-specific option scores.

**upstream signalling** (159 found)

- *Yours:* Testosterone Propionate, Testosterone Enanthate, Oxandrolone, Nandrolone Phenpropionate, Oxymetholone, Raloxifene Hydrochloride, Diethylstilbestrol Diphosphate, Estramustine Phosphate Sodium … +55 more
- Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- *Why:* Turns 'this factor has no drugs' into a real list. Most TFs are only reachable this way.
- *Cost:* Specificity falls with every step: a node three hops up reaches your gene through one branch and thousands of other genes through the rest.
- *Differs:* The only category where the drug never touches your gene or its TF — read the printed path and judge it yourself.

**literature: reported effect** (18 found)

- *Yours:* Acetaminophen, Calcifediol, Calcitriol, Cholecalciferol, Dexamethasone, Ketoconazole, Leflunomide, Rifampin … +10 more
- A chemical CTD records as raising or lowering this gene's mRNA, with no mechanism stated.
- *Why:* It is the only route that can name a drug nothing here can explain, which includes most of what CTD holds: stress, feedback and toxicity effects.
- *Cost:* Most records rest on one paper, often at a toxicology dose in one cell type. Across 36 genes only about 3% of them were reproduced in LINCS, and more papers did not raise that. Famous drugs are over-represented.
- *Differs:* Every other category predicts an effect from a mechanism. This reports an effect somebody saw, and cannot suggest a drug nobody has tested on this gene.

**direct: predicted product binder** (1 found)

- *Yours:* PARICALCITOL
- A compound predicted to bind this gene's own protein and not recorded as acting on it - either because its chemistry resembles the protein's known ligands, or because it docks into the protein's pocket.
- *Why:* It is the only route here aimed at the gene's own product rather than at a factor upstream of it, so it skips every inference about whether moving a factor moves your gene. And unlike the curated product drugs it can propose a compound no one has annotated.
- *Cost:* Predicted binding and nothing more. The model is refused below a cross-validated AUC floor, which bounds how badly it can be ranking, but a hit is a hypothesis for a binding assay. Direction is separate again: whether binding inhibits or activates the protein decides whether it serves your goal.
- *Differs:* Unlike 'direct: gene product', nothing was measured or curated for this compound at this protein - that category is annotation, this is a model. Unlike 'TF: predicted pocket binder', the pocket is on the gene's own product, so no inference about a factor controlling the gene is involved.

**promoter DNA binder** (4 found)

- *Yours:* Plicamycin, Mitoxantrone, Doxorubicin, Pentamidine
- A small molecule with a documented preference for a sequence feature this promoter carries - a quadruplex, a GC-box, an AT tract - proposed because the feature is enriched against shuffles of this same sequence.
- *Why:* It is the only route that acts on the promoter itself rather than on a protein, so it needs no druggable factor and no annotation. Where a promoter depends on one feature, occupying it is a direct intervention.
- *Cost:* Sequence preference is not selectivity. These agents bind the feature wherever it occurs and most of the genome's copies are not in your gene, so the dose that engages your promoter engages thousands of others. Nearly every approved member is cytotoxic, which is how it works rather than a side effect. And a motif is not a folded structure in a cell.
- *Differs:* Unlike every TF route, no protein is targeted - which removes the question of whether the factor controls your gene and replaces it with whether the molecule can tell your promoter from any other.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **direct: gene product** — Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- **TF: orthosteric ligand** — Binds the transcription factor's natural ligand pocket.
- **TF: inferred activity change** — A compound whose own LINCS signature moves this factor's known targets, found without any record that it binds the factor.
- **TF: predicted pocket binder** — A compound predicted by docking to fit the factor's ligand pocket, from a structure rather than from any measured activity.
- **TF: allosteric modulator** — Binds a site away from the ligand pocket, amplifying or damping the endogenous ligand's effect (PAM/NAM).
- **TF: stabiliser or glue** — Locks the factor in an active conformation, or glues it to a partner protein.
- **TF: endogenous ligand** — Changes how much natural ligand exists, by drugging the enzyme that makes or destroys it.
- **TF: degradation blocker** — Stops the factor being destroyed, so more of it accumulates.
- **off-pathway pocket binder** — A docking run against a protein that is not one of the factors found at this promoter. Usually an upstream node someone chose to dock by hand.
- **measured phenocopy** — A compound whose effect on cell viability across the DepMap panel tracks what happens when this factor is knocked out, so it behaves as though it removes the factor.
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **virtual screen hypothesis** — An approved drug predicted to bind one of your TFs, never tested against it.

### What these drugs are

Every route in this report is off-label by construction: the drug was approved for one thing and is being considered for another. The approved use is the baseline that departure is measured from, so it is listed here once per drug rather than repeated in each table above.

| Drug | What it is | Approved to treat |
|---|---|---|
| Acitretin | retinoid | severe psoriasis |
| Azacitidine | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia |
| Beclomethasone | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Bexarotene | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma |
| Calcifediol | 25-hydroxyvitamin D | vitamin D deficiency |
| Calcitriol | active vitamin D, VDR agonist | hypocalcaemia and secondary hyperparathyroidism in renal failure |
| Cholecalciferol | vitamin D3, a precursor requiring activation | vitamin D deficiency and supplementation |
| Clofibrate | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) |
| Decitabine | DNA methyltransferase inhibitor | myelodysplastic syndrome |
| Dexamethasone | potent glucocorticoid receptor agonist | inflammation, cerebral oedema, myeloma, and as an antiemetic |
| Dextrothyroxine Sodium | T4 stereoisomer | formerly hyperlipidaemia - withdrawn for cardiac mortality |
| Diethylstilbestrol Diphosphate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Doxercalciferol | vitamin D prodrug | secondary hyperparathyroidism in dialysis |
| Doxorubicin | anthracycline topoisomerase II inhibitor | many solid and haematological cancers |
| Elafibranor | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) |
| Ergocalciferol | vitamin D receptor agonist class | secondary hyperparathyroidism in renal disease, osteoporosis, and plaque psoriasis topically |
| Estradiol Acetate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Estradiol Cypionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Etretinate | retinoid | severe psoriasis - WITHDRAWN, teratogenic |
| Fenofibrate | PPARalpha agonist | hypertriglyceridaemia |
| Gemfibrozil | PPARalpha agonist | hypertriglyceridaemia |
| Hydrocortisone | glucocorticoid and mineralocorticoid | adrenal insufficiency; topical dermatitis |
| Hydrocortisone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Isotretinoin | retinoid | severe nodular acne |
| Leflunomide | broad immunosuppressant class | transplant rejection and autoimmune disease |
| Lenalidomide Hydrochloride Monohydrate | cereblon modulator | multiple myeloma, MDS |
| Levothyroxine Sodium | synthetic T4 | hypothyroidism |
| Liothyronine Sodium | synthetic T3 | hypothyroidism, myxoedema coma |
| Mometasone Furoate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Nandrolone Phenpropionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Oxandrolone | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Palovarotene | retinoid class | acne, psoriasis and other keratinisation disorders topically; acute promyelocytic leukaemia and cutaneous lymphoma systemically |
| Paricalcitol | selective VDR activator | secondary hyperparathyroidism in chronic kidney disease |
| Plicamycin | RNA synthesis inhibitor, binds GC-rich DNA | formerly hypercalcaemia of malignancy - withdrawn |
| Prednisolone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Raloxifene Hydrochloride | selective oestrogen receptor modulator | osteoporosis; breast cancer risk reduction |
| Romidepsin | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma |
| Testosterone | androgen receptor agonist | male hypogonadism |
| Testosterone Enanthate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Testosterone Propionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Tretinoin | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically |
| mitoxantrone | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis |
| vorinostat | histone deacetylase inhibitor | cutaneous T-cell lymphoma |

**No description established for 56**: Acetaminophen, Adenosine, Almotriptan Malate, Avatrombopag, Bazedoxifene Acetate, Belinostat, Bethanechol Chloride, Carboprost Tromethamine, Cevimeline Hydrochloride, Chenodiol, Chlorotrianisene, Choline Fenofibrate, Cinacalcet Hydrochloride, Copper Sulfate, Dinoprost Tromethamine, Docetaxel Anhydrous, Eletriptan Hydrobromide, Eltrombopag Olamine, Estramustine Phosphate Sodium, Etrasimod, Fenofibric Acid, Fingolimod Hydrochloride, Fingolimod Lauryl Sulfate, Givinostat Hydrochloride, Iloprost…. That means this package could not find one offline, not that the compound has no approved use. Several are environmental agents and research chemicals from the CTD tables, which have no indication because they are not drugs.

Descriptions come from: 16 class, 27 curated, 56 unknown. None of these is a label. The curated text is a summary written for orientation, ChEMBL's indication list includes uses studied but never approved, and a class line gives the kind of drug rather than its indication. Check the prescribing information before acting on any row.

### Cell-cycle compatibility

Tissue proliferation state: **unknown**

1 factor(s) here are cell-cycle restricted (E2F4). Set --proliferation to have this scored rather than just flagged.

| TF | Phase | Duty cycle | Effect multiplier |
|---|---|---|---|
| E2F4 | G1 | 46% | 1.0 |

*A partial duty cycle dilutes the population-level effect in asynchronous cells: per-cell efficacy can be much higher than a bulk RT-qPCR suggests. Synchronise or use a phase-resolved readout before calling it a failure.*

### Upstream routes to factors with no direct drug

Most TFs have no approved ligand. These paths walk upstream through signed causal edges — kinases, phosphatases, receptors — to proteins that *are* druggable. Signs multiply along the path, so an inhibitor of an inhibitor counts as activating.

Distance buys druggability and spends selectivity: a node three hops up reaches your gene through one branch and everything else through the others. Read the fan-out column.

**SOX2** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| RARG | 2 | 5 | RARG --> AKT1 --> SOX2 | activate | Acitretin, Palovarotene, Tretinoin | 31.5 |
| RARB | 2 | 6 | RARB --> AKT1 --> SOX2 | activate | Acitretin, Tretinoin | 30.9 |
| PTGER3 | 2 | 9 | PTGER3 --> EGFR --> SOX2 | activate | Misoprostol | 29.0 |
| RARA | 2 | 15 | RARA --> AKT1 --> SOX2 | activate | Acitretin, Tretinoin | 26.0 |
| DDB1 | 3 | 1 | DDB1 --> E2F1 --> BRCA1 --> SOX2 | activate | Lenalidomide Hydrochloride Monohydrate | 18.7 |
| ESR1 | 2 | 37 | ESR1 --> AKT1 --> SOX2 | activate | Estradiol Acetate, Diethylstilbestrol Diphosphate, Toremifene Citrate | 18.6 |
| MPL | 3 | 2 | MPL --> JAK2 --> EGFR --> SOX2 | activate | Avatrombopag, Avatrombopag Maleate, Lusutrombopag | 18.3 |
| CASR | 3 | 2 | CASR --> FYN --> EGFR --> SOX2 | activate | Cinacalcet Hydrochloride | 18.3 |
| AR | 2 | 38 | AR --> EGFR --> SOX2 | activate | Testosterone Propionate, Testosterone Enanthate, Oxandrolone | 18.2 |
| ESR2 | 3 | 4 | ESR2 --> RARA --> AKT1 --> SOX2 | activate | Raloxifene Hydrochloride, Diethylstilbestrol Diphosphate, Estramustine Phosphate Sodium | 17.5 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 4 downstream branches

**POU2F1** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB —\| PRKAA1 —\| POU2F1 | activate | Acetaminophen | 18.8 |
| RARG | 3 | 5 | RARG --> AKT1 —\| PRKAA1 —\| POU2F1 | activate | Acitretin, Palovarotene, Tretinoin | 17.1 |
| RARB | 3 | 6 | RARB --> AKT1 —\| PRKAA1 —\| POU2F1 | activate | Acitretin, Tretinoin | 16.7 |
| S1PR1 | 3 | 12 | S1PR1 --> LYN —\| PRKDC —\| POU2F1 | activate | Ponesimod, Siponimod, Etrasimod | 14.8 |
| RARA | 3 | 15 | RARA --> AKT1 —\| PRKAA1 —\| POU2F1 | activate | Acitretin, Tretinoin | 14.0 |
| ESR1 | 3 | 37 | ESR1 --> AKT1 —\| PRKAA1 —\| POU2F1 | activate | Estradiol Acetate, Diethylstilbestrol Diphosphate, Toremifene Citrate | 10.0 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 12 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 15 downstream branches

**FOS** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| MTNR1A | 2 | 7 | MTNR1A --> MAPK1 --> FOS | activate | Ramelteon, Tasimelteon | 30.2 |
| NOD2 | 2 | 9 | NOD2 --> MAPK14 --> FOS | activate | — | 29.2 |
| PTGER3 | 2 | 9 | PTGER3 --> SRC --> FOS | activate | Misoprostol | 28.9 |
| S1PR1 | 2 | 12 | S1PR1 --> SRC --> FOS | activate | Ponesimod, Siponimod, Etrasimod | 27.3 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK1 --> FOS | activate | Acetaminophen | 18.7 |
| ESR1 | 2 | 37 | ESR1 --> SRC --> FOS | activate | Estradiol Acetate, Diethylstilbestrol Diphosphate, Toremifene Citrate | 18.4 |
| CASR | 3 | 2 | CASR --> FYN --> MAPK14 --> FOS | activate | Cinacalcet Hydrochloride | 18.3 |
| S1PR5 | 3 | 2 | S1PR5 --> GNAI1 --> SRC --> FOS | activate | Ozanimod, Etrasimod Arginine, Fingolimod Lauryl Sulfate | 18.3 |
| ADORA3 | 3 | 2 | ADORA3 --> GNAI1 --> SRC --> FOS | activate | Adenosine | 18.3 |
| HTR1F | 3 | 2 | HTR1F --> GNAI1 --> SRC --> FOS | activate | Almotriptan Malate, Eletriptan Hydrobromide, Lasmiditan | 18.3 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ path passes through hub(s): MAPK1

**JUNB** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| NOD2 | 2 | 9 | NOD2 --> MAPK8 --> JUNB | activate | — | 30.9 |
| PTGIR | 2 | 7 | PTGIR --> STAT3 --> JUNB | activate | Treprostinil, Treprostinil Diolamine, Selexipag | 30.3 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK8 --> JUNB | activate | Acetaminophen | 18.7 |
| MPL | 3 | 2 | MPL --> JAK2 --> STAT3 --> JUNB | activate | Avatrombopag, Avatrombopag Maleate, Lusutrombopag | 18.3 |
| CASR | 3 | 2 | CASR --> FYN —\| ITCH —\| JUNB | activate | Cinacalcet Hydrochloride | 18.3 |
| GHR | 3 | 4 | GHR --> JAK2 --> STAT3 --> JUNB | activate | — | 17.5 |
| CXCR4 | 3 | 5 | CXCR4 --> JAK2 --> STAT3 --> JUNB | activate | Plerixafor | 17.1 |
| RARG | 3 | 5 | RARG --> AKT1 --> STAT3 --> JUNB | activate | Acitretin, Palovarotene, Tretinoin | 17.1 |
| CHRM3 | 3 | 5 | CHRM3 --> GNAQ —\| GSK3B —\| JUNB | activate | Pilocarpine Hydrochloride, Bethanechol Chloride | 17.1 |
| HTR2B | 3 | 5 | HTR2B --> GNAQ —\| GSK3B —\| JUNB | activate | Methylergonovine Maleate | 17.1 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 4 downstream branches

**JUND** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| MTNR1A | 2 | 7 | MTNR1A --> MAPK1 --> JUND | activate | Ramelteon, Tasimelteon | 30.4 |
| NOD2 | 2 | 9 | NOD2 --> MAPK8 --> JUND | activate | — | 30.3 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK8 --> JUND | activate | Acetaminophen | 18.8 |
| RARG | 3 | 5 | RARG --> AKT1 --> BLVRA --> JUND | activate | Acitretin, Palovarotene, Tretinoin | 17.1 |
| RARB | 3 | 6 | RARB --> AKT1 --> BLVRA --> JUND | activate | Acitretin, Tretinoin | 16.8 |
| PTGIR | 3 | 7 | PTGIR --> STAT3 --> MAPK8 --> JUND | activate | Treprostinil, Treprostinil Diolamine, Selexipag | 16.3 |
| PTGER3 | 3 | 9 | PTGER3 --> SRC --> MAPK1 --> JUND | activate | Misoprostol | 15.7 |
| S1PR1 | 3 | 12 | S1PR1 --> LYN --> MAPK3 --> JUND | activate | Ponesimod, Siponimod, Etrasimod | 14.9 |
| RARA | 3 | 15 | RARA --> AKT1 --> BLVRA --> JUND | activate | Acitretin, Tretinoin | 14.0 |
| NR3C1 | 3 | 16 | NR3C1 --> LCK --> MAPK3 --> JUND | activate | Hydrocortisone Acetate, Prednisolone Acetate, Mometasone Furoate | 14.0 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 12 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 15 downstream branches

**JUN** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| AR | 1 | 38 | AR --> JUN | activate | Testosterone Propionate, Testosterone Enanthate, Oxandrolone | 44.3 |
| ESR2 | 2 | 4 | ESR2 --> MAPK11 --> JUN | activate | Raloxifene Hydrochloride, Diethylstilbestrol Diphosphate, Estramustine Phosphate Sodium | 32.6 |
| MTNR1A | 2 | 7 | MTNR1A --> MAPK1 --> JUN | activate | Ramelteon, Tasimelteon | 30.0 |
| NOD2 | 2 | 9 | NOD2 --> MAPK8 --> JUN | activate | — | 29.2 |
| NR3C1 | 2 | 16 | NR3C1 --> SMAD3 --> JUN | activate | Hydrocortisone Acetate, Prednisolone Acetate, Mometasone Furoate | 25.7 |
| RXRA | 2 | 22 | RXRA --> AR --> JUN | activate | Bexarotene, Acitretin, Etretinate | 22.7 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK8 --> JUN | activate | Acetaminophen | 18.7 |
| NR1H4 | 3 | 1 | NR1H4 --> RXRA --> AR --> JUN | activate | Ursodiol, Chenodiol, Obeticholic Acid | 18.7 |
| ESR1 | 2 | 37 | ESR1 --> AR --> JUN | activate | Estradiol Acetate, Diethylstilbestrol Diphosphate, Toremifene Citrate | 18.5 |
| CASR | 3 | 2 | CASR --> FYN --> MAPK14 --> JUN | activate | Cinacalcet Hydrochloride | 18.3 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ path passes through hub(s): MAPK1

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

#### No stabilisation route for 8 factor(s)

FOS, FOSL2, JUN, JUNB, JUND, KLF5, POU2F1, SOX2

These have no ligand pocket, no direct drug, and **nothing curated against them or their degradation machinery**. For these factors stabilisation is a research programme, not a route you can take: making it real needs the E3 ligase that degrades the factor *and* a ligand for that ligase, and neither is published here.

The mechanism class is real — it has shipped approved drugs, listed once below — but those act on **other proteins entirely**, and none of them is evidence that any factor above is tractable.

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Cell-type reach

Three things must overlap: where the drug goes, where the gene is expressed, and where the transcription factor is. Any empty intersection kills the route regardless of motif quality or potency.

| TF | Co-expressed with gene | Shared tissues |
|---|---|---|
| SOX2 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| POU2F1 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| FOS | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| JUNB | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| JUND | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| FOSL2 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| JUN | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| KLF5 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| STAT3 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| NFYB | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |

### Human genetic evidence

Human genetic variation in CYP24A1 associates with: hypercalcemia, infantile, nephrolithiasis, atopic eczema. Read this as a preview of what modulating the gene may do, including the adverse effects — this is the reasoning that anticipated PCSK9 inhibition years before the trials.

| Phenotype | Genetic score | Overall |
|---|---|---|
| hypercalcemia, infantile | 0.902 | 0.814 |
| nephrolithiasis | 0.92 | 0.568 |
| atopic eczema | 0.913 | 0.565 |
| urolithiasis | 0.699 | 0.531 |
| ureterolithiasis | 0.846 | 0.515 |
| bladder calculus | 0.806 | 0.49 |
| hereditary disease | 0.773 | 0.472 |
| aortic stenosis | 0.677 | 0.411 |
| vitamin D deficiency | 0.593 | 0.367 |
| deep vein thrombosis | 0.59 | 0.359 |
| COVID-19 | 0.576 | 0.355 |
| allergic disease | 0.578 | 0.352 |
| Abnormal thrombosis | 0.578 | 0.352 |

#### Does variation in the factor actually move this gene? (trans-eQTL)

The strongest evidence available anywhere in this report. A motif says a factor *could* bind; a trans-eQTL says natural genetic variation near it **does** change your gene's expression, in real human tissue, with exposure randomised by meiosis.

**No factor can be supported this way for CYP24A1.** CYP24A1 has no significant trans-eQTL in eQTLGen at all, so the test returns nothing for every factor here. That is a property of the dataset's power and of the tissue (whole blood), not evidence that CYP24A1 is unregulated.

*No factor tested (SOX2, POU2F1, FOS, JUNB, JUND, FOSL2) has trans-eQTL support for this gene.*

*eQTLGen is whole blood at n=31,684, FDR<0.05. A hit is strong causal evidence. **Absence is weak** — trans effects are small and detecting them needs enormous power, so a null overwhelmingly means underpowered rather than unregulated, and a relationship operating outside blood would not appear at all.*

### Routes a promoter analysis cannot see

Everything above acts on transcription. These act after it, around it, or on a different gene entirely. For some genes that is where the whole problem is: if the transcript is destroyed, the allele is silenced, or the protein is degraded too fast, more transcription does not help and may not be possible.

#### Could a relative do the job instead?

Whether a paralog could substitute for this gene is unanswered. It is worth asking explicitly, because the question 'restore this function' has more answers than 'raise this transcript'.

- CYP24A1 is not in the curated redundancy table, which is short and certainly incomplete. Absence here is not evidence that no paralog substitutes — it means nobody has written one down in this file.

#### After transcription

A promoter-level analysis stops at transcription. These act after it, and for some genes that is where the whole problem is: if the transcript is destroyed or mistranslated, doubling transcription doubles the amount that is destroyed. 1 of 6 mechanism(s) here have an FDA-approved agent. The rest are listed so that an empty transcriptional result is not mistaken for 'nothing can be done' — they are mechanisms, not treatments.

| Mechanism | Raises | Status | Applies when |
|---|---|---|---|
| splice correction | functional protein | approved | a splice-site or splice-regulatory variant, a poison exon, |
| transcript stabilisation | transcript | investigational | the transcript is short-lived and the promoter is already  |
| microRNA de-repression | transcript | investigational | a specific miRNA is a documented repressor of your gene in |
| translational enhancement | protein | investigational | protein is low while transcript is normal — the signature  |
| nonsense readthrough | functional protein | approved elsewhere | the gene carries a premature termination codon (PTC) |
| nonsense-mediated decay (NMD) inhibition | transcript | investigational | the gene carries a PTC and the transcript is being degrade |

**Splice correction** — redirects splicing so a productive transcript is made — exon inclusion, exon skipping, or suppression of a cryptic site.

- Exon-skipping restores a partly functional, internally deleted protein, not a wild-type one — the DMD agents are judged on dystrophin that is shorter than normal.
- Accelerated approvals in this class rest on a surrogate endpoint (protein level), and the confirmatory functional benefit has been contested.

#### Raising the protein instead of the transcript

- These leave transcription alone. If your protein is made normally and degraded too fast, every transcriptional route in this report is aimed at the wrong step - and nothing in this package can tell the two apart.
- The experiment that distinguishes them: compare transcript and protein. Normal mRNA with low protein points here; low mRNA points at everything above.

| Mechanism | Specificity | Status | Why |
|---|---|---|---|
| proteasome inhibition | 0.03 | approved | essentially none: the proteasome degrades most of the proteome |
| lysosomal / autophagy blockade | 0.05 | approved | non-specific, and the two degradation routes are distinct - wh |
| pharmacological chaperone | 0.85 | approved | genuinely high: a chaperone binds one protein, and this is the |
| E3 ligase inhibition | 0.7 | clinical | potentially high - an E3 has a limited substrate set - but it  |

*Three of these four raise most of the proteome, which makes them excellent diagnostics and poor therapies. The pharmacological chaperone is the exception — and it needs a chaperone to exist for your protein.*

### Candidates that may beat the current standard

**No incumbent found.** Nothing in the curated tables or CTD currently controls this gene's expression. There is no bar to clear, so candidates below are ranked on their own merits — check the gap is real rather than a hole in the databases.

*Accounting for all 0 candidates: 0 already linked, 0 novel and non-dominated, 0 novel but beaten on every axis by another candidate, 0 with no advantage. All of them are in results.json under `improvements.all_candidates`.*

#### Compounds already reported to change this gene

These come from CTD's curated chemical–gene interactions: someone measured an expression change and a curator recorded it. That makes them the most concrete evidence in this report — but the direction is what decides whether each is a lead or a warning.

180 chemicals total: **85** move it the way you want, 48 move it the opposite way, 47 unclear.

**Expect the 'predicted by this run?' column to say no, and do not read that as a failed match.** CTD indexes MeSH, which is a vocabulary for the toxicology and pharmacology literature rather than a drug list, so what comes back is three different kinds of thing:

- **chemical** (153) — a compound. The only kind that can match a candidate.
- **exposure** (20) — formaldehyde, asbestos, bisphenol, tobacco smoke. Real measured effects on your gene, and not things anyone can be given.
- **category** (7) — MeSH headings such as "Anti-HIV Agents" or "Plant Extracts". Not compounds at all, so the column is blank rather than no.

Across all three directions, **0 of 180 matched a candidate from this run**. A low number here is the normal result: this report's drug tables are approved pharmacology, and most of CTD is not. It means the two sources are looking at different parts of chemical space, not that either is wrong.

**Reported to move CYP24A1 UP** — leads:

| Chemical | Kind | Predicted by this run? | CTD actions |
|---|---|---|---|
| Acetaminophen | chemical | no | affects^expression; increases^expression |
| Beclomethasone | chemical | no | increases^expression |
| Calcifediol | chemical | no | affects^expression; increases^expression |
| Calcitriol | chemical | no | affects^expression; increases^expression |
| Cholecalciferol | chemical | no | affects^expression; increases^expression |
| Copper Sulfate | chemical | no | increases^expression |
| Dexamethasone | chemical | no | affects^expression; increases^expression |
| Hydrocortisone | chemical | no | increases^expression |
| Ketoconazole | chemical | no | affects^expression; increases^expression |
| Leflunomide | chemical | no | increases^expression |
| Paricalcitol | chemical | no | increases^expression |
| Phenobarbital | chemical | no | increases^expression |
| Rifampin | chemical | no | increases^expression |
| Sorafenib | chemical | no | increases^expression |
| Testosterone | chemical | no | affects^expression; increases^expression |
| Tretinoin | chemical | no | affects^expression; increases^expression |
| Valproic Acid | chemical | no | affects^expression; increases^expression |
| Zinc Sulfate | chemical | no | increases^expression |

**22 compounds have a measured effect on this gene and no mechanism from this analysis.** Exposures and MeSH headings are excluded from that count, because the scan was never going to propose formaldehyde and listing it as a miss says nothing.

Most of what remains is research chemistry rather than medicine: tool compounds, probes and natural-product fractions that are measured because somebody had them on a shelf. Read them as evidence the gene IS movable, and check the few that are approved drugs individually.

Acetaminophen, Beclomethasone, Calcifediol, Calcitriol, Cholecalciferol, Copper Sulfate, Cyclosporine, Dexamethasone, Estradiol, Ethinyl Estradiol, Hydrocortisone, Ketoconazole, Leflunomide, Mirdametinib, Paricalcitol, Phenobarbital, Rifampin, Sorafenib, Testosterone, Tretinoin, Valproic Acid, Zinc Sulfate

**Reported to move it the WRONG way** — avoid, or use as negative controls:

Ampicillin, Arsenic Trioxide, Ascorbic Acid, Cisplatin, Copper, Flurbiprofen, Gemcitabine, Niclosamide, Progesterone, Raloxifene Hydrochloride, Tamoxifen, Teriflunomide

*CTD aggregates across cell types, doses and assays. A 'both' or contradictory annotation usually means the effect is context-dependent rather than that the data are wrong — check the underlying study before relying on any single row.*

#### No candidate scored above 0.5 on any axis

There is no incumbent here, so nothing was filtered out for failing to beat one — the candidates themselves are simply weak. Widen the promoter window or check whether this gene is regulated through an enhancer outside it.

### Drug–gene links with no recorded direction

66 interaction(s) across 20 factor(s) come from curated sources that record a drug–gene relationship **without saying which way it goes**. They are held back from the recommendations above, because the entire question here is direction and an undirected link cannot answer it.

They are still leads — someone curated a real relationship — but they need a literature check, not a prediction. Many are chemotherapeutics or broad stressors whose association with a transcription factor reflects cell-wide perturbation rather than a specific mechanism.

| TF | Drug | Reported by | DGIdb score |
|---|---|---|---|
| STAT5A | Tofacitinib | DTC | 2.34 |
| CTCF | Capivasertib | PRISM | 2.04 |
| KLF5 | Secnidazole | PRISM | 1.95 |
| STAT5B | Tofacitinib | DTC | 1.56 |
| ETV1 | Trametinib Dimethyl Sulfoxide | CIViC | 1.51 |
| ETV5 | Trametinib Dimethyl Sulfoxide | CIViC | 1.51 |
| NFIB | Paliperidone | ClinPGx | 1.11 |
| TFAP4 | Azacitidine | PRISM | 1.09 |
| STAT6 | Esomeprazole Sodium | ClinPGx | 0.99 |
| JUND | Clozapine | NCI | 0.94 |
| JUN | Copper Chloride | DTC | 0.9 |
| POU5F1 | Allopurinol | ClinPGx | 0.73 |
| POU5F1 | Phenytoin Sodium | ClinPGx | 0.56 |
| NHLH1 | Cisplatin | ClinPGx | 0.55 |
| JUNB | Brigatinib | PRISM | 0.53 |
| POU5F1 | Carbamazepine | ClinPGx | 0.48 |
| FOS | Bromocriptine Mesylate | NCI | 0.47 |
| E2F4 | Crizotinib | PRISM | 0.46 |
| JUN | Atomoxetine Hydrochloride | DTC | 0.45 |
| NFIB | Clozapine | ClinPGx | 0.39 |
| NFIB | Methylphenidate | ClinPGx | 0.37 |
| STAT5B | Dasatinib Anhydrous | Drug Repurposing Hub; TEND; TdgClinicalTrial | 0.34 |
| FOS | Pilocarpine Hydrochloride | NCI | 0.33 |
| FOS | Baclofen | NCI | 0.31 |
| JUN | Bupropion Hydrochloride | DTC | 0.3 |
| … 41 more | | | |

*Pass `--include-unknown-direction` to fold these into the main tables anyway — but every such row will read UNKNOWN / unclear, which is what they are.*

### Approved drugs acting on these factors

No approved drug with activating pharmacology was found for any TF matched in this promoter. That is the usual result outside the nuclear-receptor family; see the README for indirect strategies.

## Promoter of CYP24A1 (ENST00000216862)

- Transcript chosen: **distinct TSS (representing 2 within 500 bp)**
- Assembly: from GTF  |  Strand: -
- TSS: chr20:54173986
- Window: chr20:54172987-54178986 (-5000 / +1000, 6000 bp)
- Motifs scanned: 879 at p ≤ 1e-05 per position
- **The window was widened.** The requested window (-1000/+200) showed no signal above chance (16 sites), so it was rescanned at -5000/+1000, which gave 107 sites against about 60.6 expected by chance. Sites further out say less about this gene than core promoter sites do. `--no-widen` keeps the requested window.
- 311 of those motifs were not scanned: in this window their best possible match is more common than p = 1e-05, so any hit would be a chance match at a higher rate than stated. A looser `--pvalue` scans them, with more chance hits.
- Non-redundant sites: 107 across 89 TFs

> **⚠ Weak scan.** 107 sites against 61 expected by chance over 4358 scannable bases. Most of this list is noise; treat only the highest-scoring sites with ChIP support as real.

*about 57% of these 107 sites (range 46% to 72%) would be expected on shuffled sequence of this window, so roughly 46 are above chance; among the best-scoring sites the share is lower: 75% of the top 5, 100% of the top 10, 93% of the top 20. Estimated from 30 shuffled copies that keep the repeat mask and every dinucleotide count. The formula that assumes each motif matches at exactly its stated p-value would expect 49 and is not used. Nothing here says which sites are real; only ChIP evidence separates them.*

#### What else is in this window

Every other specificity check here asks whether a drug or factor affects other genes. This asks the prior question: is the sequence that was scanned uniquely this gene's?

The window overlaps 1 other transcription unit(s), including ENSG00000286587 on the opposite strand. Sites called here are not unambiguously your gene's, and for the antisense overlap the direction logic does not carry.

| Gene | Biotype | Location | Strand | Relation |
|---|---|---|---|---|
| ENSG00000286587 | lncRNA | chr20:54059120-54190210 | + | antisense overlap |

- **ENSG00000286587** — window lies inside this gene on the OPPOSITE strand. A site called here may drive that transcript instead, and an antisense RNA can stabilise or repress your gene — so the pipeline's direction logic does not reliably carry

*An antisense overlap matters more than it looks. The direction logic in this report assumes a factor acts on YOUR promoter; where the window lies inside an antisense transcript, a site may drive that instead, and antisense RNA can stabilise or repress the sense gene — so the sign does not carry.*

*This checks annotation, not sequence uniqueness. BLAT the window if you need to know it appears nowhere else in the genome.*

### Ligand-activated (nuclear) receptors with sites here

These are the directly druggable ones.

| TF | Role | Sites | Best rel. score | Best site (vs TSS) | ChIP peaks | Score |
|---|---|---|---|---|---|---|
| HNF4A | activator | 1 | 1.0 | CAAAGTCCA (-3680) | 4 | 78.3 |
| ESRRB | activator | 1 | 0.94 | ACAAGGTCAT (-3640) | 0 | 56.5 |
| HNF4G | activator | 1 | 1.0 | CAAAGTCCA (-3680) | 0 | 48.3 |

#### Why the drug tables below look the way they do

**89 factor(s) had predicted sites here. 2 of them have an approved drug, and all 2 are nuclear receptors.**

Roughly 48 of the ~1,600 human transcription factors have a ligand-binding pocket, and the steroid, retinoid, thyroid and vitamin-D families are most of them. So the drug tables are not a ranking that hormone drugs won — they are what is left after every factor without a ligand has been dropped.

**24 of the top 25 factors by score have no drug and appear in no drug table:** NFYB, KLF5, STAT3, FOS, JUNB, JUND, FOSL2, WT1, CEBPA, JUN, E2F4, ELK1 ….

That is where the gene-specific signal is. A half-site like AGGTCA occurs every few kilobases by chance, so nuclear receptors surface for almost any promoter; these factors surfaced for *this* one. The routes that can reach them are the upstream-network, stabilisation and protein-level sections, not the ligand tables.

### Top-ranked transcription factors

| Rank | TF | Role | NR | Other genes it controls | Cells need it? | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NFYB | activator |  | ~5000 | yes, most | 1 | 1.0 | -60 | 2 | 0 | 81.2 |
| 2 | HNF4A | activator | yes | ~300 | a few lines | 1 | 1.0 | -3680 | 4 | 0 | 78.3 |
| 3 | KLF5 | activator |  | ~5000 | some lines | 2 | 0.997 | +359 | 1 | 0 | 77.7 |
| 4 | STAT3 | activator |  | ~1200 | a few lines | 1 | 1.0 | +309 | 5 | 0 | 75.3 |
| 5 | FOS | activator |  | ~5000 | a few lines | 3 | 1.0 | -2778 | 1 | 0 | 75.0 |
| 6 | JUNB | activator |  | ~5000 | some lines | 5 | 1.0 | -2778 | 4 | 0 | 75.0 |
| 7 | JUND | activator |  | not counted | a few lines | 3 | 1.0 | -2778 | 3 | 0 | 75.0 |
| 8 | FOSL2 | activator |  | not counted | some lines | 4 | 0.998 | -2779 | 3 | 0 | 75.0 |
| 9 | WT1 | activator |  | not counted | a few lines | 1 | 0.921 | -320 | 1 | 0 | 72.7 |
| 10 | CEBPA | activator |  | not counted | a few lines | 1 | 0.947 | -372 | 2 | 0 | 72.7 |
| 11 | JUN | dual |  | ~5000 | some lines | 3 | 1.0 | -2778 | 5 | 0 | 69.0 |
| 12 | E2F4 | activator |  | ~5000 | a few lines | 1 | 0.916 | +675 | 2 | 0 | 68.6 |
| 13 | ELK1 | activator |  | ~5000 | a few lines | 1 | 0.904 | -652 | 1 | 0 | 68.4 |
| 14 | FOSB | activator |  | 189 | a few lines | 1 | 1.0 | -2778 | 2 | 0 | 68.3 |
| 15 | SP1 | activator |  | ~5000 | some lines | 1 | 1.0 | -4929 | 6 | 0 | 68.3 |
| 16 | ONECUT1 | activator |  | not counted | a few lines | 1 | 0.968 | -4472 | 1 | 0 | 67.4 |
| 17 | GATA1 | activator |  | ~300 | a few lines | 1 | 0.926 | -4827 | 1 | 0 | 66.1 |
| 18 | GLIS3 | dual |  | not counted | no | 1 | 0.863 | -260 | 1 | 0 | 66.1 |
| 19 | FOXO1 | dual |  | ~1200 | a few lines | 2 | 0.977 | -3548 | 2 | 0 | 65.0 |
| 20 | CTCF | repressor |  | ~5000 | yes, nearly all | 3 | 0.804 | +528 | 3 | 0 | 62.6 |
| 21 | BACH1 | dual |  | not counted | a few lines | 1 | 1.0 | -2778 | 1 | 0 | 62.3 |
| 22 | ERF | repressor |  | not counted | no | 2 | 0.856 | -432 | 1 | 0 | 61.9 |
| 23 | SP2 | repressor |  | not counted | some lines | 2 | 1.0 | -4929 | 2 | 0 | 61.2 |
| 24 | KLF1 | unknown |  | not counted | a few lines | 1 | 1.0 | -312 | 1 | 0 | 60.2 |
| 25 | FOSL1 | unknown |  | not counted | some lines | 3 | 1.0 | -2778 | 1 | 0 | 60.0 |

**Other genes it controls** is the collateral question. A drug acting through a factor moves everything that factor regulates, so this column says how wide the effect will be. Read it against the score: a strong route through a factor that controls two thousand genes is a worse starting point than a weaker route through one that controls thirty.

A count written `~N` is a bundled estimate for a factor the regulon file does not list (14 here), not a count of curated edges. Treat it as an order of magnitude.

`not counted` means the regulon table has no entry and there is no bundled estimate either (55 here). It is not a small number. Most of the roughly 1,600 human transcription factors have no curated regulon at all.

Counts come from literature curation, so they measure what has been published as much as what the factor does. SP1 has thousands of recorded targets partly because it regulates many genes and partly because thousands of people assayed it.

### What a cell loses without each factor (DepMap)

The regulon column counts how many genes a factor controls, from literature curation, which measures how much anyone has published as much as what the factor does. This is the same collateral question answered by measurement: CRISPR knockout of each factor across roughly 1,100 cancer cell lines.

It also reaches factors the curated tables cannot. A factor nobody has written a regulon for still has a knockout phenotype here.

88 of 89 factor(s) were scored in DepMap. 1 (CTCF) is needed by nearly every cell line tested, so strong inhibition is likely cytotoxic. 65 (BACH1, BACH2, CEBPA, CUX1, E2F4) show a dependency confined to a minority of lines, which is the informative case: the requirement has a context, and whether your tissue shares it is the question to ask next. 19 (ERF, ESRRB, FIGLA, GLIS1, GLIS3) were dispensable throughout. That bounds the proliferative cost of removing them and says nothing about a differentiated tissue.

| Factor | Cells need it? | Mean effect | Lines dependent | Most dependent lineage |
|---|---|---|---|---|
| CTCF | essential in nearly every line | -1.68 | 1147/1150 (100%) | Biliary Tract (-1.80) |
| NFYB | essential in most lines | -0.70 | 961/1150 (84%) | Thyroid (-0.84) |
| KLF5 | essential in a substantial minority | -0.41 | 307/1150 (27%) | Bowel (-0.99) |
| FOSL1 | essential in a substantial minority | -0.41 | 368/1150 (32%) | Pancreas (-0.64) |
| TFAP4 | essential in a minority of lines | -0.35 | 263/1150 (23%) | Myeloid (-0.54) |
| GLI1 | essential in a minority of lines | -0.30 | 85/1150 (7%) | Soft Tissue (-0.43) |
| JUN | essential in a minority of lines | -0.29 | 192/1150 (17%) | Pleura (-0.54) |
| SPIC | essential in a minority of lines | -0.28 | 66/1150 (6%) | Prostate (-0.36) |
| JUNB | essential in a minority of lines | -0.27 | 163/1150 (14%) | Pancreas (-0.43) |
| SP2 | essential in a minority of lines | -0.26 | 100/1150 (9%) | Bladder/Urinary Tract (-0.34) |
| STAT5A | essential in a handful of lines | -0.25 | 50/1150 (4%) | Peripheral Nervous System (-0.31) |
| MGA | essential in a handful of lines | -0.21 | 44/1150 (4%) | Head and Neck (-0.29) |
| STAT3 | essential in a handful of lines | -0.19 | 39/1150 (3%) | Cervix (-0.31) |
| STAT5B | essential in a handful of lines | -0.19 | 27/1150 (2%) | Myeloid (-0.49) |
| ELK1 | essential in a handful of lines | -0.17 | 3/1150 (0%) | Ovary/Fallopian Tube (-0.21) |
| FOSL2 | essential in a minority of lines | -0.17 | 60/1150 (5%) | Ovary/Fallopian Tube (-0.29) |
| JUND | essential in a handful of lines | -0.16 | 28/1150 (2%) | Thyroid (-0.25) |
| KLF1 | essential in a handful of lines | -0.15 | 17/1150 (2%) | Myeloid (-0.26) |
| PATZ1 | essential in a handful of lines | -0.15 | 9/1150 (1%) | Peripheral Nervous System (-0.29) |
| KLF10 | essential in a handful of lines | -0.14 | 33/1150 (3%) | Eye (-0.25) |
| FOXO1 | essential in a handful of lines | -0.13 | 29/1150 (2%) | Lymphoid (-0.30) |
| E2F4 | essential in a handful of lines | -0.13 | 14/1150 (1%) | Lymphoid (-0.20) |
| STAT6 | essential in a handful of lines | -0.12 | 3/1150 (0%) | Eye (-0.20) |
| ZNF281 | essential in a handful of lines | -0.11 | 29/1150 (2%) | Myeloid (-0.26) |
| NFATC4 | essential in a handful of lines | -0.11 | 2/1150 (0%) | Myeloid (-0.15) |

**Mean effect** is on DepMap's own scale: 0 means knockout changed nothing, and -1 is the median of known common-essential genes. A line counts as dependent below -0.5.

#### What this does not mean

- **These are cancer cell lines dividing in a dish.** Essentiality here is a statement about proliferating transformed cells. A factor whose job is in a post-mitotic or differentiated cell looks dispensable for reasons that have nothing to do with its importance in your tissue.
- **Knockout is not inhibition.** The protein is removed entirely, including scaffolding roles no ligand would touch. Essential genes are drugged successfully at partial occupancy all the time, so a dependency is a reason for care rather than a veto.
- **Nothing here is about your gene.** This bounds the cost of acting on a factor. Whether acting on it moves your gene is the rest of the report.
- **A dispensable result is not permission.** It bounds the proliferative cost and nothing else.

**The selective ones are the interesting ones.** BACH1, BACH2, CEBPA, CUX1, E2F4, EBF2 matter in a minority of lines rather than all of them, so the requirement has a context. The lineage column is where to look for what that context is, and whether your tissue shares it is the question this panel cannot answer.

### Drugs that behave like losing the factor

Every other route in this report starts from something somebody wrote down - a curated ligand, a motif, a structure, a signed edge. None of them can surface a drug nobody has connected to the factor. This one reads a measured phenotype instead: if a compound hurts the same cell lines that knocking the factor out hurts, it behaves as though it removes the factor.

| Factor | Drug | Partial r | Raw r | Lines | q | Behaves as |
|---|---|---|---|---|---|---|
| KLF5 | **GEFITINIB** | 0.319 | 0.3006 | 459 | 1e-06 | phenocopies |
| KLF5 | **LAPATINIB** | 0.3143 | 0.2865 | 465 | 1e-06 | phenocopies |
| KLF5 | **OSIMERTINIB** | 0.3046 | 0.273 | 466 | 1e-06 | phenocopies |
| KLF5 | **IBRUTINIB** | 0.3002 | 0.2699 | 460 | 1e-06 | phenocopies |
| KLF5 | **MOBOCERTINIB** | 0.2776 | 0.2609 | 703 | 1e-06 | phenocopies |
| KLF5 | **ALPELISIB** | 0.2622 | 0.2397 | 465 | 2e-06 | phenocopies |

**Partial r** controls for how easily each cell line dies under the whole panel. Some lines die under almost anything and those same lines suffer most knockouts, which puts a correlation between almost any drug and almost any gene. The control moves a number either way - it removes a shared component that can inflate a correlation or mask one - so where the two columns differ a lot the raw value was mostly about the cell lines. Estimates use at least 25 shared lines.

**This is not a target.** A correlated drug may act on the factor, on anything else the same lines need, or on a pathway the factor sits in, and the first is the least likely of the three for any one hit. Nothing here records any of them binding anything: the finding is a reason to run a binding assay, not a result from one.

**A proliferation readout.** Both panels score growth of transformed cells in a dish, so a factor whose role is in a differentiated or post-mitotic cell looks silent here for reasons that have nothing to do with its importance.

### Every route that acts on this gene

225 routes across 7 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

**Why established binders sit at tier 2, just under the direct agents.** A drug already measured to bind one of your factors has half its chain settled: the drug-to-receptor step has assay data behind it, so the pharmacology is de-risked and only the receptor-to-gene step is a hypothesis. A curated tier-3 ligand is annotation on *both* halves. But they stay below tier 1, because a direct agent acts on your gene itself and skips the receptor-to-gene inference entirely — a stronger claim however well-characterised the binder is. The cost of the indirect route is in the row: everything else that factor regulates moves with your gene.

- 27 × curated — mechanism annotated by a human curator
- 193 × predicted — a model says so
- 5 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 3 | TF: abundance change | teniposide | STAT3 | 1 | raises STAT3 transcript in 10 profile(s) | 1.0 | 0.45 | 0.5 | 54.6 |
| 3 | TF: abundance change | selumetinib | STAT3 | 1 | raises STAT3 transcript in 10 profile(s) | 1.0 | 0.45 | 0.5 | 54.6 |
| 3 | TF: known regulator | Cholecalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Paricalcitol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Doxercalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Ergocalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Calcifediol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: known regulator | Calcitriol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 3 | TF: abundance change | dasatinib | STAT3 | 1 | raises STAT3 transcript in 7 profile(s) | 0.7 | 0.45 | 0.5 | 45.0 |
| 3 | TF: abundance change | vorinostat | FOSL2 | 1 | raises FOSL2 transcript in 88 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | FOS | 1 | raises FOS transcript in 29 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | JUNB | 1 | raises JUNB transcript in 28 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | FOSL2 | 1 | raises FOSL2 transcript in 21 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | JUND | 1 | raises JUND transcript in 14 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | vorinostat | JUND | 1 | raises JUND transcript in 13 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | vorinostat | STAT3 | 1 | raises STAT3 transcript in 10 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | KLF5 | 1 | raises KLF5 transcript in 9 profile(s) | 0.9 | 0.1 | 0.5 | 41.6 |
| 3 | TF: abundance change | mitoxantrone | STAT3 | 1 | raises STAT3 transcript in 8 profile(s) | 0.8 | 0.1 | 0.5 | 38.4 |
| 3 | chromatin-level | Vorinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Romidepsin | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Belinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Panobinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Valproic acid | chromatin (genome-wide) | 3 | weak HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Azacitidine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Decitabine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Tazemetostat | chromatin (genome-wide) | 3 | EZH2 inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Givinostat Hydrochloride | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 4 | upstream signalling | Treprostinil | PTGIR | 2 | activate PTGIR → STAT3 (1 hop) | 0.75 | 0.851 | 0.6 | 49.9 |
| 4 | upstream signalling | Treprostinil Diolamine | PTGIR | 2 | activate PTGIR → STAT3 (1 hop) | 0.75 | 0.851 | 0.6 | 49.9 |
| 4 | upstream signalling | Selexipag | PTGIR | 2 | activate PTGIR → STAT3 (1 hop) | 0.75 | 0.851 | 0.6 | 49.9 |

Best in each category:

- **TF: abundance change** — teniposide (tier 3, score 54.6)
- **TF: known regulator** — Cholecalciferol (tier 3, score 45.2)
- **chromatin-level** — Vorinostat (tier 3, score 23.4)
- **upstream signalling** — Treprostinil (tier 4, score 49.9)
- **literature: reported effect** — Acetaminophen (tier 4, score 33.3)
- **direct: predicted product binder** — PARICALCITOL (tier 5, score 62.3)
- **promoter DNA binder** — Plicamycin (tier 5, score 31.2)

*No route found in: direct: gene or transcript, direct: gene product, TF: orthosteric ligand, TF: inferred activity change, TF: predicted pocket binder, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, TF: degradation blocker, off-pathway pocket binder, measured phenocopy, empirical signature, virtual screen hypothesis. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

| Category | What it changes | Your drugs in it | Assay that sees it |
|---|---|---|---|
| TF: abundance change | **abundance** | teniposide, selumetinib, dasatinib, vorinostat … +1 | RT-qPCR and a western for the FACTOR first, then RT-qPCR for your gene |
| TF: known regulator | **abundance** | Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol … +30 | RT-qPCR for your gene after treatment; check the factor's own activity (a reporter or a known target) moved first |
| chromatin-level | **abundance, genome-wide** | Vorinostat, Romidepsin, Belinostat, Panobinostat … +5 | RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened |
| upstream signalling | **abundance** | Treprostinil, Treprostinil Diolamine, Selexipag, Iloprost … +49 | RT-qPCR then western/ELISA; confirm the pathway node moved too, or you cannot attribute the effect |
| literature: reported effect | **abundance (reported)** | Acetaminophen, Calcifediol, Calcitriol, Cholecalciferol … +14 | RT-qPCR in your cells; read the paper for the dose, cell type and species first |
| direct: predicted product binder | **activity** | PARICALCITOL | a binding assay against the protein first, then a functional assay for direction - RT-qPCR will not see this route at all, and a western blot reads it as no effect |
| promoter DNA binder | **abundance** | Plicamycin, Mitoxantrone, Doxorubicin, Pentamidine | an EMSA or footprint on this promoter first, then a reporter, then RT-qPCR - and a counter-screen on an unrelated promoter carrying the same feature, because selectivity is the whole question here |

*Every category present here works through transcription, so all of them move your gene's product **abundance** and an RT-qPCR is the right first readout. The one exception in the taxonomy — a drug on the gene product itself, which changes activity only — did not appear in this run.*

**TF: abundance change** (13 found)

- *Yours:* teniposide, selumetinib, dasatinib, vorinostat, mitoxantrone
- A compound measured to change how much of the factor's own transcript exists, rather than what the factor does.
- *Why:* The only route here that needs the factor to be expressed rather than druggable, so it reaches factors with no ligand, no structure and no curated upstream kinase - which is most of them.
- *Cost:* More factor is not more factor ACTIVITY. For a ligand-activated receptor the two are nearly unrelated, and transcript is not protein. Broad agents that shut down transcription move every factor and are flagged.
- *Differs:* Unlike 'empirical signature', which is a compound measured against YOUR GENE, this was measured against the FACTOR - one further inferential step from what you care about.

**TF: known regulator** (36 found)

- *Yours:* Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol, Calcifediol, Calcitriol, Bexarotene, Acitretin … +26 more
- Drugs that act on a factor CollecTRI or perturbation data records as raising or lowering your gene, whether or not the motif scan ranked that factor.
- *Why:* Reaches factors the motif scan misses, such as those that bind through a partner or a distant enhancer.
- *Cost:* The link comes from a database or a knockout in other cells, so it may not hold in yours. The drug's effect on the factor is taken from its recorded action.
- *Differs:* Starts from the literature on your gene, not from its promoter sequence. CollecTRI and CTD share sources, so agreement between them is not independent.

**chromatin-level** (9 found)

- *Yours:* Vorinostat, Romidepsin, Belinostat, Panobinostat, Valproic acid, Azacitidine, Decitabine, Tazemetostat … +1 more
- HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.
- *Why:* Reliable positive control that your reporter *can* be induced.
- *Cost:* Specificity ~0.05. These move thousands of promoters. Almost never an acceptable therapeutic route on that basis alone.
- *Differs:* Not targeted at your gene in any sense; it is in the table so you can see what a non-specific option scores.

**upstream signalling** (144 found)

- *Yours:* Treprostinil, Treprostinil Diolamine, Selexipag, Iloprost, Iloprost Tromethamine, Avatrombopag, Avatrombopag Maleate, Lusutrombopag … +45 more
- Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- *Why:* Turns 'this factor has no drugs' into a real list. Most TFs are only reachable this way.
- *Cost:* Specificity falls with every step: a node three hops up reaches your gene through one branch and thousands of other genes through the rest.
- *Differs:* The only category where the drug never touches your gene or its TF — read the printed path and judge it yourself.

**literature: reported effect** (18 found)

- *Yours:* Acetaminophen, Calcifediol, Calcitriol, Cholecalciferol, Dexamethasone, Ketoconazole, Leflunomide, Rifampin … +10 more
- A chemical CTD records as raising or lowering this gene's mRNA, with no mechanism stated.
- *Why:* It is the only route that can name a drug nothing here can explain, which includes most of what CTD holds: stress, feedback and toxicity effects.
- *Cost:* Most records rest on one paper, often at a toxicology dose in one cell type. Across 36 genes only about 3% of them were reproduced in LINCS, and more papers did not raise that. Famous drugs are over-represented.
- *Differs:* Every other category predicts an effect from a mechanism. This reports an effect somebody saw, and cannot suggest a drug nobody has tested on this gene.

**direct: predicted product binder** (1 found)

- *Yours:* PARICALCITOL
- A compound predicted to bind this gene's own protein and not recorded as acting on it - either because its chemistry resembles the protein's known ligands, or because it docks into the protein's pocket.
- *Why:* It is the only route here aimed at the gene's own product rather than at a factor upstream of it, so it skips every inference about whether moving a factor moves your gene. And unlike the curated product drugs it can propose a compound no one has annotated.
- *Cost:* Predicted binding and nothing more. The model is refused below a cross-validated AUC floor, which bounds how badly it can be ranking, but a hit is a hypothesis for a binding assay. Direction is separate again: whether binding inhibits or activates the protein decides whether it serves your goal.
- *Differs:* Unlike 'direct: gene product', nothing was measured or curated for this compound at this protein - that category is annotation, this is a model. Unlike 'TF: predicted pocket binder', the pocket is on the gene's own product, so no inference about a factor controlling the gene is involved.

**promoter DNA binder** (4 found)

- *Yours:* Plicamycin, Mitoxantrone, Doxorubicin, Pentamidine
- A small molecule with a documented preference for a sequence feature this promoter carries - a quadruplex, a GC-box, an AT tract - proposed because the feature is enriched against shuffles of this same sequence.
- *Why:* It is the only route that acts on the promoter itself rather than on a protein, so it needs no druggable factor and no annotation. Where a promoter depends on one feature, occupying it is a direct intervention.
- *Cost:* Sequence preference is not selectivity. These agents bind the feature wherever it occurs and most of the genome's copies are not in your gene, so the dose that engages your promoter engages thousands of others. Nearly every approved member is cytotoxic, which is how it works rather than a side effect. And a motif is not a folded structure in a cell.
- *Differs:* Unlike every TF route, no protein is targeted - which removes the question of whether the factor controls your gene and replaces it with whether the molecule can tell your promoter from any other.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **direct: gene product** — Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- **TF: orthosteric ligand** — Binds the transcription factor's natural ligand pocket.
- **TF: inferred activity change** — A compound whose own LINCS signature moves this factor's known targets, found without any record that it binds the factor.
- **TF: predicted pocket binder** — A compound predicted by docking to fit the factor's ligand pocket, from a structure rather than from any measured activity.
- **TF: allosteric modulator** — Binds a site away from the ligand pocket, amplifying or damping the endogenous ligand's effect (PAM/NAM).
- **TF: stabiliser or glue** — Locks the factor in an active conformation, or glues it to a partner protein.
- **TF: endogenous ligand** — Changes how much natural ligand exists, by drugging the enzyme that makes or destroys it.
- **TF: degradation blocker** — Stops the factor being destroyed, so more of it accumulates.
- **off-pathway pocket binder** — A docking run against a protein that is not one of the factors found at this promoter. Usually an upstream node someone chose to dock by hand.
- **measured phenocopy** — A compound whose effect on cell viability across the DepMap panel tracks what happens when this factor is knocked out, so it behaves as though it removes the factor.
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **virtual screen hypothesis** — An approved drug predicted to bind one of your TFs, never tested against it.

### What these drugs are

Every route in this report is off-label by construction: the drug was approved for one thing and is being considered for another. The approved use is the baseline that departure is measured from, so it is listed here once per drug rather than repeated in each table above.

| Drug | What it is | Approved to treat |
|---|---|---|
| Acitretin | retinoid | severe psoriasis |
| Azacitidine | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia |
| Beclomethasone | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Bexarotene | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma |
| Calcifediol | 25-hydroxyvitamin D | vitamin D deficiency |
| Calcitriol | active vitamin D, VDR agonist | hypocalcaemia and secondary hyperparathyroidism in renal failure |
| Cholecalciferol | vitamin D3, a precursor requiring activation | vitamin D deficiency and supplementation |
| Clofibrate | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) |
| Decitabine | DNA methyltransferase inhibitor | myelodysplastic syndrome |
| Dexamethasone | potent glucocorticoid receptor agonist | inflammation, cerebral oedema, myeloma, and as an antiemetic |
| Dextrothyroxine Sodium | T4 stereoisomer | formerly hyperlipidaemia - withdrawn for cardiac mortality |
| Diethylstilbestrol Diphosphate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Doxercalciferol | vitamin D prodrug | secondary hyperparathyroidism in dialysis |
| Doxorubicin | anthracycline topoisomerase II inhibitor | many solid and haematological cancers |
| Elafibranor | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) |
| Ergocalciferol | vitamin D receptor agonist class | secondary hyperparathyroidism in renal disease, osteoporosis, and plaque psoriasis topically |
| Estradiol Acetate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Estradiol Cypionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Etretinate | retinoid | severe psoriasis - WITHDRAWN, teratogenic |
| Fenofibrate | PPARalpha agonist | hypertriglyceridaemia |
| Gemfibrozil | PPARalpha agonist | hypertriglyceridaemia |
| Hydrocortisone | glucocorticoid and mineralocorticoid | adrenal insufficiency; topical dermatitis |
| Hydrocortisone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Isotretinoin | retinoid | severe nodular acne |
| Leflunomide | broad immunosuppressant class | transplant rejection and autoimmune disease |
| Levothyroxine Sodium | synthetic T4 | hypothyroidism |
| Liothyronine Sodium | synthetic T3 | hypothyroidism, myxoedema coma |
| Mometasone Furoate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Palovarotene | retinoid class | acne, psoriasis and other keratinisation disorders topically; acute promyelocytic leukaemia and cutaneous lymphoma systemically |
| Paricalcitol | selective VDR activator | secondary hyperparathyroidism in chronic kidney disease |
| Plicamycin | RNA synthesis inhibitor, binds GC-rich DNA | formerly hypercalcaemia of malignancy - withdrawn |
| Prednisolone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Raloxifene Hydrochloride | selective oestrogen receptor modulator | osteoporosis; breast cancer risk reduction |
| Romidepsin | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma |
| Testosterone | androgen receptor agonist | male hypogonadism |
| Tretinoin | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically |
| mitoxantrone | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis |
| selumetinib | MEK1/2 inhibitor | neurofibromatosis type 1 plexiform neurofibroma |
| vorinostat | histone deacetylase inhibitor | cutaneous T-cell lymphoma |

**No description established for 56**: Acetaminophen, Adenosine, Almotriptan Malate, Avatrombopag, Bazedoxifene Acetate, Belinostat, Bethanechol Chloride, Carboprost Tromethamine, Cevimeline Hydrochloride, Chenodiol, Chlorotrianisene, Choline Fenofibrate, Cinacalcet Hydrochloride, Copper Sulfate, Dinoprost Tromethamine, Docetaxel Anhydrous, Eletriptan Hydrobromide, Eltrombopag Olamine, Estramustine Phosphate Sodium, Etrasimod, Fenofibric Acid, Fingolimod Hydrochloride, Fingolimod Lauryl Sulfate, Givinostat Hydrochloride, Iloprost…. That means this package could not find one offline, not that the compound has no approved use. Several are environmental agents and research chemicals from the CTD tables, which have no indication because they are not drugs.

Descriptions come from: 12 class, 27 curated, 56 unknown. None of these is a label. The curated text is a summary written for orientation, ChEMBL's indication list includes uses studied but never approved, and a class line gives the kind of drug rather than its indication. Check the prescribing information before acting on any row.

### Cell-cycle compatibility

Tissue proliferation state: **unknown**

1 factor(s) here are cell-cycle restricted (E2F4). Set --proliferation to have this scored rather than just flagged.

| TF | Phase | Duty cycle | Effect multiplier |
|---|---|---|---|
| E2F4 | G1 | 46% | 1.0 |

*A partial duty cycle dilutes the population-level effect in asynchronous cells: per-cell efficacy can be much higher than a bulk RT-qPCR suggests. Synchronise or use a phase-resolved readout before calling it a failure.*

### Upstream routes to factors with no direct drug

Most TFs have no approved ligand. These paths walk upstream through signed causal edges — kinases, phosphatases, receptors — to proteins that *are* druggable. Signs multiply along the path, so an inhibitor of an inhibitor counts as activating.

Distance buys druggability and spends selectivity: a node three hops up reaches your gene through one branch and everything else through the others. Read the fan-out column.

**NFYB** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| PTGER3 | 3 | 9 | PTGER3 --> SRC —\| TP53 —\| NFYB | activate | Misoprostol | 15.7 |
| S1PR1 | 3 | 12 | S1PR1 --> SRC —\| TP53 —\| NFYB | activate | Ponesimod, Siponimod, Etrasimod | 14.8 |
| NR3C1 | 3 | 16 | NR3C1 --> HDAC1 —\| TP53 —\| NFYB | activate | Hydrocortisone Acetate, Prednisolone Acetate, Mometasone Furoate | 13.7 |
| ESR1 | 3 | 37 | ESR1 --> SRC —\| TP53 —\| NFYB | activate | Estradiol Acetate, Diethylstilbestrol Diphosphate, Toremifene Citrate | 10.0 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 12 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 16 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 37 downstream branches

**KLF5** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| TUBB | 3 | 2 | TUBB —\| SMAD2 --> SMURF2 —\| KLF5 | activate | Docetaxel Anhydrous | 18.3 |
| CASR | 3 | 2 | CASR --> FYN --> PRKCD --> KLF5 | activate | Cinacalcet Hydrochloride | 18.3 |
| RARG | 3 | 5 | RARG --> AKT1 --> PRKCD --> KLF5 | activate | Acitretin, Palovarotene, Tretinoin | 17.1 |
| CHRM3 | 3 | 5 | CHRM3 --> GNAQ —\| GSK3B —\| KLF5 | activate | Pilocarpine Hydrochloride, Bethanechol Chloride | 17.1 |
| HTR2B | 3 | 5 | HTR2B --> GNAQ —\| GSK3B —\| KLF5 | activate | Methylergonovine Maleate | 17.1 |
| RARB | 3 | 6 | RARB --> AKT1 --> PRKCD --> KLF5 | activate | Acitretin, Tretinoin | 16.7 |
| HTR2A | 3 | 6 | HTR2A --> GNAQ —\| GSK3B —\| KLF5 | activate | — | 16.7 |
| CHRM1 | 3 | 6 | CHRM1 --> GNAQ —\| GSK3B —\| KLF5 | activate | Pilocarpine Hydrochloride, Xanomeline, Cevimeline Hydrochloride | 16.7 |
| MTNR1A | 3 | 7 | MTNR1A --> GNAQ —\| GSK3B —\| KLF5 | activate | Ramelteon, Tasimelteon | 16.4 |
| PTGFR | 3 | 7 | PTGFR --> GNAQ —\| GSK3B —\| KLF5 | activate | Dinoprost Tromethamine, Carboprost Tromethamine | 16.3 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 5 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 6 downstream branches

**STAT3** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| PTGIR | 1 | 7 | PTGIR --> STAT3 | activate | Treprostinil, Treprostinil Diolamine, Selexipag | 71.4 |
| MPL | 2 | 2 | MPL --> JAK2 --> STAT3 | activate | Avatrombopag, Avatrombopag Maleate, Lusutrombopag | 33.6 |
| GHR | 2 | 4 | GHR --> JAK2 --> STAT3 | activate | — | 32.1 |
| CXCR4 | 2 | 5 | CXCR4 --> JAK2 --> STAT3 | activate | Plerixafor | 31.5 |
| RARG | 2 | 5 | RARG --> AKT1 --> STAT3 | activate | Acitretin, Palovarotene, Tretinoin | 31.3 |
| RARB | 2 | 6 | RARB --> AKT1 --> STAT3 | activate | Acitretin, Tretinoin | 30.6 |
| NOD2 | 2 | 9 | NOD2 --> MAP3K14 --> STAT3 | activate | — | 28.9 |
| PTGER3 | 2 | 9 | PTGER3 --> EGFR --> STAT3 | activate | Misoprostol | 28.8 |
| PPARA | 2 | 11 | PPARA —\| FOXO1 —\| STAT3 | activate | Choline Fenofibrate, Elafibranor, Gemfibrozil | 27.6 |
| S1PR1 | 2 | 12 | S1PR1 --> SRC --> STAT3 | activate | Ponesimod, Siponimod, Etrasimod | 27.6 |

- ⚠ path passes through hub(s): AKT1
- ⚠ path passes through hub(s): EGFR
- ⚠ path passes through hub(s): SRC

**FOS** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| MTNR1A | 2 | 7 | MTNR1A --> MAPK1 --> FOS | activate | Ramelteon, Tasimelteon | 30.2 |
| NOD2 | 2 | 9 | NOD2 --> MAPK14 --> FOS | activate | — | 29.2 |
| PTGER3 | 2 | 9 | PTGER3 --> SRC --> FOS | activate | Misoprostol | 28.9 |
| S1PR1 | 2 | 12 | S1PR1 --> SRC --> FOS | activate | Ponesimod, Siponimod, Etrasimod | 27.3 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK1 --> FOS | activate | Acetaminophen | 18.7 |
| ESR1 | 2 | 37 | ESR1 --> SRC --> FOS | activate | Estradiol Acetate, Diethylstilbestrol Diphosphate, Toremifene Citrate | 18.4 |
| CASR | 3 | 2 | CASR --> FYN --> MAPK14 --> FOS | activate | Cinacalcet Hydrochloride | 18.3 |
| S1PR5 | 3 | 2 | S1PR5 --> GNAI1 --> SRC --> FOS | activate | Ozanimod, Etrasimod Arginine, Fingolimod Lauryl Sulfate | 18.3 |
| ADORA3 | 3 | 2 | ADORA3 --> GNAI1 --> SRC --> FOS | activate | Adenosine | 18.3 |
| HTR1F | 3 | 2 | HTR1F --> GNAI1 --> SRC --> FOS | activate | Almotriptan Malate, Eletriptan Hydrobromide, Lasmiditan | 18.3 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ path passes through hub(s): MAPK1

**JUNB** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| NOD2 | 2 | 9 | NOD2 --> MAPK8 --> JUNB | activate | — | 30.9 |
| PTGIR | 2 | 7 | PTGIR --> STAT3 --> JUNB | activate | Treprostinil, Treprostinil Diolamine, Selexipag | 30.3 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK8 --> JUNB | activate | Acetaminophen | 18.7 |
| MPL | 3 | 2 | MPL --> JAK2 --> STAT3 --> JUNB | activate | Avatrombopag, Avatrombopag Maleate, Lusutrombopag | 18.3 |
| CASR | 3 | 2 | CASR --> FYN —\| ITCH —\| JUNB | activate | Cinacalcet Hydrochloride | 18.3 |
| GHR | 3 | 4 | GHR --> JAK2 --> STAT3 --> JUNB | activate | — | 17.5 |
| CXCR4 | 3 | 5 | CXCR4 --> JAK2 --> STAT3 --> JUNB | activate | Plerixafor | 17.1 |
| RARG | 3 | 5 | RARG --> AKT1 --> STAT3 --> JUNB | activate | Acitretin, Palovarotene, Tretinoin | 17.1 |
| CHRM3 | 3 | 5 | CHRM3 --> GNAQ —\| GSK3B —\| JUNB | activate | Pilocarpine Hydrochloride, Bethanechol Chloride | 17.1 |
| HTR2B | 3 | 5 | HTR2B --> GNAQ —\| GSK3B —\| JUNB | activate | Methylergonovine Maleate | 17.1 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 4 downstream branches

**JUND** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| MTNR1A | 2 | 7 | MTNR1A --> MAPK1 --> JUND | activate | Ramelteon, Tasimelteon | 30.4 |
| NOD2 | 2 | 9 | NOD2 --> MAPK8 --> JUND | activate | — | 30.3 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK8 --> JUND | activate | Acetaminophen | 18.8 |
| RARG | 3 | 5 | RARG --> AKT1 --> BLVRA --> JUND | activate | Acitretin, Palovarotene, Tretinoin | 17.1 |
| RARB | 3 | 6 | RARB --> AKT1 --> BLVRA --> JUND | activate | Acitretin, Tretinoin | 16.8 |
| PTGIR | 3 | 7 | PTGIR --> STAT3 --> MAPK8 --> JUND | activate | Treprostinil, Treprostinil Diolamine, Selexipag | 16.3 |
| PTGER3 | 3 | 9 | PTGER3 --> SRC --> MAPK1 --> JUND | activate | Misoprostol | 15.7 |
| S1PR1 | 3 | 12 | S1PR1 --> LYN --> MAPK3 --> JUND | activate | Ponesimod, Siponimod, Etrasimod | 14.9 |
| RARA | 3 | 15 | RARA --> AKT1 --> BLVRA --> JUND | activate | Acitretin, Tretinoin | 14.0 |
| NR3C1 | 3 | 16 | NR3C1 --> LCK --> MAPK3 --> JUND | activate | Hydrocortisone Acetate, Prednisolone Acetate, Mometasone Furoate | 14.0 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 12 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 15 downstream branches

### Mechanism options: orthosteric vs allosteric

Allosteric sites are not constrained to bind the endogenous ligand, so they diverge across a receptor family where the orthosteric pocket cannot. That usually buys subtype selectivity — at the cost of a bounded maximal effect and a hard dependence on endogenous tone.

**HNF4A**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| orthosteric agonist | 1 | full | typically lower | no | no |  |

*No approved allosteric modulator at this target; only orthosteric pharmacology, with the usual subtype-selectivity problem.*

**Before pursuing any allosteric option, verify that the pathway has endogenous tone in your cell type.** A PAM with no ligand to modulate does nothing at all, and that is a measurement, not a prediction — nothing in this package can check it for you.

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

#### No stabilisation route for 7 factor(s)

FOS, FOSL2, JUNB, JUND, KLF5, NFYB, STAT3

These have no ligand pocket, no direct drug, and **nothing curated against them or their degradation machinery**. For these factors stabilisation is a research programme, not a route you can take: making it real needs the E3 ligase that degrades the factor *and* a ligand for that ligase, and neither is published here.

The mechanism class is real — it has shipped approved drugs, listed once below — but those act on **other proteins entirely**, and none of them is evidence that any factor above is tractable.

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Cell-type reach

Three things must overlap: where the drug goes, where the gene is expressed, and where the transcription factor is. Any empty intersection kills the route regardless of motif quality or potency.

| TF | Co-expressed with gene | Shared tissues |
|---|---|---|
| NFYB | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| HNF4A | yes (1) | kidney |
| KLF5 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| STAT3 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| FOS | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| JUNB | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| JUND | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| FOSL2 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| WT1 | yes (7) | cervix, uterine, endometrium, kidney, placenta, seminal vesi |
| CEBPA | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |

### Human genetic evidence

Human genetic variation in CYP24A1 associates with: hypercalcemia, infantile, nephrolithiasis, atopic eczema. Read this as a preview of what modulating the gene may do, including the adverse effects — this is the reasoning that anticipated PCSK9 inhibition years before the trials.

| Phenotype | Genetic score | Overall |
|---|---|---|
| hypercalcemia, infantile | 0.902 | 0.814 |
| nephrolithiasis | 0.92 | 0.568 |
| atopic eczema | 0.913 | 0.565 |
| urolithiasis | 0.699 | 0.531 |
| ureterolithiasis | 0.846 | 0.515 |
| bladder calculus | 0.806 | 0.49 |
| hereditary disease | 0.773 | 0.472 |
| aortic stenosis | 0.677 | 0.411 |
| vitamin D deficiency | 0.593 | 0.367 |
| deep vein thrombosis | 0.59 | 0.359 |
| COVID-19 | 0.576 | 0.355 |
| allergic disease | 0.578 | 0.352 |
| Abnormal thrombosis | 0.578 | 0.352 |

#### Does variation in the factor actually move this gene? (trans-eQTL)

The strongest evidence available anywhere in this report. A motif says a factor *could* bind; a trans-eQTL says natural genetic variation near it **does** change your gene's expression, in real human tissue, with exposure randomised by meiosis.

**No factor can be supported this way for CYP24A1.** CYP24A1 has no significant trans-eQTL in eQTLGen at all, so the test returns nothing for every factor here. That is a property of the dataset's power and of the tissue (whole blood), not evidence that CYP24A1 is unregulated.

*No factor tested (NFYB, HNF4A, KLF5, STAT3, FOS, JUNB) has trans-eQTL support for this gene.*

*eQTLGen is whole blood at n=31,684, FDR<0.05. A hit is strong causal evidence. **Absence is weak** — trans effects are small and detecting them needs enormous power, so a null overwhelmingly means underpowered rather than unregulated, and a relationship operating outside blood would not appear at all.*

### Routes a promoter analysis cannot see

Everything above acts on transcription. These act after it, around it, or on a different gene entirely. For some genes that is where the whole problem is: if the transcript is destroyed, the allele is silenced, or the protein is degraded too fast, more transcription does not help and may not be possible.

#### Could a relative do the job instead?

Whether a paralog could substitute for this gene is unanswered. It is worth asking explicitly, because the question 'restore this function' has more answers than 'raise this transcript'.

- CYP24A1 is not in the curated redundancy table, which is short and certainly incomplete. Absence here is not evidence that no paralog substitutes — it means nobody has written one down in this file.

#### After transcription

A promoter-level analysis stops at transcription. These act after it, and for some genes that is where the whole problem is: if the transcript is destroyed or mistranslated, doubling transcription doubles the amount that is destroyed. 1 of 6 mechanism(s) here have an FDA-approved agent. The rest are listed so that an empty transcriptional result is not mistaken for 'nothing can be done' — they are mechanisms, not treatments.

| Mechanism | Raises | Status | Applies when |
|---|---|---|---|
| splice correction | functional protein | approved | a splice-site or splice-regulatory variant, a poison exon, |
| transcript stabilisation | transcript | investigational | the transcript is short-lived and the promoter is already  |
| microRNA de-repression | transcript | investigational | a specific miRNA is a documented repressor of your gene in |
| translational enhancement | protein | investigational | protein is low while transcript is normal — the signature  |
| nonsense readthrough | functional protein | approved elsewhere | the gene carries a premature termination codon (PTC) |
| nonsense-mediated decay (NMD) inhibition | transcript | investigational | the gene carries a PTC and the transcript is being degrade |

**Splice correction** — redirects splicing so a productive transcript is made — exon inclusion, exon skipping, or suppression of a cryptic site.

- Exon-skipping restores a partly functional, internally deleted protein, not a wild-type one — the DMD agents are judged on dystrophin that is shorter than normal.
- Accelerated approvals in this class rest on a surrogate endpoint (protein level), and the confirmatory functional benefit has been contested.

#### Raising the protein instead of the transcript

- These leave transcription alone. If your protein is made normally and degraded too fast, every transcriptional route in this report is aimed at the wrong step - and nothing in this package can tell the two apart.
- The experiment that distinguishes them: compare transcript and protein. Normal mRNA with low protein points here; low mRNA points at everything above.

| Mechanism | Specificity | Status | Why |
|---|---|---|---|
| proteasome inhibition | 0.03 | approved | essentially none: the proteasome degrades most of the proteome |
| lysosomal / autophagy blockade | 0.05 | approved | non-specific, and the two degradation routes are distinct - wh |
| pharmacological chaperone | 0.85 | approved | genuinely high: a chaperone binds one protein, and this is the |
| E3 ligase inhibition | 0.7 | clinical | potentially high - an E3 has a limited substrate set - but it  |

*Three of these four raise most of the proteome, which makes them excellent diagnostics and poor therapies. The pharmacological chaperone is the exception — and it needs a chaperone to exist for your protein.*

### Candidates that may beat the current standard

**The bar: strength 0.498, specificity 0.454, efficiency 0.665** (predicted) — the best established value on *each axis*, not one drug's scores. Representative incumbent: Genistein (via TF).

1 of 3 candidates already have a reported link to this gene (CTD returned 180 chemicals with curated expression effects).

*Accounting for all 3 candidates: 1 already linked, 2 novel and non-dominated, 0 novel but beaten on every axis by another candidate, 0 with no advantage. All of them are in results.json under `improvements.all_candidates`.*

#### Compounds already reported to change this gene

These come from CTD's curated chemical–gene interactions: someone measured an expression change and a curator recorded it. That makes them the most concrete evidence in this report — but the direction is what decides whether each is a lead or a warning.

180 chemicals total: **85** move it the way you want, 48 move it the opposite way, 47 unclear.

**Expect the 'predicted by this run?' column to say no, and do not read that as a failed match.** CTD indexes MeSH, which is a vocabulary for the toxicology and pharmacology literature rather than a drug list, so what comes back is three different kinds of thing:

- **chemical** (153) — a compound. The only kind that can match a candidate.
- **exposure** (20) — formaldehyde, asbestos, bisphenol, tobacco smoke. Real measured effects on your gene, and not things anyone can be given.
- **category** (7) — MeSH headings such as "Anti-HIV Agents" or "Plant Extracts". Not compounds at all, so the column is blank rather than no.

Across all three directions, **1 of 180 matched a candidate from this run**. A low number here is the normal result: this report's drug tables are approved pharmacology, and most of CTD is not. It means the two sources are looking at different parts of chemical space, not that either is wrong.

**Reported to move CYP24A1 UP** — leads:

| Chemical | Kind | Predicted by this run? | CTD actions |
|---|---|---|---|
| Acetaminophen | chemical | no | affects^expression; increases^expression |
| Beclomethasone | chemical | no | increases^expression |
| Calcifediol | chemical | no | affects^expression; increases^expression |
| Calcitriol | chemical | no | affects^expression; increases^expression |
| Cholecalciferol | chemical | no | affects^expression; increases^expression |
| Copper Sulfate | chemical | no | increases^expression |
| Dexamethasone | chemical | no | affects^expression; increases^expression |
| Hydrocortisone | chemical | no | increases^expression |
| Ketoconazole | chemical | no | affects^expression; increases^expression |
| Leflunomide | chemical | no | increases^expression |
| Paricalcitol | chemical | no | increases^expression |
| Phenobarbital | chemical | no | increases^expression |
| Rifampin | chemical | no | increases^expression |
| Sorafenib | chemical | no | increases^expression |
| Testosterone | chemical | no | affects^expression; increases^expression |
| Tretinoin | chemical | no | affects^expression; increases^expression |
| Valproic Acid | chemical | no | affects^expression; increases^expression |
| Zinc Sulfate | chemical | no | increases^expression |

**22 compounds have a measured effect on this gene and no mechanism from this analysis.** Exposures and MeSH headings are excluded from that count, because the scan was never going to propose formaldehyde and listing it as a miss says nothing.

Most of what remains is research chemistry rather than medicine: tool compounds, probes and natural-product fractions that are measured because somebody had them on a shelf. Read them as evidence the gene IS movable, and check the few that are approved drugs individually.

Acetaminophen, Beclomethasone, Calcifediol, Calcitriol, Cholecalciferol, Copper Sulfate, Cyclosporine, Dexamethasone, Estradiol, Ethinyl Estradiol, Hydrocortisone, Ketoconazole, Leflunomide, Mirdametinib, Paricalcitol, Phenobarbital, Rifampin, Sorafenib, Testosterone, Tretinoin, Valproic Acid, Zinc Sulfate

**Reported to move it the WRONG way** — avoid, or use as negative controls:

Ampicillin, Arsenic Trioxide, Ascorbic Acid, Cisplatin, Copper, Flurbiprofen, Gemcitabine, Niclosamide, Progesterone, Raloxifene Hydrochloride, Tamoxifen, Teriflunomide

*CTD aggregates across cell types, doses and assays. A 'both' or contradictory annotation usually means the effect is context-dependent rather than that the data are wrong — check the underlying study before relying on any single row.*

#### No unlinked candidate beat the standard of care

Every candidate that scored well is already a reported regulator of this gene, or failed to clear the per-axis bar by more than the prediction error. This is the common result and it is informative: it says the obvious pharmacology has been explored.

### Drug–gene links with no recorded direction

75 interaction(s) across 24 factor(s) come from curated sources that record a drug–gene relationship **without saying which way it goes**. They are held back from the recommendations above, because the entire question here is direction and an undirected link cannot answer it.

They are still leads — someone curated a real relationship — but they need a literature check, not a prediction. Many are chemotherapeutics or broad stressors whose association with a transcription factor reflects cell-wide perturbation rather than a specific mechanism.

| TF | Drug | Reported by | DGIdb score |
|---|---|---|---|
| ELK1 | Ganciclovir | NCI | 11.71 |
| STAT5A | Tofacitinib | DTC | 2.34 |
| CTCF | Capivasertib | PRISM | 2.04 |
| KLF5 | Secnidazole | PRISM | 1.95 |
| STAT5B | Tofacitinib | DTC | 1.56 |
| ETV1 | Trametinib Dimethyl Sulfoxide | CIViC | 1.51 |
| ETV5 | Trametinib Dimethyl Sulfoxide | CIViC | 1.51 |
| TFAP4 | Azacitidine | PRISM | 1.09 |
| STAT6 | Esomeprazole Sodium | ClinPGx | 0.99 |
| JUND | Clozapine | NCI | 0.94 |
| JUN | Copper Chloride | DTC | 0.9 |
| MEF2A | Metformin | ClinPGx | 0.85 |
| NHLH1 | Cisplatin | ClinPGx | 0.55 |
| JUNB | Brigatinib | PRISM | 0.53 |
| FOXO1 | Epirubicin | ClinPGx | 0.52 |
| FOS | Bromocriptine Mesylate | NCI | 0.47 |
| E2F4 | Crizotinib | PRISM | 0.46 |
| JUN | Atomoxetine Hydrochloride | DTC | 0.45 |
| WT1 | Deferoxamine | NCI | 0.45 |
| WT1 | Dimethyl Sulfoxide | NCI | 0.42 |
| HNF4A | Urea | NCI | 0.39 |
| STAT5B | Dasatinib Anhydrous | Drug Repurposing Hub; TEND; TdgClinicalTrial | 0.34 |
| FOS | Pilocarpine Hydrochloride | NCI | 0.33 |
| HNF4A | Toremifene Citrate | DTC | 0.33 |
| FOS | Baclofen | NCI | 0.31 |
| … 50 more | | | |

*Pass `--include-unknown-direction` to fold these into the main tables anyway — but every such row will read UNKNOWN / unclear, which is what they are.*

### Approved drugs acting on these factors

No approved drug with activating pharmacology was found for any TF matched in this promoter. That is the usual result outside the nuclear-receptor family; see the README for indirect strategies.

---

### How to read this

A PWM hit means the sequence *resembles* a binding site. It is a hypothesis generator, not evidence of binding, and not evidence that binding would increase transcription of this gene in your cell type. Rank order here is a heuristic (weights printed in the JSON output), not a validated predictor. Every row needs wet-lab confirmation: promoter-reporter assay for the element, agonist dose-response with RT-qPCR/protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it.

---

# Testing any of this in the lab

Everything above is prediction. This part is experiment — how you would find out whether any of it holds, and it is placed last because it is long, not because it is optional.

**If you run one thing, run the reporter assay with the scrambled-motif control.** A drug that induces the wild-type promoter but not the mutant has an effect running through the element this report predicted. A drug that induces both means the prediction was wrong. That single control turns a hit into a mechanism.

## Designing a direct modality for CYP24A1

**Delivery to unknown: unknown** — —

specify a tissue: deliverability decides feasibility more than design quality does

**Delivery to unknown is unknown. Design quality will be the limiting factor.**

### splice-switching oligonucleotide

the main way to RAISE output of a specific gene: correct a splicing defect or force inclusion of an exon, as nusinersen does for SMN2

*requires an exploitable splicing event — a pseudoexon, a weak splice site, a skipped exon. Check the transcript structure before assuming one exists.*

### CRISPRa (dCas9-VP64/SAM)

forces transcription from the endogenous promoter; the cleanest test that raising this gene helps

*see the screen design section for guides*

### gene addition (AAV / lentiviral)

supply a working copy outright; the approved route for monogenic loss of function

*AAV cargo limit ~4.7 kb; immunogenicity and pre-existing neutralising antibodies exclude a substantial fraction of patients*

*Seed-region (guide positions 2-8) complementarity drives most siRNA off-target silencing and is a transcriptome-wide question. Run every candidate through BLAST against RefSeq, or a seed-match scan, before synthesising anything. The checks here are local sanity only.*

## Screen design

*One plan for the gene. The reporter insert below uses the ENST00000395954 window; clone a second construct for the other TSS if you want to compare them.*

The outcome chain above is almost always limited by *gene responsiveness in your cell type* — the link this package cannot see. A screen measures it instead of assuming it. Everything above prioritises this experiment; it does not replace it.

### Reporter construct

- Insert: 20:54170954-54176953 (5000 bp upstream / 1000 bp downstream of the ENST00000395954 TSS)
- Vector: promoterless luciferase (pGL4.10-style) or GFP for FACS-based screening; luciferase for arrayed compound screens, GFP for pooled CRISPRa
- Controls:
  - Promoterless vector: establishes the floor. Without it you cannot tell activation from vector background.
  - Constitutive promoter (CMV/EF1a): establishes the ceiling and reports transfection or transduction efficiency per well.
  - Scrambled-motif variant of your top predicted element: the single most informative control here. If a drug induces the wild-type reporter but not the mutant, the effect runs through that element. If it induces both, the mechanism is elsewhere and the prediction was wrong.
  - Dual-reporter normalisation (Renilla or similar) to separate real induction from changes in cell number or general translation — cytotoxic compounds otherwise read as repressors.
- An episomal reporter lacks the chromatin context of the endogenous locus. A hit should be confirmed at the native gene by RT-qPCR before it is believed.

### CRISPRa screen

- Approach: dCas9-VP64 with SAM (MS2-p65-HSF1) or dCas9-VPR for stronger activation
- Focused library: 20 TFs prioritised by this run (SOX2, POU2F1, FOS, JUNB, JUND, FOSL2, JUN, KLF5, STAT3, NFYB...)
- Unbiased: ~1,600 human TFs (Calabrese or a TF-focused sublibrary). Run this if the focused set comes back empty — an unbiased screen is the only thing that finds a factor nobody has associated with your promoter.
- Coverage: 3-5 sgRNAs per gene, >=500 cells per guide at transduction, MOI ~0.3 to keep single integrations
- Readout: FACS sort GFP-high vs GFP-low, then sequence guide abundance; or arrayed luciferase for a focused set

Positive-control sgRNAs at your promoter (12 clean of 12 shown):

| Protospacer (20 nt) | PAM | Strand | vs TSS | GC | Note |
|---|---|---|---|---|---|
| `CGGTGCTTCACAGGCTGTAT` | GG | + | -288 | 55% | ok |
| `GGCTGTATGGTCTGCTATGG` | GG | + | -276 | 55% | ok |
| `GTATGGTCTGCTATGGGGCA` | GG | + | -272 | 55% | ok |
| `CAGGCTGTATGGTCTGCTAT` | GG | + | -278 | 50% | ok |
| `AGGCTGTATGGTCTGCTATG` | GG | + | -277 | 50% | ok |
| `GATCTGAGTTGCTGGTGATG` | GG | + | -87 | 50% | ok |
| `ATCTGAGTTGCTGGTGATGG` | GG | + | -86 | 50% | ok |
| `ACAAGTGCACGATCGCAAAA` | GG | + | -245 | 45% | ok |
| `AGATCTGAGTTGCTGGTGAT` | GG | + | -88 | 45% | ok |
| `CAATTCCTAGCGGTAAAAGG` | GG | - | -51 | 45% | ok |

*Guides listed target YOUR promoter directly, for a positive-control arm proving the reporter can be activated at all. The library guides target the TFs. Run these through CRISPOR or CHOPCHOP for off-target scoring before ordering — no off-target analysis is done here.*

### Compound screen

- Prestwick Chemical Library (~1,520 approved drugs, high chemical and pharmacological diversity)
- Selleck FDA-approved Drug Library (~3,000)
- NCATS Pharmaceutical Collection (NPC) — free, and the most complete approved set

- Design: Arrayed, 384-well, 8-point dose-response on the actives. Single-concentration primary screens at 10 uM produce a cytotoxicity-dominated hit list; dose-response separates real induction from cells dying.
- Counter-screen: Run the same plate against a constitutive-promoter reporter. Anything that moves both is affecting transcription or translation globally, not your promoter, and that is the majority of primary hits.

*A full-library screen is unbiased and will find mechanisms this package cannot predict — that is its advantage, not a criticism of it. Use the cherry-pick list to get an answer in a week; run the full library when you want the answer you did not expect.*

### What to expect

- **compound screen hit rate**: 0.5-2% at a sensible threshold. Most primary hits are cytotoxic or assay-interfering; expect to lose 60-80% at counter-screen.
- **crispra hit rate**: typically a handful of TFs from ~1,600. A screen returning hundreds of hits is reporting a global effect on the reporter, not specific regulation.
- **replicates**: 3 biological replicates minimum; 2 gives you no estimate of variance and a Z' you cannot trust
- **quality metric**: Z' >= 0.5 between positive and negative controls before you screen anything. Below that the screen cannot distinguish a real hit from noise and running it wastes the library.
- **what a negative result means**: If no TF in a genome-wide CRISPRa screen activates your reporter, the promoter fragment is probably not sufficient — the regulation likely lives in an enhancer outside your window. Widen it (--upstream 50000) or use the endogenous locus.
