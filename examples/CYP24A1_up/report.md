# Upregulation strategy for CYP24A1 — promoter / TF / drug analysis

Objective: **increase** transcription of CYP24A1.

> **Only FDA approved drugs are shown.** 798 entries are left out. 262 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 42 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

## What would have to happen

**CYP24A1 dosage: tolerant.** LOEUF 1.51: the population carries about as many loss-of-function variants as chance predicts, so losing ONE copy is tolerated. That is not evidence the gene is unimportant. A recessive condition is invisible to this measure entirely - carriers are healthy, so selection never saw the phenotype. Mosaic variegated aneuploidy works this way: biallelic BUB1B, with unaffected heterozygous parents. One copy of this gene is dispensable in the population, so if your condition is recessive the relevant question is how much function the variant alleles retain - not what the population tolerates.

*Read from gnomAD constraint v4.1; ClinGen dosage curation 07 Oct,2026.*

## Evidence available in this run

**⚠ DEGRADED** — Read the ranking as a reading order, not as evidence. At least one pillar of the analysis is missing. Measured evidence that DID arrive: chip, lincs, ctd, genetics.

- the motif scan is 'no signal' — about 86% of the hit list is expected to be chance (range 46% to 100%)

| Source | State | What it returned |
|---|---|---|
| chip | ok | local BED: peaks for 145 factor(s) over the window |
| ctd | ok | 180 chemical(s) with curated expression effects on this gene; 85 move it the way you want |
| depmap | ok | 14 of 14 factor(s) scored in CRISPRGeneEffect.csv |
| direction_model | ok, nothing found | 0 of 1 grade(s) got a model call; too few labelled agonists and antagonists to train |
| expression | ok | 16 tissue(s) for CYP24A1 |
| factor_abundance | ok | 24 compound(s) move a factor's own transcript the needed way; 0 factor(s) skipped for an unrecorded role |
| genetics | ok | 13 associated phenotype(s) (Open Targets); 206 genome-wide significant GWAS Catalog association(s), mapped or nearest gene, which is weaker than fine-mapped support |
| knockout | ok | 5 of 14 factor(s) were perturbed in the knockout data; the gene moved for 2 |
| lincs | ok | 3 signature hit(s) for this gene |
| product_screen | ok | 1 unannotated candidate(s) at the CYP24A1 protein, cross-validated AUC 1.0 |
| promoter_dna | ok | 3 sequence feature(s) enriched over 200 shuffles; 9 agent(s) prefer them |
| regulators | ok | 19 recorded regulator(s) of the gene had no drug route yet; 2 have an approved drug that moves them the needed way (13 route(s)) |
| regulon | ok | 590 TFs with targets called from ChIP occupancy. A peak near a gene is binding, not regulation |
| screen | ok | no factor has an approved drug that names its target, so there was no receptor to screen |
| sensitivity | ok | 2 correlated compound(s) across 6 of 6 factor(s) scanned |
| trans_eqtl | ok, nothing found | 0 of 6 factor(s) have trans-eQTL support |
| upstream | ok | druggable routes for 4 factor(s); 2059 of 2974 reachable upstream node(s) checked across 8 factor(s) |

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
| 1 | **vorinostat** | histone deacetylase inhibitor | cutaneous T-cell lymphoma | NFIB | — | 3 | — | 54.6 |
| 2 | **vemurafenib** | raises SOX2 transcript in 5 profile(s) | not established here | SOX2 | — | 3 | — | 38.6 |
| 3 | **Vorinostat** | histone deacetylase inhibitor | cutaneous T-cell lymphoma | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 4 | **Romidepsin** | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 5 | **Belinostat** | HDAC inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 6 | **Panobinostat** | HDAC inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 7 | **Valproic acid** | weak HDAC inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 8 | **Azacitidine** | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 9 | **Decitabine** | DNA methyltransferase inhibitor | myelodysplastic syndrome | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 10 | **Tazemetostat** | EZH2 inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |

**Every row is off-label.** The drug was approved for what the fourth column says, and is listed here for something else. That column is the baseline this proposal departs from, and the distance between it and your gene is part of what you are judging - a topical acne gel and an infused cytotoxic can sit on adjacent rows with similar scores.

*Ordered by evidence tier first, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. **BBB?** is not part of the ordering — it is an applicability filter you apply yourself, since whether brain exposure is wanted or unwanted depends on where your gene needs to move. See the note under the graded table for what a yes does and does not establish.*

### Why each one

Every route below is **evidence tier 3** — a human curator annotated this mechanism.

**The direction rule**, which applies to every group below: you want to raise the gene, so you agonise an activator and block a repressor.

#### Through **NFIB**

- Role: **activator**
- Binding: 1 predicted site(s), the best scoring 0.846 -935 bp from the TSS; **no ChIP peak over this promoter**, although ChIP data was available — the motif is unconfirmed.

1. **vorinostat**
   - ⚠ measured against NFIB's own transcript, not against your gene. That NFIB controls your gene is the inference this route rests on

#### Through **SOX2**

- Role: **activator**
- Binding: 2 predicted site(s), the best scoring 0.911 +134 bp from the TSS; **1 ChIP peak(s) confirm it actually binds here** — the strongest link in this chain.

2. **vemurafenib**
   - ⚠ measured against SOX2's own transcript, not against your gene. That SOX2 controls your gene is the inference this route rests on

#### Through **chromatin (genome-wide)** — 8 candidates

- The factor has no activator/repressor annotation, so the direction of its effect on your gene is assumed, not known.

- ⚠ acts genome-wide: will move thousands of promoters including yours. Useful as a positive control that your reporter CAN be induced; almost never an acceptable therapeutic route on specificity grounds

3. **Vorinostat** — HDAC inhibitor
4. **Romidepsin** — HDAC inhibitor
5. **Belinostat** — HDAC inhibitor
6. **Panobinostat** — HDAC inhibitor
7. **Valproic acid** — weak HDAC inhibitor
8. **Azacitidine** — DNMT inhibitor
9. **Decitabine** — DNMT inhibitor
10. **Tazemetostat** — EZH2 inhibitor

*Not medical advice. "Approved" means approved for its own indication — using one of these to move an unrelated gene is off-label repurposing with the full systemic pharmacology attached. The experimental sections at the end of this report are how you would actually test any of it.*

## Promoter of CYP24A1 (ENST00000395954)

- Transcript chosen: **distinct TSS**
- Assembly: from GTF  |  Strand: -
- TSS: chr20:54171953
- Window: chr20:54171754-54172953 (-1000 / +200, 1200 bp)
- Motifs scanned: 879 at p ≤ 1e-05 per position
- 321 of those motifs were not scanned: in this window their best possible match is more common than p = 1e-05, so any hit would be a chance match at a higher rate than stated. A looser `--pvalue` scans them, with more chance hits.
- Non-redundant sites: 13 across 14 TFs

> **⚠ This scan found no signal.** 13 sites is what chance alone predicts (11 over 851 scannable bases). This scan found nothing. Tighten to --pvalue 1e-6, or narrow the window (--upstream 1000 --downstream 200), or both — and use --chip-bed so ranking rests on measured binding rather than sequence.
>
> Everything below is ranked from a list of chance matches. Fix the scan before reading any of it.

*about 86% of these 13 sites (range 46% to 100%) would be expected on shuffled sequence of this window, so roughly 2 are above chance; among the best-scoring sites the share is lower: 77% of the top 5, 62% of the top 10. Estimated from 30 shuffled copies that keep the repeat mask and every dinucleotide count. The formula that assumes each motif matches at exactly its stated p-value would expect 9 and is not used. Nothing here says which sites are real; only ChIP evidence separates them.*

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

**14 factor(s) had predicted sites here. 0 of them have an approved drug (0 of those are nuclear receptors).**

Roughly 48 of the ~1,600 human transcription factors have a ligand-binding pocket, and the steroid, retinoid, thyroid and vitamin-D families are most of them. So the drug tables are not a ranking that hormone drugs won — they are what is left after every factor without a ligand has been dropped.

**14 of the top 25 factors by score have no drug and appear in no drug table:** SOX2, POU5F1, POU2F1, ZNF530, EBF2, EBF4, MTF1, NFIB, ZNF707, ZNF320, ZNF410, KLF9 ….

That is where the gene-specific signal is. A half-site like AGGTCA occurs every few kilobases by chance, so nuclear receptors surface for almost any promoter; these factors surfaced for *this* one. The routes that can reach them are the upstream-network, stabilisation and protein-level sections, not the ligand tables.

### Top-ranked transcription factors

| Rank | TF | Role | NR | Other genes it controls | Cells need it? | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | SOX2 | activator |  | ~1200 | a few lines | 2 | 0.911 | +134 | 1 | 0 | 79.7 |
| 2 | POU5F1 | repressor |  | ~1200 | a few lines | 1 | 0.911 | +134 | 1 | 0 | 65.9 |
| 3 | POU2F1 | activator |  | not counted | some lines | 1 | 0.882 | +133 | 0 | 0 | 55.6 |
| 4 | ZNF530 | repressor |  | 253 | no | 1 | 0.922 | -16 | 0 | 0 | 49.9 |
| 5 | EBF2 | activator |  | not counted | a few lines | 1 | 0.954 | -868 | 0 | 0 | 48.7 |
| 6 | EBF4 | activator |  | not counted | no | 1 | 0.898 | -869 | 0 | 0 | 47.0 |
| 7 | MTF1 | activator |  | not counted | some lines | 1 | 0.848 | -757 | 0 | 0 | 46.0 |
| 8 | NFIB | activator |  | not counted | a few lines | 1 | 0.846 | -935 | 0 | 0 | 45.2 |
| 9 | ZNF707 | unknown |  | 348 | no | 1 | 0.836 | -15 | 0 | 0 | 42.8 |
| 10 | ZNF320 | unknown |  | 2148 | a few lines | 1 | 0.835 | -29 | 0 | 0 | 42.3 |
| 11 | ZNF410 | unknown |  | 676 | a few lines | 1 | 0.89 | -429 | 0 | 0 | 35.2 |
| 12 | KLF9 | unknown |  | not counted | a few lines | 1 | 0.937 | -716 | 0 | 0 | 34.0 |
| 13 | EBF3 | unknown |  | not counted | a few lines | 1 | 0.958 | -868 | 0 | 0 | 33.8 |
| 14 | ZNF816 | unknown |  | 590 | a few lines | 1 | 0.869 | -503 | 0 | 0 | 33.7 |

**Other genes it controls** is the collateral question. A drug acting through a factor moves everything that factor regulates, so this column says how wide the effect will be. Read it against the score: a strong route through a factor that controls two thousand genes is a worse starting point than a weaker route through one that controls thirty.

A count written `~N` is a bundled estimate for a factor the regulon file does not list (2 here), not a count of curated edges. Treat it as an order of magnitude.

`not counted` means the regulon table has no entry and there is no bundled estimate either (7 here). It is not a small number. Most of the roughly 1,600 human transcription factors have no curated regulon at all.

Counts come from literature curation, so they measure what has been published as much as what the factor does. SP1 has thousands of recorded targets partly because it regulates many genes and partly because thousands of people assayed it.

### What a cell loses without each factor (DepMap)

The regulon column counts how many genes a factor controls, from literature curation, which measures how much anyone has published as much as what the factor does. This is the same collateral question answered by measurement: CRISPR knockout of each factor across roughly 1,100 cancer cell lines.

It also reaches factors the curated tables cannot. A factor nobody has written a regulon for still has a knockout phenotype here.

14 of 14 factor(s) were scored in DepMap. 11 (EBF2, EBF3, KLF9, MTF1, NFIB) show a dependency confined to a minority of lines, which is the informative case: the requirement has a context, and whether your tissue shares it is the question to ask next. 3 (EBF4, ZNF530, ZNF707) were dispensable throughout. That bounds the proliferative cost of removing them and says nothing about a differentiated tissue.

| Factor | Cells need it? | Mean effect | Lines dependent | Most dependent lineage |
|---|---|---|---|---|
| POU2F1 | essential in a minority of lines | -0.34 | 146/1150 (13%) | Kidney (-0.44) |
| MTF1 | essential in a minority of lines | -0.21 | 110/1150 (10%) | Myeloid (-0.43) |
| POU5F1 | essential in a handful of lines | -0.11 | 17/1150 (2%) | Biliary Tract (-0.17) |
| KLF9 | essential in a handful of lines | -0.08 | 2/1150 (0%) | CNS/Brain (-0.11) |
| ZNF320 | essential in a handful of lines | -0.07 | 5/1150 (0%) | Eye (-0.15) |
| SOX2 | essential in a handful of lines | -0.06 | 53/1150 (5%) | Esophagus/Stomach (-0.18) |
| ZNF816 | essential in a handful of lines | -0.04 | 3/1150 (0%) | Bone (-0.10) |
| ZNF707 | dispensable across the panel | -0.03 | 0/1150 (0%) | — |
| ZNF410 | essential in a handful of lines | -0.02 | 1/1150 (0%) | Myeloid (-0.06) |
| EBF2 | essential in a handful of lines | -0.02 | 2/1150 (0%) | Pancreas (-0.05) |
| NFIB | essential in a handful of lines | -0.01 | 16/1150 (1%) | Prostate (-0.10) |
| ZNF530 | dispensable across the panel | +0.02 | 0/1150 (0%) | — |
| EBF4 | dispensable across the panel | +0.03 | 0/1150 (0%) | — |
| EBF3 | essential in a handful of lines | +0.06 | 4/1150 (0%) | Peripheral Nervous System (-0.06) |

**Mean effect** is on DepMap's own scale: 0 means knockout changed nothing, and -1 is the median of known common-essential genes. A line counts as dependent below -0.5.

#### What this does not mean

- **These are cancer cell lines dividing in a dish.** Essentiality here is a statement about proliferating transformed cells. A factor whose job is in a post-mitotic or differentiated cell looks dispensable for reasons that have nothing to do with its importance in your tissue.
- **Knockout is not inhibition.** The protein is removed entirely, including scaffolding roles no ligand would touch. Essential genes are drugged successfully at partial occupancy all the time, so a dependency is a reason for care rather than a veto.
- **Nothing here is about your gene.** This bounds the cost of acting on a factor. Whether acting on it moves your gene is the rest of the report.
- **A dispensable result is not permission.** It bounds the proliferative cost and nothing else.

**The selective ones are the interesting ones.** EBF2, EBF3, KLF9, MTF1, NFIB, POU2F1 matter in a minority of lines rather than all of them, so the requirement has a context. The lineage column is where to look for what that context is, and whether your tissue shares it is the question this panel cannot answer.

### Drugs that behave like losing the factor

Every other route in this report starts from something somebody wrote down - a curated ligand, a motif, a structure, a signed edge. None of them can surface a drug nobody has connected to the factor. This one reads a measured phenotype instead: if a compound hurts the same cell lines that knocking the factor out hurts, it behaves as though it removes the factor.

| Factor | Drug | Partial r | Raw r | Lines | q | Behaves as |
|---|---|---|---|---|---|---|
| ZNF530 | **GANCICLOVIR** | 0.2307 | 0.2179 | 460 | 0.003978 | phenocopies |

**Partial r** controls for how easily each cell line dies under the whole panel. Some lines die under almost anything and those same lines suffer most knockouts, which puts a correlation between almost any drug and almost any gene. The control moves a number either way - it removes a shared component that can inflate a correlation or mask one - so where the two columns differ a lot the raw value was mostly about the cell lines. Estimates use at least 25 shared lines.

**This is not a target.** A correlated drug may act on the factor, on anything else the same lines need, or on a pathway the factor sits in, and the first is the least likely of the three for any one hit. Nothing here records any of them binding anything: the finding is a reason to run a binding assay, not a result from one.

**A proliferation readout.** Both panels score growth of transformed cells in a dish, so a factor whose role is in a differentiated or post-mitotic cell looks silent here for reasons that have nothing to do with its importance.

### Every route that acts on this gene

106 routes across 7 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

**Why established binders sit at tier 2, just under the direct agents.** A drug already measured to bind one of your factors has half its chain settled: the drug-to-receptor step has assay data behind it, so the pharmacology is de-risked and only the receptor-to-gene step is a hypothesis. A curated tier-3 ligand is annotation on *both* halves. But they stay below tier 1, because a direct agent acts on your gene itself and skips the receptor-to-gene inference entirely — a stronger claim however well-characterised the binder is. The cost of the indirect route is in the row: everything else that factor regulates moves with your gene.

- 10 × curated — mechanism annotated by a human curator
- 92 × predicted — a model says so
- 4 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 3 | TF: abundance change | vorinostat | NFIB | 1 | raises NFIB transcript in 123 profile(s) | 1.0 | 0.45 | 0.5 | 54.6 |
| 3 | TF: abundance change | vemurafenib | SOX2 | 1 | raises SOX2 transcript in 5 profile(s) | 0.5 | 0.45 | 0.5 | 38.6 |
| 3 | chromatin-level | Vorinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Romidepsin | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Belinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Panobinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Valproic acid | chromatin (genome-wide) | 3 | weak HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Azacitidine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Decitabine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Tazemetostat | chromatin (genome-wide) | 3 | EZH2 inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 4 | TF: known regulator | Cholecalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Paricalcitol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Doxercalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Ergocalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Calcifediol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Calcitriol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Bexarotene | RXRA | 1 | activate RXRA, a recorded regulator of the ge | 0.5 | 0.5 | 0.6 | 42.0 |
| 4 | TF: known regulator | Acitretin | RXRA | 1 | activate RXRA, a recorded regulator of the ge | 0.5 | 0.5 | 0.6 | 42.0 |
| 4 | TF: known regulator | Etretinate | RXRA | 1 | activate RXRA, a recorded regulator of the ge | 0.5 | 0.5 | 0.6 | 42.0 |
| 4 | upstream signalling | Acitretin | RARG | 3 | activate RARG → SOX2 (2 hops) | 0.562 | 0.556 | 0.6 | 32.5 |
| 4 | upstream signalling | Palovarotene | RARG | 3 | activate RARG → SOX2 (2 hops) | 0.562 | 0.556 | 0.6 | 32.5 |
| 4 | upstream signalling | Tretinoin | RARG | 3 | activate RARG → SOX2 (2 hops) | 0.562 | 0.556 | 0.6 | 32.5 |
| 4 | upstream signalling | Acitretin | RARB | 3 | activate RARB → SOX2 (2 hops) | 0.562 | 0.543 | 0.6 | 32.3 |
| 4 | upstream signalling | Tretinoin | RARB | 3 | activate RARB → SOX2 (2 hops) | 0.562 | 0.543 | 0.6 | 32.3 |
| 4 | upstream signalling | Ramelteon | MTNR1A | 3 | activate MTNR1A → POU5F1 (2 hops) | 0.562 | 0.532 | 0.6 | 32.1 |
| 4 | upstream signalling | Tasimelteon | MTNR1A | 3 | activate MTNR1A → POU5F1 (2 hops) | 0.562 | 0.532 | 0.6 | 32.1 |
| 4 | upstream signalling | Misoprostol | PTGER3 | 3 | activate PTGER3 → SOX2 (2 hops) | 0.562 | 0.51 | 0.6 | 31.6 |
| 4 | upstream signalling | Acitretin | RARA | 3 | activate RARA → SOX2 (2 hops) | 0.562 | 0.455 | 0.6 | 30.5 |
| 4 | upstream signalling | Tretinoin | RARA | 3 | activate RARA → SOX2 (2 hops) | 0.562 | 0.455 | 0.6 | 30.5 |
| 4 | upstream signalling | Bexarotene | RXRA | 3 | activate RXRA → POU5F1 (2 hops) | 0.562 | 0.403 | 0.6 | 29.5 |

Best in each category:

- **TF: abundance change** — vorinostat (tier 3, score 54.6)
- **chromatin-level** — Vorinostat (tier 3, score 23.4)
- **TF: known regulator** — Cholecalciferol (tier 4, score 45.2)
- **upstream signalling** — Acitretin (tier 4, score 32.5)
- **measured phenocopy** — GANCICLOVIR (tier 4, score 28.6)
- **direct: predicted product binder** — PARICALCITOL (tier 5, score 62.3)
- **promoter DNA binder** — Mitoxantrone (tier 5, score 31.2)

*No route found in: direct: gene or transcript, direct: gene product, TF: orthosteric ligand, TF: inferred activity change, TF: predicted pocket binder, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, TF: degradation blocker, off-pathway pocket binder, empirical signature, virtual screen hypothesis. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

| Category | What it changes | Your drugs in it | Assay that sees it |
|---|---|---|---|
| TF: abundance change | **abundance** | vorinostat, vemurafenib | RT-qPCR and a western for the FACTOR first, then RT-qPCR for your gene |
| chromatin-level | **abundance, genome-wide** | Vorinostat, Romidepsin, Belinostat, Panobinostat … +4 | RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened |
| TF: known regulator | **abundance** | Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol … +5 | RT-qPCR for your gene after treatment; check the factor's own activity (a reporter or a known target) moved first |
| upstream signalling | **abundance** | Acitretin, Palovarotene, Tretinoin, Ramelteon … +52 | RT-qPCR then western/ELISA; confirm the pathway node moved too, or you cannot attribute the effect |
| measured phenocopy | **abundance** | GANCICLOVIR | a binding assay against the factor, then a reporter, then RT-qPCR for your gene - the correlation names a candidate and establishes none of those steps |
| direct: predicted product binder | **activity** | PARICALCITOL | a binding assay against the protein first, then a functional assay for direction - RT-qPCR will not see this route at all, and a western blot reads it as no effect |
| promoter DNA binder | **abundance** | Mitoxantrone, Doxorubicin, Pentamidine | an EMSA or footprint on this promoter first, then a reporter, then RT-qPCR - and a counter-screen on an unrelated promoter carrying the same feature, because selectivity is the whole question here |

*Every category present here works through transcription, so all of them move your gene's product **abundance** and an RT-qPCR is the right first readout. The one exception in the taxonomy — a drug on the gene product itself, which changes activity only — did not appear in this run.*

**TF: abundance change** (2 found)

- *Yours:* vorinostat, vemurafenib
- A compound measured to change how much of the factor's own transcript exists, rather than what the factor does.
- *Why:* The only route here that needs the factor to be expressed rather than druggable, so it reaches factors with no ligand, no structure and no curated upstream kinase - which is most of them.
- *Cost:* More factor is not more factor ACTIVITY. For a ligand-activated receptor the two are nearly unrelated, and transcript is not protein. Broad agents that shut down transcription move every factor and are flagged.
- *Differs:* Unlike 'empirical signature', which is a compound measured against YOUR GENE, this was measured against the FACTOR - one further inferential step from what you care about.

**chromatin-level** (8 found)

- *Yours:* Vorinostat, Romidepsin, Belinostat, Panobinostat, Valproic acid, Azacitidine, Decitabine, Tazemetostat
- HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.
- *Why:* Reliable positive control that your reporter *can* be induced.
- *Cost:* Specificity ~0.05. These move thousands of promoters. Almost never an acceptable therapeutic route on that basis alone.
- *Differs:* Not targeted at your gene in any sense; it is in the table so you can see what a non-specific option scores.

**TF: known regulator** (9 found)

- *Yours:* Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol, Calcifediol, Calcitriol, Bexarotene, Acitretin … +1 more
- Drugs that act on a factor CollecTRI or perturbation data records as raising or lowering your gene, whether or not the motif scan ranked that factor.
- *Why:* Reaches factors the motif scan misses, such as those that bind through a partner or a distant enhancer.
- *Cost:* The link comes from a database or a knockout in other cells, so it may not hold in yours. The drug's effect on the factor is taken from its recorded action.
- *Differs:* Starts from the literature on your gene, not from its promoter sequence. CollecTRI and CTD share sources, so agreement between them is not independent.

**upstream signalling** (82 found)

- *Yours:* Acitretin, Palovarotene, Tretinoin, Ramelteon, Tasimelteon, Misoprostol, Bexarotene, Etretinate … +48 more
- Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- *Why:* Turns 'this factor has no drugs' into a real list. Most TFs are only reachable this way.
- *Cost:* Specificity falls with every step: a node three hops up reaches your gene through one branch and thousands of other genes through the rest.
- *Differs:* The only category where the drug never touches your gene or its TF — read the printed path and judge it yourself.

**measured phenocopy** (1 found)

- *Yours:* GANCICLOVIR
- A compound whose effect on cell viability across the DepMap panel tracks what happens when this factor is knocked out, so it behaves as though it removes the factor.
- *Why:* It needs no annotation of any kind. Every other route here starts from something a curator, a motif or a structure already asserted, so none of them can surface a drug nobody has connected to the factor. This one can, because it reads a measured phenotype.
- *Cost:* Correlation across cell lines, not target engagement. A hit may act on the factor, on anything else the same lines need, or on a pathway the factor sits in, and the first is the least likely of the three. Both panels score proliferation of transformed cells, so a factor that matters in a differentiated cell leaves no trace.
- *Differs:* Unlike 'empirical signature', nothing was measured about your gene or about transcription at all - the readout is whether cells live.

**direct: predicted product binder** (1 found)

- *Yours:* PARICALCITOL
- A compound predicted to bind this gene's own protein and not recorded as acting on it - either because its chemistry resembles the protein's known ligands, or because it docks into the protein's pocket.
- *Why:* It is the only route here aimed at the gene's own product rather than at a factor upstream of it, so it skips every inference about whether moving a factor moves your gene. And unlike the curated product drugs it can propose a compound no one has annotated.
- *Cost:* Predicted binding and nothing more. The model is refused below a cross-validated AUC floor, which bounds how badly it can be ranking, but a hit is a hypothesis for a binding assay. Direction is separate again: whether binding inhibits or activates the protein decides whether it serves your goal.
- *Differs:* Unlike 'direct: gene product', nothing was measured or curated for this compound at this protein - that category is annotation, this is a model. Unlike 'TF: predicted pocket binder', the pocket is on the gene's own product, so no inference about a factor controlling the gene is involved.

**promoter DNA binder** (3 found)

- *Yours:* Mitoxantrone, Doxorubicin, Pentamidine
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
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **virtual screen hypothesis** — An approved drug predicted to bind one of your TFs, never tested against it.

### What these drugs are

Every route in this report is off-label by construction: the drug was approved for one thing and is being considered for another. The approved use is the baseline that departure is measured from, so it is listed here once per drug rather than repeated in each table above.

| Drug | What it is | Approved to treat |
|---|---|---|
| Acitretin | retinoid | severe psoriasis |
| Azacitidine | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia |
| Bexarotene | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma |
| Calcifediol | 25-hydroxyvitamin D | vitamin D deficiency |
| Calcitriol | active vitamin D, VDR agonist | hypocalcaemia and secondary hyperparathyroidism in renal failure |
| Cholecalciferol | vitamin D3, a precursor requiring activation | vitamin D deficiency and supplementation |
| Decitabine | DNA methyltransferase inhibitor | myelodysplastic syndrome |
| Dextrothyroxine Sodium | T4 stereoisomer | formerly hyperlipidaemia - withdrawn for cardiac mortality |
| Diethylstilbestrol Diphosphate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Doxercalciferol | vitamin D prodrug | secondary hyperparathyroidism in dialysis |
| Doxorubicin | anthracycline topoisomerase II inhibitor | many solid and haematological cancers |
| Elafibranor | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) |
| Ergocalciferol | vitamin D receptor agonist class | secondary hyperparathyroidism in renal disease, osteoporosis, and plaque psoriasis topically |
| Estradiol Acetate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Estradiol Cypionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Etretinate | retinoid | severe psoriasis - WITHDRAWN, teratogenic |
| Lenalidomide Hydrochloride Monohydrate | cereblon modulator | multiple myeloma, MDS |
| Levothyroxine Sodium | synthetic T4 | hypothyroidism |
| Liothyronine Sodium | synthetic T3 | hypothyroidism, myxoedema coma |
| Mitoxantrone | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis |
| Nandrolone Phenpropionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Oxandrolone | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Palovarotene | retinoid class | acne, psoriasis and other keratinisation disorders topically; acute promyelocytic leukaemia and cutaneous lymphoma systemically |
| Paricalcitol | selective VDR activator | secondary hyperparathyroidism in chronic kidney disease |
| Raloxifene Hydrochloride | selective oestrogen receptor modulator | osteoporosis; breast cancer risk reduction |
| Romidepsin | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma |
| Seladelpar | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) |
| Testosterone Enanthate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Testosterone Propionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Tretinoin | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically |
| vorinostat | histone deacetylase inhibitor | cutaneous T-cell lymphoma |

**No description established for 34**: Acetaminophen, Avatrombopag, Belinostat, Chenodiol, Chlorotrianisene, Cinacalcet Hydrochloride, Dronabinol, Eltrombopag Olamine, Estramustine Phosphate Sodium, Etrasimod, Fingolimod Hydrochloride, Fingolimod Lauryl Sulfate, GANCICLOVIR, Lusutrombopag, Midodrine Hydrochloride, Misoprostol, Nabilone, Obeticholic Acid, Oxymetholone, Ozanimod Hydrochloride, Panobinostat, Pentamidine, Phenylephrine Hydrochloride, Phenylpropanolamine Hydrochloride, Ponesimod…. That means this package could not find one offline, not that the compound has no approved use. Several are environmental agents and research chemicals from the CTD tables, which have no indication because they are not drugs.

Descriptions come from: 11 class, 20 curated, 34 unknown. None of these is a label. The curated text is a summary written for orientation, ChEMBL's indication list includes uses studied but never approved, and a class line gives the kind of drug rather than its indication. Check the prescribing information before acting on any row.

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

**POU5F1** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| MTNR1A | 2 | 7 | MTNR1A --> MAPK1 —\| POU5F1 | inhibit | Ramelteon, Tasimelteon | 30.3 |
| NOD2 | 2 | 9 | NOD2 --> MAPK1 —\| POU5F1 | inhibit | — | 29.2 |
| RXRA | 2 | 22 | RXRA —\| CTNNB1 --> POU5F1 | inhibit | Bexarotene, Acitretin, Etretinate | 22.8 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK3 —\| POU5F1 | inhibit | Acetaminophen | 18.7 |
| NR1H4 | 3 | 1 | NR1H4 --> RXRA —\| CTNNB1 --> POU5F1 | inhibit | Ursodiol, Chenodiol, Obeticholic Acid | 18.7 |
| PPARD | 3 | 4 | PPARD --> RXRA —\| CTNNB1 --> POU5F1 | inhibit | Seladelpar, Elafibranor, Seladelpar Lysine | 17.4 |
| THRB | 3 | 6 | THRB --> RXRA —\| CTNNB1 --> POU5F1 | inhibit | Liothyronine Sodium, Liothyronine, Levothyroxine Sodium | 16.7 |
| RBX1 | 3 | 7 | RBX1 --> CUL1 —\| CTNNB1 --> POU5F1 | inhibit | Lenalidomide Hydrochloride Monohydrate | 16.3 |
| THRA | 3 | 8 | THRA --> RXRA —\| CTNNB1 --> POU5F1 | inhibit | Liothyronine Sodium, Liothyronine, Levothyroxine Sodium | 16.0 |
| VDR | 3 | 9 | VDR --> RXRA —\| CTNNB1 --> POU5F1 | inhibit | Cholecalciferol, Paricalcitol, Doxercalciferol | 15.7 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 4 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 6 downstream branches

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

**MTF1** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| P2RY4 | 3 | 6 | P2RY4 --> GNA13 —\| LATS1 —\| MTF1 | activate | — | 16.7 |
| S1PR2 | 3 | 8 | S1PR2 --> GNA13 —\| LATS1 —\| MTF1 | activate | Fingolimod Lauryl Sulfate, Fingolimod Hydrochloride | 16.1 |
| P2RY2 | 3 | 9 | P2RY2 --> GNA13 —\| LATS1 —\| MTF1 | activate | — | 15.7 |
| ADRA1A | 3 | 10 | ADRA1A --> GNA13 —\| LATS1 —\| MTF1 | activate | Phenylephrine Hydrochloride, Midodrine Hydrochloride, Phenylpropanolamine Hydrochloride | 15.4 |
| P2RY6 | 3 | 10 | P2RY6 --> GNA13 —\| LATS1 —\| MTF1 | activate | — | 15.4 |
| CNR1 | 3 | 10 | CNR1 --> GNA13 —\| LATS1 —\| MTF1 | activate | Nabilone, Dronabinol | 15.4 |
| P2RY1 | 3 | 11 | P2RY1 --> GNA13 —\| LATS1 —\| MTF1 | activate | — | 15.1 |
| HTR6 | 3 | 11 | HTR6 --> GNA13 —\| LATS1 —\| MTF1 | activate | — | 15.1 |
| S1PR3 | 3 | 12 | S1PR3 --> GNA13 —\| LATS1 —\| MTF1 | activate | Fingolimod Lauryl Sulfate, Fingolimod Hydrochloride | 14.8 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 10 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 11 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 12 downstream branches

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

#### No stabilisation route for 8 factor(s)

EBF2, EBF4, MTF1, NFIB, POU2F1, POU5F1, SOX2, ZNF530

These have no ligand pocket, no direct drug, and **nothing curated against them or their degradation machinery**. For these factors stabilisation is a research programme, not a route you can take: making it real needs the E3 ligase that degrades the factor *and* a ligand for that ligase, and neither is published here.

The mechanism class is real — it has shipped approved drugs, listed once below — but those act on **other proteins entirely**, and none of them is evidence that any factor above is tractable.

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Cell-type reach

Three things must overlap: where the drug goes, where the gene is expressed, and where the transcription factor is. Any empty intersection kills the route regardless of motif quality or potency.

| TF | Co-expressed with gene | Shared tissues |
|---|---|---|
| SOX2 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| POU5F1 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| POU2F1 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| ZNF530 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| EBF2 | yes (11) | breast, cervix, uterine, endometrium, esophagus, lung, prost |
| EBF4 | yes (14) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| MTF1 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| NFIB | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| ZNF707 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| ZNF320 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |

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

*No factor tested (SOX2, POU5F1, POU2F1, ZNF530, EBF2, EBF4) has trans-eQTL support for this gene.*

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

9 interaction(s) across 2 factor(s) come from curated sources that record a drug–gene relationship **without saying which way it goes**. They are held back from the recommendations above, because the entire question here is direction and an undirected link cannot answer it.

They are still leads — someone curated a real relationship — but they need a literature check, not a prediction. Many are chemotherapeutics or broad stressors whose association with a transcription factor reflects cell-wide perturbation rather than a specific mechanism.

| TF | Drug | Reported by | DGIdb score |
|---|---|---|---|
| NFIB | Paliperidone | ClinPGx | 1.11 |
| POU5F1 | Allopurinol | ClinPGx | 0.73 |
| POU5F1 | Phenytoin Sodium | ClinPGx | 0.56 |
| POU5F1 | Carbamazepine | ClinPGx | 0.48 |
| NFIB | Clozapine | ClinPGx | 0.39 |
| NFIB | Methylphenidate | ClinPGx | 0.37 |
| POU5F1 | Dehydrated Alcohol | NCI | 0.3 |
| NFIB | Risperidone | ClinPGx | 0.22 |
| NFIB | Vincristine | ClinPGx | 0.13 |

*Pass `--include-unknown-direction` to fold these into the main tables anyway — but every such row will read UNKNOWN / unclear, which is what they are.*

### Approved drugs acting on these factors

No approved drug with activating pharmacology was found for any TF matched in this promoter. That is the usual result outside the nuclear-receptor family; see the README for indirect strategies.

## Promoter of CYP24A1 (ENST00000216862)

- Transcript chosen: **distinct TSS (representing 2 within 500 bp)**
- Assembly: from GTF  |  Strand: -
- TSS: chr20:54173986
- Window: chr20:54173787-54174986 (-1000 / +200, 1200 bp)
- Motifs scanned: 879 at p ≤ 1e-05 per position
- 243 of those motifs were not scanned: in this window their best possible match is more common than p = 1e-05, so any hit would be a chance match at a higher rate than stated. A looser `--pvalue` scans them, with more chance hits.
- Non-redundant sites: 16 across 17 TFs

> **⚠ This scan found no signal.** 16 sites is what chance alone predicts (14 over 1077 scannable bases). This scan found nothing. Tighten to --pvalue 1e-6, or narrow the window (--upstream 1000 --downstream 200), or both — and use --chip-bed so ranking rests on measured binding rather than sequence.
>
> Everything below is ranked from a list of chance matches. Fix the scan before reading any of it.

*about 89% of these 16 sites (range 50% to 100%) would be expected on shuffled sequence of this window, so roughly 2 are above chance; among the best-scoring sites the share is lower: 100% of the top 5, 100% of the top 10. Estimated from 30 shuffled copies that keep the repeat mask and every dinucleotide count. The formula that assumes each motif matches at exactly its stated p-value would expect 14 and is not used. Nothing here says which sites are real; only ChIP evidence separates them.*

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

**17 factor(s) had predicted sites here. 1 of them have an approved drug (0 of those are nuclear receptors).**

Roughly 48 of the ~1,600 human transcription factors have a ligand-binding pocket, and the steroid, retinoid, thyroid and vitamin-D families are most of them. So the drug tables are not a ranking that hormone drugs won — they are what is left after every factor without a ligand has been dropped.

**16 of the top 25 factors by score have no drug and appear in no drug table:** NFYB, CEBPD, CEBPA, TBX2, ELK1, GLIS3, FOXO1, TBX3, STAT1, GLIS1, PRDM4, PAX6 ….

That is where the gene-specific signal is. A half-site like AGGTCA occurs every few kilobases by chance, so nuclear receptors surface for almost any promoter; these factors surfaced for *this* one. The routes that can reach them are the upstream-network, stabilisation and protein-level sections, not the ligand tables.

### Top-ranked transcription factors

| Rank | TF | Role | NR | Other genes it controls | Cells need it? | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | NFYB | activator |  | ~5000 | yes, most | 1 | 1.0 | -60 | 1 | 0 | 81.2 |
| 2 | CEBPD | activator |  | not counted | a few lines | 2 | 1.0 | -886 | 1 | 0 | 73.3 |
| 3 | CEBPA | activator |  | not counted | a few lines | 1 | 0.938 | -372 | 2 | 0 | 72.4 |
| 4 | TBX2 | repressor |  | not counted | a few lines | 1 | 0.969 | +53 | 1 | 0 | 70.0 |
| 5 | ELK1 | activator |  | ~5000 | a few lines | 1 | 0.875 | -652 | 1 | 0 | 67.5 |
| 6 | GLIS3 | dual |  | not counted | no | 1 | 0.877 | -260 | 1 | 0 | 66.5 |
| 7 | FOXO1 | dual |  | ~1200 | a few lines | 1 | 0.875 | -652 | 1 | 0 | 61.5 |
| 8 | NFE2L2 | activator |  | ~1200 | some lines | 1 | 0.895 | -656 | 0 | 1 | 58.1 |
| 9 | TBX3 | dual |  | 3222 | a few lines | 1 | 0.957 | +53 | 0 | 0 | 54.2 |
| 10 | STAT1 | unknown |  | ~1200 | no | 1 | 0.869 | -544 | 2 | 0 | 53.3 |
| 11 | GLIS1 | activator |  | not counted | no | 1 | 0.851 | -260 | 0 | 0 | 51.7 |
| 12 | PRDM4 | activator |  | not counted | a few lines | 1 | 0.907 | -913 | 0 | 0 | 47.1 |
| 13 | PAX6 | dual |  | 1892 | a few lines | 1 | 0.853 | -817 | 0 | 0 | 39.9 |
| 14 | STAT2 | unknown |  | 1536 | a few lines | 1 | 0.869 | -544 | 0 | 0 | 33.3 |
| 15 | IRF7 | unknown |  | not counted | a few lines | 1 | 0.859 | -544 | 0 | 0 | 33.0 |
| 16 | ZNF677 | unknown |  | 823 | no | 1 | 0.894 | -867 | 0 | 0 | 31.9 |
| 17 | PRDM9 | unknown |  | not counted | a few lines | 1 | 0.873 | -848 | 0 | 0 | 31.3 |

**Other genes it controls** is the collateral question. A drug acting through a factor moves everything that factor regulates, so this column says how wide the effect will be. Read it against the score: a strong route through a factor that controls two thousand genes is a worse starting point than a weaker route through one that controls thirty.

A count written `~N` is a bundled estimate for a factor the regulon file does not list (5 here), not a count of curated edges. Treat it as an order of magnitude.

`not counted` means the regulon table has no entry and there is no bundled estimate either (8 here). It is not a small number. Most of the roughly 1,600 human transcription factors have no curated regulon at all.

Counts come from literature curation, so they measure what has been published as much as what the factor does. SP1 has thousands of recorded targets partly because it regulates many genes and partly because thousands of people assayed it.

### What a cell loses without each factor (DepMap)

The regulon column counts how many genes a factor controls, from literature curation, which measures how much anyone has published as much as what the factor does. This is the same collateral question answered by measurement: CRISPR knockout of each factor across roughly 1,100 cancer cell lines.

It also reaches factors the curated tables cannot. A factor nobody has written a regulon for still has a knockout phenotype here.

17 of 17 factor(s) were scored in DepMap. 12 (CEBPA, CEBPD, ELK1, FOXO1, IRF7) show a dependency confined to a minority of lines, which is the informative case: the requirement has a context, and whether your tissue shares it is the question to ask next. 4 (GLIS1, GLIS3, STAT1, ZNF677) were dispensable throughout. That bounds the proliferative cost of removing them and says nothing about a differentiated tissue.

| Factor | Cells need it? | Mean effect | Lines dependent | Most dependent lineage |
|---|---|---|---|---|
| NFYB | essential in most lines | -0.70 | 961/1150 (84%) | Thyroid (-0.84) |
| NFE2L2 | essential in a minority of lines | -0.24 | 109/1150 (10%) | Pleura (-0.38) |
| ELK1 | essential in a handful of lines | -0.17 | 3/1150 (0%) | Ovary/Fallopian Tube (-0.21) |
| TBX3 | essential in a handful of lines | -0.16 | 50/1150 (4%) | Prostate (-0.53) |
| CEBPD | essential in a handful of lines | -0.15 | 16/1150 (1%) | Prostate (-0.25) |
| FOXO1 | essential in a handful of lines | -0.13 | 29/1150 (2%) | Lymphoid (-0.30) |
| STAT2 | essential in a handful of lines | -0.12 | 8/1150 (1%) | Uterus (-0.17) |
| CEBPA | essential in a handful of lines | -0.07 | 22/1150 (2%) | Myeloid (-0.36) |
| TBX2 | essential in a handful of lines | -0.06 | 17/1150 (2%) | Skin (-0.15) |
| STAT1 | dispensable across the panel | -0.01 | 0/1150 (0%) | — |
| PRDM4 | essential in a handful of lines | -0.00 | 1/1150 (0%) | Prostate (-0.11) |
| IRF7 | essential in a handful of lines | +0.01 | 1/1150 (0%) | Peripheral Nervous System (-0.03) |
| GLIS1 | dispensable across the panel | +0.03 | 0/1150 (0%) | — |
| PAX6 | essential in a handful of lines | +0.05 | 1/1150 (0%) | Eye (+0.01) |
| ZNF677 | dispensable across the panel | +0.08 | 0/1150 (0%) | — |
| PRDM9 | essential in a handful of lines | +0.16 | 1/1150 (0%) | Soft Tissue (+0.12) |
| GLIS3 | dispensable across the panel | +0.17 | 0/1150 (0%) | — |

**Mean effect** is on DepMap's own scale: 0 means knockout changed nothing, and -1 is the median of known common-essential genes. A line counts as dependent below -0.5.

#### What this does not mean

- **These are cancer cell lines dividing in a dish.** Essentiality here is a statement about proliferating transformed cells. A factor whose job is in a post-mitotic or differentiated cell looks dispensable for reasons that have nothing to do with its importance in your tissue.
- **Knockout is not inhibition.** The protein is removed entirely, including scaffolding roles no ligand would touch. Essential genes are drugged successfully at partial occupancy all the time, so a dependency is a reason for care rather than a veto.
- **Nothing here is about your gene.** This bounds the cost of acting on a factor. Whether acting on it moves your gene is the rest of the report.
- **A dispensable result is not permission.** It bounds the proliferative cost and nothing else.

**The selective ones are the interesting ones.** CEBPA, CEBPD, ELK1, FOXO1, IRF7, NFE2L2 matter in a minority of lines rather than all of them, so the requirement has a context. The lineage column is where to look for what that context is, and whether your tissue shares it is the question this panel cannot answer.

### Drugs that behave like losing the factor

Every other route in this report starts from something somebody wrote down - a curated ligand, a motif, a structure, a signed edge. None of them can surface a drug nobody has connected to the factor. This one reads a measured phenotype instead: if a compound hurts the same cell lines that knocking the factor out hurts, it behaves as though it removes the factor.

| Factor | Drug | Partial r | Raw r | Lines | q | Behaves as |
|---|---|---|---|---|---|---|
| CEBPD | **ABEMACICLIB** | -0.2171 | -0.1749 | 462 | 0.017208 | opposes |

**Partial r** controls for how easily each cell line dies under the whole panel. Some lines die under almost anything and those same lines suffer most knockouts, which puts a correlation between almost any drug and almost any gene. The control moves a number either way - it removes a shared component that can inflate a correlation or mask one - so where the two columns differ a lot the raw value was mostly about the cell lines. Estimates use at least 25 shared lines.

**This is not a target.** A correlated drug may act on the factor, on anything else the same lines need, or on a pathway the factor sits in, and the first is the least likely of the three for any one hit. Nothing here records any of them binding anything: the finding is a reason to run a binding assay, not a result from one.

**A proliferation readout.** Both panels score growth of transformed cells in a dish, so a factor whose role is in a differentiated or post-mitotic cell looks silent here for reasons that have nothing to do with its importance.

### Every route that acts on this gene

128 routes across 9 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

**Why established binders sit at tier 2, just under the direct agents.** A drug already measured to bind one of your factors has half its chain settled: the drug-to-receptor step has assay data behind it, so the pharmacology is de-risked and only the receptor-to-gene step is a hypothesis. A curated tier-3 ligand is annotation on *both* halves. But they stay below tier 1, because a direct agent acts on your gene itself and skips the receptor-to-gene inference entirely — a stronger claim however well-characterised the binder is. The cost of the indirect route is in the row: everything else that factor regulates moves with your gene.

- 19 × curated — mechanism annotated by a human curator
- 105 × predicted — a model says so
- 4 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 3 | TF: degradation blocker | Dimethyl fumarate | KEAP1 | 1 | degradation blocker | 0.75 | 0.65 | 0.7 | 56.2 |
| 3 | TF: degradation blocker | Omaveloxolone | KEAP1 | 1 | degradation blocker | 0.75 | 0.65 | 0.7 | 56.2 |
| 3 | TF: degradation blocker | Monomethyl fumarate | KEAP1 | 1 | degradation blocker | 0.75 | 0.65 | 0.7 | 56.2 |
| 3 | TF: orthosteric ligand | Omaveloxolone | NFE2L2 | 1 | ligand at the transcription factor | 0.595 | 0.42 | 0.725 | 45.3 |
| 3 | TF: abundance change | vorinostat | ELK1 | 1 | raises ELK1 transcript in 139 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | dasatinib | CEBPD | 1 | raises CEBPD transcript in 51 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | mitoxantrone | CEBPA | 1 | raises CEBPA transcript in 10 profile(s) | 1.0 | 0.1 | 0.5 | 44.8 |
| 3 | TF: abundance change | sirolimus | TBX2 | 1 | lowers TBX2 transcript in 3 profile(s) | 0.3 | 0.45 | 0.5 | 32.2 |
| 3 | TF: abundance change | dasatinib | TBX2 | 1 | lowers TBX2 transcript in 6 profile(s) | 0.6 | 0.1 | 0.5 | 32.0 |
| 3 | TF: abundance change | vorinostat | NFE2L2 | 1 | raises NFE2L2 transcript in 4 profile(s) | 0.4 | 0.1 | 0.5 | 25.6 |
| 3 | TF: abundance change | selumetinib | TBX2 | 1 | lowers TBX2 transcript in 4 profile(s) | 0.4 | 0.1 | 0.5 | 25.6 |
| 3 | chromatin-level | Vorinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Romidepsin | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Belinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Panobinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Valproic acid | chromatin (genome-wide) | 3 | weak HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Azacitidine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Decitabine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Tazemetostat | chromatin (genome-wide) | 3 | EZH2 inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 4 | upstream signalling | Progesterone | PGR | 2 | activate PGR → FOXO1 (1 hop) | 0.75 | 0.976 | 0.6 | 52.8 |
| 4 | upstream signalling | Norethindrone | PGR | 2 | activate PGR → FOXO1 (1 hop) | 0.75 | 0.976 | 0.6 | 52.8 |
| 4 | upstream signalling | Ethynodiol Diacetate | PGR | 2 | activate PGR → FOXO1 (1 hop) | 0.75 | 0.976 | 0.6 | 52.8 |
| 4 | upstream signalling | Hydroxyprogesterone Caproate | PGR | 2 | activate PGR → FOXO1 (1 hop) | 0.75 | 0.976 | 0.6 | 52.8 |
| 4 | TF: known regulator | Cholecalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Paricalcitol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Doxercalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Ergocalciferol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Calcifediol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Calcitriol | VDR | 1 | activate VDR, a recorded regulator of the gen | 0.6 | 0.5 | 0.6 | 45.2 |
| 4 | TF: known regulator | Bexarotene | RXRA | 1 | activate RXRA, a recorded regulator of the ge | 0.5 | 0.5 | 0.6 | 42.0 |

Best in each category:

- **TF: degradation blocker** — Dimethyl fumarate (tier 3, score 56.2)
- **TF: orthosteric ligand** — Omaveloxolone (tier 3, score 45.3)
- **TF: abundance change** — vorinostat (tier 3, score 44.8)
- **chromatin-level** — Vorinostat (tier 3, score 23.4)
- **upstream signalling** — Progesterone (tier 4, score 52.8)
- **TF: known regulator** — Cholecalciferol (tier 4, score 45.2)
- **measured phenocopy** — ABEMACICLIB (tier 4, score 28.1)
- **direct: predicted product binder** — Paricalcitol (tier 5, score 63.8)
- **promoter DNA binder** — Plicamycin (tier 5, score 31.2)

*No route found in: direct: gene or transcript, direct: gene product, TF: inferred activity change, TF: predicted pocket binder, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, off-pathway pocket binder, empirical signature, virtual screen hypothesis. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

| Category | What it changes | Your drugs in it | Assay that sees it |
|---|---|---|---|
| TF: degradation blocker | **abundance** | Dimethyl fumarate, Omaveloxolone, Monomethyl fumarate | RT-qPCR for your gene; western for the FACTOR to confirm it accumulated |
| TF: orthosteric ligand | **abundance** | Omaveloxolone | RT-qPCR then western/ELISA |
| TF: abundance change | **abundance** | vorinostat, dasatinib, mitoxantrone, sirolimus … +1 | RT-qPCR and a western for the FACTOR first, then RT-qPCR for your gene |
| chromatin-level | **abundance, genome-wide** | Vorinostat, Romidepsin, Belinostat, Panobinostat … +4 | RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened |
| upstream signalling | **abundance** | Progesterone, Norethindrone, Ethynodiol Diacetate, Hydroxyprogesterone Caproate … +56 | RT-qPCR then western/ELISA; confirm the pathway node moved too, or you cannot attribute the effect |
| TF: known regulator | **abundance** | Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol … +5 | RT-qPCR for your gene after treatment; check the factor's own activity (a reporter or a known target) moved first |
| measured phenocopy | **abundance** | ABEMACICLIB | a binding assay against the factor, then a reporter, then RT-qPCR for your gene - the correlation names a candidate and establishes none of those steps |
| direct: predicted product binder | **activity** | Paricalcitol | a binding assay against the protein first, then a functional assay for direction - RT-qPCR will not see this route at all, and a western blot reads it as no effect |
| promoter DNA binder | **abundance** | Plicamycin, Mitoxantrone, Doxorubicin | an EMSA or footprint on this promoter first, then a reporter, then RT-qPCR - and a counter-screen on an unrelated promoter carrying the same feature, because selectivity is the whole question here |

*Every category present here works through transcription, so all of them move your gene's product **abundance** and an RT-qPCR is the right first readout. The one exception in the taxonomy — a drug on the gene product itself, which changes activity only — did not appear in this run.*

**TF: degradation blocker** (3 found)

- *Yours:* Dimethyl fumarate, Omaveloxolone, Monomethyl fumarate
- Stops the factor being destroyed, so more of it accumulates.
- *Why:* The class that actually shipped: NRF2 via KEAP1, HIF via PHD inhibitors.
- *Cost:* Effect is capped by the synthesis rate, and stabilising a factor stabilises it everywhere.
- *Differs:* Raises protein *abundance*, visible on a western blot — unlike a conformational stabiliser, which does not.

**TF: orthosteric ligand** (1 found)

- *Yours:* Omaveloxolone
- Binds the transcription factor's natural ligand pocket.
- *Why:* The best-characterised route, and where nearly all approved TF pharmacology lives.
- *Cost:* Only ~48 human TFs (the nuclear receptors) have such a pocket, and the factor's whole regulon moves with your gene.
- *Differs:* Clamps the receptor on or off regardless of endogenous signalling — contrast the allosteric and ligand-axis routes, which scale what is already happening.

**TF: abundance change** (7 found)

- *Yours:* vorinostat, dasatinib, mitoxantrone, sirolimus, selumetinib
- A compound measured to change how much of the factor's own transcript exists, rather than what the factor does.
- *Why:* The only route here that needs the factor to be expressed rather than druggable, so it reaches factors with no ligand, no structure and no curated upstream kinase - which is most of them.
- *Cost:* More factor is not more factor ACTIVITY. For a ligand-activated receptor the two are nearly unrelated, and transcript is not protein. Broad agents that shut down transcription move every factor and are flagged.
- *Differs:* Unlike 'empirical signature', which is a compound measured against YOUR GENE, this was measured against the FACTOR - one further inferential step from what you care about.

**chromatin-level** (8 found)

- *Yours:* Vorinostat, Romidepsin, Belinostat, Panobinostat, Valproic acid, Azacitidine, Decitabine, Tazemetostat
- HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.
- *Why:* Reliable positive control that your reporter *can* be induced.
- *Cost:* Specificity ~0.05. These move thousands of promoters. Almost never an acceptable therapeutic route on that basis alone.
- *Differs:* Not targeted at your gene in any sense; it is in the table so you can see what a non-specific option scores.

**upstream signalling** (95 found)

- *Yours:* Progesterone, Norethindrone, Ethynodiol Diacetate, Hydroxyprogesterone Caproate, Raloxifene Hydrochloride, Diethylstilbestrol Diphosphate, Estramustine Phosphate Sodium, Toremifene Citrate … +52 more
- Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- *Why:* Turns 'this factor has no drugs' into a real list. Most TFs are only reachable this way.
- *Cost:* Specificity falls with every step: a node three hops up reaches your gene through one branch and thousands of other genes through the rest.
- *Differs:* The only category where the drug never touches your gene or its TF — read the printed path and judge it yourself.

**TF: known regulator** (9 found)

- *Yours:* Cholecalciferol, Paricalcitol, Doxercalciferol, Ergocalciferol, Calcifediol, Calcitriol, Bexarotene, Acitretin … +1 more
- Drugs that act on a factor CollecTRI or perturbation data records as raising or lowering your gene, whether or not the motif scan ranked that factor.
- *Why:* Reaches factors the motif scan misses, such as those that bind through a partner or a distant enhancer.
- *Cost:* The link comes from a database or a knockout in other cells, so it may not hold in yours. The drug's effect on the factor is taken from its recorded action.
- *Differs:* Starts from the literature on your gene, not from its promoter sequence. CollecTRI and CTD share sources, so agreement between them is not independent.

**measured phenocopy** (1 found)

- *Yours:* ABEMACICLIB
- A compound whose effect on cell viability across the DepMap panel tracks what happens when this factor is knocked out, so it behaves as though it removes the factor.
- *Why:* It needs no annotation of any kind. Every other route here starts from something a curator, a motif or a structure already asserted, so none of them can surface a drug nobody has connected to the factor. This one can, because it reads a measured phenotype.
- *Cost:* Correlation across cell lines, not target engagement. A hit may act on the factor, on anything else the same lines need, or on a pathway the factor sits in, and the first is the least likely of the three. Both panels score proliferation of transformed cells, so a factor that matters in a differentiated cell leaves no trace.
- *Differs:* Unlike 'empirical signature', nothing was measured about your gene or about transcription at all - the readout is whether cells live.

**direct: predicted product binder** (1 found)

- *Yours:* Paricalcitol
- A compound predicted to bind this gene's own protein and not recorded as acting on it - either because its chemistry resembles the protein's known ligands, or because it docks into the protein's pocket.
- *Why:* It is the only route here aimed at the gene's own product rather than at a factor upstream of it, so it skips every inference about whether moving a factor moves your gene. And unlike the curated product drugs it can propose a compound no one has annotated.
- *Cost:* Predicted binding and nothing more. The model is refused below a cross-validated AUC floor, which bounds how badly it can be ranking, but a hit is a hypothesis for a binding assay. Direction is separate again: whether binding inhibits or activates the protein decides whether it serves your goal.
- *Differs:* Unlike 'direct: gene product', nothing was measured or curated for this compound at this protein - that category is annotation, this is a model. Unlike 'TF: predicted pocket binder', the pocket is on the gene's own product, so no inference about a factor controlling the gene is involved.

**promoter DNA binder** (3 found)

- *Yours:* Plicamycin, Mitoxantrone, Doxorubicin
- A small molecule with a documented preference for a sequence feature this promoter carries - a quadruplex, a GC-box, an AT tract - proposed because the feature is enriched against shuffles of this same sequence.
- *Why:* It is the only route that acts on the promoter itself rather than on a protein, so it needs no druggable factor and no annotation. Where a promoter depends on one feature, occupying it is a direct intervention.
- *Cost:* Sequence preference is not selectivity. These agents bind the feature wherever it occurs and most of the genome's copies are not in your gene, so the dose that engages your promoter engages thousands of others. Nearly every approved member is cytotoxic, which is how it works rather than a side effect. And a motif is not a folded structure in a cell.
- *Differs:* Unlike every TF route, no protein is targeted - which removes the question of whether the factor controls your gene and replaces it with whether the molecule can tell your promoter from any other.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **direct: gene product** — Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- **TF: inferred activity change** — A compound whose own LINCS signature moves this factor's known targets, found without any record that it binds the factor.
- **TF: predicted pocket binder** — A compound predicted by docking to fit the factor's ligand pocket, from a structure rather than from any measured activity.
- **TF: allosteric modulator** — Binds a site away from the ligand pocket, amplifying or damping the endogenous ligand's effect (PAM/NAM).
- **TF: stabiliser or glue** — Locks the factor in an active conformation, or glues it to a partner protein.
- **TF: endogenous ligand** — Changes how much natural ligand exists, by drugging the enzyme that makes or destroys it.
- **off-pathway pocket binder** — A docking run against a protein that is not one of the factors found at this promoter. Usually an upstream node someone chose to dock by hand.
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **virtual screen hypothesis** — An approved drug predicted to bind one of your TFs, never tested against it.

### What these drugs are

Every route in this report is off-label by construction: the drug was approved for one thing and is being considered for another. The approved use is the baseline that departure is measured from, so it is listed here once per drug rather than repeated in each table above.

| Drug | What it is | Approved to treat |
|---|---|---|
| Acitretin | retinoid | severe psoriasis |
| Azacitidine | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia |
| Bexarotene | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma |
| Calcifediol | 25-hydroxyvitamin D | vitamin D deficiency |
| Calcitriol | active vitamin D, VDR agonist | hypocalcaemia and secondary hyperparathyroidism in renal failure |
| Cholecalciferol | vitamin D3, a precursor requiring activation | vitamin D deficiency and supplementation |
| Decitabine | DNA methyltransferase inhibitor | myelodysplastic syndrome |
| Desoxycorticosterone Acetate | mineralocorticoid class | adrenal insufficiency and salt-wasting states |
| Desoxycorticosterone Pivalate | mineralocorticoid class | adrenal insufficiency and salt-wasting states |
| Diethylstilbestrol Diphosphate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Doxercalciferol | vitamin D prodrug | secondary hyperparathyroidism in dialysis |
| Doxorubicin | anthracycline topoisomerase II inhibitor | many solid and haematological cancers |
| Ergocalciferol | vitamin D receptor agonist class | secondary hyperparathyroidism in renal disease, osteoporosis, and plaque psoriasis topically |
| Estradiol Acetate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Estradiol Cypionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Etretinate | retinoid | severe psoriasis - WITHDRAWN, teratogenic |
| Fludrocortisone Acetate | mineralocorticoid class | adrenal insufficiency and salt-wasting states |
| Hydrocortisone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Lenalidomide Hydrochloride Monohydrate | cereblon modulator | multiple myeloma, MDS |
| Mometasone Furoate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Norethindrone | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Palovarotene | retinoid class | acne, psoriasis and other keratinisation disorders topically; acute promyelocytic leukaemia and cutaneous lymphoma systemically |
| Paricalcitol | selective VDR activator | secondary hyperparathyroidism in chronic kidney disease |
| Plicamycin | RNA synthesis inhibitor, binds GC-rich DNA | formerly hypercalcaemia of malignancy - withdrawn |
| Prednisolone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Progesterone | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Raloxifene Hydrochloride | selective oestrogen receptor modulator | osteoporosis; breast cancer risk reduction |
| Romidepsin | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma |
| Tretinoin | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically |
| mitoxantrone | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis |
| selumetinib | MEK1/2 inhibitor | neurofibromatosis type 1 plexiform neurofibroma |
| sirolimus | mTOR inhibitor | transplant rejection prophylaxis, lymphangioleiomyomatosis |
| vorinostat | histone deacetylase inhibitor | cutaneous T-cell lymphoma |

**No description established for 41**: ABEMACICLIB, Acetaminophen, Alprostadil, Avatrombopag, Belinostat, Bethanechol Chloride, Carboprost Tromethamine, Cevimeline Hydrochloride, Chlorotrianisene, Cinacalcet Hydrochloride, Dimethyl fumarate, Dinoprost Tromethamine, Eltrombopag Olamine, Estramustine Phosphate Sodium, Ethynodiol Diacetate, Etrasimod, Fenfluramine Hydrochloride, Hydroxyprogesterone Caproate, Iloprost, Lusutrombopag, Methylergonovine Maleate, Midodrine Hydrochloride, Misoprostol, Monomethyl fumarate, Omaveloxolone…. That means this package could not find one offline, not that the compound has no approved use. Several are environmental agents and research chemicals from the CTD tables, which have no indication because they are not drugs.

Descriptions come from: 13 class, 20 curated, 41 unknown. None of these is a label. The curated text is a summary written for orientation, ChEMBL's indication list includes uses studied but never approved, and a class line gives the kind of drug rather than its indication. Check the prescribing information before acting on any row.

### Graded candidates

Combined = geometric mean of TF confidence and drug quality. A weak TF prediction caps the grade no matter how good the drug.

| Grade | Combined | Drug | BBB? | TF | TF conf. | Drug qual. | Potency | Protein select. | Gene select. | Evidence | Flags |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **B** | 60.4 | Omaveloxolone | no | NFE2L2 | 0.581 | 0.66 | 0.5 | 1.0 | 0.42 | 0.5 | no potency data; scored neutral; thin evidence base (n=0 potency measu |

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

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

**CEBPD** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| CHRM3 | 3 | 5 | CHRM3 --> GNAQ —\| GSK3B —\| CEBPD | activate | Pilocarpine Hydrochloride, Bethanechol Chloride | 17.1 |
| HTR2B | 3 | 5 | HTR2B --> GNAQ —\| GSK3B —\| CEBPD | activate | Methylergonovine Maleate | 17.1 |
| HTR2A | 3 | 6 | HTR2A --> GNAQ —\| GSK3B —\| CEBPD | activate | — | 16.7 |
| CHRM1 | 3 | 6 | CHRM1 --> GNAQ —\| GSK3B —\| CEBPD | activate | Pilocarpine Hydrochloride, Xanomeline, Cevimeline Hydrochloride | 16.7 |
| MTNR1A | 3 | 7 | MTNR1A --> GNAQ —\| GSK3B —\| CEBPD | activate | Ramelteon, Tasimelteon | 16.4 |
| PTGFR | 3 | 7 | PTGFR --> GNAQ —\| GSK3B —\| CEBPD | activate | Dinoprost Tromethamine, Carboprost Tromethamine | 16.3 |
| PTGIR | 3 | 7 | PTGIR --> GNAQ —\| GSK3B —\| CEBPD | activate | Treprostinil, Treprostinil Diolamine, Selexipag | 16.3 |
| ADRA1D | 3 | 7 | ADRA1D --> GNAQ —\| GSK3B —\| CEBPD | activate | Phenylephrine Hydrochloride, Midodrine Hydrochloride | 16.3 |
| PTGER1 | 3 | 8 | PTGER1 --> GNAQ —\| GSK3B —\| CEBPD | activate | Alprostadil | 16.0 |
| HTR2C | 3 | 8 | HTR2C --> GNAQ —\| GSK3B —\| CEBPD | activate | Fenfluramine Hydrochloride | 16.0 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 5 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 6 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 7 downstream branches

**CEBPA** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| RBX1 | 3 | 7 | RBX1 —\| ESR1 --> NR0B2 —\| CEBPA | activate | Lenalidomide Hydrochloride Monohydrate | 16.4 |
| NOD2 | 3 | 9 | NOD2 --> MAPK14 —\| MAPK3 —\| CEBPA | activate | — | 15.8 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 7 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 9 downstream branches
- ⚠ path passes through hub(s): MAPK14, MAPK3

**ELK1** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| ESR2 | 2 | 4 | ESR2 --> MAPK11 --> ELK1 | activate | Raloxifene Hydrochloride, Diethylstilbestrol Diphosphate, Estramustine Phosphate Sodium | 34.8 |
| MTNR1A | 2 | 7 | MTNR1A --> MAPK1 --> ELK1 | activate | Ramelteon, Tasimelteon | 30.2 |
| NOD2 | 2 | 9 | NOD2 --> MAPK1 --> ELK1 | activate | — | 29.7 |
| TRPV1 | 3 | 1 | TRPV1 --> PRKCB --> MAPK1 --> ELK1 | activate | Acetaminophen | 18.7 |
| CASR | 3 | 2 | CASR --> FYN --> MAPK14 --> ELK1 | activate | Cinacalcet Hydrochloride | 18.3 |
| RARG | 3 | 5 | RARG --> AKT1 --> MAPK1 --> ELK1 | activate | Acitretin, Palovarotene, Tretinoin | 17.1 |
| RARB | 3 | 6 | RARB --> AKT1 --> MAPK1 --> ELK1 | activate | Acitretin, Tretinoin | 16.7 |
| PTGIR | 3 | 7 | PTGIR --> STAT3 --> MAPK8 --> ELK1 | activate | Treprostinil, Treprostinil Diolamine, Selexipag | 16.3 |
| RBX1 | 3 | 7 | RBX1 --> NFKB1 --> CREBBP --> ELK1 | activate | Lenalidomide Hydrochloride Monohydrate | 16.3 |
| PTGER3 | 3 | 9 | PTGER3 --> SRC --> MAPK1 --> ELK1 | activate | Misoprostol | 15.7 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 5 downstream branches

**GLIS3** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| CASR | 3 | 2 | CASR --> FYN —\| ITCH —\| GLIS3 | activate | Cinacalcet Hydrochloride | 18.3 |
| S1PR1 | 3 | 12 | S1PR1 --> FYN —\| ITCH —\| GLIS3 | activate | Ponesimod, Siponimod, Etrasimod | 14.8 |
| NR3C1 | 3 | 16 | NR3C1 —\| MAPK8 --> ITCH —\| GLIS3 | activate | Hydrocortisone Acetate, Prednisolone Acetate, Mometasone Furoate | 13.7 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 12 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 16 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches

**FOXO1** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| PGR | 1 | 1 | PGR --> FOXO1 | activate | Progesterone, Norethindrone, Ethynodiol Diacetate | 98.5 |
| RXRA | 2 | 22 | RXRA --> NR1I2 --> FOXO1 | activate | Bexarotene, Acitretin, Etretinate | 33.6 |
| RBX1 | 2 | 7 | RBX1 —\| ESR1 —\| FOXO1 | activate | Lenalidomide Hydrochloride Monohydrate | 30.0 |
| VDR | 2 | 9 | VDR --> CTNNB1 --> FOXO1 | activate | Cholecalciferol, Paricalcitol, Doxercalciferol | 29.2 |
| NR3C1 | 2 | 16 | NR3C1 --> CEBPB --> FOXO1 | activate | Hydrocortisone Acetate, Prednisolone Acetate, Mometasone Furoate | 25.8 |
| ESR1 | 2 | 37 | ESR1 --> PGR --> FOXO1 | activate | Estradiol Acetate, Diethylstilbestrol Diphosphate, Toremifene Citrate | 23.1 |
| NR3C2 | 3 | 1 | NR3C2 --> NR3C1 --> CEBPB --> FOXO1 | activate | Desoxycorticosterone Acetate, Desoxycorticosterone Pivalate, Fludrocortisone Acetate | 18.7 |
| MPL | 3 | 2 | MPL --> JAK2 —\| ESR1 —\| FOXO1 | activate | Avatrombopag, Avatrombopag Maleate, Lusutrombopag | 18.3 |
| CASR | 3 | 2 | CASR --> FYN —\| ITCH —\| FOXO1 | activate | Cinacalcet Hydrochloride | 18.3 |
| GHR | 3 | 4 | GHR --> JAK2 —\| ESR1 —\| FOXO1 | activate | — | 17.4 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 1 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 4 downstream branches

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

**NFE2L2** — **no ligand pocket** · **a curated route exists**

| Drug | Mechanism | Binds | Raises | Status |
|---|---|---|---|---|
| Dimethyl fumarate (2013) | degradation blocker | KEAP1 | abundance | approved |
| Omaveloxolone (2023) | degradation blocker | KEAP1 | abundance | approved |
| Monomethyl fumarate (2020) | degradation blocker | KEAP1 | abundance | approved |

- Degradation blockers raise protein LEVELS; conformational stabilisers raise activity with levels unchanged. Pick your assay accordingly — a western blot cannot see the second kind.

#### No stabilisation route for 7 factor(s)

CEBPA, CEBPD, ELK1, FOXO1, GLIS3, NFYB, TBX2

These have no ligand pocket, no direct drug, and **nothing curated against them or their degradation machinery**. For these factors stabilisation is a research programme, not a route you can take: making it real needs the E3 ligase that degrades the factor *and* a ligand for that ligase, and neither is published here.

The mechanism class is real — it has shipped approved drugs, listed once below — but those act on **other proteins entirely**, and none of them is evidence that any factor above is tractable.

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Cell-type reach

Three things must overlap: where the drug goes, where the gene is expressed, and where the transcription factor is. Any empty intersection kills the route regardless of motif quality or potency.

| TF | Co-expressed with gene | Shared tissues |
|---|---|---|
| NFYB | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| CEBPD | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| CEBPA | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| TBX2 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| ELK1 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| GLIS3 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| FOXO1 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| NFE2L2 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| TBX3 | yes (15) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |
| STAT1 | yes (16) | breast, cervix, uterine, endometrium, esophagus, kidney, lun |

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

*No factor tested (NFYB, CEBPD, CEBPA, TBX2, ELK1, GLIS3) has trans-eQTL support for this gene.*

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

*Accounting for all 1 candidates: 0 already linked, 1 novel and non-dominated, 0 novel but beaten on every axis by another candidate, 0 with no advantage. All of them are in results.json under `improvements.all_candidates`.*

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

#### Not previously linked (no incumbent to beat)

| Candidate | Route | TF | Strength | Specificity | Efficiency | Wins on | Trade-offs |
|---|---|---|---|---|---|---|---|
| Omaveloxolone | via TF | NFE2L2 | 0.512 | 0.42 | 0.725 | stronger activation, superior efficiency | none |

*No incumbent exists, so these cleared no bar — they are simply the candidates that scored above 0.5 on at least one axis. Hypotheses to test.*

#### Non-dominated set (nothing beats these on every axis)

- **Omaveloxolone** (NEW link) — no incumbent for this gene; nothing to beat

### Strength of evidence for each factor route

Two links between a drug and your gene come from evidence and differ from row to row: whether the drug engages its target at a dose that can be reached, and whether the factor is likely to act on this promoter. Each is a 0 to 1 heuristic and not a calibrated probability. They are shown as they are, and they are not multiplied into a chance that the gene moves.

| Drug | TF | Drug engagement | Promoter activation | Weaker link | Knockout evidence |
|---|---|---|---|---|---|
| Omaveloxolone ★ | NFE2L2 | 0.775 | 0.557 | promoter activation | gene rises when the factor is lost, 1 of 27 |

★ = named in the comparison against the standard of care above, pinned here so the two sections cannot disagree.

**Knockout evidence** says whether your gene is among the genes that changed when the factor was knocked out, knocked down or overexpressed in published experiments. A response shows the factor can move the gene somewhere, not that it will in your cells or that a drug acting on the factor will. "Not among top changes" is weak, because each list holds only the few hundred genes that changed most.

Where the knockout data and the annotated role disagree, one of them is wrong for this gene, and it matters whether the factor is to be activated or blocked:

- NFE2L2: annotated as an activator, yet the gene rises when it is lost in the knockout data

Order candidates by tier and score, as the tables above do. These columns say where the evidence behind a route is thin.

## Virtual screening

A screen needs a model, and a model needs data. This one trains on measured binding at the receptor, so it runs only where ChEMBL holds at least 20 ligands at pChEMBL 6 or better, and only when the trained model separates known actives from decoys at CV AUC 0.7 or above.

**Screened:** NFE2L2.

**This is why the same receptor families appear here whatever gene you run.** The filter selects for how much medicinal chemistry a receptor has attracted, not for anything about your promoter. Nuclear receptors, and the retinoid and steroid receptors above all, have ligand series running to thousands of compounds; most transcription factors have none at all. So a receptor reaching this section means it is well studied, and the drugs listed under it are that receptor's pharmacology, not a result about your gene.

Read these rows as *candidate ligands for a factor that has a site in your promoter*. Whether moving that factor moves your gene is a separate question, and the activity and evidence sections above speak to it. A long list here is evidence about the receptor's chemistry; it is not evidence that the gene will respond.

### Virtual screen of approved drugs against NFE2L2

Binding model, CV AUC 0.994 — trained on 249 measured actives vs 747 decoys, screened 3311 approved drugs.

*Direction not annotated: only 1 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

**Screen rejected:** no known actives appear among the top hits — the model is not ranking them above the library, so the novel list cannot be trusted No candidates are reported from it.

### Drug–gene links with no recorded direction

133 interaction(s) across 4 factor(s) come from curated sources that record a drug–gene relationship **without saying which way it goes**. They are held back from the recommendations above, because the entire question here is direction and an undirected link cannot answer it.

They are still leads — someone curated a real relationship — but they need a literature check, not a prediction. Many are chemotherapeutics or broad stressors whose association with a transcription factor reflects cell-wide perturbation rather than a specific mechanism.

| TF | Drug | Reported by | DGIdb score |
|---|---|---|---|
| ELK1 | Ganciclovir | NCI | 11.71 |
| TBX2 | Atenolol | ClinPGx | 2.6 |
| FOXO1 | Epirubicin | ClinPGx | 0.52 |
| FOXO1 | Cyclophosphamide Anhydrous | ClinPGx | 0.21 |
| FOXO1 | Fluorouracil | ClinPGx | 0.14 |
| NFE2L2 | Dimethyl Fumarate | DTC; Drug Repurposing Hub | 0.06 |
| NFE2L2 | Atorvastatin Calcium | DTC | 0.03 |
| NFE2L2 | Atovaquone | DTC | 0.03 |
| NFE2L2 | Furazolidone | DTC | 0.03 |
| NFE2L2 | Hydrocortisone Acetate | DTC | 0.03 |
| NFE2L2 | Nandrolone Decanoate | DTC | 0.03 |
| NFE2L2 | Neostigmine Methylsulfate | DTC | 0.03 |
| NFE2L2 | Dicumarol | DTC; NCI | 0.02 |
| NFE2L2 | Adefovir Dipivoxil | DTC | 0.02 |
| NFE2L2 | Deferoxamine Mesylate | DTC | 0.02 |
| NFE2L2 | Dydrogesterone | DTC | 0.02 |
| NFE2L2 | Hydrocortisone Sodium Succinate | DTC | 0.02 |
| NFE2L2 | Methyclothiazide | DTC | 0.02 |
| NFE2L2 | Methylprednisolone Acetate | DTC | 0.02 |
| NFE2L2 | Oxyphenonium Bromide | DTC | 0.02 |
| NFE2L2 | Pioglitazone Hydrochloride | DTC | 0.02 |
| NFE2L2 | Pipobroman | DTC | 0.02 |
| NFE2L2 | Prednisolone Acetate | DTC | 0.02 |
| NFE2L2 | Rosiglitazone Maleate | DTC | 0.02 |
| NFE2L2 | Stanozolol | DTC | 0.02 |
| … 108 more | | | |

*Pass `--include-unknown-direction` to fold these into the main tables anyway — but every such row will read UNKNOWN / unclear, which is what they are.*

### Approved drugs acting on these factors

| TF | TF role | Do what to this TF | Drug | Action | Drug direction | First approval | Source |
|---|---|---|---|---|---|---|---|
| NFE2L2 | activator | increase | Omaveloxolone | ACTIVATOR | up | 2023 | ChEMBL |

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

- Insert: 20:54171754-54172953 (1000 bp upstream / 200 bp downstream of the ENST00000395954 TSS)
- Vector: promoterless luciferase (pGL4.10-style) or GFP for FACS-based screening; luciferase for arrayed compound screens, GFP for pooled CRISPRa
- Controls:
  - Promoterless vector: establishes the floor. Without it you cannot tell activation from vector background.
  - Constitutive promoter (CMV/EF1a): establishes the ceiling and reports transfection or transduction efficiency per well.
  - Scrambled-motif variant of your top predicted element: the single most informative control here. If a drug induces the wild-type reporter but not the mutant, the effect runs through that element. If it induces both, the mechanism is elsewhere and the prediction was wrong.
  - Dual-reporter normalisation (Renilla or similar) to separate real induction from changes in cell number or general translation — cytotoxic compounds otherwise read as repressors.
- An episomal reporter lacks the chromatin context of the endogenous locus. A hit should be confirmed at the native gene by RT-qPCR before it is believed.

### CRISPRa screen

- Approach: dCas9-VP64 with SAM (MS2-p65-HSF1) or dCas9-VPR for stronger activation
- Focused library: 14 TFs prioritised by this run (SOX2, POU5F1, POU2F1, ZNF530, EBF2, EBF4, MTF1, NFIB, ZNF707, ZNF320...)
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
