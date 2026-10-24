# Upregulation strategy for BUB1B — promoter / TF / drug analysis

Objective: **increase** transcription of BUB1B.

> **Only FDA approved drugs are shown.** 16 entries are left out. 11 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 4 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

## What would have to happen

**BUB1B dosage: unknown.**

*Read from gnomAD constraint v4.1; ClinGen dosage curation 07 Oct,2026.*

## Evidence available in this run

**✓ INTERPRETABLE** — Scan is enriched and measured evidence is available (chip). The ranking rests on more than annotation.

| Source | State | What it returned |
|---|---|---|
| chip | ok | local BED: peaks for 534 factor(s) over the window |
| depmap | ok | 17 of 17 factor(s) scored in CRISPRGeneEffect.csv |
| factor_abundance | ok, nothing found | no compound moved a factor's own transcript the needed way in enough contexts to report (8 factor(s) asked, 0 of them recorded regulators the motif scan did not rank) |
| fame | ok | 4 candidate(s) found in the fame table over 11 gene(s); 1 are recorded for most of them |
| knockout | ok | 7 of 17 factor(s) were perturbed in the knockout data; the gene moved for 3 |
| proliferation | ok, nothing found | 5 candidate(s) checked against the mitotic programme; 0 move it along with the gene |
| protein_abundance | ok, nothing found | BUB1B is not in the table, so no compound was measured against its protein. That is an absent measurement, not a negative one |
| regulon | ok | 590 TFs with targets called from ChIP occupancy. A peak near a gene is binding, not regulation |
| utr | ok | 5'UTR 152 bp, 3'UTR 364 bp; 1 feature(s) beat a shuffle of the same sequence |

## A mitotic gene

BUB1B is a mitotic gene. Its transcript varies with the cell-cycle phase, so a compound that changes the phase mix of a culture can move it without acting on it: read any RT-qPCR from synchronised or phase-resolved cells. Compounds suspected to hold cells in mitosis are marked "halts mitosis", and their routes score half, because changing the cell's mitosis moves this gene for reasons unrelated to treating it.

## How much of the problem this run examined

**5 of the 6 steps between a gene and a working protein were searched for this gene; 1 are named without being searched.** The ranking holds 9 transcriptional route(s) against 21 from every other step combined. That ratio is a property of this package, not of the biology: it searches transcription hardest because that is what it was built to do. If the pathogenic variant is in the coding sequence or a splice site, the limiting step is one of the others, and a full table above is not evidence that it was considered.

| Step | This run | Candidates | What is missing |
|---|---|---|---|
| transcription | searched | 9 | — |
| splicing and isoform choice | searched | 20 | no drug was searched against them. Splice modulation has an approved small molecule (risdiplam), but nothing here looks for one against this gene's exons, and no public corpus of… |
| mRNA stability | searched | 1 | the element is annotated and no drug was matched to it. No corpus of transcript half-life under drug is read, and the microRNA that represses it is not identified. |
| translation | searched | 0 | a uORF is a lever and no drug was matched to it. Nothing reads ribosome profiling, so whether protein is low while message is normal - the signature that would point here - still needs the… |
| folding | listed only | — |  |
| degradation | searched | 0 | a drug-to-protein-abundance table is on disk and does not quantify BUB1B at all, so nothing in it bears on this protein. That is an absent measurement, not a negative one. |

What would close each gap:

- **folding** — measured ligands for this protein; the screen refuses below the number needed to train a model, rather than guessing.

*"Listed only" means the mechanism and its drug classes are named in this report and nothing was looked up for this gene. It is a menu, not a search.*

### Measured on the protein

**BUB1B is not in the drug-to-protein table, so no compound was checked against its protein.** That is an absent measurement and not a negative one: the table quantifies a few thousand proteins in one cell line, and this is not one of them.

### Moving the gene, or moving the cell cycle

**No candidate compound moved the mitotic programme with BUB1B unusually often.** Each was checked against 25 other cell-cycle transcripts across its LINCS signatures, and compared with the 1.9% of all signatures that move the programme anyway. None reached the 95th percentile, so none of these rows stands out as a change in how many cells are dividing. That is not a clean bill: the test needs the compound to be in the files at all.

### How much of this order is fame?

**1 of 4 checked candidate(s) are recorded for most of the 11 genes measured here.** This is the control that matters most for reading the order below. Tested across these runs, a baseline that only counts how many other genes record a drug predicts CTD better (AUC 0.93) than every route in this package combined (0.53), and among drugs of equal fame the routes fall to chance. A compound recorded for every gene is telling you about the literature, not about yours. Promiscuous rows lose a little score within their tier and say so; none changes tier.

| Compound | Genes it is recorded for |
|---|---|
| Valproic acid | 8 of 11 (73%) |

## The plate: every compound, in no order

**9 distinct compounds, in no order.** Chemical families are not folded together here and no compound is demoted for being well studied, because a plate measures what a compound does to cells and neither of those is about that. Rows are alphabetical inside each group; the groups say where the evidence places a compound, and for 0 of the 9 it places them nowhere - those rest on the gene having moved, with no mechanism attached. 0 further compound(s) were found and left off on safety grounds, listed below: that judgement is clinical, and some of them are reasonable mechanistic controls on a plate.

### transcription (9)

| Compound | Acts on | Site evidence | Routes | Found by |
|---|---|---|---|---|
| Azacitidine | chromatin (genome-wide) | — | 1 | catalogue |
| Belinostat | chromatin (genome-wide) | — | 1 | catalogue |
| Decitabine | chromatin (genome-wide) | — | 1 | catalogue |
| Givinostat Hydrochloride | chromatin (genome-wide) | — | 1 | catalogue |
| Panobinostat | chromatin (genome-wide) | — | 1 | catalogue |
| Romidepsin | chromatin (genome-wide) | — | 1 | catalogue |
| Tazemetostat | chromatin (genome-wide) | — | 1 | catalogue |
| Valproic acid | chromatin (genome-wide) | — | 1 | catalogue |
| Vorinostat | chromatin (genome-wide) | — | 1 | catalogue |

**No compound is placed at acts on the gene or its product, splicing and isoform choice, mRNA stability, translation, folding, degradation.** Those steps were searched and returned annotations rather than compounds, so they contribute nothing to a plate however the groups above are weighted. Giving them equal weight here would mean weighting empty sets.

`screen_list.tsv` holds the same rows for a pipetting robot or a plate map.

## What the BUB1B transcript itself does to its own output

Two of the six steps between a gene and a working protein happen after the message is made: how readily it is translated, and how long it survives. Both are partly set by the untranslated regions, which are read here against a shuffle of the same sequence, so a feature counts only when there is more of it than that base composition gives by chance.

No upstream ORF in the 152 bp 5'UTR, so translational repression by that route is not available here. The 364 bp 3'UTR is enriched for the AU-rich pentamer: 4 against 1.22 expected (p 0.0196). That marks the transcript for rapid decay, so its half-life is a lever - and again, no compound was matched to it.

**149 named regulator(s) of the BUB1B transcript** from `data\transcript_sites.tsv`. RBP: AARS, AATF, ABCF1, AGGF1, AKAP1, AKAP8L, APOBEC3C, AQR, BCCIP, BCLAF1 (148 of 148 measured rather than predicted) m6A: m6A site (1 of 1 measured rather than predicted)

Read that count the right way round. An eCLIP peak is matched across the whole gene body, introns included, so a well expressed gene matches most of the proteins the atlas assays at all - BUB1B is not unusual in matching this many. Every row is a real measured footprint, and the list ranks nothing: it is the set of proteins seen on this transcript, not the ones that set its level. Narrowing it needs the peak positions against the UTRs, which this table does not carry.

**No compound was searched against any of this.** A uORF, a quadruplex or a decay element is a place an intervention could act, not an intervention; and a motif in the annotated UTR is not proof the element is used in your cells.

## Exons that may brake BUB1B output

Some genes make part of their transcript with an extra exon that puts a stop codon early, so the message is destroyed. An antisense oligonucleotide that makes the cell skip that exon turns those transcripts into working ones and raises the protein, from the working copy only. 12 decay isoform(s) of BUB1B in Ensembl, 11 of them with an exon the coding isoforms do not use: 2 novel exon(s) and 18 extended or shifted coding exon(s).

No decay isoform of this gene carries measurable message in Cells_Cultured_fibroblasts. The candidate brake exons below sit on transcripts that are not detectably made in this tissue, so there is nothing for an oligonucleotide to redirect. The non-productive message that does exist is mostly retained intron (13.2%), which a brake-exon oligonucleotide does not address - a different intervention would be needed to recover it. 10 decay isoform(s) Ensembl lists now are newer than the GENCODE build GTEx was quantified on, so they have no measurement at all - which is not a measurement of zero. Measured on 11 of the 31 isoforms Ensembl lists, 8.6 TPM in total.

| Kind of transcript | Share of this gene's message |
|---|---|
| productive | 86.8% |
| retained intron | 13.2% |
| no defined coding sequence | 0.0% |
| decay | 0.0% |

Measured in Cells_Cultured_fibroblasts (gtex_v8). An isoform Ensembl lists that this dataset never quantified is counted as unmeasured and left out of the shares above, rather than treated as zero.

| Exon (genomic) | Length | Kind | Decay isoforms using it | Message they carry |
|---|---|---|---|---|
| 40166359–40166482 | 124 bp | novel exon | 1 | 0.0% |
| 40176383–40176436 | 54 bp | novel exon | 1 | 0.0% (1 unmeasured) |
| 40161026–40161255 | 230 bp | extended or shifted coding exon | 2 | 0.0% (2 unmeasured) |
| 40200933–40200980 | 48 bp | extended or shifted coding exon | 2 | 0.0% (2 unmeasured) |
| 40220564–40221134 | 571 bp | extended or shifted coding exon | 2 | 0.0% (2 unmeasured) |
| 40160984–40161255 | 272 bp | extended or shifted coding exon | 1 | 0.0% (1 unmeasured) |
| 40161006–40161255 | 250 bp | extended or shifted coding exon | 1 | 0.0% (1 unmeasured) |

13 more extended or shifted coding exons are not shown; most are first-exon start variants, which are not brakes.

**This is a list of places to look, not a design.** An annotated exon is one Ensembl has seen in some tissue, and the share beside it is a median across GTEx donors in one tissue rather than a measurement in your cells. Whether the exon is included in them, how often, and whether skipping it raises protein is a question for RNA sequencing and a splicing assay. Nothing here picks an oligonucleotide sequence or checks it against the rest of the transcriptome.

## Model-predicted medications for BUB1B

Compounds this analysis predicts would **raise** BUB1B, with the reasoning for each. Everything here is a **hypothesis generated from annotation**, not a recommendation — the tier on each row tells you how far the claim is reaching.

| # | Medication | What it is | Approved to treat | Acts on | BBB? | Tier | Grade | Score |
|---|---|---|---|---|---|---|---|---|
| 1 | **Vorinostat** | histone deacetylase inhibitor | cutaneous T-cell lymphoma | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 2 | **Romidepsin** | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 3 | **Belinostat** | HDAC inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 4 | **Panobinostat** | HDAC inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 5 | **Azacitidine** | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 6 | **Decitabine** | DNA methyltransferase inhibitor | myelodysplastic syndrome | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 7 | **Tazemetostat** | EZH2 inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 8 | **Givinostat Hydrochloride** | HDAC inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 23.4 |
| 9 | **Valproic acid** | weak HDAC inhibitor | not established here | chromatin (genome-wide) | — | 3 | — | 19.9 |

**Every row is off-label.** The drug was approved for what the fourth column says, and is listed here for something else. That column is the baseline this proposal departs from, and the distance between it and your gene is part of what you are judging - a topical acne gel and an infused cytotoxic can sit on adjacent rows with similar scores.

*Ordered by evidence tier first, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. **BBB?** is not part of the ordering — it is an applicability filter you apply yourself, since whether brain exposure is wanted or unwanted depends on where your gene needs to move. See the note under the graded table for what a yes does and does not establish.*

### Why each one

Every route below is **evidence tier 3** — a human curator annotated this mechanism.

#### Through **chromatin (genome-wide)** — 9 candidates

- The factor has no activator/repressor annotation, so the direction of its effect on your gene is assumed, not known.

- ⚠ acts genome-wide: will move thousands of promoters including yours. Useful as a positive control that your reporter CAN be induced; almost never an acceptable therapeutic route on specificity grounds

1. **Vorinostat** — HDAC inhibitor
2. **Romidepsin** — HDAC inhibitor
3. **Belinostat** — HDAC inhibitor
4. **Panobinostat** — HDAC inhibitor
5. **Azacitidine** — DNMT inhibitor
6. **Decitabine** — DNMT inhibitor
7. **Tazemetostat** — EZH2 inhibitor
8. **Givinostat Hydrochloride** — HDAC inhibitor
9. **Valproic acid** — weak HDAC inhibitor

*Not medical advice. "Approved" means approved for its own indication — using one of these to move an unrelated gene is off-label repurposing with the full systemic pharmacology attached. The experimental sections at the end of this report are how you would actually test any of it.*

## Promoter of BUB1B (ENST00000287598)

- Transcript chosen: **MANE Select**
- Assembly: from GTF  |  Strand: +
- TSS: chr15:40161069
- Window: chr15:40160069-40161268 (-1000 / +200, 1200 bp)
- Motifs scanned: 879 at p ≤ 1e-05 per position
- 254 of those motifs were not scanned: in this window their best possible match is more common than p = 1e-05, so any hit would be a chance match at a higher rate than stated. A looser `--pvalue` scans them, with more chance hits.
- Non-redundant sites: 17 across 17 TFs
- **Measured occupancy splits the scan.** 10 of 17 factors have a measured peak over this window and 7 do not. A shuffled copy of the window cannot produce a peak, so the share of this scan the null expects by chance - about 8 of its 17 factors - is attributed to those 7: about 100% of them are consistent with chance, against roughly 7% of the measured ones. This matters for the ranking because NR1H4 has a ligand and no measured occupancy: every drug reaching the gene through it rests on a predicted site. RXRA has both a ligand and a measured peak, which is the stronger footing. A peak says something bound here in the cell type the atlas used. It does not say this factor is what bound, and it is not a measurement in your cells.

*Scan quality: 17 sites against 8 expected by chance over 739 scannable bases (2.2x).*

*about 45% of these 17 sites (range 24% to 65%) would be expected on shuffled sequence of this window, so roughly 9 are above chance; among the best-scoring sites the share is lower: 31% of the top 5, 39% of the top 10. Estimated from 30 shuffled copies that keep the repeat mask and every dinucleotide count. The formula that assumes each motif matches at exactly its stated p-value would expect 9 and is not used. Nothing here says which sites are real; only ChIP evidence separates them.*

### Ligand-activated (nuclear) receptors with sites here

These are the directly druggable ones.

| TF | Role | Sites | Best rel. score | Best site (vs TSS) | ChIP peaks | Score |
|---|---|---|---|---|---|---|
| RXRA | unknown | 2 | 0.907 | AGGCCACTGACCT (-413) | 1 | 69.2 |
| NR1H4 | unknown | 1 | 0.882 | AGGTCAGTGGCCT (-413) | 0 | 51.3 |
| NR4A2 | unknown | 1 | 0.907 | AGGCCACTGACCT (-413) | 0 | 42.2 |
| NR2F6 | unknown | 1 | 0.876 | TAGGTCAGTGGCCTA (-414) | 0 | 41.1 |
| NR2F1 | unknown | 1 | 0.872 | TAGGCCACTGACCTA (-414) | 0 | 40.9 |

#### Why the drug tables below look the way they do

**17 factor(s) had predicted sites here. 2 of them have an approved drug, and all 2 are nuclear receptors.**

Roughly 48 of the ~1,600 human transcription factors have a ligand-binding pocket, and the steroid, retinoid, thyroid and vitamin-D families are most of them. So the drug tables are not a ranking that hormone drugs won — they are what is left after every factor without a ligand has been dropped.

**15 of the top 25 factors by score have no drug and appear in no drug table:** OSR2, THAP11, KLF3, SP3, KLF16, ZNF76, ZNF143, FLI1, TFAP4, FOXD2, NR4A2, NR2F6 ….

That is where the gene-specific signal is. A half-site like AGGTCA occurs every few kilobases by chance, so nuclear receptors surface for almost any promoter; these factors surfaced for *this* one. The routes that can reach them are the upstream-network, stabilisation and protein-level sections, not the ligand tables.

### Top-ranked transcription factors

| Rank | TF | Role | NR | Other genes it controls | Cells need it? | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | RXRA | unknown | yes | ~1200 | a few lines | 2 | 0.907 | -413 | 1 | 4 | 69.2 |
| 2 | OSR2 | unknown |  | not counted | a few lines | 1 | 1.0 | +51 | 1 | 0 | 66.5 |
| 3 | THAP11 | unknown |  | 1588 | some lines | 2 | 0.98 | -160 | 1 | 0 | 66.1 |
| 4 | KLF3 | unknown |  | not counted | a few lines | 1 | 0.989 | -58 | 1 | 0 | 66.0 |
| 5 | SP3 | unknown |  | ~5000 | a few lines | 1 | 0.979 | -60 | 1 | 0 | 65.6 |
| 6 | KLF16 | unknown |  | not counted | some lines | 1 | 0.973 | -60 | 1 | 0 | 65.4 |
| 7 | ZNF76 | unknown |  | not counted | a few lines | 2 | 0.927 | -152 | 1 | 0 | 64.7 |
| 8 | ZNF143 | unknown |  | not counted | yes, most | 1 | 0.923 | -153 | 1 | 0 | 61.2 |
| 9 | FLI1 | unknown |  | not counted | a few lines | 1 | 0.866 | -426 | 1 | 0 | 54.5 |
| 10 | TFAP4 | unknown |  | not counted | some lines | 1 | 0.866 | -426 | 1 | 0 | 54.5 |
| 11 | NR1H4 | unknown | yes | ~300 | a few lines | 1 | 0.882 | -413 | 0 | 3 | 51.3 |
| 12 | FOXD2 | unknown |  | not counted | no | 1 | 0.931 | -211 | 0 | 0 | 46.6 |
| 13 | NR4A2 | unknown | yes | not counted | a few lines | 1 | 0.907 | -413 | 0 | 0 | 42.2 |
| 14 | NR2F6 | unknown | yes | not counted | a few lines | 1 | 0.876 | -414 | 0 | 0 | 41.1 |
| 15 | NR2F1 | unknown | yes | not counted | a few lines | 1 | 0.872 | -414 | 0 | 0 | 40.9 |
| 16 | PRDM9 | unknown |  | not counted | a few lines | 1 | 0.889 | -997 | 0 | 0 | 37.5 |
| 17 | ZNF768 | unknown |  | 826 | a few lines | 1 | 0.979 | -992 | 0 | 0 | 34.0 |

**Other genes it controls** is the collateral question. A drug acting through a factor moves everything that factor regulates, so this column says how wide the effect will be. Read it against the score: a strong route through a factor that controls two thousand genes is a worse starting point than a weaker route through one that controls thirty.

A count written `~N` is a bundled estimate for a factor the regulon file does not list (3 here), not a count of curated edges. Treat it as an order of magnitude.

`not counted` means the regulon table has no entry and there is no bundled estimate either (12 here). It is not a small number. Most of the roughly 1,600 human transcription factors have no curated regulon at all.

Counts come from literature curation, so they measure what has been published as much as what the factor does. SP1 has thousands of recorded targets partly because it regulates many genes and partly because thousands of people assayed it.

### What a cell loses without each factor (DepMap)

The regulon column counts how many genes a factor controls, from literature curation, which measures how much anyone has published as much as what the factor does. This is the same collateral question answered by measurement: CRISPR knockout of each factor across roughly 1,100 cancer cell lines.

It also reaches factors the curated tables cannot. A factor nobody has written a regulon for still has a knockout phenotype here.

17 of 17 factor(s) were scored in DepMap. 13 (FLI1, KLF3, NR1H4, NR2F1, NR2F6) show a dependency confined to a minority of lines, which is the informative case: the requirement has a context, and whether your tissue shares it is the question to ask next. 1 (FOXD2) was dispensable throughout. That bounds the proliferative cost of removing them and says nothing about a differentiated tissue.

| Factor | Cells need it? | Mean effect | Lines dependent | Most dependent lineage |
|---|---|---|---|---|
| ZNF143 | essential in most lines | -0.52 | 632/1150 (55%) | Lymphoid (-0.64) |
| THAP11 | essential in a substantial minority | -0.43 | 428/1150 (37%) | Lymphoid (-0.69) |
| KLF16 | essential in a substantial minority | -0.36 | 292/1150 (25%) | CNS/Brain (-0.49) |
| TFAP4 | essential in a minority of lines | -0.35 | 263/1150 (23%) | Myeloid (-0.54) |
| SP3 | essential in a handful of lines | -0.13 | 17/1150 (2%) | Myeloid (-0.20) |
| NR1H4 | essential in a handful of lines | -0.11 | 7/1150 (1%) | Bladder/Urinary Tract (-0.16) |
| NR2F6 | essential in a handful of lines | -0.07 | 2/1150 (0%) | Prostate (-0.15) |
| NR4A2 | essential in a handful of lines | -0.07 | 4/1150 (0%) | Ovary/Fallopian Tube (-0.10) |
| RXRA | essential in a handful of lines | -0.07 | 31/1150 (3%) | Bladder/Urinary Tract (-0.27) |
| ZNF76 | essential in a handful of lines | -0.07 | 1/1150 (0%) | Head and Neck (-0.12) |
| ZNF768 | essential in a handful of lines | -0.06 | 1/1150 (0%) | Myeloid (-0.11) |
| KLF3 | essential in a handful of lines | -0.04 | 12/1150 (1%) | Pancreas (-0.24) |
| NR2F1 | essential in a handful of lines | -0.04 | 4/1150 (0%) | CNS/Brain (-0.09) |
| FOXD2 | dispensable across the panel | -0.02 | 0/1150 (0%) | — |
| OSR2 | essential in a handful of lines | -0.02 | 1/1150 (0%) | Liver (-0.06) |
| FLI1 | essential in a handful of lines | +0.00 | 39/1150 (3%) | Bone (-0.58) |
| PRDM9 | essential in a handful of lines | +0.16 | 1/1150 (0%) | Soft Tissue (+0.12) |

**Mean effect** is on DepMap's own scale: 0 means knockout changed nothing, and -1 is the median of known common-essential genes. A line counts as dependent below -0.5.

#### What this does not mean

- **These are cancer cell lines dividing in a dish.** Essentiality here is a statement about proliferating transformed cells. A factor whose job is in a post-mitotic or differentiated cell looks dispensable for reasons that have nothing to do with its importance in your tissue.
- **Knockout is not inhibition.** The protein is removed entirely, including scaffolding roles no ligand would touch. Essential genes are drugged successfully at partial occupancy all the time, so a dependency is a reason for care rather than a veto.
- **Nothing here is about your gene.** This bounds the cost of acting on a factor. Whether acting on it moves your gene is the rest of the report.
- **A dispensable result is not permission.** It bounds the proliferative cost and nothing else.

**The selective ones are the interesting ones.** FLI1, KLF3, NR1H4, NR2F1, NR2F6, NR4A2 matter in a minority of lines rather than all of them, so the requirement has a context. The lineage column is where to look for what that context is, and whether your tissue shares it is the question this panel cannot answer.

### Every route that acts on this gene

9 routes across 1 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

**Why established binders sit at tier 2, just under the direct agents.** A drug already measured to bind one of your factors has half its chain settled: the drug-to-receptor step has assay data behind it, so the pharmacology is de-risked and only the receptor-to-gene step is a hypothesis. A curated tier-3 ligand is annotation on *both* halves. But they stay below tier 1, because a direct agent acts on your gene itself and skips the receptor-to-gene inference entirely — a stronger claim however well-characterised the binder is. The cost of the indirect route is in the row: everything else that factor regulates moves with your gene.

- 9 × curated — mechanism annotated by a human curator

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 3 | chromatin-level | Vorinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Romidepsin | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Belinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Panobinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Azacitidine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Decitabine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Tazemetostat | chromatin (genome-wide) | 3 | EZH2 inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Givinostat Hydrochloride | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Valproic acid | chromatin (genome-wide) | 3 | weak HDAC inhibitor | 0.7 | 0.05 | 0.45 | 19.9 |

Best in each category:

- **chromatin-level** — Vorinostat (tier 3, score 23.4)

*No route found in: direct: gene or transcript, direct: gene product, direct: predicted product binder, promoter DNA binder, TF: orthosteric ligand, TF: inferred activity change, TF: abundance change, protein abundance change, TF: predicted pocket binder, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, TF: degradation blocker, TF: known regulator, upstream signalling, literature: reported effect, off-pathway pocket binder, measured phenocopy, empirical signature, virtual screen hypothesis. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

| Category | What it changes | Your drugs in it | Assay that sees it |
|---|---|---|---|
| chromatin-level | **abundance, genome-wide** | Vorinostat, Romidepsin, Belinostat, Panobinostat … +5 | RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened |

*Every category present here works through transcription, so all of them move your gene's product **abundance** and an RT-qPCR is the right first readout. The one exception in the taxonomy — a drug on the gene product itself, which changes activity only — did not appear in this run.*

**chromatin-level** (9 found)

- *Yours:* Vorinostat, Romidepsin, Belinostat, Panobinostat, Azacitidine, Decitabine, Tazemetostat, Givinostat Hydrochloride … +1 more
- HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.
- *Why:* Reliable positive control that your reporter *can* be induced.
- *Cost:* Specificity ~0.05. These move thousands of promoters. Almost never an acceptable therapeutic route on that basis alone.
- *Differs:* Not targeted at your gene in any sense; it is in the table so you can see what a non-specific option scores.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **direct: gene product** — Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- **direct: predicted product binder** — A compound predicted to bind this gene's own protein and not recorded as acting on it - either because its chemistry resembles the protein's known ligands, or because it docks into the protein's pocket.
- **promoter DNA binder** — A small molecule with a documented preference for a sequence feature this promoter carries - a quadruplex, a GC-box, an AT tract - proposed because the feature is enriched against shuffles of this same sequence.
- **TF: orthosteric ligand** — Binds the transcription factor's natural ligand pocket.
- **TF: inferred activity change** — A compound whose own LINCS signature moves this factor's known targets, found without any record that it binds the factor.
- **TF: abundance change** — A compound measured to change how much of the factor's own transcript exists, rather than what the factor does.
- **protein abundance change** — A compound measured to change how much of your gene's own PROTEIN exists, in a proteomics experiment under drug
- **TF: predicted pocket binder** — A compound predicted by docking to fit the factor's ligand pocket, from a structure rather than from any measured activity.
- **TF: allosteric modulator** — Binds a site away from the ligand pocket, amplifying or damping the endogenous ligand's effect (PAM/NAM).
- **TF: stabiliser or glue** — Locks the factor in an active conformation, or glues it to a partner protein.
- **TF: endogenous ligand** — Changes how much natural ligand exists, by drugging the enzyme that makes or destroys it.
- **TF: degradation blocker** — Stops the factor being destroyed, so more of it accumulates.
- **TF: known regulator** — Drugs that act on a factor CollecTRI or perturbation data records as raising or lowering your gene, whether or not the motif scan ranked that factor.
- **upstream signalling** — Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- **literature: reported effect** — A chemical CTD records as raising or lowering this gene's mRNA, with no mechanism stated.
- **off-pathway pocket binder** — A docking run against a protein that is not one of the factors found at this promoter. Usually an upstream node someone chose to dock by hand.
- **measured phenocopy** — A compound whose effect on cell viability across the DepMap panel tracks what happens when this factor is knocked out, so it behaves as though it removes the factor.
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **virtual screen hypothesis** — An approved drug predicted to bind one of your TFs, never tested against it.

### What these drugs are

Every route in this report is off-label by construction: the drug was approved for one thing and is being considered for another. The approved use is the baseline that departure is measured from, so it is listed here once per drug rather than repeated in each table above.

| Drug | What it is | Approved to treat |
|---|---|---|
| Azacitidine | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia |
| Decitabine | DNA methyltransferase inhibitor | myelodysplastic syndrome |
| Romidepsin | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma |
| Vorinostat | histone deacetylase inhibitor | cutaneous T-cell lymphoma |

**No description established for 5**: Belinostat, Givinostat Hydrochloride, Panobinostat, Tazemetostat, Valproic acid. That means this package could not find one offline, not that the compound has no approved use. Several are environmental agents and research chemicals from the CTD tables, which have no indication because they are not drugs.

Descriptions come from: 4 curated, 5 unknown. None of these is a label. The curated text is a summary written for orientation, ChEMBL's indication list includes uses studied but never approved, and a class line gives the kind of drug rather than its indication. Check the prescribing information before acting on any row.

### Routes a promoter analysis cannot see

Everything above acts on transcription. These act after it, around it, or on a different gene entirely. For some genes that is where the whole problem is: if the transcript is destroyed, the allele is silenced, or the protein is degraded too fast, more transcription does not help and may not be possible.

#### Could a relative do the job instead?

Whether a paralog could substitute for this gene is unanswered. It is worth asking explicitly, because the question 'restore this function' has more answers than 'raise this transcript'.

- BUB1B is not in the curated redundancy table, which is short and certainly incomplete. Absence here is not evidence that no paralog substitutes — it means nobody has written one down in this file.

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

**BUB1B is destroyed by APC/C with Cdc20**, through a KEN box and a D-box, at mitotic exit, once the checkpoint is satisfied. Blocking that specific step raises the protein without touching transcription, and without the whole-proteome cost of a proteasome inhibitor. No approved drug in this package's tables acts on CDC20 or ANAPC2 or ANAPC11 or FZR1, so this route needs chemistry that does not yet exist here. The index is built from annotated mechanisms, so a covalent agent could be missing from it. It is not only the ligase: acetylation of the protein is reported to block its recognition by Cdc20, so the balance - not the ligase alone - sets the level. Verify in the primary literature before building on it; this is a pointer, not a result of this run. BUBR1 is a checkpoint protein that the checkpoint's own target complex destroys, so its level is tied to the cycle it controls. Blocking that destruction also blocks mitotic exit, which is why this route is harder than it looks for a gene of this kind.

**No deubiquitinase is recorded here as rescuing BUB1B.** Degradation is a balance, and only one side of it is known for this protein: the ligase that destroys it is named above, and the enzyme that would reverse that is not. That is a gap in the biology rather than in this package, and it is worth knowing, because the counter-regulation is where a selective lever would sit if one existed. No approved drug activates a deubiquitinase. Essentially all DUB chemistry is inhibitors, and inhibiting the DUB that rescues a protein LOWERS that protein - the opposite of what raising it needs. Treat any DUB compound found below as the wrong direction unless its label says otherwise.

### Drug–gene links with no recorded direction

26 interaction(s) across 3 factor(s) come from curated sources that record a drug–gene relationship **without saying which way it goes**. They are held back from the recommendations above, because the entire question here is direction and an undirected link cannot answer it.

They are still leads — someone curated a real relationship — but they need a literature check, not a prediction. Many are chemotherapeutics or broad stressors whose association with a transcription factor reflects cell-wide perturbation rather than a specific mechanism.

| TF | Drug | Reported by | DGIdb score |
|---|---|---|---|
| NR4A2 | Gemcitabine | CIViC | 0.22 |
| NR1H4 | Imipramine Hydrochloride | DTC | 0.09 |
| RXRA | Etodolac | TdgClinicalTrial | 0.09 |
| RXRA | Imipramine Hydrochloride | DTC | 0.09 |
| NR1H4 | Saquinavir Mesylate | DTC | 0.07 |
| RXRA | Saquinavir Mesylate | DTC | 0.07 |
| NR1H4 | Methylene Blue | DTC | 0.06 |
| NR1H4 | Nimodipine | DTC | 0.06 |
| RXRA | Rosuvastatin | ClinPGx | 0.04 |
| NR1H4 | Phenazopyridine Hydrochloride | DTC | 0.03 |
| NR1H4 | Levoketoconazole | DTC | 0.03 |
| NR1H4 | Raloxifene Hydrochloride | DTC | 0.03 |
| RXRA | Raloxifene Hydrochloride | DTC | 0.02 |
| NR1H4 | Melphalan | DTC | 0.02 |
| RXRA | Vitamin A | Drug Repurposing Hub | 0.02 |
| NR1H4 | Diethylstilbestrol | DTC | 0.02 |
| RXRA | Niclosamide | DTC | 0.02 |
| RXRA | Pimozide | DTC | 0.02 |
| RXRA | Daunorubicin Hydrochloride | DTC | 0.02 |
| RXRA | Dactinomycin | DTC | 0.01 |
| RXRA | Testosterone | DTC | 0.01 |
| RXRA | Docetaxel Anhydrous | ClinPGx | 0.01 |
| RXRA | Progesterone | DTC | 0.01 |
| RXRA | Tretinoin | DTC | 0.01 |
| NR1H4 | Doxorubicin Hydrochloride | DTC | 0.0 |
| … 1 more | | | |

*Pass `--include-unknown-direction` to fold these into the main tables anyway — but every such row will read UNKNOWN / unclear, which is what they are.*

### Approved drugs acting on these factors

| TF | TF role | Do what to this TF | Drug | Action | Drug direction | First approval | Source |
|---|---|---|---|---|---|---|---|
| RXRA | unknown | either | Acitretin | AGONIST | up | 1996 | ChEMBL |
| RXRA | unknown | either | Bexarotene | AGONIST | up | 1999 | ChEMBL |
| RXRA | unknown | either | Sulindac | INHIBITOR | down | — | DGIdb |
| RXRA | unknown | either | Etretinate ⚠withdrawn | AGONIST | up | 1986 | ChEMBL |
| NR1H4 | unknown | either | Chenodiol | AGONIST | up | 1983 | ChEMBL |
| NR1H4 | unknown | either | Ursodiol | AGONIST | up | 1987 | ChEMBL |
| NR1H4 | unknown | either | Obeticholic Acid ⚠withdrawn | AGONIST | up | 2016 | ChEMBL |

---

### How to read this

A PWM hit means the sequence *resembles* a binding site. It is a hypothesis generator, not evidence of binding, and not evidence that binding would increase transcription of this gene in your cell type. Rank order here is a heuristic (weights printed in the JSON output), not a validated predictor. Every row needs wet-lab confirmation: promoter-reporter assay for the element, agonist dose-response with RT-qPCR/protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it.
