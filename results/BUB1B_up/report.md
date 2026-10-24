# Upregulation strategy for BUB1B — promoter / TF / drug analysis

Objective: **increase** transcription of BUB1B.

> **Only FDA approved drugs are shown.** 490 entries are left out. 223 compounds are not in the FDA Orange Book: research compounds, drugs approved only outside the United States, and biologics, which the Orange Book does not list. A further 14 are FDA approved only for local use (topical, eye, inhaled and the like) and are hidden unless the target tissue is skin; `--tissue skin` or `--include-local-only` shows them. `results.json` keeps all of them, and `--include-non-fda` shows everything.

## What would have to happen

**What BUB1B has to do: about 2x from the weak allele.** Only one copy makes anything, and it is impaired. Raising transcription raises the amount of that impaired product; it does not improve it. Whether more of it is enough depends on how much function the allele keeps, which no database here holds - so the fold above is what transcription could deliver, not a dose that is known to be sufficient. What the dead allele produces is not declared (--broken-allele-makes). If it makes a truncated protein rather than nothing, raising transcription amplifies it too. If the problem is the quality of the protein rather than its amount, the routes that address it are stabilisation, readthrough and the protein-level sections, not the promoter. **No route carries a measured effect on this gene.** All 253 rest on annotation, so nothing here has a size to compare against the requirement.

**BUB1B dosage: unknown.**

*Read from gnomAD constraint v4.1; ClinGen dosage curation 07 Oct,2026.*

## Evidence available in this run

**⚠ DEGRADED** — Read the ranking as a reading order, not as evidence. At least one pillar of the analysis is missing. Measured evidence that DID arrive: chip, lincs, ctd, genetics.

- these tests could not be completed: product_screen — which is not the same as their having returned a negative result

| Source | State | What it returned |
|---|---|---|
| chip | ok | local BED: peaks for 534 factor(s) over the window |
| ctd | ok | 184 chemical(s) with curated expression effects on this gene; 71 move it the way you want |
| depmap | ok | 17 of 17 factor(s) scored in CRISPRGeneEffect.csv |
| direction_model | ok, nothing found | 0 of 12 grade(s) got a model call; too few labelled agonists and antagonists to train |
| docking | ok, nothing found | no pocket validated, so no docking hit is reported |
| expression | ok | 28 tissue(s) for BUB1B |
| factor_abundance | ok | 46 compound(s) move a factor's own transcript the needed way; 19 factor(s) asked, 11 of them recorded regulators the motif scan did not rank; 8 factor(s) skipped for an unrecorded role |
| fame | ok | 84 candidate(s) found in the fame table over 11 gene(s); 22 are recorded for most of them |
| genetics | ok | 6 associated phenotype(s) (Open Targets); 12 genome-wide significant GWAS Catalog association(s), mapped or nearest gene, which is weaker than fine-mapped support |
| knockout | ok | 7 of 17 factor(s) were perturbed in the knockout data; the gene moved for 3 |
| lincs | ok | 20 signature hit(s) for this gene |
| product_screen | **FAILED** | BUB1B's own product has 1 measured ligand(s) available, below the 20 needed to train a binding model. A screen without them would rank the library on nothing. The source that answered had only these;  |
| proliferation | ok | 92 candidate(s) checked against the mitotic programme; 2 move it along with the gene |
| promoter_dna | ok, nothing found | no sequence feature beat the shuffled null, so no DNA-directed agent is proposed |
| protein_abundance | ok, nothing found | BUB1B is not in the table, so no compound was measured against its protein. That is an absent measurement, not a negative one |
| regulators | ok | 227 recorded regulator(s) of the gene had no drug route yet; 8 have an approved drug that moves them the needed way (38 route(s)); 219 have none, and are listed rather than dropped |
| regulon | ok | 590 TFs with targets called from ChIP occupancy. A peak near a gene is binding, not regulation |
| screen | ok | 2 receptor(s) screened |
| sensitivity | ok | 13 correlated compound(s) across 6 of 6 factor(s) scanned |
| trans_eqtl | ok, nothing found | 0 of 6 factor(s) have trans-eQTL support |
| upstream | ok | druggable routes for 3 factor(s); 1468 of 1468 reachable upstream node(s) checked across 8 factor(s) |
| utr | ok | 5'UTR 152 bp, 3'UTR 364 bp; 1 feature(s) beat a shuffle of the same sequence |

**These tests could not be completed.** Asking and getting no answer is not the same as asking and being told no, so nothing below rests on them. What each one would have told you:

- **product_screen** — The gene's own protein was not screened for binders. A compound that acts on the protein rather than on its transcription is absent rather than ruled out - and that is often the shortest route to changing what the gene DOES.

## Drugs targeting BUB1B directly

None found in the curated table or in ChEMBL/DGIdb for this direction. The curated list is hand-maintained and goes out of date; check FDA/EMA and clinicaltrials.gov before concluding nothing exists.

## A mitotic gene

BUB1B is a mitotic gene. Its transcript varies with the cell-cycle phase, so a compound that changes the phase mix of a culture can move it without acting on it: read any RT-qPCR from synchronised or phase-resolved cells. Compounds suspected to hold cells in mitosis are marked "halts mitosis", and their routes score half, because changing the cell's mitosis moves this gene for reasons unrelated to treating it.

## How much of the problem this run examined

**6 of the 6 steps between a gene and a working protein were searched for this gene; 0 are named without being searched.** The ranking holds 105 transcriptional route(s) against 21 from every other step combined. That ratio is a property of this package, not of the biology: it searches transcription hardest because that is what it was built to do. If the pathogenic variant is in the coding sequence or a splice site, the limiting step is one of the others, and a full table above is not evidence that it was considered.

| Step | This run | Candidates | What is missing |
|---|---|---|---|
| transcription | searched | 105 | — |
| splicing and isoform choice | searched | 20 | no drug was searched against them. Splice modulation has an approved small molecule (risdiplam), but nothing here looks for one against this gene's exons, and no public corpus of… |
| mRNA stability | searched | 1 | the element is annotated and no drug was matched to it. No corpus of transcript half-life under drug is read, and the microRNA that represses it is not identified. |
| translation | searched | 0 | a uORF is a lever and no drug was matched to it. Nothing reads ribosome profiling, so whether protein is low while message is normal - the signature that would point here - still needs the… |
| folding | searched | 0 | BUB1B's own product has 1 measured ligand(s) available, below the 20 needed to train a binding model. A screen without them would rank the library on nothing. The source that answered had… |
| degradation | searched | 0 | a drug-to-protein-abundance table is on disk and does not quantify BUB1B at all, so nothing in it bears on this protein. That is an absent measurement, not a negative one. |

*"Listed only" means the mechanism and its drug classes are named in this report and nothing was looked up for this gene. It is a menu, not a search.*

## The problem may be the protein, not the amount

**You declared that the only working copy is itself impaired, and the transcript carries a premature stop.** Everything in the promoter sections below changes how MUCH BUB1B is made. None of it changes what is made. Where the pathogenic allele yields an impaired product, more transcription yields more of the same impairment — and because transcription raises both copies, it raises the broken one too.

The routes that address the product itself are computed in this run and sit lower down:

- **Stabilisation** — 7 factor(s) assessed for a molecule that holds the protein in its working shape instead of making more of it.
- **Product stabilisation** — whether the gene's own protein can be kept from being degraded.
- **After transcription** — 6 route(s) acting on the transcript or the protein rather than the promoter. Readthrough and NMD inhibition are the ones that address a premature stop, and for a stop codon they address the problem where extra transcription does not.
- **Brake exons** — 20 exon(s) used only by decay-prone isoforms, candidates for an oligonucleotide that raises useful output rather than total output.
- **Designed modality** — a direct construct for this gene.

*This is a pointer, not a ranking. Nothing here establishes that a protein-level route works either; it establishes that it is aimed at the right thing.*

### Measured on the protein

**BUB1B is not in the drug-to-protein table, so no compound was checked against its protein.** That is an absent measurement and not a negative one: the table quantifies a few thousand proteins in one cell line, and this is not one of them.

## What the literature says moves this gene

**Most of what is recorded as raising this gene is not a medicine.** Of 59 chemicals CTD records as raising it, 37 fall into a carcinogen, tumour promoter, mitogen, inflammogen or industrial class, and 7 of them are FDA approved drugs (Choline, Estradiol, Lidocaine, Methimazole, Phenobarbital, Silver Nitrate …). That pattern is worth pausing on: when the published way to raise a gene is to damage or inflame the tissue, the gene is most likely responding to that insult rather than being regulated in a way a drug could copy. Check whether the transcriptional objective is the right objective before reading the ranking below.

| Kind of chemical | How many |
|---|---|
| unclassified | 17 |
| carcinogen or genotoxin | 9 |
| pesticide or industrial chemical | 7 |
| inflammogen or irritant | 6 |
| tumour promoter | 6 |
| hormone or mitogen | 5 |
| nutrient or dietary | 5 |
| particulate or metal salt | 4 |

*The classes come from a keyword vocabulary over the CTD chemical names, not from a toxicology database. An unclassified chemical was not examined and is not thereby safe, and a chemical can be both a medicine and a hazard - phenobarbital is an approved barbiturate and a classic rodent tumour promoter.*

### Moving the gene, or moving the cell cycle

**2 of 92 checked compound(s) move the mitotic programme along with BUB1B.** BUB1B's transcript tracks the cell-cycle phase, so a compound that moves the rest of the programme the same way is most simply read as changing how many cells are dividing rather than as acting on this gene. Their scores are halved and each says so. The bar is the 95th percentile of all compounds, against a background of 1.9% of signatures moving the programme anyway. This is a measurement over LINCS set membership, not a list of known antimitotics; a compound the files do not cover gets no verdict, which is not a pass.

| Compound | Signatures moving the programme | Background |
|---|---|---|
| clobetasolpropionate | 2 of 13 (15%) | 1.9% |
| ingenol320dibenzoate | 5 of 40 (12%) | 1.9% |

### How much of this order is fame?

**22 of 84 checked candidate(s) are recorded for most of the 11 genes measured here.** This is the control that matters most for reading the order below. Tested across these runs, a baseline that only counts how many other genes record a drug predicts CTD better (AUC 0.93) than every route in this package combined (0.53), and among drugs of equal fame the routes fall to chance. A compound recorded for every gene is telling you about the literature, not about yours. Promiscuous rows lose a little score within their tier and say so; none changes tier.

| Compound | Genes it is recorded for |
|---|---|
| Estradiol | 10 of 11 (91%) |
| Pirinixic Acid | 9 of 11 (82%) |
| Lipopolysaccharides | 8 of 11 (73%) |
| Valproic acid | 8 of 11 (73%) |
| Acrylamide | 7 of 11 (64%) |
| Dietary Sucrose | 7 of 11 (64%) |
| Nanotubes, Carbon | 7 of 11 (64%) |
| Nickel | 7 of 11 (64%) |
| Okadaic Acid | 7 of 11 (64%) |
| Phenobarbital | 7 of 11 (64%) |
| Thioacetamide | 7 of 11 (64%) |
| Tretinoin | 7 of 11 (64%) |

## Where the sources disagree

**5 ranked drugs are recorded by CTD as moving BUB1B the other way** (the goal is to move it up). A route predicts the drug should help; CTD says someone measured the opposite. A route whose own direction at its factor is not recorded is lowered to tier 4; any other is only flagged, since CTD mixes cell types and doses. Check it before testing the drug.

| Drug | Tier | Acts on | Caution |
|---|---|---|---|
| dasatinib | 3 | TFAP2C | CTD records the opposite direction |
| Valproic acid | 3 | chromatin (genome-wide) | CTD records the opposite direction |
| Fulvestrant | 4 | NR1H4 | CTD records the opposite direction and this route's own direction is unknown |
| Tretinoin | 4 | RXRA | CTD records the opposite direction and this route's own direction is unknown |
| Troglitazone | 4 | PPARG | CTD records the opposite direction |

**CTD and LINCS** both give a direction for 1 drug and agree on 1 (100%): Dasatinib.

Both sources agree that Dasatinib moves the gene the OPPOSITE way to the one this run wants. Two independent sources agreeing against a candidate is the strongest negative here, and it is worth more than the agreement rate above.

## The plate: every compound, in no order

**58 distinct compounds, in no order.** Chemical families are not folded together here and no compound is demoted for being well studied, because a plate measures what a compound does to cells and neither of those is about that. Rows are alphabetical inside each group; the groups say where the evidence places a compound, and for 8 of the 58 it places them nowhere - those rest on the gene having moved, with no mechanism attached. 6 further compound(s) were found and left off on safety grounds, listed below: that judgement is clinical, and some of them are reasonable mechanistic controls on a plate.

### transcription (50)

| Compound | Acts on | Site evidence | Routes | Found by |
|---|---|---|---|---|
| Acitretin | RARA, RXRA, RXRB, RXRG | measured | 6 | drug_scoring, network, regulators |
| Balsalazide Disodium | PPARG | — | 1 | network |
| Bazedoxifene Acetate | ESR2 | — | 1 | regulators |
| belinostat | HDAC1, HDAC3, HDAC9, SOX17, chromatin (genome-wide) | — | 5 | abundance, catalogue, regulators |
| Bexarotene | RXRA, RXRB, RXRG | measured | 6 | drug_scoring, network, regulators, screening |
| Chenodiol | NR1H4 | predicted only | 1 | drug_scoring |
| Chlorotrianisene | ESR2 | — | 1 | regulators |
| Cinacalcet Hydrochloride | CASR | — | 1 | network |
| Dextrothyroxine Sodium | THRB | — | 1 | network |
| Estramustine Phosphate Sodium | ESR2 | — | 1 | regulators |
| Etrasimod | S1PR1 | — | 2 | network |
| Etretinate | RARA, RXRA, RXRB, RXRG | measured | 6 | drug_scoring, network, regulators |
| Fulvestrant | NR1H4 | predicted only | 1 | screening |
| Givinostat Hydrochloride | HDAC1, HDAC3, HDAC9, chromatin (genome-wide) | — | 4 | catalogue, regulators |
| Hydrocortisone Acetate | NR3C1 | — | 1 | network |
| Isotretinoin | RARA | — | 1 | regulators |
| Lenalidomide Hydrochloride Monohydrate | RBX1 | — | 1 | network |
| Levothyroxine | RXRA, THRB | measured | 2 | network, screening |
| Liothyronine Sodium | THRB | — | 2 | network |
| Misoprostol | PTGER3 | — | 1 | network |
| Mometasone Furoate | NR3C1 | — | 1 | network |
| Nandrolone Phenpropionate | AR | — | 1 | network |
| Obeticholic Acid | NR1H4 | predicted only | 2 | drug_scoring, screening |
| Odevixibat | NR1H4 | predicted only | 1 | screening |
| Olsalazine Sodium | PPARG | — | 1 | network |
| Ospemifene | ESR2 | — | 1 | regulators |
| Oxandrolone | AR | — | 1 | network |
| Oxymetholone | AR | — | 1 | network |
| Ozanimod Hydrochloride | S1PR1 | — | 1 | network |
| Panobinostat | HDAC1, HDAC3, HDAC9, chromatin (genome-wide) | — | 4 | catalogue, regulators |
| Paramethasone Acetate | NR3C1 | — | 1 | network |
| Pioglitazone Hydrochloride | PPARG | — | 1 | network |
| Ponesimod | S1PR1 | — | 1 | network |
| Prednisolone Acetate | NR3C1 | — | 1 | network |
| Raloxifene Hydrochloride | ESR2 | — | 1 | regulators |
| Resmetirom | THRB | — | 1 | network |
| Romidepsin | HDAC1, HDAC3, HDAC9, chromatin (genome-wide) | — | 4 | catalogue, regulators |
| Rosiglitazone | PPARG | — | 1 | network |
| Simvastatin | NR1H4 | predicted only | 1 | screening |
| Siponimod | S1PR1 | — | 1 | network |
| Sulindac | RXRA | measured | 1 | drug_scoring |
| Tazemetostat | EZH2, chromatin (genome-wide) | — | 2 | catalogue, regulators |
| Testosterone Enanthate | AR | — | 1 | network |
| Testosterone Propionate | AR | — | 1 | network |
| Toremifene Citrate | ESR2 | — | 1 | regulators |
| Tretinoin | RARA, RXRA | measured | 2 | regulators, screening |
| Troglitazone | PPARG | — | 1 | network |
| Ursodiol | NR1H4 | predicted only | 1 | drug_scoring |
| Valproic acid | chromatin (genome-wide) | — | 1 | catalogue |
| vorinostat | HDAC1, HDAC3, NRL, chromatin (genome-wide) | — | 4 | abundance, catalogue, regulators |

### mechanism not established (8)

| Compound | Acts on | Site evidence | Routes | Found by |
|---|---|---|---|---|
| BIVALIRUDIN | KLF16 | measured | 1 | sensitivity |
| Estradiol | not stated | — | 1 | ctd |
| Lidocaine | not stated | — | 1 | ctd |
| Methimazole | not stated | — | 1 | ctd |
| PAROMOMYCIN | KLF16 | measured | 1 | sensitivity |
| Phenobarbital | not stated | — | 1 | ctd |
| PITAVASTATIN | KLF16 | measured | 1 | sensitivity |
| TIPIRACIL | THAP11 | measured | 1 | sensitivity |

**No compound is placed at acts on the gene or its product, splicing and isoform choice, mRNA stability, translation, folding, degradation.** Those steps were searched and returned annotations rather than compounds, so they contribute nothing to a plate however the groups above are weighted. Giving them equal weight here would mean weighting empty sets.

### Found and left off (6)

The reason in each case is clinical. On a plate several of these are reasonable mechanistic controls, and adding one back is a decision this table exists to inform.

| Compound | Why it is not on the list |
|---|---|
| Azacitidine | works by damaging DNA or the spindle, which is the defect being treated |
| dasatinib | ETS2 is a proto-oncogene and this route raises it |
| Decitabine | works by damaging DNA or the spindle, which is the defect being treated |
| mitoxantrone | works by damaging DNA or the spindle, which is the defect being treated |
| Paclitaxel | works by damaging DNA or the spindle, which is the defect being treated |
| thioridazine hydrochloride | FOXM1 is a proto-oncogene and this route raises it |

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

> **⚠ Read as a reading order, not as evidence.** The run's evidence verdict is *degraded*: at least one pillar of the analysis is missing. See the evidence table above.

| # | Medication | What it is | Approved to treat | Acts on | BBB? | Tier | Grade | Score |
|---|---|---|---|---|---|---|---|---|
| 1 | **Obeticholic Acid** | binds NR1H4 (measured), and NR1H4 is predicted to act on this promoter - so the drug would reach your gene THROUGH NR1H4, not by touching it directly | not established here | NR1H4 | no | 3 | D | 51.5 |
| 2 | **Odevixibat** | binds NR1H4 (measured), and NR1H4 is predicted to act on this promoter - so the drug would reach your gene THROUGH NR1H4, not by touching it directly | not established here | NR1H4 | — | 3 | — | 49.3 |
| 3 | **Bexarotene** | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma | RXRA | yes? | 3 | B | 48.7 |
| 4 | **dasatinib** | raises TFAP2C transcript in 7 profile(s) | not established here | TFAP2C | — | 3 | — | 45.0 |
| 5 | **Chenodiol** | ligand at the transcription factor | not established here | NR1H4 | no | 3 | B | 44.4 |
| 6 | **Ursodiol** | ligand at the transcription factor | not established here | NR1H4 | no | 3 | B | 44.4 |
| 7 | **Acitretin** | retinoid | severe psoriasis | RXRA | yes? | 3 | B | 40.7 |
| 8 | **Etretinate** | retinoid | severe psoriasis - WITHDRAWN, teratogenic | RXRA | yes? | 3 | C | 40.7 |
| 9 | **Levothyroxine** | synthetic T4 | hypothyroidism | RXRA | — | 3 | — | 40.1 |
| 10 | **Sulindac** | ligand at the transcription factor | not established here | RXRA | — | 3 | B | 39.7 |

**Every row is off-label.** The drug was approved for what the fourth column says, and is listed here for something else. That column is the baseline this proposal departs from, and the distance between it and your gene is part of what you are judging - a topical acne gel and an infused cytotoxic can sit on adjacent rows with similar scores.

*Ordered by evidence tier first, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. **BBB?** is not part of the ordering — it is an applicability filter you apply yourself, since whether brain exposure is wanted or unwanted depends on where your gene needs to move. See the note under the graded table for what a yes does and does not establish.*

**Engage** and **promoter** below are the two links that come from evidence: the drug engaging its target at a dose that can be reached, and the factor acting on this promoter. They differ between drugs and are 0 to 1 heuristics, not probabilities. Order on tier and score.

### Why each one

Every route below is **evidence tier 3** — a human curator annotated this mechanism.

#### Through **NR1H4** — 4 candidates

- The factor has no activator/repressor annotation, so the direction of its effect on your gene is assumed, not known.
- Binding: 1 predicted site(s), the best scoring 0.882 -413 bp from the TSS; no ChIP peak here, but the ChIP table never covers this factor, so its binding was not measured either way.
- Collateral: the factor controls roughly 300 gene(s), which all move with yours; clinical burden **prohibitive**.

1. **Obeticholic Acid** — grade D, engage 0.74, promoter 0.51, knockout: — — binds NR1H4 (measured), and NR1H4 is predicted to act on this promoter - so the drug would reach your gene THROUGH NR1H4, not by touching it directly
   - ⚠ half-measured: the Obeticholic Acid-NR1H4 interaction has assay data behind it (pChEMBL 8.0), so the pharmacology is de-risked. What is untested is whether moving NR1H4 moves your gene - the one step a binding assay cannot tell you
   - *Would falsify it:* Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.
2. **Odevixibat** — binds NR1H4 (measured), and NR1H4 is predicted to act on this promoter - so the drug would reach your gene THROUGH NR1H4, not by touching it directly
   - ⚠ half-measured: the Odevixibat-NR1H4 interaction has assay data behind it (pChEMBL 9.8), so the pharmacology is de-risked. What is untested is whether moving NR1H4 moves your gene - the one step a binding assay cannot tell you
3. **Chenodiol** — grade B, engage 0.74, promoter 0.51, knockout: — — ligand at the transcription factor
   - ⚠ measured: pChEMBL 6.65 against NR1H4 (n=3183), with a recorded action
   - *Would falsify it:* Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.
4. **Ursodiol** — grade B, engage 0.74, promoter 0.51, knockout: — — ligand at the transcription factor
   - ⚠ measured: pChEMBL 6.65 against NR1H4 (n=3183), with a recorded action
   - *Would falsify it:* Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.

#### Through **RXRA** — 5 candidates

- The factor has no activator/repressor annotation, so the direction of its effect on your gene is assumed, not known.
- Binding: 2 predicted site(s), the best scoring 0.907 -413 bp from the TSS; **1 ChIP peak(s) confirm it actually binds here** — the strongest link in this chain.
- Collateral: the factor controls roughly 1200 gene(s), which all move with yours; clinical burden **high** (retinoid).
- Family: RXRA shares its 9-cis retinoic acid pocket with RXRB, RXRG, so a ligand here will tend to engage them too, and no subtype-selective agent is recorded for this family.

5. **Bexarotene** — grade B, engage 0.70, promoter 0.79, knockout: not among top changes, 0 of 2 — binds RXRA (measured), and RXRA is predicted to act on this promoter - so the drug would reach your gene THROUGH RXRA, not by touching it directly
   - ⚠ half-measured: the Bexarotene-RXRA interaction has assay data behind it (pChEMBL 6.82), so the pharmacology is de-risked. What is untested is whether moving RXRA moves your gene - the one step a binding assay cannot tell you
   - *Would falsify it:* Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.
6. **Acitretin** — grade B, engage 0.70, promoter 0.79, knockout: not among top changes, 0 of 2 — ligand at the transcription factor
   - ⚠ measured: pChEMBL 6.37 against RXRA (n=675), with a recorded action
   - *Would falsify it:* Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.
7. **Etretinate** — grade C, engage 0.70, promoter 0.79, knockout: not among top changes, 0 of 2 — ligand at the transcription factor
   - ⚠ measured: pChEMBL 6.37 against RXRA (n=675), with a recorded action
   - *Would falsify it:* Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.
8. **Levothyroxine** — binds RXRA (measured), and RXRA is predicted to act on this promoter - so the drug would reach your gene THROUGH RXRA, not by touching it directly
   - ⚠ half-measured: the Levothyroxine-RXRA interaction has assay data behind it (pChEMBL 6.03), so the pharmacology is de-risked. What is untested is whether moving RXRA moves your gene - the one step a binding assay cannot tell you
9. **Sulindac** — grade B, engage 0.70, promoter 0.79, knockout: not among top changes, 0 of 2 — ligand at the transcription factor
   - ⚠ measured: pChEMBL 6.37 against RXRA (n=675), with a recorded action
   - *Would falsify it:* Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.

#### Through **TFAP2C**

- The factor has no activator/repressor annotation, so the direction of its effect on your gene is assumed, not known.

10. **dasatinib**
   - ⚠ measured against TFAP2C's own transcript, not against your gene. That TFAP2C controls your gene is the inference this route rests on

*Not medical advice. "Approved" means approved for its own indication — using one of these to move an unrelated gene is off-label repurposing with the full systemic pharmacology attached. The experimental sections at the end of this report are how you would actually test any of it.*

## Promoter of BUB1B (ENST00000287598)

- Transcript chosen: **distinct TSS (representing 2 within 500 bp)**
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

#### Recorded regulators with no approved drug

These are not motif matches. Each is named by a curated table or by a knockout experiment as regulating this gene, and none has an approved drug that moves it the way this run needs. They generate no route and are listed so the regulator section is not read as the full set of regulators. This is the list a chemistry effort would aim at.

| Factor | Effect on the gene | Would need to be | Recorded by |
|---|---|---|---|
| ETS2 | raises it | activated | knockout data |
| FLI1 | raises it | activated | knockout data |
| FOXA1 | raises it | activated | knockout data |
| FOXM1 | raises it | activated | CollecTRI |
| FOXP1 | lowers it | inhibited | knockout data |
| MYC | raises it | activated | CollecTRI |
| NRL | lowers it | inhibited | knockout data |
| SOX17 | lowers it | inhibited | knockout data |
| SRF | lowers it | inhibited | knockout data |
| TFAP2C | raises it | activated | knockout data |
| TP53 | raises it | activated | CollecTRI |
| ZNF143 | raises it | activated | CollecTRI |

A further **207** factor(s) reach the gene only through one of the regulators above, so each is a longer way round with a sign that is the product of two edges. They are in `results.json` rather than here, because at this distance the list is a property of the network and not of this gene.

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

### Drugs that behave like losing the factor

Every other route in this report starts from something somebody wrote down - a curated ligand, a motif, a structure, a signed edge. None of them can surface a drug nobody has connected to the factor. This one reads a measured phenotype instead: if a compound hurts the same cell lines that knocking the factor out hurts, it behaves as though it removes the factor.

| Factor | Drug | Partial r | Raw r | Lines | q | Behaves as |
|---|---|---|---|---|---|---|
| KLF16 | **PITAVASTATIN** | 0.2436 | 0.203 | 463 | 0.00039 | phenocopies |
| KLF16 | **BIVALIRUDIN** | 0.2179 | 0.2074 | 451 | 0.006996 | phenocopies |
| THAP11 | **TIPIRACIL** | -0.2073 | -0.1459 | 444 | 0.036841 | opposes |
| KLF16 | **PAROMOMYCIN** | 0.2022 | 0.2015 | 451 | 0.014411 | phenocopies |

**Partial r** controls for how easily each cell line dies under the whole panel. Some lines die under almost anything and those same lines suffer most knockouts, which puts a correlation between almost any drug and almost any gene. The control moves a number either way - it removes a shared component that can inflate a correlation or mask one - so where the two columns differ a lot the raw value was mostly about the cell lines. Estimates use at least 25 shared lines.

**This is not a target.** A correlated drug may act on the factor, on anything else the same lines need, or on a pathway the factor sits in, and the first is the least likely of the three for any one hit. Nothing here records any of them binding anything: the finding is a reason to run a binding assay, not a result from one.

**A proliferation readout.** Both panels score growth of transformed cells in a dish, so a factor whose role is in a differentiated or post-mitotic cell looks silent here for reasons that have nothing to do with its importance.

### Structure-based docking

Docking asks whether a molecule fits the pocket. It says nothing about direction: a compound that binds a factor may activate it, block it, or do neither to transcription of this gene.

| Pocket | Docked | Failed | ROC AUC | EF 1% | Validated? |
|---|---|---|---|---|---|
| BUB1B | 0 | 0 | not tested | not tested | **no** |

Validated means the same scoring function separated known actives from decoys in this pocket, at AUC 0.70 or better and 2x or better enrichment in the top 1%. Below that the geometry may still be right and the score is not evidence.

**1 pocket did not validate** (BUB1B). Their hits are listed for the record and should not be ranked against hits from a pocket that did.

### Every route that acts on this gene

105 routes across 7 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

**Why established binders sit at tier 2, just under the direct agents.** A drug already measured to bind one of your factors has half its chain settled: the drug-to-receptor step has assay data behind it, so the pharmacology is de-risked and only the receptor-to-gene step is a hypothesis. A curated tier-3 ligand is annotation on *both* halves. But they stay below tier 1, because a direct agent acts on your gene itself and skips the receptor-to-gene inference entirely — a stronger claim however well-characterised the binder is. The cost of the indirect route is in the row: everything else that factor regulates moves with your gene.

- 24 × curated — mechanism annotated by a human curator
- 74 × predicted — a model says so
- 7 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 3 | TF: orthosteric ligand | Obeticholic Acid | NR1H4 | 1 | binds NR1H4 (measured), and NR1H4 is predicte | 0.766 | 0.5 | 0.65 | 51.5 |
| 3 | TF: orthosteric ligand | Odevixibat | NR1H4 | 1 | binds NR1H4 (measured), and NR1H4 is predicte | 0.698 | 0.5 | 0.65 | 49.3 |
| 3 | TF: orthosteric ligand | Bexarotene | RXRA | 1 | binds RXRA (measured), and RXRA is predicted  | 0.678 | 0.5 | 0.65 | 48.7 |
| 3 | TF: abundance change | dasatinib | TFAP2C | 1 | raises TFAP2C transcript in 7 profile(s) | 0.7 | 0.45 | 0.5 | 45.0 |
| 3 | TF: orthosteric ligand | Chenodiol | NR1H4 | 1 | ligand at the transcription factor | 0.614 | 0.337 | 0.765 | 44.4 |
| 3 | TF: orthosteric ligand | Ursodiol | NR1H4 | 1 | ligand at the transcription factor | 0.614 | 0.337 | 0.765 | 44.4 |
| 3 | TF: orthosteric ligand | Obeticholic Acid | NR1H4 | 1 | ligand at the transcription factor | 0.614 | 0.337 | 0.765 | 44.4 |
| 3 | TF: orthosteric ligand | Acitretin | RXRA | 1 | ligand at the transcription factor | 0.627 | 0.21 | 0.737 | 40.7 |
| 3 | TF: orthosteric ligand | Bexarotene | RXRA | 1 | ligand at the transcription factor | 0.627 | 0.21 | 0.737 | 40.7 |
| 3 | TF: orthosteric ligand | Etretinate | RXRA | 1 | ligand at the transcription factor | 0.627 | 0.21 | 0.737 | 40.7 |
| 3 | TF: orthosteric ligand | Levothyroxine | RXRA | 1 | binds RXRA (measured), and RXRA is predicted  | 0.409 | 0.5 | 0.65 | 40.1 |
| 3 | TF: orthosteric ligand | Sulindac | RXRA | 1 | ligand at the transcription factor | 0.597 | 0.21 | 0.737 | 39.7 |
| 3 | TF: orthosteric ligand | Simvastatin | NR1H4 | 1 | binds NR1H4 (measured), and NR1H4 is predicte | 0.388 | 0.5 | 0.65 | 39.4 |
| 3 | TF: abundance change | dasatinib | FOXA1 | 1 | raises FOXA1 transcript in 2 profile(s) | 0.2 | 0.45 | 0.5 | 29.0 |
| 3 | TF: abundance change | vorinostat | NRL | 1 | lowers NRL transcript in 2 profile(s) | 0.2 | 0.45 | 0.5 | 29.0 |
| 3 | TF: abundance change | dasatinib | SOX17 | 1 | lowers SOX17 transcript in 2 profile(s) | 0.2 | 0.45 | 0.5 | 29.0 |
| 3 | TF: abundance change | belinostat | SOX17 | 1 | lowers SOX17 transcript in 2 profile(s) | 0.2 | 0.45 | 0.5 | 29.0 |
| 3 | chromatin-level | Vorinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Romidepsin | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Belinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Panobinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Tazemetostat | chromatin (genome-wide) | 3 | EZH2 inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Givinostat Hydrochloride | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Valproic acid | chromatin (genome-wide) | 3 | weak HDAC inhibitor | 0.7 | 0.05 | 0.45 | 19.9 |
| 4 | TF: orthosteric ligand | Fulvestrant | NR1H4 | 1 | binds NR1H4 (measured), and NR1H4 is predicte | 0.372 | 0.5 | 0.65 | 38.9 |
| 4 | TF: orthosteric ligand | Tretinoin | RXRA | 1 | binds RXRA (measured), and RXRA is predicted  | 0.46 | 0.5 | 0.65 | 35.4 |
| 4 | TF: known regulator | Raloxifene Hydrochloride | ESR2 | 2 | inhibit ESR2, which regulates FOXM1, MYC, a r | 0.4 | 0.5 | 0.6 | 32.3 |
| 4 | TF: known regulator | Estramustine Phosphate Sodium | ESR2 | 2 | inhibit ESR2, which regulates FOXM1, MYC, a r | 0.4 | 0.5 | 0.6 | 32.3 |
| 4 | TF: known regulator | Toremifene Citrate | ESR2 | 2 | inhibit ESR2, which regulates FOXM1, MYC, a r | 0.4 | 0.5 | 0.6 | 32.3 |
| 4 | TF: known regulator | Chlorotrianisene | ESR2 | 2 | inhibit ESR2, which regulates FOXM1, MYC, a r | 0.4 | 0.5 | 0.6 | 32.3 |

Best in each category:

- **TF: orthosteric ligand** — Obeticholic Acid (tier 3, score 51.5)
- **TF: abundance change** — dasatinib (tier 3, score 45.0)
- **chromatin-level** — Vorinostat (tier 3, score 23.4)
- **TF: known regulator** — Raloxifene Hydrochloride (tier 4, score 32.3)
- **upstream signalling** — Balsalazide Disodium (tier 4, score 29.8)
- **literature: reported effect** — Lidocaine (tier 4, score 29.3)
- **measured phenocopy** — PITAVASTATIN (tier 4, score 29.0)

*No route found in: direct: gene or transcript, direct: gene product, direct: predicted product binder, promoter DNA binder, TF: inferred activity change, protein abundance change, TF: predicted pocket binder, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, TF: degradation blocker, off-pathway pocket binder, empirical signature, virtual screen hypothesis. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

| Category | What it changes | Your drugs in it | Assay that sees it |
|---|---|---|---|
| TF: orthosteric ligand | **abundance** | Obeticholic Acid, Odevixibat, Bexarotene, Chenodiol … +8 | RT-qPCR then western/ELISA |
| TF: abundance change | **abundance** | dasatinib, vorinostat, belinostat, mitoxantrone … +1 | RT-qPCR and a western for the FACTOR first, then RT-qPCR for your gene |
| chromatin-level | **abundance, genome-wide** | Vorinostat, Romidepsin, Belinostat, Panobinostat … +5 | RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened |
| TF: known regulator | **abundance** | Raloxifene Hydrochloride, Estramustine Phosphate Sodium, Toremifene Citrate, Chlorotrianisene … +13 | RT-qPCR for your gene after treatment; check the factor's own activity (a reporter or a known target) moved first |
| upstream signalling | **abundance** | Balsalazide Disodium, Olsalazine Sodium, Rosiglitazone, Pioglitazone Hydrochloride … +26 | RT-qPCR then western/ELISA; confirm the pathway node moved too, or you cannot attribute the effect |
| literature: reported effect | **abundance (reported)** | Lidocaine, Methimazole, Estradiol, Phenobarbital … +1 | RT-qPCR in your cells; read the paper for the dose, cell type and species first |
| measured phenocopy | **abundance** | PITAVASTATIN, BIVALIRUDIN, TIPIRACIL, PAROMOMYCIN | a binding assay against the factor, then a reporter, then RT-qPCR for your gene - the correlation names a candidate and establishes none of those steps |

*Every category present here works through transcription, so all of them move your gene's product **abundance** and an RT-qPCR is the right first readout. The one exception in the taxonomy — a drug on the gene product itself, which changes activity only — did not appear in this run.*

**TF: orthosteric ligand** (14 found)

- *Yours:* Obeticholic Acid, Odevixibat, Bexarotene, Chenodiol, Ursodiol, Acitretin, Etretinate, Levothyroxine … +4 more
- Binds the transcription factor's natural ligand pocket.
- *Why:* The best-characterised route, and where nearly all approved TF pharmacology lives.
- *Cost:* Only ~48 human TFs (the nuclear receptors) have such a pocket, and the factor's whole regulon moves with your gene.
- *Differs:* Clamps the receptor on or off regardless of endogenous signalling — contrast the allosteric and ligand-axis routes, which scale what is already happening.

**TF: abundance change** (9 found)

- *Yours:* dasatinib, vorinostat, belinostat, mitoxantrone, thioridazine hydrochloride
- A compound measured to change how much of the factor's own transcript exists, rather than what the factor does.
- *Why:* The only route here that needs the factor to be expressed rather than druggable, so it reaches factors with no ligand, no structure and no curated upstream kinase - which is most of them.
- *Cost:* More factor is not more factor ACTIVITY. For a ligand-activated receptor the two are nearly unrelated, and transcript is not protein. Broad agents that shut down transcription move every factor and are flagged.
- *Differs:* Unlike 'empirical signature', which is a compound measured against YOUR GENE, this was measured against the FACTOR - one further inferential step from what you care about.

**chromatin-level** (9 found)

- *Yours:* Vorinostat, Romidepsin, Belinostat, Panobinostat, Tazemetostat, Givinostat Hydrochloride, Valproic acid, Azacitidine … +1 more
- HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.
- *Why:* Reliable positive control that your reporter *can* be induced.
- *Cost:* Specificity ~0.05. These move thousands of promoters. Almost never an acceptable therapeutic route on that basis alone.
- *Differs:* Not targeted at your gene in any sense; it is in the table so you can see what a non-specific option scores.

**TF: known regulator** (28 found)

- *Yours:* Raloxifene Hydrochloride, Estramustine Phosphate Sodium, Toremifene Citrate, Chlorotrianisene, Ospemifene, Bazedoxifene Acetate, Tazemetostat Hydrobromide, Romidepsin … +9 more
- Drugs that act on a factor CollecTRI or perturbation data records as raising or lowering your gene, whether or not the motif scan ranked that factor.
- *Why:* Reaches factors the motif scan misses, such as those that bind through a partner or a distant enhancer.
- *Cost:* The link comes from a database or a knockout in other cells, so it may not hold in yours. The drug's effect on the factor is taken from its recorded action.
- *Differs:* Starts from the literature on your gene, not from its promoter sequence. CollecTRI and CTD share sources, so agreement between them is not independent.

**upstream signalling** (36 found)

- *Yours:* Balsalazide Disodium, Olsalazine Sodium, Rosiglitazone, Pioglitazone Hydrochloride, Troglitazone, Cinacalcet Hydrochloride, Liothyronine Sodium, Liothyronine … +22 more
- Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- *Why:* Turns 'this factor has no drugs' into a real list. Most TFs are only reachable this way.
- *Cost:* Specificity falls with every step: a node three hops up reaches your gene through one branch and thousands of other genes through the rest.
- *Differs:* The only category where the drug never touches your gene or its TF — read the printed path and judge it yourself.

**literature: reported effect** (5 found)

- *Yours:* Lidocaine, Methimazole, Estradiol, Phenobarbital, Paclitaxel
- A chemical CTD records as raising or lowering this gene's mRNA, with no mechanism stated.
- *Why:* It is the only route that can name a drug nothing here can explain, which includes most of what CTD holds: stress, feedback and toxicity effects.
- *Cost:* Most records rest on one paper, often at a toxicology dose in one cell type. Across 36 genes only about 3% of them were reproduced in LINCS, and more papers did not raise that. Famous drugs are over-represented.
- *Differs:* Every other category predicts an effect from a mechanism. This reports an effect somebody saw, and cannot suggest a drug nobody has tested on this gene.

**measured phenocopy** (4 found)

- *Yours:* PITAVASTATIN, BIVALIRUDIN, TIPIRACIL, PAROMOMYCIN
- A compound whose effect on cell viability across the DepMap panel tracks what happens when this factor is knocked out, so it behaves as though it removes the factor.
- *Why:* It needs no annotation of any kind. Every other route here starts from something a curator, a motif or a structure already asserted, so none of them can surface a drug nobody has connected to the factor. This one can, because it reads a measured phenotype.
- *Cost:* Correlation across cell lines, not target engagement. A hit may act on the factor, on anything else the same lines need, or on a pathway the factor sits in, and the first is the least likely of the three. Both panels score proliferation of transformed cells, so a factor that matters in a differentiated cell leaves no trace.
- *Differs:* Unlike 'empirical signature', nothing was measured about your gene or about transcription at all - the readout is whether cells live.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **direct: gene product** — Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- **direct: predicted product binder** — A compound predicted to bind this gene's own protein and not recorded as acting on it - either because its chemistry resembles the protein's known ligands, or because it docks into the protein's pocket.
- **promoter DNA binder** — A small molecule with a documented preference for a sequence feature this promoter carries - a quadruplex, a GC-box, an AT tract - proposed because the feature is enriched against shuffles of this same sequence.
- **TF: inferred activity change** — A compound whose own LINCS signature moves this factor's known targets, found without any record that it binds the factor.
- **protein abundance change** — A compound measured to change how much of your gene's own PROTEIN exists, in a proteomics experiment under drug
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
| Decitabine | DNA methyltransferase inhibitor | myelodysplastic syndrome |
| Dextrothyroxine Sodium | T4 stereoisomer | formerly hyperlipidaemia - withdrawn for cardiac mortality |
| Estradiol | oestrogen receptor agonist | menopausal hormone replacement |
| Etretinate | retinoid | severe psoriasis - WITHDRAWN, teratogenic |
| Fulvestrant | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Hydrocortisone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Isotretinoin | retinoid | severe nodular acne |
| Lenalidomide Hydrochloride Monohydrate | cereblon modulator | multiple myeloma, MDS |
| Levothyroxine | synthetic T4 | hypothyroidism |
| Liothyronine Sodium | synthetic T3 | hypothyroidism, myxoedema coma |
| Mometasone Furoate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Nandrolone Phenpropionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Oxandrolone | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Pioglitazone Hydrochloride | PPARgamma agonist | type 2 diabetes |
| Prednisolone Acetate | glucocorticoid class | inflammation and immune suppression - asthma, autoimmune disease, transplant rejection, some leukaemias |
| Raloxifene Hydrochloride | selective oestrogen receptor modulator | osteoporosis; breast cancer risk reduction |
| Romidepsin | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma |
| Rosiglitazone | PPARgamma agonist (thiazolidinedione) | type 2 diabetes; restricted for cardiovascular risk |
| Simvastatin | HMG-CoA reductase inhibitor | hypercholesterolaemia |
| Testosterone Enanthate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Testosterone Propionate | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers |
| Tretinoin | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically |
| mitoxantrone | topoisomerase II inhibitor, DNA intercalator | acute myeloid leukaemia, prostate cancer, multiple sclerosis |
| vorinostat | histone deacetylase inhibitor | cutaneous T-cell lymphoma |

**No description established for 37**: BIVALIRUDIN, Balsalazide Disodium, Bazedoxifene Acetate, Chenodiol, Chlorotrianisene, Cinacalcet Hydrochloride, Estramustine Phosphate Sodium, Etrasimod, Givinostat Hydrochloride, Lidocaine, Methimazole, Misoprostol, Obeticholic Acid, Odevixibat, Olsalazine Sodium, Ospemifene, Oxymetholone, Ozanimod Hydrochloride, PAROMOMYCIN, PITAVASTATIN, Paclitaxel, Panobinostat, Paramethasone Acetate, Phenobarbital, Ponesimod…. That means this package could not find one offline, not that the compound has no approved use. Several are environmental agents and research chemicals from the CTD tables, which have no indication because they are not drugs.

Descriptions come from: 8 class, 19 curated, 37 unknown. None of these is a label. The curated text is a summary written for orientation, ChEMBL's indication list includes uses studied but never approved, and a class line gives the kind of drug rather than its indication. Check the prescribing information before acting on any row.

### Graded candidates

Combined = geometric mean of TF confidence and drug quality. A weak TF prediction caps the grade no matter how good the drug.

| Grade | Combined | Drug | BBB? | TF | TF conf. | Drug qual. | Potency | Protein select. | Gene select. | Evidence | Flags |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **B** | 66.5 | Acitretin | yes? | RXRA | 0.692 | 0.605 | 0.343 | 0.5 | 0.42 | 1.0 | RXRA controls ~1200 genes; expect broad transcriptional effects; direc |
| **B** | 66.5 | Bexarotene | yes? | RXRA | 0.692 | 0.605 | 0.343 | 0.5 | 0.42 | 1.0 | RXRA controls ~1200 genes; expect broad transcriptional effects; direc |
| **B** | 66.2 | Sulindac | — | RXRA | 0.692 | 0.596 | 0.343 | 0.5 | 0.42 | 1.0 | RXRA controls ~1200 genes; expect broad transcriptional effects; direc |
| **B** | 55.8 | Chenodiol | no | NR1H4 | 0.513 | 0.681 | 0.413 | 0.5 | 0.673 | 1.0 | direction not independently predicted |
| **B** | 55.8 | Ursodiol | no | NR1H4 | 0.513 | 0.681 | 0.413 | 0.5 | 0.673 | 1.0 | direction not independently predicted |
| **C** | 41.0 | Etretinate | yes? | RXRA | 0.692 | 0.121 | 0.343 | 0.5 | 0.42 | 1.0 | RXRA controls ~1200 genes; expect broad transcriptional effects; WITHD |
| **D** | 34.5 | Obeticholic Acid | no | NR1H4 | 0.513 | 0.136 | 0.413 | 0.5 | 0.673 | 1.0 | WITHDRAWN from market; direction not independently predicted |

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

A **yes** predicts *passive diffusion only*, and is the weaker of the two answers. Active efflux by P-glycoprotein at the endothelium pumps an otherwise suitable molecule straight back out; known substrates are checked here and a hit overrides a pass outright.

**A question mark means undecided** — the score landed within noise of the threshold, and 3 of 7 row(s) here did. This is not pedantry. **Atenolol clears a complete five-parameter MPO by 0.04 and is a textbook poor penetrant**, so a margin that thin carries no information. Read a question-marked cell as "look this one up", never as the answer it leans toward.

**— is not a no.** 1 of 7 row(s) returned no physicochemical properties — common for DGIdb- and Open Targets-sourced rows, which carry a name and a mechanism but no structure. Nothing was measured and nothing is being claimed; look the compound up before excluding it on this column.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

### Upstream routes to factors with no direct drug

Most TFs have no approved ligand. These paths walk upstream through signed causal edges — kinases, phosphatases, receptors — to proteins that *are* druggable. Signs multiply along the path, so an inhibitor of an inhibitor counts as activating.

Distance buys druggability and spends selectivity: a node three hops up reaches your gene through one branch and everything else through the others. Read the fan-out column.

**SP3** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| NOD2 | 3 | 9 | NOD2 --> MAPK14 —\| MAPK3 --> SP3 | inhibit | — | 16.3 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 9 downstream branches
- ⚠ path passes through hub(s): MAPK14, MAPK3

**KLF16** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| RBX1 | 3 | 7 | RBX1 —\| ESR1 --> SRC --> KLF16 | inhibit | Lenalidomide Hydrochloride Monohydrate | 16.3 |
| NR3C1 | 3 | 16 | NR3C1 —\| MAPK8 --> SRC --> KLF16 | inhibit | Hydrocortisone Acetate, Prednisolone Acetate, Mometasone Furoate | 13.8 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 16 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 7 downstream branches
- ⚠ path passes through hub(s): SRC

**FLI1** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| PPARG | 2 | 20 | PPARG —\| HDAC1 --> FLI1 | inhibit | Balsalazide Disodium, Olsalazine Sodium, Rosiglitazone | 24.4 |
| CASR | 3 | 2 | CASR --> FYN --> PRKCD —\| FLI1 | inhibit | Cinacalcet Hydrochloride | 18.3 |
| THRB | 3 | 6 | THRB --> PPARG —\| HDAC1 --> FLI1 | inhibit | Liothyronine Sodium, Liothyronine, Levothyroxine Sodium | 16.7 |
| RXRG | 3 | 7 | RXRG --> PPARG —\| HDAC1 --> FLI1 | inhibit | Bexarotene, Acitretin, Etretinate | 16.4 |
| PTGER3 | 3 | 9 | PTGER3 --> SRC --> PRKCD —\| FLI1 | inhibit | Misoprostol | 15.7 |
| S1PR1 | 3 | 12 | S1PR1 --> LCK --> PRKCD —\| FLI1 | inhibit | Ponesimod, Siponimod, Etrasimod | 14.9 |
| RXRB | 3 | 14 | RXRB --> PPARG —\| HDAC1 --> FLI1 | inhibit | Bexarotene, Acitretin, Etretinate | 14.3 |
| RXRA | 3 | 22 | RXRA —\| CTNNB1 --> EP300 --> FLI1 | inhibit | Bexarotene, Acitretin, Etretinate | 12.4 |
| AR | 3 | 38 | AR --> NKX3-1 —\| HDAC1 --> FLI1 | inhibit | Testosterone Propionate, Testosterone Enanthate, Oxandrolone | 9.9 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 12 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 14 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 2 downstream branches

### Mechanism options: orthosteric vs allosteric

Allosteric sites are not constrained to bind the endogenous ligand, so they diverge across a receptor family where the orthosteric pocket cannot. That usually buys subtype selectivity — at the cost of a bounded maximal effect and a hard dependence on endogenous tone.

**RXRA**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| orthosteric agonist | 4 | full | typically lower | no | no | Bexarotene, Acitretin |

*No approved allosteric modulator at this target; only orthosteric pharmacology, with the usual subtype-selectivity problem.*

**Before pursuing any allosteric option, verify that the pathway has endogenous tone in your cell type.** A PAM with no ligand to modulate does nothing at all, and that is a measurement, not a prediction — nothing in this package can check it for you.

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

#### No stabilisation route for 7 factor(s)

KLF16, KLF3, OSR2, SP3, THAP11, ZNF143, ZNF76

These have no ligand pocket, no direct drug, and **nothing curated against them or their degradation machinery**. For these factors stabilisation is a research programme, not a route you can take: making it real needs the E3 ligase that degrades the factor *and* a ligand for that ligase, and neither is published here.

The mechanism class is real — it has shipped approved drugs, listed once below — but those act on **other proteins entirely**, and none of them is evidence that any factor above is tractable.

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Cell-type reach

Three things must overlap: where the drug goes, where the gene is expressed, and where the transcription factor is. Any empty intersection kills the route regardless of motif quality or potency.

| TF | Co-expressed with gene | Shared tissues |
|---|---|---|
| RXRA | yes (27) | adipose tissue, appendix, bone marrow, breast, cervix, uteri |
| OSR2 | yes (21) | adipose tissue, appendix, breast, cervix, uterine, colon, en |
| THAP11 | yes (28) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |
| KLF3 | yes (28) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |
| SP3 | yes (28) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |
| KLF16 | yes (28) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |
| ZNF76 | yes (28) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |
| ZNF143 | yes (28) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |
| FLI1 | yes (28) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |
| TFAP4 | yes (27) | adipose tissue, appendix, b-cells, bone marrow, breast, cerv |

### Human genetic evidence

Human genetic variation in BUB1B associates with: mosaic variegated aneuploidy syndrome, mosaic variegated aneuploidy syndrome 1, premature chromatid separation trait. Read this as a preview of what modulating the gene may do, including the adverse effects — this is the reasoning that anticipated PCSK9 inhibition years before the trials.

| Phenotype | Genetic score | Overall |
|---|---|---|
| mosaic variegated aneuploidy syndrome | 0.821 | 0.815 |
| mosaic variegated aneuploidy syndrome 1 | 0.941 | 0.794 |
| premature chromatid separation trait | 0.771 | 0.631 |
| colorectal cancer | 0.662 | 0.503 |
| hereditary disease | 0.776 | 0.472 |
| colon carcinoma | 0.294 | 0.38 |

#### Does variation in the factor actually move this gene? (trans-eQTL)

The strongest evidence available anywhere in this report. A motif says a factor *could* bind; a trans-eQTL says natural genetic variation near it **does** change your gene's expression, in real human tissue, with exposure randomised by meiosis.

**No factor can be supported this way for BUB1B.** BUB1B has no significant trans-eQTL in eQTLGen at all, so the test returns nothing for every factor here. That is a property of the dataset's power and of the tissue (whole blood), not evidence that BUB1B is unregulated.

*No factor tested (RXRA, OSR2, THAP11, KLF3, SP3, KLF16) has trans-eQTL support for this gene.*

*eQTLGen is whole blood at n=31,684, FDR<0.05. A hit is strong causal evidence. **Absence is weak** — trans effects are small and detecting them needs enormous power, so a null overwhelmingly means underpowered rather than unregulated, and a relationship operating outside blood would not appear at all.*

### Routes a promoter analysis cannot see

Everything above acts on transcription. These act after it, around it, or on a different gene entirely. For some genes that is where the whole problem is: if the transcript is destroyed, the allele is silenced, or the protein is degraded too fast, more transcription does not help and may not be possible.

#### Could a relative do the job instead?

Whether a paralog could substitute for this gene is unanswered. It is worth asking explicitly, because the question 'restore this function' has more answers than 'raise this transcript'.

- BUB1B is not in the curated redundancy table, which is short and certainly incomplete. Absence here is not evidence that no paralog substitutes — it means nobody has written one down in this file.

#### After transcription

A promoter-level analysis stops at transcription. These act after it, and for some genes that is where the whole problem is: if the transcript is destroyed or mistranslated, doubling transcription doubles the amount that is destroyed. **Your gene has a premature stop codon**, so readthrough and NMD inhibition are listed first — they address that allele directly, and no amount of extra transcription does. 1 of 6 mechanism(s) here have an FDA-approved agent. The rest are listed so that an empty transcriptional result is not mistaken for 'nothing can be done' — they are mechanisms, not treatments.

| Mechanism | Raises | Status | Applies when |
|---|---|---|---|
| nonsense readthrough | functional protein | approved elsewhere | the gene carries a premature termination codon (PTC) |
| nonsense-mediated decay (NMD) inhibition | transcript | investigational | the gene carries a PTC and the transcript is being degrade |
| splice correction | functional protein | approved | a splice-site or splice-regulatory variant, a poison exon, |
| transcript stabilisation | transcript | investigational | the transcript is short-lived and the promoter is already  |
| microRNA de-repression | transcript | investigational | a specific miRNA is a documented repressor of your gene in |
| translational enhancement | protein | investigational | protein is low while transcript is normal — the signature  |

**Nonsense readthrough** — makes the ribosome read through a premature stop codon, so a full-length protein is produced from an allele that otherwise yields a truncated one.

- Only helps a nonsense allele. A missense, frameshift or deletion allele gets nothing from readthrough.
- Aminoglycoside readthrough is weak and the doses implied are nephrotoxic and ototoxic; this is why the class has not become a therapy despite thirty years of interest.

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

### Candidates that may beat the current standard

**The bar: strength 0.578, specificity 0.21, efficiency 0.737** (predicted) — the best established value on *each axis*, not one drug's scores. Representative incumbent: Alitretinoin (via TF).

2 of 12 candidates already have a reported link to this gene (CTD returned 184 chemicals with curated expression effects).

*Accounting for all 12 candidates: 2 already linked, 4 novel and non-dominated, 2 novel but beaten on every axis by another candidate, 4 with no advantage. All of them are in results.json under `improvements.all_candidates`.*

#### Compounds already reported to change this gene

These come from CTD's curated chemical–gene interactions: someone measured an expression change and a curator recorded it. That makes them the most concrete evidence in this report — but the direction is what decides whether each is a lead or a warning.

184 chemicals total: **71** move it the way you want, 81 move it the opposite way, 32 unclear.

**Expect the 'predicted by this run?' column to say no, and do not read that as a failed match.** CTD indexes MeSH, which is a vocabulary for the toxicology and pharmacology literature rather than a drug list, so what comes back is three different kinds of thing:

- **chemical** (155) — a compound. The only kind that can match a candidate.
- **exposure** (25) — formaldehyde, asbestos, bisphenol, tobacco smoke. Real measured effects on your gene, and not things anyone can be given.
- **category** (4) — MeSH headings such as "Anti-HIV Agents" or "Plant Extracts". Not compounds at all, so the column is blank rather than no.

Across all three directions, **2 of 184 matched a candidate from this run**. A low number here is the normal result: this report's drug tables are approved pharmacology, and most of CTD is not. It means the two sources are looking at different parts of chemical space, not that either is wrong.

**Reported to move BUB1B UP** — leads:

| Chemical | Kind | Predicted by this run? | CTD actions |
|---|---|---|---|
| Estradiol | chemical | no | affects^expression; increases^expression |
| Lidocaine | chemical | no | increases^expression |
| Methimazole | chemical | no | increases^expression |
| Phenobarbital | chemical | no | increases^expression |

**Reported to move BUB1B protein UP** (1 compounds of 15 with a protein report). The protein level can change through stability, translation or secretion with no change at the promoter, so these are not evidence that the promoter is involved. A compound that also appears above is reported at both levels.

| Chemical | Predicted by this run? | CTD actions |
|---|---|---|
| Paclitaxel ⚠halts mitosis | no | increases^expression |

⚠halts mitosis: the compound is suspected to hold cells in mitosis instead of letting them proceed, and BUB1B is a mitotic gene. A report that it raises BUB1B (mRNA or protein) can reflect cells held in mitosis, where such genes are at their highest, and not a switch at the promoter.

**10 compounds have a measured effect on this gene and no mechanism from this analysis.** Exposures and MeSH headings are excluded from that count, because the scan was never going to propose formaldehyde and listing it as a miss says nothing.

Most of what remains is research chemistry rather than medicine: tool compounds, probes and natural-product fractions that are measured because somebody had them on a shelf. Read them as evidence the gene IS movable, and check the few that are approved drugs individually.

Acetaminophen, Cisplatin, Estradiol, Ethinyl Estradiol, Folic Acid, Lidocaine, Methimazole, Phenobarbital, Propylthiouracil, Testosterone

**Reported to move it the WRONG way** — avoid, or use as negative controls:

Amphetamine, Azathioprine, Bicalutamide, Calcitriol, Cannabidiol, Copper Sulfate, Cyclosporine, Dasatinib, Doxorubicin, Etoposide, Fluorouracil, Fulvestrant, Ketamine, Leflunomide, Methamphetamine, Methotrexate, Oxaliplatin, Piroxicam, Pravastatin, Progesterone, Sirolimus, Sunitinib, Topotecan, Tretinoin, Troglitazone …

*CTD aggregates across cell types, doses and assays. A 'both' or contradictory annotation usually means the effect is context-dependent rather than that the data are wrong — check the underlying study before relying on any single row.*

#### Not previously linked, and better on at least one axis

| Candidate | Route | TF | Strength | Specificity | Efficiency | Wins on | Trade-offs |
|---|---|---|---|---|---|---|---|
| Chenodiol | via TF | NR1H4 | 0.503 | 0.337 | 0.765 | higher specificity, superior efficiency | stronger activation |
| Ursodiol | via TF | NR1H4 | 0.503 | 0.337 | 0.765 | higher specificity, superior efficiency | stronger activation |

*These are predicted scores compared against measured ones, so a margin was required before counting a win. They are hypotheses to test, not improvements demonstrated.*

Novel but beaten on every axis by another candidate (shown so nothing is hidden): Obeticholic Acid

#### Non-dominated set (nothing beats these on every axis)

- **Chenodiol** (NEW link) — trade-off: better on higher specificity, superior efficiency; worse on stronger activation
- **Ursodiol** (NEW link) — trade-off: better on higher specificity, superior efficiency; worse on stronger activation

### Strength of evidence for each factor route

Two links between a drug and your gene come from evidence and differ from row to row: whether the drug engages its target at a dose that can be reached, and whether the factor is likely to act on this promoter. Each is a 0 to 1 heuristic and not a calibrated probability. They are shown as they are, and they are not multiplied into a chance that the gene moves.

| Drug | TF | Drug engagement | Promoter activation | Weaker link | Knockout evidence |
|---|---|---|---|---|---|
| Acitretin ★ | RXRA | 0.704 | 0.793 | drug engagement | not among top changes, 0 of 2 |
| Bexarotene ★ | RXRA | 0.704 | 0.793 | drug engagement | not among top changes, 0 of 2 |
| Chenodiol ★ | NR1H4 | 0.736 | 0.509 | promoter activation | — |
| Ursodiol ★ | NR1H4 | 0.736 | 0.509 | promoter activation | — |
| Sulindac | RXRA | 0.704 | 0.793 | drug engagement | not among top changes, 0 of 2 |
| Etretinate | RXRA | 0.704 | 0.793 | drug engagement | not among top changes, 0 of 2 |
| Obeticholic Acid | NR1H4 | 0.736 | 0.509 | promoter activation | — |

★ = named in the comparison against the standard of care above, pinned here so the two sections cannot disagree.

**Knockout evidence** says whether your gene is among the genes that changed when the factor was knocked out, knocked down or overexpressed in published experiments. A response shows the factor can move the gene somewhere, not that it will in your cells or that a drug acting on the factor will. "Not among top changes" is weak, because each list holds only the few hundred genes that changed most.

Order candidates by tier and score, as the tables above do. These columns say where the evidence behind a route is thin.

## Virtual screening

A screen needs a model, and a model needs data. This one trains on measured binding at the receptor, so it runs only where ChEMBL holds at least 20 ligands at pChEMBL 6 or better, and only when the trained model separates known actives from decoys at CV AUC 0.7 or above.

**Screened:** NR1H4, RXRA.

**This is why the same receptor families appear here whatever gene you run.** The filter selects for how much medicinal chemistry a receptor has attracted, not for anything about your promoter. Nuclear receptors, and the retinoid and steroid receptors above all, have ligand series running to thousands of compounds; most transcription factors have none at all. So a receptor reaching this section means it is well studied, and the drugs listed under it are that receptor's pharmacology, not a result about your gene.

Read these rows as *candidate ligands for a factor that has a site in your promoter*. Whether moving that factor moves your gene is a separate question, and the activity and evidence sections above speak to it. A long list here is evidence about the receptor's chemistry; it is not evidence that the gene will respond.

### Virtual screen of approved drugs against RXRA

Binding model, CV AUC 0.967 — trained on None measured actives vs 102 decoys, screened 3311 approved drugs.

| Quadrant | Drug | What it is | Approved to treat | Measured pChEMBL | P(binds), model | Novelty |
|---|---|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma | 6.82 | 0.792 (in-sample) | — |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically | 6.48 | 0.375 (in-sample) | — |
| ESTABLISHED LIGAND, untested on this gene | Levothyroxine | synthetic T4 | hypothyroidism | 6.03 | 0.161 (in-sample) | — |
| no call | Lumacaftor | — | not established here | — | 0.441 | 0.703 |
| no call | Isotretinoin | retinoid | severe nodular acne | — | 0.375 | 0.794 |
| no call | Eluxadoline | — | not established here | — | 0.338 | 0.781 |
| no call | Seladelpar | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) | — | 0.322 | 0.738 |
| no call | Cobicistat | — | not established here | — | 0.321 | 0.766 |
| no call | Iptacopan | — | not established here | — | 0.321 | 0.741 |
| no call | Zaleplon | — | not established here | — | 0.317 | 0.743 |
| no call | Iptacopan Hydrochloride | — | not established here | — | 0.315 | 0.747 |
| no call | Seladelpar Lysine | PPAR agonist class | type 2 diabetes (thiazolidinediones) and hypertriglyceridaemia (fibrates) | — | 0.31 | 0.725 |
| no call | Telmisartan | — | not established here | — | 0.303 | 0.69 |
| no call | Methotrexate | — | not established here | — | 0.294 | 0.788 |
| no call | Pafolacianine | — | not established here | — | 0.294 | 0.854 |
| no call | Eprosartan | — | not established here | — | 0.291 | 0.736 |

Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts, novel chemistry and an untested gene link, so expect most to be wrong and treat each as needing a binding assay first.

**Read the indication column before the scores.** A receptor only reaches this section when it has a large ligand series, and such series are rarely homogeneous: a topical gel, an oral oncology agent and a withdrawn compound can sit on adjacent rows with near-identical numbers. The scores rank binding to the receptor and say nothing about whether the drug could be given for your purpose.

**The two score columns are not interchangeable, and each quadrant has only one real number.** An established ligand was a training molecule for this model, so its model score is an in-sample fit and its novelty is structurally zero, because the nearest known ligand the lookup finds is itself. Those rows are ranked by measured pChEMBL, which is data. Predicted rows have no measurement, so they are ranked by the model score, which is a prediction. Comparing a pChEMBL to a probability compares a measurement to a guess.

### Virtual screen of approved drugs against NR1H4

Binding model, CV AUC 0.999 — trained on 2176 measured actives vs 3307 decoys, screened 3311 approved drugs.

*Direction not annotated: only 9 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | What it is | Approved to treat | Measured pChEMBL | P(binds), model | Novelty |
|---|---|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Odevixibat | — | not established here | 9.8 | 0.846 (in-sample) | — |
| ESTABLISHED LIGAND, untested on this gene | Obeticholic Acid | — | not established here | 8.0 | 0.993 (in-sample) | — |
| ESTABLISHED LIGAND, untested on this gene | Simvastatin | HMG-CoA reductase inhibitor | hypercholesterolaemia | 6.24 | 0.369 (in-sample) | — |
| ESTABLISHED LIGAND, untested on this gene | Fulvestrant | sex steroid class | contraception, hormone replacement, and hormone-sensitive cancers | 6.1 | 0.397 (in-sample) | — |
| predicted binder, close to known chemistry | Odevixibat Sesquihydrate | — | not established here | — | 0.84 | 0.0 |

Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts, novel chemistry and an untested gene link, so expect most to be wrong and treat each as needing a binding assay first.

**Read the indication column before the scores.** A receptor only reaches this section when it has a large ligand series, and such series are rarely homogeneous: a topical gel, an oral oncology agent and a withdrawn compound can sit on adjacent rows with near-identical numbers. The scores rank binding to the receptor and say nothing about whether the drug could be given for your purpose.

**The two score columns are not interchangeable, and each quadrant has only one real number.** An established ligand was a training molecule for this model, so its model score is an in-sample fit and its novelty is structurally zero, because the nearest known ligand the lookup finds is itself. Those rows are ranked by measured pChEMBL, which is data. Predicted rows have no measurement, so they are ranked by the model score, which is a prediction. Comparing a pChEMBL to a probability compares a measurement to a guess.

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

---

# Testing any of this in the lab

Everything above is prediction. This part is experiment — how you would find out whether any of it holds, and it is placed last because it is long, not because it is optional.

**If you run one thing, run the reporter assay with the scrambled-motif control.** A drug that induces the wild-type promoter but not the mutant has an effect running through the element this report predicted. A drug that induces both means the prediction was wrong. That single control turns a hit into a mechanism.

## Designing a direct modality for BUB1B

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

The outcome chain above is almost always limited by *gene responsiveness in your cell type* — the link this package cannot see. A screen measures it instead of assuming it. Everything above prioritises this experiment; it does not replace it.

### Primary screen: BUB1B protein per cell, in the patient's own cells

The endpoint is the protein and not a reporter, because the only working copy is itself impaired, and the transcript carries a premature stop. A promoter reporter measures transcription from a cloned fragment, so it cannot see a stabiliser, a splice modifier or a readthrough agent - and those are the mechanisms in play here. This screen is also mechanism-agnostic, which means it does not depend on any ranking in this report being right.

**Cells**

- Patient-derived primary fibroblasts (skin biopsy) or EBV transformed lymphoblastoid lines. Primary fibroblasts are closer to the biology; LCLs scale better to 384-well.
- On the same plate: a carrier parent if available, and at least two unrelated healthy donors. The carrier is the informative control, because the target level is what a healthy heterozygote already demonstrates is tolerated.
- Keep every well sub-confluent and proliferating. This gene's protein tracks the cell cycle, so a density difference between plates is indistinguishable from a drug effect.

**Readout**

- Primary: high-content immunofluorescence, per-cell protein intensity. Imaging is worth the extra effort here because it gives per-cell data, and per-cell data is what lets you gate on cell-cycle phase.
- Alternatives at 384-well: in-cell western (two-colour, normalised to a total-protein stain) or HTRF/AlphaLISA. Both give a well average, so neither can separate 'more protein per cell' from 'more cells in G2/M'.
- Confirmation: capillary or conventional western on the confirmed hits, which is the only readout that shows you the band is the right size - relevant when one allele may make a truncated product.
- Validate the antibody against a knockdown before screening. An antibody that does not drop with siRNA is measuring something else, and the whole screen then measures it too.

**Normalisation, and the control that matters most**

- Total protein stain, or a housekeeping protein that is NOT cell-cycle regulated. Not PCNA, not a cyclin, not Ki-67.
- Co-stain DNA content or phospho-histone H3 in the same well. This is the control that matters most: it turns 'the signal went up' into 'the signal went up at matched cell-cycle phase', and without it the screen's strongest hits will be compounds that arrest cells in mitosis.

**Controls**

- Positive: a proteasome inhibitor (MG132, or bortezomib). It raises almost any protein, which makes it a poor therapy and an excellent control - it proves the assay can detect an increase at all, and it sets the ceiling.
- Negative: vehicle (DMSO) at the screening concentration, many wells, spread across the plate rather than in one column.
- Cytotoxicity in parallel (CellTiter-Glo or equivalent) on a sister plate. A dying well reads as whatever the normalisation does under stress, and most primary hits in screens like this are simply toxic.
- Cell-cycle counter-screen on the primary hits: anything that arrests cells in G2/M raises this protein for a reason that is not a treatment.

**Thresholds**

- Z' >= 0.5 between the proteasome inhibitor and vehicle before screening anything. Below that the plate cannot tell a hit from noise.
- Primary hit: >= 1.5x vehicle, after excluding wells with more than 30% viability loss or a mitotic index more than twice vehicle.
- Confirm on an 8-point dose response in triplicate. A real stabiliser gives a dose-response; a plate artefact does not.

**Triage, in order**

- 1. Dose-response confirms, cytotoxicity excluded.
- 2. Western confirms the increase and the band size.
- 3. The increase survives cell-cycle matching.
- 4. RT-qPCR says whether the transcript moved too. Protein up with transcript flat points at stability or translation; both up points at transcription, and that route is testable against the rest of this report.
- 5. Read the factor the route runs through, in the same wells, and then test whether the effect survives losing it. See 'Testing the route' below: this is the step that separates a compound that works for the reason this report gives from one that works for another reason, and it is also what makes a negative interpretable.
- 6. Only then is it a candidate, and it is a candidate for a cell assay of function, not for a patient.

**Testing the route**

Everything above establishes that a hit is real. This establishes whether it works for the reason this report gives, which is a different question and the one the ranking rests on. Run it on confirmed hits only, not in the primary screen.

*What to measure*

- The factor's own level, in the same wells as the BUB1B readout (RXRA, OSR2, THAP11, KLF3). A drug can engage a factor without changing how much of it there is, so for a ligand read activity rather than abundance: a target-gene panel the factor is already known to control, or its occupancy at the predicted site by ChIP-qPCR.
- BUB1B protein and BUB1B transcript, which step 4 has already produced for this compound.

*What each result means*

| Result | What to do |
|---|---|
| BUB1B up, factor engaged | the route is supported, and the next compound to try is the one with the best potency at that factor |
| BUB1B up, factor NOT engaged | a real hit by another mechanism. Do not optimise on the factor, which would be tuning the wrong variable; follow the transcript reading instead, because protein up with transcript flat points at stabilisation and almost nothing in this report covers that axis |
| BUB1B flat, factor engaged | the factor-to-gene link is false in this cell type. This is the most economical result in the whole screen: it retires every compound that reaches the gene through that factor at once, not just the one tested |
| BUB1B flat, factor NOT engaged | nothing was tested. The compound did not reach its target at this dose or time, so this says nothing about the route, and discarding the branch on it would be an error |

*The decisive test*

- Knock the factor down (siRNA or a degron line), then add the compound. If BUB1B still rises without the factor, the effect does not run through it and the route is wrong whatever the correlations say.
- This is the practical form of the site-directed promoter mutant named elsewhere in this report. The mutant is the stronger evidence, because it shows the effect runs through one element rather than through the factor by any route; the knockdown is faster, cheaper and answers for the whole family of compounds at once.
- Include the knockdown alone as its own arm. A factor whose loss moves the gene by itself tells you the direction to push it, which is the question the computational part of this report could not settle.

*A factor that partners promiscuously - a nuclear receptor that heterodimerises with many others, for one - makes its abundance a poor readout, because the compound may reach the gene through a partner instead. For those, occupancy at the predicted element or the promoter mutant is the readout that discriminates, and a level alone is not.*

**Library**

- an FDA-approved collection, so anything that hits is already through safety
  - Broad Drug Repurposing Hub, approved subset (~2,400 compounds)
  - Prestwick Chemical Library (~1,520, all off-patent and approved)
  - Selleck FDA-approved Drug Library (~2,900)
  - NCATS Pharmaceutical Collection (NPC, ~2,800)
- Format: single concentration primary, typically 10 uM, in duplicate; two exposures (24 h and 48 h), because a stabiliser and a transcriptional effect separate in time

*A screen in fibroblasts or LCLs is not the affected tissue. A compound that raises the protein there may not reach or work in the tissue that matters, and the exposure check in `exposure.py` is the next question for anything that survives.*

### CRISPRa screen

- Approach: dCas9-VP64 with SAM (MS2-p65-HSF1) or dCas9-VPR for stronger activation
- Focused library: 17 TFs prioritised by this run (RXRA, OSR2, THAP11, KLF3, SP3, KLF16, ZNF76, ZNF143, FLI1, TFAP4...)
- Unbiased: ~1,600 human TFs (Calabrese or a TF-focused sublibrary). Run this if the focused set comes back empty — an unbiased screen is the only thing that finds a factor nobody has associated with your promoter.
- Coverage: 3-5 sgRNAs per gene, >=500 cells per guide at transduction, MOI ~0.3 to keep single integrations
- Readout: FACS sort GFP-high vs GFP-low, then sequence guide abundance; or arrayed luciferase for a focused set

Positive-control sgRNAs at your promoter (12 clean of 12 shown):

| Protospacer (20 nt) | PAM | Strand | vs TSS | GC | Note |
|---|---|---|---|---|---|
| `ACTCCGGACGGGTGAGATTT` | GG | + | -396 | 55% | ok |
| `GAGATTTGGGGCAGACAGCA` | GG | + | -383 | 55% | ok |
| `AGATTTGGGGCAGACAGCAG` | GG | + | -382 | 55% | ok |
| `TAAGCCTGCTGCACTTCCAC` | GG | + | -198 | 55% | ok |
| `ATTATGCACCGCGAGTTCGC` | GG | + | -148 | 55% | ok |
| `GCGAGTTCGCGGACTAAACA` | GG | + | -138 | 55% | ok |
| `ACAAGCCCCAGAATGCCTTG` | GG | + | -108 | 55% | ok |
| `TGTTTAGTCCGCGAACTCGC` | GG | - | -138 | 55% | ok |
| `GCGAACTCGCGGTGCATAAT` | GG | - | -148 | 55% | ok |
| `CGAACTCGCGGTGCATAATG` | GG | - | -149 | 55% | ok |

*Guides listed target YOUR promoter directly, for a positive-control arm proving the reporter can be activated at all. The library guides target the TFs. Run these through CRISPOR or CHOPCHOP for off-target scoring before ordering — no off-target analysis is done here.*

### Compound screen

- Prestwick Chemical Library (~1,520 approved drugs, high chemical and pharmacological diversity)
- Selleck FDA-approved Drug Library (~3,000)
- NCATS Pharmaceutical Collection (NPC) — free, and the most complete approved set

- Design: Arrayed, 384-well, 8-point dose-response on the actives. Single-concentration primary screens at 10 uM produce a cytotoxicity-dominated hit list; dose-response separates real induction from cells dying.
- Counter-screen: Run the same plate against a constitutive-promoter reporter. Anything that moves both is affecting transcription or translation globally, not your promoter, and that is the majority of primary hits.

#### Cherry-pick list

**What this is.** "Cherry-picking" means having the vendor pull these named compounds into a custom plate rather than buying a whole library. It is a cost and time decision: ~48 wells you can run in a week, against ~1,520 compounds for Prestwick or ~3,000 for the Selleck set.

**Where these came from.** Every compound below was proposed earlier in *this run* — nothing external was consulted to build the list:

- 0 direct gene-targeted agent(s) — priority 1, because they act on the gene or its product and skip every inferential step
- 12 graded TF/drug route(s) — priority 2 if graded A/B, else 3

Priority means:

- **1** — acts directly on your gene, transcript or protein
- **2** — TF route the grading believes (A or B)
- **3** — TF route with a weaker grade — included to fill the plate, and the first rows to cut if wells are scarce

**What it cannot do.** Every compound here was proposed by this run. A cherry-pick therefore tests the model and nothing beyond it — it cannot discover a mechanism the pipeline never predicted. If the whole plate is negative, the honest reading is that the promoter/TF model is wrong for this gene, not that no drug moves it.

Showing 7 of 7 selected (from 12 candidates in this run), in priority order:

| Priority | Compound | Why it is on the plate | Pulled from |
|---|---|---|---|
| 2 | Acitretin | RXRA route, grade B | graded TF route via RXRA |
| 2 | Bexarotene | RXRA route, grade B | graded TF route via RXRA |
| 2 | Sulindac | RXRA route, grade B | graded TF route via RXRA |
| 2 | Chenodiol | NR1H4 route, grade B | graded TF route via NR1H4 |
| 2 | Ursodiol | NR1H4 route, grade B | graded TF route via NR1H4 |
| 3 | Etretinate | RXRA route, grade C | graded TF route via RXRA |
| 3 | Obeticholic Acid | NR1H4 route, grade D | graded TF route via NR1H4 |

*A full-library screen is unbiased and will find mechanisms this package cannot predict — that is its advantage, not a criticism of it. Use the cherry-pick list to get an answer in a week; run the full library when you want the answer you did not expect.*

### What to expect

- **compound screen hit rate**: 0.5-2% primary at a sensible threshold, and most of those are cytotoxic or cell-cycle effects. Expect to lose 60-80% at counter-screen and to finish with single digits of confirmed compounds, possibly none.
- **replicates**: duplicate primary, triplicate for confirmation. A single primary well cannot be distinguished from a pipetting error.
- **quality metric**: Z' >= 0.5 against the proteasome-inhibitor control, per plate, recorded for every plate.
- **what a negative result means**: No approved compound raises this protein in this cell type at this exposure. That is an informative result and a cheap one: it says repurposing is unlikely to solve it, and redirects effort to chemistry aimed at the specific mechanism - the degradation route or the splice defect - rather than to more screening.
- **why this screen and not the reporter**: A promoter reporter measures transcription from a cloned fragment. It cannot see a stabiliser, a splice modifier or a readthrough agent, which are the mechanisms in play when the protein is the problem.
