# Downregulation strategy for CCR5 — promoter / TF / drug analysis

Objective: **decrease** transcription of CCR5.

## Drugs targeting CCR5 directly

These act on the gene, transcript or protein itself, skipping the transcription-factor route entirely. Where one of these exists it is almost always the better starting point: sequence-directed agents affect your gene and essentially nothing else.

| Drug | Modality | Specificity | Raises protein? | Mechanism | Approved | Source |
|---|---|---|---|---|---|---|
| Maraviroc | small molecule | 0.55 | activity only | C-C chemokine receptor type 5 antagonist | 2007 | ChEMBL |

*"activity only" means the drug raises the function of protein already present without increasing its abundance. If your readout is a western blot rather than an activity assay, that distinction decides whether the experiment appears to work.*

## Promoter of CCR5 (ENST00000292303)

- Transcript chosen: **distinct TSS (representing 2 within 500 bp)**
- Assembly: GRCh38  |  Strand: +
- TSS: chr3:46370946
- Window: chr3:46365946-46371945 (-5000 / +1000, 6000 bp)
- Motifs scanned: 879 at p ≤ 0.0001 per position
- Non-redundant sites: 900 across 410 TFs

> **⚠ Weak scan.** 900 sites against 744 expected by chance over 4230 scannable bases. Most of this list is noise; treat only the highest-scoring sites with ChIP support as real.

### Ligand-activated (nuclear) receptors with sites here

These are the directly druggable ones.

| TF | Role | Sites | Best rel. score | Best site (vs TSS) | ChIP peaks | Score |
|---|---|---|---|---|---|---|
| ESRRG | activator | 3 | 0.987 | CCAAGGTCA (-999) | 0 | 65.8 |
| NR1I3 | activator | 2 | 1.0 | TGAACTTT (-2691) | 0 | 61.7 |
| PPARD | repressor | 1 | 0.839 | AGGGTCACAGTTTG (+184) | 0 | 61.5 |
| ESRRB | activator | 2 | 0.95 | CCAAGGTCAC (-999) | 0 | 61.4 |
| PGR | activator | 2 | 0.971 | ACATCATGT (-3492) | 0 | 60.8 |
| RARA | repressor | 6 | 0.839 | GTGATCAAAAGTTCATTT (-2698) | 0 | 58.7 |
| NR2C2 | activator | 1 | 1.0 | GGGTCA (+185) | 0 | 57.8 |
| PPARG | dual | 6 | 0.872 | AGAGCCAAGGTCA (-1003) | 0 | 56.4 |
| NR1H3 | activator | 4 | 0.873 | TGAACTTAAATAGACT (+684) | 0 | 53.9 |
| ESR1 | activator | 1 | 0.778 | GGGCTTCCGTGACCT (-996) | 0 | 52.9 |
| NR2F1 | repressor | 4 | 0.923 | CAAGGTCACGGA (-998) | 0 | 52.4 |
| NR2F2 | repressor | 2 | 1.0 | AAGGTCA (+890) | 0 | 51.8 |
| RXRB | activator | 1 | 0.788 | AGGGTCACAGTTTG (+184) | 0 | 51.5 |
| PPARA | dual | 4 | 1.0 | AAGGTCA (-997) | 0 | 50.2 |
| NR1D2 | repressor | 3 | 0.839 | TGGGTTAGTCTGCCT (-1100) | 0 | 49.6 |
| RXRA | unknown | 12 | 0.819 | GGGACACAGGGTTAA (-505) | 0 | 48.8 |
| NR2F6 | repressor | 4 | 0.837 | GTGATCAAAAGTTCA (-2698) | 0 | 48.6 |
| NR5A1 | activator | 2 | 0.854 | GAGCCAAGGTCA (-1002) | 0 | 48.5 |
| NR4A2 | activator | 1 | 0.948 | GAAGGTCA (+890) | 0 | 48.4 |
| NR4A1 | activator | 1 | 0.928 | GAAGGTCA (+890) | 0 | 47.8 |
| NR1I2 | activator | 1 | 0.873 | CAAACTCCCTGCACC (+625) | 0 | 47.7 |
| NR6A1 | dual | 3 | 0.887 | TATGGTCAAGTTCA (-848) | 0 | 47.4 |
| THRB | unknown | 7 | 0.854 | TGGGCTTCCGTGACCTTG (-998) | 0 | 46.8 |
| VDR | unknown | 2 | 1.0 | GAGTTCA (-1881) | 0 | 46.8 |
| RARG | unknown | 5 | 0.851 | AAGGTCACGGAAGCCCA (-997) | 0 | 46.8 |

### Top-ranked transcription factors

| Rank | TF | Role | NR | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|
| 1 | NHLH1 | activator |  | 7 | 0.928 | +55 | 0 | 0 | 65.9 |
| 2 | ESRRG | activator | yes | 3 | 0.987 | -999 | 0 | 1 | 65.8 |
| 3 | GATA2 | activator |  | 3 | 1.0 | +174 | 0 | 0 | 64.7 |
| 4 | ETV2 | activator |  | 4 | 0.922 | +104 | 0 | 0 | 64.2 |
| 5 | MAFG | activator |  | 2 | 0.869 | +6 | 0 | 0 | 62.5 |
| 6 | OSR2 | activator |  | 2 | 0.949 | +79 | 0 | 0 | 62.4 |
| 7 | NR1I3 | activator | yes | 2 | 1.0 | -2691 | 0 | 2 | 61.7 |
| 8 | PPARD | repressor | yes | 1 | 0.839 | +184 | 0 | 5 | 61.5 |
| 9 | PBX1 | activator |  | 1 | 0.971 | -28 | 0 | 0 | 61.5 |
| 10 | ESRRB | activator | yes | 2 | 0.95 | -999 | 0 | 1 | 61.4 |
| 11 | GATA4 | activator |  | 2 | 1.0 | +174 | 0 | 0 | 61.4 |
| 12 | ZNF784 | activator |  | 1 | 0.972 | +42 | 0 | 0 | 61.0 |
| 13 | MAFK | repressor |  | 2 | 0.862 | +6 | 0 | 0 | 60.8 |
| 14 | PGR | activator | yes | 2 | 0.971 | -3492 | 0 | 1 | 60.8 |
| 15 | MEIS1 | activator |  | 1 | 0.943 | -27 | 0 | 0 | 60.7 |
| 16 | DUX4 | activator |  | 1 | 0.941 | -27 | 0 | 0 | 60.6 |
| 17 | ONECUT1 | activator |  | 1 | 0.939 | -26 | 0 | 0 | 60.6 |
| 18 | BARX2 | activator |  | 1 | 0.983 | -89 | 0 | 0 | 59.8 |
| 19 | RARA | repressor | yes | 6 | 0.839 | -2698 | 0 | 9 | 58.7 |
| 20 | MAF | activator |  | 1 | 0.851 | +6 | 0 | 0 | 58.6 |
| 21 | VAX2 | repressor |  | 5 | 1.0 | -440 | 0 | 0 | 58.5 |
| 22 | NR2C2 | activator | yes | 1 | 1.0 | +185 | 0 | 0 | 57.8 |
| 23 | NKX3-2 | repressor |  | 1 | 0.949 | -96 | 0 | 0 | 57.1 |
| 24 | SMAD2 | activator |  | 2 | 1.0 | +430 | 0 | 0 | 56.8 |
| 25 | BHLHA15 | activator |  | 3 | 0.969 | +707 | 0 | 0 | 56.6 |

### Every route that acts on this gene

216 routes across 4 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

- 1 × clinical — approved for this gene, documented in humans
- 201 × curated — mechanism annotated by a human curator
- 13 × predicted — a model says so
- 1 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 1 | direct: gene product | Maraviroc | CCR5 | 0 | C-C chemokine receptor type 5 antagonist | 0.55 | 0.55 | 0.75 | 60.0 |
| 3 | TF: orthosteric ligand | Calcitriol | VDR | 1 | ligand at the transcription factor | 0.856 | 1.0 | 0.99 | 75.2 |
| 3 | TF: orthosteric ligand | Seladelpar | PPARD | 1 | ligand at the transcription factor | 0.87 | 1.0 | 0.933 | 74.5 |
| 3 | TF: orthosteric ligand | Calcipotriene | VDR | 1 | ligand at the transcription factor | 0.867 | 1.0 | 0.925 | 74.2 |
| 3 | TF: orthosteric ligand | Fulvestrant | ESR1 | 1 | ligand at the transcription factor | 0.845 | 0.926 | 0.967 | 72.3 |
| 3 | TF: orthosteric ligand | Trifarotene | RARG | 1 | ligand at the transcription factor | 0.8 | 1.0 | 0.865 | 70.9 |
| 3 | TF: orthosteric ligand | Cholecalciferol | VDR | 1 | ligand at the transcription factor | 0.796 | 1.0 | 0.862 | 70.7 |
| 3 | TF: orthosteric ligand | Mifepristone | PGR | 1 | ligand at the transcription factor | 0.902 | 0.741 | 1.0 | 69.6 |
| 3 | TF: orthosteric ligand | Belzutifan | HLF | 1 | ligand at the transcription factor | 0.747 | 1.0 | 0.863 | 69.2 |
| 3 | TF: orthosteric ligand | Seladelpar Lysine | PPARD | 1 | ligand at the transcription factor | 0.679 | 1.0 | 0.725 | 64.2 |
| 3 | TF: orthosteric ligand | Elacestrant | ESR1 | 1 | ligand at the transcription factor | 0.678 | 0.926 | 0.794 | 63.5 |
| 3 | TF: orthosteric ligand | Danthron | RXRA | 1 | ligand at the transcription factor | 0.647 | 1.0 | 0.725 | 63.2 |
| 3 | TF: orthosteric ligand | Dextrothyroxine | THRB | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Dextrothyroxine Sodium | THRB | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Liothyronine | THRB | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Alfacalcidol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Calcifediol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Dihydrotachysterol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Ergocalciferol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Paricalcitol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Levothyroxine Sodium | THRB | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Liothyronine Sodium | THRB | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Palovarotene | RARG | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Dextrothyroxine | THRA | 1 | ligand at the transcription factor | 0.621 | 1.0 | 0.725 | 62.4 |
| 3 | TF: orthosteric ligand | Dextrothyroxine Sodium | THRA | 1 | ligand at the transcription factor | 0.621 | 1.0 | 0.725 | 62.4 |
| 3 | TF: orthosteric ligand | Liothyronine | THRA | 1 | ligand at the transcription factor | 0.621 | 1.0 | 0.725 | 62.4 |
| 3 | TF: orthosteric ligand | Levothyroxine Sodium | THRA | 1 | ligand at the transcription factor | 0.621 | 1.0 | 0.725 | 62.4 |
| 3 | TF: orthosteric ligand | Liothyronine Sodium | THRA | 1 | ligand at the transcription factor | 0.621 | 1.0 | 0.725 | 62.4 |
| 3 | TF: orthosteric ligand | Resmetirom | THRB | 1 | ligand at the transcription factor | 0.606 | 1.0 | 0.731 | 62.0 |
| 3 | TF: orthosteric ligand | Elacestrant Hydrochloride | ESR1 | 1 | ligand at the transcription factor | 0.642 | 0.926 | 0.725 | 61.0 |

Best in each category:

- **direct: gene product** — Maraviroc (tier 1, score 60.0)
- **TF: orthosteric ligand** — Calcitriol (tier 3, score 75.2)
- **upstream signalling** — Ibrutinib (tier 4, score 21.3)
- **virtual screen hypothesis** — Palovarotene (tier 5, score 51.5)

*No route found in: direct: gene or transcript, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, TF: degradation blocker, empirical signature, chromatin-level. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

**direct: gene product** (1 found)

- Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- *Why:* Immediate, well-understood pharmacology; the protein is right there rather than three inferential steps away.
- *Cost:* Usually changes the protein's *activity*, not its abundance. If your readout is a western blot, it can look like nothing happened.
- *Differs:* Does not touch expression at all — it is here because it changes functional output, which is often what you actually want.

**TF: orthosteric ligand** (201 found)

- Binds the transcription factor's natural ligand pocket.
- *Why:* The best-characterised route, and where nearly all approved TF pharmacology lives.
- *Cost:* Only ~48 human TFs (the nuclear receptors) have such a pocket, and the factor's whole regulon moves with your gene.
- *Differs:* Clamps the receptor on or off regardless of endogenous signalling — contrast the allosteric and ligand-axis routes, which scale what is already happening.

**upstream signalling** (13 found)

- Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- *Why:* Turns 'this factor has no drugs' into a real list. Most TFs are only reachable this way.
- *Cost:* Specificity falls with every step: a node three hops up reaches your gene through one branch and thousands of other genes through the rest.
- *Differs:* The only category where the drug never touches your gene or its TF — read the printed path and judge it yourself.

**virtual screen hypothesis** (1 found)

- An approved drug predicted to bind one of your TFs, never tested against it.
- *Why:* Real repurposing space — most approved drugs have been assayed against a tiny fraction of the proteome.
- *Cost:* An extrapolation beyond the training chemistry. Expect most to be wrong; each is a binding assay, not a finding.
- *Differs:* The only category with no prior evidence of the drug–target pair existing at all.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **TF: allosteric modulator** — Binds a site away from the ligand pocket, amplifying or damping the endogenous ligand's effect (PAM/NAM).
- **TF: stabiliser or glue** — Locks the factor in an active conformation, or glues it to a partner protein.
- **TF: endogenous ligand** — Changes how much natural ligand exists, by drugging the enzyme that makes or destroys it.
- **TF: degradation blocker** — Stops the factor being destroyed, so more of it accumulates.
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **chromatin-level** — HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.

### Graded candidates

Combined = geometric mean of TF confidence and drug quality. A weak TF prediction caps the grade no matter how good the drug.

| Grade | Combined | Drug | TF | TF conf. | Drug qual. | Potency | Protein select. | Gene select. | Evidence | Flags |
|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 70.1 | Seladelpar | PPARD | 0.615 | 0.953 | 0.925 | 1.0 | 1.0 | 0.75 | direction not independently predicted |
| **B** | 66.1 | Seladelpar Lysine | PPARD | 0.615 | 0.785 | 0.5 | 1.0 | 1.0 | 0.5 | no potency data; scored neutral; thin evidence base (n=0 potency measu |
| **B** | 65.6 | Diethylstilbestrol | ESRRG | 0.658 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 64.9 | Elafibranor | PPARD | 0.615 | 0.739 | 0.4 | 0.741 | 1.0 | 0.75 | direction not independently predicted |
| **B** | 64.7 | Tamibarotene | RARA | 0.587 | 0.811 | 0.581 | 0.741 | 1.0 | 0.75 | direction not independently predicted |
| **B** | 63.6 | Elacestrant | ESR1 | 0.529 | 0.832 | 0.58 | 1.0 | 0.926 | 0.75 | structure supports antagonist (p=0.273, moderate confidence) |
| **B** | 63.4 | Fulvestrant | ESR1 | 0.529 | 0.967 | 0.917 | 1.0 | 0.926 | 1.0 | structure supports antagonist (p=0.458, low confidence) |
| **B** | 62.8 | Bezafibrate | PPARD | 0.615 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.8 | Tretinoin | PPARD | 0.615 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.6 | Clotrimazole | NR1I3 | 0.617 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.6 | Meclizine | NR1I3 | 0.617 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.4 | Diethylstilbestrol | ESRRB | 0.614 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.2 | Elacestrant Hydrochloride | ESR1 | 0.529 | 0.77 | 0.5 | 1.0 | 0.926 | 0.5 | no potency data; scored neutral; thin evidence base (n=0 potency measu |
| **B** | 61.2 | Mifepristone | PGR | 0.608 | 0.621 | 1.0 | 0.741 | 1.0 | 1.0 | black-box warning; direction not independently predicted |
| **B** | 60.8 | Acitretin | RARA | 0.587 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 60.8 | Etretinate | RARA | 0.587 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 60.8 | Isotretinoin | RARA | 0.587 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 60.8 | Trifarotene | RARA | 0.587 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.7 | Bexarotene | PPARG | 0.564 | 0.644 | 0.5 | 0.5 | 0.949 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.7 | Diclofenac Sodium | PPARG | 0.564 | 0.644 | 0.5 | 0.5 | 0.949 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.7 | Resveratrol | PPARG | 0.564 | 0.644 | 0.5 | 0.5 | 0.949 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.7 | Calcitriol | VDR | 0.468 | 0.994 | 0.975 | 1.0 | 1.0 | 1.0 | direction not independently predicted |
| **B** | 58.3 | Adapalene | RARA | 0.587 | 0.575 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.3 | Alitretinoin | RARA | 0.587 | 0.575 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.3 | Tazarotene | RARA | 0.587 | 0.575 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |

### Cell-cycle compatibility

Tissue proliferation state: **unknown**

3 factor(s) here are cell-cycle restricted (MYC, E2F4, E2F3). Set --proliferation to have this scored rather than just flagged.

| TF | Phase | Duty cycle | Effect multiplier |
|---|---|---|---|
| MYC | G1 | 46% | 1.0 |
| E2F4 | G1 | 46% | 1.0 |
| E2F3 | G1/S | 25% | 1.0 |

*A partial duty cycle dilutes the population-level effect in asynchronous cells: per-cell efficacy can be much higher than a bulk RT-qPCR suggests. Synchronise or use a phase-resolved readout before calling it a failure.*

### Upstream routes to factors with no direct drug

Most TFs have no approved ligand. These paths walk upstream through signed causal edges — kinases, phosphatases, receptors — to proteins that *are* druggable. Signs multiply along the path, so an inhibitor of an inhibitor counts as activating.

Distance buys druggability and spends selectivity: a node three hops up reaches your gene through one branch and everything else through the others. Read the fan-out column.

**MAFG** (no direct pharmacology)

| Upstream target | Hops | Fan-out | Path | Act on it | Drugs | Score |
|---|---|---|---|---|---|---|
| BTK | 3 | 19 | BTK --> GTF2I --> PRRX1 —| MAFG | inhibit | Ibrutinib, Tirabrutinib, Acalabrutinib | 20.0 |
| JAK2 | 3 | 58 | JAK2 --> GTF2I --> PRRX1 —| MAFG | inhibit | Ruxolitinib Phosphate, Momelotinib, Baricitinib | 11.4 |
| SRC | 3 | 280 | SRC --> GTF2I --> PRRX1 —| MAFG | inhibit | Bosutinib, Vandetanib, Tirbanibulin | 3.4 |

- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 19 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 280 downstream branches
- ⚠ 3 hops upstream: druggable, but it reaches your gene through one of 58 downstream branches

### Mechanism options: orthosteric vs allosteric

Allosteric sites are not constrained to bind the endogenous ligand, so they diverge across a receptor family where the orthosteric pocket cannot. That usually buys subtype selectivity — at the cost of a bounded maximal effect and a hard dependence on endogenous tone.

**ESRRG**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| inhibitor | 1 | unknown | unknown | no | no | Diethylstilbestrol |

*No approved allosteric modulator at this target; only orthosteric pharmacology, with the usual subtype-selectivity problem.*

**Before pursuing any allosteric option, verify that the pathway has endogenous tone in your cell type.** A PAM with no ligand to modulate does nothing at all, and that is a measurement, not a prediction — nothing in this package can check it for you.

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

**NHLH1** — **no ligand pocket** · stabilisation recommended

- NHLH1 has no ligand-binding pocket and no direct drug. This is exactly the case stabilisation exists for: its protein-protein interfaces and its degradation machinery are both druggable even though its surface is not.
- Nothing curated for this factor specifically. The tractable questions are: which E3 ligase degrades it (block that to raise abundance), and which partner complex sequesters or activates it (glue or disrupt that). Both are established, approved mechanisms — see the exemplars.

  Approved proof the mechanism class works, on other targets:
  - **Tafamidis** (2019) — conformational stabiliser: Binds the transthyretin tetramer interface and prevents dissociation — the cleanest approved example of locking a protein into a stable shape
  - **Migalastat** (2016) — conformational stabiliser: Pharmacological chaperone; stabilises mutant alpha-galactosidase A into a foldable, traffickable form
  - **Sirolimus** (1999) — molecular glue: Glues FKBP12 to mTOR; neither interaction occurs without the drug
  - **Tacrolimus** (1994) — molecular glue: Glues FKBP12 to calcineurin

**GATA2** — **no ligand pocket** · stabilisation recommended

- GATA2 has no ligand-binding pocket and no direct drug. This is exactly the case stabilisation exists for: its protein-protein interfaces and its degradation machinery are both druggable even though its surface is not.
- Nothing curated for this factor specifically. The tractable questions are: which E3 ligase degrades it (block that to raise abundance), and which partner complex sequesters or activates it (glue or disrupt that). Both are established, approved mechanisms — see the exemplars.

  Approved proof the mechanism class works, on other targets:
  - **Tafamidis** (2019) — conformational stabiliser: Binds the transthyretin tetramer interface and prevents dissociation — the cleanest approved example of locking a protein into a stable shape
  - **Migalastat** (2016) — conformational stabiliser: Pharmacological chaperone; stabilises mutant alpha-galactosidase A into a foldable, traffickable form
  - **Sirolimus** (1999) — molecular glue: Glues FKBP12 to mTOR; neither interaction occurs without the drug
  - **Tacrolimus** (1994) — molecular glue: Glues FKBP12 to calcineurin

**ETV2** — **no ligand pocket** · stabilisation recommended

- ETV2 has no ligand-binding pocket and no direct drug. This is exactly the case stabilisation exists for: its protein-protein interfaces and its degradation machinery are both druggable even though its surface is not.
- Nothing curated for this factor specifically. The tractable questions are: which E3 ligase degrades it (block that to raise abundance), and which partner complex sequesters or activates it (glue or disrupt that). Both are established, approved mechanisms — see the exemplars.

  Approved proof the mechanism class works, on other targets:
  - **Tafamidis** (2019) — conformational stabiliser: Binds the transthyretin tetramer interface and prevents dissociation — the cleanest approved example of locking a protein into a stable shape
  - **Migalastat** (2016) — conformational stabiliser: Pharmacological chaperone; stabilises mutant alpha-galactosidase A into a foldable, traffickable form
  - **Sirolimus** (1999) — molecular glue: Glues FKBP12 to mTOR; neither interaction occurs without the drug
  - **Tacrolimus** (1994) — molecular glue: Glues FKBP12 to calcineurin

**MAFG** — **no ligand pocket** · stabilisation recommended

- MAFG has no ligand-binding pocket and no direct drug. This is exactly the case stabilisation exists for: its protein-protein interfaces and its degradation machinery are both druggable even though its surface is not.
- Nothing curated for this factor specifically. The tractable questions are: which E3 ligase degrades it (block that to raise abundance), and which partner complex sequesters or activates it (glue or disrupt that). Both are established, approved mechanisms — see the exemplars.

  Approved proof the mechanism class works, on other targets:
  - **Tafamidis** (2019) — conformational stabiliser: Binds the transthyretin tetramer interface and prevents dissociation — the cleanest approved example of locking a protein into a stable shape
  - **Migalastat** (2016) — conformational stabiliser: Pharmacological chaperone; stabilises mutant alpha-galactosidase A into a foldable, traffickable form
  - **Sirolimus** (1999) — molecular glue: Glues FKBP12 to mTOR; neither interaction occurs without the drug
  - **Tacrolimus** (1994) — molecular glue: Glues FKBP12 to calcineurin

**OSR2** — **no ligand pocket** · stabilisation recommended

- OSR2 has no ligand-binding pocket and no direct drug. This is exactly the case stabilisation exists for: its protein-protein interfaces and its degradation machinery are both druggable even though its surface is not.
- Nothing curated for this factor specifically. The tractable questions are: which E3 ligase degrades it (block that to raise abundance), and which partner complex sequesters or activates it (glue or disrupt that). Both are established, approved mechanisms — see the exemplars.

  Approved proof the mechanism class works, on other targets:
  - **Tafamidis** (2019) — conformational stabiliser: Binds the transthyretin tetramer interface and prevents dissociation — the cleanest approved example of locking a protein into a stable shape
  - **Migalastat** (2016) — conformational stabiliser: Pharmacological chaperone; stabilises mutant alpha-galactosidase A into a foldable, traffickable form
  - **Sirolimus** (1999) — molecular glue: Glues FKBP12 to mTOR; neither interaction occurs without the drug
  - **Tacrolimus** (1994) — molecular glue: Glues FKBP12 to calcineurin

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Cell-type reach

Three things must overlap: where the drug goes, where the gene is expressed, and where the transcription factor is. Any empty intersection kills the route regardless of motif quality or potency.

| TF | Co-expressed with gene | Shared tissues |
|---|---|---|
| NHLH1 | unknown | — |
| ESRRG | unknown | — |
| GATA2 | unknown | — |
| ETV2 | unknown | — |
| MAFG | unknown | — |
| OSR2 | unknown | — |
| NR1I3 | unknown | — |
| PPARD | unknown | — |
| PBX1 | unknown | — |
| ESRRB | unknown | — |

### Human genetic evidence

Human genetic variation in CCR5 associates with: HIV-1 infection, COVID-19, type 1 diabetes mellitus. Read this as a preview of what modulating the gene may do, including the adverse effects — this is the reasoning that anticipated PCSK9 inhibition years before the trials.

| Phenotype | Genetic score | Overall |
|---|---|---|
| HIV-1 infection | 0.635 | 0.629 |
| COVID-19 | 0.238 | 0.414 |
| type 1 diabetes mellitus | 0.517 | 0.339 |
| Erythema nodosum | 0.551 | 0.335 |
| inflammatory bowel disease | 0.484 | 0.309 |
| ulcerative colitis | 0.416 | 0.266 |
| juvenile idiopathic arthritis | 0.421 | 0.265 |
| celiac disease | 0.361 | 0.222 |

- **NHLH1 → CCR5**: To test this link genetically: look for trans-eQTLs of NHLH1 affecting CCR5 in eQTLGen (blood, n>30,000) or GTEx. A hit is strong causal evidence for the regulatory relationship. Absence is weak evidence against it — trans-eQTL detection needs very large samples and most real regulatory links are undetectable at current power.
- **ESRRG → CCR5**: To test this link genetically: look for trans-eQTLs of ESRRG affecting CCR5 in eQTLGen (blood, n>30,000) or GTEx. A hit is strong causal evidence for the regulatory relationship. Absence is weak evidence against it — trans-eQTL detection needs very large samples and most real regulatory links are undetectable at current power.
- **GATA2 → CCR5**: To test this link genetically: look for trans-eQTLs of GATA2 affecting CCR5 in eQTLGen (blood, n>30,000) or GTEx. A hit is strong causal evidence for the regulatory relationship. Absence is weak evidence against it — trans-eQTL detection needs very large samples and most real regulatory links are undetectable at current power.

### Candidates that may beat the current standard

**The bar: strength 0.55, specificity 0.55, efficiency 0.75** (measured) — the best established value on *each axis*, not one drug's scores. Representative incumbent: Maraviroc (direct).

5 of 102 candidates already have a reported link to this gene (CTD returned 152 chemicals with curated expression effects).

*Accounting for all 102 candidates: 5 already linked, 1 novel and non-dominated, 32 novel but beaten on every axis by another candidate, 64 with no advantage. All of them are in results.json under `improvements.all_candidates`.*

#### Compounds already reported to change this gene

These come from CTD's curated chemical–gene interactions: someone measured an expression change and a curator recorded it. That makes them the most concrete evidence in this report — but the direction is what decides whether each is a lead or a warning.

152 chemicals total: **39** move it the way you want, 89 move it the opposite way, 24 unclear.

**Reported to move CCR5 DOWN** — leads:

| Chemical | Predicted by this run? | CTD actions |
|---|---|---|
| (+)-Jq1 Compound | **no** | decreases^expression |
| 1,3-Butadiene | **no** | decreases^expression |
| 2,2-Dichloro-1-(Triphenylphosphonio)Vinyl Formamide Perchlorate | **no** | decreases^expression |
| 3-Chlorophenol | **no** | decreases^expression |
| A-83-01 | **no** | decreases^expression |
| Ammonium 2,3,3,3-Tetrafluoro-2-(Heptafluoropropoxy)-Propanoate | **no** | decreases^expression |
| Ancriviroc | **no** | decreases^expression |
| Anti-Hiv Agents | **no** | decreases^expression |
| Antirheumatic Agents | **no** | decreases^expression |
| Azoxymethane | **no** | decreases^expression |
| Bisphenol S | **no** | decreases^expression |
| Celecoxib | **no** | decreases^expression |
| Chloroprene | **no** | decreases^expression |
| Cirsimarin | **no** | decreases^expression |
| Dextran Sulfate | **no** | decreases^expression |
| Dinoprostone | **no** | decreases^expression |
| Diuron | **no** | decreases^expression |
| Doxorubicin | **no** | decreases^expression |
| Enniatins | **no** | decreases^expression |
| Ethyl Methanesulfonate | **no** | decreases^expression |
| Fenretinide | **no** | decreases^expression |
| Formaldehyde | **no** | decreases^expression |
| Gsk-J4 | **no** | decreases^expression |
| Hydrolyzable Tannins | **no** | decreases^expression |
| Isopentenyl Pyrophosphate | **no** | decreases^expression |
| … 14 more | | |

**53 of these were NOT predicted by this run at all.** They have a measured effect on your gene and no mechanism from this analysis — which makes them worth reading closely, since an empirical hit the model missed is either a route outside its scope or a clue that the scan is wrong.

(+)-Jq1 Compound, 1,3-Butadiene, 2,2-Dichloro-1-(Triphenylphosphonio)Vinyl Formamide Perchlorate, 2,4,5,2',4',5'-Hexachlorobiphenyl, 3-Chlorophenol, A-83-01, Aflatoxin B1, Ammonium 2,3,3,3-Tetrafluoro-2-(Heptafluoropropoxy)-Propanoate, Ancriviroc, Anti-Hiv Agents, Antirheumatic Agents, Azoxymethane, Benzo(A)Pyrene, Bisphenol A, Bisphenol S, Celecoxib, Chloroprene, Cirsimarin, Dextran Sulfate, Dinoprostone, Diuron, Doxorubicin, Drugs, Chinese Herbal, Enniatins, Ethyl Methanesulfonate …

**Reported to move it the WRONG way** — avoid, or use as negative controls:

2,4-Dichlorophenoxyacetic Acid, 2,4-Dinitrofluorobenzene Sulfonic Acid, 2-(2-Amino-3-Methoxyphenyl)-4H-1-Benzopyran-4-One, 2-Chloroethyl Ethyl Sulfide, 2-Naphthylamine, 4-Hydroxyphenyl 4-Isopropoxyphenylsulfone, Acetamide, Air Pollutants, Amphetamine, Amphotericin B, Deoxycholate Drug Combination, Arsenic, Asbestos, Crocidolite, Atorvastatin, Bardoxolone Methyl, Benzalkonium Compounds, Benzene, Biological Factors, Bisphenol F, Bleomycin, Bromine, Cadmium, Cadmium Chloride, Cannabidiol, Captax, Carbon Tetrachloride …

*CTD aggregates across cell types, doses and assays. A 'both' or contradictory annotation usually means the effect is context-dependent rather than that the data are wrong — check the underlying study before relying on any single row.*

#### Not previously linked, and better on at least one axis

| Candidate | Route | TF | Strength | Specificity | Efficiency | Wins on | Trade-offs |
|---|---|---|---|---|---|---|---|
| Calcitriol | via TF | VDR | 0.677 | 1.0 | 0.99 | higher specificity, superior efficiency | none |

*These are predicted scores compared against measured ones, so a margin was required before counting a win. They are hypotheses to test, not improvements demonstrated.*

Novel but beaten on every axis by another candidate (shown so nothing is hidden): Seladelpar, Seladelpar Lysine, Elafibranor, Tamibarotene, Elacestrant, Fulvestrant, Elacestrant Hydrochloride, Mifepristone, Calcipotriene, Trifarotene, Cholecalciferol, Belzutifan ...

#### Non-dominated set (nothing beats these on every axis)

- **Calcitriol** (NEW link) — beats the incumbent on higher specificity, superior efficiency

### Predicted transcriptional outcome

P(gene moves) = P(engages) x P(direction) x P(activates here) x P(responsive). Four uncertain links multiply, so these numbers are small by construction — a chain at 0.7 per link lands near 0.24. Anything reporting 0.9 here has hidden a link.

| Drug | TF | P(engage) | P(direction) | P(promoter) | P(responsive) | **P(gene moves)** | Verdict | Limiting step |
|---|---|---|---|---|---|---|---|---|
| Calcitriol ★ | VDR | 0.95 | 0.75 | 0.478 | 0.4 | **0.136** (0.048–0.225) | unlikely via this route | gene responsiveness in your cells |
| Seladelpar | PPARD | 0.95 | 0.75 | 0.58 | 0.4 | **0.165** (0.083–0.248) | speculative | gene responsiveness in your cells |
| Mifepristone | PGR | 0.95 | 0.75 | 0.576 | 0.4 | **0.164** (0.082–0.246) | speculative | gene responsiveness in your cells |
| Fulvestrant | ESR1 | 0.95 | 0.75 | 0.52 | 0.4 | **0.148** (0.052–0.245) | unlikely via this route | gene responsiveness in your cells |
| Belzutifan | HLF | 0.887 | 0.75 | 0.54 | 0.4 | **0.144** (0.05–0.237) | unlikely via this route | gene responsiveness in your cells |
| Diethylstilbestrol | ESRRG | 0.775 | 0.75 | 0.611 | 0.4 | **0.142** (0.071–0.213) | unlikely via this route | gene responsiveness in your cells |
| Tamibarotene | RARA | 0.811 | 0.75 | 0.561 | 0.4 | **0.137** (0.068–0.205) | unlikely via this route | gene responsiveness in your cells |
| Calcipotriene | VDR | 0.95 | 0.75 | 0.478 | 0.4 | **0.136** (0.048–0.225) | unlikely via this route | gene responsiveness in your cells |
| Seladelpar Lysine | PPARD | 0.775 | 0.75 | 0.58 | 0.4 | **0.135** (0.067–0.202) | unlikely via this route | gene responsiveness in your cells |
| Bezafibrate | PPARD | 0.775 | 0.75 | 0.58 | 0.4 | **0.135** (0.067–0.202) | unlikely via this route | gene responsiveness in your cells |
| Tretinoin | PPARD | 0.775 | 0.75 | 0.58 | 0.4 | **0.135** (0.067–0.202) | unlikely via this route | gene responsiveness in your cells |
| Clotrimazole | NR1I3 | 0.775 | 0.75 | 0.582 | 0.4 | **0.135** (0.068–0.203) | unlikely via this route | gene responsiveness in your cells |
| Meclizine | NR1I3 | 0.775 | 0.75 | 0.582 | 0.4 | **0.135** (0.068–0.203) | unlikely via this route | gene responsiveness in your cells |
| Diethylstilbestrol | ESRRB | 0.775 | 0.75 | 0.58 | 0.4 | **0.135** (0.067–0.202) | unlikely via this route | gene responsiveness in your cells |
| Trifarotene | RARG | 0.933 | 0.75 | 0.477 | 0.4 | **0.134** (0.047–0.221) | unlikely via this route | gene responsiveness in your cells |
| Cholecalciferol | VDR | 0.929 | 0.75 | 0.478 | 0.4 | **0.133** (0.047–0.22) | unlikely via this route | gene responsiveness in your cells |
| Acitretin | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Etretinate | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Isotretinoin | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Trifarotene | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Adapalene | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Alitretinoin | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Tazarotene | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Tretinoin | RARA | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.196) | unlikely via this route | gene responsiveness in your cells |
| Elafibranor | PPARD | 0.73 | 0.75 | 0.58 | 0.4 | **0.127** (0.064–0.191) | unlikely via this route | gene responsiveness in your cells |

★ = named in the comparison against the standard of care above, pinned here so the two sections cannot disagree.

**Most informative next experiment** (for Calcitriol): Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.

### Virtual screen of approved drugs against PPARD

Binding model, CV AUC 0.999 — trained on 1264 measured actives vs 3309 decoys, screened 3311 approved drugs.

*Direction not annotated: only 5 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Seladelpar | 0.988 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Elafibranor | 0.898 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Pemafibrate | 0.798 | 0.0 |  |
| predicted binder, close to known chemistry | Seladelpar Lysine | 0.73 | 0.0 |  |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against PGR

Binding model, CV AUC 0.995 — trained on 1647 measured actives vs 3239 decoys, screened 3311 approved drugs.

*Direction not annotated: only 27 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Mifepristone | 0.994 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Drospirenone | 0.946 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Fluocinolone Acetonide | 0.945 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Ulipristal Acetate | 0.924 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Norethindrone | 0.906 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Fluticasone Propionate | 0.886 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Halobetasol Propionate | 0.885 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Chlormadinone Acetate | 0.881 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Bicalutamide | 0.874 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Fluticasone Furoate | 0.872 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Triamcinolone Acetonide | 0.869 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Desonide | 0.867 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Norethindrone Acetate | 0.866 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Canrenone | 0.863 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Medroxyprogesterone Acetate | 0.858 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Megestrol Acetate | 0.848 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Flurandrenolide | 0.847 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Clobetasol Propionate | 0.845 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Flunisolide | 0.843 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Halcinonide | 0.84 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Beclomethasone Dipropionate | 0.834 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Etonogestrel | 0.832 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Testosterone | 0.83 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Mometasone Furoate | 0.82 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Mestranol | 0.82 | 0.0 | binds this receptor with measured activity, and has no repor |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against RARA

Binding model, CV AUC 0.997 — trained on 211 measured actives vs 633 decoys, screened 3311 approved drugs.

*Direction not annotated: only 2 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.988 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamibarotene | 0.984 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Trifarotene | 0.956 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.937 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.937 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Adapalene | 0.902 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tazarotene | 0.819 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.937 | 0.0 |  |
| predicted binder, close to known chemistry | Palovarotene | 0.864 | 0.429 |  |
| no call | Acoramidis | 0.563 | 0.681 | outside applicability domain; score withheld |
| no call | Tafamidis Meglumine | 0.559 | 0.727 | outside applicability domain; score withheld |
| no call | Deferasirox | 0.552 | 0.673 | outside applicability domain; score withheld |
| no call | Acoramidis Hydrochloride | 0.547 | 0.686 | outside applicability domain; score withheld |
| no call | Iptacopan | 0.51 | 0.726 | outside applicability domain; score withheld |
| no call | Iptacopan Hydrochloride | 0.491 | 0.733 | outside applicability domain; score withheld |
| no call | Eprosartan | 0.489 | 0.741 | outside applicability domain; score withheld |
| no call | Probenecid | 0.458 | 0.69 | outside applicability domain; score withheld |
| no call | Bentiromide | 0.453 | 0.661 | outside applicability domain; score withheld |
| no call | Eprosartan Mesylate | 0.445 | 0.756 | outside applicability domain; score withheld |
| no call | Eltrombopag | 0.393 | 0.736 | outside applicability domain; score withheld |
| no call | Eltrombopag Olamine | 0.368 | 0.755 | outside applicability domain; score withheld |
| no call | Pranlukast | 0.366 | 0.727 | outside applicability domain; score withheld |
| no call | Novobiocin | 0.353 | 0.76 | outside applicability domain; score withheld |
| no call | Opicapone | 0.351 | 0.733 | outside applicability domain; score withheld |
| no call | Febuxostat | 0.343 | 0.689 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against HLF

Binding model, CV AUC 0.997 — trained on 273 measured actives vs 819 decoys, screened 3311 approved drugs.

*Direction not annotated: only 1 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Belzutifan | 0.995 | 0.0 | binds this receptor with measured activity, and has no repor |
| no call | Doravirine | 0.459 | 0.694 | outside applicability domain; score withheld |
| no call | Nimesulide | 0.367 | 0.7 | outside applicability domain; score withheld |
| no call | Tucatinib | 0.361 | 0.722 | outside applicability domain; score withheld |
| no call | Lenacapavir | 0.324 | 0.81 | outside applicability domain; score withheld |
| no call | Lenacapavir Sodium | 0.321 | 0.827 | outside applicability domain; score withheld |
| no call | Crisaborole | 0.319 | 0.739 | outside applicability domain; score withheld |
| no call | Laropiprant | 0.307 | 0.658 | outside applicability domain; score withheld |
| no call | Oteseconazole | 0.293 | 0.795 | outside applicability domain; score withheld |
| no call | Bicalutamide | 0.286 | 0.662 | outside applicability domain; score withheld |
| no call | Etravirine | 0.27 | 0.746 | outside applicability domain; score withheld |
| no call | Venetoclax | 0.261 | 0.775 | outside applicability domain; score withheld |
| no call | Lapatinib | 0.259 | 0.683 | outside applicability domain; score withheld |
| no call | Rucaparib | 0.251 | 0.788 | outside applicability domain; score withheld |
| no call | Tavaborole | 0.25 | 0.803 | outside applicability domain; score withheld |
| no call | Rucaparib Camsylate | 0.246 | 0.709 | outside applicability domain; score withheld |
| no call | Levocabastine | 0.242 | 0.722 | outside applicability domain; score withheld |
| no call | Levocabastine Hydrochloride | 0.238 | 0.708 | outside applicability domain; score withheld |
| no call | Cilostazol | 0.228 | 0.742 | outside applicability domain; score withheld |
| no call | Regorafenib | 0.228 | 0.686 | outside applicability domain; score withheld |
| no call | Pretomanid | 0.224 | 0.727 | outside applicability domain; score withheld |
| no call | Ivosidenib | 0.223 | 0.777 | outside applicability domain; score withheld |
| no call | Triclosan | 0.219 | 0.708 | outside applicability domain; score withheld |
| no call | Lapatinib Ditosylate | 0.216 | 0.688 | outside applicability domain; score withheld |
| no call | Omarigliptin | 0.208 | 0.742 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against ESR1

Binding model, CV AUC 0.998 — trained on 2014 measured actives vs 3294 decoys, screened 3311 approved drugs.

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Raloxifene | 0.995 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Bazedoxifene | 0.987 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Raloxifene Hydrochloride | 0.967 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Lasofoxifene | 0.959 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Ethinyl Estradiol | 0.953 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Diethylstilbestrol | 0.927 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estradiol | 0.919 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Fulvestrant | 0.912 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estrone | 0.896 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estriol | 0.829 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamoxifen | 0.802 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estetrol Anhydrous | 0.785 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Mestranol | 0.543 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tibolone | 0.49 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Norethindrone | 0.437 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Ethylestrenol | 0.42 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Stanozolol | 0.385 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Hydroquinone | 0.886 | 0.0 | predicted agonist (p=0.635, low confidence) |
| predicted binder, close to known chemistry | Bazedoxifene Acetate | 0.886 | 0.0 |  |
| predicted binder, close to known chemistry | Dienestrol | 0.821 | 0.0 | predicted agonist (p=0.635, low confidence) |
| predicted binder, close to known chemistry | Hexestrol | 0.805 | 0.0 | predicted agonist (p=0.599, low confidence) |
| predicted binder, close to known chemistry | Estetrol | 0.77 | 0.0 | predicted agonist (p=0.741, moderate confidence) |
| predicted binder, close to known chemistry | Mequinol | 0.743 | 0.0 | predicted agonist (p=0.606, low confidence) |
| predicted binder, close to known chemistry | Monobenzone | 0.73 | 0.0 |  |
| predicted binder, close to known chemistry | Lasofoxifene Tartrate | 0.728 | 0.0 |  |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against RXRA

Binding model, CV AUC 0.998 — trained on 1021 measured actives vs 3063 decoys, screened 3311 approved drugs.

*Direction not annotated: only 1 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.999 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.896 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.896 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Levothyroxine | 0.37 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.896 | 0.0 |  |
| predicted binder, close to known chemistry | Tamibarotene | 0.805 | 0.0 |  |
| predicted binder, close to known chemistry | Oxaprozin | 0.802 | 0.0 |  |
| NOVEL HYPOTHESIS | Palovarotene | 0.859 | 0.516 | extrapolation: nearest known ligand Tanimoto 0.48. The model |
| no call | Celecoxib | 0.19 | 0.684 | outside applicability domain; score withheld |
| no call | Lurbinectedin | 0.174 | 0.778 | outside applicability domain; score withheld |
| no call | Omidenepag Isopropyl | 0.173 | 0.758 | outside applicability domain; score withheld |
| no call | Trabectedin | 0.155 | 0.763 | outside applicability domain; score withheld |
| no call | Pralsetinib | 0.145 | 0.778 | outside applicability domain; score withheld |
| no call | Gemfibrozil | 0.123 | 0.651 | outside applicability domain; score withheld |
| no call | Pemafibrate | 0.118 | 0.683 | outside applicability domain; score withheld |
| no call | Axitinib | 0.093 | 0.78 | outside applicability domain; score withheld |
| no call | Doravirine | 0.093 | 0.761 | outside applicability domain; score withheld |
| no call | Laropiprant | 0.082 | 0.725 | outside applicability domain; score withheld |
| no call | Entrectinib | 0.055 | 0.704 | outside applicability domain; score withheld |
| no call | Oxamniquine | 0.051 | 0.738 | outside applicability domain; score withheld |
| no call | Ceritinib | 0.049 | 0.653 | outside applicability domain; score withheld |
| no call | Ecamsule | 0.049 | 0.709 | outside applicability domain; score withheld |
| no call | Chlorzoxazone | 0.03 | 0.729 | outside applicability domain; score withheld |
| no call | Molsidomine | 0.03 | 0.737 | outside applicability domain; score withheld |
| no call | Lorlatinib | 0.029 | 0.717 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against THRB

Binding model, CV AUC 1.0 — trained on 33 measured actives vs 99 decoys, screened 3311 approved drugs.

*Direction not annotated: only 6 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

**Screen rejected:** no known actives appear among the top hits — the model is not ranking them above the library, so the novel list cannot be trusted No candidates are reported from it.

### Virtual screen of approved drugs against VDR

Binding model, CV AUC 0.996 — trained on 370 measured actives vs 1110 decoys, screened 3311 approved drugs.

*Direction not annotated: only 18 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Calcipotriene | 0.995 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Calcitriol | 0.991 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Doxercalciferol | 0.989 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tacalcitol | 0.982 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Cholecalciferol | 0.94 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Calcifediol Anhydrous | 0.989 | 0.183 |  |
| predicted binder, close to known chemistry | Calcifediol | 0.984 | 0.197 |  |
| predicted binder, close to known chemistry | Paricalcitol | 0.974 | 0.317 |  |
| predicted binder, close to known chemistry | Ergocalciferol | 0.949 | 0.233 |  |
| predicted binder, close to known chemistry | Alfacalcidol | 0.932 | 0.0 |  |
| predicted binder, close to known chemistry | Dihydrotachysterol | 0.873 | 0.339 |  |
| no call | Omacetaxine Mepesuccinate | 0.214 | 0.792 | outside applicability domain; score withheld |
| no call | Paritaprevir | 0.164 | 0.699 | outside applicability domain; score withheld |
| no call | Pafolacianine | 0.154 | 0.707 | outside applicability domain; score withheld |
| no call | Sildenafil Citrate | 0.149 | 0.714 | outside applicability domain; score withheld |
| no call | Ponesimod | 0.147 | 0.695 | outside applicability domain; score withheld |
| no call | Dalfopristin | 0.146 | 0.829 | outside applicability domain; score withheld |
| no call | Pafolacianine Sodium | 0.146 | 0.699 | outside applicability domain; score withheld |
| no call | Fusidic Acid | 0.145 | 0.753 | outside applicability domain; score withheld |
| no call | Sunitinib Malate | 0.142 | 0.714 | outside applicability domain; score withheld |
| no call | Elagolix | 0.142 | 0.763 | outside applicability domain; score withheld |
| no call | Crinecerfont | 0.141 | 0.8 | outside applicability domain; score withheld |
| no call | Sunitinib | 0.14 | 0.667 | outside applicability domain; score withheld |
| no call | Carbenoxolone | 0.139 | 0.737 | outside applicability domain; score withheld |
| no call | Elagolix Sodium | 0.138 | 0.769 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against RARG

Binding model, CV AUC 0.998 — trained on 246 measured actives vs 738 decoys, screened 3311 approved drugs.

*Direction not annotated: only 2 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.983 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamibarotene | 0.978 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.961 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.961 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Adapalene | 0.959 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Trifarotene | 0.95 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tazarotene | 0.804 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.961 | 0.0 |  |
| predicted binder, close to known chemistry | Palovarotene | 0.89 | 0.429 |  |
| no call | Tafamidis Meglumine | 0.571 | 0.727 | outside applicability domain; score withheld |
| no call | Deferasirox | 0.562 | 0.672 | outside applicability domain; score withheld |
| no call | Acoramidis | 0.549 | 0.681 | outside applicability domain; score withheld |
| no call | Acoramidis Hydrochloride | 0.535 | 0.686 | outside applicability domain; score withheld |
| no call | Iptacopan | 0.522 | 0.726 | outside applicability domain; score withheld |
| no call | Eprosartan | 0.521 | 0.735 | outside applicability domain; score withheld |
| no call | Iptacopan Hydrochloride | 0.504 | 0.733 | outside applicability domain; score withheld |
| no call | Eprosartan Mesylate | 0.443 | 0.721 | outside applicability domain; score withheld |
| no call | Probenecid | 0.442 | 0.69 | outside applicability domain; score withheld |
| no call | Febuxostat | 0.415 | 0.689 | outside applicability domain; score withheld |
| no call | Eltrombopag | 0.383 | 0.736 | outside applicability domain; score withheld |
| no call | Bentiromide | 0.375 | 0.661 | outside applicability domain; score withheld |
| no call | Eltrombopag Olamine | 0.351 | 0.755 | outside applicability domain; score withheld |
| no call | Telmisartan | 0.35 | 0.737 | outside applicability domain; score withheld |
| no call | Opicapone | 0.312 | 0.733 | outside applicability domain; score withheld |
| no call | Ledipasvir | 0.294 | 0.663 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against RARB

Binding model, CV AUC 0.998 — trained on 272 measured actives vs 816 decoys, screened 3311 approved drugs.

*Direction not annotated: only 1 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.99 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamibarotene | 0.982 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Adapalene | 0.96 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Trifarotene | 0.952 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.936 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.936 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tazarotene | 0.849 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.936 | 0.0 |  |
| predicted binder, close to known chemistry | Palovarotene | 0.897 | 0.429 |  |
| no call | Acoramidis | 0.597 | 0.672 | outside applicability domain; score withheld |
| no call | Acoramidis Hydrochloride | 0.574 | 0.677 | outside applicability domain; score withheld |
| no call | Iptacopan | 0.507 | 0.716 | outside applicability domain; score withheld |
| no call | Deferasirox | 0.506 | 0.673 | outside applicability domain; score withheld |
| no call | Iptacopan Hydrochloride | 0.481 | 0.707 | outside applicability domain; score withheld |
| no call | Eprosartan Mesylate | 0.466 | 0.672 | outside applicability domain; score withheld |
| no call | Bentiromide | 0.415 | 0.661 | outside applicability domain; score withheld |
| no call | Telmisartan | 0.409 | 0.685 | outside applicability domain; score withheld |
| no call | Eltrombopag | 0.379 | 0.736 | outside applicability domain; score withheld |
| no call | Febuxostat | 0.376 | 0.689 | outside applicability domain; score withheld |
| no call | Pranlukast | 0.338 | 0.727 | outside applicability domain; score withheld |
| no call | Eltrombopag Olamine | 0.335 | 0.753 | outside applicability domain; score withheld |
| no call | Ledipasvir | 0.329 | 0.663 | outside applicability domain; score withheld |
| no call | Lusutrombopag | 0.293 | 0.762 | outside applicability domain; score withheld |
| no call | Opicapone | 0.288 | 0.733 | outside applicability domain; score withheld |
| no call | Alogliptin Benzoate | 0.277 | 0.75 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against THRA

Binding model, CV AUC 1.0 — trained on 33 measured actives vs 99 decoys, screened 3311 approved drugs.

*Direction not annotated: only 6 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

**Screen rejected:** no known actives appear among the top hits — the model is not ranking them above the library, so the novel list cannot be trusted No candidates are reported from it.

### Drug–gene links with no recorded direction

747 interaction(s) across 73 factor(s) come from curated sources that record a drug–gene relationship **without saying which way it goes**. They are held back from the recommendations above, because the entire question here is direction and an undirected link cannot answer it.

They are still leads — someone curated a real relationship — but they need a literature check, not a prediction. Many are chemotherapeutics or broad stressors whose association with a transcription factor reflects cell-wide perturbation rather than a specific mechanism.

| TF | Drug | Reported by | DGIdb score |
|---|---|---|---|
| ELK1 | Ganciclovir | NCI | 13.23 |
| SMAD4 | Lysine | CIViC | 10.59 |
| NR1I2 | Rifaximin | DTC | 4.41 |
| TFAP2A | Amprenavir | NCI | 4.41 |
| EHF | Ziprasidone | PharmGKB | 4.23 |
| GABPA | Omaveloxolone | TTD | 3.78 |
| MAFK | Rifampin | PharmGKB | 3.65 |
| GATA4 | Acenocoumarol | PharmGKB | 3.53 |
| PGR | Ulipristal Acetate | ChEMBL; GuideToPharmacology; TEND; TTD; TdgClinicalTrial | 3.48 |
| NFATC4 | Anhydrous Tacrolimus | PharmGKB | 3.41 |
| SOX10 | Vemurafenib | CIViC | 3.11 |
| PROX1 | Atenolol | PharmGKB | 2.94 |
| SMAD2 | Amifostine Anhydrous | NCI | 2.65 |
| NEUROD1 | Repaglinide | PharmGKB | 2.52 |
| NEUROD1 | Deferoxamine | NCI | 2.35 |
| ETV1 | Trametinib Dimethyl Sulfoxide | CIViC | 2.3 |
| ETV5 | Trametinib Dimethyl Sulfoxide | CIViC | 2.3 |
| KLF4 | Hydroxyurea | PharmGKB | 2.3 |
| TCF7L2 | Linagliptin | PharmGKB | 2.21 |
| NR1I3 | Oxymorphone | DTC | 2.08 |
| HNF1A | Rosuvastatin | PharmGKB | 2.04 |
| SREBF1 | Insulin, Regular, Human | NCI | 2.04 |
| STAT5B | Tofacitinib | DTC | 1.92 |
| ATF3 | Progesterone | NCI | 1.74 |
| PLAGL1 | Olanzapine | PharmGKB | 1.56 |
| … 722 more | | | |

*Pass `--include-unknown-direction` to fold these into the main tables anyway — but every such row will read UNKNOWN / unclear, which is what they are.*

### Approved drugs acting on these factors

| TF | TF role | Do what to this TF | Drug | Action | Drug direction | First approval | Source |
|---|---|---|---|---|---|---|---|
| ESRRG | activator | decrease | Diethylstilbestrol | INHIBITOR | down | — | DGIdb |
| NR1I3 | activator | decrease | Clotrimazole | INHIBITOR | down | — | DGIdb |
| NR1I3 | activator | decrease | Meclizine | INHIBITOR | down | — | DGIdb |
| PPARD | repressor | increase | Bezafibrate | AGONIST | up | — | DGIdb |
| PPARD | repressor | increase | Elafibranor | AGONIST | up | 2024 | ChEMBL |
| PPARD | repressor | increase | Seladelpar | AGONIST | up | 2024 | ChEMBL |
| PPARD | repressor | increase | Seladelpar Lysine | AGONIST | up | 2024 | ChEMBL |
| PPARD | repressor | increase | Tretinoin | AGONIST | up | — | DGIdb |
| ESRRB | activator | decrease | Diethylstilbestrol | INHIBITOR | down | — | DGIdb |
| PGR | activator | decrease | Mifepristone | ANTAGONIST | down | 2000 | ChEMBL |
| RARA | repressor | increase | Acitretin | AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Adapalene | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Alitretinoin | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Etretinate | AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Isotretinoin | AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Tamibarotene | AGONIST | up | — | ChEMBL |
| RARA | repressor | increase | Tazarotene | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Tretinoin | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Trifarotene | AGONIST | up | — | DGIdb |
| PPARG | dual | decrease | Bexarotene | INHIBITOR | down | — | DGIdb |
| PPARG | dual | decrease | Diclofenac Sodium | INHIBITOR | down | — | DGIdb |
| PPARG | dual | decrease | Resveratrol | INHIBITOR | down | — | DGIdb |
| HLF | activator | decrease | Belzutifan | INHIBITOR | down | 2021 | ChEMBL |
| ESR1 | activator | decrease | Bazedoxifene | INHIBITOR | down | — | DGIdb |
| ESR1 | activator | decrease | Clomiphene | INHIBITOR | down | — | DGIdb |
| … 76 more, see tf_summary.tsv | | | | | | | |

## Promoter of CCR5 (ENST00000445772)

- Transcript chosen: **distinct TSS**
- Assembly: GRCh38  |  Strand: +
- TSS: chr3:46372670
- Window: chr3:46367670-46373669 (-5000 / +1000, 6000 bp)
- Motifs scanned: 879 at p ≤ 0.0001 per position
- Non-redundant sites: 953 across 420 TFs

> **⚠ Weak scan.** 953 sites against 774 expected by chance over 4401 scannable bases. Most of this list is noise; treat only the highest-scoring sites with ChIP support as real.

### Ligand-activated (nuclear) receptors with sites here

These are the directly druggable ones.

| TF | Role | Sites | Best rel. score | Best site (vs TSS) | ChIP peaks | Score |
|---|---|---|---|---|---|---|
| RARA | repressor | 9 | 0.907 | GGGTCAAGCATGGTGCA (-108) | 0 | 72.2 |
| ESRRG | activator | 3 | 0.987 | CCAAGGTCA (-2723) | 0 | 64.6 |
| ESRRB | activator | 3 | 0.95 | CCAAGGTCAC (-2723) | 0 | 63.5 |
| RORB | activator | 2 | 0.974 | AACTGGGTCA (-97) | 0 | 62.7 |
| NR1I3 | activator | 2 | 1.0 | TGAACTTT (-4415) | 0 | 61.7 |
| RORC | repressor | 2 | 0.988 | AACTGGGTCA (-97) | 0 | 61.6 |
| RORA | unknown | 3 | 0.966 | AAACTGGGTCA (-97) | 0 | 60.8 |
| RXRG | unknown | 4 | 0.907 | GGGTCAAGCATGGTGCA (-108) | 0 | 58.7 |
| PPARD | repressor | 2 | 0.838 | AGGGTCACAGTTTG (-1541) | 0 | 55.6 |
| PPARG | dual | 4 | 0.871 | AGAGCCAAGGTCA (-2727) | 0 | 55.1 |
| NR1H2 | activator | 3 | 1.0 | AGGTCA (-2720) | 0 | 55.0 |
| NR2E3 | activator | 3 | 1.0 | AAGCTT (-2865) | 0 | 55.0 |
| NR1H3 | activator | 8 | 0.873 | TGAACTTAAATAGACT (-1041) | 0 | 52.3 |
| NR2C2 | activator | 2 | 0.859 | CAGGTCAGAGATGG (+451) | 0 | 52.3 |
| NR4A2 | activator | 2 | 0.948 | GAAGGTCA (-835) | 0 | 52.0 |
| ESR1 | activator | 1 | 0.778 | GGGCTTCCGTGACCT (-2720) | 0 | 51.7 |
| NR2F1 | repressor | 5 | 0.921 | CAAGGTCACGGA (-2722) | 0 | 51.1 |
| HNF4G | activator | 1 | 0.936 | CAAAGTCCC (+514) | 0 | 50.6 |
| HNF4A | activator | 1 | 0.93 | CAAAGTCCC (+514) | 0 | 50.4 |
| NR2F2 | repressor | 2 | 1.0 | AAGGTCA (-2721) | 0 | 50.2 |
| PPARA | dual | 5 | 1.0 | AAGGTCA (-2721) | 0 | 49.0 |
| NR1D2 | repressor | 3 | 0.837 | TGGGTTAGTCTGCCT (-2824) | 0 | 48.6 |
| NR2F6 | repressor | 4 | 0.836 | GTGATCAAAAGTTCA (-4422) | 0 | 48.6 |
| NR4A1 | activator | 1 | 0.929 | GAAGGTCA (-835) | 0 | 48.1 |
| NR5A1 | activator | 2 | 0.856 | GAGCCAAGGTCA (-2726) | 0 | 47.4 |

### Top-ranked transcription factors

| Rank | TF | Role | NR | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|
| 1 | RARA | repressor | yes | 9 | 0.907 | -108 | 0 | 9 | 72.2 |
| 2 | ESRRG | activator | yes | 3 | 0.987 | -2723 | 0 | 1 | 64.6 |
| 3 | ESRRB | activator | yes | 3 | 0.95 | -2723 | 0 | 1 | 63.5 |
| 4 | RORB | activator | yes | 2 | 0.974 | -97 | 0 | 0 | 62.7 |
| 5 | NFIA | activator |  | 3 | 1.0 | +291 | 0 | 0 | 62.2 |
| 6 | NR1I3 | activator | yes | 2 | 1.0 | -4415 | 0 | 2 | 61.7 |
| 7 | RORC | repressor | yes | 2 | 0.988 | -97 | 0 | 0 | 61.6 |
| 8 | CTCF | repressor |  | 7 | 0.782 | +18 | 0 | 0 | 61.3 |
| 9 | RORA | unknown | yes | 3 | 0.966 | -97 | 0 | 1 | 60.8 |
| 10 | PRDM5 | repressor |  | 1 | 0.954 | +16 | 0 | 0 | 59.9 |
| 11 | SOX10 | activator |  | 4 | 1.0 | -469 | 0 | 0 | 59.6 |
| 12 | THAP1 | repressor |  | 2 | 0.962 | +164 | 0 | 0 | 59.0 |
| 13 | RXRG | unknown | yes | 4 | 0.907 | -108 | 0 | 4 | 58.7 |
| 14 | HOXA10 | activator |  | 1 | 0.97 | +122 | 0 | 0 | 58.5 |
| 15 | STAT6 | activator |  | 3 | 0.993 | +570 | 0 | 0 | 58.4 |
| 16 | TFAP4 | activator |  | 8 | 0.88 | +327 | 0 | 0 | 58.0 |
| 17 | DLX5 | activator |  | 1 | 0.94 | +122 | 0 | 0 | 57.6 |
| 18 | ETV5 | activator |  | 4 | 0.863 | +327 | 0 | 0 | 57.5 |
| 19 | FIGLA | activator |  | 3 | 0.863 | +327 | 0 | 0 | 57.5 |
| 20 | DLX2 | activator |  | 1 | 0.934 | +122 | 0 | 0 | 57.4 |
| 21 | ONECUT1 | activator |  | 2 | 0.953 | +298 | 0 | 0 | 57.4 |
| 22 | TCF12 | activator |  | 6 | 1.0 | +761 | 0 | 0 | 57.2 |
| 23 | NFAT5 | activator |  | 3 | 1.0 | +776 | 0 | 0 | 57.2 |
| 24 | HOXD9 | repressor |  | 1 | 0.975 | +122 | 0 | 0 | 57.1 |
| 25 | CDX1 | activator |  | 1 | 0.924 | +121 | 0 | 0 | 57.1 |

### Every route that acts on this gene

150 routes across 3 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

- 1 × clinical — approved for this gene, documented in humans
- 148 × curated — mechanism annotated by a human curator
- 1 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 1 | direct: gene product | Maraviroc | CCR5 | 0 | C-C chemokine receptor type 5 antagonist | 0.55 | 0.55 | 0.75 | 60.0 |
| 3 | TF: orthosteric ligand | Calcitriol | VDR | 1 | ligand at the transcription factor | 0.856 | 1.0 | 0.99 | 75.2 |
| 3 | TF: orthosteric ligand | Calcipotriene | VDR | 1 | ligand at the transcription factor | 0.867 | 1.0 | 0.925 | 74.2 |
| 3 | TF: orthosteric ligand | Seladelpar | PPARD | 1 | ligand at the transcription factor | 0.855 | 1.0 | 0.933 | 74.0 |
| 3 | TF: orthosteric ligand | Fulvestrant | ESR1 | 1 | ligand at the transcription factor | 0.842 | 0.926 | 0.967 | 72.2 |
| 3 | TF: orthosteric ligand | Trifarotene | RARG | 1 | ligand at the transcription factor | 0.797 | 1.0 | 0.865 | 70.8 |
| 3 | TF: orthosteric ligand | Cholecalciferol | VDR | 1 | ligand at the transcription factor | 0.796 | 1.0 | 0.862 | 70.7 |
| 3 | TF: orthosteric ligand | Seladelpar Lysine | PPARD | 1 | ligand at the transcription factor | 0.664 | 1.0 | 0.725 | 63.7 |
| 3 | TF: orthosteric ligand | Elacestrant | ESR1 | 1 | ligand at the transcription factor | 0.675 | 0.926 | 0.794 | 63.4 |
| 3 | TF: orthosteric ligand | Alfacalcidol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Calcifediol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Dihydrotachysterol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Ergocalciferol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Paricalcitol | VDR | 1 | ligand at the transcription factor | 0.642 | 1.0 | 0.725 | 63.0 |
| 3 | TF: orthosteric ligand | Dextrothyroxine | THRB | 1 | ligand at the transcription factor | 0.639 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Dextrothyroxine Sodium | THRB | 1 | ligand at the transcription factor | 0.639 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Liothyronine | THRB | 1 | ligand at the transcription factor | 0.639 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Dextrothyroxine | THRA | 1 | ligand at the transcription factor | 0.637 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Dextrothyroxine Sodium | THRA | 1 | ligand at the transcription factor | 0.637 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Liothyronine | THRA | 1 | ligand at the transcription factor | 0.637 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Levothyroxine Sodium | THRB | 1 | ligand at the transcription factor | 0.639 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Liothyronine Sodium | THRB | 1 | ligand at the transcription factor | 0.639 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Palovarotene | RARG | 1 | ligand at the transcription factor | 0.639 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Levothyroxine Sodium | THRA | 1 | ligand at the transcription factor | 0.637 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Liothyronine Sodium | THRA | 1 | ligand at the transcription factor | 0.637 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Danthron | RXRA | 1 | ligand at the transcription factor | 0.637 | 1.0 | 0.725 | 62.9 |
| 3 | TF: orthosteric ligand | Resmetirom | THRB | 1 | ligand at the transcription factor | 0.603 | 1.0 | 0.731 | 61.9 |
| 3 | TF: orthosteric ligand | Elacestrant Hydrochloride | ESR1 | 1 | ligand at the transcription factor | 0.639 | 0.926 | 0.725 | 60.9 |
| 3 | TF: orthosteric ligand | Tamibarotene | RARA | 1 | ligand at the transcription factor | 0.742 | 0.741 | 0.795 | 60.4 |
| 3 | TF: orthosteric ligand | Doxercalciferol | VDR | 1 | ligand at the transcription factor | 0.531 | 1.0 | 0.627 | 57.5 |

Best in each category:

- **direct: gene product** — Maraviroc (tier 1, score 60.0)
- **TF: orthosteric ligand** — Calcitriol (tier 3, score 75.2)
- **virtual screen hypothesis** — Palovarotene (tier 5, score 51.5)

*No route found in: direct: gene or transcript, TF: allosteric modulator, TF: stabiliser or glue, TF: endogenous ligand, TF: degradation blocker, upstream signalling, empirical signature, chromatin-level. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

**direct: gene product** (1 found)

- Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- *Why:* Immediate, well-understood pharmacology; the protein is right there rather than three inferential steps away.
- *Cost:* Usually changes the protein's *activity*, not its abundance. If your readout is a western blot, it can look like nothing happened.
- *Differs:* Does not touch expression at all — it is here because it changes functional output, which is often what you actually want.

**TF: orthosteric ligand** (148 found)

- Binds the transcription factor's natural ligand pocket.
- *Why:* The best-characterised route, and where nearly all approved TF pharmacology lives.
- *Cost:* Only ~48 human TFs (the nuclear receptors) have such a pocket, and the factor's whole regulon moves with your gene.
- *Differs:* Clamps the receptor on or off regardless of endogenous signalling — contrast the allosteric and ligand-axis routes, which scale what is already happening.

**virtual screen hypothesis** (1 found)

- An approved drug predicted to bind one of your TFs, never tested against it.
- *Why:* Real repurposing space — most approved drugs have been assayed against a tiny fraction of the proteome.
- *Cost:* An extrapolation beyond the training chemistry. Expect most to be wrong; each is a binding assay, not a finding.
- *Differs:* The only category with no prior evidence of the drug–target pair existing at all.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **TF: allosteric modulator** — Binds a site away from the ligand pocket, amplifying or damping the endogenous ligand's effect (PAM/NAM).
- **TF: stabiliser or glue** — Locks the factor in an active conformation, or glues it to a partner protein.
- **TF: endogenous ligand** — Changes how much natural ligand exists, by drugging the enzyme that makes or destroys it.
- **TF: degradation blocker** — Stops the factor being destroyed, so more of it accumulates.
- **upstream signalling** — Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **chromatin-level** — HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.

### Graded candidates

Combined = geometric mean of TF confidence and drug quality. A weak TF prediction caps the grade no matter how good the drug.

| Grade | Combined | Drug | TF | TF conf. | Drug qual. | Potency | Protein select. | Gene select. | Evidence | Flags |
|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 74.7 | Tamibarotene | RARA | 0.722 | 0.811 | 0.581 | 0.741 | 1.0 | 0.75 | direction not independently predicted |
| **A** | 70.2 | Acitretin | RARA | 0.722 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **A** | 70.2 | Etretinate | RARA | 0.722 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **A** | 70.2 | Isotretinoin | RARA | 0.722 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **A** | 70.2 | Trifarotene | RARA | 0.722 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 67.4 | Adapalene | RARA | 0.722 | 0.575 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 67.4 | Alitretinoin | RARA | 0.722 | 0.575 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 67.4 | Tazarotene | RARA | 0.722 | 0.575 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 67.4 | Tretinoin | RARA | 0.722 | 0.575 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 65.4 | Seladelpar | PPARD | 0.556 | 0.953 | 0.925 | 1.0 | 1.0 | 0.75 | direction not independently predicted |
| **B** | 64.7 | Diethylstilbestrol | ESRRG | 0.646 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 63.9 | Diethylstilbestrol | ESRRB | 0.635 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.6 | Clotrimazole | NR1I3 | 0.617 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.6 | Meclizine | NR1I3 | 0.617 | 0.649 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 62.6 | Elacestrant | ESR1 | 0.517 | 0.832 | 0.58 | 1.0 | 0.926 | 0.75 | structure supports antagonist (p=0.273, moderate confidence) |
| **B** | 62.4 | Fulvestrant | ESR1 | 0.517 | 0.967 | 0.917 | 1.0 | 0.926 | 1.0 | structure supports antagonist (p=0.458, low confidence) |
| **B** | 62.3 | Cholesterol | RORA | 0.608 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 61.7 | Seladelpar Lysine | PPARD | 0.556 | 0.785 | 0.5 | 1.0 | 1.0 | 0.5 | no potency data; scored neutral; thin evidence base (n=0 potency measu |
| **B** | 61.1 | Elacestrant Hydrochloride | ESR1 | 0.517 | 0.77 | 0.5 | 1.0 | 0.926 | 0.5 | no potency data; scored neutral; thin evidence base (n=0 potency measu |
| **B** | 60.6 | Elafibranor | PPARD | 0.556 | 0.739 | 0.4 | 0.741 | 1.0 | 0.75 | direction not independently predicted |
| **B** | 59.1 | Acitretin | RXRG | 0.587 | 0.603 | 0.5 | 0.5 | 0.55 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 59.1 | Etretinate | RXRG | 0.587 | 0.603 | 0.5 | 0.5 | 0.55 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.5 | Bezafibrate | PPARD | 0.556 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.5 | Tretinoin | PPARD | 0.556 | 0.66 | 0.5 | 0.5 | 1.0 | 0.25 | no potency data; scored neutral; barely characterised: no potency or m |
| **B** | 58.5 | Calcitriol | VDR | 0.467 | 0.994 | 0.975 | 1.0 | 1.0 | 1.0 | direction not independently predicted |

### Cell-cycle compatibility

Tissue proliferation state: **unknown**

3 factor(s) here are cell-cycle restricted (E2F6, E2F3, MYC). Set --proliferation to have this scored rather than just flagged.

| TF | Phase | Duty cycle | Effect multiplier |
|---|---|---|---|
| E2F6 | G1/S | 25% | 1.0 |
| E2F3 | G1/S | 25% | 1.0 |
| MYC | G1 | 46% | 1.0 |

*A partial duty cycle dilutes the population-level effect in asynchronous cells: per-cell efficacy can be much higher than a bulk RT-qPCR suggests. Synchronise or use a phase-resolved readout before calling it a failure.*

### Mechanism options: orthosteric vs allosteric

Allosteric sites are not constrained to bind the endogenous ligand, so they diverge across a receptor family where the orthosteric pocket cannot. That usually buys subtype selectivity — at the cost of a bounded maximal effect and a hard dependence on endogenous tone.

**ESRRG**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| inhibitor | 1 | unknown | unknown | no | no | Diethylstilbestrol |

*No approved allosteric modulator at this target; only orthosteric pharmacology, with the usual subtype-selectivity problem.*

**ESRRB**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| inhibitor | 1 | unknown | unknown | no | no | Diethylstilbestrol |

*No approved allosteric modulator at this target; only orthosteric pharmacology, with the usual subtype-selectivity problem.*

**Before pursuing any allosteric option, verify that the pathway has endogenous tone in your cell type.** A PAM with no ligand to modulate does nothing at all, and that is a measurement, not a prediction — nothing in this package can check it for you.

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

**NFIA** — **no ligand pocket** · stabilisation recommended

- NFIA has no ligand-binding pocket and no direct drug. This is exactly the case stabilisation exists for: its protein-protein interfaces and its degradation machinery are both druggable even though its surface is not.
- Nothing curated for this factor specifically. The tractable questions are: which E3 ligase degrades it (block that to raise abundance), and which partner complex sequesters or activates it (glue or disrupt that). Both are established, approved mechanisms — see the exemplars.

  Approved proof the mechanism class works, on other targets:
  - **Tafamidis** (2019) — conformational stabiliser: Binds the transthyretin tetramer interface and prevents dissociation — the cleanest approved example of locking a protein into a stable shape
  - **Migalastat** (2016) — conformational stabiliser: Pharmacological chaperone; stabilises mutant alpha-galactosidase A into a foldable, traffickable form
  - **Sirolimus** (1999) — molecular glue: Glues FKBP12 to mTOR; neither interaction occurs without the drug
  - **Tacrolimus** (1994) — molecular glue: Glues FKBP12 to calcineurin

**CTCF** — **no ligand pocket** · stabilisation recommended

- CTCF has no ligand-binding pocket and no direct drug. This is exactly the case stabilisation exists for: its protein-protein interfaces and its degradation machinery are both druggable even though its surface is not.
- Nothing curated for this factor specifically. The tractable questions are: which E3 ligase degrades it (block that to raise abundance), and which partner complex sequesters or activates it (glue or disrupt that). Both are established, approved mechanisms — see the exemplars.

  Approved proof the mechanism class works, on other targets:
  - **Tafamidis** (2019) — conformational stabiliser: Binds the transthyretin tetramer interface and prevents dissociation — the cleanest approved example of locking a protein into a stable shape
  - **Migalastat** (2016) — conformational stabiliser: Pharmacological chaperone; stabilises mutant alpha-galactosidase A into a foldable, traffickable form
  - **Sirolimus** (1999) — molecular glue: Glues FKBP12 to mTOR; neither interaction occurs without the drug
  - **Tacrolimus** (1994) — molecular glue: Glues FKBP12 to calcineurin

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Cell-type reach

Three things must overlap: where the drug goes, where the gene is expressed, and where the transcription factor is. Any empty intersection kills the route regardless of motif quality or potency.

| TF | Co-expressed with gene | Shared tissues |
|---|---|---|
| RARA | unknown | — |
| ESRRG | unknown | — |
| ESRRB | unknown | — |
| RORB | unknown | — |
| NFIA | unknown | — |
| NR1I3 | unknown | — |
| RORC | unknown | — |
| CTCF | unknown | — |
| RORA | unknown | — |
| PRDM5 | unknown | — |

### Human genetic evidence

Human genetic variation in CCR5 associates with: HIV-1 infection, COVID-19, type 1 diabetes mellitus. Read this as a preview of what modulating the gene may do, including the adverse effects — this is the reasoning that anticipated PCSK9 inhibition years before the trials.

| Phenotype | Genetic score | Overall |
|---|---|---|
| HIV-1 infection | 0.635 | 0.629 |
| COVID-19 | 0.238 | 0.414 |
| type 1 diabetes mellitus | 0.517 | 0.339 |
| Erythema nodosum | 0.551 | 0.335 |
| inflammatory bowel disease | 0.484 | 0.309 |
| ulcerative colitis | 0.416 | 0.266 |
| juvenile idiopathic arthritis | 0.421 | 0.265 |
| celiac disease | 0.361 | 0.222 |

- **RARA → CCR5**: To test this link genetically: look for trans-eQTLs of RARA affecting CCR5 in eQTLGen (blood, n>30,000) or GTEx. A hit is strong causal evidence for the regulatory relationship. Absence is weak evidence against it — trans-eQTL detection needs very large samples and most real regulatory links are undetectable at current power.
- **ESRRG → CCR5**: To test this link genetically: look for trans-eQTLs of ESRRG affecting CCR5 in eQTLGen (blood, n>30,000) or GTEx. A hit is strong causal evidence for the regulatory relationship. Absence is weak evidence against it — trans-eQTL detection needs very large samples and most real regulatory links are undetectable at current power.
- **ESRRB → CCR5**: To test this link genetically: look for trans-eQTLs of ESRRB affecting CCR5 in eQTLGen (blood, n>30,000) or GTEx. A hit is strong causal evidence for the regulatory relationship. Absence is weak evidence against it — trans-eQTL detection needs very large samples and most real regulatory links are undetectable at current power.

### Candidates that may beat the current standard

**The bar: strength 0.55, specificity 0.55, efficiency 0.75** (measured) — the best established value on *each axis*, not one drug's scores. Representative incumbent: Maraviroc (direct).

5 of 100 candidates already have a reported link to this gene (CTD returned 152 chemicals with curated expression effects).

*Accounting for all 100 candidates: 5 already linked, 1 novel and non-dominated, 30 novel but beaten on every axis by another candidate, 64 with no advantage. All of them are in results.json under `improvements.all_candidates`.*

#### Compounds already reported to change this gene

These come from CTD's curated chemical–gene interactions: someone measured an expression change and a curator recorded it. That makes them the most concrete evidence in this report — but the direction is what decides whether each is a lead or a warning.

152 chemicals total: **39** move it the way you want, 89 move it the opposite way, 24 unclear.

**Reported to move CCR5 DOWN** — leads:

| Chemical | Predicted by this run? | CTD actions |
|---|---|---|
| (+)-Jq1 Compound | **no** | decreases^expression |
| 1,3-Butadiene | **no** | decreases^expression |
| 2,2-Dichloro-1-(Triphenylphosphonio)Vinyl Formamide Perchlorate | **no** | decreases^expression |
| 3-Chlorophenol | **no** | decreases^expression |
| A-83-01 | **no** | decreases^expression |
| Ammonium 2,3,3,3-Tetrafluoro-2-(Heptafluoropropoxy)-Propanoate | **no** | decreases^expression |
| Ancriviroc | **no** | decreases^expression |
| Anti-Hiv Agents | **no** | decreases^expression |
| Antirheumatic Agents | **no** | decreases^expression |
| Azoxymethane | **no** | decreases^expression |
| Bisphenol S | **no** | decreases^expression |
| Celecoxib | **no** | decreases^expression |
| Chloroprene | **no** | decreases^expression |
| Cirsimarin | **no** | decreases^expression |
| Dextran Sulfate | **no** | decreases^expression |
| Dinoprostone | **no** | decreases^expression |
| Diuron | **no** | decreases^expression |
| Doxorubicin | **no** | decreases^expression |
| Enniatins | **no** | decreases^expression |
| Ethyl Methanesulfonate | **no** | decreases^expression |
| Fenretinide | **no** | decreases^expression |
| Formaldehyde | **no** | decreases^expression |
| Gsk-J4 | **no** | decreases^expression |
| Hydrolyzable Tannins | **no** | decreases^expression |
| Isopentenyl Pyrophosphate | **no** | decreases^expression |
| … 14 more | | |

**53 of these were NOT predicted by this run at all.** They have a measured effect on your gene and no mechanism from this analysis — which makes them worth reading closely, since an empirical hit the model missed is either a route outside its scope or a clue that the scan is wrong.

(+)-Jq1 Compound, 1,3-Butadiene, 2,2-Dichloro-1-(Triphenylphosphonio)Vinyl Formamide Perchlorate, 2,4,5,2',4',5'-Hexachlorobiphenyl, 3-Chlorophenol, A-83-01, Aflatoxin B1, Ammonium 2,3,3,3-Tetrafluoro-2-(Heptafluoropropoxy)-Propanoate, Ancriviroc, Anti-Hiv Agents, Antirheumatic Agents, Azoxymethane, Benzo(A)Pyrene, Bisphenol A, Bisphenol S, Celecoxib, Chloroprene, Cirsimarin, Dextran Sulfate, Dinoprostone, Diuron, Doxorubicin, Drugs, Chinese Herbal, Enniatins, Ethyl Methanesulfonate …

**Reported to move it the WRONG way** — avoid, or use as negative controls:

2,4-Dichlorophenoxyacetic Acid, 2,4-Dinitrofluorobenzene Sulfonic Acid, 2-(2-Amino-3-Methoxyphenyl)-4H-1-Benzopyran-4-One, 2-Chloroethyl Ethyl Sulfide, 2-Naphthylamine, 4-Hydroxyphenyl 4-Isopropoxyphenylsulfone, Acetamide, Air Pollutants, Amphetamine, Amphotericin B, Deoxycholate Drug Combination, Arsenic, Asbestos, Crocidolite, Atorvastatin, Bardoxolone Methyl, Benzalkonium Compounds, Benzene, Biological Factors, Bisphenol F, Bleomycin, Bromine, Cadmium, Cadmium Chloride, Cannabidiol, Captax, Carbon Tetrachloride …

*CTD aggregates across cell types, doses and assays. A 'both' or contradictory annotation usually means the effect is context-dependent rather than that the data are wrong — check the underlying study before relying on any single row.*

#### Not previously linked, and better on at least one axis

| Candidate | Route | TF | Strength | Specificity | Efficiency | Wins on | Trade-offs |
|---|---|---|---|---|---|---|---|
| Calcitriol | via TF | VDR | 0.676 | 1.0 | 0.99 | higher specificity, superior efficiency | none |

*These are predicted scores compared against measured ones, so a margin was required before counting a win. They are hypotheses to test, not improvements demonstrated.*

Novel but beaten on every axis by another candidate (shown so nothing is hidden): Tamibarotene, Seladelpar, Elacestrant, Fulvestrant, Seladelpar Lysine, Elacestrant Hydrochloride, Elafibranor, Calcipotriene, Cholecalciferol, Trifarotene, Alfacalcidol, Calcifediol ...

#### Non-dominated set (nothing beats these on every axis)

- **Calcitriol** (NEW link) — beats the incumbent on higher specificity, superior efficiency

### Predicted transcriptional outcome

P(gene moves) = P(engages) x P(direction) x P(activates here) x P(responsive). Four uncertain links multiply, so these numbers are small by construction — a chain at 0.7 per link lands near 0.24. Anything reporting 0.9 here has hidden a link.

| Drug | TF | P(engage) | P(direction) | P(promoter) | P(responsive) | **P(gene moves)** | Verdict | Limiting step |
|---|---|---|---|---|---|---|---|---|
| Calcitriol ★ | VDR | 0.95 | 0.75 | 0.477 | 0.4 | **0.136** (0.048–0.224) | unlikely via this route | gene responsiveness in your cells |
| Tamibarotene | RARA | 0.811 | 0.75 | 0.655 | 0.4 | **0.159** (0.08–0.239) | speculative | gene responsiveness in your cells |
| Seladelpar | PPARD | 0.95 | 0.75 | 0.539 | 0.4 | **0.154** (0.054–0.254) | speculative | gene responsiveness in your cells |
| Acitretin | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Etretinate | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Isotretinoin | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Trifarotene | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Adapalene | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Alitretinoin | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Tazarotene | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Tretinoin | RARA | 0.775 | 0.75 | 0.655 | 0.4 | **0.152** (0.076–0.228) | speculative | gene responsiveness in your cells |
| Fulvestrant | ESR1 | 0.95 | 0.75 | 0.512 | 0.4 | **0.146** (0.051–0.241) | unlikely via this route | gene responsiveness in your cells |
| Diethylstilbestrol | ESRRG | 0.775 | 0.75 | 0.602 | 0.4 | **0.14** (0.07–0.21) | unlikely via this route | gene responsiveness in your cells |
| Diethylstilbestrol | ESRRB | 0.775 | 0.75 | 0.595 | 0.4 | **0.138** (0.069–0.207) | unlikely via this route | gene responsiveness in your cells |
| Calcipotriene | VDR | 0.95 | 0.75 | 0.477 | 0.4 | **0.136** (0.048–0.224) | unlikely via this route | gene responsiveness in your cells |
| Clotrimazole | NR1I3 | 0.775 | 0.75 | 0.582 | 0.4 | **0.135** (0.068–0.203) | unlikely via this route | gene responsiveness in your cells |
| Meclizine | NR1I3 | 0.775 | 0.75 | 0.582 | 0.4 | **0.135** (0.068–0.203) | unlikely via this route | gene responsiveness in your cells |
| Cholesterol | RORA | 0.775 | 0.75 | 0.575 | 0.4 | **0.134** (0.067–0.201) | unlikely via this route | gene responsiveness in your cells |
| Cholecalciferol | VDR | 0.929 | 0.75 | 0.477 | 0.4 | **0.133** (0.046–0.219) | unlikely via this route | gene responsiveness in your cells |
| Trifarotene | RARG | 0.933 | 0.75 | 0.469 | 0.4 | **0.131** (0.046–0.217) | unlikely via this route | gene responsiveness in your cells |
| Acitretin | RXRG | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.195) | unlikely via this route | gene responsiveness in your cells |
| Etretinate | RXRG | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.195) | unlikely via this route | gene responsiveness in your cells |
| Alitretinoin | RXRG | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.195) | unlikely via this route | gene responsiveness in your cells |
| Bexarotene | RXRG | 0.775 | 0.75 | 0.561 | 0.4 | **0.13** (0.065–0.195) | unlikely via this route | gene responsiveness in your cells |
| Elacestrant | ESR1 | 0.811 | 0.75 | 0.512 | 0.4 | **0.125** (0.044–0.205) | unlikely via this route | gene responsiveness in your cells |

★ = named in the comparison against the standard of care above, pinned here so the two sections cannot disagree.

**Most informative next experiment** (for Calcitriol): Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.

### Virtual screen of approved drugs against RARA

Binding model, CV AUC 0.997 — trained on 211 measured actives vs 633 decoys, screened 3311 approved drugs.

*Direction not annotated: only 2 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.988 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamibarotene | 0.984 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Trifarotene | 0.956 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.937 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.937 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Adapalene | 0.902 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tazarotene | 0.819 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.937 | 0.0 |  |
| predicted binder, close to known chemistry | Palovarotene | 0.864 | 0.429 |  |
| no call | Acoramidis | 0.563 | 0.681 | outside applicability domain; score withheld |
| no call | Tafamidis Meglumine | 0.559 | 0.727 | outside applicability domain; score withheld |
| no call | Deferasirox | 0.552 | 0.673 | outside applicability domain; score withheld |
| no call | Acoramidis Hydrochloride | 0.547 | 0.686 | outside applicability domain; score withheld |
| no call | Iptacopan | 0.51 | 0.726 | outside applicability domain; score withheld |
| no call | Iptacopan Hydrochloride | 0.491 | 0.733 | outside applicability domain; score withheld |
| no call | Eprosartan | 0.489 | 0.741 | outside applicability domain; score withheld |
| no call | Probenecid | 0.458 | 0.69 | outside applicability domain; score withheld |
| no call | Bentiromide | 0.453 | 0.661 | outside applicability domain; score withheld |
| no call | Eprosartan Mesylate | 0.445 | 0.756 | outside applicability domain; score withheld |
| no call | Eltrombopag | 0.393 | 0.736 | outside applicability domain; score withheld |
| no call | Eltrombopag Olamine | 0.368 | 0.755 | outside applicability domain; score withheld |
| no call | Pranlukast | 0.366 | 0.727 | outside applicability domain; score withheld |
| no call | Novobiocin | 0.353 | 0.76 | outside applicability domain; score withheld |
| no call | Opicapone | 0.351 | 0.733 | outside applicability domain; score withheld |
| no call | Febuxostat | 0.343 | 0.689 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against PPARD

Binding model, CV AUC 0.999 — trained on 1264 measured actives vs 3309 decoys, screened 3311 approved drugs.

*Direction not annotated: only 5 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Seladelpar | 0.988 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Elafibranor | 0.898 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Pemafibrate | 0.798 | 0.0 |  |
| predicted binder, close to known chemistry | Seladelpar Lysine | 0.73 | 0.0 |  |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against ESR1

Binding model, CV AUC 0.998 — trained on 2014 measured actives vs 3294 decoys, screened 3311 approved drugs.

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Raloxifene | 0.995 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Bazedoxifene | 0.987 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Raloxifene Hydrochloride | 0.967 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Lasofoxifene | 0.959 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Ethinyl Estradiol | 0.953 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Diethylstilbestrol | 0.927 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estradiol | 0.919 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Fulvestrant | 0.912 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estrone | 0.896 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estriol | 0.829 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamoxifen | 0.802 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Estetrol Anhydrous | 0.785 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Mestranol | 0.543 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tibolone | 0.49 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Norethindrone | 0.437 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Ethylestrenol | 0.42 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Stanozolol | 0.385 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Hydroquinone | 0.886 | 0.0 | predicted agonist (p=0.635, low confidence) |
| predicted binder, close to known chemistry | Bazedoxifene Acetate | 0.886 | 0.0 |  |
| predicted binder, close to known chemistry | Dienestrol | 0.821 | 0.0 | predicted agonist (p=0.635, low confidence) |
| predicted binder, close to known chemistry | Hexestrol | 0.805 | 0.0 | predicted agonist (p=0.599, low confidence) |
| predicted binder, close to known chemistry | Estetrol | 0.77 | 0.0 | predicted agonist (p=0.741, moderate confidence) |
| predicted binder, close to known chemistry | Mequinol | 0.743 | 0.0 | predicted agonist (p=0.606, low confidence) |
| predicted binder, close to known chemistry | Monobenzone | 0.73 | 0.0 |  |
| predicted binder, close to known chemistry | Lasofoxifene Tartrate | 0.728 | 0.0 |  |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against VDR

Binding model, CV AUC 0.996 — trained on 370 measured actives vs 1110 decoys, screened 3311 approved drugs.

*Direction not annotated: only 18 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Calcipotriene | 0.995 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Calcitriol | 0.991 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Doxercalciferol | 0.989 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tacalcitol | 0.982 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Cholecalciferol | 0.94 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Calcifediol Anhydrous | 0.989 | 0.183 |  |
| predicted binder, close to known chemistry | Calcifediol | 0.984 | 0.197 |  |
| predicted binder, close to known chemistry | Paricalcitol | 0.974 | 0.317 |  |
| predicted binder, close to known chemistry | Ergocalciferol | 0.949 | 0.233 |  |
| predicted binder, close to known chemistry | Alfacalcidol | 0.932 | 0.0 |  |
| predicted binder, close to known chemistry | Dihydrotachysterol | 0.873 | 0.339 |  |
| no call | Omacetaxine Mepesuccinate | 0.214 | 0.792 | outside applicability domain; score withheld |
| no call | Paritaprevir | 0.164 | 0.699 | outside applicability domain; score withheld |
| no call | Pafolacianine | 0.154 | 0.707 | outside applicability domain; score withheld |
| no call | Sildenafil Citrate | 0.149 | 0.714 | outside applicability domain; score withheld |
| no call | Ponesimod | 0.147 | 0.695 | outside applicability domain; score withheld |
| no call | Dalfopristin | 0.146 | 0.829 | outside applicability domain; score withheld |
| no call | Pafolacianine Sodium | 0.146 | 0.699 | outside applicability domain; score withheld |
| no call | Fusidic Acid | 0.145 | 0.753 | outside applicability domain; score withheld |
| no call | Sunitinib Malate | 0.142 | 0.714 | outside applicability domain; score withheld |
| no call | Elagolix | 0.142 | 0.763 | outside applicability domain; score withheld |
| no call | Crinecerfont | 0.141 | 0.8 | outside applicability domain; score withheld |
| no call | Sunitinib | 0.14 | 0.667 | outside applicability domain; score withheld |
| no call | Carbenoxolone | 0.139 | 0.737 | outside applicability domain; score withheld |
| no call | Elagolix Sodium | 0.138 | 0.769 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against THRB

Binding model, CV AUC 1.0 — trained on 33 measured actives vs 99 decoys, screened 3311 approved drugs.

*Direction not annotated: only 6 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

**Screen rejected:** no known actives appear among the top hits — the model is not ranking them above the library, so the novel list cannot be trusted No candidates are reported from it.

### Virtual screen of approved drugs against RARG

Binding model, CV AUC 0.998 — trained on 246 measured actives vs 738 decoys, screened 3311 approved drugs.

*Direction not annotated: only 2 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.983 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamibarotene | 0.978 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.961 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.961 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Adapalene | 0.959 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Trifarotene | 0.95 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tazarotene | 0.804 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.961 | 0.0 |  |
| predicted binder, close to known chemistry | Palovarotene | 0.89 | 0.429 |  |
| no call | Tafamidis Meglumine | 0.571 | 0.727 | outside applicability domain; score withheld |
| no call | Deferasirox | 0.562 | 0.672 | outside applicability domain; score withheld |
| no call | Acoramidis | 0.549 | 0.681 | outside applicability domain; score withheld |
| no call | Acoramidis Hydrochloride | 0.535 | 0.686 | outside applicability domain; score withheld |
| no call | Iptacopan | 0.522 | 0.726 | outside applicability domain; score withheld |
| no call | Eprosartan | 0.521 | 0.735 | outside applicability domain; score withheld |
| no call | Iptacopan Hydrochloride | 0.504 | 0.733 | outside applicability domain; score withheld |
| no call | Eprosartan Mesylate | 0.443 | 0.721 | outside applicability domain; score withheld |
| no call | Probenecid | 0.442 | 0.69 | outside applicability domain; score withheld |
| no call | Febuxostat | 0.415 | 0.689 | outside applicability domain; score withheld |
| no call | Eltrombopag | 0.383 | 0.736 | outside applicability domain; score withheld |
| no call | Bentiromide | 0.375 | 0.661 | outside applicability domain; score withheld |
| no call | Eltrombopag Olamine | 0.351 | 0.755 | outside applicability domain; score withheld |
| no call | Telmisartan | 0.35 | 0.737 | outside applicability domain; score withheld |
| no call | Opicapone | 0.312 | 0.733 | outside applicability domain; score withheld |
| no call | Ledipasvir | 0.294 | 0.663 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against THRA

Binding model, CV AUC 1.0 — trained on 33 measured actives vs 99 decoys, screened 3311 approved drugs.

*Direction not annotated: only 6 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

**Screen rejected:** no known actives appear among the top hits — the model is not ranking them above the library, so the novel list cannot be trusted No candidates are reported from it.

### Virtual screen of approved drugs against RXRA

Binding model, CV AUC 0.998 — trained on 1021 measured actives vs 3063 decoys, screened 3311 approved drugs.

*Direction not annotated: only 1 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.999 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.896 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.896 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Levothyroxine | 0.37 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.896 | 0.0 |  |
| predicted binder, close to known chemistry | Tamibarotene | 0.805 | 0.0 |  |
| predicted binder, close to known chemistry | Oxaprozin | 0.802 | 0.0 |  |
| NOVEL HYPOTHESIS | Palovarotene | 0.859 | 0.516 | extrapolation: nearest known ligand Tanimoto 0.48. The model |
| no call | Celecoxib | 0.19 | 0.684 | outside applicability domain; score withheld |
| no call | Lurbinectedin | 0.174 | 0.778 | outside applicability domain; score withheld |
| no call | Omidenepag Isopropyl | 0.173 | 0.758 | outside applicability domain; score withheld |
| no call | Trabectedin | 0.155 | 0.763 | outside applicability domain; score withheld |
| no call | Pralsetinib | 0.145 | 0.778 | outside applicability domain; score withheld |
| no call | Gemfibrozil | 0.123 | 0.651 | outside applicability domain; score withheld |
| no call | Pemafibrate | 0.118 | 0.683 | outside applicability domain; score withheld |
| no call | Axitinib | 0.093 | 0.78 | outside applicability domain; score withheld |
| no call | Doravirine | 0.093 | 0.761 | outside applicability domain; score withheld |
| no call | Laropiprant | 0.082 | 0.725 | outside applicability domain; score withheld |
| no call | Entrectinib | 0.055 | 0.704 | outside applicability domain; score withheld |
| no call | Oxamniquine | 0.051 | 0.738 | outside applicability domain; score withheld |
| no call | Ceritinib | 0.049 | 0.653 | outside applicability domain; score withheld |
| no call | Ecamsule | 0.049 | 0.709 | outside applicability domain; score withheld |
| no call | Chlorzoxazone | 0.03 | 0.729 | outside applicability domain; score withheld |
| no call | Molsidomine | 0.03 | 0.737 | outside applicability domain; score withheld |
| no call | Lorlatinib | 0.029 | 0.717 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Virtual screen of approved drugs against RARB

Binding model, CV AUC 0.998 — trained on 272 measured actives vs 816 decoys, screened 3311 approved drugs.

*Direction not annotated: only 1 drugs at this receptor carry an agonist/antagonist label. The hits below are predicted BINDERS — whether each would activate or inhibit is unknown and needs a functional assay.*

| Quadrant | Drug | P(wanted direction) | Novelty | Note |
|---|---|---|---|---|
| ESTABLISHED LIGAND, untested on this gene | Bexarotene | 0.99 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tamibarotene | 0.982 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Adapalene | 0.96 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Trifarotene | 0.952 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tretinoin | 0.936 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Alitretinoin | 0.936 | 0.0 | binds this receptor with measured activity, and has no repor |
| ESTABLISHED LIGAND, untested on this gene | Tazarotene | 0.849 | 0.0 | binds this receptor with measured activity, and has no repor |
| predicted binder, close to known chemistry | Isotretinoin | 0.936 | 0.0 |  |
| predicted binder, close to known chemistry | Palovarotene | 0.897 | 0.429 |  |
| no call | Acoramidis | 0.597 | 0.672 | outside applicability domain; score withheld |
| no call | Acoramidis Hydrochloride | 0.574 | 0.677 | outside applicability domain; score withheld |
| no call | Iptacopan | 0.507 | 0.716 | outside applicability domain; score withheld |
| no call | Deferasirox | 0.506 | 0.673 | outside applicability domain; score withheld |
| no call | Iptacopan Hydrochloride | 0.481 | 0.707 | outside applicability domain; score withheld |
| no call | Eprosartan Mesylate | 0.466 | 0.672 | outside applicability domain; score withheld |
| no call | Bentiromide | 0.415 | 0.661 | outside applicability domain; score withheld |
| no call | Telmisartan | 0.409 | 0.685 | outside applicability domain; score withheld |
| no call | Eltrombopag | 0.379 | 0.736 | outside applicability domain; score withheld |
| no call | Febuxostat | 0.376 | 0.689 | outside applicability domain; score withheld |
| no call | Pranlukast | 0.338 | 0.727 | outside applicability domain; score withheld |
| no call | Eltrombopag Olamine | 0.335 | 0.753 | outside applicability domain; score withheld |
| no call | Ledipasvir | 0.329 | 0.663 | outside applicability domain; score withheld |
| no call | Lusutrombopag | 0.293 | 0.762 | outside applicability domain; score withheld |
| no call | Opicapone | 0.288 | 0.733 | outside applicability domain; score withheld |
| no call | Alogliptin Benzoate | 0.277 | 0.75 | outside applicability domain; score withheld |

*Two independent axes here. **ESTABLISHED LIGAND, untested on this gene** is the most actionable class: the drug reliably binds this receptor and has no reported effect on your gene, so the pharmacology is already de-risked and only the receptor-to-gene step is a hypothesis. NOVEL HYPOTHESIS rows are hypotheses on both counts — novel chemistry AND an untested gene link — so expect most to be wrong and treat each as requiring a binding assay first.*

### Drug–gene links with no recorded direction

829 interaction(s) across 73 factor(s) come from curated sources that record a drug–gene relationship **without saying which way it goes**. They are held back from the recommendations above, because the entire question here is direction and an undirected link cannot answer it.

They are still leads — someone curated a real relationship — but they need a literature check, not a prediction. Many are chemotherapeutics or broad stressors whose association with a transcription factor reflects cell-wide perturbation rather than a specific mechanism.

| TF | Drug | Reported by | DGIdb score |
|---|---|---|---|
| EOMES | Glatiramer | PharmGKB | 15.12 |
| ELK1 | Ganciclovir | NCI | 13.23 |
| SMAD4 | Lysine | CIViC | 10.59 |
| NR1I2 | Rifaximin | DTC | 4.41 |
| TFAP2A | Amprenavir | NCI | 4.41 |
| CUX2 | Rifampin | PharmGKB | 3.65 |
| MAFK | Rifampin | PharmGKB | 3.65 |
| STAT5A | Tofacitinib | DTC | 3.61 |
| GATA4 | Acenocoumarol | PharmGKB | 3.53 |
| NFATC4 | Anhydrous Tacrolimus | PharmGKB | 3.41 |
| SOX10 | Vemurafenib | CIViC | 3.11 |
| TBXT | Flunisolide Anhydrous | PharmGKB | 2.65 |
| FOXC2 | Oxaliplatin | CIViC | 2.58 |
| NEUROD1 | Repaglinide | PharmGKB | 2.52 |
| NEUROD1 | Deferoxamine | NCI | 2.35 |
| ETV1 | Trametinib Dimethyl Sulfoxide | CIViC | 2.3 |
| ETV5 | Trametinib Dimethyl Sulfoxide | CIViC | 2.3 |
| KLF4 | Hydroxyurea | PharmGKB | 2.3 |
| TCF7L2 | Linagliptin | PharmGKB | 2.21 |
| NR1I3 | Oxymorphone | DTC | 2.08 |
| HNF1A | Rosuvastatin | PharmGKB | 2.04 |
| SREBF1 | Insulin, Regular, Human | NCI | 2.04 |
| STAT5B | Tofacitinib | DTC | 1.92 |
| PLAGL1 | Olanzapine | PharmGKB | 1.56 |
| STAT6 | Esomeprazole Sodium | PharmGKB | 1.56 |
| … 804 more | | | |

*Pass `--include-unknown-direction` to fold these into the main tables anyway — but every such row will read UNKNOWN / unclear, which is what they are.*

### Approved drugs acting on these factors

| TF | TF role | Do what to this TF | Drug | Action | Drug direction | First approval | Source |
|---|---|---|---|---|---|---|---|
| RARA | repressor | increase | Acitretin | AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Adapalene | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Alitretinoin | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Etretinate | AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Isotretinoin | AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Tamibarotene | AGONIST | up | — | ChEMBL |
| RARA | repressor | increase | Tazarotene | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Tretinoin | AGONIST/AGONIST | up | — | DGIdb |
| RARA | repressor | increase | Trifarotene | AGONIST | up | — | DGIdb |
| ESRRG | activator | decrease | Diethylstilbestrol | INHIBITOR | down | — | DGIdb |
| ESRRB | activator | decrease | Diethylstilbestrol | INHIBITOR | down | — | DGIdb |
| NR1I3 | activator | decrease | Clotrimazole | INHIBITOR | down | — | DGIdb |
| NR1I3 | activator | decrease | Meclizine | INHIBITOR | down | — | DGIdb |
| RORA | unknown | either | Cholesterol | AGONIST | up | — | DGIdb |
| RXRG | unknown | either | Acitretin | AGONIST | up | — | DGIdb |
| RXRG | unknown | either | Alitretinoin | AGONIST/AGONIST | up | — | DGIdb |
| RXRG | unknown | either | Bexarotene | AGONIST/AGONIST | up | — | DGIdb |
| RXRG | unknown | either | Etretinate | AGONIST | up | — | DGIdb |
| PPARD | repressor | increase | Bezafibrate | AGONIST | up | — | DGIdb |
| PPARD | repressor | increase | Elafibranor | AGONIST | up | 2024 | ChEMBL |
| PPARD | repressor | increase | Seladelpar | AGONIST | up | 2024 | ChEMBL |
| PPARD | repressor | increase | Seladelpar Lysine | AGONIST | up | 2024 | ChEMBL |
| PPARD | repressor | increase | Tretinoin | AGONIST | up | — | DGIdb |
| PPARG | dual | decrease | Bexarotene | INHIBITOR | down | — | DGIdb |
| PPARG | dual | decrease | Diclofenac Sodium | INHIBITOR | down | — | DGIdb |
| … 74 more, see tf_summary.tsv | | | | | | | |

## Screen design

*One plan for the gene. The reporter insert below uses the ENST00000292303 window; clone a second construct for the other TSS if you want to compare them.*

The outcome chain above is almost always limited by *gene responsiveness in your cell type* — the link this package cannot see. A screen measures it instead of assuming it. Everything above prioritises this experiment; it does not replace it.

### Reporter construct

- Insert: 3:46365946-46371945 (5000 bp upstream / 1000 bp downstream of the ENST00000292303 TSS)
- Vector: promoterless luciferase (pGL4.10-style) or GFP for FACS-based screening; luciferase for arrayed compound screens, GFP for pooled CRISPRa
- Controls:
  - Promoterless vector: establishes the floor. Without it you cannot tell activation from vector background.
  - Constitutive promoter (CMV/EF1a): establishes the ceiling and reports transfection or transduction efficiency per well.
  - Scrambled-motif variant of your top predicted element: the single most informative control here. If a drug induces the wild-type reporter but not the mutant, the effect runs through that element. If it induces both, the mechanism is elsewhere and the prediction was wrong.
  - Dual-reporter normalisation (Renilla or similar) to separate real induction from changes in cell number or general translation — cytotoxic compounds otherwise read as repressors.
- An episomal reporter lacks the chromatin context of the endogenous locus. A hit should be confirmed at the native gene by RT-qPCR before it is believed.

### CRISPRa screen

- Approach: dCas9-VP64 with SAM (MS2-p65-HSF1) or dCas9-VPR for stronger activation
- Focused library: 20 TFs prioritised by this run (NHLH1, ESRRG, GATA2, ETV2, MAFG, OSR2, NR1I3, PPARD, PBX1, ESRRB...)
- Unbiased: ~1,600 human TFs (Calabrese or a TF-focused sublibrary). Run this if the focused set comes back empty — an unbiased screen is the only thing that finds a factor nobody has associated with your promoter.
- Coverage: 3-5 sgRNAs per gene, >=500 cells per guide at transduction, MOI ~0.3 to keep single integrations
- Readout: FACS sort GFP-high vs GFP-low, then sequence guide abundance; or arrayed luciferase for a focused set

Positive-control sgRNAs at your promoter (7 clean of 12 shown):

| Protospacer (20 nt) | PAM | Strand | vs TSS | GC | Note |
|---|---|---|---|---|---|
| `CTAACAGATTCTGTGTAGTG` | GG | + | -253 | 40% | ok |
| `GAATAGATCTCTGGTCTGAA` | GG | - | -177 | 40% | ok |
| `GGTCTGAAGGTTTATTTACG` | GG | - | -189 | 40% | ok |
| `TCTAACAGATTCTGTGTAGT` | GG | + | -254 | 35% | ok |
| `AAGCTAGAGAATAGATCTCT` | GG | - | -169 | 35% | ok |
| `TGGTCTGAAGGTTTATTTAC` | GG | - | -188 | 35% | ok |
| `TGATTTGTACAAGATCACAG` | GG | - | -374 | 35% | ok |
| `TGCTGTTTCTTTTGAAGGAG` | GG | - | -71 | 40% | contains TTTT: terminates U6 transcription, t |
| `TGTTTCTTTTGAAGGAGGGT` | GG | - | -74 | 40% | contains TTTT: terminates U6 transcription, t |
| `ATTTACGGGCTTTTCTCACT` | GG | - | -202 | 40% | contains TTTT: terminates U6 transcription, t |

*Guides listed target YOUR promoter directly, for a positive-control arm proving the reporter can be activated at all. The library guides target the TFs. Run these through CRISPOR or CHOPCHOP for off-target scoring before ordering — no off-target analysis is done here.*

### Compound screen

- Prestwick Chemical Library (~1,520 approved drugs, high chemical and pharmacological diversity)
- Selleck FDA-approved Drug Library (~3,000)
- NCATS Pharmaceutical Collection (NPC) — free, and the most complete approved set

- Design: Arrayed, 384-well, 8-point dose-response on the actives. Single-concentration primary screens at 10 uM produce a cytotoxicity-dominated hit list; dose-response separates real induction from cells dying.
- Counter-screen: Run the same plate against a constitutive-promoter reporter. Anything that moves both is affecting transcription or translation globally, not your promoter, and that is the majority of primary hits.

Cherry-pick, in priority order:

| Compound | Rationale |
|---|---|
| Maraviroc | direct gene-targeted agent |
| Seladelpar | PPARD route, grade A |
| Seladelpar Lysine | PPARD route, grade B |
| Diethylstilbestrol | ESRRG route, grade B |
| Elafibranor | PPARD route, grade B |
| Tamibarotene | RARA route, grade B |
| Elacestrant | ESR1 route, grade B |
| Fulvestrant | ESR1 route, grade B |
| Bezafibrate | PPARD route, grade B |
| Tretinoin | PPARD route, grade B |
| Clotrimazole | NR1I3 route, grade B |
| Meclizine | NR1I3 route, grade B |
| Diethylstilbestrol | ESRRB route, grade B |
| Elacestrant Hydrochloride | ESR1 route, grade B |
| Mifepristone | PGR route, grade B |
| Acitretin | RARA route, grade B |
| Etretinate | RARA route, grade B |
| Isotretinoin | RARA route, grade B |
| Trifarotene | RARA route, grade B |
| Bexarotene | PPARG route, grade B |
| Diclofenac Sodium | PPARG route, grade B |
| Resveratrol | PPARG route, grade B |
| Calcitriol | VDR route, grade B |
| Adapalene | RARA route, grade B |

*A full-library screen is unbiased and will find mechanisms this package cannot predict — that is its advantage, not a criticism of it. Use the cherry-pick list to get an answer in a week; run the full library when you want the answer you did not expect.*

### What to expect

- **compound screen hit rate**: 0.5-2% at a sensible threshold. Most primary hits are cytotoxic or assay-interfering; expect to lose 60-80% at counter-screen.
- **crispra hit rate**: typically a handful of TFs from ~1,600. A screen returning hundreds of hits is reporting a global effect on the reporter, not specific regulation.
- **replicates**: 3 biological replicates minimum; 2 gives you no estimate of variance and a Z' you cannot trust
- **quality metric**: Z' >= 0.5 between positive and negative controls before you screen anything. Below that the screen cannot distinguish a real hit from noise and running it wastes the library.
- **what a negative result means**: If no TF in a genome-wide CRISPRa screen activates your reporter, the promoter fragment is probably not sufficient — the regulation likely lives in an enhancer outside your window. Widen it (--upstream 50000) or use the endogenous locus.

---

### How to read this

**For downregulation, check the direct modalities first.** Approved siRNA (patisiran, givosiran, inclisiran, vutrisiran) and antisense oligonucleotides (inotersen, nusinersen, tofersen) knock down a named transcript with sequence-level specificity. A TF-level intervention hits every other gene that factor controls. Use this report when a direct modality is not an option, or to understand *why* a drug you already have lowers your gene.

A PWM hit means the sequence *resembles* a binding site. It is a hypothesis generator, not evidence of binding, and not evidence that binding would increase transcription of this gene in your cell type. Rank order here is a heuristic (weights printed in the JSON output), not a validated predictor. Every row needs wet-lab confirmation: promoter-reporter assay for the element, agonist dose-response with RT-qPCR/protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it.