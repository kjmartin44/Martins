# Upregulation strategy for DEMOGENE — promoter / TF / drug analysis

Objective: **increase** transcription of DEMOGENE.

> **FDA filter not applied:** not applied to demo data, whose drugs are synthetic

> **DEMO MODE.** Synthetic sequence, illustrative matrices and a hand-entered drug table. Not usable for real conclusions.

## Evidence available in this run

**⚠ DEGRADED** — Read the ranking as a reading order, not as evidence. At least one pillar of the analysis is missing.

- no source of measured evidence is available (ChIP, LINCS, genetics all absent or failed)

| Source | State | What it returned |
|---|---|---|
| chip | not requested | no ChIP source given |
| regulon | not requested | no regulon file; coarse bundled tiers in use |

Not requested (enable to strengthen the run): chip (`--chip-bed (preferred) or --chip-api`), regulon (`--regulon-file`).

## Drugs targeting DEMOGENE directly

None found in the curated table or in ChEMBL/DGIdb for this direction. The curated list is hand-maintained and goes out of date; check FDA/EMA and clinicaltrials.gov before concluding nothing exists.

## Model-predicted medications for DEMOGENE

Compounds this analysis predicts would **raise** DEMOGENE, with the reasoning for each. Everything here is a **hypothesis generated from annotation**, not a recommendation — the tier on each row tells you how far the claim is reaching.

> **⚠ Read as a reading order, not as evidence.** The run's evidence verdict is *degraded*: at least one pillar of the analysis is missing. See the evidence table above.

| # | Medication | What it is | Approved to treat | Acts on | BBB? | Tier | Grade | Score |
|---|---|---|---|---|---|---|---|---|
| 1 | **Calcitriol** | active vitamin D, VDR agonist | hypocalcaemia and secondary hyperparathyroidism in renal failure | VDR | no | 2 | A | 67.1 |
| 2 | **Paricalcitol** | selective VDR activator | secondary hyperparathyroidism in chronic kidney disease | VDR | — | 2 | A | 64.9 |
| 3 | **Dexamethasone** | potent glucocorticoid receptor agonist | inflammation, cerebral oedema, myeloma, and as an antiemetic | NR3C1 | no | 2 | B | 51.2 |
| 4 | **Rosiglitazone** | PPARgamma agonist (thiazolidinedione) | type 2 diabetes; restricted for cardiovascular risk | PPARG | — | 2 | B | 51.0 |
| 5 | **Tretinoin** | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically | RARA | — | 2 | B | 49.7 |
| 6 | **Bexarotene** | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma | RXRA | — | 2 | B | 48.9 |
| 7 | **Prednisolone** | glucocorticoid | asthma, autoimmune and inflammatory disease | NR3C1 | — | 2 | B | 47.8 |
| 8 | **Alitretinoin** | RAR/RXR pan-agonist | Kaposi sarcoma skin lesions; chronic hand eczema | RARA | — | 2 | B | 47.4 |
| 9 | **Hydrocortisone** | glucocorticoid and mineralocorticoid | adrenal insufficiency; topical dermatitis | NR3C1 | — | 2 | B | 44.3 |
| 10 | **Pioglitazone** | PPARgamma agonist | type 2 diabetes | PPARG | — | 2 | C | 43.0 |

**Every row is off-label.** The drug was approved for what the fourth column says, and is listed here for something else. That column is the baseline this proposal departs from, and the distance between it and your gene is part of what you are judging - a topical acne gel and an infused cytotoxic can sit on adjacent rows with similar scores.

*Ordered by evidence tier first, then score. A tier-1 route with a mediocre score beats a tier-5 route with a great one. **BBB?** is not part of the ordering — it is an applicability filter you apply yourself, since whether brain exposure is wanted or unwanted depends on where your gene needs to move. See the note under the graded table for what a yes does and does not establish.*

**Engage** and **promoter** below are the two links that come from evidence: the drug engaging its target at a dose that can be reached, and the factor acting on this promoter. They differ between drugs and are 0 to 1 heuristics, not probabilities. Order on tier and score.

### Why each one

Every route below is **evidence tier 2** — measured experimentally.

**The direction rule**, which applies to every group below: you want to raise the gene, so you agonise an activator and block a repressor — and a dual-role factor settles neither way, because it activates at some promoters and represses at others.

**One experiment falsifies all of them:** Knock the TF down (siRNA/CRISPRi) in your cell type and measure the target gene. If it does not move, no drug acting through that factor will help, and everything upstream is moot.

#### Through **VDR** — 2 candidates

- Role: **activator**
- Binding: 1 predicted site(s), the best scoring 0.785 -439 bp from the TSS; no ChIP data in this run, so the site is a sequence prediction only.
- Collateral: the factor controls roughly 300 gene(s), which all move with yours; clinical burden **moderate** (vitamin D receptor agonist).
- Family: VDR shares its structurally diverse lipophiles pocket with NR1I2, NR1I3, so a ligand here will tend to engage them too — though subtype selectivity has been achieved (calcitriol (VDR, does not meaningfully engage PXR/CAR)).

1. **Calcitriol** — grade A, engage 0.95, promoter 0.59
   - ⚠ measured: pChEMBL 9.3 against VDR (n=6), with a recorded action
2. **Paricalcitol** — grade A, engage 0.95, promoter 0.59
   - ⚠ measured: pChEMBL 8.6 against VDR (n=6), with a recorded action

#### Through **NR3C1** — 3 candidates

- Role: **dual-role**
- Binding: 1 predicted site(s), the best scoring 0.972 -850 bp from the TSS; no ChIP data in this run, so the site is a sequence prediction only.
- Collateral: the factor controls roughly 1200 gene(s), which all move with yours; clinical burden **high** (glucocorticoid).
- Family: NR3C1 shares its 3-ketosteroids (cortisol, aldosterone, progesterone, DHT) pocket with AR, NR3C2, PGR, so a ligand here will tend to engage them too — though subtype selectivity has been achieved (eplerenone (MR-selective vs spironolactone)).

- ⚠ NR3C1 controls ~1200 genes; expect broad transcriptional effects

3. **Dexamethasone** — grade B, engage 0.89, promoter 0.57
   - ⚠ measured: pChEMBL 8.0 against NR3C1 (n=6), with a recorded action
4. **Prednisolone** — grade B, engage 0.82, promoter 0.57
   - ⚠ measured: pChEMBL 7.4 against NR3C1 (n=6), with a recorded action
5. **Hydrocortisone** — grade B, engage 0.79, promoter 0.57
   - ⚠ measured: pChEMBL 7.1 against NR3C1 (n=6), with a recorded action

#### Through **PPARG** — 2 candidates

- Role: **activator**
- Binding: 2 predicted site(s), the best scoring 1.0 -610 bp from the TSS; no ChIP data in this run, so the site is a sequence prediction only.
- Collateral: the factor controls roughly 1200 gene(s), which all move with yours; clinical burden **high** (PPAR agonist).
- Family: PPARG shares its fatty acids and eicosanoids pocket with PPARA, PPARD, so a ligand here will tend to engage them too — though subtype selectivity has been achieved (seladelpar (PPARD)).

- ⚠ PPARG controls ~1200 genes; expect broad transcriptional effects
- ⚠ black-box warning

6. **Rosiglitazone** — grade B, engage 0.86, promoter 0.65
   - ⚠ measured: pChEMBL 7.8 against PPARG (n=6), with a recorded action
7. **Pioglitazone** — grade C, engage 0.74, promoter 0.65
   - ⚠ measured: pChEMBL 6.7 against PPARG (n=6), with a recorded action

#### Through **RARA** — 2 candidates

- Role: **activator**
- Binding: 3 predicted site(s), the best scoring 0.956 -430 bp from the TSS; no ChIP data in this run, so the site is a sequence prediction only.
- Collateral: the factor controls roughly 1200 gene(s), which all move with yours; clinical burden **high** (retinoid).
- Family: RARA shares its all-trans retinoic acid pocket with RARB, RARG, so a ligand here will tend to engage them too — though subtype selectivity has been achieved (trifarotene (RARG-selective)).

- ⚠ promiscuous: 6 annotated targets

8. **Tretinoin** — grade B, engage 0.92, promoter 0.68
   - ⚠ measured: pChEMBL 8.3 against RARA (n=6), with a recorded action
9. **Alitretinoin** — grade B, engage 0.88, promoter 0.68
   - ⚠ measured: pChEMBL 7.9 against RARA (n=6), with a recorded action

#### Through **RXRA**

- Role: **activator**
- Binding: 6 predicted site(s), the best scoring 1.0 -610 bp from the TSS; no ChIP data in this run, so the site is a sequence prediction only.
- Collateral: the factor controls roughly 1200 gene(s), which all move with yours; clinical burden **high** (retinoid).
- Family: RXRA shares its 9-cis retinoic acid pocket with RXRB, RXRG, so a ligand here will tend to engage them too, and no subtype-selective agent is recorded for this family.

10. **Bexarotene** — grade B, engage 0.85, promoter 0.68
   - ⚠ measured: pChEMBL 7.7 against RXRA (n=6), with a recorded action

*Not medical advice. "Approved" means approved for its own indication — using one of these to move an unrelated gene is off-label repurposing with the full systemic pharmacology attached. The experimental sections at the end of this report are how you would actually test any of it.*

## Promoter of DEMOGENE (ENST00000000000)

- Transcript chosen: **demo fixture**
- Assembly: DEMO  |  Strand: +
- TSS: chr1:1000000
- Window: chr1:999000-1000199 (-1000 / +200, 1200 bp)
- Motifs scanned: 10 at p ≤ 0.0001 per position
- Non-redundant sites: 14 across 10 TFs

*Scan quality: 14 sites against 2 expected by chance over 1200 scannable bases (7.1x).*

*about 14% of these 14 sites (range 0% to 29%) would be expected on shuffled sequence of this window, so roughly 12 are above chance; among the best-scoring sites the share is lower: 0% of the top 5, 1% of the top 10. Estimated from 30 shuffled copies that keep the repeat mask and every dinucleotide count. The formula that assumes each motif matches at exactly its stated p-value would expect 2 and is not used. Nothing here says which sites are real; only ChIP evidence separates them.*

### Ligand-activated (nuclear) receptors with sites here

These are the directly druggable ones.

| TF | Role | Sites | Best rel. score | Best site (vs TSS) | ChIP peaks | Score |
|---|---|---|---|---|---|---|
| RARA | activator | 3 | 0.956 | AGGTCAGGCTTAGGTCA (-430) | 0 | 75.5 |
| RXRA | activator | 6 | 1.0 | AGGTCAAAGGTCA (-610) | 0 | 75.3 |
| PPARG | activator | 2 | 1.0 | AGGTCAAAGGTCA (-610) | 0 | 71.9 |
| VDR | activator | 1 | 0.785 | ACTTCTTAAAGGTCA (-439) | 0 | 62.4 |
| NR3C1 | dual | 1 | 0.972 | AGAACAGTGTGTTCT (-850) | 0 | 60.1 |
| HNF4A | activator | 1 | 0.874 | CAGGTCAAAGGTCA (-611) | 0 | 53.9 |

#### Why the drug tables below look the way they do

**10 factor(s) had predicted sites here. 6 of them have an approved drug (5 of those are nuclear receptors).**

Roughly 48 of the ~1,600 human transcription factors have a ligand-binding pocket, and the steroid, retinoid, thyroid and vitamin-D families are most of them. So the drug tables are not a ranking that hormone drugs won — they are what is left after every factor without a ligand has been dropped.

**4 of the top 25 factors by score have no drug and appear in no drug table:** SP1, CREB1, NFKB1, HNF4A.

That is where the gene-specific signal is. A half-site like AGGTCA occurs every few kilobases by chance, so nuclear receptors surface for almost any promoter; these factors surfaced for *this* one. The routes that can reach them are the upstream-network, stabilisation and protein-level sections, not the ligand tables.

### Top-ranked transcription factors

| Rank | TF | Role | NR | Other genes it controls | Sites | Best rel. score | Offset | ChIP | Approved drugs | Score |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | RARA | activator | yes | ~1200 | 3 | 0.956 | -430 | 0 | 2 | 75.5 |
| 2 | RXRA | activator | yes | ~1200 | 6 | 1.0 | -610 | 0 | 1 | 75.3 |
| 3 | PPARG | activator | yes | ~1200 | 2 | 1.0 | -610 | 0 | 3 | 71.9 |
| 4 | SP1 | activator |  | ~5000 | 2 | 1.0 | -95 | 0 | 0 | 70.5 |
| 5 | CREB1 | activator |  | ~5000 | 1 | 1.0 | -140 | 0 | 0 | 65.9 |
| 6 | VDR | activator | yes | ~300 | 1 | 0.785 | -439 | 0 | 3 | 62.4 |
| 7 | REST | repressor |  | ~300 | 2 | 1.0 | -720 | 0 | 1 | 60.6 |
| 8 | NR3C1 | dual | yes | ~1200 | 1 | 0.972 | -850 | 0 | 3 | 60.1 |
| 9 | NFKB1 | dual |  | ~5000 | 1 | 1.0 | -310 | 0 | 0 | 56.2 |
| 10 | HNF4A | activator | yes | ~300 | 1 | 0.874 | -611 | 0 | 0 | 53.9 |

**Other genes it controls** is the collateral question. A drug acting through a factor moves everything that factor regulates, so this column says how wide the effect will be. Read it against the score: a strong route through a factor that controls two thousand genes is a worse starting point than a weaker route through one that controls thirty.

A count written `~N` is a bundled estimate for a factor the regulon file does not list (10 here), not a count of curated edges. Treat it as an order of magnitude.

Counts come from literature curation, so they measure what has been published as much as what the factor does. SP1 has thousands of recorded targets partly because it regulates many genes and partly because thousands of people assayed it.

### Every route that acts on this gene

31 routes across 5 mechanism categories, pooled from every method in this run and sorted by **evidence tier first**, score second.

A tier-1 route with a mediocre score beats a tier-5 route with a great one. These are not the same kind of claim and averaging them would destroy the information you need.

**Why established binders sit at tier 2, just under the direct agents.** A drug already measured to bind one of your factors has half its chain settled: the drug-to-receptor step has assay data behind it, so the pharmacology is de-risked and only the receptor-to-gene step is a hypothesis. A curated tier-3 ligand is annotation on *both* halves. But they stay below tier 1, because a direct agent acts on your gene itself and skips the receptor-to-gene inference entirely — a stronger claim however well-characterised the binder is. The cost of the indirect route is in the row: everything else that factor regulates moves with your gene.

- 10 × measured — the drug's target engagement is observed (LINCS, ChIP, or a measured receptor binder)
- 18 × curated — mechanism annotated by a human curator
- 3 × hypothesis — generated, untested

| Tier | Route | Drug | Acts on | Steps | Mechanism | Str | Spec | Eff | Score |
|---|---|---|---|---|---|---|---|---|---|
| 2 | TF: orthosteric ligand | Calcitriol | VDR | 1 | ligand at the transcription factor | 0.906 | 0.673 | 0.963 | 67.1 |
| 2 | TF: orthosteric ligand | Paricalcitol | VDR | 1 | ligand at the transcription factor | 0.861 | 0.673 | 0.923 | 64.9 |
| 2 | TF: orthosteric ligand | Dexamethasone | NR3C1 | 1 | ligand at the transcription factor | 0.788 | 0.311 | 0.863 | 51.2 |
| 2 | TF: orthosteric ligand | Rosiglitazone | PPARG | 1 | ligand at the transcription factor | 0.795 | 0.311 | 0.843 | 51.0 |
| 2 | TF: orthosteric ligand | Tretinoin | RARA | 1 | ligand at the transcription factor | 0.86 | 0.153 | 0.893 | 49.7 |
| 2 | TF: orthosteric ligand | Bexarotene | RXRA | 1 | ligand at the transcription factor | 0.792 | 0.247 | 0.833 | 48.9 |
| 2 | TF: orthosteric ligand | Prednisolone | NR3C1 | 1 | ligand at the transcription factor | 0.72 | 0.311 | 0.802 | 47.8 |
| 2 | TF: orthosteric ligand | Alitretinoin | RARA | 1 | ligand at the transcription factor | 0.815 | 0.153 | 0.853 | 47.4 |
| 2 | TF: orthosteric ligand | Hydrocortisone | NR3C1 | 1 | ligand at the transcription factor | 0.686 | 0.247 | 0.772 | 44.3 |
| 2 | TF: orthosteric ligand | Pioglitazone | PPARG | 1 | ligand at the transcription factor | 0.671 | 0.247 | 0.733 | 43.0 |
| 3 | TF: orthosteric ligand | Demo-Pam-2 | VDR | 1 | ligand at the transcription factor | 0.576 | 0.437 | 0.688 | 44.4 |
| 3 | TF: endogenous ligand | Carbenoxolone | HSD11B2 | 2 | inhibit catabolism (cortisol) via HSD11B2 | 0.7 | 0.6 | 0.65 | 43.5 |
| 3 | TF: allosteric modulator | Demo-Pam-1 | PPARG | 1 | PAM | 0.6 | 0.273 | 0.688 | 40.6 |
| 3 | TF: endogenous ligand | Tretinoin | — | 2 | supplement ligand (all-trans retinoic acid) v | 0.7 | 0.3 | 0.65 | 36.5 |
| 3 | TF: endogenous ligand | Calcifediol | — | 2 | supplement ligand (calcitriol) via — | 0.7 | 0.3 | 0.65 | 36.5 |
| 3 | TF: endogenous ligand | Ketoconazole | CYP24A1 | 2 | inhibit catabolism (calcitriol) via CYP24A1 | 0.7 | 0.3 | 0.65 | 36.5 |
| 3 | TF: endogenous ligand | Alfacalcidol | — | 2 | supplement ligand (calcitriol) via — | 0.7 | 0.3 | 0.65 | 36.5 |
| 3 | TF: endogenous ligand | Cholecalciferol | — | 2 | supplement ligand (calcitriol) via — | 0.7 | 0.3 | 0.65 | 36.5 |
| 3 | TF: endogenous ligand | Hydrocortisone | — | 2 | supplement ligand (cortisol) via — | 0.7 | 0.3 | 0.65 | 36.5 |
| 3 | TF: orthosteric ligand | Valproic Acid | REST | 1 | ligand at the transcription factor | 0.421 | 0.195 | 0.562 | 30.2 |
| 3 | chromatin-level | Vorinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Romidepsin | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Belinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Panobinostat | chromatin (genome-wide) | 3 | HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Valproic acid | chromatin (genome-wide) | 3 | weak HDAC inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Azacitidine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Decitabine | chromatin (genome-wide) | 3 | DNMT inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 3 | chromatin-level | Tazemetostat | chromatin (genome-wide) | 3 | EZH2 inhibitor | 0.7 | 0.05 | 0.45 | 23.4 |
| 5 | TF: stabiliser or glue | MRL-24 | PPARG | 1 | conformational stabiliser | 0.6 | 0.65 | 0.4 | 45.4 |
| 5 | TF: endogenous ligand | Talarozole | CYP26A1/B1 | 2 | inhibit catabolism (all-trans retinoic acid)  | 0.5 | 0.6 | 0.35 | 33.2 |

Best in each category:

- **TF: orthosteric ligand** — Calcitriol (tier 2, score 67.1)
- **TF: endogenous ligand** — Carbenoxolone (tier 3, score 43.5)
- **TF: allosteric modulator** — Demo-Pam-1 (tier 3, score 40.6)
- **chromatin-level** — Vorinostat (tier 3, score 23.4)
- **TF: stabiliser or glue** — MRL-24 (tier 5, score 45.4)

*No route found in: direct: gene or transcript, direct: gene product, direct: predicted product binder, promoter DNA binder, TF: inferred activity change, TF: abundance change, TF: predicted pocket binder, TF: degradation blocker, TF: known regulator, upstream signalling, off-pathway pocket binder, measured phenocopy, empirical signature, virtual screen hypothesis. An absent category is usually a real gap in the pharmacology, but check you enabled the relevant flag (--lincs, --screen, --upstream-network).*

*Steps = how far the drug acts from your gene: 0 on the gene or transcript, 1 on a TF that binds the promoter, 2+ upstream of that TF. Each step costs specificity.*

#### What these categories mean

Several of these sound alike and behave very differently. The distinction that most often decides an experiment is whether a route changes protein **abundance** or protein **activity** — a western blot only sees the first.

| Category | What it changes | Your drugs in it | Assay that sees it |
|---|---|---|---|
| TF: orthosteric ligand | **abundance** | Calcitriol, Paricalcitol, Dexamethasone, Rosiglitazone … +8 | RT-qPCR then western/ELISA |
| TF: endogenous ligand | **abundance** | Carbenoxolone, Tretinoin, Calcifediol, Ketoconazole … +5 | RT-qPCR then western/ELISA; also measure the endogenous ligand, which is what actually moved |
| TF: allosteric modulator | **abundance** | Demo-Pam-1 | RT-qPCR then western/ELISA |
| chromatin-level | **abundance, genome-wide** | Vorinostat, Romidepsin, Belinostat, Panobinostat … +4 | RT-qPCR plus a housekeeping panel; if everything moved, nothing specific happened |
| TF: stabiliser or glue | **abundance** | MRL-24 | RT-qPCR for your gene; note the DRUG acts on the factor's activity or level, not your gene's |

*Every category present here works through transcription, so all of them move your gene's product **abundance** and an RT-qPCR is the right first readout. The one exception in the taxonomy — a drug on the gene product itself, which changes activity only — did not appear in this run.*

**TF: orthosteric ligand** (12 found)

- *Yours:* Calcitriol, Paricalcitol, Dexamethasone, Rosiglitazone, Tretinoin, Bexarotene, Prednisolone, Alitretinoin … +4 more
- Binds the transcription factor's natural ligand pocket.
- *Why:* The best-characterised route, and where nearly all approved TF pharmacology lives.
- *Cost:* Only ~48 human TFs (the nuclear receptors) have such a pocket, and the factor's whole regulon moves with your gene.
- *Differs:* Clamps the receptor on or off regardless of endogenous signalling — contrast the allosteric and ligand-axis routes, which scale what is already happening.

**TF: endogenous ligand** (9 found)

- *Yours:* Carbenoxolone, Tretinoin, Calcifediol, Ketoconazole, Alfacalcidol, Cholecalciferol, Hydrocortisone, Talarozole … +1 more
- Changes how much natural ligand exists, by drugging the enzyme that makes or destroys it.
- *Why:* The enzyme is often far more druggable than the receptor, and some controlling enzymes are tissue-restricted, giving selectivity the receptor cannot.
- *Cost:* The ligand acts on every receptor it binds, everywhere — a systemic hormone change has a systemic phenotype. Homeostatic feedback also partly opposes you.
- *Differs:* Never touches the receptor. Anastrozole and finasteride work this way and appear in no other category.

**TF: allosteric modulator** (1 found)

- *Yours:* Demo-Pam-1
- Binds a site away from the ligand pocket, amplifying or damping the endogenous ligand's effect (PAM/NAM).
- *Why:* Allosteric sites are not conserved across a receptor family, so subtype selectivity is achievable where orthosteric chemistry has failed. Circadian and tissue patterning survive.
- *Cost:* Bounded maximal effect, and it does nothing without endogenous ligand present.
- *Differs:* Same receptor as the orthosteric route, opposite trade-off: less maximal effect, more selectivity.

**chromatin-level** (8 found)

- *Yours:* Vorinostat, Romidepsin, Belinostat, Panobinostat, Valproic acid, Azacitidine, Decitabine, Tazemetostat
- HDAC, DNMT and EZH2 inhibitors — genome-wide transcriptional effects.
- *Why:* Reliable positive control that your reporter *can* be induced.
- *Cost:* Specificity ~0.05. These move thousands of promoters. Almost never an acceptable therapeutic route on that basis alone.
- *Differs:* Not targeted at your gene in any sense; it is in the table so you can see what a non-specific option scores.

**TF: stabiliser or glue** (1 found)

- *Yours:* MRL-24
- Locks the factor in an active conformation, or glues it to a partner protein.
- *Why:* Needs no ligand pocket, so it reaches factors that orthosteric chemistry cannot touch at all.
- *Cost:* Mostly preclinical for TF targets; the approved examples are on other protein classes.
- *Differs:* Raises activity per molecule with protein levels unchanged — the opposite of a degradation blocker.

Not found here, and what you would be missing:

- **direct: gene or transcript** — Acts on your gene's own DNA or RNA — siRNA, antisense, splice-switching oligos, gene therapy.
- **direct: gene product** — Binds the protein your gene encodes — an antibody, enzyme inhibitor or receptor ligand acting on the product itself.
- **direct: predicted product binder** — A compound predicted to bind this gene's own protein and not recorded as acting on it - either because its chemistry resembles the protein's known ligands, or because it docks into the protein's pocket.
- **promoter DNA binder** — A small molecule with a documented preference for a sequence feature this promoter carries - a quadruplex, a GC-box, an AT tract - proposed because the feature is enriched against shuffles of this same sequence.
- **TF: inferred activity change** — A compound whose own LINCS signature moves this factor's known targets, found without any record that it binds the factor.
- **TF: abundance change** — A compound measured to change how much of the factor's own transcript exists, rather than what the factor does.
- **TF: predicted pocket binder** — A compound predicted by docking to fit the factor's ligand pocket, from a structure rather than from any measured activity.
- **TF: degradation blocker** — Stops the factor being destroyed, so more of it accumulates.
- **TF: known regulator** — Drugs that act on a factor CollecTRI or perturbation data records as raising or lowering your gene, whether or not the motif scan ranked that factor.
- **upstream signalling** — Drugs a kinase, phosphatase or receptor upstream of the transcription factor.
- **off-pathway pocket binder** — A docking run against a protein that is not one of the factors found at this promoter. Usually an upstream node someone chose to dock by hand.
- **measured phenocopy** — A compound whose effect on cell viability across the DepMap panel tracks what happens when this factor is knocked out, so it behaves as though it removes the factor.
- **empirical signature** — A compound observed to move your gene in a cell-based transcriptomic assay, mechanism unknown.
- **virtual screen hypothesis** — An approved drug predicted to bind one of your TFs, never tested against it.

### What these drugs are

Every route in this report is off-label by construction: the drug was approved for one thing and is being considered for another. The approved use is the baseline that departure is measured from, so it is listed here once per drug rather than repeated in each table above.

| Drug | What it is | Approved to treat |
|---|---|---|
| Alfacalcidol | vitamin D prodrug | renal osteodystrophy, rickets |
| Alitretinoin | RAR/RXR pan-agonist | Kaposi sarcoma skin lesions; chronic hand eczema |
| Azacitidine | DNA methyltransferase inhibitor | myelodysplastic syndrome, acute myeloid leukaemia |
| Bexarotene | RXR-selective agonist (rexinoid) | cutaneous T-cell lymphoma |
| Calcifediol | 25-hydroxyvitamin D | vitamin D deficiency |
| Calcitriol | active vitamin D, VDR agonist | hypocalcaemia and secondary hyperparathyroidism in renal failure |
| Cholecalciferol | vitamin D3, a precursor requiring activation | vitamin D deficiency and supplementation |
| Decitabine | DNA methyltransferase inhibitor | myelodysplastic syndrome |
| Dexamethasone | potent glucocorticoid receptor agonist | inflammation, cerebral oedema, myeloma, and as an antiemetic |
| Hydrocortisone | glucocorticoid and mineralocorticoid | adrenal insufficiency; topical dermatitis |
| Paricalcitol | selective VDR activator | secondary hyperparathyroidism in chronic kidney disease |
| Pioglitazone | PPARgamma agonist | type 2 diabetes |
| Prednisolone | glucocorticoid | asthma, autoimmune and inflammatory disease |
| Romidepsin | histone deacetylase inhibitor | cutaneous and peripheral T-cell lymphoma |
| Rosiglitazone | PPARgamma agonist (thiazolidinedione) | type 2 diabetes; restricted for cardiovascular risk |
| Tretinoin | RAR agonist, all-trans retinoic acid | acute promyelocytic leukaemia; acne topically |
| Vorinostat | histone deacetylase inhibitor | cutaneous T-cell lymphoma |

**No description established for 11**: Belinostat, Carbenoxolone, Demo-Pam-1, Demo-Pam-2, Ketoconazole, Liarozole, MRL-24, Panobinostat, Talarozole, Tazemetostat, Valproic Acid. That means this package could not find one offline, not that the compound has no approved use. Several are environmental agents and research chemicals from the CTD tables, which have no indication because they are not drugs.

Descriptions come from: 17 curated, 11 unknown. None of these is a label. The curated text is a summary written for orientation, ChEMBL's indication list includes uses studied but never approved, and a class line gives the kind of drug rather than its indication. Check the prescribing information before acting on any row.

### Graded candidates

Combined = geometric mean of TF confidence and drug quality. A weak TF prediction caps the grade no matter how good the drug.

| Grade | Combined | Drug | BBB? | TF | TF conf. | Drug qual. | Potency | Protein select. | Gene select. | Evidence | Flags |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 75.6 | Calcitriol | no | VDR | 0.624 | 0.916 | 1.0 | 1.0 | 0.673 | 0.75 | — |
| **A** | 74.6 | Paricalcitol | — | VDR | 0.624 | 0.892 | 0.9 | 1.0 | 0.673 | 0.75 | — |
| **B** | 67.5 | Dexamethasone | no | NR3C1 | 0.601 | 0.759 | 0.75 | 0.741 | 0.42 | 0.75 | NR3C1 controls ~1200 genes; expect broad transcriptional effects |
| **B** | 65.7 | Prednisolone | — | NR3C1 | 0.601 | 0.718 | 0.6 | 0.741 | 0.42 | 0.75 | NR3C1 controls ~1200 genes; expect broad transcriptional effects |
| **B** | 63.7 | Demo-Pam-1 | — | PPARG | 0.719 | 0.565 | 0.5 | 0.65 | 0.42 | 0.25 | no potency data; scored neutral; PAM: subtype selectivity typically be |
| **B** | 63.5 | Hydrocortisone | — | NR3C1 | 0.601 | 0.671 | 0.525 | 0.588 | 0.42 | 0.75 | NR3C1 controls ~1200 genes; expect broad transcriptional effects |
| **B** | 61.5 | Demo-Pam-2 | — | VDR | 0.624 | 0.606 | 0.5 | 0.65 | 0.673 | 0.25 | no potency data; scored neutral; PAM: subtype selectivity typically be |
| **B** | 59.1 | Bexarotene | — | RXRA | 0.753 | 0.464 | 0.675 | 0.588 | 0.42 | 0.75 | RXRA controls ~1200 genes; expect broad transcriptional effects; black |
| **B** | 59.1 | Rosiglitazone | — | PPARG | 0.719 | 0.485 | 0.7 | 0.741 | 0.42 | 0.75 | PPARG controls ~1200 genes; expect broad transcriptional effects; blac |
| **B** | 58.5 | Tretinoin | — | RARA | 0.755 | 0.454 | 0.825 | 0.364 | 0.42 | 0.75 | promiscuous: 6 annotated targets; RARA controls ~1200 genes; expect br |
| **B** | 57.6 | Alitretinoin | — | RARA | 0.755 | 0.44 | 0.725 | 0.364 | 0.42 | 0.75 | promiscuous: 6 annotated targets; RARA controls ~1200 genes; expect br |
| **C** | 54.5 | Pioglitazone | — | PPARG | 0.719 | 0.413 | 0.425 | 0.588 | 0.42 | 0.75 | PPARG controls ~1200 genes; expect broad transcriptional effects; blac |
| **F** | 22.9 | Valproic Acid | yes | REST | 0.606 | 0.087 | 0.0 | 0.29 | 0.673 | 0.75 | weak affinity (pChEMBL 4.2); promiscuous: 8 annotated targets; black-b |

**BBB?** is CNS MPO (Wager et al.) against the same scaled threshold the delivery section uses, so a drug cannot read "yes" here and "unlikely" there. It answers whether the molecule can reach brain at all: for a CNS target that decides whether a candidate is usable, and for a peripheral target it decides whether brain exposure is an unwanted side effect rather than the goal.

A **yes** predicts *passive diffusion only*, and is the weaker of the two answers. Active efflux by P-glycoprotein at the endothelium pumps an otherwise suitable molecule straight back out; known substrates are checked here and a hit overrides a pass outright.

**— is not a no.** 10 of 13 row(s) returned no physicochemical properties — common for DGIdb- and Open Targets-sourced rows, which carry a name and a mechanism but no structure. Nothing was measured and nothing is being claimed; look the compound up before excluding it on this column.

A **no** is the more reliable direction. A molecule failing on size, polarity or lipophilicity does not cross by diffusion, and no formulation fixes that short of bypassing the barrier — intrathecal delivery, convection-enhanced delivery, or focused ultrasound.

### Mechanism options: orthosteric vs allosteric

Allosteric sites are not constrained to bind the endogenous ligand, so they diverge across a receptor family where the orthosteric pocket cannot. That usually buys subtype selectivity — at the cost of a bounded maximal effect and a hard dependence on endogenous tone.

**RARA**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| orthosteric agonist | 2 | full | typically lower | no | no | Alitretinoin, Tretinoin |

*No approved allosteric modulator at this target; only orthosteric pharmacology, with the usual subtype-selectivity problem.*

**RXRA**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| orthosteric agonist | 1 | full | typically lower | no | no | Bexarotene |

*No approved allosteric modulator at this target; only orthosteric pharmacology, with the usual subtype-selectivity problem.*

**PPARG**

| Mechanism | Drugs | Max effect | Subtype selectivity | Outcompeted by ligand? | Needs endogenous ligand | Examples |
|---|---|---|---|---|---|---|
| PAM ⭐ | 1 | bounded (ceiling) | typically higher | no | YES | Demo-Pam-1 |
| orthosteric agonist | 2 | full | typically lower | no | no | Pioglitazone, Rosiglitazone |

*Both allosteric and orthosteric options exist at this target. The allosteric route usually buys subtype selectivity and preserved signalling dynamics at the cost of a bounded maximal effect and a hard dependence on endogenous ligand. Which matters more depends on whether you need a large fold change or a clean one.*

**Before pursuing any allosteric option, verify that the pathway has endogenous tone in your cell type.** A PAM with no ligand to modulate does nothing at all, and that is a measurement, not a prediction — nothing in this package can check it for you.

### Changing the natural ligand instead of the receptor

A nuclear receptor's activity is set by how much of its endogenous ligand is present, and the enzymes that make or destroy that ligand are often far more druggable than the receptor. Anastrozole, finasteride and metyrapone all work this way — none of them bind a receptor at all.

| Strategy | Effect on the receptor |
|---|---|
| supplement ligand | up |
| inhibit catabolism | up |
| inhibit synthesis | down |
| sequester ligand | down |

- Endocrine axes are homeostatic. Feedback will partly oppose you, so the effect is smaller and slower than the enzyme inhibition suggests, and it can rebound on withdrawal.
- Standing caveat: the ligand acts on every receptor it binds, in every tissue. Changing a systemic hormone level has a systemic phenotype — aromatase inhibitors cause bone loss because estradiol does far more than drive your gene.
- Supplementation only works if ligand supply is limiting. If the receptor is already saturated in your tissue, adding more does nothing — check the axis has headroom first.

#### RARA — ligand: all-trans retinoic acid

Controlling enzymes: CYP26, CYP26A1/B1

| Drug | Strategy | Acts on | Local? | Status |
|---|---|---|---|---|
| Talarozole | inhibit catabolism | CYP26A1/B1 | **yes** | clinical |
| Tretinoin (1971) | supplement ligand | — | no | approved |
| Liarozole | inhibit catabolism | CYP26 | no | clinical |

The chain, link by link:

- Talarozole blocks CYP26A1/B1 → more all-trans retinoic acid → more RARA activity → up at your gene
- Tretinoin supplies all-trans retinoic acid → more RARA activity → up at your gene
- Liarozole blocks CYP26 → more all-trans retinoic acid → more RARA activity → up at your gene

Drug-specific caveats:

- **Talarozole** — CYP26 expression is induced by retinoic acid itself, so the system compensates
- **Liarozole** — NOT approved (clinical)

- At least one route is tissue-selective: the controlling enzyme is not expressed everywhere, so the ligand change is local rather than systemic. That is selectivity the receptor itself cannot offer, and the main reason to prefer this approach.

#### VDR — ligand: calcitriol

Controlling enzymes: CYP24A1

| Drug | Strategy | Acts on | Local? | Status |
|---|---|---|---|---|
| Calcifediol (2016) | supplement ligand | — | no | approved |
| Ketoconazole (1981) | inhibit catabolism | CYP24A1 | no | approved |
| Alfacalcidol (1978) | supplement ligand | — | no | approved |
| Cholecalciferol (1940) | supplement ligand | — | no | approved |

The chain, link by link:

- Calcifediol supplies calcitriol → more VDR activity → up at your gene
- Ketoconazole blocks CYP24A1 → more calcitriol → more VDR activity → up at your gene
- Alfacalcidol supplies calcitriol → more VDR activity → up at your gene
- Cholecalciferol supplies calcitriol → more VDR activity → up at your gene

Drug-specific caveats:

- **Ketoconazole** — profoundly non-selective across CYPs; hepatotoxicity

#### NR3C1 — ligand: cortisol

Controlling enzymes: CYP11B1, CYP17A1/CYP11B1, HSD11B2

| Drug | Strategy | Acts on | Local? | Status |
|---|---|---|---|---|
| Carbenoxolone (1962) | inhibit catabolism | HSD11B2 | **yes** | approved |
| Hydrocortisone (1952) | supplement ligand | — | no | approved |

The chain, link by link:

- Carbenoxolone blocks HSD11B2 → more cortisol → more NR3C1 activity → up at your gene
- Hydrocortisone supplies cortisol → more NR3C1 activity → up at your gene

Drug-specific caveats:

- **Carbenoxolone** — the classic tissue-selective example, but HSD11B2 inhibition in kidney causes sodium retention and hypertension

*Same axis, opposite direction — do not confuse these with the list above: Levoketoconazole (inhibit synthesis), Osilodrostat (inhibit synthesis), Metyrapone (inhibit synthesis).*

- At least one route is tissue-selective: the controlling enzyme is not expressed everywhere, so the ligand change is local rather than systemic. That is selectivity the receptor itself cannot offer, and the main reason to prefer this approach.

**No ligand axis for:** RXRA, PPARG, SP1, CREB1, REST, NFKB1, HNF4A

Most transcription factors are ligand-independent, so this is usually a real absence rather than a gap in the table. Receptors with a curated axis here: AHR (kynurenine), AR (dihydrotestosterone), ESR1 (estradiol), NR1H4 (bile acids), NR3C1 (cortisol), RARA (all-trans retinoic acid), THRB (triiodothyronine), VDR (calcitriol).

### Locking the factor into an active state

Nuclear receptors have a ligand-binding pocket; almost no other TF does, which is why the graded tables keep returning NR agonists. Stabilisation sidesteps the pocket: every factor has protein-protein interfaces and degradation machinery, and both are druggable.

**PPARG** — has a ligand pocket · **a curated route exists**

| Drug | Mechanism | Binds | Raises | Status |
|---|---|---|---|---|
| MRL-24 | conformational stabiliser | PPARG | activity | preclinical |

- PPARG has a ligand-binding pocket, so conventional agonists are available and are usually the simpler route. Stabilisation is worth considering where you need activity without the full agonist phenotype — the PPARG Ser273 story is the example.

#### No stabilisation route for 3 factor(s)

CREB1, REST, SP1

These have no ligand pocket, no direct drug, and **nothing curated against them or their degradation machinery**. For these factors stabilisation is a research programme, not a route you can take: making it real needs the E3 ligase that degrades the factor *and* a ligand for that ligase, and neither is published here.

The mechanism class is real — it has shipped approved drugs, listed once below — but those act on **other proteins entirely**, and none of them is evidence that any factor above is tractable.

*Conformational stabilisers raise activity with protein levels unchanged; degradation blockers raise levels. A western blot cannot see the first kind — pick the assay to match the mechanism.*

### Routes a promoter analysis cannot see

Everything above acts on transcription. These act after it, around it, or on a different gene entirely. For some genes that is where the whole problem is: if the transcript is destroyed, the allele is silenced, or the protein is degraded too fast, more transcription does not help and may not be possible.

#### Could a relative do the job instead?

Whether a paralog could substitute for this gene is unanswered. It is worth asking explicitly, because the question 'restore this function' has more answers than 'raise this transcript'.

- DEMOGENE is not in the curated redundancy table, which is short and certainly incomplete. Absence here is not evidence that no paralog substitutes — it means nobody has written one down in this file.

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

*Accounting for all 13 candidates: 0 already linked, 1 novel and non-dominated, 11 novel but beaten on every axis by another candidate, 1 with no advantage. All of them are in results.json under `improvements.all_candidates`.*

#### Not previously linked (no incumbent to beat)

| Candidate | Route | TF | Strength | Specificity | Efficiency | Wins on | Trade-offs |
|---|---|---|---|---|---|---|---|
| Calcitriol | via TF | VDR | 0.744 | 0.673 | 0.963 | stronger activation, higher specificity, superior efficiency | none |

*No incumbent exists, so these cleared no bar — they are simply the candidates that scored above 0.5 on at least one axis. Hypotheses to test.*

Novel but beaten on every axis by another candidate (shown so nothing is hidden): Paricalcitol, Dexamethasone, Prednisolone, Demo-Pam-1, Hydrocortisone, Demo-Pam-2, Bexarotene, Rosiglitazone, Tretinoin, Alitretinoin, Pioglitazone

#### Non-dominated set (nothing beats these on every axis)

- **Calcitriol** (NEW link) — no incumbent for this gene; nothing to beat

### Strength of evidence for each factor route

Two links between a drug and your gene come from evidence and differ from row to row: whether the drug engages its target at a dose that can be reached, and whether the factor is likely to act on this promoter. Each is a 0 to 1 heuristic and not a calibrated probability. They are shown as they are, and they are not multiplied into a chance that the gene moves.

| Drug | TF | Drug engagement | Promoter activation | Weaker link |
|---|---|---|---|---|
| Calcitriol ★ | VDR | 0.95 | 0.587 | promoter activation |
| Tretinoin | RARA | 0.921 | 0.678 | promoter activation |
| Alitretinoin | RARA | 0.876 | 0.678 | promoter activation |
| Bexarotene | RXRA | 0.854 | 0.677 | promoter activation |
| Rosiglitazone | PPARG | 0.865 | 0.654 | promoter activation |
| Paricalcitol | VDR | 0.95 | 0.587 | promoter activation |
| Dexamethasone | NR3C1 | 0.887 | 0.571 | promoter activation |
| Demo-Pam-1 | PPARG | 0.775 | 0.654 | promoter activation |
| Pioglitazone | PPARG | 0.741 | 0.654 | promoter activation |
| Prednisolone | NR3C1 | 0.82 | 0.571 | promoter activation |
| Demo-Pam-2 | VDR | 0.775 | 0.587 | promoter activation |
| Hydrocortisone | NR3C1 | 0.786 | 0.571 | promoter activation |
| Valproic Acid | REST | 0.55 | 0.575 | drug engagement |

★ = named in the comparison against the standard of care above, pinned here so the two sections cannot disagree.

Order candidates by tier and score, as the tables above do. These columns say where the evidence behind a route is thin.

### Approved drugs acting on these factors

| TF | TF role | Do what to this TF | Drug | Action | Drug direction | First approval | Source |
|---|---|---|---|---|---|---|---|
| RARA | activator | increase | Alitretinoin | AGONIST | up | 1999 | DEMO |
| RARA | activator | increase | Tretinoin | AGONIST | up | 1971 | DEMO |
| RXRA | activator | increase | Bexarotene | AGONIST | up | 1999 | DEMO |
| PPARG | activator | increase | Demo-Pam-1 | POSITIVE ALLOSTERIC MODULATOR | up | 2015 | DEMO |
| PPARG | activator | increase | Pioglitazone | AGONIST | up | 1999 | DEMO |
| PPARG | activator | increase | Rosiglitazone | AGONIST | up | 1999 | DEMO |
| VDR | activator | increase | Calcitriol | AGONIST | up | 1978 | DEMO |
| VDR | activator | increase | Demo-Pam-2 | POSITIVE ALLOSTERIC MODULATOR | up | 2018 | DEMO |
| VDR | activator | increase | Paricalcitol | AGONIST | up | 1998 | DEMO |
| REST | repressor | decrease | Valproic Acid | INHIBITOR | down | 1978 | DEMO |
| NR3C1 | dual | increase | Dexamethasone | AGONIST | up | 1958 | DEMO |
| NR3C1 | dual | increase | Hydrocortisone | AGONIST | up | 1952 | DEMO |
| NR3C1 | dual | increase | Prednisolone | AGONIST | up | 1955 | DEMO |

---

### How to read this

A PWM hit means the sequence *resembles* a binding site. It is a hypothesis generator, not evidence of binding, and not evidence that binding would increase transcription of this gene in your cell type. Rank order here is a heuristic (weights printed in the JSON output), not a validated predictor. Every row needs wet-lab confirmation: promoter-reporter assay for the element, agonist dose-response with RT-qPCR/protein readout for the drug, and a site-directed mutant of the element to show the drug effect runs through it.
