# CYP24A1, raise it

A positive control. CYP24A1 is switched on by vitamin D acting through the vitamin
D receptor (VDR), so calcitriol and its analogues are the textbook answer. The
question is whether the pipeline finds them and how it ranks them.

Command (the launcher's options without docking, on 20 October 2026):

```sh
python -m promoter_tf_drugs --gene CYP24A1 --direction up \
  --pfm-file data/jaspar2024.txt --regulon-file data/chip_regulon.tsv \
  --chip-bed data/remap2022.bed.gz --upstream 1000 --downstream 200 --pvalue 1e-5 \
  --all-tss --tf-weight 0.7 --cell-type-check --predict-direction \
  --upstream-network --upstream-depth 3 --regulator-route --lincs --genetics \
  --screen --design-screen --design-modality --promoter-dna --screen-product \
  --find-analogs --sensitivity-scan --out results/CYP24A1_up_final
```

What the run shows:

- **The verdict is "degraded".** The motif scan found 13 sites where shuffled copies
  of the window give about 11, so about 86% of them are expected to be chance. The
  report says so at the top and does not read the ranking as evidence.
- **The textbook answer is found, but not first.** Calcitriol, paricalcitol,
  cholecalciferol, doxercalciferol, ergocalciferol and calcifediol appear at ranks
  10 to 15 of the 65 approved drugs shown, through the known regulators route
  (CollecTRI and a knockout experiment name VDR as a regulator of the gene). The
  promoter scan did not rank VDR at all. The upstream search also lists some of
  these drugs, but only as acting on VDR upstream of other factors, at a lower
  score.
- **They rank behind nine tier-3 routes.** Those are degradation blockers such as
  dimethyl fumarate, and LINCS-derived effects on factor transcripts such as
  dasatinib at CEBPD. Tier decides the order first, and a tier-4 route cannot beat a
  tier-3 one. The retrospective check in `docs/checks.md` found that the routes rank
  drugs no better than chance once a drug's fame is held equal, so it gives no reason
  to prefer those routes. This run is one more case where the textbook answer ranks
  below them.
- **Nothing here shows that any listed drug raises CYP24A1.** The routes are
  hypotheses generated from annotation.
