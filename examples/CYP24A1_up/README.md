# CYP24A1, raise it

A positive control. CYP24A1 is switched on by vitamin D acting through the vitamin
D receptor (VDR), so calcitriol and its analogues are the textbook answer. The
question is whether the pipeline finds them and how it ranks them.

Command (the launcher's options without docking, on 21 October 2026):

```sh
python -m promoter_tf_drugs --gene CYP24A1 --direction up \
  --pfm-file data/jaspar2024.txt --regulon-file data/chip_regulon.tsv \
  --chip-bed data/remap2022.bed.gz --upstream 1000 --downstream 200 --pvalue 1e-5 \
  --all-tss --tf-weight 0.7 --cell-type-check --predict-direction \
  --upstream-network --upstream-depth 3 --regulator-route --regulator-steps 2 \
  --lincs --genetics --screen --design-screen --design-modality --promoter-dna \
  --screen-product --find-analogs --sensitivity-scan --out results/CYP24A1_up_v2
```

What the run shows:

- **The verdict is "degraded".** The default window showed no signal above chance (13
  sites), so the scan was retried at -5000/+1000, which shows a weak signal (89 sites
  against about 65 expected by chance). The report says both at the top and does not read
  the ranking as evidence.
- **The textbook answer is found and ranks first.** Cholecalciferol, paricalcitol,
  doxercalciferol, ergocalciferol, calcifediol and calcitriol are rows 1 to 6 of the 99
  approved drugs shown, through the known regulators route. VDR is a regulator of CYP24A1
  in both CollecTRI and the knockout data, and a regulator that both name is tier 3,
  provisionally.
- **Read that ranking with care.** The tier 3 rule for known regulators was chosen after
  looking at 36 genes, and CYP24A1 is one of them, so this run is not a test of it. In an
  earlier snapshot of this same gene, before the rule, the same drugs were tier 4 and ranked
  10 to 15 of 65, behind nine tier 3 routes built from LINCS signatures and degradation
  blockers. The retrospective check in `docs/checks.md` found that the routes rank drugs no
  better than chance once a drug's fame is held equal, and held-out genes are still needed.
- **The CTD route adds rows.** 13 of the approved drugs shown are reported effects from
  CTD, marked as such, with no mechanism.
- **Nothing here shows that any listed drug raises CYP24A1.** The routes are
  hypotheses generated from annotation.
