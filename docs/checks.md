# Checks on candidate drugs

[Documentation index](index.md)

Two commands read saved results and ask whether the candidates deserve
belief. Neither changes a rank or a score, and neither runs the pipeline.

- `exposure` asks whether a person taking the approved dose reaches the
  concentration at which the evidence was seen.
- `retrospective` asks whether the pipeline's routes point at drugs already
  known to move the gene.

## Exposure against need

```bash
python -m promoter_tf_drugs.exposure results/MYGENE_up/results.json \
    --cmax-table data/human_cmax.tsv --top 60
```

It takes the top candidate drugs from the result and sets two numbers against
each other, then writes `exposure.json` and `exposure.md` next to the result.

**Needed** is the lowest LINCS dose at which the drug moved this gene. When
there is none, it is the drug's potency at the factor the route goes through,
from its pChEMBL. The second is weaker, because it assumes the gene responds
when the factor does, and the report says which was used.

**Reached** is a published human peak blood level (Cmax). It comes from your
table first, then from human plasma records in ChEMBL.

The ratio gives one of four labels: reaches (at or above), within tenfold,
far below, or not assessable.

### What the label means

Most values are total drug in plasma, and only the free fraction acts. So
"far below" is firm: total already under the target means free is lower
still. "Reaches" is the weakest claim, because free drug can be a few percent
of the total. Put an unbound fraction (`fu`) in your table and the free level
is used instead, and labelled as such.

A drug with no exposure value, or with no known effective concentration, is
not assessable. That is not the same as unreachable, and it is never reported
as far below.

### The Cmax table

A tab or comma separated file with `name` and either `cmax_nM`, or `cmax`
and `unit` (and `mw` for a mass unit such as ng/mL). Optional columns are
`fu` (unbound fraction between 0 and 1), `dose`, `source` and `note`.

```text
name        cmax    unit    mw      fu      dose        source
Imatinib    2.6     ug/mL   493.6   0.05    400 mg      (your source)
```

The package ships no values, because none could be checked against a label
here. Use a source you trust and name it in `source`. A drug listed twice
keeps its larger value. A row that cannot be read is skipped with a warning
and never read as zero.

### Where ChEMBL fills the gap, and how much to trust it

Without `--no-chembl`, drugs the table lacks are looked up in ChEMBL's human
Cmax records. A record is kept only when its description says human and names
no tissue, fluid other than plasma, animal, drug interaction or special
population. Records for one drug that differ by more than thirty fold are not
averaged: the drug is left not assessable and the range is shown.

Those filters exist because of what ChEMBL returned. Calcitriol's "human"
records came from ketoconazole and ritonavir interaction studies and read as
4.5 and 18 uM, thousands of times the real level, which produced a false
"reaches". Coverage is also thin. On saved results for BUB1B, CCR5 and HTT,
about one drug in thirty had both a usable exposure and a known effective
concentration. Read a run mostly as a list of what could not be checked.

### Only FDA approved drugs are assessed

By default the dose check assesses drugs in the FDA Orange Book and leaves out
the rest (`--include-non-fda` keeps them). A research compound has no human
dose to compare with, and a drug approved only outside the United States has
no FDA label to read one from. Biologics are not in the Orange Book, so they
are left out too. With no Orange Book on disk the filter is not applied.

### Why so many drugs are not assessable

The report lists the reasons with counts. On the BUB1B results, before the
FDA filter, 56 of the 60 candidates could not be compared, for these reasons:

- About 23 are not approved drugs: LINCS research compounds such as
  roscovitine, triptolide and BRD-K ids. They have no human dose, and no
  ChEMBL record under that name. The FDA filter now removes them.
- About 23 are approved drugs for which ChEMBL holds no Cmax record at all,
  among them tretinoin, ursodiol, levothyroxine and selumetinib. ChEMBL stores
  a Cmax only when a paper deposited it, so coverage is uneven.
- About 8 have Cmax records, and none is from people in plasma. Most are mouse
  or rat studies from medicinal chemistry papers, as for obeticholic acid,
  fulvestrant and vorinostat.
- One has human records more than thirty fold apart (belinostat, whose oral
  and intravenous doses are mixed), so none is used.
- 30 have no concentration to compare with, which is a separate gap.

Three changes came from reading these. A record that says "patient with solid
tumor" or "lung cancer" names who was dosed and no longer counts as a tissue
measurement; that rule had dropped bexarotene's only human record. A receptor
screen's measured potency now counts as the concentration to compare with,
where before only a LINCS dose or a graded drug did, so tretinoin and
odevixibat have a target concentration. And a drug the approved table lacks
is looked up by the ChEMBL id the screen found, or by an exact name match.
Together these moved the BUB1B run from 2 drugs with a verdict to 4. The rest
are limited by the data, not the code: more values need a source of human
exposure beyond ChEMBL, which is what `--cmax-table` is for.

Label text from openFDA was tried and not used. The pharmacokinetics sections
are flattened tables where a Cmax can belong to an interaction study or a
different dose, and a wrong number read silently is worse than none.

### Other limits

- LINCS doses are nominal concentrations in a dish, and plasma is not a tissue.
- A single-dose Cmax understates steady state for a drug that accumulates.
- Potency at a factor is not the concentration at which the gene responds.

## Retrospective check

```bash
python -m promoter_tf_drugs.retrospective results/GENE1_up results/GENE2_down \
    --out results/retrospective
```

Give it result folders or `results.json` files. It reads only saved results
and the local LINCS files, and writes `retrospective.json` and
`retrospective.md`. When a gene has several saved results in the same
direction, the newest is used and the rest are listed, so a rerun is not
counted as another gene.

### The two answer keys

**CTD** is read from the saved result. It records chemical effects on gene
expression from the literature, at any dose, in any cell line. It plays no
part in route scores. Most of what it holds for a gene is a toxicant, not a
drug anyone would repurpose, so few of its entries are in the approved library.

**LINCS** is read from the local signature files. A drug counts as wanted when
it moved the gene the wanted way in at least two signatures and at least twice
as often as the other way, and as the opposite when the reverse holds. Mixed
signs and single sightings are set aside. Only drugs the files contain are
compared, because a drug never tested has not been seen not moving the gene.
For many genes this leaves too few positives to score, and the report says so.

### Which routes are graded

Routes through a transcription factor, an upstream signalling node, a
chromatin class, or promoter DNA are the predictions. Routes built from LINCS
(a factor whose own transcript moves, or whose activity is inferred from
signatures) are left out of the group graded against LINCS, so nothing is
scored on its own input. LINCS signature routes are scored against CTD only.
Routes that bind the gene's own product, and virtual screen hypotheses, are in
no group, because the first says nothing about expression and the second is
labelled using CTD.

### The two tests

Both are areas under the ROC curve, where 0.5 is chance. A permutation test
gives how often shuffled labels do as well.

- Reported or not: among approved drugs, do the ones the key records as moving
  the gene the wanted way score higher than the ones it has no record of? A
  drug with no route scores zero, so reaching a drug at all counts.
- Right direction: among drugs the key records as moving the gene at all, do
  the wanted ones score higher than the opposite ones? Both groups were seen
  moving the gene, so being well studied does not help.

The first test also reports a baseline that knows nothing about the gene, the
number of routes a drug has. A predicted score that does not beat it may only
reflect drugs that appear everywhere. Beside every score is how many of the
positives have a route at all. A low score with few positives reached is a
coverage result, which is different from reaching them and ranking them badly.

### Positive controls

The module holds a short table of genes whose regulation by an approved drug
is textbook pharmacology (CYP3A4 and rifampicin, FKBP5 and glucocorticoids,
and so on). When a run is of one of those genes in the listed direction, the
report says where the expected drugs landed among the candidates, by which
route, or that they are absent. The table is the author's judgement of
established pharmacology, written before the runs were made, with no
citations. Read it critically and change it in `POSITIVE_CONTROLS`.

### What it cannot show

- The keys are noisy. A perfect predictor would not score 1.0.
- It tests expression, not benefit. Moving a gene is not treating a disease.
- It tests the past. The value of the pipeline is in the drugs no record
  covers, which this cannot score.
- A handful of genes carries wide uncertainty. The range and the count above
  chance are shown beside the mean, and a run with fewer than five positives or
  five negatives is not scored.
- Routes that remain in the LINCS-graded group may still use scores fitted
  with signature data somewhere upstream. That is not traced.

## First results

Run on 17 October 2026 against an uncommitted working tree, so treat the
numbers as a snapshot. Thirty-six genes were run with the standard options, 6
chosen for textbook drug regulation and 30 others (xenobiotic response,
inflammation, lipid and stress genes, plus CCR5, HTT, BUB1B, SCN1A, GJC1 and
ABCA4). The outputs are in `results/retrospective_36`, and the runs behind
them are in `results/validation` and `results`.

**Recall is low.** Of 1,666 approved drugs that CTD records as moving these
genes the wanted way, 142 (9 percent) appear in any predicted route. Of the six
positive controls, two have a known drug among the candidates at all:

| Gene | Known answer | Where it landed |
|---|---|---|
| CYP24A1 | vitamin D analogues | rank 2 of 163, via VDR |
| ANGPTL4 | fibrates, glitazones | rank 63 of 114, via PPARA |
| CYP3A4 | rifampicin | absent from 111 |
| FKBP5 | glucocorticoids | absent from 145 |
| LDLR | statins | absent from 185 |
| CYP1A1 | omeprazole and relatives | absent from 161 |

**The ranking cannot separate known regulators once drug fame is held equal.**
Against CTD the predicted routes average an AUC of 0.53 on the
reported-or-not test over the 34 genes that could be scored, with 26 of 34
above chance. A baseline that only counts how many other genes record a drug
averages 0.93. Among only the drugs recorded for at least three other genes,
where fame no longer separates the groups, the routes average 0.53 (chance).
On the right-direction test, where both groups were seen moving the gene, they
average 0.50. The apparent skill is popularity. The LINCS key could be scored
for 15 genes and showed the same.

**One cause of misses is the scan window.** For CYP3A4 and FKBP5 the receptor
that mediates the drug (PXR, the glucocorticoid receptor) has no site in the
default window of 1,000 bases upstream and 200 downstream, so no route can
reach the drug through it. Rerunning CYP3A4 with 8,500 bases upstream put the
receptor in the scan and rifampicin in the candidates, through NR1I2, at rank
131 of 219. Widening the window fixed reach and not rank. That is one gene's
evidence.

**Another is routes the model does not draw.** Statins raise LDLR by lowering
cholesterol, which activates SREBP2, and no step of that reaches a drug in the
network the pipeline reads.

## What the ranking review found

The review looked at why the lists are ordered as they are. It did not find a
way to improve the order against these keys, and it changed no weight to chase
them. It found and fixed three defects, listed in the compatibility note in
[Inputs and outputs](data.md).

- **The top of every list is a tie.** In the 12 genes first reviewed, 20 to 30
  of the top 30 drugs shared an exact tier and score with at least one other,
  and the 12 lists had between 4 and 43 distinct levels for 40 to 220 drugs.
  After the fixes the share of tied top-30 drugs fell from 0.95 to about 0.78
  on the same saved results. Order within a tie is arbitrary, so a rank such as
  131 of 219 can mean a position inside a block of identical routes.
- **Half the factor drugs skipped the potency lookup.** Their ids were RxNorm
  or NCIt, and the lookups need a ChEMBL id, so they scored as average and had
  their safety lookups skipped. For CYP3A4 the share with a measured potency
  rose from 27 to 49 of 83 once ids and targets were found.
- **Receptor-screen ligands ignored potency and direction.** All scored 0.7
  whatever their measured affinity, and a drug known to move the factor the
  wrong way ranked with one known to move it the right way.

Fitting a rule to CTD does not help. A category weight fitted on half the genes
and scored on the other half raised the candidate AUC from 0.62 to 0.67, a gain
whose interval excluded zero. That gain disappears once fame is held equal (0.47
to 0.49, both chance) and does not appear on the right-direction test (a
difference of +0.01, interval -0.07 to +0.09). It was preferring well-known drug
classes, so it was not adopted.

One finding is left open. Rifampicin's potency at PXR is real (pChEMBL 5.9 over
16 measurements) but moderate, as it is for most ligands of that receptor, and
the potency component reads anything under pChEMBL 9 as less than strong. A
drug whose blood level exceeds its potency is a better candidate than one that
is potent and never reaches it, which needs the exposure data the dose check
lacks. Until that exists, potency is scored on an absolute scale.

None of this was used to change a score or default to make an example look
stronger. As it stands, a candidate list from this pipeline is a reading order,
a drug whose factor sits outside the window cannot appear, and a known
regulator can be present and still ranked well below the middle.
