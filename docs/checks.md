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

### The FDA label table

`python scripts/build_fda_cmax.py` builds `data/fda_label_cmax.tsv` from FDA
labels, and the dose check reads it when it exists, below any table you name.
With no arguments it looks up the FDA approved drugs in the `results.json` files
under `results/`; `--from-results`, `--drugs` and `--all-orange-book` choose
others. The launcher runs it on the finished results before the dose check.
openFDA allows 1,000 requests a day without a key and 120,000 with one (put it in
`OPENFDA_API_KEY`).

It reads sentences in the pharmacokinetics section that state a peak:
"Cmax", "peak plasma concentration", "maximum concentration". It takes each
number with a concentration unit, with the dose, whether it is single dose or
steady state, and the sentence. It refuses a sentence or heading about renal or
hepatic impairment, elderly, children, sex or race, a drug interaction, food, or
a metabolite. A drug is written to the table when its remaining figures agree
within thirty fold, and the row holds their median with the range, the label, the
dose and the sentence in `note` so it can be checked. A protein binding statement
gives the unbound fraction, from which the dose check reports a free level. A
range such as "6% to 12% bound" gives none.

Anything else goes to `data/fda_label_cmax_review.tsv` with the reason: figures
that disagree, a mass unit with no molecular weight, or every sentence refused.

Coverage is partial, and that is the result and not a defect. On 170 FDA
approved drugs from three saved runs, 38 got a value, 14 were set aside and 118
had nothing to read. The labels state a number for the peak in prose for about
one drug in four. The rest give only a time to peak ("C max achieved within 2 to 4
hours" for imatinib, with no value), put the figures in a table (not read: the
tables openFDA carries for these sections are mostly interaction ratios, the
exact trap the rules above avoid), or say nothing. A value that is read is still
one figure from one population at one dose, so treat a row as a lead to check
against the label and not as a measured blood level. The text reader can also
be wrong in ways the tests do not cover.

Label text from openFDA was tried first and dropped, because a wrong number read
silently is worse than none. The reader above is that attempt with the refusals
added.

### Other limits

- LINCS doses are nominal concentrations in a dish, and plasma is not a tissue.
- A single-dose Cmax understates steady state for a drug that accumulates.
- Potency at a factor is not the concentration at which the gene responds.

## Where the literature's drugs are lost

```bash
python scripts/positive_fate.py <results.json files or folders> --out results/retrospective_fate
```

The retrospective check finds that routes recover few of the approved drugs CTD
records as moving a gene. This script follows each such drug that has no route
and says where it stops. For each, it looks up the drug's ChEMBL mechanism
targets (every library row of the drug, since a mechanism sits on one salt
form) and places it at the first stage that explains the miss.

On the 36 saved runs, 1,677 approved drugs are recorded by CTD as moving their
gene the wanted way:

| Stage | Drugs | Share |
|---|---|---|
| found by a predicted route | 144 | 8.6% |
| acts on a factor the run ranked, yet has no route | 2 | 0.1% |
| a target is a transcription factor the scan did not rank | 206 | 12.3% |
| non-factor target within three steps upstream of a factor the search covered | 335 | 20.0% |
| the same, but only of a ranked factor the search did not cover | 59 | 3.5% |
| target more than three steps upstream | 291 | 17.4% |
| target not in the signalling network | 265 | 15.8% |
| no mechanism target in ChEMBL | 375 | 22.4% |

A drug with several targets takes the stage nearest to being found. Reading it:

- **Drug lookup is not the bottleneck.** Once a factor is ranked, its drugs are
  found: only 2 drugs act on a ranked factor and lack a route. Misses come from
  which factors are ranked and how far upstream the search goes.
- **About 38% cannot be reached by a target-based method at all.** Those are
  drugs with no ChEMBL target (urethane, hydrogen peroxide, zinc salts, DMSO)
  and drugs whose target is not in the signalling network (RNA polymerase,
  DHODH, ERG11). CTD records them for toxicology or stress effects. Recall
  cannot approach 100% by improving the routes.
- **The upstream groups (335 and 59) are an upper bound.** They count a
  non-factor target within three steps of a ranked factor without checking that
  the path has the right sign or that the drug acts in the right direction. The
  targets are PTGS2, HTR2A, DRD2, PTGS1 and KIT. The upstream search covers only
  the top eight factors, which is why 59 of these reach a ranked factor outside
  it.
- **Transcription factors that were never ranked.** 208 missed drugs have a
  transcription factor among their targets, mostly NR3C1, ESR1, PGR, AR, PPARA,
  ESR2 and VDR. For 97 of them (47%) the gene is a known CollecTRI target of
  that factor, and for 42 (20%) the gene responded when it was knocked out; 110
  (53%) have one or the other. That is 6.6% of all 1,677 drugs, which would lift
  recall from 8.6% to about 15% through a route that starts from known
  regulators and not from motifs at the promoter. CollecTRI and CTD both draw
  on the published literature, so some of that overlap is not independent, and
  the 6.6% is an estimate of what such a route could reach, not a measured gain.

### Why the upstream group is missed

`python scripts/upstream_misses.py <results.json files>` replays the filters of
the upstream search on the 598 drugs whose target is within three steps of a
ranked factor (the two upstream rows above, plus transcription factors that
are also upstream nodes). The first column is the search as it runs, with eight
factors. The second asks every ranked factor, with or without a drug of its own
and with no cap on how many:

| First filter that removes the drug | Eight factors searched | Every ranked factor |
|---|---|---|
| reaches only a factor outside the eight searched | 82 (13.7%) | 0 |
| no searched factor within three steps | 16 (2.7%) | 16 (2.7%) |
| no path whose signs give the effect the goal needs | 65 (10.9%) | 32 (5.4%) |
| no recorded mechanism at that target in the approved-drug index | 26 (4.3%) | 46 (7.7%) |
| the drug moves the node the wrong way for that path | 252 (42.1%) | 306 (51.2%) |
| passes every filter, so the search's limits left it out | 157 (26.3%) | 198 (33.1%) |

- **Protein families were missing from the index, and now are included.**
  ChEMBL records aspirin and ibuprofen at a "PROTEIN FAMILY" target (PTGS1 and
  PTGS2), and haloperidol at the dopamine receptor family. The index used to
  keep only single proteins and complexes, which left 156 drugs without an entry
  (26.1%). It now keeps families too, lists the drug at every member, and marks
  it with `target_type = "PROTEIN FAMILY"` and a source ending in "family
  target". The drug is not known to act on one member alone, so read such a
  route as "acts on this protein family".
- **Many of the drugs that moved into the wrong way column were hidden before.**
  With families included, aspirin and ibuprofen are found at PTGS2, and their
  direction there is down, so the search correctly asks whether lowering PTGS2
  helps. For most of these drugs it does not, so the count of wrong way drugs
  rose from 156 to 252.
- **Widening the search adds at most 41 drugs.** Asking all ranked factors moves
  the drugs that pass every filter from 157 to 198, which is 6.9% of the 598. It
  is an upper bound, because the per-factor limits are not replayed. Most of
  the 198 are at NR3C1, ESR1, AR, ESR2 and VDR, which have dozens of approved
  drugs each, so the limits would cut most of them. Searching more factors is
  therefore a small gain compared with the limits, and with starting from known
  regulators of the gene.
- **The limits hide drugs at busy nodes.** The search keeps ten druggable nodes
  per factor and five drugs per node. This is a count of drugs the filters did
  not remove, not a measurement of which limit did.
- **Direction and sign filters remove about half.** Those are the search
  working as designed. They can still be wrong when a factor's role is unknown,
  because a role that is neither activator nor dual is treated as a repressor.

**Compatibility:** `Drug` has a new field, `target_type`, which is `None` unless
the target is a protein family. The approved-drug index returns more drugs, at
more genes, than before, so upstream routes can list drugs such as aspirin that
an older run did not. Saved results are unchanged and still read correctly.

### What starting from regulators adds

`python scripts/regulator_gain.py <results.json files>` builds the
`--regulator-route` routes for each of the 36 saved runs, adds them to the run's
routes, and scores the result against CTD the way the retrospective does. It
compares each source alone. The check covers the runs that have a CTD record
(1,666 approved drugs that CTD records as moving a gene the wanted way).

| Source | Runs with new routes | New drugs | CTD wanted | CTD opposite | Wanted drugs with a route |
|---|---|---|---|---|---|
| none (routes already in the runs) | 36 | 3,109 | 144 | 91 | 143 of 1,666 (8.6%) |
| CollecTRI | 32 | 1,119 | 76 | 62 | 190 of 1,666 (11.4%) |
| knockout data, transcription factors | 14 | 113 | 14 | 7 | 152 of 1,666 (9.1%) |
| both | 32 | 1,159 | 79 | 64 | 192 of 1,666 (11.5%) |

- **It finds more.** CollecTRI regulators add 47 of the drugs CTD records, a third
  more than the routes already had. Knockout data adds 9, in 14 of the 36 runs.
- **Its direction rate is not separated from chance.** Among the new drugs that
  CTD has an opinion on, 55% (CollecTRI) and 67% (knockout, 21 drugs) had the
  direction the goal asked for. An approved drug picked at random has 48%, and
  the routes the runs already had reach 61%. Resampling the 36 genes 5,000 times
  puts the CollecTRI gain over chance at -1% to +14% (95% interval) and the
  knockout gain at -16% to +49%. Both intervals include zero. The routes the runs
  already had do separate from chance (+4% to +20%).
- **The two-experiment rule for knockout data was chosen before the check, and
  the check does not favour a different one.** Asking for one experiment finds
  more (162 of 1,666 against 152) at 64% right direction (-1% to +32% over
  chance). Asking for three finds fewer (148) at 75% on 12 drugs, with an
  interval from -29% to +58%. None of the three is separated from chance.
- **It does not rank better.** The mean AUC for the right-direction test stays
  at 0.51 with the new routes added. They widen the list and do not sort it.
- **The CollecTRI gain is not independent.** CollecTRI and CTD draw on the same
  literature, so a drug recorded in both was probably found through the same
  papers. The knockout rows use measured experiments and are the fairer test.
  Counting 9 more drugs is the part of the gain that does not depend on shared
  curation.
- **Most new drugs are not in CTD for the gene.** 981 of 1,119 are neither
  wanted nor opposite there. That says CTD has no record, not that the drug does
  nothing.

For these reasons the routes are tier 4, and the route is off unless
`--regulator-route` is given.

The classification describes each drug's targets. It does not prove why a route
was missed, and an early version of this script counted failed ChEMBL requests
as "no target" and put 40% there. It now retries and keeps failures out.

## Rare disease genes whose drugs are known

Three genes where an existing drug is known to act through a factor, raised with all
the options except docking: ABCD2 (X-linked adrenoleukodystrophy, where ABCD2
compensates for ABCD1), LPL (lipoprotein lipase deficiency) and FXN (Friedreich
ataxia). Where a drug that should be there was not, `scripts/drug_trace.py GENE up DRUG`
follows it through the chain and stops at the first missing link: in the library, on the
Orange Book, a mechanism in ChEMBL, a recorded regulator of the gene, the right
direction, within reach upstream, and in CTD.

| Gene | Expected drug | Result |
|---|---|---|
| ABCD2 up | levothyroxine and liothyronine (THRA) | found, ranks 1 and 3 of 66 approved drugs |
| ABCD2 up | fenofibrate, gemfibrozil (PPARA) | found, ranks 18 to 20, through known regulators |
| ABCD2 up | valproic acid | found, rank 11, chromatin-level |
| LPL up | fenofibrate, gemfibrozil, bexarotene | found, ranks 7, 14 and 5 |
| FXN up | pioglitazone, rosiglitazone (PPARG) | found, ranks 19 to 20, through known regulators |
| FXN up | valproic acid and vorinostat | found, ranks 8 and 1 |
| FXN up | givinostat | **missing at first**, rank 12 after the fix below |
| FXN up | omaveloxolone | **still missing**, see below |

The scan verdict was "degraded" in all three: no signal for ABCD2 (9 sites against 11.6
by chance) and FXN weak (9 against 4.7), enriched for LPL (31 against 12.8). Bezafibrate
and pemafibrate are not FDA approved and are hidden by the filter, as intended.

Where the misses came from, and what was done:

| Cause | Example | Change |
|---|---|---|
| The approved-drug index dropped three ChEMBL target types that name human proteins: complex groups, fusion proteins and selectivity groups (120, 19 and 9 human mechanisms of approved drugs) | hydroxyurea at ribonucleotide reductase, imatinib at BCR and ABL1 | Included, listed at every member, and marked on the route |
| The curated chromatin list was written by hand | givinostat, approved in 2024 for Duchenne | Extended from ChEMBL's mechanisms at the HDACs, DNMTs and EZH2. It adds one drug today |
| The gene's own regulators are few | a factor two steps away | `--regulator-steps 2`, below |
| A drug with a measured activity at a receptor and no curated mechanism | rifampin at PXR: its only mechanism is the bacterial polymerase | `--regulator-measured`, below |
| No recorded link between the factor and the gene | omaveloxolone acts through NRF2, and no table links NRF2 to FXN | Not fixed. None of CollecTRI, TRRUST, the knockout data or the ChIP occupancy table names it, and nothing is within three signed steps. ReMap has one NRF2 peak among the 594 factors with peaks over the FXN window, so occupancy alone is too unspecific to add |

**Two steps.** `--regulator-steps 2` also lists drugs at the factors that regulate the
gene's recorded regulators, in CollecTRI, with the two signs multiplied. A factor whose
two routes disagree on the sign is left out. Against CTD on the 36 saved runs:

| | Wanted drugs reached (of 1,666) | Right-direction share | Gain over chance, 95% |
|---|---|---|---|
| Routes the runs already had | 143 | 61% | +4% to +20% |
| Known regulators, one step | 190 | 55% | -1% to +14% |
| Known regulators, two steps | 216 | 59% | +2% to +18% |
| Two steps plus measured activities | 378 | 58% | +4% to +13% |

The two-step interval excludes zero where the one-step interval did not, and the drugs
it adds (41 wanted, 20 opposite) have a 67% right-direction share. This was chosen
after looking at the same 36 genes, so it is not a held-out result. The launcher uses
two steps and the command line defaults to one.

**Measured activity.** `--regulator-measured` asks ChEMBL's activity records for a
regulator that the mechanism table gives no drug for, and lists approved drugs with a
pChEMBL value of 5 or better, the ten best. It finds rifampin at PXR on its own: PXR has
110 approved drugs with such a record and rifampin ranks fourth (6.7). The same list
puts droxidopa first (8.4) and includes ceftazidime, so the records are noise as well as
signal, and they give no direction, so every drug is listed with direction "unclear".
Against CTD on the 36 saved runs it reached 378 of the 1,666 wanted drugs (22.7%),
against 216 (13.0%) with two steps alone and 143 (8.6%) for the routes the runs already
had. It adds 4,748 - 1,849 = 2,899 new drugs to the 1,849 from two steps, of which 199
are recorded as raising the gene the wanted way and 148 the other way, a 57% share
against 48% for a random approved drug. The mean AUC for the reported-or-not test rises
from 0.553 to 0.591. The drugs have no direction, so this says that drugs with a
measured activity at a regulator are somewhat more likely than chance to be among those
CTD records as moving the gene the wanted way, and says nothing about which way a
particular drug acts.

The cost is one lookup per regulator with no curated drug, roughly 250 per gene on a
first run (the check took more than an hour for 700 distinct regulators), cached after
that, and it makes the list much longer: about 80 more drugs per gene on average. It is off
by default and the launcher does not pass it. Like the two-step option it was tried on
the same 36 genes it was scored on, so it is not a held-out result.

## Wider promoter windows, and the fate numbers after the fixes

**Five genes, four windows.** FXN, ABCD2, LPL, CYP24A1 and BUB1B were scanned at four
window sizes with the corrected chance estimate (20 shuffled copies). Sites above chance
rise with the window for most genes (ABCD2 0, 19, 20, 83; CYP24A1 2, 12, 48, 68; BUB1B 9,
13, 40, 75 at -1000/+200, -2000/+500, -5000/+1000 and -10000/+1000), but the enrichment
ratio stays at 1.4 to 2.1 and the chance share at 40% to 55%. A wider window finds more real
sites and more chance ones in about the same proportion.

**All 36 genes, two windows.** Each gene was run at -1000/+200 and at -5000/+1000 with
the options that drive drug recovery (scan, ChIP, factor drugs, upstream network, two-step
known regulators, LINCS) and none of the slow extras, then scored against CTD with the
cleaner key (1,489 wanted drugs). `scripts/window_comparison.py` runs and scores it.

| Arm | Wanted drugs with a route | Mean AUC, right direction | Mean AUC, reported or not |
|---|---|---|---|
| -1000/+200 | 214 (14.4%) | 0.531 | 0.572 |
| -5000/+1000 | 241 (16.2%) | 0.536 | 0.571 |
| Adaptive: the narrow window, the wide one where the narrow scan shows no signal | 224 (15.0%) | 0.534 | 0.573 |

The wide window reaches 27 more wanted drugs and does not move the AUC. By scan verdict,
19 of the 36 genes show no signal at -1000/+200, 11 weak and 6 enriched; at -5000/+1000
it is 18, 16 and 2. The adaptive rule switched 7 genes (ABCB1, CDKN1A, CYP24A1, CYP3A4,
GCLC, HMGCR and MMP9). That rule is on by default (`--no-widen` turns it off). The
evidence for it is a small gain in drugs reached and an honest scan verdict, not a better
ranking.

**The fate funnel, before and after.** The same funnel (`scripts/positive_fate.py`) on the
new runs at the default window:

| Stage | Before: 1,677 drugs | After: 1,494 drugs |
|---|---|---|
| found | 144 (8.6%) | 215 (14.4%) |
| acts on a found factor | 3 (0.2%) | 0 |
| in reach upstream | 598 (35.7%) | 445 (29.8%) |
| factor not found | 1 (0.1%) | 7 (0.5%) |
| too far upstream | 291 (17.4%) | 278 (18.6%) |
| not in the network | 265 (15.8%) | 227 (15.2%) |
| no target recorded | 375 (22.4%) | 322 (21.6%) |

Found rose from 8.6% to 14.4%. Part of that is the cleaner CTD key: its 1,489 wanted drugs
are fewer than the old 1,666, and the old routes reach 146 of them (9.8%). The rest, 9.8% to
14.4%, comes from the changes made to the routes. At the wide window found is 241 (16.1%).
The runs differ in two ways besides the fixes: the new ones used the lean option set, and
every fix was tried on these same genes, so it is not a held-out improvement. The share
with no target recorded fell by under a point. Many of those drugs act on bacterial,
viral or non-protein targets, or have no recorded mechanism, which no route here can reach.

## The motif scan's chance estimate

The reports say what share of the motif sites called would be expected by chance.
That number used to come from a formula: matrices times two strands times scannable
bases times the p-value. Checking it against shuffled copies of a real promoter showed
it was wrong, and showed two faults in the scan that made it wrong.

The window was BUB1B, 1,200 bases, 739 scannable after repeat masking, 879 JASPAR
2024 matrices, p of 1e-5.

| | Old scan | Corrected scan |
|---|---|---|
| Matrices scanned | 879 | 625 (254 could not reach p = 1e-5) |
| Sites called | 37 | 17 |
| Expected by the formula | 13 | 9 |
| Found in shuffled copies | 27 +/- 10 (40 copies) | 7.7 (30 copies) |
| Share expected by chance | 35% by formula, about 70% by shuffles | 46%, range 24% to 65% |
| Best 5 and best 10 sites | not given | 31% and 39% |

- **Some matrices could not reach the p-value.** A short or weakly specific matrix
  can have a best possible match that is itself more common than 1e-5 in the window.
  The old scan then used the best score as its threshold, so every consensus match
  was called, at a rate of 4.5e-5 for NKX2-5 in this window. 254 of 879 matrices were
  in this state. They are now left out, counted, and named in the report, and every
  site called has a p-value at or below the one requested (the highest is 9.9e-6).
- **The strands shared one threshold.** The reverse strand scores the reverse
  complement of the matrix, whose scores are distributed differently when A and T
  are not equal. The thresholds differ by more than 0.3 bits for 354 of 570 reachable
  matrices in this window. Each strand now has its own.
- **The formula undercounted by half.** With both faults the shuffled copies gave 27
  sites where the formula gave 13. With them fixed the two agree (7.7 and 9) and the
  formula is a fair cross-check, slightly high because scores are discrete.
- **The estimate is now measured.** The run scans shuffled copies that keep the repeat
  mask and every two-base count, collapses overlapping hits the same way, and uses
  their mean. The range is the 10th to 90th percentile of the copies. A shuffle that
  moved repeat bases into the scanned region gave a different answer, so each run of
  unmasked bases is shuffled on its own. `--null-shuffles 0` returns to the formula.
- **What it does not say.** A site is a sequence match. Which sites are real is a
  question for ChIP, and the shuffles give no way to tell them apart. A shuffle that
  keeps two-base counts is a fair null for the base composition and a poor one for
  anything longer, such as a repeat family the mask missed.

Compatibility: sites that were called only through the faults are gone, so a run now
lists fewer sites and fewer factors, and the factor ranking can differ from an earlier
run of the same gene. `scan_quality` in `results.json` keeps its fields and adds
`expected_by_formula`, `null`, `false_discovery_range`, `top_site_fdp`,
`motifs_loaded`, `motifs_scanned` and `motifs_unreachable`. `expected_by_chance` is now
the shuffled mean when shuffles were run. Each hit has a `pvalue`. The verdict
thresholds are unchanged, so a window can now be called "weak" or "no signal" where
the formula called it enriched, and the reverse.

## What CTD's expression record contains

CTD's answer for one gene is a list of sentences, and its `InteractionActions` field
lists every action in a sentence, not only the ones about this gene. Reading that field
turned any chemical with an `increases^expression` action into a chemical that raises
the gene. `python scripts/ctd_check.py <results.json files>` sorts every row that has an
expression action by what its sentence says, across the 36 saved genes
(100,979 rows):

| What the sentence says | Share of rows |
|---|---|
| The chemical changes this gene's mRNA | 39% |
| Expression of this gene's protein | 26% |
| The chemical modulates another's effect ("inhibits the reaction [X results in increased expression of G]") | 24% |
| Another gene's expression | 8% |
| The effect comes through another gene | 3% |
| A shape that is not read | 1% |

A chemical that blocks an inducer was read as raising the gene, and one that raises
another gene was read as raising this one. The fetch now reads the sentence. It counts a
chemical that is the subject of an expression clause about this gene's mRNA, or about
the gene with no form stated. A chemical that modulates another's effect is kept as a
known link with no direction. Protein expression and indirect effects are left out.

Scored against the saved routes, the cleaner key has fewer approved drugs recorded as
raising a gene (1,489 against 1,666) and the routes reach slightly more of them (146
against 143, so 9.8% against 8.6%). The mean AUC for the reported-or-not test moves from
0.535 to 0.543. Neither change is large, so the cleaner key does not rescue the
ranking, and the 8.6% figure is not mostly an artefact of the key.

Compatibility: `improvements.reported_links` in a new run has fewer chemicals and
different directions for some, so counts such as "77 move it the way you want" fall.
Saved results keep what they recorded. `discovery.ctd_expression_links(client, gene,
strict=False)` returns the earlier reading.

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
- **The 36 saved runs predate two corrections.** They were made before the
  motif scan stopped calling matrices that cannot reach the requested p-value and
  before CTD's sentences were read, so the figures in this section describe the
  earlier code. The 36 genes were rerun with the corrected scan, the cleaner key
  and the later route fixes; see "Wider promoter windows, and the fate numbers after
  the fixes" below, where found rose from 8.6% to 14.4%.

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
