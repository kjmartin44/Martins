# promoter-tf-drugs

Given a gene, this finds its promoter, scans for transcription-factor binding
sites, works out which of those factors are **activators** (and which are
ligand-activated **nuclear receptors**), and looks up **FDA/EMA-approved drugs**
whose pharmacology pushes those factors in the direction that would raise
output from your target gene.

Requires Python 3.10+ and numpy. Nothing else.

```bash
# runs fully offline on bundled fixtures - start here to see the output format
python -m promoter_tf_drugs --demo --gene DEMOGENE --out results/demo

# the other direction: drugs predicted to LOWER the gene
python -m promoter_tf_drugs --demo --gene DEMOGENE --direction down --out results/demo_down

# real run
python -m promoter_tf_drugs --gene PAH --upstream 2000 --downstream 500 \
    --pfm-file data/jaspar2024.txt \
    --chip-bed data/remap2022.bed.gz \
    --out results/PAH -v

python -m pytest tests/
```

Outputs: `report.md`, `results.json`, `tf_summary.tsv`, `sites.bed`
(load in IGV/UCSC alongside the gene), `promoter.fa`.
The real run above uses an optional ReMap file. On Windows, run
`scripts/download_remap.bat` once before using `--chip-bed`.

## Repository layout

- `promoter_tf_drugs/` contains the Python package.
- `tests/` contains the test suite.
- `data/` contains the included JASPAR, TRRUST, and HPA reference files.
  `scripts/download_remap.bat` downloads the optional ReMap file there.
- `scripts/` contains Windows launchers. Run them from any folder; they switch
  to the repository root before starting the pipeline.
- `examples/CCR5_down/` contains a saved analysis. New runs go in `results/`,
  which Git ignores.

---

## Cell-type reach (`--cell-type-check`, `--tissue cns`)

Three things must overlap, and without this flag the pipeline checks none of
them:

```
where the drug goes  ∩  where the gene is expressed  ∩  where the TF is
```

Any empty intersection kills the route regardless of motif quality, potency or
grade. **The third term is the one usually skipped**: a beautifully scored
hepatocyte transcription factor is irrelevant if you are treating neurons.
`--cell-type-check` pulls tissue expression for your gene and each top TF and
heavily down-weights factors with no shared tissue.

### CNS penetration: MPO, not a rule of thumb

`--tissue cns` now uses the CNS Multiparameter Optimisation score (Wager et
al., *ACS Chem Neurosci* 2010) rather than a MW/TPSA/logP cutoff: six
properties, each through a desirability function, summed to 0–6, with ≥4
conventionally favourable. Note the TPSA term is a **hump** — too apolar is
penalised as well as too polar, which a one-sided cutoff misses entirely.

Sanity check on known drugs: caffeine 4.0/4, diazepam 3.56/4, plicamycin 1.0/4.

ClogD and pKa need tooling beyond RDKit, so scores are usually computed over 4
of 6 parameters and **reported as such** rather than silently rescaled. A
worked limitation, encoded as a test: atenolol scores 3.17/4, reading as
favourable, yet it is a textbook poor penetrant — the missing pKa term (~9.6,
strongly basic) is exactly what would catch it. Pass `pka` and `logd` when you
have them, and treat a borderline pass on a basic amine with suspicion.

### The bigger reason CNS drugs fail

Physicochemistry predicts *passive diffusion*. The commonest reason a suitable
molecule never reaches brain is **active efflux by P-glycoprotein** at the
endothelium — a compound can have an excellent MPO and be pumped straight back
out. Known substrates override a good score entirely. Substrate prediction is
genuinely poor, so a negative flag is **not** reassurance, and the report says
so rather than implying a clean result.

### Other compartments, and one hard stop

Be honest about limits: **BBB is the only barrier with a decent predictive
model.** For liver, kidney, muscle or tumour, distribution depends on
transporters, protein binding and perfusion in ways physicochemistry does not
capture — the module says so instead of inventing a number.

One absolute constraint that is easy to miss: **antibodies and other biologics
cannot reach intracellular targets.** Every transcription factor is
intracellular. An antibody against a TF is not a weak option, it is not an
option, and it is rejected rather than scored.

---

## Screen design (`--design-screen`)

Look at what limits the outcome chain: almost always *gene responsiveness in
your cell type*, the link this package cannot see. A screen measures it instead
of assuming it. **Everything else here prioritises this experiment; none of it
replaces it.**

A CRISPRa screen across the ~1,600 human TFs answers "which factor actually
drives my promoter" without a prior. An arrayed approved-compound library
answers the drug question the same way, and is exactly what surfaces the
"anti-parasitic that massively activates your promoter" case — because it
assumes nothing about mechanism.

Produced: a reporter plan with the controls that make a hit interpretable,
CRISPRa sgRNAs in the productive window, a cherry-pick compound list ordered by
the pipeline's grades, and realistic hit rates.

Three things in the output that matter more than the guide sequences:

**The scrambled-motif control.** A mutant reporter with your top predicted
element destroyed. If a drug induces wild-type but not mutant, the effect runs
through that element. If it induces both, your prediction was wrong. This one
control converts a screen hit into a mechanism.

**The counter-screen.** Run the same plate against a constitutive promoter.
Anything moving both is affecting transcription or translation globally rather
than your promoter — and that is the majority of primary hits.

**What a negative result means.** If nothing in a genome-wide CRISPRa screen
activates your reporter, the promoter fragment is probably not sufficient and
the regulation lives in an enhancer outside your window. Widen it or go to the
endogenous locus. That is a real result, not a failed experiment.

sgRNAs are placed in the -400 to -50 window where dCas9-VP64/SAM activation
actually works; guides containing TTTT are flagged because U6 will terminate
them. No off-target analysis is done — run candidates through CRISPOR or
CHOPCHOP before ordering.

### AlphaFold docking, and the trap specific to TFs

`docking.py` can now use AlphaFold models, gated on pLDDT. Two reasons
predicted structures dock badly, both of which apply here:

**Apo side chains.** Models are unliganded; side chains have not reorganised
around a ligand and cryptic pockets are closed. Even a 95-pLDDT model carries
this warning.

**Low pLDDT is a disorder prediction, not a pocket.** This is the
TF-specific trap. Activation domains are frequently intrinsically disordered;
AlphaFold reports that as low pLDDT, and a docking program sees an open greasy
cavity and scores compounds into it with excellent affinities. Every such hit
is an artefact. Models below 50% confidently-folded residues are **refused**,
partially-ordered models restrict docking to folded regions, and poses
contacting low-pLDDT residues are rejected.

Use it to triage which pockets deserve an experiment, never as evidence of
binding.

---

## Measured expression changes: LINCS L1000 (`--lincs`)

**The only empirical evidence in this package.** Everything else infers: a
motif resembles a site, a factor is annotated as an activator, a drug is
annotated as an agonist. L1000 measured what happened to the transcriptome when
cells were dosed. Where it agrees with a predicted route, you have a real
candidate rather than a hypothesis.

```
python -m promoter_tf_drugs --gene YOURGENE --lincs --out results\YOURGENE
```

Two uses. **Direct**: which approved drugs moved your gene? **Connectivity**
(`--signature-file`, exported from CLUE): which drugs produce a whole-signature
resembling a known activator of your gene, *regardless of annotated target*?
The second is what CMap was designed for and how you find a compound whose
official mechanism has nothing to do with your pathway but whose cellular
effect does.

### Four things that decide whether a hit is real

**Landmark vs inferred.** L1000 physically measures ~978 landmark transcripts
and *computationally infers* the other ~11,350. If your gene is not a landmark,
its value is a model prediction from the landmarks — not a measurement. This is
the commonest misuse of the dataset and it is almost never stated. Hits are
weighted 1.0 vs 0.4 and the report prints `measured` or `INFERRED` per row.
Unknown status is treated as inferred, deliberately.

**Cell line.** Signatures come from a fixed panel (MCF7, PC3, A375, HA1E, VCAP,
HEPG2...). If your gene is not expressed in the line that produced the
signature, there is nothing to move and absence of signal means nothing. The
cell line is printed on every row for that reason.

**Non-specific stress dominates.** A large share of strong L1000 signatures are
a version of "this cell is dying" — HSP induction, ATF4/CHOP, cell-cycle
collapse. Cytotoxic compounds therefore connect to almost everything, and a
naive connectivity search returns the same two dozen drugs for any query. Those
are flagged and demoted rather than reported as hits.

**Dose and time.** 6 h and 24 h signatures of one compound differ
substantially, and high doses drift toward the stress response above.

---

## Human genetic evidence (`--genetics`)

**A correction first.** GWAS and PheWAS do not record whether a gene's
expression changes in patients taking a drug — that relation is not in those
datasets. GWAS links germline variants to traits; PheWAS scans many traits
against one variant. Neither observes drug exposure and expression together.
For drug-induced expression change you want LINCS, above.

What genetics does better than the thing being asked for:

**Genetic support for the TF → gene link (trans-eQTL).** If variants near a
transcription factor associate with your target gene's expression in humans,
that is in-vivo causal evidence that the factor regulates your gene — real
tissue, real people, exposure randomised by meiosis. A motif says a protein
*could* bind; a trans-eQTL says variation in it *does* change your gene's
output. This is the strongest regulatory evidence available anywhere in the
pipeline. Check eQTLGen (blood, n>30,000) and the eQTL Catalogue; the report
tells you exactly what to look up.

**What happens if you succeed (cis-eQTL → phenotype).** Drug-target Mendelian
randomisation run backwards: variants raising your gene's expression are
lifelong natural experiments for the intervention you are designing, including
its on-target adverse effects. This is the reasoning that anticipated PCSK9
inhibition years before the trials read out.

Caveats that decide applicability: eQTLs are **tissue-specific** (a liver eQTL
says nothing about brain); variants act **lifelong** while drugs act for
months, so developmental compensation confounds the comparison in a direction
you cannot predict; **horizontal pleiotropy** can break the instrument; and
**absence is weak** — no trans-eQTL usually means insufficient power, not
absence of regulation.

---

## Cell-cycle compatibility (`--proliferation`)

```
--proliferation cycling | post-mitotic | mixed | unknown
```

**This is a different clock from what a PAM preserves.** A PAM keeps the
endogenous ligand's signalling dynamics — pulsatile, circadian,
tissue-restricted — on a timescale of seconds to hours. Cell-cycle dependence
asks a blunter question: *is this factor active in the cells you are treating
at all?* A PAM at a G1/S-restricted factor still does nothing in a cell sitting
in G2.

**Post-mitotic tissue: routes are removed, not demoted.** In neurons,
cardiomyocytes and other terminally differentiated cells, E2F-dependent
promoters are held down by the DREAM complex (E2F4, DP, p130, MuvB). Those
targets are not quiescent-but-available — they are *actively repressed*. A
route through E2F1 or FOXM1 there does not work at any dose, so ranking it
lower would still invite someone to try it. It is deleted, and the report says
which ones and why. If you are working on a neurological or cardiac
indication, set this flag.

The exception is **E2F4**, the repressive arm of DREAM, which persists in
quiescence and is kept.

**Cycling cells: duty cycle dilutes the measurement.** E2F1 is active across
roughly the late-G1-to-early-S window, about 25% of the cycle. In an
asynchronous population only that fraction of cells can respond at any moment,
so the population-level effect is scaled accordingly. This matters practically:
**per-cell efficacy can be much higher than a bulk RT-qPCR suggests.** Before
concluding a drug does not work, synchronise or use a phase-resolved readout —
the report says so for any factor below 25% duty.

**Mitotic bookmarking** is flagged where relevant. Most TFs are evicted from
chromatin during mitosis; a minority (FOXA1, GATA1, ESRRB, CTCF, RUNX2) stay
bound and reseed transcription in the daughters, making them less
phase-sensitive.

Unannotated factors are assumed cycle-independent — the common case, and the
one that avoids an unearned penalty. With `--proliferation unknown` (the
default) restricted factors are flagged but not scored, so you are told the
issue exists without the tool guessing at your biology.

---

## Multi-hop upstream search (`--upstream`)

The rest of the package is one hop: a drug that binds your TF. That ceiling is
severe — most of the ~1,600 human TFs have no approved ligand at all. CREB1 has
none, but PKA does, adenylyl cyclase does, and the GPCRs above them do. This
walks upstream through signed causal edges (OmniPath/SIGNOR) to find druggable
kinases, phosphatases and receptors, and reports the path that explains each
one.

It runs only for TFs with **no** direct pharmacology, since that is the gap it
exists to fill.

```
python -m promoter_tf_drugs --gene YOURGENE --pfm-file data/jaspar2024.txt \
    --upstream --upstream-depth 3 --out results\YOURGENE
```

### Four ways this goes wrong, all handled

**Signs must multiply.** Inhibiting an inhibitor is activation. A path's net
direction is the product of its edge signs, so a 4-hop path with two inhibitory
edges is net activating. Unsigned networks — plain PPI, STRING-style resources —
cannot support this analysis at all, which is why signed causal edges are
required rather than preferred. Tested in both directions, including recovery
of double-negative routes.

**Specificity collapses with distance**, and this works *against* the
specificity axis the rest of the package optimises. A node three hops up
reaches your gene through one branch and everything else through the others.
Hop count carries an explicit penalty, fan-out is printed per path, and paths
at ≥3 hops get a warning. Distance buys druggability and spends selectivity;
the trade is shown rather than hidden in a rank.

**Random walks find hubs, not biology.** Node degree in curated networks tracks
how much a protein has been *studied*. AKT1, TP53, EGFR and MAPK1 have enormous
edge counts, so an unnormalised walk lands on them from any seed. The RWR here
is degree-normalised and hub nodes are flagged. If your top hit is AKT1 for
every gene you try, that is the artefact, not a finding.

**Confidence decays multiplicatively** — 0.75 per hop, so 3 hops is 0.42 before
the drug is considered.

### Direction inverts twice

To raise a gene through an *activating* TF you must raise the TF; through a
*repressor* you must lower it; and a `--direction down` goal flips it again.
Then the path's own net sign decides whether you activate or inhibit the
upstream node. All four combinations are tested, and the report prints "act on
it: activate / inhibit" so you do not have to work it out.

Every ranking comes with its path printed (`RECEPTOR2 —| PHOSPHATASE —| TF`).
A rank without a mechanism is not actionable — you should be able to read the
chain and judge whether it is plausible in your cell type, which is a judgement
no network can make for you.

---

## Changing the natural ligand instead of the receptor

A nuclear receptor's activity is set by how much of its endogenous ligand is
around, and that concentration is controlled by enzymes. Those enzymes are
frequently far more druggable than the receptor — and several of the most
successful drugs in the pharmacopoeia work this way without binding the
receptor at all:

| drug | acts on | result |
|---|---|---|
| anastrozole | aromatase | less estradiol → **less** ER activity |
| finasteride | 5α-reductase | less DHT → **less** AR activity |
| metyrapone | CYP11B1 | less cortisol → **less** GR activity |
| carbenoxolone | HSD11B2 | cortisol not inactivated → **more** GR activity |
| levothyroxine | — | supplies T4 → **more** TR activity |

None of these appear in any other module, because every other one looks for
ligands *of the receptor*. This looks one step sideways.

Four strategies and the direction each produces:

| strategy | receptor activity |
|---|---|
| supplement ligand | up |
| inhibit catabolism | up |
| inhibit synthesis | down |
| sequester ligand | down |

**Two reasons to try it.** The enzyme is often druggable where the receptor is
not. And some routes are **tissue-selective in a way the receptor cannot be**:
HSD11B2 converts cortisol to inert cortisone only where it is expressed, so
inhibiting it raises cortisol *locally*. Tissue-selective routes are ranked
first and score higher on specificity for exactly this reason.

**Two reasons it fails.** The ligand acts on every receptor it binds,
everywhere — aromatase inhibitors cause bone loss because estradiol does far
more than drive your gene, and that systemic phenotype is usually why the drug
is hard to use. And endocrine axes are homeostatic, so feedback partly opposes
you: the effect is smaller and slower than the enzyme inhibition suggests, and
can rebound on withdrawal.

For supplementation specifically, check the axis has headroom. If the receptor
is already saturated in your tissue, more ligand does nothing.

---

## Orthosteric vs allosteric (PAM/NAM)

Every run compares the mechanism classes available at your top receptors,
because a positive allosteric modulator is not just another way to activate a
target — it trades on four axes and they do not all point the same way.

| | orthosteric agonist | PAM |
|---|---|---|
| max effect | full | **bounded (ceiling)** |
| subtype selectivity | typically lower | **typically higher** |
| preserves signalling dynamics | no | **yes** |
| needs endogenous ligand | no | **YES** |
| desensitisation risk | higher | lower |

**Why PAMs are selective.** The orthosteric pocket is under selection to bind
the same endogenous ligand across a receptor family, so it is conserved and
hard to discriminate. Allosteric sites carry no such constraint and diverge
freely — which is why subtype-selective PAMs exist for targets where decades of
orthosteric chemistry produced nothing selective. For a pipeline whose main
failure mode is "this drug moves far more than your gene", that is the
strongest argument for looking at them. The specificity axis gets a 1.30x
multiplier for allosteric mechanisms, justified by that conservation argument
rather than by fitting.

**Why they preserve dynamics.** An agonist clamps the receptor on wherever the
drug reaches. A PAM amplifies whatever the endogenous ligand is already doing,
so tissue-restricted, pulsatile and circadian patterning survives. For a TF
under feedback or circadian control that can be the difference between a
physiological increase and a flat, dysregulated one.

**The catch, and it is a hard one.** No endogenous tone, no effect. If the
pathway is silent in your cell type, a PAM does nothing at all. **Verify tone
before pursuing any allosteric option** — that is a measurement in your system,
not something this package can predict. PAMs also keep their lower action
confidence (0.65) in the grade, so the trade is not free.

Reality check for nuclear receptors: approved PAMs cluster on GPCRs and ion
channels (benzodiazepines at GABA-A, cinacalcet at CaSR, ivacaftor at CFTR).
True clinical allosteric modulators of nuclear receptors are rare, though
allosteric sites are established — the PPARgamma alternate site, the RXR
tetramer interface, Nurr1. Expect most PAM hits to sit on upstream signalling
nodes rather than on the TF itself; they are labelled as such.

`--no-mechanisms` to skip.

---

## Running offline

Public APIs go down. Every stage has a local-file route except gene lookup:

| stage | offline flag | where to get the file |
|---|---|---|
| motifs | `--pfm-file` | JASPAR bulk PFM download |
| ChIP evidence | `--chip-bed` | ReMap 2022 / ChIP-Atlas BED |
| regulon breadth | `--regulon-file` | TRRUST v2 TSV |
| tissue expression | `--expression-file` | HPA `rna_tissue_consensus.tsv`, or a GTEx median-TPM matrix |
| LINCS signatures | `--lincs-file` | CLUE / iLINCS export, `drug,gene,score[,cell_line,is_landmark]` |
| everything | `--offline` | replays the on-disk cache only |

`--expression-file` accepts both shapes people actually have: HPA's long format
(Gene, Gene name, Tissue, nTPM) and a wide gene-by-tissue matrix like GTEx. It
takes **priority over the API**, so a broken Open Targets endpoint stops
mattering. Downloads: HPA at https://www.proteinatlas.org/about/download
("RNA consensus tissue gene data"), GTEx at
https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression.

`--lincs-file` applies the same landmark/inferred weighting as the API route:
if your table has an `is_landmark` column it is honoured, otherwise values are
scored as inferred, which is the conservative reading.

Only Ensembl (gene coordinates and sequence) has no offline substitute — but
its results cache, so a gene you have fetched once works offline afterwards.

---

## Expanding one compound (`explain`)

A run scatters evidence about a single drug across eight sections — the route
catalogue, the grade table and its components, the outcome chain, the mechanism
comparison, the ligand axis, the LINCS signature, the comparison against the
standard of care, and the screen list. To judge one candidate you would
otherwise have to find it in each.

```
python -m promoter_tf_drugs.explain results\PAH_up\results.json calcitriol
```

Partial names work; where several compounds match, all of them are reported
rather than the first being picked silently.

Output pulls together the route and its axes, the full grade breakdown with
the underlying pChEMBL and regulon evidence, the mechanism class and its
trade-offs, any ligand-axis route, the measured LINCS signature if there is
one, the outcome chain with its **limiting step and the experiment that would
narrow it**, and how the compound compares against the per-axis bar — including
whether it has a previously reported link to your gene.

---

## Narrowing the list (`shortlist`)

A full run takes tens of minutes; deciding what to order should take a second.
This filters an existing `results.json` with no network and no re-run, so you
can tighten and loosen criteria interactively.

```
python -m promoter_tf_drugs.shortlist results\PAH_up\results.json --preset strict
```

Or double-click `scripts/shortlist.bat`.

| preset | the question it answers |
|---|---|
| `strict` | what would I actually order next week? |
| `mechanism` | what do I understand well enough to design an experiment around? |
| `empirical` | what has actually been measured? |
| `specific` | what will not move half the transcriptome? |
| `novel` | what has nobody tried on this gene? |

Presets are starting points; individual flags override them
(`--max-tier`, `--min-grade`, `--max-steps`, `--min-specificity`,
`--max-regulon`, `--require-chip`, `--exclude-conflicts`, `--only-novel`).

**An empty `strict` result is informative, not a failure.** It says no approved
drug clears the bar for your gene on current evidence — which is worth knowing
before you spend months on one. Loosen one criterion at a time to find which is
binding; the output tells you where to start, and `--require-chip` is usually
it.

Two things the filters get right that a naive version would not. **Missing
grades pass**: letter grades only exist for TF-mediated routes, so a direct
siRNA or approved gene therapy has none, and dropping a tier-1 clinical agent
for lacking an inapplicable field is backwards. And **`strict` explicitly
excludes genome-wide agents**, because chromatin drugs have no grade either and
would otherwise slip past the grade filter on a technicality despite
specificity around 0.05.

---

## Every route in one table

Each module answers its own question and prints its own section. That left the
obvious question unanswered: *what are ALL the ways to move this gene, and
which is best?* The report now opens each promoter with a pooled ranking of
every route found, across ten mechanism categories:

| | |
|---|---|
| direct: gene or transcript | ASO, siRNA, splice modifiers, gene therapy |
| direct: gene product | drugs on the protein itself |
| TF: orthosteric ligand | classic agonists/antagonists |
| TF: allosteric modulator | PAM/NAM |
| TF: stabiliser or glue | conformational locks, molecular glues |
| TF: degradation blocker | NRF2/HIF-style stabilisation |
| upstream signalling | kinases and receptors N hops up |
| empirical signature | LINCS hits, mechanism unknown |
| chromatin-level | HDAC/DNMT/EZH2 inhibitors |
| virtual screen hypothesis | predicted, never tested |

### Sorted by evidence tier first, score second

Collapsing these into one number would be the wrong fix. "Approved, documented
in humans" and "a docking program liked this pose" differ by orders of
magnitude in how much they should move you, and averaging destroys exactly the
information you need. So:

| tier | meaning |
|---|---|
| 1 | clinical — approved for this gene, documented in humans |
| 2 | measured — observed experimentally (LINCS landmark, ChIP) |
| 3 | curated — a human curator annotated the mechanism |
| 4 | predicted — a model says so |
| 5 | hypothesis — generated, untested |

**A tier-1 route with a mediocre score beats a tier-5 route with a great one,
every time.** Anything at tier 4 or 5 is an experiment you have not done. Seeing
them side by side shows how far you are reaching, rather than hiding it.

Two things the table makes visible that the separate sections did not.
**Chromatin drugs are scored as remote and unspecific** (specificity 0.05,
directness 3): they will move almost any promoter, which makes them an
excellent positive control that your reporter *can* be induced, and almost
never an acceptable therapeutic route. And **absent categories are named**, so
an empty "upstream signalling" row tells you either that no route exists or
that you forgot `--upstream-network`.

---

## Finding something BETTER than what exists (the main goal)

Ranking tells you what to read first. "Stronger activation, higher specificity,
or superior efficiency" is a **comparative** claim, which needs three things a
ranking does not give you. All three are implemented:

**1. A bar, not a single incumbent.** The comparison is against the **per-axis
best across every established option** — the strongest value anyone has
achieved on strength, on specificity and on efficiency, whichever drug holds
each. Picking one incumbent by summing its axes would contradict the Pareto
design: a challenger could "beat the incumbent" on specificity while a
different established drug is more specific still.

If nothing exists, the report says so and ranks candidates on their own merits
rather than quietly comparing against zero.

**Nothing is withheld.** Every candidate is accounted for in the report:
already-linked, novel and non-dominated, novel but beaten on every axis
(these are named), and novel with no advantage. The full list is always in
`results.json` under `improvements.all_candidates`.

**2. Three axes, kept separate.** A weighted sum lets a candidate that is
slightly stronger and drastically less specific come out ahead — exactly the
trade this package exists to prevent. So candidates are compared by **Pareto
dominance**: interesting means beating the incumbent on at least one axis
*without* being meaningfully worse on another, and the winning axis is named.

| axis | what it means |
|---|---|
| strength | potency x unambiguous pharmacology x how well the TF drives this promoter |
| specificity | protein selectivity **x** regulon breadth, multiplied |
| efficiency | effect per unit dose, exposure and risk |

Trade-offs are reported *as* trade-offs: "better on efficiency; worse on
specificity" — never collapsed into a single verdict.

**3. A novelty test.** CTD curates chemical-gene interactions with explicit
"increases/decreases expression of" actions — the exact relation this pipeline
predicts. Candidates already in CTD for your gene are labelled known links
(still useful, as confirmation); the rest are flagged NEW. Matching is
deliberately generous, because presenting a known relationship as a discovery
is the one genuinely embarrassing failure mode, while over-calling a link as
known merely costs you a candidate.

### The margin, and why it matters

An incumbent's numbers are **measured**. A challenger's are **predicted** by the
models here. Comparing them directly is biased toward the challenger, and not
for the obvious reason: you are looking at that candidate *because* its
predicted score came out high, so the same noise that made you notice it has
already inflated it. Selection is not symmetric even when error is.

So a predicted challenger must clear a measured incumbent by 0.15 on an axis
before the difference counts. This is why some candidates are reported as
"within noise of the incumbent" rather than as improvements — and why an empty
result is informative rather than a failure: it says the obvious pharmacology
has been explored.

`--no-improvements` to skip.

---

## Predicted transcriptional outcome

```
P(gene moves) = P(engages) x P(direction) x P(activates here) x P(responsive)
```

Receptor-level agonism is **not** a transcriptional outcome, and conflating the
two is the most common error in this kind of analysis. Tamoxifen antagonises
ERalpha in breast and activates it in endometrium and bone — same drug, same
receptor, same pose. What differs is which coregulators the cell has: helix 12's
position creates a surface, and whether that surface recruits a coactivator
(LXXLL) or a corepressor (CoRNR) depends on their relative abundance in that
cell type.

So the pipeline reports the chain rather than a single number, and names the
weakest link plus the experiment that would narrow it. Typical output is
P ≈ 0.15, limited by "gene responsiveness in your cells" — the link nothing in
this package can see. Four honest links multiply to small numbers; a chain at
0.7 per link lands near 0.24. **Anything reporting 0.9 here has hidden a link.**

Ligands predicted between 0.35 and 0.65 are flagged **tissue-dependent**: that
is the SERM regime, where direction is a property of your cell type rather than
the molecule, and no structural model can settle it. A coactivator recruitment
assay can. `--no-outcome` to skip.

---

## Virtual screening (`--screen`)

Scores every approved drug with a structure against your receptor's classifier,
to find pairs nobody has tested.

**What this cannot do:** discover a mechanism nobody has characterised. Every
model here is trained on published measurements; a classifier cannot know about
a binding mode never observed. That is a property of the training data, not of
the method, and no amount of compute changes it.

**What it can do:** there are ~4,000 approved drugs and each has been assayed
against a tiny fraction of the proteome. A 1974 drug for an unrelated
indication may well engage your receptor — if nobody ran the assay, that pair
is in no database. Surfacing those is real repurposing.

The two are in tension: predictions are reliable in proportion to similarity to
the training set, and similar molecules are exactly the ones already tested. So
the interesting quadrant is also the least reliable. Every hit is labelled
rather than hiding that behind one score:

| quadrant | meaning |
|---|---|
| known ligand | already annotated here. A positive control, not a finding |
| confident prediction | close to known ligands. Most likely to replicate, least likely to be new |
| **NOVEL HYPOTHESIS** | dissimilar and predicted active. The repurposing quadrant, and an **extrapolation** — requires a binding assay before it means anything |
| no call | outside the applicability domain; score withheld |

Before reporting anything, `enrichment_check` asks whether known ligands rank
above the rest. **A screen that cannot recover its own training ligands is
declared unusable and its novel list suppressed.** This is the check people
skip, and without it a broken model still produces a confident-looking ranking.

Expect most NOVEL HYPOTHESIS rows to be wrong. Single-digit hit rates are
considered a success in virtual screening.

---

## Designing a direct modality (`--design-modality`)

`direct.py` *looks up* approved agents already targeting your gene, which for
most genes returns nothing. This **designs** one instead: siRNA and gapmer ASO
candidates from the transcript, CRISPRi guides at the promoter, and — the part
that actually matters — a deliverability verdict.

```
python -m promoter_tf_drugs --gene CCR5 --direction down --tissue cns \
    --design-modality --pfm-file data\jaspar2024.txt --out results\CCR5
```

### Delivery is the constraint, not design

Designing a potent siRNA is close to solved: the sequence rules are well
established and a competent design usually knocks down 80-95% in transfected
cells. Programs do not fail there. **They fail because the molecule never
reaches the cell type that matters.** So deliverability is scored first:

| compartment | status | route |
|---|---|---|
| liver | **solved** | GalNAc, subcutaneous — inclisiran, givosiran, lumasiran, vutrisiran |
| CNS | solved with a caveat | intrathecal — nusinersen, tofersen; means repeated LPs |
| eye | solved | intravitreal |
| lung, kidney, muscle, tumour | partial | chemistry- and compartment-dependent |
| **T cells, macrophages, microglia** | **NOT solved in vivo** | ex vivo editing of HSCs or T cells, then reinfusion |

That last row is the one to read carefully for immune-cell targets. There is no
routine systemic oligo delivery to T cells or macrophages, which is exactly why
CCR5 work went the ex vivo route rather than systemic knockdown — and ex vivo
editing is a transplant programme with conditioning, not a drug. When delivery
scores below 0.5 the report says plainly that **delivery, not design, is the
blocker**, and that a potent oligo for that target is a cell-culture reagent.

### Design rules applied

**siRNA**: GC 30-52%, no run of 4+ identical bases, thermodynamic asymmetry
favouring antisense loading into RISC, first ~75 nt skipped (ribosome-loaded),
immunostimulatory motifs flagged.

**Gapmer ASO** (5-10-5 LNA/DNA/LNA): GC 40-60%, CpG flagged as
immunostimulatory in phosphorothioate backbones, G-quadruplex and homopolymer
runs flagged. Works on nuclear RNA, so it reaches pre-mRNA and nuclear-retained
transcripts that siRNA cannot.

**CRISPRi** uses a different window from CRISPRa — −50 to +300 around the TSS
rather than −400 to −50 — because KRAB silencing is most effective at the TSS
while VP64 activation works upstream of it. Both windows are asserted distinct
in the tests.

**Off-target is not checked here.** Seed-region (guide positions 2-8)
complementarity drives most siRNA off-target silencing and is a
transcriptome-wide question. Run every candidate through BLAST against RefSeq
before synthesising anything; the checks here are local sanity only.

---

## Drugs that skip the TF route entirely

Every run first checks whether anything targets your gene **directly** — its
own protein, its transcript, or the gene itself. This section is printed at the
top of the report, because when such a drug exists it is almost always the
better starting point, and it is the first thing a reviewer will ask whether
you checked.

The TF route is indirect and lossy: predict a site, hope the factor binds it in
your cell type, hope it activates rather than represses, hope the drug moves it
enough, and accept that everything else that factor controls moves too. A
sequence-directed agent skips all of that.

Covered: approved siRNA, antisense and splice-switching oligos, gene
therapies, antibodies, and small molecules whose ChEMBL/DGIdb target is your
gene's own product. `--no-direct` turns it off.

**The upregulating members of this class are the interesting ones**, because
there is no general "siRNA for turning genes up". Nusinersen and risdiplam
raise SMN protein by correcting SMN2 splicing; the DMD exon-skippers restore a
partly functional dystrophin; the AAV therapies supply a working copy outright.

One distinction the databases will not make for you: an **activator raises the
activity of protein already present, not its abundance**. Ivacaftor potentiates
CFTR channel gating; CFTR protein levels barely move. If your readout is a
western blot rather than an activity assay, that difference decides whether the
experiment appears to work. Those rows are flagged "activity only".

The curated table is hand-maintained and **will go out of date** — approvals
happen constantly. Treat it as a prompt to check FDA/EMA labelling, not as an
authority.

---

## Direction: `--direction up` (default) or `down`

The logic inverts cleanly, because what you want from a promoter-bound factor
depends on both its role and your goal:

| | `--direction up` | `--direction down` |
|---|---|---|
| bound **activator** | agonise it | **block it** |
| bound **repressor** | block it | agonise it |
| example | dexamethasone → GR → GR-driven genes up | mifepristone → GR → GR-driven genes down |

Scoring shifts too: in `down` mode a proximal activator outranks everything,
because blocking is the easy half of pharmacology.

**Downregulation is much better served by the drug catalogue.** Antagonists,
inverse agonists and degraders hugely outnumber agonists, and they exist for
receptors with no useful agonist at all — tamoxifen and fulvestrant for ERα,
enzalutamide for AR, mifepristone for GR. Expect the `down` run to return more
rows than the `up` run on the same gene. There is also more non-NR reach:
plicamycin displaces SP1 from GC boxes, and bromodomain/CDK-family inhibitors
shut down transcription broadly.

**But check the direct modalities first if you want a gene knocked down.**
Approved siRNA (patisiran, givosiran, inclisiran, vutrisiran) and antisense
oligonucleotides (inotersen, nusinersen, tofersen) silence a *named transcript*
with sequence-level specificity. A TF-level intervention hits every other gene
that factor controls — for GR or SP1 that is thousands. Use this pipeline for
knockdown when a direct modality isn't available, or to explain *why* a drug
you already have lowers your gene. The same caution applies in reverse: there
is no upregulating equivalent of siRNA, which is exactly why the `up` direction
is harder and why this analysis is more useful there.

---

## Drug grading

Every TF/drug pair gets a 0-100 combined grade (A-F). On by default;
`--no-grade` to skip, `--tissue cns` to filter on exposure, `--tf-weight` to
shift the balance.

```
combined = 100 x TFconfidence^w  x  drugQuality^(1-w)        w = --tf-weight
```

**Geometric, not additive — that is the point.** The route only works if the
TF really binds your promoter **and** the drug really moves that TF. Add the
two and a superb drug rescues a worthless TF prediction; multiply them and a
weak TF hit caps the grade no matter how good the compound is, which is the
honest behaviour. Two tests assert exactly this in both directions.

Drug quality is a weighted geometric mean of:

| component | weight | source |
|---|---|---|
| potency | 0.25 | median pChEMBL for this drug **against this target** (9 → 1.0, 5 → 0.0) |
| protein selectivity | 0.15 | distinct ChEMBL mechanism targets; 10 targets scores < 0.35 |
| **gene selectivity** | 0.15 | **how many genes the TF controls** — see below |
| action confidence | 0.15 | full agonist/antagonist 1.0, partial 0.7, "modulator" 0.6, unlabelled 0.3 |
| maturity | 0.10 | max_phase; 4 → 1.0 |
| tissue access | 0.10 | only with `--tissue cns`; MW/TPSA/logP heuristic |
| **evidence depth** | 0.10 | **how much was actually measured** — see below |

### Two kinds of specificity

Molecular selectivity asks how many *proteins* a drug binds. That is the wrong
denominator here. A compound binding exactly one target, 1 nM, zero
off-targets — a perfect selectivity score — will still move tens of thousands
of promoters if that target is SP1.

So specificity is scored on two independent layers and a drug needs both:

| | question | component |
|---|---|---|
| molecular | how many proteins does the drug bind? | `selectivity` |
| transcriptional | how many genes does that protein control? | `regulon` |

SP1 scores 0.16; VDR scores 0.67. Getting your gene via a factor with 40
targets is a cleaner experiment *and* a more plausible therapy than getting it
via a global regulator, even when the motif score is worse.

Supply real numbers with `--regulon-file` (TRRUST v2 from
https://www.grnpedia.org/trrust/, or a DoRothEA export). Without one, coarse
bundled tiers are used — deliberately buckets rather than invented counts.

### Evidence depth, and the scrutiny trap

`selectivity` counts *annotated* targets. A drug studied for forty years
accumulates a long target list and looks promiscuous; an obscure compound with
three papers looks pristine. A low target count can mean selective, or can mean
nobody looked — and without a correction the score rewards ignorance. The
`evidence_depth` component separates the two, so "clean" and "unstudied" stop
looking identical.

then multiplied by a safety factor: withdrawn x0.20, black-box x0.65. Penalised
drugs are still **shown**, flagged, never silently dropped — you should see
plicamycin and decide about it yourself.

Every component is printed per row, so you can see which one is carrying or
sinking a given drug. Raise `--tf-weight` when you have real binding evidence
(`--chip-bed`); lower it when the promoter call is speculative.

These weights are **not fitted to outcome data.** Nothing here predicts whether
a drug will actually move your protein in your cells; it is a defensible
reading order for a list that is otherwise unsorted.

---

## Read this before you trust a result

**"Transcriptional activator receptors" is two different things.** Almost
certainly you want the second one:

| | activating TFs | ligand-activated (nuclear) receptors |
|---|---|---|
| how many | ~1,000 of the ~1,600 human TFs | 48 |
| examples | SP1, CREB1, NFY, ELK1, STAT3 | GR, PPARγ, RAR/RXR, VDR, TRβ, ERα, FXR, LXR |
| drugged how | not directly — via upstream signalling | small-molecule agonist binds the receptor's LBD |
| approved agonists | essentially none | dozens |

The pipeline annotates both (`role` and `nuclear_receptor` columns), but the
drug table will be dominated by nuclear receptors, because that is where the
druggable biology actually is. Outside the NR family the short list of
pharmacologically reachable activators is roughly: **NFE2L2** (KEAP1
inhibitors — dimethyl fumarate, omaveloxolone), **HIF1A/2A** (PHD inhibitors —
roxadustat and relatives), **AHR** (tapinarof), **HSF1** (indirect), and
transcription raised by broad chromatin drugs (HDAC inhibitors such as
valproate, vorinostat; DNMT inhibitors such as decitabine, azacitidine).

**A motif hit is not a binding site.** An 8–12 bp PWM matches by chance every
few kilobases per strand. A 1.2 kb promoter will "contain" sites for hundreds
of TFs — in real sequence, in shuffled sequence, in any sequence. Three things
in this package exist to fight that, and you should use all of them:

- p-value thresholding rather than the usual "80% of max score", so motifs of
  different lengths are comparable (`--pvalue`, default 1e-4)
- background composition taken from the promoter itself, so GC-rich motifs
  don't win automatically in a CpG island
- repeat masking on by default (Ensembl `mask=soft`); `--no-mask-repeats` to
  disable

None of that makes a prediction true. **`--chip-bed` is the step that does.**
Intersecting with ReMap/ChIP-Atlas peaks converts "could bind" into "was
observed bound", and it is worth more than every scoring tweak combined. The
scoring function applies a hard discount when you don't supply it.

**Binding ≠ activation ≠ more protein.** NR3C1 activates at a GRE and represses
at an nGRE in the same cell. A drug that doubles TF activity may not move your
gene at all if the limiting step is elsewhere — and for many disease genes the
bottleneck is splicing (SMN2), transcript stability, or translation, none of
which a promoter-level intervention touches.

**Not medical advice.** "Approved" means approved for its own indication. Using
an NR agonist to nudge an unrelated gene is off-label repurposing with the full
systemic pharmacology attached — glucocorticoids and retinoids in particular
are pleiotropic, and their nuclear receptors regulate thousands of genes
simultaneously. Any of this belongs in a lab, then a protocol, not in a
prescription.

---

## Pipeline

```
gene symbol
  │  Ensembl REST  lookup/symbol -> MANE Select transcript -> strand-aware TSS
  ▼
promoter window  (-1000/+200 by default, soft-masked)
  │  JASPAR CORE PFMs -> log-odds PSSM -> exact p-value threshold -> both strands
  ▼
binding-site predictions  (collapsed across overlapping windows)
  │  MyGene.info/GO: GO:0001228 activator, GO:0001227 repressor,
  │  GO:0004879 nuclear receptor + a bundled 48-NR list
  ▼
activators / repressors / nuclear receptors
  │  ChEMBL mechanism (action_type + max_phase 4) and DGIdb
  │  direction filter: activator -> agonist,  repressor -> antagonist
  ▼
ranked TFs + approved drugs
```

External services, all keyless: Ensembl REST, JASPAR, MyGene.info, ChEMBL,
DGIdb. Every response is cached under `~/.cache/promoter_tf_drugs`, so re-runs
are free and `--offline` replays them. Public endpoints do move; if one 404s,
that's the first thing to check.

### Scoring

A transparent weighted sum, printed in `results.json` so you can change it.
It is a triage heuristic, **not** a validated predictor — treat the ranking as
a reading order, not a probability.

| component | weight | rationale |
|---|---|---|
| best relative PWM score | 30 | match quality |
| number of sites | 10 | homotypic clusters are more often functional |
| proximity to TSS | 15 | `exp(-|offset|/400)` |
| annotated activator role | 15 | GO evidence, 60% credit for dual-role |
| ChIP-seq peak present | 20 | measured binding; heavily discounted without `--chip-bed` |
| approved drug available | 10 | actionability |

---

## Key options

| flag | why |
|---|---|
| `--upstream / --downstream` | 1000/200 is a core-promoter default. Many real regulatory elements are enhancers tens of kb away — see below. |
| `--all-tss` | genes with alternative promoters; scans each distinct TSS |
| `--pfm-file` | bulk JASPAR file; avoids ~900 API calls and pins the release |
| `--chip-bed` | ReMap 2022 / ChIP-Atlas peaks — the evidence step |
| `--pvalue` | 1e-4 default; 1e-5 for a short list, 1e-3 to cast wide |
| `--include-investigational` | drops the max_phase 4 filter |
| `--demo`, `--offline` | no network |

---

## Known limitations, in rough order of how much they'll bite you

1. **Promoters aren't where most regulation happens.** A ±1 kb window misses
   enhancers, which do most of the tissue-specific work. Extend with your own
   regions of interest (Hi-C/ABC-model links, or an ENCODE cCRE set) — the
   scanner takes any sequence.
2. **No chromatin context.** Pass a DNase/ATAC BED for your cell type; a site
   in closed chromatin is not available for binding.
3. **No TF expression filter.** A hepatocyte TF is irrelevant in neurons.
   Filter the output against GTEx or your own RNA-seq.
4. **GO annotations are incomplete**, so `role: unknown` is weak evidence of
   nothing, and activator/repressor assignment is context-free.
5. **No cooperativity or spacing constraints.** Nuclear receptor half-site
   spacing (DR1/DR3/DR4/DR5) determines which dimer binds; the PWMs partly
   capture this, the scoring does not.
6. **No conservation filter.** phastCons/phyloP over the site would help.

## The complementary approach you should also run

Rather than reasoning promoter → TF → drug, go straight from data:
query **LINCS L1000 / CLUE** for approved compounds whose expression signature
*already* upregulates your gene, in a relevant cell line. It skips every
inferential step above. The two approaches are strongest together — L1000 tells
you *which* drug works, this pipeline tells you *why*, and the "why" is what
lets you predict whether it will hold in your cell type. Signature data:
https://clue.io/ or https://maayanlab.cloud/sigcom-lincs/.

## Validating a hit

1. Clone the promoter window (`promoter.fa`) into a luciferase reporter;
   confirm the drug induces it.
2. Mutate the specific predicted element; the drug response should be lost.
3. ChIP-qPCR for the TF at that element, ± drug.
4. Endogenous readout in the right cell type: RT-qPCR then western/ELISA.
   Transcript up without protein up means your bottleneck was never here.
