# How the pipeline works, start to finish

[Documentation index](index.md)

This page follows one run from the question to the files it writes, and lists
every database it reads and why. It describes what the code does. It does not
say the answers are right: the last section says what has been checked.

## The question it asks

Given a gene and a direction (raise it or lower it), which approved drugs
plausibly move it, by what route, and how much would you have to believe to
test one? The output is a ranked list of hypotheses with the evidence behind
each. A drug in the list is not shown to work, and a run that finishes does not
validate anything.

You also say what is wrong with the gene, if you know. There are seven
dosage scenarios: haploinsufficiency (one copy lost, the other intact),
biallelic hypomorph (one copy weak, the other dead - the recessive
arrangement), hypomorph (both copies weak), toxic gain, null, silenced
allele, and modifier (no mutation). The
scenario fixes how far and which way the gene has to move, so the ranking can be
read against a requirement. The package reads no variant, so this is declared by
you and never inferred.

## Before the first run

The reference files are git-ignored, so a fresh checkout has none of them.
`scripts\setup_data.bat` lists what is missing and fetches it: motifs, ChIP peaks,
GENCODE, expression, LINCS, the knockout libraries, DepMap, PRISM, the GWAS
Catalog, gnomAD, ClinGen, the Orange Book, CollecTRI and BindingDB, plus the
approved drug library and the regulon built from the ChIP peaks. It leaves files
you already have alone. See [data setup](data-setup.md).

## The stages, in order

**1. Resolve the gene and pull its promoter.** The symbol goes to Ensembl (or
to a local GENCODE file and genome) and becomes a transcript, a start site and a
sequence. The transcript is the MANE Select one if there is one, then the
Ensembl canonical, then the longest protein-coding one, and the choice is
recorded. For a minus-strand gene the start is the transcript's end. With
`--all-tss` each distinct start is its own window. The default window is 1,000
bases before the start and 200 after.

**2. Scan the promoter for motifs.** Position-weight matrices from JASPAR are
slid along the sequence and a site is called at a p-value of 1e-4 (the command
default; the launcher uses 1e-5). Each strand has its own threshold, and a
matrix whose best possible match is more common than that p-value in the window
is not scanned, because every hit it gave would be a chance match at a higher
rate than stated. Overlapping calls are collapsed. The run reports how many
calls to expect by chance for that window by scanning 30 shuffled copies of it,
which keep the repeat mask and every two-base count. At these thresholds many
calls are noise, and the report says what share. When the requested window shows no
signal above chance, the same transcript is rescanned at 5,000 bases before the start
and 1,000 after, and the wider window is kept only if it shows some; the report says
what was tried either way, and `--no-widen` keeps the requested window. Sites far from
the start say less about the gene, so a widened scan is more of a lead than a core one.

**3. Say what each factor does.** Gene Ontology annotations decide whether a
factor is an activator, a repressor, both, or unknown, and a curated list marks
nuclear receptors, the one class with many small-molecule ligands. The role
matters later: whether to activate or block a factor depends on it.

**4. Look for measured evidence that the site is real.** ReMap ChIP-seq peaks
are intersected with the window. A motif says a factor could bind; a peak says it
was seen bound in some cell type. The window is also checked for other
transcription units (an antisense gene, a neighbour), which would make every
site call ambiguous.

**5. Rank the factors.** Each factor gets a score out of 100: motif quality (30),
number of sites (10), closeness to the start (15), annotated activator role
(15), a ChIP peak (20, discounted when no ChIP data exists), and whether an
approved drug exists (10). These weights are a judgement, not a fitted model. Factors that cannot work in
the target tissue's cell-cycle state (cycling, post-mitotic, mixed) are removed
at this point.

**6. Find approved drugs at each factor.** ChEMBL's mechanism table gives a
curated action (agonist, antagonist, and so on) and the approval flag
(max phase 4). DGIdb and Open Targets add recall. A drug is kept only when it
moves the factor the way the goal needs: an agonist for an activator when
raising the gene, an antagonist for it when lowering. A drug whose direction is
unrecorded is held apart as unusable, not mixed in. Drugs that arrive with an
RxNorm or NCIt id have their ChEMBL id found by name so their potency and
safety can be looked up.

**7. Grade each factor and drug pair.** The drug side combines measured potency
against that factor (pChEMBL from ChEMBL or BindingDB), selectivity (how many
targets it has), how many genes the factor controls (regulon size from TRRUST,
CollecTRI or the ChIP regulon, with DepMap knockout effects as a second view),
how well the action is defined, how mature the drug is, tissue reach (a CNS
score when asked), and how much evidence exists at all. These are combined as a
weighted geometric mean, so one very poor component pulls the whole grade down,
then multiplied down for a withdrawn drug or a boxed warning. The drug score and
the factor score are combined with `--tf-weight` (0.5, or 0.7 when ChIP data is
present).

**8. Look for routes that skip the factor.** Thirteen further searches run next
to the factor route, among them the upstream network, the known regulators and the
CTD literature record. Most produce routes that go into the same ranked list;
three (genetics, CTD and the dosage constraint) produce evidence that is
attached to the list instead, and the designed modalities produce a separate
section. [The next section](#the-routes-in-detail) takes each in turn.

**9. Put every route on one scale.** The catalogue gives each route a tier, which
says what kind of claim it is: 1 clinical, 2 measured, 3 curated, 4 predicted, 5
hypothesis. It also gives three axes (strength, specificity, efficiency) and a
distance from the gene. Within a tier the score is
`100 x (0.4 strength + 0.35 specificity + 0.25 efficiency) / (1 + 0.25 x distance)`.
Tier decides the order first, so a weak tier-1 route beats a strong tier-5 one.
Compounds that move a large share of the transcriptome in LINCS, and compounds
suspected of holding cells in mitosis when the gene is mitotic, have their score halved within
the tier. A factor's ligand found by a screen is scored like a graded one, and a
graded drug with measured binding of 1 uM or better and a recorded action is tier
2.

**10. Record what actually ran.** Every source is logged as it ran: found
something, found nothing, or failed. A failed lookup is never read as a zero.
The verdict at the top of the report is "interpretable" or "degraded", depending
on whether a pillar of the analysis is missing, and the motif scan's own quality
is part of it.

**11. Write the outputs.** `report.md` is the full account, `summary.md` the one
page (a top list with one row per mechanism family and a sentence on how
each would act), `results.json` everything, `tf_summary.tsv` the factor table,
`sites.bed` the sites for a genome browser.

The reports list only FDA approved drugs. A drug counts as approved when its
name matches an active ingredient in the FDA Orange Book, with salt forms
ignored, or when it is one of the curated gene therapies. Research compounds
from the LINCS library, drugs approved only outside the United States such as
pemafibrate, and biologics are left out of `report.md`, `summary.md`,
`tf_summary.tsv`, the panel funnel, the combined list and the dose check. The
Orange Book does not list biologics (they are in the Purple Book, which the
package does not read), so an approved antibody is hidden too. Drugs the Orange
Book lists only for local routes (topical, eye, inhaled, nasal, rectal and the
like) are hidden as well, unless the target tissue is skin: `--tissue skin`, or
an `--activity-tissue` that names skin, keeps them, and `--include-local-only`
keeps them for any tissue. A drug with any systemic route, such as tretinoin by
mouth, stays. This hides an eye drop for a retinal gene too; use
`--include-local-only` for that case. `results.json`
keeps every route, and each report says how many entries were left out. With no
Orange Book on disk the filter is not applied and the report says so. Run with
`--include-non-fda` to show everything.

## The routes in detail

Each search below asks a different question, reads different data, and enters
the ranked list in a different way. The list puts every route on one scale
(stage 9): a tier for the kind of claim, three axes (strength, specificity,
efficiency), and a distance from the gene. This table says where each search
lands.

The promoter motif route is in the table with the rest, and not above it.
The package was called `promoter-tf-drugs` after that route because it was
written first, which is a fact about the history and not a ranking. It is
the only route whose
factor-to-gene step is a sequence resemblance rather than a record of
anything, so it is the one most weakened by a noisy scan, and above a false
discovery proportion of 25% its tiers are capped (see
[methods](methods.md#gates-that-demote-a-route)). The routes that rest on a
measurement of the gene itself - the CTD record, the LINCS signature, the
PRISM phenocopy - are stronger evidence about the gene, though weaker about
mechanism.

| Search | Category in the list | Tier | Strength is | Distance |
|---|---|---|---|---|
| Drugs on the gene's own protein or transcript | direct: gene product, or gene or transcript | 1 if approved, else 3 | 0.85 if it raises the protein, 0.55 if only activity | 0 |
| Promoter motif scan, then drugs at the factors it finds | TF: orthosteric ligand, and the other TF categories | 2 with measured binding at the factor, else 3; capped at 3 when the scan is noisy | potency, action and factor confidence | 1 |
| Upstream network | upstream signalling | 4 | confidence of the best path | 1 plus the hops |
| CTD literature record | literature: reported effect | 4 | 0.5 for two or more papers, 0.35 for one | 2 |
| Known regulators (`--regulator-route`) | TF: known regulator | 4, or 3 for a direct regulator that CollecTRI and the knockout data both name with the same sign | 0.5 for a CollecTRI edge, 0.6 with knockout support, 0.35 from knockout data alone | 1 |
| LINCS signature on the gene | empirical signature | 2 if the gene was measured, 3 if imputed, 4 otherwise | size of the signature score | 0 |
| LINCS effect on a factor's own transcript | TF: abundance change | 3 if the factor's value was measured or well imputed, else 4 | how many signatures agree | 1 |
| Inferred factor activity | TF: inferred activity change | 3 | median z of the regulon | 1 |
| Receptor ligand screen, established ligand | TF: orthosteric ligand | 2 | potency, action and factor confidence | 1 |
| Receptor ligand screen, new hypothesis | virtual screen hypothesis | 5 | the model's probability, or its measured transfer AUC | 1 |
| Docking | predicted product or pocket binder | 4 | the pocket's validated AUC, not the drug's score | 0, 1 or 2 |
| PRISM phenocopy | measured phenocopy | 4 | partial correlation | 1 |
| Product ligand screen, new hypothesis | direct: predicted product binder | 5 | the model's probability | 0 |

Genetics, CTD, DepMap's knockout effects and the dosage constraint add no
route. They are described last.

### Drugs on the gene's own protein or transcript

This is the shortest route and the one a reviewer asks about first, because it
skips every inference in stages 2 to 7: no motif, no factor, no hope that a
factor binds in your cells.

Two things are collected. For the protein, approved drugs whose ChEMBL or DGIdb
target is the gene itself, kept when their action points the way the goal needs
(an inhibitor to lower it, an agonist or activator to raise it). For the
transcript or the gene, a hand-written table of approved antisense,
splice-switching and RNA-interference drugs and gene therapies (nusinersen and
risdiplam for SMN2, the exon-skipping drugs for DMD, voretigene for RPE65, and
so on). That table is maintained by hand and goes out of date, so it is a prompt
to check, not an authority.

In the list these are tier 1 when approved, since they are approved for this
gene and documented in people, and tier 3 when not. A withdrawn drug is dropped
to tier 3 and marked. Specificity comes from the modality (gene therapy 1.0,
siRNA 0.98, antisense 0.97, antibody 0.85, small molecule 0.55). One
distinction decides how the result reads: an enzyme activator raises the
activity of protein that is already there and not the amount, so a western blot
will show nothing even when it works, and the route says so.

A second part screens the protein's own known ligands with the same classifier
described under receptor screens, and keeps only compounds nobody has annotated
(tier 5). It needs at least 20 measured ligands and a classifier that separates
them from decoys, and refuses otherwise. BUB1B has one, so it was refused, and
the evidence ledger recorded that as a failure and not as an empty result.

### An upstream network search

Most of the roughly 1,600 human transcription factors have no drug of their own.
CREB1 has none, but protein kinase A does, and so does adenylyl cyclase, and so
do the receptors above them. This search walks upstream from each factor through
signed, directed edges from OmniPath (kinase to substrate, receptor to effector,
mostly from SIGNOR), then attaches approved drugs to the proteins it reaches. It
checks every reachable node against an index of approved drugs and their
targets, not only the best few, because the best-connected nodes are mostly
chromatin and repair proteins that no approved drug touches.

Four rules keep it honest.

- **Signs multiply.** Inhibiting an inhibitor activates. A path's net effect is
  the product of its edge signs, and a drug is kept only if moving its node the
  way the path needs moves the factor the way the goal needs.
- **Every hop costs.** Each edge carries a confidence of 0.75, so three hops
  leave 0.42 before the drug is considered. Specificity falls with every hop too,
  because a receptor several steps up reaches your gene through one branch and
  thousands of others through the rest.
- **Hubs are flagged.** The search is a random walk with restart, normalised for
  node degree, and any node with more than 150 edges is marked. Degree tracks how
  much a protein has been studied. If AKT1 is the top hit for every gene you try,
  that is the artefact.
- **Direction is checked at the factor.** Whether to activate or block a factor
  depends on its role at the promoter, and a factor of unknown role cannot be
  judged.

In the list these are category upstream signalling, tier 4, with a distance of
one plus the number of hops. The path is printed beside the route, so you can
judge whether it is plausible in your cell type. They are real paths in
literature-curated networks and are often three steps from the factor.

Two details change what the search can return. The index of approved drugs keeps
drugs ChEMBL records at a protein family (aspirin at cyclooxygenase, haloperidol at
the dopamine receptors), a group of complexes (hydroxyurea at ribonucleotide
reductase), a fusion protein (imatinib at BCR and ABL1) or a selectivity group, and
lists them at every member, marked as acting on the group. And the search covers the first eight ranked factors that have no drug of
their own and lists five drugs per node, which `--upstream-factors` and
`--upstream-drugs-per-node` change. Replaying it on 36 saved runs showed that asking
every ranked factor adds at most 41 drugs, while the limits hide more at busy nodes
such as NR3C1; see [checks](checks.md#why-the-upstream-group-is-missed).

### Known regulators

The promoter scan starts from the sequence and asks which factors might bind. This
route starts from the literature and asks which factors are recorded as regulating
the gene, then what approved drug acts on each. It reaches factors the scan did not
rank, such as those that bind through a partner or a distant enhancer. It is off
unless `--regulator-route` is given.

The factors come from CollecTRI, which signs each edge, and from the knockout files,
where a transcription factor that was lost or raised and moved the gene the same way
in at least two experiments counts. A factor that already has a drug route from the
scan is skipped. The drug must move the factor the way the sign and the goal need:
switch on an activator or block a repressor to raise the gene, and the reverse to
lower it. A drug recorded with an unclear action is listed and marked.

`--regulator-steps 2` also lists drugs at the factors that regulate those regulators
in CollecTRI, with the two signs multiplied (a factor whose routes disagree on the
sign is left out). Those routes sit at a distance of 2 and score lower. The launcher
uses 2 and the command line defaults to 1. `--regulator-measured` adds, for a
regulator that ChEMBL's mechanism table gives no drug for, approved drugs with a
measured activity at it. The records give no direction and name many unrelated
drugs, so it is off by default. Against CTD it reached 22.7% of the wanted drugs in the 36-gene check, at about 80 more drugs per gene.

In the list these are category TF: known regulator, at a distance of 1. Most are tier 4.
A direct regulator that CollecTRI and the knockout data both name with the same sign is
tier 3, provisionally: in the later 36-gene runs those drugs had the right direction 78%
of the time (32 against 9, so a small group), against 69% for the other predicted routes
and 59% for a random approved drug. The promotion was found on the genes it was tuned on
and has not been tested on held-out genes. Checked against CTD on 36 genes, 55% of the
drugs the route found had the right direction
(CollecTRI) and 67% (knockout data, 21 drugs), where an approved drug picked at random
has 48%. Resampling the genes gives intervals for the gain that include zero, so it
widens the list and has not been shown to sort it. CollecTRI and CTD share
literature, so agreement between them is not independent. On CYP24A1 the promoter scan
did not rank VDR, and this route found the vitamin D analogues at VDR itself (the
upstream search listed some of them only as acting on VDR upstream of other factors).
It cannot reach a drug that
ChEMBL does not record at the factor, which is why rifampicin at PXR is missed. See
[checks](checks.md#what-starting-from-regulators-adds).

### LINCS signatures

The L1000 experiment measured what happened to cells' transcripts when they were
dosed with a compound. That makes it the only link between a drug and expression
that was observed and not inferred, and it is used twice.

The first use asks which approved compounds moved this gene the way you want. The
local files are Enrichr's consensus signatures: for each compound, cell line,
dose and time, the set of genes that went up and the set that went down. They
carry membership and no magnitude, so a hit is scored at a nominal value and
labelled as set membership. The web service version carries z scores.

Four things decide how far to trust a hit.

- **Measured or imputed.** L1000 physically measures about 978 landmark genes and
  infers the remaining roughly 11,350 from them. A hit on a landmark gene is a
  measurement (tier 2). A hit on a gene the experiment predicted with high
  confidence is imputed (tier 3). Anything else is tier 4.
- **Cell line.** If the gene is not expressed in the line, nothing can move and
  silence means nothing. The cell line is printed.
- **Stress.** Many strong signatures are some version of "this cell is dying".
  Compounds that connect to a large share of unrelated queries are marked
  non-specific and left out. In addition, a compound whose signatures are among
  the longest 10 percent in the files has its score halved within its tier.
- **Dose and time.** The dose is usually 10 micromolar, which matters for whether
  a person could reach it (see `exposure`).

The second use asks which compounds change a factor's own transcript the needed
way. That gives the TF: abundance change routes (tier 3 or 4), where the drug
changes how much of the factor the cell makes and not what it does.

### Factor activity inferred from target expression

A motif scan says a factor could bind. This asks whether the factor is
demonstrably active, from the opposite data: if a factor is active in a tissue,
the genes it controls leave a footprint in that tissue's expression. For each
factor with at least 10 measured targets, the Human Protein Atlas tissue
signature is rank-normalised and tested against the factor's signed regulon.
CollecTRI is used when present because it signs every edge, falling back to
TRRUST, which leaves about half its edges unsigned. Two statistics are
reported: a directional one over signed edges (positive means the factor looks
activated), and a magnitude one over all edges. A call is directional
only when at least about a third of its edges carry a sign. The p-value assumes the
targets are independent, which they are not, so the launcher uses 2,000
permutations for an empirical null. Two factors whose regulons overlap by half or
more cannot be told apart, and the call says so.

It does two jobs. It corrects the motif scan: a factor with a strong site and a
flat regulon is probably scan noise, and one with a weak site and a coherently
moved regulon is probably real. And it is a drug source that needs no annotated
binding: the same statistic is run on each compound's LINCS signatures, and a
compound counts as moving a factor when it does so the same way across enough
profiles and cell lines. Those are the TF: inferred activity change routes, at
tier 3.

The limit that matters: the footprint is a reference tissue average, not your
sample. A receptor can be transcribed everywhere and active only where its
ligand reaches, so a quiet regulon does not mean a nuclear receptor is absent.

### Receptor ligand screens

For the top-scoring factors that are receptors, the package asks which approved
drugs bind them, including ones nobody has annotated there. Known ligands (ChEMBL
activity at pChEMBL 6 or better, or a local BindingDB table) train a random
forest on molecular fingerprints. The model is used only if it has at least 20
actives and separates them from decoys at an AUC of 0.70 or better. It then
scores every approved drug that has a structure and puts each in one of four
groups.

- **Established ligand.** Already measured at this receptor. Not a discovery, but
  a positive control: a model that ranks these low is not working.
- **Confident prediction.** Close to known ligands and predicted active.
- **New hypothesis.** Unlike any known ligand (nearest Tanimoto below 0.55) and
  predicted active at 0.70 or better. This is the repurposing group, and it is an
  extrapolation into chemistry the model has not seen.
- **No call.** Too far outside the model's domain (nearest ligand below 0.35).

A receptor with too few ligands of its own borrows its relatives', and the
transfer is measured against ligands the model never saw before anything is
reported.

An established ligand enters at tier 2 as TF: orthosteric ligand, with a strength
built from its measured potency, whether its action at the factor is known to
point the right way, and the factor's own confidence, the same formula as a
graded drug. A new hypothesis enters at tier 5 and states that it was never
tested at that receptor. Most new hypotheses will be wrong, which is normal for
virtual screening. Unless the receptor has enough drugs with an agonist or
antagonist label, the screen finds binders only and says that direction is not
known.

### Docking

Docking asks whether a compound fits a pocket. Its score is a poor estimate of
affinity even when the pocket is good, so the package first asks whether the
pocket can be trusted.

The structure comes from the PDB and must carry a co-crystallised ligand, which
places the search box. For a nuclear receptor the conformation (activating or
inactivating) can be read from whether a coactivator or corepressor peptide is
bound; for a kinase, a GPCR or an enzyme it cannot, and the run says the
conformation is unverified. Then the gate: the pocket is asked to recover its
own known ligands. At least 5 actives (up to 60) and 20 decoys must dock, the
AUC must be 0.70 or better and the enrichment in the top one percent at least 2.
A pocket that fails is refused, and the refusal is a result with its reason. Only
after it passes is the whole approved library docked, at exhaustiveness 16, with
gnina (which also gives a neural-network pose score) or smina.

Two kinds of pocket are docked. The gene's own protein is docked with the
conformation left open, since a binder there changes what the protein does. The
five best-ranked factors are docked in the conformation the goal needs: lowering
a gene through an activator means docking the activator's inactive pocket.

In the list a docking hit is tier 4. Its strength is the pocket's validated AUC
and not the drug's score, because the score separates binders from non-binders
and ranks binders badly. A hit on the gene's own protein is a direct product
binder (distance 0), a hit on a factor found at the promoter is a pocket binder
(distance 1), and a hit on a protein that is not at the promoter is marked
off-pathway (distance 2). Docking says nothing about whether the compound
activates or blocks. Where a receptor was both docked and screened by ligand
model, the two are intersected, since they share no input and agreement is worth
more than either alone. The combined list fuses docking with the other routes by
rank, so a validated pocket is not buried under tier-3 routes.

### Genetics

Three separate lookups, none of which makes a route.

- **The gene.** Open Targets and, if present, a local GWAS Catalog table say
  which traits human variation in this gene is linked to. Variants that change a
  gene's expression are lifelong natural experiments for the intervention you
  propose, unwanted effects included, so this is evidence about what moving the
  gene does to a person.
- **The factor.** Whether variation near each of the top three factors is linked
  to the gene.
- **Trans-eQTLs.** eQTLGen's significant trans-eQTL table (59,786 of them over
  6,297 genes) is checked for variants near each of the top six factors that move
  this gene's expression in blood. This is the strongest link the package can
  show: natural genetic variation at the factor changes the gene's expression in
  real people, with exposure randomised by meiosis. It is also usually empty. In
  the ABCA4, BUB1B and CCR5 runs none of the six factors had support, and the
  ledger says "nothing found" and not "failed".

GWAS and PheWAS do not record whether a drug changes a gene's expression, so
they cannot be used for that, and the report says so.

### CTD

The Comparative Toxicogenomics Database curates, from papers, which chemicals
raise or lower a gene's expression. It is asked once per gene. The answer does
three things. Each candidate is marked previously linked when CTD already records
it moving this gene, and novel otherwise. Chemicals CTD records that the
pipeline did not predict are listed. And the discovery step compares the novel
candidates with the best already-linked one on the three axes to say which beat
what exists. The lookup can be answered "nothing known" or can fail, and those
are kept apart: a failed lookup means nothing can be called novel.

CTD's own action list names every action in a sentence, not only the ones about the
gene, so the package reads the sentence. Across 36 genes only 39% of the rows with an
expression action describe the chemical changing this gene's mRNA. The rest are protein
expression (26%), a chemical that blocks or promotes another chemical's effect (24%),
another gene's expression (8%) and effects through another gene (3%). Only the first
counts as a direction. A chemical that modulates another's effect is kept as a known
link with no direction. See [checks](checks.md#what-ctds-expression-record-contains).

CTD also supplies a route of its own, **literature: reported effect**: one row for
each chemical CTD records as moving this gene the wanted way (exposures and MeSH headings
left out), at tier 4 and a distance of 2 because the path from the drug to the gene is
unknown. Its strength is 0.5 for two or more papers and 0.35 for one, which is a fixed
rule and not a calibration. It says what somebody saw and cannot name a drug nobody
has tested on this gene. It is left out of the retrospective check, since CTD is the
answer key there, and `--no-ctd-route` turns it off. CTD's reports on the gene's
protein level ("increased expression of the protein") are read too, as a separate
stream. A chemical reported at both levels is one row, a chemical reported only for
the protein is a row that says so, and it is not read as a sign that the promoter is
involved, since protein level can change through stability or translation. See
[checks](checks.md#using-ctd-directly). The run then compares CTD with LINCS and with the
routes: a ranked drug that CTD records moving the gene the other way is marked, and drugs on
which the two sources give opposite directions are listed. A route with no recorded
direction at its factor for such a drug is lowered to tier 4, and any other is only flagged.
See [checks](checks.md#when-ctd-records-the-opposite-direction).

CTD does not change another route's score. It does set the label of receptor-screen
hits: an established ligand that CTD already links to the gene is labelled as
confirming a reported link and is not offered as a new route. Its records cover
toxicology doses and exposures as readily as clinical findings, so a link here
is weaker evidence than it sounds.

### DepMap and PRISM

Both come from cancer cell lines growing in a dish, and both say nothing about a
differentiated or non-dividing tissue.

DepMap's CRISPR screens say what a cell loses when a factor is removed. For each
candidate factor the package reports the mean effect, how many lines depend on
it (below a score of -0.5), and whether that dependence is narrow or everywhere.
A factor essential everywhere is a poor target for a chronic drug. This is a
cost and not a benefit, and it feeds the burden stage of the funnel.

PRISM is the other half of the same idea. It measures how roughly 6,000 drugs
affect the viability of the same lines. If a drug and a factor's knockout hurt
the same lines, the drug behaves as though it removes the factor, with no record
of it binding anything. The number reported is a partial correlation that
controls for how easily each line dies under the whole panel, because cell-line
fragility otherwise creates a large correlation between almost any drug and
almost any gene. The raw value is printed beside it, and a large gap means it
was mostly about the lines. Comparisons across too few lines are refused. These
become the measured phenocopy routes at tier 4. A drug that phenocopies a
knockout serves a goal that needs the factor lowered and works against one that
needs it raised, and routes pointing the wrong way are dropped.

### Knockout signatures

The step from a factor to your gene was the one the package could not see. A
drug may bind a factor, and the factor may sit on your promoter, and your gene
may still not respond to it. The knockout libraries from Enrichr say what
experiments found when a factor was lost or raised: the genes that changed
most. Two are read, from local files. `LINCS_L1000_CRISPR_KO_Consensus_Sigs`
has the consensus up and down genes for 5,212 CRISPR knockouts in the L1000
cell lines. `TF_Perturbations_Followed_by_Expression` has knockouts, knockdowns
and overexpression of transcription factors from GEO, in human and mouse.

For each factor on the list, the report says whether your gene is among the
changed genes, and which way. A gene that falls when the factor is lost is one
the factor activates, and one that rises is one it represses. Overexpression
reads the other way round. Mutation, truncation and fusion experiments are not
read, since they can raise or lower activity. The answer is one of three:
a measured response (with how many experiments and which way), perturbed with
the gene not among the top changes, or factor not tested. The second is weak,
because each list holds only a few hundred genes, so a gene can respond and
not be on it. The cells are cancer lines, embryonic cells and mouse tissue. A
response says the factor can move the gene somewhere. It does not say that a
drug acting on the factor will move it in your cells. Where the measured
direction contradicts the factor's annotated role, the report says so.

This replaces a fixed placeholder for responsiveness (0.40 for every factor).
The placeholder is still saved in `results.json` and is no longer shown.

The same files are read the other way for the known regulators route: for one gene,
which factors changed it when they were lost or raised.

### Designed modalities

For a gene that has no approved gene-directed drug, `--design-modality` writes
a design instead of looking one up, and puts the answer behind deliverability.
For lowering, it proposes siRNA (21-mers chosen by GC content and strand
asymmetry), RNase H gapmers (16-mers), CRISPRi guides in the window 50 bases
before to 300 after the start site, and nuclease knockout or base editing. For
raising it proposes a splice-switching oligonucleotide, CRISPRa and gene
addition by AAV, with notes on what each requires and no sequences, since there
is no general antisense for turning a gene up.

Delivery comes first. Liver is solved (GalNAc conjugates, under the skin), the
central nervous system is reachable by injection into the spinal canal, the eye
by injection into the vitreous and the lung by inhalation. T cells, macrophages
and microglia have no routine systemic route, so a potent oligo for them is a
cell-culture reagent until one exists. The design says whether delivery or
design is the limit. Local checks only are done for off-target effects, and every
candidate has to be searched against the transcriptome before it means anything.
These are not ranked routes. They appear as their own section, with the screen
that would test them.

### The dosage constraint

Every route ends with "and your gene moves". Two things nothing else checks are
whether there is room to move it and how far it must move.

gnomAD counts loss-of-function variants across about 800,000 people against the
number expected, and ClinGen's panel curates whether losing or gaining a copy is
harmful. Together they say whether the gene's dose has been under selection. A
gene can be tolerant when one copy is lost and still cause severe disease when
both are lost (mosaic variegated aneuploidy is two broken BUB1B copies, and
carriers are healthy), so a tolerant score is never read as "this gene does not
matter".

The scenario you declare then sets the requirement from the mechanism, not from a
prediction. One working copy of two means roughly doubling the working copy to
restore the total. A toxic product to be halved is a stated target. A null allele
cannot be helped by raising transcription at all, and the report says so and
treats every route as context. The requirement is stated in words and is not
compared with any route in the same units. The package reports which routes carry
a measured effect (a LINCS hit with a z score) and which carry none, and never
converts a z score to a fold change, since that needs the plate's spread, which
the data does not contain. A route with a good score is therefore never claimed
to deliver "enough".

## What the launcher adds

`scripts/run_everything.bat` wraps this in five stages. It docks the gene's own
protein, docks the five best-ranked factors, runs the analysis with the docking
results folded in (the upstream network, the known regulators, LINCS, genetics, the
virtual screens, the DNA-directed agents and the sensitivity scan all on), narrows
the routes through a funnel against a control gene (ACTB) in `panel`, and builds one
list from the panel's survivors and every validated docking result in `combined`.
Before any of that it checks ten key reference files and offers to run
`setup_data.bat` when one is missing.

Two optional steps follow the five stages. The dose check runs by default and reads
the finished results: it builds the FDA label table of peak blood levels for the
candidates and sets them against the concentration the evidence was seen at. The
exploratory dock does not run by default. The launcher asks, with an answer of no
unless you type y, whether to run it when the gene's own pocket is refused. If you
decline, it prints the command for running it later. When it runs, it finds its own
template structure by searching the PDB, so it needs no input, and it docks up to
three distinct ligand sites that the search finds, each as a separate run. The result
is not validated and takes hours.

Before stage 1 the launcher prints a line for each route that its
data files rule out: the PRISM and DepMap matrices for the measured phenocopy
route, the GWAS Catalog table, and a tissue for factor activity. It offers to
download the GWAS Catalog (about 70 MB) and lists the tissues in the expression
table so you can pick one by number. It also offers the DepMap 24Q2 and PRISM 24Q2
tables (about 420 MB and 70 MB), fetched from figshare because the DepMap
portal asks for a browser check. Each file is kept only if its published MD5
checksum matches. Both are needed for the measured phenocopy route. The funnel removes routes in stages
(evidence floor, genome-wide agents, large regulons, wide tissue footprint,
unachievable selectivity, heavy clinical burden) and reports what each stage
removed.

## Checks run on saved results

- `retrospective` asks whether the routes point at drugs already known to move a
  gene, against CTD and LINCS, with a drug-fame control.
- `scripts/positive_fate.py` follows each drug CTD records as raising a gene through
  the pipeline and says where it was lost. `scripts/upstream_misses.py` replays the
  upstream search to say which filter removed the drugs within reach.
  `scripts/regulator_gain.py` scores the known regulators route against CTD, with
  intervals. `scripts/ctd_check.py` sorts CTD's expression rows by what the sentence
  says. `scripts/drug_trace.py GENE up DRUG` follows one drug through the chain and
  stops at the first missing link.
- `exposure` sets a published human blood level against the concentration the
  evidence was seen at.
- `exploratory` docks into a pocket transplanted onto an AlphaFold model, from a
  template it finds by sequence search, labelled unvalidated, with a class check that
  has to pass. It is run on request and never by default.

## Databases it reads, and why

"Live" means a request during the run, cached for 30 days. "Local" means a file
under `data/` that you download once.

| Source | Kind | Used for |
|---|---|---|
| Ensembl REST | live | transcripts, start sites, promoter sequence, neighbouring genes, isoforms for brake exons |
| GENCODE (GTF, hg38 genome) | local | the same without the network, and the neighbour check |
| JASPAR 2024 | local or live | the binding motifs that are scanned |
| ReMap 2022 | local or live | measured ChIP occupancy, and a ChIP-derived regulon |
| Gene Ontology (GOA), mygene.info, HGNC | local or live | factor role, nuclear receptor flag, symbol and UniProt resolution |
| ChEMBL | live or local | drug mechanisms and approval, potency against a target, drug warnings, ligand sets for screens and the docking gate, the approved library |
| DGIdb | live | broader drug-gene interactions, for recall |
| Open Targets | live | drug and target links, tissue expression, and the traits human variation in the gene is linked to |
| BindingDB | local | measured affinities and known ligands for the docking gate |
| Broad Drug Repurposing Hub | local | the screening library, clinical phase and targets |
| FDA Orange Book | local | route of administration |
| TRRUST, CollecTRI (via OmniPath) | local | how many genes a factor controls, and signed edges for the activity test |
| OmniPath | live | signed upstream network edges |
| LINCS L1000 (Enrichr or SigCom signatures, CLUE gene info) | local or live | measured drug effects on the gene, and how broadly a compound acts |
| Human Protein Atlas (local table), Open Targets expression (live) | local or live | where the gene and the factor are expressed, and the tissue signature; a GTEx table is also accepted |
| DepMap and PRISM | local | cell dependency on a factor and drug sensitivity patterns |
| gnomAD, ClinGen | local | how much dosage change the gene tolerates, and curated dosage sensitivity |
| eQTLGen, GWAS Catalog | local or live | human genetic evidence that a factor moves the gene |
| Enrichr knockout libraries (LINCS CRISPR KO, TF perturbations) | local | whether your gene changes when a factor is lost |
| CTD | live | curated chemical effects on the gene's expression, to mark what is already known |
| RCSB PDB | live | experimental structures with ligands for docking |
| AlphaFold DB | live | predicted structures, used only when the exploratory dock is asked for |

Not every source is read in every run. The command `--check` lists which local
files are present, and the evidence ledger in each report says which sources
ran. Several of these sources set their own terms (CTD, ReMap, TRRUST and CLUE
restrict commercial use), so read them before sharing `data/` or a results
archive.

## What has and has not been checked

The package has tests, and tests show the code does what it says, not that the
science is right. Against saved results for 36 genes, the retrospective check
found that the routes recover 9 percent of the approved drugs CTD records as
raising a gene, rank them no better than chance once a drug's fame is held equal,
and place the known drug among the candidates for 2 of 6 textbook cases. One
cause is the scan window, which misses factors that bind far from the start.
Read a candidate list as a reading order. [Checks on candidate
drugs](checks.md) has the numbers and their limits.
