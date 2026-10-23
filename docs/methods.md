# Methods and limitations

[Documentation index](index.md)

This page describes the implementation and the assumptions encoded in it. It is
not an external validation of its scientific models or curated claims.

## From gene to candidate routes

1. Ensembl resolves a gene and selects a promoter. Transcript selection prefers
   MANE Select, then canonical, then a fallback transcript. `--all-tss` considers
   distinct clusters of starts rather than every transcript separately.
2. JASPAR position frequency matrices become log-odds motifs. The scanner checks
   both strands and normally skips repeats and ambiguous bases. Overlapping hits
   are collapsed. Its threshold is per position, not a study-wide false discovery
   rate. The scan-quality summary compares site counts with the sites found in 30
   shuffled copies of the window, and gives a range. A matrix that cannot reach the
   requested p-value in the window is left out and counted.
3. Annotation identifies TF roles and nuclear receptors. Optional ChIP overlap
   provides another source of evidence. It counts peaks for each factor over the
   whole promoter window, so it says the factor was seen somewhere in the window,
   not at the motif site. A motif match alone is not observed binding.
4. ChEMBL and optional DGIdb queries associate drugs with factors. The requested
   direction and the annotated TF role determine which actions are retained.
5. Scoring modules rank TF/drug pairs. Additional modules propose direct,
   upstream, ligand-axis, stabilization, and empirical routes. The catalogue
   combines them into a table ordered by evidence tier and score.

## Scores

The TF score allocates 30 points to motif quality, 10 to site count, 15 to
proximity, 15 to annotated role, 20 to ChIP evidence, and 10 to druggability.
Without ChIP input, the code grants discounted motif-derived credit in that
component. These are chosen weights, not fitted clinical estimates.

Drug grading combines evidence and selectivity with TF confidence. `--tf-weight`
controls the latter's share. Regulon size estimates how broadly a TF may act;
missing local evidence can fall back to curated tiers. Protein selectivity and
gene-expression specificity are different quantities.

The outcome module multiplies heuristic terms for engagement, direction,
promoter regulation, and responsiveness. Numbers labeled as probabilities in
reports are not calibrated success probabilities. Use the limiting-step notes
to identify missing evidence, not to claim a numerical chance of efficacy.

## Optional analyses

| Module or option | What it contributes | Main limit |
| --- | --- | --- |
| Direct routes | Drugs and modalities associated with the gene itself. | Acting on its protein does not necessarily change its transcription. |
| `--upstream-network` | Signed paths to druggable upstream targets. | Sign and confidence depend on curated edges; longer paths incur penalties. Every reachable node is checked against an index of approved-drug targets built from the local drug table (about 80 ChEMBL requests the first time). Without the table, or if the index fails, each node is asked about separately and the search stops after 100 per factor; the evidence ledger says which. |
| `--regulator-route` | Approved drugs at factors CollecTRI and the knockout data record as regulating the gene, including factors the motif scan did not rank. | The links come from a database and from knockouts in other cells. CollecTRI and CTD share literature, so agreement between them is not independent. Checked against CTD on 36 genes, 55% of these drugs had the right direction, where 48% would by chance; that gain is within the noise (95% interval -1% to +14%). |
| Ligand axes | Curated synthesis, catabolism, sequestration, and supplementation routes. | Context and feedback are not quantitatively modeled. |
| Mechanism comparison | Orthosteric, allosteric, and related mechanism classes. | Class labels do not establish a useful effect in the target tissue. |
| Stabilization | Curated stabilizers, glues, and degradation-related routes. | Coverage is incomplete and target effects need confirmation. |
| `--cell-type-check` | Expression overlap between the gene and TFs. | Tissue-level expression is not cell-specific regulation or exposure. |
| `--tissue cns` | Partial CNS MPO properties and curated delivery flags. | Missing properties, transport, dose, and active efflux limit interpretation. |
| `--proliferation` | Curated rules that flag, penalize, or remove routes. | They are assumptions about tissue state, not observations from the sample. |
| `--lincs` / `--lincs-file` | Signed expression-effect records and stress flags. | Cell line, dose, time, and measured versus inferred status matter. |
| `--genetics` | Trait associations and guidance for regulatory support. | Associations alone do not establish a causal drug effect. |
| `--predict-direction` | Structure-based direction models when training is possible. | Training coverage and chemical similarity constrain applicability. |
| `--screen` | Candidate ranking from approved compound structures. | Model scores and similarity are hypotheses, not binding measurements. |
| `--design-screen` | Reporter, guide, control, and compound suggestions. | Plans require experimental review; guide off-target analysis is absent. |
| `--design-modality` | Sequence-based direct modality suggestions. | Delivery and off-target validation remain separate work. |
| `docking.py` helpers | Optional structure and docking utilities. | Docking scores do not establish binding or agonism; external tools are needed. |

The signature module includes connectivity helpers, but the main pipeline does
not currently call them for `--signature-file`. The genetics module also includes
manual lookup guidance; do not read that as a completed trans-eQTL analysis.

## Gates that demote a route

Six checks can lower a route after it is built, and a seventh warns
without touching one. None deletes anything, and each says on the row what
it did and why. They were added after a review of a
BUB1B run that produced confident tier-2 candidates on evidence that did not
support them.

**A noisy scan caps a motif-derived tier.** `motifs.scan_quality` reports the
false discovery proportion: the share of the site list expected on shuffled
copies of the same window. Above 25%, routes whose factor came from that scan
are capped at tier 3. Tier 2 means the route's evidence is observed, and for
these routes only the drug-to-factor half is observed; the factor-to-gene half
is a sequence resemblance. A scan that is half noise cannot support it. Routes
that do not rest on the scan — CTD, LINCS, the gene's own product, the
sensitivity correlation — are untouched.

**A compound that moves the cell cycle is halved.** For a gene whose
transcript tracks the cell cycle (`cellcycle.MITOTIC_GENES`),
`proliferation.assess` asks whether each candidate also moves the rest of the
mitotic programme the same way across its LINCS signatures. The bar is the
95th percentile of all compounds, not a flat share, because the background is
not zero: about 1.9% of signatures raise three or more mitotic markers on
their own, and about 10.9% lower them. A compound above the bar is most simply
read as changing how many cells are dividing. This is a measurement and it
takes precedence over the hand-written `systemic.MITOTIC_DISRUPTORS` list. A
compound the files do not cover gets no verdict, which is not a pass.

**An agent that causes the defect being treated sinks.** With
`--avoid-genotoxic`, a compound that works by damaging DNA or the spindle —
alkylators, platinums, topoisomerase poisons, azacitidine, decitabine, the
taxanes — drops to tier 5 with its score cut tenfold. It is for a patient who
already has chromosomal instability or is cancer-predisposed, where such an
agent causes the problem it is being offered for. The rows stay in the report
so the search remains visible.

**A compound the literature records for every gene loses a little
score.** This is the control the retrospective found mattered most. Across
the saved runs, a baseline that counts only how many OTHER genes record a
drug predicts CTD with an AUC of 0.93, while the routes score 0.53, and
among drugs of equal fame the routes fall to chance (0.53 on
reported-or-not, 0.50 on the right direction). A reader of one gene's
report never saw that finding.
`scripts/build_fame_table.py` counts, over your saved results, how many
distinct genes name each chemical; a candidate in the top decile loses 15%
of its score and keeps its tier. Slight on purpose: being well studied is
not evidence of being wrong, and the claim is only that such a row says
little about your gene in particular. Tetrachlorodibenzodioxin,
benzo(a)pyrene and bisphenol A are recorded for every gene in the runs
here, which is what the confound looks like. The table's resolution is
however many genes you have run, every message states the denominator, and
a compound the table does not hold gets no verdict.

**An agent that raises an oncogene sinks too.** The same
`--avoid-genotoxic` flag reaches a second way: not a drug that damages
DNA, but a TARGET whose dose the route moves in the direction that causes
cancer - raising a proto-oncogene, or lowering a tumour suppressor. The
abundance scan offered to raise MYC, FOXM1 and ETS2 for a
cancer-predisposed patient and nothing said so. Raising a tumour
suppressor is left alone, because more of it is not the hazard, and a
route whose direction is unrecorded is left alone rather than guessed at.

**A declared scenario is checked against the population data.** This one
warns rather than demoting, because the declaration is the clinician's.
`--scenario haploinsufficiency` asserts one allele is lost and the other is
intact. When the gene's pLI says losing one copy is not selected against, that
assertion disagrees with the frequency data: healthy people carry a broken
copy, so the condition is probably recessive and the second allele is probably
not wild type either. The report says so and points at
`--scenario biallelic_hypomorph`, which does not assume a wild-type allele.
The check is reported rather than refused, because the constraint tables are
read during the analysis and a gene can be haploinsufficient in one tissue and
tolerated in the population.

## The six steps between a gene and a protein

How much working protein a gene yields is set by six steps in series:
transcription, splicing and isoform choice, mRNA stability, translation,
folding, and degradation. This package was built around the first one, so
a report can fill two hundred rows about transcription and say nothing
about the other five without anything looking wrong.

Every report therefore opens with an audit of all six, which marks each
one as **searched**, **listed only** or **not covered**:

- *searched* — a gene-specific lookup ran and returned a count. The
  count may be zero, and zero is then a result.
- *listed only* — the mechanism and the drug classes that work that way
  are named, and nothing was looked up for this gene. A menu, not a
  search. Most of the post-transcriptional section is this.
- *not covered* — neither.

On a real BUB1B run the audit reports several hundred transcriptional
routes against a couple of dozen from every other step combined. That
ratio is a property of the package rather than of the biology, and the
audit says so, because a full table is not evidence that the other steps
were considered. Each gap names what would close it, so the limitation is
a concrete next step rather than a disclaimer.

This is the evidence ledger's idea applied to biology instead of to data
sources: the thing most likely to mislead is not a wrong number, it is a
question nobody asked.

### What each of the other five steps actually looks at

| Step | What runs | Source |
|---|---|---|
| Splicing and isoform choice | exons that only decay-prone isoforms use, each weighed by the share of the gene's message its isoforms carry | Ensembl isoforms; GTEx v8 median transcript TPM |
| mRNA stability | the 3'UTR scanned for AU-rich elements, Pumilio response elements and DRACH, each against a shuffle of the same sequence | the transcript's own sequence |
| Translation | the 5'UTR scanned for upstream ORFs with Kozak strength, G-quadruplexes, and a 5'TOP motif at the cap | the transcript's own sequence |

`--promoter-dna` scans the promoter for the features a DNA-binding agent
prefers. It reads the same sequence as the motif scan and honours the same
repeat mask: soft-masked (lowercase) bases are skipped, and each unmasked
stretch is scanned on its own so no feature is called across a join. That
matters because Alu and simple repeats are exactly what a GC-rich or
homopolymer pattern matches best, and the shuffled null cannot catch it —
the null is built from the same sequence, so a repeat is in both.
| Folding | a ligand screen against the gene's own protein, which is what a pharmacological chaperone would have to be | ChEMBL, and docking where a pocket passes its control |
| Degradation | the ligase that destroys the protein, the enzyme that reverses it, and compounds measured on the protein itself | curated degron and DUB tables; a drug-to-protein table |

Every sequence feature above is tested against a dinucleotide shuffle of
the same sequence, so a hit means more of the element than this base
composition gives by chance. None of them is a measurement in your cells,
and none has a compound matched to it.

### Weighing a brake exon

Listing exons that only decay-prone isoforms use is the easy half. The
list is drawn from the annotation, and the annotation says an isoform has
been *seen*, not that it carries any appreciable share of the message —
so twenty candidate exons read as twenty opportunities whether or not any
of them matters.

`isoforms.py` supplies the missing number from GTEx v8 median transcript
TPM, by default in cultured fibroblasts. Each candidate exon then carries
the share of the gene's message its isoforms hold, and the step as a
whole gets a verdict: decay isoforms carrying a real share make splice
correction a lever, and decay isoforms carrying none mean there is
nothing to redirect.

Two cautions are built into the output rather than left to the reader.
GTEx v8 is built on GENCODE v26, so isoforms annotated since then have
**no measurement at all** — they are counted separately and never folded
in as zeros. And a transcript's share is not its flux: a decay isoform is
destroyed quickly, so a low steady-state share is partly a consequence of
the decay it is named for. The share bounds what redirecting it could add
to the pool, which is the question, and is not the rate through it.

recount3 would answer the same question over far more samples. It is
distributed as R objects, so it belongs in a table you supply at
`data/transcript_expression.tsv` rather than in a live lookup; a file
there is read in preference to the API.

### Naming the factor, where a motif cannot

The UTR scans say an element is in excess. They cannot say which
microRNA or which RNA-binding protein uses it, and that name is what
turns an annotation into something to act on. The datasets that carry the
names are bulk downloads under their own licences, so none ships with the
package; all of them reshape into one table at
`data/transcript_sites.tsv` with columns `gene, factor, kind [, region,
evidence, source]`.

| Source | What it adds | Where |
|---|---|---|
| TargetScan | predicted microRNA sites, conserved families | targetscan.org, `Predicted_Targets_Info.txt` |
| miRTarBase | microRNA targets with experimental support | mirtarbase.cuhk.edu.cn, `hsa_MTI.csv` |
| ENCODE eCLIP | measured footprints for ~150 RNA-binding proteins | encodeproject.org, eCLIP IDR peaks (BED) |
| m6A atlases (REPIC, m6A-Atlas) | measured methylation sites rather than DRACH motifs | the atlas's own download |

`scripts/build_transcript_sites.py` does the reshaping, and
`scripts/setup_data.py --only transcript_sites` fetches the three
downloads and runs it. A miRTarBase row carries that entry's own
experiment list, so a reporter assay is not reported as the same thing as
a CLIP hit. Two miRTarBase classes are left out. `Non-Functional MTI`
records that the pair does *not* act, so listing it as a regulator would
say the opposite of the measurement. `Functional MTI (Weak)` is a hit
from a high-throughput screen rather than a direct assay, and in release
11 it is 7,309,069 of the 7,339,072 human rows — keeping it makes the
table 880 MB, which every run then reads to answer about one gene.
`--mirtarbase-weak` puts that class back with every row labelled weak.

An eCLIP peak is an interval, mapped to a gene by overlap on
the same strand with the gene body, introns included, so those rows say
`gene body` rather than claiming the protein sits in a UTR. That is a
broad criterion. A well expressed gene carries peaks for most of the
assayed proteins — BUB1B matches 148 of the 168 — so the list says which
proteins have been seen on the transcript, and does not rank them or
argue that any one of them sets its level.

Rows marked `predicted` are reported separately from rows with
experimental evidence, because a TargetScan prediction and an eCLIP peak
are not the same claim. Ribosome profiling (GWIPS-viz, Trips-Viz) would
replace the 5'UTR sequence proxies with observed translation efficiency;
nothing reads it yet, and the coverage audit names it as the gap.

## Measured on the protein, not the transcript

Every other abundance claim here is about mRNA. `abundance.py` asks
whether a compound moved a factor's own transcript; `signatures.py` asks
the same of the target gene. Where the pathogenic allele makes an
impaired or unstable protein, that gap is the whole question, because a
stabiliser or a degradation blocker raises protein with the transcript
flat and no other route in this package can see one.

`proteomics.py` reads a drug-to-protein-abundance table and builds tier-2
routes from it — the strongest non-clinical evidence here, because the
measurement is on the thing the patient needs with nothing inferred in
between. It says nothing about mechanism: a stabiliser, a translation
effect and a transcriptional one are indistinguishable in the table and
separate on an RT-qPCR.

Which datasets answer this question matters more than the reader, because
two that are routinely named do not:

| Resource | What it measures | Answers this? |
|---|---|---|
| decryptE | dose-resolved proteome under ~144 drugs, Jurkat | yes |
| ProTargetMiner | proteome under ~50 anticancer compounds, A549 | yes |
| PRIDE / ProteomeXchange, MassIVE | thousands of individual drug-treatment studies | per study |
| LINCS P100 | ~90 phosphopeptides, signalling state | no |
| LINCS GCP | ~60 histone marks | no |
| CPTAC | tumour proteogenomics, no drug arm | no |
| DepMap proteomics | baseline protein levels, no drug arm | no |
| Thermal proteome profiling (Meltome) | melting shift, target engagement | no |
| Kinobeads chemoproteomics | binding of ~240 kinase inhibitors | no |

Coverage is the honest limit. decryptE and ProTargetMiner together reach
a low hundreds of compounds against an approved library of several
thousand, and each is a single cell line, so a compound absent from them
is unexamined and not negative. The repositories in the third row cover
far more compounds between them, but each study has its own design, dose
and normalisation, so they are the way to answer a question about one
compound you already care about rather than a corpus to screen against.

The last two rows answer a neighbouring question that is worth not
confusing with this one. Binding is not abundance: a drug can engage a
protein without changing how much of it exists. They are still worth
consulting when the gene is itself a kinase, because their coverage is
far better than the abundance datasets and a measured interaction is a
real starting point.

Whatever table you supply, the report states how much of it touches your
protein — how many compounds have any reading for it, and how those
split between raising it, lowering it and moving it too little to count.
That is what separates a measured negative from a protein the table never
mentions, and the two read identically without it.

## Which screen the run designs

`--design-screen` writes the experiment the analysis prioritises. The
primary assay depends on what was declared about the alleles, because the
two possible assays measure different things.

By default it designs a **promoter reporter**: the scanned window cloned
upstream of luciferase or GFP, with a scrambled-motif variant of the top
element as the control that decides whether a drug effect runs through
that element at all.

When the caller declares that the PROTEIN is the problem — `--scenario
hypomorph` or `biallelic_hypomorph`, `--has-ptc`, or
`--broken-allele-makes truncated` — it designs a different screen: the
**FDA-approved drug library against the gene's protein in the patient's
own cells**. A reporter measures transcription from a cloned fragment, so
it cannot see a stabiliser, a splice modifier or a readthrough agent, and
those are the mechanisms in play when an allele makes an impaired
product. The protein screen is also mechanism-agnostic, so unlike every
ranked table in the report it does not depend on this package's
inferences being right.

The controls in that design are not boilerplate. For a gene whose protein
tracks the cell cycle, a compound that arrests cells raises the signal
without touching the gene, so the plan insists on a cell-cycle co-stain
in the same well, a housekeeping normaliser that is not itself
cycle-regulated, and a cytotoxicity counter-screen in parallel. Those are
the same confounds the proliferation gate measures in the LINCS data, met
at the bench instead.

The endpoint is chosen from a declaration and never inferred, the same
rule the dosage scenarios follow: nothing here can read a variant.

## What the literature says moves the gene

CTD records which chemicals have been published as changing a gene's
expression. The run already read that list to decide which candidates were
novel. It now also says what KIND of chemical is on it: carcinogens, tumour
promoters, mitogens, inflammogens, industrial chemicals, nutrients, or
approved drugs. When most of what raises a gene is not a medicine, the gene is
more likely responding to injury than being regulated in a way a drug could
copy, and that is worth knowing before reading a ranking. The classes come
from a keyword vocabulary over the chemical names, not from a toxicology
database: an unmatched chemical is unclassified, which is not the same as
safe, and a chemical can be both a medicine and a hazard.

## Regulators with no drug

The regulator route lists approved drugs at factors recorded as regulating the
gene. It used to drop a recorded regulator that no approved drug reaches,
which left the section listing whatever happened to be druggable rather than
what is recorded as running the gene. Those factors are now listed separately.
They generate no route — there is no molecule to offer — and they are the list
of what a chemistry effort would aim at.

## Interpreting results

Approval filtering depends on source fields and curated entries. For example,
ChEMBL's maximum phase is used as an approval proxy. The package does not perform
a current FDA/EMA authorization check for each indication or jurisdiction.

The catalogue contains curated priors and broad mechanisms in addition to
gene-specific evidence. A populated table can coexist with no convincing promoter
signal. A “novel” label means no link was found in the sources queried, not that
the literature contains none. A comparison with an incumbent uses selected
heuristic axes and does not demonstrate superiority.

Optional-service failures can resemble absence of evidence in the JSON. Read
the logs. Match input assemblies, record dataset versions, inspect the selected
transcript, and assess whether evidence applies to the intended cell type before
using any candidate to plan further research.
