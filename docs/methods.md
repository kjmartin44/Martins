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

Four checks can lower a route after it is built. None deletes anything, and
each says on the row what it did and why. They were added after a review of a
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

**A declared scenario is checked against the population data.**
`--scenario haploinsufficiency` asserts one allele is lost and the other is
intact. When the gene's pLI says losing one copy is not selected against, that
assertion disagrees with the frequency data: healthy people carry a broken
copy, so the condition is probably recessive and the second allele is probably
not wild type either. The report says so and points at
`--scenario biallelic_hypomorph`, which does not assume a wild-type allele.
The check is reported rather than refused, because the constraint tables are
read during the analysis and a gene can be haploinsufficient in one tissue and
tolerated in the population.

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
