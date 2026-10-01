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
   rate. The scan-quality summary compares site counts with a rough null estimate.
3. Annotation identifies TF roles and nuclear receptors. Optional ChIP overlap
   provides another source of evidence. A motif match alone is not observed binding.
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
| `--upstream-network` | Signed paths to druggable upstream targets. | Sign and confidence depend on curated edges; longer paths incur penalties. |
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
