# Docking

[Documentation index](index.md)

Docking asks one question: can this molecule occupy that pocket. It is
worth running because it is the only step in the package that looks at a
structure rather than at an annotation, and it is worth distrusting for
the same reason — a score is a prediction, and most of this page is about
telling the two apart.

## What it can and cannot say

**It says a molecule fits.** Nothing more. It does not say the molecule
activates or inhibits, and that direction is what decides which way your
gene moves. A compound that fits the agonist-bound conformation of a
receptor is as likely to be an antagonist of it.

**A score is not an affinity.** Even a pocket that passes the validation
below pulls its known binders toward the top of a list without ordering
them usefully once they are there.

**An unvalidated pocket means nothing at all.** Which is why the package
refuses to report one.

## The validation gate

Before any library is screened, the pocket is tested on compounds whose
answer is already known: the receptor's own measured ligands, against
decoys drawn from the approved library. If docking cannot separate those,
its scores on untested compounds carry no information, and the run
reports the refusal instead of a ranking.

The bars are an AUC floor and a top-1% enrichment floor. A target with
more known ligands than the gate needs is sampled — sixty against ten
times as many decoys puts the standard error on the AUC near 0.02, and a
target with thousands would otherwise spend more docking on the gate than
on the library it is gating. The sample is drawn at random under the run's
seed, never by potency: taking the tightest binders would make the actives
easier to recover and raise the AUC, which is tuning the gate to pass.

## Which program

Two tiers, and the difference is not speed.

**gnina** docks the way smina does and then rescores every pose with a
convolutional network trained on crystal structures. That is a second and
different judgement about whether a pose looks real, and it produces two
numbers the empirical function cannot:

| number | what it says |
| --- | --- |
| `CNNscore` | how pose-like the geometry is, 0 to 1 |
| `CNNaffinity` | predicted affinity in pK, comparable between ligands |

The Vina-style `minimizedAffinity` is still produced and still reported.
The three disagree often enough to be worth seeing side by side.

**smina** is the fallback and the easier install: one static file,
empirical scoring only. **AutoDock Vina** also works and needs PDBQT for
both receptor and ligand.

Installation for all three, including the CUDA libraries gnina needs and
how it is reached from Windows through WSL, is in [inputs and
outputs](data.md#docking).

## Was it measured?

Yes, on a complex whose answer is in the file. PDB 1DB1 is the vitamin D
receptor with calcitriol bound. Redocking the ligand from a conformer
embedded from its own SMILES, so no crystal coordinate survives into the
input:

| rank | RMSD from crystal | affinity | CNNscore | CNNaffinity |
| --- | --- | --- | --- | --- |
| 1 | 1.04 Å | -12.79 | 0.959 | 8.56 |
| 2 | 1.77 Å | -11.79 | 0.693 | 7.95 |
| 3 | 9.65 Å | -10.02 | 0.628 | 7.91 |

The top-ranked pose is correct by the usual 2 Å convention and both
scoring functions chose it. The CNN separated it far more sharply: ranks
1 and 2 differ by 1.0 kcal/mol of affinity but by 0.27 of CNNscore, and
the CNN collapses to 0.63 for a pose 9.7 Å out of place that the
empirical score still rates -10.0.

That is a pose-accuracy measurement on one complex. It is not a
measurement of whether a screen finds real drugs — the AUC and enrichment
floors are for that — and it says nothing about direction.

## Which score should rank the hits?

Measured per pocket rather than assumed. The gate scores its known
ligands under every metric the engine produced, which costs nothing
because the poses are already scored, and reports what each was worth.

For CCR5 against its maraviroc-bound structure:

| metric | AUC | enrichment |
| --- | --- | --- |
| CNNaffinity | 0.967 | 11.31 |
| affinity | 0.937 | 9.69 |
| CNNscore | 0.704 | 9.69 |

The CNN's affinity prediction separated real binders better than the
empirical function. CNNscore did much worse, which is what it should do:
it rates whether a pose looks real, not how tightly a molecule binds.

Ranking still uses the affinity, so there is a fixed thing to compare
against, and the note in the report says plainly that this is one
pocket's ligands rather than a general finding about scoring functions.

## Which pocket gets screened

The box is placed on a co-crystallised ligand, because that is the only
reliable marker of where the orthosteric site is. Choosing that ligand
has three traps. All three were hit by CCR5 and none of them fails
loudly — the run completes and the scores describe another site.

**Modified amino acids are not ligands.** Sulfotyrosine, phosphoserine
and selenomethionine appear as HETATM records because they are
non-standard, not because anything is bound. CCR5's N-terminal structures
carry sulfated tyrosines and one was selected as the ligand.

**NMR entries hold many models of the same atoms.** Grouping by residue
pooled all of them, so a sixteen-atom residue in a ten-model file counted
as a hundred and sixty atoms with a centre averaged across ten
conformations — a point no model occupies. Largest-group-wins then puts
every multi-model structure ahead of a real crystal ligand.

**Crystallisation lipids are usually the largest group.** A membrane
protein solved in a lipidic cubic phase is full of monoolein, which lies
against the outside of the protein. These are ranked last rather than
excluded, because for a lipid receptor the fatty acid is the genuine
ligand; where nothing else exists the choice is made and flagged.

## Conformation, and when it cannot be read

`--state` picks the agonist-bound or antagonist-bound structure, because
docking the wrong one answers a different question. The conformation is
read from bound coactivator and corepressor peptides: an LXXLL peptide
can only sit in the AF-2 groove when helix 12 is in its active position.

That signal is specific to nuclear receptors. A GPCR, a kinase or an
enzyme carries no such peptide, so every structure comes back `unknown`
and a run asking for a particular conformation is refused. The refusal
distinguishes the two cases, because they are not the same claim:

- a structure **classified as the other conformation** is a finding
- a structure that **could not be classified** is not a finding at all,
  and reporting it as one would turn a measurement that was never made
  into a negative result

`--state any` docks whichever pocket exists and records the conformation
as unverified. The result carries `conformation_verified: false` and a
caveat saying that binding and direction are different questions and only
the first was addressed. For a GPCR this is the only way to screen at all.

## Running it

A screen is thousands of independent jobs, and a docking program spreads
a single job across cores with heavy diminishing returns. So `--jobs` is
the speed lever, not `--cpu`:

| setting | per compound | 3,315 compounds |
| --- | --- | --- |
| `--cpu 16` | 11.4 s | 10.5 h |
| `--jobs 8 --cpu 4` | 2.0 s | 1.8 h |

Scores are identical either way — same receptor, same box, same seed, and
the only change is how many run at once. `--jobs` defaults to 1.

```sh
python -m promoter_tf_drugs.docking --gene VDR --state active \
    --jobs 8 --cpu 4 --library data/approved_library.tsv \
    --out results/vdr/dock.json
```

The gene's own protein is a target too, with `--product`: the
conformation then follows `--direction` directly, because no factor sits
between the product and the gene.

```sh
python -m promoter_tf_drugs.docking --product --gene CCR5 \
    --direction down --state any --jobs 8 --cpu 4 \
    --library data/approved_library.tsv --out results/ccr5/dock.json
```

Screens report progress and a rate as they go. A run that printed nothing
until it finished looked exactly like a wedged one.

## Folding it back in

Docking usually runs on a different machine from the one with the network,
so the two halves are separable: `--prepare` does every networked step and
leaves a folder, `--bundle` docks it with no network at all, and
`--docking-results` brings the JSON back into an analysis.

```sh
python -m promoter_tf_drugs --gene CCR5 --direction down \
    --docking-results results/ccr5/dock.json --out results/ccr5
```

Hits join the route tables at tier 4, as `TF: predicted pocket binder`, or
as a direct route when the target was the gene's own product. A route
whose structure had an unreadable conformation says so rather than
implying the question was never asked.

## What gets dropped, and why

Compounds absent from a ranking are counted by cause rather than
explained by the likeliest one. A molecule that could not be given a 3D
structure and one whose docking crashed are different problems and send
you to different places. For the CCR5 library, 20 of 3,315 failed, all of
them at conformer generation — macrocycles and peptides such as
daptomycin, oritavancin and pibrentasvir, which ETKDG cannot build.

## Salts and counter-ions

About a third of the approved library is stored as a salt: the drug with
chloride, sodium, calcium, sulfate, maleate or another counter-ion, and
sometimes two copies of the drug. By default each structure is docked exactly
as stored, so the counter-ion is scored as part of the ligand. The salt
dissociates in solution and does not bind the receptor, and the extra atoms
add contacts, charge and size to the score. A bare ion such as `Ca+2` also
draws RDKit's `UFFTYPER` warning, because it has no force field type.

`--library-limit` takes the first N compounds in the order the table is
written, and the bundled table is alphabetical. A limit of 600 docks the drugs
from 2-mercaptoethanesulfonic acid to chlorpropamide and nothing after, so it
is the start of the alphabet, not a sample. A CCR5 screen run that way never
docked maraviroc, which sits at position 1,782 and is both the approved CCR5
drug and the ligand in the pocket being docked. The launcher now docks the
whole table, and a limit that drops compounds logs a warning saying so.

`--strip-counter-ions` docks the drug alone, with one copy of each remaining
component. A piece is removed when it has no carbon, when RDKit's default
salt list recognises it, or when it is on a short list of organic counter-ions
in `docking.py` (citrate, succinate, besylate, meglumine and similar). Charges
are left as stored. A piece that can also be the drug, such as piperazine,
lysine or caffeine, is never removed on its own: piperazine citrate becomes
piperazine, and an entry with two such components stays whole and is counted
as ambiguous. On the bundled library this changes about 1,130 of 3,315
entries and leaves about 40 ambiguous.

The option applies to the validation gate as well as the library, so the
known ligands and decoys are prepared the same way as the drugs they are
compared with. The saved result gains a `ligand_prep` block with the number of
entries changed and ambiguous; it is absent when the option is off. Scores
from the two modes are not interchangeable, so compare runs only within one.

## A caution about size

Docking scores grow with molecular size: a larger molecule makes more
contacts with any pocket. In the CCR5 screen the top hits by affinity
included ergotamine, plicamycin, digoxin and deslanoside, with a median
of 55 heavy atoms against maraviroc's 37. Both the empirical score and the
CNN show the correlation.

No size correction is applied. Normalising a score to make a ranking look
better is tuning, and the package does not do it — the caution belongs in
how you read the order.

## Exploratory docking on an AlphaFold model

The normal stage refuses a target whose pocket cannot recover its own known
ligands. For a gene with an AlphaFold model, a crystallised relative, and no
structure of its own with a ligand, `exploratory.py` asks a weaker question:
does anything in the library fit a model of the pocket.

```bash
python -m promoter_tf_drugs.exploratory --gene BUB1B --uniprot O60566     --template-pdb 4QPM --template-chain A --template-ligand ADP     --jobs 8 --cpu 4 --reference
```

It downloads the current model, aligns the template to it on sequence,
superposes the aligned core, and moves the template ligand into the model's
frame to place the docking box. It refuses when fewer than 100 residues align,
identity is under 25 percent, the superposition core is more than 3 A off, or
the pocket residues average below pLDDT 70. Residues below pLDDT 50 are left
out of the receptor. The model as a whole may be mostly disordered, as BUB1B is
(49 percent confidently folded); only the pocket is judged.

Compounds are scored, then ranked by how far each score sits from what its size
predicts. Only compounds between half and 1.6 times the template ligand's heavy
atoms are ranked, because a line fitted through every size is pulled around by
molecules that cannot occupy the pocket. Everything keeps its raw score in
`exploratory.json`.

The class check is the safeguard. Approved kinase inhibitors in the library
should rank above the rest. If their size corrected AUC is under 0.6, the
output says the check failed and the ranking is not a list of what fits.
`--reference` docks the same compounds into the template's own pocket, which
shows whether a failure comes from the model. Neither check says anything about
the target specifically.

On BUB1B with 110 sampled compounds, kinase inhibitors scored worse than chance
in the AlphaFold pocket (AUC 0.34) and well above chance in the open BUB1
pocket (0.68 size corrected, 0.80 raw). Eleven of the 27 ADP atoms land inside
model atoms. The model's pocket is closed to these compounds, and its motif
reads DFS where BUB1 has DLG. A ranking from it is not usable.

Output is `exploratory.json` and `exploratory.md` under `--out`. Every result
carries the banner that it is unvalidated. The main pipeline does not call this
module, and nothing from it enters a candidate table.
