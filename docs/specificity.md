# Specificity, and what counts as evidence

[Documentation index](index.md)

Two questions decide whether a result here is worth acting on. *Did the
evidence actually arrive?* And *how much does this route move besides your
gene?* This page covers both.

## Did the evidence arrive?

Every optional analysis degrades gracefully — a dead API costs you that
section, not the run. That is right for robustness and wrong for inference,
because the scoring layer cannot distinguish these:

```
the lookup ran and found nothing       -> evidence of absence
the lookup broke and returned nothing  -> absence of evidence
```

Both arrive as an empty dict. A run where ChIP, tissue expression and LINCS
have all silently died still prints grades, Pareto fronts and a tier-1 row,
and is **structurally identical** to a run where all three worked.

So every report opens with an **Evidence available in this run** table:

| state | meaning |
|---|---|
| `ok` | ran, and its result may be used as evidence |
| `FAILED` | requested and did not deliver. Named, with what it cost you |
| `not requested` | you did not ask for it. A choice, not a fault |

Two places this changes a number rather than just a message:

**A broken ChIP lookup no longer scores as measured absence.** ReMap's API
response shape has changed between releases; when it does, a promoter comes
back with one or two factors instead of dozens. That used to award every
factor a hard zero on the 20-point evidence axis — presenting a failed
request as though binding had been measured and ruled out. Below three
factors over a promoter-sized region the result is now refused as scoring
input. A local `--chip-bed` that parses and finds nothing is the opposite
case: a real negative, and it is kept as one.

**An unavailable expression source is separated from a tissue-restricted
gene.** Both produced "unknown" for every factor before.

### The run verdict

Two things must hold before a ranking is evidence rather than a reading
order: the scan found signal above chance, **and** at least one source of
measurement worked. Neither suffices alone — an enriched scan with no ChIP is
still unconfirmed sequence, and perfect ChIP over a noise scan confirms
binding at sites you had no reason to look at. When both fail the report says
**NOT INTERPRETABLE** at the top and explains why.

### False discovery proportion

`scan_quality` reports the share of the hit list expected to be chance,
because the enrichment ratio badly understates it. **1.2x enrichment sounds
like mild signal; it means roughly 83% of the rows below are noise.** The
report prints both, plus how many sites are expected to be real — while
stating that nothing in the sequence says *which*. Only ChIP separates them.

A weak scan is usually fixed by the window rather than the threshold: a
±1 kb core-promoter window at `--pvalue 1e-5` is interpretable, whereas
5 kb/1 kb at `1e-4` expects ~740 chance hits across ~880 matrices.

## Four axes of collateral

A drug can be perfectly selective and still move far more than your gene.
These are scored separately because they fail independently.

### 1. Regulon breadth — how many genes the factor controls

A compound binding exactly one target, 1 nM, zero off-targets still moves
tens of thousands of promoters if that target is SP1. Supply real counts with
`--regulon-file` (TRRUST v2).

**The trap: regulon size is confounded with study effort.** TRRUST is
literature-curated, so a factor with 8 recorded targets is far likelier
under-assayed than genuinely narrow, and a hard cutoff would select for
obscurity and call it specificity. This is the same correction
`drug_scoring.evidence_depth` already applies to protein selectivity.
`regulon_evidence()` therefore marks counts below 15 as **understudied**, and
surfaces the provenance — `regulon file`, `bundled tier`, `unknown` — that a
single float collapses. `--trustworthy-regulon` requires a real count.

### 2. TF tissue breadth — how body-wide the factor is

**There is no drug-to-tissue map in this package and one must not be
invented.** The proxy is the *transcription factor's* expression breadth, and
the stage is named "TF tissue footprint" so the limitation is visible where
it is read.

A count needs a denominator: HPA has ~50 tissues, GTEx ~54, and the Open
Targets vocabulary varies per gene, so `n_tf = 30` means different things in
different runs. Without a local table `tissue_universe_size()` returns None
and the axis refuses to score rather than fabricating a fraction.

**For a nuclear receptor the proxy is withheld entirely.** mRNA presence is
not activity: VDR is transcribed almost everywhere while 1,25-D signalling is
not. Ligand availability sets the footprint of a ligand-gated factor, and
nuclear receptors were 201 of 216 routes in a real CCR5 run — so a bad proxy
there would bias the whole output.

### 3. Clinical burden — how many systems move

Dexamethasone is a clean, potent, highly selective NR3C1 agonist. It passes
*every* other filter here, and it acts on essentially every tissue, because
that is what the glucocorticoid receptor does. `selectivity` counts
**proteins**; `regulon` counts **genes**; neither counts **systems**.

Three signals combine: a curated pleiotropic drug-class table; the regulatory
flags (`WITHDRAWN`, `black-box`) read from the same strings `shortlist`
filters on; and — the one that catches dexamethasone — whether the *target*
is a global regulator, reusing `regulon.GLOBAL_REGULATORS` rather than
duplicating it. A drug absent from the table is reported **unassessed**,
never clean.

### 4. Receptor family — can a selective ligand exist at all?

An orthosteric pocket is under selection to bind the same endogenous ligand
across a family, so it is **conserved and hard to discriminate**. Thyroid
hormone binds THRA and THRB because both evolved to bind it.

Visible in a real CCR5 run: THRA *and* THRB both surfaced with
dextrothyroxine, liothyronine and levothyroxine scoring against both; RARA
*and* RARG both surfaced with tretinoin, trifarotene and alitretinoin
crossing between them. Two family members scoring well is **one bet, not
two**.

Why family size and not ligand count: **ligand count measures study effort**,
the same trap as regulon size. A factor with two known ligands usually means
nobody did the chemistry — and the causation runs backwards, since nuclear
receptors have many ligands *because* they have a pocket. Few ligands
correlates with *no pocket*, the hardest case rather than the cleanest.

A large family is a **risk, not a verdict**: subtype-selective agents exist
and were often developed precisely because the pan-agonist was unusable —
resmetirom for THRB, trifarotene for RARG. Where one exists it is named and
the drug is exempt. This axis is left off unless you ask for it, because
unlike the other three it is a structural prior rather than a measurement;
the test is a binding panel.

## Trans-eQTL: the strongest evidence available

Everything else infers. A motif says a factor *could* bind. A trans-eQTL says
natural genetic variation near it **does** change your gene's expression, in
human tissue, with exposure randomised by meiosis.

The eQTLGen significant-results file is 2 MB, so the lookup runs rather than
being recommended. A worked result: **CCR5 has no significant trans-eQTL at
all** — it is not among the 6,297 genes with one, across 31,684 whole-blood
samples at FDR<0.05. So the test returns nothing for every factor, and
checking them one at a time is wasted effort.

**Absence is weak here.** Trans effects are small and detecting them needs
enormous power, so a null overwhelmingly means underpowered rather than
unregulated. eQTLGen is also whole blood: a relationship operating in
hepatocytes or neurons will not appear however real it is.

## What else is in the promoter window?

Every check above asks whether a drug or factor affects other genes. This
asks the prior question: **is the sequence you scanned uniquely your gene's?**

Found in a real CCR5 run. The window `chr3:46,365,946-46,371,945` sits
entirely inside **CCR5AS**, an antisense lncRNA spanning
`chr3:46,363,984-46,407,184` on the opposite strand. **CCR2**, the
co-regulated paralog, is 5 kb beyond the edge — `--upstream 10000` would have
swallowed it.

An antisense overlap matters more than it looks: the direction logic
("activator bound here, so block it") assumes the factor acts on *your*
promoter. Where the window lies inside an antisense transcript, a site may
drive that instead, and antisense RNA can stabilise or repress the sense
gene — so the sign does not carry.

This checks **annotation, not sequence uniqueness**. BLAT the window if you
need to know it appears nowhere else. For CCR5 that distinction came out
well: CCR5 and CCR2 promoters share almost no sequence (12-mer Jaccard
0.0005, 6 shared 12-mers against 0 for a shuffled control) despite the genes
being paralogs — the risk was neighbouring transcription units, not
duplication.
