# Improvement review

[Documentation index](index.md)

This review was written against the repository as inspected on 2026-09-30.
Items that have since been done are listed at the bottom with what changed;
everything under the headings below is still outstanding.

## Highest priority

1. **Record the source release each lookup answered from.** Most of this item
   is done. Stages report four states (not requested, FAILED, empty, ok); the
   result carries the package version, a timestamp, every argument and a hash
   of every input file; and the two dosage tables now record their release in
   the result and in the report, read from the file rather than configured -
   gnomAD from its filename, falling back to the column spelling, and ClinGen
   from the date on its header. Where a file does not state a release the field
   says so instead of guessing, because a wrong release is worse than none.

   What remains is the networked sources. A ChEMBL or Open Targets response
   carries no version, so two runs a year apart are still indistinguishable in
   the saved result even though the underlying data moved. That needs either a
   per-service version endpoint where one exists or a recorded request date per
   source, and it is a larger change than the local files were.
2. **Audit curated biological claims and scoring.** Source-specific approval
   filters, tissue rules, clinical comparison labels, and outcome “probabilities”
   need a documented evidence trail and independent benchmarks. Successful
   software tests cannot establish scientific validity.

## Maintainability and release work

- Run the full test suite on Linux before pushing. The two CI failures above
  were found without access to the CI logs, by testing a clean clone on
  Windows under Python 3.10 and 3.12, compiling every file under 3.10 and
  3.14, and running individual stdlib-only modules under WSL. That found
  everything that was visible, and it cannot show a Linux-only failure in a
  test that needs NumPy. The development machine is Windows and its WSL has no
  pip, so those tests have so far only run in CI. A virtual environment in WSL
  with the dev and chemistry extras would close the gap.
- Finish splitting `pipeline.py`. The report formatting moved to `report.py`
  and the file went from 4,649 lines to 1,493, but `run` is still 1,185 of
  those and does the fetching, the scoring and the assembly in one function.
  The test module is still one file covering every domain.
- Add a versioned result schema and transcript identifiers to multi-promoter TSV
  rows. Treat schema changes as compatibility changes for shortlist and explain.
- Add input validation to the remaining helper CLIs. The exported tables and
  the Windows launchers are done. TSV is written through `csv.writer`, so a
  name containing a tab or a newline no longer splits the row, and the BED
  name field is sanitised because that format cannot be quoted at all. Every
  launcher prompt now validates, and the validation itself was the hole —
  see below.
- Rename the bundled HPA lookup, and decide the redistribution question.
  Provenance and hashes are done: the result carries the package version, a
  UTC timestamp and, for every file argument, its size, modification time and
  sha256, plus the release for the two dosage tables. `datasets.py` is the
  inventory of all twenty-four tables with their sources, and the data guide
  states that this repository ships none of them and that a derived table
  inherits from its inputs.

  What is left needs a person rather than code. The terms are not uniform
  across these sources and change between releases, so a written inventory
  of them would need checking and re-checking — and a licence stated wrongly
  is acted on, where one left unstated is looked up. And `data/hpa.tsv` is
  still named as though it were expression evidence when it is a tissue
  lookup; renaming it is a deliberate data-maintenance change, not a
  refactor.
- Choose a software license before distributing a public release. No project
  license was present during this review; selecting one requires the owner's
  decision and is not inferred from the dependency licenses.

- ~~Add periodic opt-in service contract checks.~~ **Done**, and it found
  a live bug on its first run — see below. `python -m
  promoter_tf_drugs.services` probes each web service with a query whose
  answer cannot disappear (TP53, aspirin) and checks that the paths the
  parsers actually walk are present and populated. Four outcomes, never
  two: `ok`, `CONTRACT CHANGED`, `unreachable`, and `PROBE BROKEN` — an
  outage is a reason to try later, a schema change is a reason to edit a
  parser, and a broken check says nothing about the service at all.
  Probes call the package's own request functions rather than composing
  their own, because a probe that builds its own parameters checks a
  fiction.

  The other two parts of this item were already done. Progress reporting prints a checkpoint per stage with
  elapsed time and what it found, and `--max-requests N` caps how many
  requests a run may send. The budget counts requests that leave rather
  than lookups, so a cache hit stays free and a rerun resumes; retries of
  one request count once; and running out raises `HttpError`, the same
  exception a dead host raises, so every `_optional` path degrades as it
  already does for an outage rather than meeting a new failure mode. The
  report says when the ceiling bound the run, above the evidence table
  rather than below it — a cap that silently truncates evidence is worse
  than no cap, because every table underneath looks identical whether a
  source had nothing to say or was never asked.

  What is left is the service contract check: nothing verifies that a
  response still has the shape the parser expects, so a renamed field
  degrades into an empty result rather than an error.


## Tried and dropped: a cell-cycle-shift check from the LINCS files

The idea was to warn when a compound raises a mitotic gene only because it
arrests cells in G2/M, by asking how much of a G2/M gene panel its signatures
move the same way. It cannot work with the local `.gmt` files: even for
nocodazole, paclitaxel, vincristine and tozasertib the median share of a
54-gene panel in the up lists is zero, and the highest of 2,207 compounds is
two genes. The lists are short and most of the genes are inferred rather than
measured. Checking the premise directly for BUB1B also found no pattern: 1 of
468 mitotic-disruptor signatures has BUB1B in its up list (0.21%), against 19
of 32,665 others (0.06%), and only 20 compounds in all of LINCS raise it. What
was kept instead is the principled part: a curated tag for compounds that
disturb mitosis, scored down for mitotic genes, and a plain caution that a
mitotic transcript moves with the phase mix of the culture.

## Proposed: reading the conformation beyond nuclear receptors

This is a scope, not a change. It was written on 2026-10-08 from queries
against GPCRdb and KLIFS, and the numbers below are what they returned that
day.

**The gap.** `--state` picks the agonist-bound or antagonist-bound structure,
but the package can only tell which is which for nuclear receptors, from a
bound coactivator peptide. A GPCR, a kinase or an enzyme carries no such
peptide, so every structure comes back `unknown`, a request for a particular
state is refused, and the only option is `--state any` with the direction left
unclaimed.

**Reading the state is a lookup, not geometry.** Two public databases already
classify each structure, and both answered from this network.

- GPCRdb labels every GPCR structure Inactive, Active or Intermediate. CCR5
  has 11: seven inactive and four active.
- KLIFS labels every kinase structure by DFG (in, out, out-like) and by the
  position of the αC helix. EGFR has 565, covering all the usual combinations.

So the adapter is small: ask the database for the structure's state, record
the source, and treat a structure the database does not know as `unknown`, as
now. The state is then read from a curated source rather than from the file,
and the result should say so.

**What the lookup does not solve.**

- *The box needs a ligand.* The package places the box on a co-crystallised
  ligand. CCR5's four active-state structures carry a chemokine or nothing, no
  small molecule, so there is nothing to place a box on. Docking into the
  active state would need the box carried over from an inactive structure by
  aligning the two pockets, which is new work with its own error. CCR5 also has
  no small-molecule agonists to validate against: maraviroc is an antagonist.
- *Coverage varies a lot.* BUB1B has one human structure in KLIFS (6TLJ, DFG
  out-like, αC out) and it holds no ligand, so docking would refuse for the
  same reason it does today. BUB1B's kinase domain is also reported to be
  catalytically inactive (it is a pseudokinase), which makes "active versus
  inactive" a weak idea for it. Check that before relying on it.
- *Direction cannot be derived.* For a nuclear receptor the wanted state comes
  from the factor's role and your goal. A GPCR has no such role annotation, so
  the wanted state would have to be stated by the user.
- *The gate needs known ligands per state.* A pocket is validated against
  agonists for the active state and antagonists for the inactive one, with at
  least five of each. Active-state GPCR pockets will often refuse.

**By class.**

| Class | What to build | Worth it? |
|---|---|---|
| Kinases | KLIFS adapter; dock into DFG-in and DFG-out structures separately | Best: hundreds of structures for many kinases, and the DFG state separates type I from type II inhibitors |
| GPCRs | GPCRdb adapter; box transfer for active-state structures without a ligand | Second: inactive-state docking already works, the gain is the active state, which is the part with the weakest data |
| Enzymes | Nothing new | There is no on/off pocket; an inhibitor fits the active site, and the existing gate already tests whether docking recovers known inhibitors |

**Recommendation.** Do kinases first, GPCRs second, and leave enzymes alone.
Each adapter should return the existing `Conformation` record with a `source`
field, keep the refusal for a structure that cannot be classified, and add a
test for each of: state known, state unknown, database unreachable (which must
read as unreachable, not as "unknown state").

## Done since this review was written

**Global expression state removed.** `delivery._LOCAL_EXPRESSION` survived
between `run` calls, so a second analysis in one process that asked for no
expression file silently used the first one's table and reported the source
as present. `reset_expression` is now called at the start of every run, and
the "disabled for this run" flag resets with it rather than outliving the run
it names. The provider-per-run shape the item asked for would still be
cleaner; the leak it existed to prevent is closed.

**Direction prediction wired into the docking screen.** A screen returned
compounds that fit a pocket and said nothing about whether any of them turned
it on or off. The per-receptor classifier in `agonism.py` is now pointed at
the hits, and reads its training labels from a local drug table's
mechanism-of-action column as well as from ChEMBL, so it works offline. It
annotates and never ranks, and refuses on molecules outside its applicability
domain rather than calling them weakly.

**Report formatting split out of `pipeline.py`** into `report.py`, one way:
`report` imports nothing from `pipeline`.

**Several reporting faults found by running the thing.** A failed evidence
source printed a bullet with nothing after the dash, because `IMPACT` covered
seven of the fourteen sources a real run records. The panel compared the top
route against a runner-up from a different evidence tier and reported the
resulting negative margin as a tie. The funnel counted removals in distinct
drug/target pairs and unique kills in rows, so a stage could report more
unique kills than kills. The docking report asserted a cause for dropped
compounds it had not measured. Each has a regression test.

**Impossible flag combinations are refused before the run starts**, rather
than twelve minutes in, after the scan and every lookup have already run.

**The panel forwarded twenty-one of the fifty-seven arguments it parses.**
`--genetics`, `--lincs`, `--upstream-network`, `--tf-activity` and
`--docking-results` were accepted, listed in `--help`, and dropped. Both CLIs
now derive the arguments from `run`'s own signature.

**Static guards for three bug classes**, each checked by reintroducing the
bug: a local read before assignment, a CLI flag parsed but never read, and a
documented command using a flag its parser does not accept.

## Done in this round

**The result now says what produced it.** The query block records all
fifty-eight arguments rather than twenty-four, and a `provenance` block
carries the package version, a UTC timestamp, the Python version and,
for every file argument, its size, modification time and sha256. Hashes
are skipped above 64 MB - the ReMap BED is 1.4 GB - and the record says
so, because a blank reads like a file nobody looked at. The result also
carries a schema version, and `shortlist` and `explain` warn when they
are handed a file written by a newer package instead of reading it
blind.

**A failed CTD lookup no longer reads as a novelty claim.** The outage
in the validation record is gone - the endpoint returns 200 and 153
curated chemicals for CCR5 - but checking it exposed something worse
than an outage. `ctd_expression_links` returned `{}` when the request
failed, which is the value it also returns when CTD answers and knows of
nothing, so every candidate came back `previously_linked: False` and the
report printed a bold "no". The function's own docstring calls that "the
one genuinely embarrassing failure mode". Failure now returns `None`,
and the ledger distinguishes a search that found nothing from one that
never happened.

**Advertised evidence matches the queries.** `--signature-file` says in
its own help that it does not run a connectivity search, and the
trans-eQTL lookup that used to be described as advice is now performed.

**Windows launchers meet all four criteria** - local environment
preferred, `curl -fL` so an HTTP error cannot be saved under the name of
the file you wanted, content checks rather than size, and validated
input. Writing them turned up a bug of its own: four were saved with
Unix line endings, which `.gitattributes` forbids for exactly the reason
that bit - `cmd.exe` mis-parses an LF-only batch file, and a valid path
was reported missing while an invalid preset was accepted.

**Reports say what an incomplete test cost, in plain words.** "requested
source(s) failed: promoter_dna" named a module and told the reader
nothing. Every one of the fifteen sources now has a sentence a person
can act on, and the heading says these are tests that could not be
completed - which is not the same as their having returned a negative
result.

**A run says what it is doing while it does it.** One checkpoint per
stage, to stderr so a piped report is untouched, with elapsed time and
what the stage found rather than what it was about to try.

**How much the gene has to move is now stated.** `magnitude.py` takes
the shape of the problem - haploinsufficiency, a hypomorph, a toxic
gain, a null, a silenced allele, or a gene with no mutation at all - and
turns it into a requirement. Two of those shapes cannot be solved by
anything this package ranks, and the report says so above the tables
rather than offering two hundred routes for a problem none of them
addresses. `constraint.py` adds the companion question, from gnomAD and
ClinGen: is this gene's dose under selection, and is too much of it
independently harmful.

## Found by running it, not by reading it

None of these was in the review. Each has a regression test.

- The docking screen could not reach its own binary: `_docking_binary`
  returns a WSL path and two call sites handed it to a Windows
  subprocess, which the caller caught at DEBUG and recorded as `None` -
  so a whole library would have scored nothing and reported an empty
  screen rather than a broken one.
- Pocket selection chose a sulfated tyrosine in a ten-model NMR peptide
  over a maraviroc-bound crystal structure, for three independent
  reasons.
- The panel parsed thirty-six of fifty-seven arguments and forwarded
  none of them, recording `use_genetics: False` after being asked for
  `True`.
- The panel compared the top route against a runner-up from a different
  evidence tier and reported the resulting negative margin as a tie.
- The funnel counted removals in drug/target pairs and unique kills in
  rows, so a stage could report more unique kills than kills.
- An evidence source that failed printed a bullet with nothing after the
  dash, because the impact table covered seven of fourteen sources.
- The Orange Book notice printed once per drug rather than once per run.
- A refusal that happened before the screen started wrote no file, so a
  target that was tried and refused looked identical to one never tried.
- The launcher's own input check was the injection it was meant to
  prevent. `echo !GENE!| findstr /r /c:"..."` does not protect the
  pipeline: with a gene of `CCR5 & del *` the line runs as `echo CCR5`
  and then `del *| findstr ...`, so findstr reads del's output,
  errorlevel comes back 0, and the value is accepted **because** it
  contained a separator. Delayed expansion does not cover a pipeline;
  quoting the echoed value and the pattern does. Four checks across two
  launchers, and the gene is the one that matters — both take it from
  `%~1`, so it can arrive from a shortcut or a dropped file.
- ClinGen dosage curation loaded zero genes from the real file. The
  loader dropped every line starting with `#` before looking for a
  header, and in the real file the header *is* a `#` line. So the header
  was discarded, the first gene was read as one, no column matched, and
  it returned `{}` for all 1,532 curated genes without a word. The
  fixture it was tested against had been written to match the parser
  rather than the file.
- ClinGen's score was read as a severity scale. Two of its values are
  category codes above the grade's range: `40` means "dosage sensitivity
  unlikely" and `v >= 2` made it True, so the thirty-five genes most
  confidently curated as *insensitive* were reported as
  haploinsufficient — PCSK9 among them, which is the README's own worked
  example. And `0` means "no evidence available", 1,306 genes, which was
  read as a measured negative. Both errors ran the same way: turning an
  absence into a statement.
- gnomAD was read one transcript at a time and the lowest LOEUF kept.
  That is a minimum over about a dozen draws per gene, so it lands below
  the typical one and every multi-transcript gene drifts toward
  "constrained". 9,207 of 18,204 genes differ from their MANE Select
  value, 922 by more than 0.2 LOEUF — wider than the gap between this
  package's own sensitive and tolerant lines, so it changed verdicts and
  not just numbers.
- Where gnomAD and ClinGen disagree, only one was reported. 44 genes are
  curated haploinsufficient while carrying about as many loss-of-function
  variants as chance predicts, and the reason is usually informative:
  late onset, incomplete penetrance, or a phenotype that is not harmful.
  PCSK9 is the clean example — its loss has a phenotype, which is what
  ClinGen records, but the phenotype is low LDL, so selection never
  removed the carriers.
- "No variance" meant "a standard deviation of exactly zero", which is
  a narrower thing. A column constant to every digit anyone would print
  can still have a std near 4e-17, and `_pairwise_r` divided by it and
  reported r = -0.134 for a profile carrying no information — breaking
  its own docstring, which promises that a correlation that could not be
  computed is not a correlation of none. The same rule was written a
  second time thirty lines below as `.std() > 0`, around a copy of the
  same formula, which is how the two spellings came to disagree. Now one
  `_pearson`, one flatness test, both callers. The tolerance is relative
  to the values' own magnitude and sits about nine orders of magnitude
  below the spread of any real log fold change column, because a looser
  one would start discarding small genuine effects — a change to results
  rather than to robustness.
- The Python 3.10 jobs failed from 2026-10-01, on both operating systems,
  and nothing on a development machine could show it. `predictions.py` put a
  multi-line expression inside a single-quoted f-string, which is legal from
  3.12 (PEP 701) and a SyntaxError on 3.10, so the package did not import
  there and every test errored before it ran. `ast.parse(feature_version=...)`
  cannot find this, because it does not change how f-strings are tokenised;
  only a real 3.10 interpreter rejects the file. Found by compiling every
  file under 3.10. The workflow now does that as its first step, with the
  file and line in the output, and a tokenizer-based test catches the same
  constructs on 3.12 and later.
- A second CI failure affected every Linux job whatever the Python version.
  `wsldock.to_wsl_path` turned `C:\Users\x\y.pdb` into `/mnt/c/Users/x/y.pdb` by
  calling `Path(...).resolve()` and reading the drive off the result, which
  only works where `Path` is a WindowsPath. On Linux a backslash is an ordinary
  filename character, so the whole string is one relative name, `resolve()`
  prefixes the current directory, and the function returned
  `/home/runner/work/.../C:/Users/x/y.pdb`. In production this never mattered,
  because the function only runs on Windows to reach a Linux binary through
  WSL, but it is a pure string conversion and its test failed on Linux. Found
  by running the real function under WSL's Python, which needs only the
  standard library for this module. The conversion for other hosts is now its
  own function, `_from_windows_path`, so a test can call it on any machine and
  not only on the one it was broken on. Windows behaviour is unchanged.
- The test that checks the 3.10 floor could not run on 3.10. It read
  `pyproject.toml` with `tomllib`, which exists from 3.11, so on the
  interpreter it was written about it failed with ModuleNotFoundError.
- The command reference depended on the shell it was generated in. argparse
  wraps help to `$COLUMNS`, and the generator never set it: unset or 80
  reproduced the committed file exactly, 120 rewrote 320 lines, and 200
  rewrote 208. The CI step that fails on `git diff` would have reported stale
  docs that were not stale. The generators now pin the width and write LF.
- Four failures in the first attempt to reproduce CI were not CI failures.
  They were `WinError 3` from cache paths over Windows' 260 character limit,
  caused by running the test from a deeply nested scratch directory. The
  development guide now says so, because the next person to reproduce CI on
  Windows will meet the same four.
- `_is_mane` had never once returned True against the live API. Found
  by the new contract check on its first run. `select_transcripts`
  prefers MANE Select, falling back to Ensembl canonical, and `_is_mane`
  tests four spellings of the flag — but the lookup never sent `mane=1`,
  so Ensembl never returned any of them. Every gene silently took the
  canonical branch. Canonical is a reasonable answer and it is not
  always the same answer: where the two differ, the package scanned a
  different promoter than it believed it was scanning. Both halves are
  fixed — the request asks, and the parser reads the fifth spelling
  Ensembl actually sends.
- `run_everything.bat` could not reach the newest work. It never passed
  `--scenario`, so the requirement section was never computed, and it
  downloaded neither gnomAD nor ClinGen, so dosage was never assessed.
  Every command it issued parsed and every flag was accepted — the batch
  file was correct and simply predated the feature, which is the failure
  mode a launcher has that a library does not.
- A workup merged the previous workup's factor docking into this one.
  The launcher clears `product.json` and `factors.txt` before each run
  and not the per-factor results, which are named `<FACTOR>.json`. The
  set of factors changes with the scan and `merge_docking.py` globs the
  directory, so a stale `STAT3.json` joins today's `SP1` and `NFKB1` and
  reaches the analysis as one result. It matters because a docking run
  does not mean the same thing every time: the conformation is derived
  from the direction and the factor's role, so lowering a gene through
  an activator docks that activator's INACTIVE pocket, and a previous
  `up` run leaves a screen of the opposite pocket under exactly the
  filename a `down` run would write. Found by noticing two results in an
  output directory dated three hours before the run writing beside them.
- Nothing in a docking result said what the run had been asked for, so
  the above was undetectable after the fact. A result recorded its
  target, its pocket and whether the conformation could be read, but not
  which conformation was REQUESTED and not when the run happened.


## Changes completed in this update

Explicit package discovery and dependency extras; command entry points for saved
result tools; input validation before requests; consistent report limits; correct
BED strands for negative-strand genes; zero-valued shortlist overrides; cache
separation by request identity; concurrent cache writes; bounded POST redirects;
robust retry delays; regression tests; documentation and contributor instructions.
Local PFM loading now honors `--matrix-id` and rejects missing requested IDs.

Cache keys changed. Existing entries remain untouched but need to be populated
under the new keys before they can support offline runs.
