# Improvement review

[Documentation index](index.md)

This review was written against the repository as inspected on 2026-09-30.
Items that have since been done are listed at the bottom with what changed;
everything under the headings below is still outstanding.

## Highest priority

1. **Record stage status and provenance in results.** Optional-stage exceptions
   often become empty lists or dictionaries through `_optional`. Add a structured
   status per stage that distinguishes disabled, unavailable, empty, and complete.
   Record package version, inputs and hashes, source releases, parameters, and
   timestamps. This would make partial analyses and cached runs easier to audit.
2. **Align advertised evidence with implemented queries.** The main pipeline
   counts genes from `--signature-file` but does not invoke connectivity search.
   Genetics includes lookup guidance rather than completed regulatory validation.
   Implement those paths with fixtures and service tests if they are needed.
   The CLI help and new guides now state the signature behavior explicitly.
3. **Audit curated biological claims and scoring.** Source-specific approval
   filters, tissue rules, clinical comparison labels, and outcome “probabilities”
   need a documented evidence trail and independent benchmarks. Successful
   software tests cannot establish scientific validity.

## Maintainability and release work

- Finish splitting `pipeline.py`. The report formatting moved to `report.py`
  and the file went from 4,649 lines to 1,493, but `run` is still 1,185 of
  those and does the fetching, the scoring and the assembly in one function.
  The test module is still one file covering every domain.
- Add a versioned result schema and transcript identifiers to multi-promoter TSV
  rows. Treat schema changes as compatibility changes for shortlist and explain.
- Add input and output validation for helper CLIs and quoted TSV fields. Current
  export joins values directly; embedded tabs or newlines can break columns.
- Add explicit dataset provenance, hashes, and redistribution terms. The included
  HPA file is a lookup, not expression evidence, and should eventually be renamed
  or replaced through a deliberate data-maintenance change.
- Choose a software license before distributing a public release. No project
  license was present during this review; selecting one requires the owner's
  decision and is not inferred from the dependency licenses.
- Harden Windows launchers: prefer the local environment, quote and validate
  inputs, check download HTTP status, and validate file content instead of size.
- Add periodic opt-in service contract checks, bounded request budgets, and clear
  progress reporting. Some live analyses issue many sequential annotation and
  pharmacology requests. Cached responses help but do not eliminate that cost.
- Investigate the CTD HTTP 302 failures observed in all three requested live
  checks. Until that source works, absence of a retrieved link cannot support a
  novelty claim. See the [validation record](validation.md).

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

## Changes completed in this update

Explicit package discovery and dependency extras; command entry points for saved
result tools; input validation before requests; consistent report limits; correct
BED strands for negative-strand genes; zero-valued shortlist overrides; cache
separation by request identity; concurrent cache writes; bounded POST redirects;
robust retry delays; regression tests; documentation and contributor instructions.
Local PFM loading now honors `--matrix-id` and rejects missing requested IDs.

Cache keys changed. Existing entries remain untouched but need to be populated
under the new keys before they can support offline runs.
