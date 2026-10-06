# Improvement review

[Documentation index](index.md)

This review is based on the repository as inspected on 2026-09-30. The items
below remain follow-up work; they are not claims of completed implementation.

## Highest priority

1. **Record stage status and provenance in results.** Optional-stage exceptions
   often become empty lists or dictionaries through `_optional`. Add a structured
   status per stage that distinguishes disabled, unavailable, empty, and complete.
   Record package version, inputs and hashes, source releases, parameters, and
   timestamps. This would make partial analyses and cached runs easier to audit.
2. **Remove global expression state.** `delivery._LOCAL_EXPRESSION` can survive
   between `run` calls in one Python process. Pass an expression provider per run
   so a prior file cannot influence a later analysis unexpectedly. Until then,
   use separate processes for analyses with different expression inputs.
3. **Align advertised evidence with implemented queries.** The main pipeline
   counts genes from `--signature-file` but does not invoke connectivity search.
   Genetics includes lookup guidance rather than completed regulatory validation.
   Implement those paths with fixtures and service tests if they are needed.
   The CLI help and new guides now state the signature behavior explicitly.
4. **Audit curated biological claims and scoring.** Source-specific approval
   filters, tissue rules, clinical comparison labels, and outcome “probabilities”
   need a documented evidence trail and independent benchmarks. Successful
   software tests cannot establish scientific validity.
5. **Consider other routes to identify similarities to FDA approved Drugs**
    in the docking strategy could we identify potential ligands/agonists/
    antagonists depending on the goal of the screen then compare similarity to
    known FDA approved ligands to guess on action
6. **guess on if a molecule will act as an agonist or antagonist based on structure**
   if the mechanism of a certain ligand is unknown, compare structures to structures with
   known interactions to make an educated guess if it would improve or dampen gene activity
   or expression.

## Maintainability and release work

- Split orchestration, report formatting, and export out of the large
  `pipeline.py`; split the large test module by domain while preserving coverage.
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

## Changes completed in this update

Explicit package discovery and dependency extras; command entry points for saved
result tools; input validation before requests; consistent report limits; correct
BED strands for negative-strand genes; zero-valued shortlist overrides; cache
separation by request identity; concurrent cache writes; bounded POST redirects;
robust retry delays; regression tests; documentation and contributor instructions.
Local PFM loading now honors `--matrix-id` and rejects missing requested IDs.

Cache keys changed. Existing entries remain untouched but need to be populated
under the new keys before they can support offline runs.
