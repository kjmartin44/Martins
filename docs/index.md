# Documentation

promoter-tf-drugs builds candidate routes from a gene's promoter to transcription
factors and drugs. Start with the synthetic demo, then run a gene with explicit
inputs and inspect the evidence behind each route.

- [Installation](installation.md): environments, extras, and Windows launchers.
- [Usage](usage.md): analysis, directions, offline runs, and saved-result tools.
- [Inputs and outputs](data.md): file formats, coordinates, and cache behavior.
- [Methods and limitations](methods.md): how routes are constructed and ranked.
- [Specificity and evidence](specificity.md): the four collateral axes, and
  how a failed lookup is kept from reading as a measured negative.
- [Comparing several genes](panel.md): the narrowing funnel and its control.
- [Docking](docking.md): what a docking score is worth, the validation
  gate it has to pass first, and how the pocket is chosen.
- [Command reference](cli.md): complete help for all five commands.
- [Python API](api.md): calling the package from Python.
- [How the package is laid out](architecture.md): the chain from gene to
  ranked table, and which of the sixty-one modules owns each step.
- [Development](development.md): tests, builds, and repository structure.
- [Validation record](validation.md): the requested live gene runs and caveats.
- [Improvement review](improvements.md): remaining work and its priority.

The command and API documentation describes the current implementation. Scores
and route labels are research heuristics; they do not establish efficacy.
