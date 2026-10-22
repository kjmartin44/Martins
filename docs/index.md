# Documentation

generoutes collects candidate ways to move one gene's expression and
ranks them by the strength of the evidence behind each. The routes are
independent: the gene's own product, its recorded regulators, compounds
measured to move it, its promoter's predicted binding sites, and others.
Start with the synthetic demo, then run a gene with explicit inputs and
inspect the evidence behind each route.

- [How the pipeline works](how-it-works.md): one run from the question to the
  output files, and every database it reads and why.
- [Pipeline overview](pipeline.mmd): one run as a short Mermaid
  flowchart, with the routes in parallel rather than one chain.
  [The detailed flowchart](pipeline_detail.mmd) shows the branches, the
  shuffled-sequence null and the four collateral axes.
- [A full workup](workflow.md): the five steps, end to end, with the
  commands and what would make you stop at each one.
- [Installation](installation.md): environments, extras, and Windows launchers.
- [Running it from the command line](cli-guide.md): what is required, what has a
  default, and what each switch costs.
- [Usage](usage.md): analysis, directions, offline runs, and saved-result tools.
- [Setting up the data](data-setup.md): a walkthrough from an empty `data/`
  folder to a full run, with every download command and a check for each file.
- [Getting the analysis tables](datasets.md): every table the analysis can
  read, where to obtain it, and the command that says what you are missing.
- [Inputs and outputs](data.md): file formats, coordinates, and cache behavior.
- [Running with no network](offline.md): every data file, where it comes
  from, where to put it, and the three that need building first.
- [Methods and limitations](methods.md): how routes are constructed and ranked.
- [Specificity and evidence](specificity.md): the four collateral axes, and
  how a failed lookup is kept from reading as a measured negative.
- [Comparing several genes](panel.md): the narrowing funnel and its control,
  and the combined list that fuses docking with the rest.
- [Docking](docking.md): what a docking score is worth, the validation
  gate it has to pass first, and how the pocket is chosen.
- [Command reference](cli.md): complete help for all five commands.
- [Python API](api.md): calling the package from Python.
- [How the package is laid out](architecture.md): the chain from gene to
  ranked table, and which module owns each step.
- [Development](development.md): tests, builds, and repository structure.
- [Validation record](validation.md): what was run and checked, and what that
  does not show.
- [Example outputs](../examples/README.md): saved reports, one synthetic and one
  real run.
- [Checks on candidate drugs](checks.md): whether a person reaches the
  effective concentration, and whether the routes recover drugs known to move
  a gene.
- [Improvement review](improvements.md): remaining work and its priority.

The command and API documentation describes the current implementation. Scores
and route labels are research heuristics; they do not establish efficacy.
