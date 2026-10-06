"""Command-line interface.

    python -m promoter_tf_drugs --demo --gene DEMOGENE --out results/
    python -m promoter_tf_drugs --gene SMN2 --upstream 2000 --downstream 500 \
        --pfm-file data/jaspar2024.txt \
        --chip-bed data/remap2022.bed.gz --out results/SMN2
"""

from __future__ import annotations

import argparse
import logging
import sys

from .pipeline import run, write_outputs, markdown_report
from .http import HttpError
from . import progress as _progress


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="promoter-tf-drugs",
        description="Find a gene's promoter, the TFs predicted to bind it, and "
                    "approved drugs that push those TFs in whichever direction "
                    "moves your target gene (--direction up|down).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--gene", required=True, help="gene symbol or Ensembl gene ID")
    p.add_argument("--out", default="results", help="output directory")
    p.add_argument("--top", type=int, default=25, help="rows in the report tables")
    add_analysis_args(p)
    return p


def add_analysis_args(p: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Every flag that configures the analysis itself, shared by both CLIs.

    `panel` needs the same ~40 options as the single-gene entry point, and it
    used to borrow them by reaching into the other parser's private
    `_actions` and re-registering the same Action objects in two parsers.
    That shares mutable state between parsers, which is exactly the kind of
    thing that works until it silently does not. Declared once here instead.

    `--gene` and `--out` stay with their own CLIs, because the panel takes a
    gene LIST and writes a different set of files.
    """
    p.add_argument("--direction", choices=["up", "down"], default="up",
                   help="up: find drugs that should raise target-gene expression "
                        "(agonise activators, block repressors). "
                        "down: the inverse (block activators, agonise repressors)")
    p.add_argument("--species", default="homo_sapiens")
    p.add_argument("--upstream", type=int, default=1000, help="bp 5' of the TSS")
    p.add_argument("--downstream", type=int, default=200, help="bp 3' of the TSS")
    p.add_argument("--all-tss", action="store_true",
                   help="scan every distinct TSS, not just the MANE/canonical one")

    p.add_argument("--pvalue", type=float, default=1e-4,
                   help="per-position false-positive rate for motif hits")
    p.add_argument("--min-ic", type=float, default=6.0,
                   help="drop matrices below this information content (bits)")
    p.add_argument("--pfm-file", help="local JASPAR-format PFM file (recommended)")
    p.add_argument("--matrix-id", action="append", dest="matrix_ids",
                   help="restrict to specific JASPAR matrix IDs (repeatable)")

    p.add_argument("--no-mask-repeats", action="store_true",
                   help="scan repeat-masked (lowercase) sequence too; off by default "
                        "because Alu/simple repeats generate large numbers of "
                        "spurious hits")
    p.add_argument("--no-grade", action="store_true",
                   help="skip drug grading (faster; grading costs ~3 ChEMBL calls "
                        "per drug on the first run)")
    p.add_argument("--predict-direction", action="store_true",
                   help="independently predict agonist vs antagonist from ligand "
                        "structure (needs: pip install rdkit scikit-learn). Trains a "
                        "per-receptor classifier on ChEMBL labels and flags rows where "
                        "it disagrees with the database annotation")
    p.add_argument("--regulon-file",
                   help="TRRUST/DoRothEA TF-target TSV. Scores how many GENES each "
                        "TF controls - a drug hitting one protein that runs 5000 "
                        "genes is not specific. Without this, coarse bundled tiers "
                        "are used. https://www.grnpedia.org/trrust/")
    p.add_argument("--no-direct", action="store_true",
                   help="skip the search for drugs targeting the gene itself")
    p.add_argument("--design-modality", action="store_true",
                   help="DESIGN a direct modality for this gene rather than "
                        "looking up existing ones: siRNA and ASO candidates from "
                        "the transcript, CRISPRi guides, and a deliverability "
                        "verdict for the tissue you named in --tissue")
    p.add_argument("--design-screen", action="store_true",
                   help="design the wet-lab screen this analysis prioritises: "
                        "reporter construct with controls, CRISPRa sgRNAs, and a "
                        "cherry-pick compound list with realistic hit rates")
    p.add_argument("--no-stabilizers", action="store_true",
                   help="skip the stabilisation section (conformational "
                        "stabilisers, molecular glues, degradation blockers) — "
                        "the routes that work on TFs with no ligand pocket")
    p.add_argument("--cell-type-check", action="store_true",
                   help="check that each TF is expressed in the same tissues as "
                        "your gene, and that the drug reaches them. Heavily "
                        "down-weights factors with no tissue overlap")
    p.add_argument("--lincs", action="store_true",
                   help="query LINCS L1000 drug/gene expression signatures; "
                        "distinguishes measured landmarks from inferred values")
    p.add_argument("--gtf",
                   help="GENCODE/Ensembl GTF for the neighbouring-gene check. "
                        "Replaces the Ensembl overlap call, which 500s under "
                        "load. data/gencode.basic.gtf.gz is picked up "
                        "automatically if present (34 MB, from "
                        "ftp.ebi.ac.uk/pub/databases/gencode)")
    p.add_argument("--has-ptc", action="store_true",
                   help="the gene carries a premature termination codon in "
                        "the alleles you care about. Moves nonsense "
                        "readthrough and NMD inhibition to the front of the "
                        "post-transcriptional section, because for a PTC "
                        "allele those address the problem and extra "
                        "transcription does not")
    p.add_argument("--genome",
                   help="soft-masked genome FASTA, for fully offline runs. "
                        "With --gtf this removes the last Ensembl "
                        "dependency: gene models come from the GTF and "
                        "sequence from here. A UCSC hg38.fa is already "
                        "soft-masked, so repeat masking is unchanged. "
                        "data/hg38.fa is found automatically; a .fai index "
                        "is built on first use")
    p.add_argument("--promoter-fasta",
                   help="read promoter sequences from a local FASTA instead "
                        "of Ensembl, matched on exact coordinates. This "
                        "package's own promoter.fa works, so a previous run "
                        "is its own offline source — Ensembl 500s then cost "
                        "you nothing")
    p.add_argument("--expression-file",
                   help="local tissue-expression table (HPA "
                        "rna_tissue_consensus.tsv or a GTEx median-TPM matrix). "
                        "Replaces the Open Targets lookup entirely — faster, "
                        "reproducible, and works with no network")
    p.add_argument("--lincs-gmt", nargs=2, metavar=("UP_GMT", "DOWN_GMT"),
                   help="Enrichr LINCS_L1000_Chem_Pert .gmt files. Replaces "
                        "SigCom, which returns HTTP 500 for days at a time. "
                        "data/lincs_chem_pert_{up,down}.gmt are picked up "
                        "automatically if present (16 MB + 14 MB)")
    p.add_argument("--lincs-geneinfo", metavar="TSV",
                   help="CLUE geneinfo table, which records whether each gene "
                        "was MEASURED on the L1000 array (978 landmarks) or "
                        "imputed from those that were. Without it every LINCS "
                        "hit is tiered as predicted, so measured evidence "
                        "ranks BELOW a curated ligand. "
                        "data/lincs_geneinfo.txt is picked up automatically "
                        "if present (1.1 MB)")
    p.add_argument("--lincs-file",
                   help="local drug/gene signature table (drug,gene,score[,"
                        "cell_line,is_landmark]) exported from CLUE or iLINCS. "
                        "Offline replacement for --lincs")
    p.add_argument("--signature-file",
                   help="reference gene,score signature; records the gene count "
                        "when --lincs or --lincs-file is active. Does not run "
                        "a connectivity search")
    p.add_argument("--genetics", action="store_true",
                   help="human genetic evidence via Open Targets: what "
                        "phenotypes associate with variation in this gene "
                        "(drug-target MR in reverse)")
    p.add_argument("--proliferation", default="unknown",
                   choices=["cycling", "post-mitotic", "mixed", "unknown"],
                   help="proliferation state of your target tissue. post-mitotic "
                        "(neurons, cardiomyocytes) REMOVES cell-cycle-restricted "
                        "TFs such as E2F/FOXM1, whose targets are actively "
                        "repressed there; cycling down-weights them by duty cycle")
    p.add_argument("--upstream-network", action="store_true",
                   help="walk upstream through signed causal edges (OmniPath) to "
                        "find druggable kinases/receptors controlling TFs that have "
                        "no direct pharmacology. Signs propagate, so inhibiting an "
                        "inhibitor counts as activation")
    p.add_argument("--upstream-depth", type=int, default=3,
                   help="max hops upstream. Each hop costs selectivity; 4 is "
                        "usually where paths stop being interpretable")
    p.add_argument("--include-unknown-direction", action="store_true",
                   help="include curated drug-gene links that record no "
                        "direction. Off by default: they cannot answer the "
                        "question this pipeline asks, and they fill the tables "
                        "with UNKNOWN/unclear rows")
    p.add_argument("--no-mechanisms", action="store_true",
                   help="skip the orthosteric vs allosteric (PAM/NAM) comparison")
    p.add_argument("--no-improvements", action="store_true",
                   help="skip the comparison against the current standard of care")
    p.add_argument("--no-outcome", action="store_true",
                   help="skip the transcriptional-outcome probability chain")
    p.add_argument("--screen", action="store_true",
                   help="virtual-screen ALL approved drugs against the top "
                        "receptors to find untested repurposing pairs. Slow "
                        "(downloads ~4000 structures once) and needs rdkit")
    p.add_argument("--tissue", choices=["cns"],
                   help="down-weight drugs unlikely to reach this compartment")
    p.add_argument("--tf-activity", action="store_true",
                   help="infer which factors are actually ACTIVE from their "
                        "targets' expression (VIPER-style), instead of only "
                        "guessing from motifs. Uses --regulon-file, or "
                        "data/collectri.tsv or data/trrust.tsv if present, "
                        "plus --expression-file and --activity-tissue for "
                        "the signature. CollecTRI is the better table: it "
                        "signs every edge and covers about three times as "
                        "many factors as TRRUST. This is the "
                        "independent check on the motif scan: a site with a "
                        "flat regulon is probably scan noise")
    p.add_argument("--activity-tissue", metavar="TISSUE",
                   help="tissue name as spelled in the expression table "
                        "(e.g. liver, brain, kidney) to infer TF activity "
                        "in. Required by --tf-activity")
    p.add_argument("--docking-results", metavar="JSON",
                   help="fold in a docking run performed elsewhere. "
                        "Written by `python -m promoter_tf_drugs.docking "
                        "--bundle DIR --out JSON`, which is the step "
                        "that needs Linux. Prepare and report here, "
                        "dock there")
    p.add_argument("--promoter-dna", action="store_true",
                   help="look for sequence features in the promoter that "
                        "DNA-binding small molecules prefer - quadruplexes, "
                        "GC-boxes, AT tracts - and propose agents only "
                        "where a feature beats shuffles of the same "
                        "sequence. Nearly every approved agent of this "
                        "kind is cytotoxic")
    p.add_argument("--screen-product", action="store_true",
                   help="virtual-screen the approved library against the "
                        "GENE'S OWN protein, not just the transcription "
                        "factors. Finds compounds whose chemistry fits the "
                        "product's known ligands and which nobody has "
                        "annotated against it. This changes what the "
                        "existing protein does rather than how much of it "
                        "there is, so RT-qPCR will not see it")
    p.add_argument("--sensitivity-scan", action="store_true",
                   help="correlate every approved compound in the PRISM "
                        "panel against each factor's CRISPR knockout "
                        "profile. Finds drugs that behave as though they "
                        "remove the factor, with no annotation of any kind "
                        "behind them - the one route here that can surface "
                        "a compound nobody has linked to the gene. Needs "
                        "both the DepMap and PRISM matrices on disk")
    p.add_argument("--prism-file", metavar="CSV",
                   help="PRISM Repurposing LFC matrix. "
                        "data/Repurposing_Public_24Q2_LFC_COLLAPSED.csv is "
                        "picked up automatically if present, and the "
                        "compound list beside it turns BRD ids into drug "
                        "names")
    p.add_argument("--depmap-file", metavar="CSV",
                   help="DepMap CRISPRGeneEffect.csv. Measures what happens "
                        "to a cell when each factor is knocked out, which is "
                        "the collateral question the regulon count stands in "
                        "for. data/CRISPRGeneEffect.csv is picked up "
                        "automatically if present; Model.csv beside it adds "
                        "lineage")
    p.add_argument("--activity-permutations", type=int, default=0,
                   metavar="N",
                   help="build an empirical null from N resampled target "
                        "sets instead of the analytic one. The analytic "
                        "p-value assumes targets are independent, which "
                        "co-expression makes false, so use this (2000+) "
                        "before relying on a number")
    p.add_argument("--tf-weight", type=float, default=0.5,
                   help="0..1. Share of the combined grade carried by TF confidence "
                        "rather than drug quality. 0.5 weights them equally; raise it "
                        "when you trust the promoter evidence (e.g. with --chip-bed)")
    p.add_argument("--chip-bed", help="ChIP-seq peak BED (ReMap/ChIP-Atlas) for "
                                      "orthogonal binding evidence")
    p.add_argument("--chip-api", action="store_true",
                   help="get ChIP evidence from the ReMap REST API instead of a "
                        "local BED. No multi-GB download and far faster per gene; "
                        "needs network access. --chip-bed wins if both are given")

    p.add_argument("--include-investigational", action="store_true",
                   help="do not restrict to approved (max_phase 4) drugs")
    p.add_argument("--no-dgidb", action="store_true", help="ChEMBL only")
    p.add_argument("--find-analogs", action="store_true",
                   help="search ChEMBL for approved drugs chemically similar "
                        "to this report's own best candidates, subtracting "
                        "everything already annotated at the factor. Yields "
                        "one or two untested pairs per factor; costs a few "
                        "requests")
    p.add_argument("--no-opentargets", action="store_true",
                   help="skip the Open Targets drug lookup. It is a third "
                        "independent curation and typically adds one or two "
                        "approved drugs per target that ChEMBL and DGIdb "
                        "miss, at the cost of one request per factor")

    p.add_argument("--demo", action="store_true",
                   help="run entirely offline on bundled fixtures")
    p.add_argument("--offline", action="store_true",
                   help="fail rather than hit the network; use cached responses only")
    p.add_argument("--cache-dir")
    p.add_argument("--quiet", action="store_true",
                   help="no progress checkpoints. They go to stderr, so "
                        "they do not touch a piped report; this is for "
                        "when stderr itself is being read")
    p.add_argument("-v", "--verbose", action="count", default=0)
    return p


def _check_flag_combinations(parser, args) -> None:
    """Reject flag combinations that cannot work, before anything runs.

    A full run is a quarter of an hour of scanning and lookups. Learning
    at the end that a requested source could never have produced
    anything wastes all of it, and the condition is knowable before the
    first request goes out.

    This is the CLI's check, not the pipeline's: `run()` still records a
    failed source and carries on, because a library caller may want the
    rest of the run regardless.
    """
    if getattr(args, "tf_activity", False) and not getattr(
            args, "activity_tissue", None):
        parser.error(
            "--tf-activity needs --activity-tissue: it infers which "
            "factors are active from their targets' expression in one "
            "tissue, and without a tissue there is no expression profile "
            "to infer from. Name a tissue as the expression table spells "
            "it, e.g. --activity-tissue liver.")
    if getattr(args, "activity_permutations", 0) and not getattr(
            args, "tf_activity", False):
        parser.error(
            "--activity-permutations only applies to --tf-activity, "
            "which was not requested. It would do nothing here.")


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _check_flag_combinations(parser, args)
    # Checkpoints to stderr. The report goes to stdout, so a
    # piped run is unaffected.
    if not getattr(args, "quiet", False):
        _progress.enable()
    if args.top < 1:
        parser.error("--top must be a positive integer")

    # Windows consoles default to cp1252, which cannot encode the characters
    # used in the report tables. Print UTF-8 where possible; degrade to "?"
    # rather than crashing on terminals that genuinely cannot display it.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    logging.basicConfig(
        level=[logging.WARNING, logging.INFO, logging.DEBUG][min(args.verbose, 2)],
        format="%(levelname)s %(name)s: %(message)s",
    )

    try:
        return _execute(args)
    except (ValueError, OSError, HttpError) as exc:
        if args.verbose:
            logging.getLogger(__name__).exception("analysis failed")
        else:
            print(f"{parser.prog}: error: {exc}", file=sys.stderr)
        return 1


# argparse dest -> run() parameter, where the two names differ.
_ARG_TO_RUN = {
    "direction": "goal",
    "lincs": "use_lincs",
    "genetics": "use_genetics",
    "screen": "screen_receptors",
    "include_unknown_direction": "allow_unknown_direction",
}

# Flags phrased as the negation of the parameter they set.
_NEGATED_ARGS = {
    "no_mask_repeats": "mask_repeats",
    "no_grade": "grade_drugs",
    "no_direct": "include_direct",
    "no_outcome": "estimate_outcome",
    "no_improvements": "find_improvements",
    "no_mechanisms": "compare_mechanisms",
    "no_stabilizers": "include_stabilizers",
    "include_investigational": "approved_only",
    "no_dgidb": "use_dgidb",
    "no_opentargets": "use_opentargets",
}


def run_kwargs_from_args(args) -> dict:
    """Every `pipeline.run` argument a parsed namespace supplies.

    Derived from run()'s own signature rather than listed by hand, and
    shared by the single-gene CLI and the panel so the two cannot
    diverge. The panel used to keep a tuple of twenty-one keys against
    the fifty run() takes, so --genetics, --lincs, --upstream-network,
    --tf-activity and --docking-results were accepted and dropped
    without a word.

    A flag absent from one parser is simply absent here; run()'s own
    default then applies.
    """
    import inspect

    from .pipeline import run as _run

    params = {n for n in inspect.signature(_run).parameters if n != "gene"}
    out: dict = {}
    for dest, value in vars(args).items():
        if dest in _NEGATED_ARGS:
            out[_NEGATED_ARGS[dest]] = not value
        elif dest in _ARG_TO_RUN:
            out[_ARG_TO_RUN[dest]] = value
        elif dest in params:
            out[dest] = value
    return {k: v for k, v in out.items() if k in params}


def _execute(args) -> int:

    result = run(args.gene, **run_kwargs_from_args(args))

    written = write_outputs(result, args.out, top_n=args.top)
    print(markdown_report(result, top_n=args.top))
    print("\nWrote:", ", ".join(str(p) for p in written), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
