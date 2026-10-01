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

    p.add_argument("--demo", action="store_true",
                   help="run entirely offline on bundled fixtures")
    p.add_argument("--offline", action="store_true",
                   help="fail rather than hit the network; use cached responses only")
    p.add_argument("--cache-dir")
    p.add_argument("-v", "--verbose", action="count", default=0)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
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


def _execute(args) -> int:

    result = run(
        args.gene,
        upstream=args.upstream, downstream=args.downstream,
        pvalue=args.pvalue, species=args.species,
        pfm_file=args.pfm_file, chip_bed=args.chip_bed, chip_api=args.chip_api,
        matrix_ids=args.matrix_ids, all_tss=args.all_tss,
        mask_repeats=not args.no_mask_repeats, goal=args.direction,
        grade_drugs=not args.no_grade, tissue=args.tissue,
        tf_weight=args.tf_weight, predict_direction=args.predict_direction,
        regulon_file=args.regulon_file, include_direct=not args.no_direct,
        estimate_outcome=not args.no_outcome, screen_receptors=args.screen,
        find_improvements=not args.no_improvements,
        compare_mechanisms=not args.no_mechanisms,
        proliferation=args.proliferation,
        allow_unknown_direction=args.include_unknown_direction, cell_type_check=args.cell_type_check,
        include_stabilizers=not args.no_stabilizers,
        design_screen=args.design_screen,
        design_modality=args.design_modality,
        use_lincs=args.lincs,
        signature_file=args.signature_file, lincs_file=args.lincs_file,
        expression_file=args.expression_file,
        promoter_fasta=args.promoter_fasta, use_genetics=args.genetics,
        upstream_network=args.upstream_network,
        upstream_depth=args.upstream_depth,
        approved_only=not args.include_investigational,
        use_dgidb=not args.no_dgidb,
        demo=args.demo, offline=args.offline,
        cache_dir=args.cache_dir, min_ic=args.min_ic,
    )

    written = write_outputs(result, args.out, top_n=args.top)
    print(markdown_report(result, top_n=args.top))
    print("\nWrote:", ", ".join(str(p) for p in written), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
