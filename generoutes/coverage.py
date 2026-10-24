"""Which of the six steps to a protein this run actually looked at.

The gap this exists to show
---------------------------
How much functional protein a gene yields is set by six steps in series:
transcription, splicing and isoform choice, mRNA stability, translation,
folding, and degradation. A report can fill two hundred rows about the
first one and say nothing about the other five, and nothing in it would
look wrong - the tables would be full, the evidence ledger would be
green, and a reader would have no way to see that five sixths of the
problem was never examined.

That is not hypothetical. A promoter-first tool searches transcription
hardest, and for a disease caused by a coding or splice variant the
limiting step is almost certainly elsewhere.

Searched, listed, and absent
----------------------------
The distinction this module draws is the one that matters:

- **searched** - a gene-specific lookup ran and returned a count. The
  answer may be zero, and zero is then a result.
- **listed** - the report names the mechanism and the drug classes that
  work that way, but nothing was looked up for THIS gene. A menu, not a
  search. Most of the post-transcriptional section is this.
- **absent** - neither.

A step that is only listed is not covered, and saying so is the whole
point. The ledger does this for data sources; this does it for biology.
"""

from __future__ import annotations

SEARCHED = "searched"
LISTED = "listed only"
ABSENT = "not covered"

# The categories in `catalogue.CATEGORIES` that act on transcription.
_TRANSCRIPTIONAL = (
    "TF: ", "promoter DNA binder", "upstream signalling", "chromatin-level",
    "literature: reported effect", "empirical signature",
    "measured phenocopy", "off-pathway pocket binder",
    "virtual screen hypothesis",
)


def _routes(block: dict) -> list:
    return ((block.get("all_routes") or {}).get("routes") or [])


def _n_transcriptional(block: dict) -> int:
    return sum(1 for r in _routes(block)
               if any(str(r.get("category", "")).startswith(p)
                      for p in _TRANSCRIPTIONAL))


def audit(result: dict) -> list[dict]:
    """One row per step: what ran, what it found, and what would fix a gap."""
    block = next(iter(result.get("results") or []), {}) or {}
    beyond = block.get("beyond_transcription") or {}
    brake = result.get("brake_exons") or {}
    weight = brake.get("isoform_weight") or {}
    weighed = bool(weight.get("n_measured")) and not weight.get("error")
    ptx = (beyond.get("post_transcriptional") or {}).get("routes") or []
    stability = beyond.get("product_stability") or {}
    product = block.get("product_screen") or {}
    protein = block.get("protein_abundance") or []
    degron = stability.get("degradation") or {}
    u = block.get("utr") or {}
    u_tr = u.get("translation") or {}
    u_st = u.get("stability") or {}
    # A UTR scan that RAN is a search, whatever it found: "no uORF in a
    # 200 bp 5'UTR" is a result about this transcript.
    tr_searched = bool(u_tr) and not u_tr.get("error")
    st_searched = bool(u_st) and not u_st.get("error")

    rows = [
        {
            "step": "transcription",
            "state": SEARCHED,
            "found": _n_transcriptional(block),
            "what_ran": "the promoter scan, the factor ligands, the known "
                        "regulators, the upstream network, the literature "
                        "record and the expression signatures",
            "gap": "",
            "enable": "",
        },
        {
            "step": "splicing and isoform choice",
            "state": SEARCHED if brake.get("candidates") else (
                LISTED if ptx else ABSENT),
            "found": len(brake.get("candidates") or []),
            "what_ran": "the decay-isoform scan: exons used only by "
                        "transcripts the cell destroys, which an "
                        "oligonucleotide could redirect"
                        + (", each weighed by the share of the gene's "
                           "message its isoforms carry in GTEx"
                           if weighed else ""),
            "gap": "no drug was searched against them. Splice modulation "
                   "has an approved small molecule (risdiplam), but "
                   "nothing here looks for one against this gene's exons, "
                   "and no public corpus of compound-by-splicing effects "
                   "is read."
                   + ("" if weighed else " How much message the decay "
                      "isoforms carry is also unknown, so the exons are "
                      "unranked."),
            "enable": "a splice-modulator screen"
                      + ("" if weighed else ", or isoform-resolved "
                         "expression to show how much message the decay "
                         "isoforms actually carry"),
        },
        {
            "step": "mRNA stability",
            "state": SEARCHED if st_searched else (LISTED if ptx else ABSENT),
            "found": sum(1 for k in ("pentamer", "nonamer")
                         if (u_st.get(k) or {}).get("enriched")),
            "what_ran": "the 3'UTR scanned for AU-rich decay elements, "
                        "each against a shuffle of the same sequence",
            "gap": "the element is annotated and no drug was matched to "
                   "it. No corpus of transcript half-life under drug is "
                   "read, and the microRNA that represses it is not "
                   "identified.",
            "enable": "a half-life dataset (SLAM-seq or actinomycin-D "
                      "chase under compound), or a curated miRNA-target "
                      "table for this transcript",
        },
        {
            "step": "translation",
            "state": SEARCHED if tr_searched else (LISTED if ptx else ABSENT),
            "found": u_tr.get("n_observed", 0) if tr_searched else 0,
            "what_ran": "the 5'UTR scanned for upstream ORFs against a "
                        "shuffle of the same sequence, with the Kozak "
                        "strength of each",
            "gap": "a uORF is a lever and no drug was matched to it. "
                   "Nothing reads ribosome profiling, so whether protein "
                   "is low while message is normal - the signature that "
                   "would point here - still needs the two assays in "
                   "patient cells.",
            "enable": "ribosome profiling under compound, or simply the "
                      "transcript-versus-protein comparison in patient "
                      "cells, which costs two assays",
        },
        {
            "step": "folding",
            "state": SEARCHED if product else (LISTED if stability else ABSENT),
            "found": len(product.get("hits") or []),
            "what_ran": "the ligand screen against the gene's own protein, "
                        "which is what a pharmacological chaperone would "
                        "have to be",
            "gap": product.get("error") or "",
            "enable": "measured ligands for this protein; the screen "
                      "refuses below the number needed to train a model, "
                      "rather than guessing",
        },
        {
            "step": "degradation",
            "state": SEARCHED if (protein or degron) else (
                LISTED if stability else ABSENT),
            "found": len(protein),
            "what_ran": "the named degradation route for this protein, and "
                        "the drug-to-protein-abundance table",
            "gap": "",
            "enable": "",
        },
    ]

    if degron:
        rows[-1]["what_ran"] += f" ({degron.get('ligase')})"
    if not protein:
        # An absent table and a table that does not quantify this protein
        # both give no hits and need opposite things done about them.
        from . import proteomics as _pr

        gene = ((result.get("query") or {}).get("gene") or "this gene")
        if not _pr.available():
            rows[-1]["gap"] = ("no drug-to-protein-abundance table on disk, "
                               "so no compound was checked against the "
                               "protein itself - only against transcripts")
            rows[-1]["enable"] = ("decryptE or ProTargetMiner at "
                                  "data/drug_protein_abundance.tsv")
        elif not (_pr.describe(gene) or {}).get("n_compounds"):
            rows[-1]["gap"] = (f"a drug-to-protein-abundance table is on "
                               f"disk and does not quantify {gene} at all, "
                               f"so nothing in it bears on this protein. "
                               f"That is an absent measurement, not a "
                               f"negative one.")
            rows[-1]["enable"] = ("proteomics that covers this protein - "
                                  "a targeted assay, or a dataset whose "
                                  "depth reaches it")
        else:
            rows[-1]["gap"] = ("the table quantifies this protein and no "
                               "compound in it moved the protein the "
                               "wanted way. That is a measured negative "
                               "over the compounds it covers, and says "
                               "nothing about the rest of the library.")
            rows[-1]["enable"] = ("proteomics over more compounds than the "
                                  "current table holds")
    return rows


def summarise(rows: list) -> str:
    """One paragraph: how much of the problem this run actually examined."""
    if not rows:
        return ""
    searched = [r for r in rows if r["state"] == SEARCHED]
    listed = [r for r in rows if r["state"] == LISTED]
    n_t = next((r["found"] for r in rows if r["step"] == "transcription"), 0)
    others = sum(r["found"] for r in rows if r["step"] != "transcription")
    return (
        f"**{len(searched)} of the {len(rows)} steps between a gene and a "
        f"working protein were searched for this gene; {len(listed)} are "
        f"named without being searched.** The ranking holds {n_t} "
        f"transcriptional route(s) against {others} from every other step "
        f"combined. That ratio is a property of this package, not of the "
        f"biology: it searches transcription hardest because that is what "
        f"it was built to do. If the pathogenic variant is in the coding "
        f"sequence or a splice site, the limiting step is one of the "
        f"others, and a full table above is not evidence that it was "
        f"considered.")
