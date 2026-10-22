"""How much message each isoform actually carries.

Why this exists
---------------
`poisonexon` lists exons that only decay-prone isoforms use, as places an
antisense oligonucleotide could act. That list is drawn from the
annotation, and the annotation says an isoform has been *seen*, not that
it carries any appreciable share of the gene's message. Both ends matter
and they are very different results:

    the decay isoforms carry a third of the message  -> a real lever
    the decay isoforms carry none of it              -> nothing to win

Without a number the report shows a table of candidate exons and leaves
the reader to assume the first case. This supplies the number.

Where the number comes from
---------------------------
GTEx v8, median transcript TPM, through the portal's public API. GTEx is
the only resource that measures every isoform of a gene across hundreds
of donors in a named tissue, which is exactly the quantity wanted. The
default tissue is cultured fibroblasts because that is the cell the
screen in this package is designed around; any GTEx tissue id works.

recount3 would answer the same question over far more samples, but it is
distributed as R objects of a size no stdlib client should be fetching,
so it belongs in a table you supply rather than in a live lookup. A local
table is read in preference to the API for that reason.

What the number can and cannot say
----------------------------------
**It is a different annotation.** GTEx v8 is built on GENCODE v26.
Current Ensembl has isoforms that did not exist then, and those have no
measurement at all - not a measurement of zero. The two are reported
separately, because collapsing them would turn "never measured" into
"measured as absent", which is the error this whole module exists to
prevent.

**It is a median across donors, in one tissue.** A transcript at zero
median can still be present in some donors, and splicing in a patient's
cells is not obliged to match a GTEx tissue.

**TPM is not protein.** A decay isoform is destroyed quickly, so a low
steady-state share is partly a consequence of the decay it is named for.
The share bounds what redirecting it could add to the pool, which is the
question here, but it is not the flux through that isoform.
"""

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)

API = "https://gtexportal.org/api/v2"
DATASET = "gtex_v8"
GENCODE_VERSION = "v26"
GENOME = "GRCh38/hg38"

# Cultured fibroblasts: the cell the protein screen in this package is
# written around, and one of the few GTEx tissues that is a cell line you
# can actually obtain from a patient.
DEFAULT_TISSUE = "Cells_Cultured_fibroblasts"

DEFAULT_PATHS = ("data/transcript_expression.tsv",)

# Below this share of the gene's message, redirecting an isoform cannot
# add enough to matter even if every molecule were recovered. 5% of a
# transcript that is already at half normal is a 2.5% gain.
NEGLIGIBLE_SHARE = 0.05

# Biotypes whose message is non-productive, and so in principle
# recoverable. `retained_intron` is included deliberately: it is not a
# poison exon and `poisonexon` does not list it, but it is the other way
# a gene parks message in a form that makes no protein, and when it holds
# the share the report should say so rather than being silent.
NON_PRODUCTIVE = {
    "nonsense_mediated_decay": "decay",
    "retained_intron": "retained intron",
    "protein_coding_CDS_not_defined": "no defined coding sequence",
    "processed_transcript": "no defined coding sequence",
}

SOURCE = ("GTEx v8 median transcript TPM (GENCODE v26) "
          "https://gtexportal.org/api/v2 - or a local table at "
          "data/transcript_expression.tsv with columns "
          "transcript, tpm [, tissue].")


def _resolve(path=None) -> Path | None:
    cands = [path] if path else list(DEFAULT_PATHS)
    return next((Path(c) for c in cands if c and Path(c).exists()), None)


def _bare(tid: str) -> str:
    """Versions differ between sources; ENST123.10 and ENST123 are one transcript."""
    return (tid or "").strip().upper().split(".")[0]


def from_table(gene: str, *, tissue: str | None = None, path=None) -> dict:
    """Per-transcript TPM from a local table, if one is present.

    Columns, header required, order free, tab or comma:

        transcript, tpm [, tissue] [, gene] [, source]

    A `tissue` column is filtered on when `tissue` is given; a table with
    no such column is taken as already being the tissue you want.
    """
    import csv

    found = _resolve(path)
    if found is None:
        return {}
    want_tissue = (tissue or "").strip().lower()
    tpm: dict = {}
    try:
        with open(found, encoding="utf-8", errors="replace") as fh:
            sample = fh.read(8192)
            fh.seek(0)
            delim = "\t" if sample.count("\t") >= sample.count(",") else ","
            for row in csv.DictReader(fh, delimiter=delim):
                row = {(k or "").strip().lower(): v for k, v in row.items()}
                if row.get("gene") and gene and \
                        row["gene"].strip().upper() != gene.strip().upper():
                    continue
                if want_tissue and row.get("tissue") and \
                        row["tissue"].strip().lower() != want_tissue:
                    continue
                tid = _bare(row.get("transcript") or row.get("transcript_id"))
                try:
                    val = float(row.get("tpm") or row.get("median"))
                except (TypeError, ValueError):
                    continue
                if tid:
                    tpm[tid] = val
    except OSError as exc:
        log.warning("could not read the transcript table %s (%s)", found, exc)
        return {}
    if tpm:
        log.info("transcript expression for %s: %d isoform(s) from %s",
                 gene, len(tpm), found.name)
    return {"tpm": tpm, "tissue": tissue or "as supplied",
            "source": str(found), "dataset": "local table"}


def fetch(client, gene: str, *, tissue: str | None = None) -> dict:
    """Per-transcript median TPM for `gene` in one GTEx tissue.

    Prefers a local table. Returns `{}` when neither is available, and a
    dict carrying `error` when the lookup was attempted and failed, so a
    failed request is never reported as an isoform carrying no message.
    """
    tissue = tissue or DEFAULT_TISSUE
    local = from_table(gene, tissue=tissue)
    if local:
        return local
    if client is None:
        return {}
    try:
        ref = client.get(f"{API}/reference/gene",
                         {"geneId": gene, "gencodeVersion": GENCODE_VERSION,
                          "genomeBuild": GENOME})
        rows = (ref or {}).get("data") or []
        hit = next((r for r in rows
                    if (r.get("geneSymbolUpper") or "").upper()
                    == gene.strip().upper()), None)
        if hit is None:
            return {"tpm": {}, "tissue": tissue, "dataset": DATASET,
                    "error": f"GTEx (GENCODE {GENCODE_VERSION}) has no gene "
                             f"called {gene}"}
        expr = client.get(f"{API}/expression/medianTranscriptExpression",
                          {"gencodeId": hit["gencodeId"],
                           "datasetId": DATASET,
                           "tissueSiteDetailId": tissue,
                           "itemsPerPage": 200})
    except Exception as exc:                                # noqa: BLE001
        log.warning("GTEx transcript expression for %s failed: %s",
                    gene, str(exc)[:100])
        return {"tpm": {}, "tissue": tissue, "dataset": DATASET,
                "error": f"the GTEx lookup failed ({str(exc)[:80]})"}

    tpm = {}
    for r in (expr or {}).get("data") or []:
        try:
            tpm[_bare(r.get("transcriptId"))] = float(r.get("median"))
        except (TypeError, ValueError):
            continue
    tpm.pop("", None)
    if not tpm:
        return {"tpm": {}, "tissue": tissue, "dataset": DATASET,
                "error": f"GTEx returned no transcript rows for {gene} in "
                         f"{tissue}; the tissue id may be wrong"}
    log.info("transcript expression for %s: %d isoform(s) from GTEx %s",
             gene, len(tpm), tissue)
    return {"tpm": tpm, "tissue": tissue, "dataset": DATASET,
            "gencode_id": hit["gencodeId"], "source": API}


def weigh(gene_record: dict, expression: dict) -> dict:
    """Split the gene's message between productive and non-productive forms.

    `gene_record` is the Ensembl record `poisonexon.find` reads, so the
    biotypes come from current Ensembl while the TPM comes from GTEx.
    Transcripts Ensembl lists and GTEx never measured are counted
    separately and never treated as zero.
    """
    transcripts = (gene_record or {}).get("Transcript") or []
    tpm = (expression or {}).get("tpm") or {}
    out = {
        "tissue": (expression or {}).get("tissue"),
        "dataset": (expression or {}).get("dataset"),
        "error": (expression or {}).get("error"),
        "total_tpm": 0.0, "coding_tpm": 0.0,
        "by_class": {}, "per_transcript": {},
        "n_annotated": len(transcripts), "n_measured": 0,
        "unmeasured_by_class": {}, "verdict": "",
    }
    if not transcripts or not tpm:
        return out

    for t in transcripts:
        tid = _bare(t.get("id"))
        bio = t.get("biotype") or "?"
        label = NON_PRODUCTIVE.get(bio, "productive"
                                   if bio == "protein_coding" else bio)
        if tid not in tpm:
            out["unmeasured_by_class"][label] = \
                out["unmeasured_by_class"].get(label, 0) + 1
            continue
        val = tpm[tid]
        out["n_measured"] += 1
        out["total_tpm"] += val
        out["per_transcript"][tid] = {"tpm": round(val, 3), "biotype": bio,
                                      "class": label}
        out["by_class"][label] = round(
            out["by_class"].get(label, 0.0) + val, 3)
        if bio == "protein_coding":
            out["coding_tpm"] += val

    total = out["total_tpm"]
    out["total_tpm"] = round(total, 3)
    out["coding_tpm"] = round(out["coding_tpm"], 3)
    for tid, rec in out["per_transcript"].items():
        rec["share"] = round(rec["tpm"] / total, 4) if total else 0.0
    out["shares"] = {k: (round(v / total, 4) if total else 0.0)
                     for k, v in out["by_class"].items()}
    decay = out["shares"].get("decay", 0.0)
    recoverable = sum(v for k, v in out["shares"].items()
                      if k != "productive")
    out["decay_share"] = decay
    out["non_productive_share"] = round(recoverable, 4)
    out["verdict"] = _verdict(out, decay, recoverable)
    return out


def _verdict(out: dict, decay: float, recoverable: float) -> str:
    unmeasured_decay = out["unmeasured_by_class"].get("decay", 0)
    tissue = out["tissue"]
    tail = ""
    if unmeasured_decay:
        tail = (f" {unmeasured_decay} decay isoform(s) Ensembl lists now are "
                f"newer than the GENCODE build GTEx was quantified on, so "
                f"they have no measurement at all - which is not a "
                f"measurement of zero.")
    if decay >= NEGLIGIBLE_SHARE:
        head = (f"Decay isoforms carry {decay:.1%} of this gene's message in "
                f"{tissue}, so redirecting them is worth something: that is "
                f"the ceiling on what skipping a brake exon could add, "
                f"before any question of whether an oligonucleotide can be "
                f"made to do it.")
    elif decay > 0:
        head = (f"Decay isoforms carry only {decay:.1%} of this gene's "
                f"message in {tissue}. Recovering all of it would add less "
                f"than a twentieth to the pool, so splice correction is a "
                f"poor lever here on this evidence.")
    else:
        head = (f"No decay isoform of this gene carries measurable message "
                f"in {tissue}. The candidate brake exons below sit on "
                f"transcripts that are not detectably made in this tissue, "
                f"so there is nothing for an oligonucleotide to redirect.")
    other = recoverable - decay
    if other >= NEGLIGIBLE_SHARE:
        worst = max(((k, v) for k, v in out["shares"].items()
                     if k not in ("productive", "decay")),
                    key=lambda kv: kv[1], default=(None, 0))
        head += (f" The non-productive message that does exist is mostly "
                 f"{worst[0]} ({worst[1]:.1%}), which a brake-exon "
                 f"oligonucleotide does not address - a different "
                 f"intervention would be needed to recover it.")
    return head + tail


def weigh_candidates(candidates: list, weighed: dict) -> list:
    """Attach to each candidate exon the share of message it sits on.

    An exon used by three decay isoforms that are each at zero is not a
    better target than one used by a single isoform at zero. What ranks
    them is how much message they carry between them, so each candidate
    gets `share` and, separately, how many of its isoforms were never
    measured - the two must not be added together.
    """
    per = (weighed or {}).get("per_transcript") or {}
    for c in candidates or []:
        ids = [_bare(t) for t in (c.get("transcripts") or [])]
        measured = [per[i]["share"] for i in ids if i in per]
        c["share"] = round(sum(measured), 4) if measured else 0.0
        c["n_measured"] = len(measured)
        c["n_unmeasured"] = len(ids) - len(measured)
    return candidates


def summarise(weighed: dict, gene: str) -> str:
    """One paragraph for the report."""
    if not weighed or weighed.get("error"):
        why = (weighed or {}).get("error") or "no isoform-level expression " \
                                              "was available"
        return (f"**The isoforms of {gene} were not weighed** ({why}), so "
                f"any brake exon below is a place to look and not a "
                f"quantified opportunity. {SOURCE}")
    if not weighed.get("n_measured"):
        return (f"**No isoform of {gene} was measured** in "
                f"{weighed.get('tissue')}, so the brake exons below carry no "
                f"weight either way.")
    return (f"{weighed['verdict']} Measured on {weighed['n_measured']} of the "
            f"{weighed['n_annotated']} isoforms Ensembl lists, "
            f"{weighed['total_tpm']:.1f} TPM in total.")
