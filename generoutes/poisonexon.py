"""Exons that only decay-prone isoforms use: candidate brakes on output.

Why this exists
---------------
A gene with one working copy makes about half the normal amount of
transcript. One way to raise the output of that copy without touching
transcription is to remove a brake in the RNA. Many genes make some of their
transcripts with an extra ("poison") exon that puts a premature stop in the
message, so the transcript is destroyed by nonsense-mediated decay. An
antisense oligonucleotide that makes the cell skip that exon turns those
transcripts into productive ones and raises the protein, from the working copy
only, in the cells that express the gene. This is the idea behind targeted
augmentation of nuclear gene output.

What this module does
---------------------
It reads the gene's isoforms from Ensembl, which the run has already
fetched, and lists the exons that nonsense-mediated-decay isoforms contain and
coding isoforms do not. Exons are grouped by their coordinates, so an exon
used by three decay isoforms appears once with a count of three.

What it does not do
-------------------
**An annotated exon is not an included exon.** Ensembl lists isoforms that
have been seen at some level in some tissue. Whether the exon is included in
the cells you care about, how often, and whether skipping it raises protein is
a question for RNA sequencing and a splicing assay. This lists where to look.

**It is not an oligonucleotide design.** Nothing here picks a target sequence,
checks it against the rest of the transcriptome, or says the exon boundary is
reachable.

**A decay isoform without a candidate exon** is decayed for another reason (an
extended or shortened coding exon, or a splice pattern that leaves a stop
early). Those are counted and left out of the list.

**A local GTF cannot answer this.** Its records carry no exon coordinates, and
the "basic" GENCODE set leaves out most decay isoforms, so the run reads the
isoform list from Ensembl even when the gene models came from a GTF. A gene
record with no decay isoform is not the same as the gene having none.
"""

from __future__ import annotations

NMD = "nonsense_mediated_decay"
CODING = "protein_coding"


def _span(exon: dict) -> tuple:
    return int(exon["start"]), int(exon["end"])


def find(gene: dict) -> dict:
    """Candidate brake exons from an Ensembl gene record.

    `gene` is the record `genes.get_promoters` reads: `Transcript` entries,
    each with `biotype` and an `Exon` list. Returns a dict with the counts,
    the `candidates` (most-supported first) and a one-line `note`. A record
    without exon coordinates gives `error`, not an empty list.
    """
    symbol = gene.get("display_name") or gene.get("id") or "?"
    transcripts = gene.get("Transcript") or []
    out = {"gene": symbol, "n_transcripts": len(transcripts),
           "n_coding": 0, "n_nmd": 0, "n_nmd_without_candidate": 0,
           "candidates": [], "error": None, "note": ""}
    if not transcripts:
        out["error"] = "no transcripts in the gene record"
        return out
    if not any(t.get("Exon") for t in transcripts):
        out["error"] = ("the gene record carries no exon coordinates, so "
                        "isoforms cannot be compared")
        return out

    coding = [t for t in transcripts if t.get("biotype") == CODING]
    nmd = [t for t in transcripts if t.get("biotype") == NMD]
    out["n_coding"], out["n_nmd"] = len(coding), len(nmd)
    if not nmd:
        out["note"] = (f"Ensembl lists no nonsense-mediated-decay isoform of "
                       f"{symbol}, so there is no annotated brake exon to "
                       f"target. That depends on the annotation: a local "
                       f"'basic' GTF omits most of them.")
        return out

    coding_spans = {_span(e) for t in coding for e in t.get("Exon") or []}

    def overlaps_coding(span) -> bool:
        s, e = span
        return any(s <= ce and e >= cs for cs, ce in coding_spans)

    found: dict = {}
    for t in nmd:
        had = False
        for ex in t.get("Exon") or []:
            span = _span(ex)
            if span in coding_spans:
                continue                       # an exon the coding forms use
            kind = ("novel exon" if not overlaps_coding(span)
                    else "extended or shifted coding exon")
            rec = found.setdefault((span, kind), {
                "start": span[0], "end": span[1],
                "length": span[1] - span[0] + 1, "kind": kind,
                "transcripts": []})
            rec["transcripts"].append(t.get("id"))
            had = True
        if not had:
            out["n_nmd_without_candidate"] += 1

    cands = sorted(found.values(),
                   key=lambda c: (c["kind"] != "novel exon",
                                  -len(c["transcripts"]), c["start"]))
    for c in cands:
        c["n_transcripts"] = len(c["transcripts"])
    out["candidates"] = cands
    novel = sum(1 for c in cands if c["kind"] == "novel exon")
    out["note"] = (
        f"{len(nmd)} decay isoform(s) of {symbol} in Ensembl, "
        f"{len(nmd) - out['n_nmd_without_candidate']} of them with an exon "
        f"the coding isoforms do not use: {novel} novel exon(s) and "
        f"{len(cands) - novel} extended or shifted coding exon(s).")
    return out
