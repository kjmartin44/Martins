"""What else lives in the promoter window you just scanned?

The gap this closes
-------------------
Every other specificity check in this package asks "does this drug, or this
transcription factor, affect other genes". None of them asks the prior
question: **is the sequence I scanned uniquely my gene's?**

Found by inspecting a real CCR5 run. The scanned window
`chr3:46,365,946-46,371,945` sits entirely inside `CCR5AS`, an antisense
lncRNA spanning `chr3:46,363,984-46,407,184` on the opposite strand. `CCR2`,
CCR5's co-regulated paralog, is 5 kb beyond the edge — a `--upstream 10000`
would have swallowed it whole. Nothing in the pipeline said so.

Why it changes the interpretation
---------------------------------
**Antisense overlap breaks the direction logic.** The pipeline reasons
"activator bound here, so block it to lower the gene". That chain assumes the
factor acts on *your* promoter. Where a window lies inside an antisense
transcription unit, a site may drive that transcript instead, and an antisense
RNA can stabilise or repress the sense gene — so the sign does not carry.

**A neighbouring paralog is usually co-regulated.** Overlapping a paralog's
regulatory region means a shared TF is expected, not evidence of specificity.

What this is not
----------------
It checks *annotation*, not sequence uniqueness. A promoter can be unique in
the genome and still sit inside another gene, and it can be free of
neighbours while sharing sequence with a duplicated locus elsewhere. For true
uniqueness, BLAT the window — that is a different question and this module
does not answer it.

For CCR5 that distinction mattered and came out well: CCR5 and CCR2 promoters
share almost no sequence (12-mer Jaccard 0.0005, 6 shared 12-mers against 0
for a shuffled control) despite the genes being paralogs. The risk there was
neighbouring transcription units, not duplicated sequence.
"""

from __future__ import annotations

import logging

log = logging.getLogger(__name__)

ENSEMBL = "https://rest.ensembl.org"

# How far beyond the window to look. A neighbour just outside is worth
# knowing about, because widening the window is the first thing people try
# when a scan comes back weak.
DEFAULT_PADDING = 10_000


def _relation(gene: dict, start: int, end: int, strand: int) -> tuple[str, str]:
    """How this gene sits relative to the scanned window."""
    gs, ge, gstrand = gene["start"], gene["end"], gene["strand"]
    overlaps = not (ge < start or gs > end)
    antisense = gstrand != strand

    if not overlaps:
        gap = start - ge if ge < start else gs - end
        return "nearby", (f"{gap / 1000:.1f} kb outside the window — a wider "
                          f"--upstream/--downstream would include it")
    if gs <= start and ge >= end:
        base = "window lies inside this gene"
    else:
        base = "partially overlaps the window"
    if antisense:
        return "antisense overlap", (
            base + " on the OPPOSITE strand. A site called here may drive "
            "that transcript instead, and an antisense RNA can stabilise or "
            "repress your gene — so the pipeline's direction logic does not "
            "reliably carry")
    return "same-strand overlap", (
        base + " on the same strand; sites here may belong to its regulation "
        "rather than your gene's")


def genes_in_window(client, chrom: str, start: int, end: int, *,
                    strand: int = 1, species: str = "homo_sapiens",
                    padding: int = DEFAULT_PADDING,
                    exclude: str | None = None) -> dict:
    """Every annotated gene overlapping or adjacent to the scanned window.

    One cached Ensembl call per promoter, through the existing Client, so a
    re-run costs nothing.
    """
    c = str(chrom).replace("chr", "")
    lo, hi = max(1, start - padding), end + padding
    url = f"{ENSEMBL}/overlap/region/{species}/{c}:{lo}-{hi}"
    try:
        rows = client.get(url, params={"feature": "gene"})
    except Exception as exc:                                # noqa: BLE001
        log.warning("neighbour lookup failed for %s: %s", url, exc)
        return {"checked": False, "reason": str(exc)[:160], "genes": []}

    if not isinstance(rows, list):
        return {"checked": False, "reason": "unexpected response shape",
                "genes": []}

    out = []
    for g in rows:
        name = g.get("external_name") or g.get("id")
        if exclude and name and name.upper() == exclude.upper():
            continue
        kind, note = _relation(g, start, end, strand)
        out.append({
            "gene": name, "biotype": g.get("biotype"),
            "chrom": c, "start": g["start"], "end": g["end"],
            "strand": g["strand"], "relation": kind, "note": note,
        })
    out.sort(key=lambda r: (r["relation"] != "antisense overlap", r["start"]))

    inside = [r for r in out if r["relation"] != "nearby"]
    antisense = [r for r in out if r["relation"] == "antisense overlap"]
    return {
        "checked": True,
        "window": f"chr{c}:{start}-{end}",
        "padding": padding,
        "genes": out,
        "n_inside": len(inside),
        "n_antisense": len(antisense),
        "clean": not inside,
        "summary": _summarise(inside, antisense, out),
    }


def _summarise(inside, antisense, all_rows) -> str:
    if not all_rows:
        return ("No other annotated gene within the window or its "
                "surroundings — the scanned sequence is this gene's alone, "
                "as far as annotation goes.")
    if not inside:
        near = ", ".join(r["gene"] for r in all_rows[:3] if r["gene"])
        return (f"No other gene inside the window. Nearby: {near}. Widening "
                f"the window would start to include them.")
    if antisense:
        names = ", ".join(r["gene"] for r in antisense if r["gene"])
        return (f"The window overlaps {len(inside)} other transcription "
                f"unit(s), including {names} on the opposite strand. Sites "
                f"called here are not unambiguously your gene's, and for the "
                f"antisense overlap the direction logic does not carry.")
    names = ", ".join(r["gene"] for r in inside if r["gene"])
    return (f"The window overlaps {len(inside)} other gene(s): {names}. "
            f"Sites called here may belong to their regulation.")


def window_overlap(promoters) -> list:
    """How much the --all-tss windows of one gene duplicate each other.

    CCR5's two TSS windows overlapped by 4,275 of 6,000 bp — 71%. Treating
    them as two independent promoters double-counts one analysis, which
    inflates every per-gene route count downstream.
    """
    out = []
    for i in range(len(promoters)):
        for j in range(i + 1, len(promoters)):
            a, b = promoters[i], promoters[j]
            if a.get("chrom") != b.get("chrom"):
                continue
            ov = min(a["end"], b["end"]) - max(a["start"], b["start"])
            if ov <= 0:
                continue
            span = min(a["end"] - a["start"], b["end"] - b["start"]) or 1
            out.append({
                "a": a.get("transcript_id"), "b": b.get("transcript_id"),
                "overlap_bp": ov, "fraction": round(ov / span, 3),
                "note": (f"these two TSS windows share {ov} bp "
                         f"({ov / span:.0%}) — largely one analysis counted "
                         f"twice, which inflates this gene's route count"),
            })
    return out
