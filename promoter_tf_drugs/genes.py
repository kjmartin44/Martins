"""Resolve a gene symbol to a transcription start site and pull the promoter.

Uses the Ensembl REST API (rest.ensembl.org), which needs no key.

Two details that are easy to get wrong and that silently ruin everything
downstream:

1. Strand. For a minus-strand gene the TSS is the transcript *end* coordinate,
   not the start. Get this backwards and you scan the 3' end of the gene and
   find nothing meaningful.
2. Which transcript. Genes have many TSSs. Defaulting to "the first transcript
   in the list" is arbitrary. This module prefers the MANE Select transcript,
   then the Ensembl canonical, then the longest protein-coding one, and it
   records the choice so it appears in the report. If your gene has genuine
   alternative promoters, scan each of them (`--all-tss`).
"""

from __future__ import annotations

import logging
from pathlib import Path
from dataclasses import dataclass, field, asdict

from .http import Client, HttpError

log = logging.getLogger(__name__)

ENSEMBL = "https://rest.ensembl.org"


@dataclass
class Promoter:
    gene_symbol: str
    gene_id: str
    transcript_id: str
    transcript_selection: str          # why this transcript was picked
    chrom: str
    strand: int                        # +1 / -1
    tss: int                           # 1-based genomic TSS
    start: int                         # 1-based inclusive, genomic
    end: int                           # 1-based inclusive, genomic
    upstream: int
    downstream: int
    assembly: str
    sequence: str = field(repr=False, default="")

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    def offset_to_tss(self, index: int) -> int:
        """0-based index into `sequence` -> position relative to TSS.

        The sequence is always returned 5'->3' along the gene, so index 0 is
        `-upstream` and the TSS itself sits at index `upstream` (= position +1
        by convention; there is no position 0).
        """
        rel = index - self.upstream
        return rel + 1 if rel >= 0 else rel

    def to_dict(self) -> dict:
        d = asdict(self)
        d["sequence_length"] = len(self.sequence)
        d.pop("sequence")
        return d


def _lookup_symbol(client: Client, symbol: str, species: str) -> dict:
    url = f"{ENSEMBL}/lookup/symbol/{species}/{symbol}"
    try:
        return client.get(url, {"expand": 1, "content-type": "application/json"})
    except HttpError as exc:
        raise HttpError(
            f"Ensembl could not resolve '{symbol}' in {species}. "
            "Check the symbol, or pass a stable ID such as ENSG00000205571."
        ) from exc


def _lookup_id(client: Client, gene_id: str) -> dict:
    url = f"{ENSEMBL}/lookup/id/{gene_id}"
    return client.get(url, {"expand": 1, "content-type": "application/json"})


def _is_mane(tr: dict) -> bool:
    if tr.get("mane_select") or tr.get("is_mane_select"):
        return True
    tags = tr.get("tag") or tr.get("tags") or []
    if isinstance(tags, str):
        tags = [tags]
    return any("MANE_Select" in str(t) for t in tags)


# Two TSSs closer than this give windows that overlap almost completely, so
# analysing both just duplicates the whole report.
TSS_CLUSTER_BP = 500


def select_transcripts(gene: dict, all_tss: bool = False,
                       cluster_bp: int = TSS_CLUSTER_BP) -> list[tuple[dict, str]]:
    """Return [(transcript, reason)] ordered by preference."""
    transcripts = gene.get("Transcript") or []
    if not transcripts:
        raise ValueError(f"no transcripts returned for {gene.get('id')}")

    if all_tss:
        # Collapse to *distinct* TSS positions. Deduplicating on the exact
        # coordinate is not enough: CCR5 has TSSs 7 bp apart, which produced
        # three promoter windows overlapping 100% and tripled the report for
        # no information. Cluster anything within `cluster_bp`.
        # key on position only: two transcripts can share a TSS, and tuple
        # comparison would then try to order the dicts themselves
        positions = sorted(
            (((tr["start"] if gene["strand"] == 1 else tr["end"]), tr)
             for tr in transcripts),
            key=lambda pair: pair[0])

        out, clusters = [], []
        for pos, tr in positions:
            if clusters and pos - clusters[-1][0] <= cluster_bp:
                clusters[-1][2].append(tr)      # same cluster
                continue
            clusters.append((pos, tr, [tr]))

        for pos, rep, members in clusters:
            # Prefer MANE/canonical as the cluster representative.
            best = next((t for t in members if _is_mane(t)), None) or \
                next((t for t in members if t.get("is_canonical")), None) or rep
            reason = ("distinct TSS" if len(members) == 1
                      else f"distinct TSS (representing {len(members)} within "
                           f"{cluster_bp} bp)")
            out.append((best, reason))

        log.info("%d transcripts collapsed to %d distinct TSS cluster(s)",
                 len(transcripts), len(out))
        return out

    for tr in transcripts:
        if _is_mane(tr):
            return [(tr, "MANE Select")]
    for tr in transcripts:
        if tr.get("is_canonical"):
            return [(tr, "Ensembl canonical")]
    coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
    pool = coding or transcripts
    best = max(pool, key=lambda t: t["end"] - t["start"])
    return [(best, "longest protein-coding" if coding else "longest transcript")]


def fetch_sequence(client: Client, species: str, chrom: str, start: int, end: int,
                   strand: int, mask: str | None = None) -> str:
    """Genomic sequence, already oriented 5'->3' for the requested strand.

    mask="soft" returns repeats in lowercase, which the scanner can exclude.
    """
    url = f"{ENSEMBL}/sequence/region/{species}/{chrom}:{start}..{end}:{strand}"
    params = {"content-type": "application/json"}
    if mask:
        params["mask"] = mask

    try:
        data = client.get(url, params)
    except HttpError as exc:
        if not mask:
            raise HttpError(
                f"Ensembl could not return {chrom}:{start}-{end}. Either the "
                f"service is having problems (try again in a few minutes, or "
                f"check https://rest.ensembl.org) or the region is invalid. "
                f"Original error: {exc}"
            ) from exc
        # Repeat masking is an optimisation, not a requirement. Ensembl
        # intermittently 500s on mask=soft for some regions; losing the whole
        # run over it is the wrong trade.
        log.warning(
            "Ensembl rejected mask=%s for %s:%d-%d (%s). Retrying unmasked — "
            "repeat masking is OFF for this region, so expect more spurious "
            "hits in Alu and simple repeats.", mask, chrom, start, end, exc)
        params.pop("mask", None)
        data = client.get(url, params)
        mask = None

    seq = data["seq"] if mask == "soft" else data["seq"].upper()
    expected = end - start + 1
    if len(seq) != expected:
        log.warning("expected %d bp, Ensembl returned %d", expected, len(seq))
    return seq


def get_promoters(client: Client, symbol: str, *, species: str = "homo_sapiens",
                  upstream: int = 1000, downstream: int = 200,
                  all_tss: bool = False, mask: str | None = None,
                  tss_cluster_bp: int = TSS_CLUSTER_BP) -> list[Promoter]:
    """Fetch promoter window(s) for a gene symbol or Ensembl gene ID."""
    gene = (_lookup_id(client, symbol) if symbol.upper().startswith("ENS")
            else _lookup_symbol(client, symbol, species))

    chrom = gene["seq_region_name"]
    strand = int(gene["strand"])
    assembly = gene.get("assembly_name", "unknown")

    out: list[Promoter] = []
    failures: list[tuple] = []
    for tr, reason in select_transcripts(gene, all_tss=all_tss,
                                         cluster_bp=tss_cluster_bp):
        if strand == 1:
            tss = int(tr["start"])
            start, end = tss - upstream, tss + downstream - 1
        else:
            tss = int(tr["end"])
            start, end = tss - downstream + 1, tss + upstream
        start = max(start, 1)

        # One unreachable window must not cost the whole run. With --all-tss
        # a gene has several, and a transient Ensembl 500 on one of them used
        # to abort after the transcript selection, regulon and expression
        # tables had already been built. Same principle as
        # `pipeline._optional`: lose the section, not the run.
        try:
            # A local FASTA wins outright: faster, reproducible, and immune
            # to the endpoint being down.
            seq = local_sequence(chrom, start, end)
            if seq is None:
                seq = fetch_sequence(client, species, chrom, start, end,
                                     strand, mask=mask)
            else:
                log.info("using the local FASTA for %s:%d-%d (no Ensembl "
                         "request)", chrom, start, end)
        except Exception as exc:                            # noqa: BLE001
            failures.append((tr["id"], f"{chrom}:{start}-{end}",
                             f"{type(exc).__name__}: {str(exc)[:120]}"))
            log.warning(
                "could not fetch the promoter sequence for %s (%s:%d-%d): "
                "%s — skipping this TSS and continuing with the others",
                tr["id"], chrom, start, end, str(exc)[:160])
            continue

        out.append(Promoter(
            gene_symbol=gene.get("display_name", symbol),
            gene_id=gene["id"],
            transcript_id=tr["id"],
            transcript_selection=reason,
            chrom=chrom,
            strand=strand,
            tss=tss,
            start=start,
            end=end,
            upstream=upstream,
            downstream=downstream,
            assembly=assembly,
            sequence=seq,
        ))

    if failures and not out:
        # Every window failed: there is genuinely nothing to analyse, so this
        # one does raise rather than returning an empty promoter list that
        # would read downstream as "this gene has no promoter".
        detail = "; ".join(f"{t} ({r}): {e}" for t, r, e in failures[:3])
        raise HttpError(
            f"could not fetch any promoter sequence for {symbol} "
            f"({len(failures)} window(s) attempted). {detail}. Ensembl "
            f"returns HTTP 500 intermittently — everything already fetched "
            f"is cached, so re-running resumes quickly.")
    if failures:
        log.warning(
            "%d of %d promoter window(s) for %s could not be fetched and "
            "were SKIPPED: %s. The analysis below covers only the %d window(s) "
            "that succeeded.", len(failures), len(failures) + len(out),
            symbol, ", ".join(r for _, r, _ in failures), len(out))
    return out


# ------------------------------------------------- offline sequence source
# Ensembl was the one stage with no local route, so a 500 there took the
# whole run down. But the package already WRITES every window it fetched to
# promoter.fa, so a previous run is its own offline source — and any FASTA
# with matching coordinates works.
_LOCAL_SEQ: dict[str, str] = {}


def load_promoter_fasta(path) -> int:
    """Index a promoter FASTA so sequence fetches can run with no network.

    Accepts this package's own `promoter.fa` output, whose headers look like

        >GENE|ENST00000292303|chr3:46365946-46371945|strand=1|-5000/+1000

    and any FASTA whose header contains a `chrN:start-end` span. Keyed on the
    span, so a window is reused only when its coordinates match exactly —
    a near-miss would silently analyse the wrong sequence.
    """
    import re
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")

    span = re.compile(r"(?:chr)?([\w.]+):(\d+)-(\d+)")
    name, buf, n = None, [], 0
    _LOCAL_SEQ.clear()

    def _store(hdr, chunks):
        nonlocal n
        if not hdr or not chunks:
            return
        m = span.search(hdr)
        if not m:
            return
        key = f"{m.group(1)}:{m.group(2)}-{m.group(3)}"
        _LOCAL_SEQ[key] = "".join(chunks)
        n += 1

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            _store(name, buf)
            name, buf = line[1:], []
        elif line.strip():
            buf.append(line.strip())
    _store(name, buf)

    log.info("indexed %d promoter sequence(s) from %s — these will be used "
             "instead of Ensembl where the coordinates match", n, path.name)
    return n


def local_sequence(chrom: str, start: int, end: int) -> str | None:
    if not _LOCAL_SEQ:
        return None
    c = str(chrom).replace("chr", "")
    return _LOCAL_SEQ.get(f"{c}:{start}-{end}")


def fetch_cdna(client: Client, transcript_id: str) -> str | None:
    """Spliced transcript sequence, for siRNA/ASO design."""
    try:
        data = client.get(f"{ENSEMBL}/sequence/id/{transcript_id}",
                          {"type": "cdna", "content-type": "application/json"})
    except HttpError as exc:
        log.warning("cDNA fetch failed for %s: %s", transcript_id, exc)
        return None
    return (data.get("seq") or "").upper() or None
