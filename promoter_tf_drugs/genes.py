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


# Ensembl returns MANE only when asked. Without `mane=1` the transcripts
# carry `is_canonical`, `ens_canon_extended` and `gencode_primary` and
# nothing MANE-shaped, so `_is_mane` below can never match and the MANE
# branch of `select_transcripts` never runs - every gene quietly takes
# the canonical fallback instead. Canonical is a reasonable answer and
# it is not always the same answer.
_LOOKUP_PARAMS = {"expand": 1, "mane": 1,
                  "content-type": "application/json"}


def _lookup_symbol(client: Client, symbol: str, species: str) -> dict:
    url = f"{ENSEMBL}/lookup/symbol/{species}/{symbol}"
    try:
        return client.get(url, dict(_LOOKUP_PARAMS))
    except HttpError as exc:
        raise HttpError(
            f"Ensembl could not resolve '{symbol}' in {species}. "
            "Check the symbol, or pass a stable ID such as ENSG00000205571."
        ) from exc


def _lookup_id(client: Client, gene_id: str) -> dict:
    url = f"{ENSEMBL}/lookup/id/{gene_id}"
    return client.get(url, dict(_LOOKUP_PARAMS))


def _is_mane(tr: dict) -> bool:
    """Is this the MANE Select transcript?

    Five spellings, because Ensembl has used several and the earlier
    four were added without the request parameter that makes any of
    them arrive. The one the live API actually sends, given `mane=1`,
    is the last: a `MANE` list whose entries carry
    `{"type": "MANE_Select"}`. MANE Plus Clinical is deliberately not
    accepted - it is a second transcript for variant reporting, not the
    representative one.
    """
    if tr.get("mane_select") or tr.get("is_mane_select"):
        return True
    mane = tr.get("MANE") or []
    if isinstance(mane, dict):
        mane = [mane]
    for entry in mane:
        if isinstance(entry, dict) and entry.get("type") == "MANE_Select":
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
    # A local GTF answers this with no network. Ensembl's REST API returns
    # HTTP 500 under load, and this is the call that cannot degrade - without
    # a gene model there is no promoter to analyse.
    gene = local_gene(symbol) if _TX_BY_GENE else None
    if gene is None:
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
                seq = local_sequence_from_genome(chrom, start, end)
                if seq is not None:
                    log.info("using the local genome for %s:%d-%d (no "
                             "Ensembl request)", chrom, start, end)
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

# ====================================================================== #
# Fully offline: gene models from a GTF, sequence from a genome FASTA.    #
#                                                                         #
# Ensembl was the last hard network dependency, and its REST API returns   #
# HTTP 500 under load like every other public service. `--promoter-fasta`  #
# only helps for windows already fetched; a new gene still needed the      #
# network. These two routes remove that entirely.                          #
# ====================================================================== #
_TX_BY_GENE: dict = {}
_TX_SOURCE: dict = {"path": None}
_GENOME: dict = {"path": None, "index": {}}

DEFAULT_GENOME_PATHS = ("data/hg38.fa", "data/genome.fa")


def load_transcripts_gtf(path, *, force: bool = False) -> int:
    """Gene models from a GENCODE/Ensembl GTF, replacing the Ensembl lookup.

    Idempotent: a 34 MB parse takes seconds, and doing it once per
    `pipeline.run()` call took the test suite from 32s to 7 minutes.

    Reads `transcript` rows and keeps the tags transcript selection depends
    on, so `--all-tss` and MANE-Select behaviour are unchanged offline.
    """
    import gzip
    import re

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")
    if _TX_BY_GENE and not force and _TX_SOURCE.get("path") == str(path):
        return len(_TX_BY_GENE)

    opener = gzip.open if path.suffix == ".gz" else open
    gene_re = re.compile(r'gene_name "([^"]+)"')
    tx_re = re.compile(r'transcript_id "([^"]+)"')
    bt_re = re.compile(r'transcript_type "([^"]+)"')
    TAB = chr(9)

    by_gene: dict = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(TAB)
            if len(f) < 9 or f[2] != "transcript":
                continue
            g = gene_re.search(f[8])
            t = tx_re.search(f[8])
            if not (g and t):
                continue
            bt = bt_re.search(f[8])
            by_gene.setdefault(g.group(1).upper(), []).append({
                "id": t.group(1).split(".")[0],
                "start": int(f[3]), "end": int(f[4]),
                "strand": 1 if f[6] == "+" else -1,
                "biotype": bt.group(1) if bt else "unknown",
                "chrom": f[0].replace("chr", ""),
                # the flags select_transcripts relies on
                "is_canonical": 1 if 'tag "Ensembl_canonical"' in f[8] else 0,
                "is_mane_select": 'tag "MANE_Select"' in f[8],
            })

    _TX_BY_GENE.clear()
    _TX_BY_GENE.update(by_gene)
    _TX_SOURCE["path"] = str(path)
    log.info("indexed transcripts for %d genes from %s - gene lookup now "
             "runs with no network", len(by_gene), path.name)
    return len(by_gene)


def local_gene(symbol: str) -> dict | None:
    """An Ensembl-shaped gene record, built from the local GTF."""
    tx = _TX_BY_GENE.get((symbol or "").strip().upper())
    if not tx:
        return None
    return {
        "display_name": symbol.upper(),
        "id": "LOCAL:" + symbol.upper(),
        "seq_region_name": tx[0]["chrom"],
        "strand": tx[0]["strand"],
        "assembly_name": "from GTF",
        "Transcript": tx,
    }


# ------------------------------------------------------------- genome FASTA
def index_genome(fasta) -> dict:
    """Build (or read) a .fai index so sequence can be fetched by seek.

    The .fai format is five columns - name, length, byte offset of the first
    base, bases per line, bytes per line - and generating it here avoids
    requiring samtools just to read a FASTA.
    """
    fasta = Path(fasta)
    fai = fasta.with_suffix(fasta.suffix + ".fai")
    TAB = chr(9)

    if fai.exists():
        index = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            f = line.split(TAB)
            if len(f) >= 5:
                index[f[0].replace("chr", "")] = (int(f[1]), int(f[2]),
                                                  int(f[3]), int(f[4]))
        return index

    log.info("building a FASTA index for %s (one time, a few minutes for a "
             "whole genome)", fasta.name)
    index, rows = {}, []
    with open(fasta, "rb") as fh:
        name, offset, length, linebases, linewidth = None, 0, 0, 0, 0
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith(b">"):
                if name is not None:
                    index[name.replace("chr", "")] = (length, offset,
                                                      linebases, linewidth)
                    rows.append(TAB.join([name, str(length), str(offset),
                                          str(linebases), str(linewidth)]))
                name = line[1:].split()[0].decode()
                offset, length, linebases, linewidth = fh.tell(), 0, 0, 0
            else:
                if linebases == 0:
                    linebases = len(line.rstrip(b"\r\n"))
                    linewidth = len(line)
                length += len(line.rstrip(b"\r\n"))
        if name is not None:
            index[name.replace("chr", "")] = (length, offset, linebases,
                                              linewidth)
            rows.append(TAB.join([name, str(length), str(offset),
                                  str(linebases), str(linewidth)]))
    try:
        fai.write_text(chr(10).join(rows) + chr(10), encoding="utf-8")
    except OSError:
        pass
    return index


def load_genome(fasta) -> int:
    """Make a local genome the sequence source. Idempotent."""
    fasta = Path(fasta)
    if not fasta.exists():
        raise FileNotFoundError(str(fasta) + " not found")
    if _GENOME["path"] == str(fasta) and _GENOME["index"]:
        return len(_GENOME["index"])
    _GENOME["index"] = index_genome(fasta)
    _GENOME["path"] = str(fasta)
    log.info("genome loaded: %d sequence(s) from %s", len(_GENOME["index"]),
             fasta.name)
    return len(_GENOME["index"])


_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def local_sequence_from_genome(chrom, start: int, end: int,
                               strand: int = 1) -> str | None:
    """One window, by seek. Coordinates are 1-based inclusive, as Ensembl's.

    Soft-masking is whatever the FASTA carries: a UCSC hg38.fa is already
    soft-masked, so repeat masking keeps working exactly as it did against
    `mask=soft`.
    """
    idx = _GENOME["index"]
    if not idx or not _GENOME["path"]:
        return None
    key = str(chrom).replace("chr", "")
    meta = idx.get(key)
    if meta is None:
        return None
    length, offset, linebases, linewidth = meta
    if start < 1 or end > length or end < start:
        return None

    begin = offset + (start - 1) // linebases * linewidth \
        + (start - 1) % linebases
    stop = offset + (end - 1) // linebases * linewidth + (end - 1) % linebases

    with open(_GENOME["path"], "rb") as fh:
        fh.seek(begin)
        raw = fh.read(stop - begin + 1)
    seq = raw.decode("ascii", "replace")
    for ch in (chr(10), chr(13)):
        seq = seq.replace(ch, "")

    if strand == -1:
        seq = seq.translate(_COMPLEMENT)[::-1]
    return seq


def autoload_offline(gtf=None, genome=None) -> dict:
    """Pick up a GTF and genome sitting at conventional paths."""
    from . import neighbours as _nbr

    out = {"gtf": False, "genome": False}
    for cand in ([gtf] if gtf else list(_nbr.DEFAULT_GTF_PATHS)):
        if cand and Path(cand).exists():
            try:
                load_transcripts_gtf(cand)
                out["gtf"] = True
                break
            except Exception as exc:                        # noqa: BLE001
                log.warning("could not read %s: %s", cand, exc)
    for cand in ([genome] if genome else list(DEFAULT_GENOME_PATHS)):
        if cand and Path(cand).exists():
            try:
                load_genome(cand)
                out["genome"] = True
                break
            except Exception as exc:                        # noqa: BLE001
                log.warning("could not read %s: %s", cand, exc)
    return out
