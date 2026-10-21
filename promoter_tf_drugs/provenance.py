"""Did each evidence source actually produce evidence?

The failure this exists to prevent
-----------------------------------
Every optional analysis in this package degrades gracefully: a dead API costs
you that section, not the run (`pipeline._optional`). That is the right
behaviour for *robustness* and exactly the wrong behaviour for *inference*,
because the scoring layer downstream cannot tell these two apart:

    the lookup ran and found nothing        → evidence of absence
    the lookup broke and returned nothing   → absence of evidence

Both arrive as an empty dict. A run where ChIP, tissue expression and LINCS all
failed silently still prints grades, Pareto fronts and a tier-1 row, and reads
exactly like a run where all three succeeded. That report is more dangerous
than no report, because its confident structure is intact while the empirical
content underneath it is gone.

So every source records what it was *asked* to do and what it *managed* to do.
A source that was requested and failed is a loud, named failure at the top of
the report, and — where the distinction changes a number — it is refused as
scoring input rather than passed through as a neutral value.

`requested=False` is not a failure. It means you did not ask for that source,
which is a choice, not a fault, and it is reported as "not requested" so an
absent section is never mistaken for a broken one.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


# Each source names the flag that turns it on, so a "not requested" line tells
# you what to add rather than merely what is missing.
FLAGS = {
    "chip": "--chip-bed (preferred) or --chip-api",
    "expression": "--cell-type-check, ideally with --expression-file",
    "lincs": "--lincs or --lincs-file",
    "ctd": "on by default; --no-improvements turns it off",
    "genetics": "--genetics",
    "regulon": "--regulon-file",
    "upstream": "--upstream-network",
    "screen": "--screen",
    "direction_model": "--predict-direction (needs rdkit + scikit-learn)",
    "knockout": "--knockout-gmt, or data/knockout_*.gmt from Enrichr",
    "regulators": "--regulator-route (needs CollecTRI in data/)",
}

# What you lose when a source is unavailable. Stated per source because
# "degraded" on its own does not tell anyone whether to rerun.
IMPACT = {
    "chip": "Nothing measured whether these factors actually sit on this "
            "promoter, so the ranking rests on sequence resemblance "
            "alone. This is the single largest loss of evidence in a "
            "run - most of a motif hit list is chance.",
    "expression": "Factors were not checked for being switched on in the "
                  "same tissue as your gene, so a regulator that is "
                  "absent from your cells can still rank at the top.",
    "lincs": "No measured transcriptional response for this gene, so no "
             "row here is backed by a compound having actually moved it "
             "in an experiment. Check the CTD row before concluding "
             "nothing was ever measured.",
    "ctd": "The curated literature on chemicals that change this gene's "
           "expression could not be read. Two consequences: the report "
           "loses its only literature-derived evidence, and nothing can "
           "be called novel, because a link that was not retrieved is "
           "not a link that does not exist.",
    "genetics": "No human genetic evidence that the factor-to-gene link "
                "is real, so the connection rests on annotation.",
    "regulon": "How many genes each factor controls is unknown, so a "
               "factor that runs half the transcriptome is not being "
               "penalised for it and may rank alongside a specific one.",
    "upstream": "Factors with no drug acting on them directly return no "
                "route at all, so a reachable target two steps away is "
                "simply missing rather than rejected.",
    "screen": "No search for untested drug and receptor pairs, so "
              "anything not already annotated is absent rather than "
              "ruled out.",
    "direction_model": "Whether each drug activates or blocks its target "
                       "was taken from the database label alone, with no "
                       "independent structural check - so a mislabelled "
                       "entry passes through, and direction is what "
                       "decides which way your gene moves.",
    "tf_activity": "Nothing checked whether these factors are actually "
                   "DOING anything in your cells. A motif match says a "
                   "sequence resembles a binding site; this is the step "
                   "that would have asked whether the factor that "
                   "matches it is active.",
    "docking": "No structure was examined, so no row is supported by a "
               "compound physically fitting a pocket. Everything here "
               "rests on what databases say rather than on shape.",
    "product_screen": "The gene's own protein was not screened for "
                      "binders. A compound that acts on the protein "
                      "rather than on its transcription is absent "
                      "rather than ruled out - and that is often the "
                      "shortest route to changing what the gene DOES.",
    "promoter_dna": "The promoter sequence was not examined for "
                    "features a small molecule could bind directly, so "
                    "DNA-directed agents are absent rather than "
                    "excluded. Where this does run and finds nothing, it "
                    "means a molecule aimed at such a feature would have "
                    "no specific handle on this promoter - the same "
                    "feature occurs at the same density in any sequence "
                    "of that base composition.",
    "factor_abundance": "No check on whether a compound changes how MUCH "
                        "of the factor exists, as against what it does. "
                        "That is the one mechanism here a western blot "
                        "would see.",
    "trans_eqtl": "No human genetic evidence that natural variation in a "
                  "factor moves this gene, so the factor-to-gene step "
                  "rests on annotation rather than on people.",
    "regulators": "Factors recorded as regulating your gene were not "
                  "checked for drugs, so a factor the motif scan did not "
                  "rank cannot supply a route.",
    "knockout": "Nothing measured whether your gene changes when a factor "
                "is lost, so responsiveness stays unmeasured for every "
                "factor route.",
    "activity_drugs": "No compound was tested against inferred factor "
                      "activity, so that evidence generates no routes of "
                      "its own.",
    "depmap": "What a cell loses without each factor was not measured, "
              "so a factor that every cell line needs is not marked as "
              "a costly target and can rank beside a dispensable one. "
              "Cancer cell lines say little about a post-mitotic tissue "
              "even when this does run.",
    "sensitivity": "No drug was matched against a factor's knockout "
                   "profile, so the measured phenocopy route is absent: a "
                   "compound that behaves as though it removes the factor "
                   "but has no annotation anywhere stays unfound.",
}


@dataclass
class Source:
    """One evidence source's outcome for one run."""
    name: str
    requested: bool = False
    ok: bool = False
    detail: str = ""
    yielded: int = 0

    @property
    def state(self) -> str:
        """One of: not requested, FAILED, empty, ok.

        `empty` is the state the review asked for and the one most
        easily lost: the source was reached, it answered, and the answer
        was nothing. That is a measurement and it is not a find, and
        reporting it as plain `ok` puts it in the same column as a
        source that returned a hundred records.

        It remains `ok` for every decision that turns on the source
        having worked - `usable` and `empirical_sources_ok` both read
        `self.ok`, not this - because a measured nothing outranks a
        lookup that broke, which is the distinction the whole ledger
        exists to keep.
        """
        if not self.requested:
            return "not requested"
        if not self.ok:
            return "FAILED"
        if not self.yielded:
            return "empty"
        return "ok"

    @property
    def usable(self) -> bool:
        """May this source's result be used as scoring input?

        A requested-but-failed source is not usable. Refusing it here is what
        stops an empty result being scored as a measured zero.
        """
        return self.requested and self.ok

    def to_dict(self) -> dict:
        d = {
            "name": self.name,
            "requested": self.requested,
            "ok": self.ok,
            "state": self.state,
            "detail": self.detail,
            "yielded": self.yielded,
            "enable_with": FLAGS.get(self.name),
        }
        if self.requested and not self.ok:
            d["impact"] = IMPACT.get(self.name, "")
        return d


@dataclass
class Ledger:
    """Every evidence source's outcome, collected across a run."""
    sources: dict = field(default_factory=dict)

    def record(self, name: str, *, requested: bool, ok: bool,
               detail: str = "", yielded: int = 0) -> Source:
        s = Source(name=name, requested=requested, ok=ok, detail=detail,
                   yielded=yielded)
        # Across several TSSs the same source is recorded repeatedly. Keep the
        # best outcome: one promoter failing a lookup that succeeded elsewhere
        # is a per-region result, not a broken source.
        prev = self.sources.get(name)
        if prev is None or (s.usable and not prev.usable):
            self.sources[name] = s
        elif prev.usable and s.usable:
            prev.yielded = max(prev.yielded, s.yielded)
        return self.sources[name]

    def get(self, name: str) -> Source:
        return self.sources.get(name) or Source(name=name)

    def usable(self, name: str) -> bool:
        return self.get(name).usable

    def failures(self) -> list:
        return [s for s in self.sources.values() if s.requested and not s.ok]

    # Sources carrying MEASUREMENT rather than annotation. CTD belongs here:
    # a curator recorded that someone dosed cells and the gene's expression
    # moved. That is an observation, not an inference from a motif — and it
    # is frequently the only empirical source still standing when LINCS is
    # unreachable, which it often is.
    EMPIRICAL = ("chip", "lincs", "ctd", "genetics")

    def empirical_sources_ok(self) -> int:
        return sum(1 for n in self.EMPIRICAL if self.usable(n))

    def empirical_names(self) -> list:
        return [n for n in self.EMPIRICAL if self.usable(n)]

    def to_dict(self) -> dict:
        return {
            "sources": {n: s.to_dict() for n, s in sorted(self.sources.items())},
            "n_failed": len(self.failures()),
            "failed": sorted(s.name for s in self.failures()),
            "empirical_sources_ok": self.empirical_sources_ok(),
            "empirical_available": self.empirical_names(),
        }


def verdict(ledger: Ledger, scan_quality: dict | None = None) -> dict:
    """Can this run's ranking be read as evidence, or only as a reading order?

    Two independent things have to hold. The scan has to have found signal
    above chance, and at least one source of measurement has to have worked.
    Neither is sufficient alone: an enriched scan with no ChIP is still
    unconfirmed sequence, and perfect ChIP over a noise scan confirms binding
    at sites you had no reason to look at.
    """
    fdp = (scan_quality or {}).get("false_discovery_proportion")
    fdp_range = (scan_quality or {}).get("false_discovery_range")
    n_sites = (scan_quality or {}).get("observed")
    scan_verdict = (scan_quality or {}).get("verdict")
    n_emp = ledger.empirical_sources_ok()
    failed = [s.name for s in ledger.failures()]

    scan_bad = scan_verdict in ("no signal", "weak")
    reasons = []
    if scan_bad:
        extra = (f" — about {fdp:.0%} of the hit list is expected to be chance"
                 if isinstance(fdp, (int, float)) else "")
        if extra and fdp_range:
            extra += f" (range {fdp_range[0]:.0%} to {fdp_range[1]:.0%})"
        reasons.append(f"the motif scan is '{scan_verdict}'{extra}")
    if n_emp == 0:
        reasons.append("no source of measured evidence is available "
                       "(ChIP, LINCS, genetics all absent or failed)")
    if failed:
        reasons.append("these tests could not be completed: "
                       + ", ".join(failed)
                       + " — which is not the same as their having "
                         "returned a negative result")

    have = ledger.empirical_names()
    if scan_bad and n_emp == 0:
        level = "not interpretable"
        headline = (
            "This run cannot identify a target. The candidate list is a "
            "ranking of chance motif matches with nothing measured to test "
            "it against.")
    elif reasons:
        level = "degraded"
        headline = ("Read the ranking as a reading order, not as evidence. "
                    "At least one pillar of the analysis is missing."
                    + (f" Measured evidence that DID arrive: "
                       f"{', '.join(have)}." if have else ""))
    else:
        level = "interpretable"
        headline = (f"Scan is enriched and measured evidence is available "
                    f"({', '.join(have)}). The ranking rests on more than "
                    f"annotation.")

    return {"level": level, "headline": headline, "reasons": reasons,
            "n_empirical_sources": n_emp,
            "empirical_available": have,
            "scan_verdict": scan_verdict,
            "scan_sites": n_sites,
            "false_discovery_proportion": fdp,
            "false_discovery_range": fdp_range}


# ===================================================================== #
# What produced this result.                                            #
# ===================================================================== #
# The result schema's version. Bump it when the shape changes in a way
# that something reading an older result would get wrong - a renamed
# field, a changed unit, a different meaning for the same key. Adding a
# key is not such a change.
RESULT_SCHEMA = 1


def check_schema(data: dict, tool: str) -> str | None:
    """Warn when a saved result was written by a newer package.

    Returns the warning text, and prints it to stderr. None when the
    file is one this version understands - including a file with no
    schema at all, which predates the field and is read quietly rather
    than warned about every time.
    """
    import sys

    found = (data.get("provenance") or {}).get("schema")
    if found is None or found <= RESULT_SCHEMA:
        return None
    version = (data.get("provenance") or {}).get("package_version")
    msg = (f"{tool}: this result uses schema {found} and this package "
           f"understands {RESULT_SCHEMA}"
           + (f" (written by version {version})" if version else "")
           + ". Reading it anyway; fields added since are ignored, and "
             f"any whose meaning changed are read with the older one. "
             f"Upgrade the package to read it as it was written.")
    print(msg, file=sys.stderr)
    return msg

# Above this, a file's contents are not hashed. The ReMap BED is 1.4 GB
# and the full BindingDB release is 8.6 GB, and hashing either on every
# run costs real time for a number that is almost never compared.
HASH_LIMIT_BYTES = 64 * 1024 * 1024


def describe_input(path) -> dict:
    """Identify one input file: size, modification time, and a hash.

    The hash is skipped above `HASH_LIMIT_BYTES`, and the record says
    so. A missing hash with no explanation reads like a file nobody
    looked at, which is a different claim from one too large to bother
    with.
    """
    import datetime as _dt
    import hashlib

    p = Path(path)
    out: dict = {"path": str(path)}
    try:
        st = p.stat()
    except OSError as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    out["bytes"] = st.st_size
    out["modified"] = _dt.datetime.fromtimestamp(
        st.st_mtime, _dt.timezone.utc).isoformat(timespec="seconds")

    if st.st_size > HASH_LIMIT_BYTES:
        out["sha256"] = None
        out["hash_skipped"] = (
            f"{st.st_size / 1048576:.0f} MB is over the "
            f"{HASH_LIMIT_BYTES // 1048576} MB limit, so size and "
            f"modification time identify it here")
        return out

    try:
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        out["sha256"] = h.hexdigest()
    except OSError as exc:
        out["sha256"] = None
        out["hash_skipped"] = f"unreadable: {type(exc).__name__}"
    return out


def run_provenance(inputs: dict | None = None) -> dict:
    """Package version, when it ran, and which files it read.

    `inputs` maps an argument name to the path it was given. Arguments
    that were not given are left out rather than recorded as null: the
    absence of a file is already in the query block, and repeating it
    here as an input would imply one was read.
    """
    import datetime as _dt
    import platform
    import sys

    from . import __version__

    described = {}
    for name, path in (inputs or {}).items():
        if path:
            described[name] = describe_input(path)

    return {
        "schema": RESULT_SCHEMA,
        "package_version": __version__,
        "created": _dt.datetime.now(_dt.timezone.utc).isoformat(
            timespec="seconds"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "inputs": described,
    }
