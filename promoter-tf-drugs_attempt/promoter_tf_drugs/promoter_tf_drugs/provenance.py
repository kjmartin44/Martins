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


# Each source names the flag that turns it on, so a "not requested" line tells
# you what to add rather than merely what is missing.
FLAGS = {
    "chip": "--chip-bed (preferred) or --chip-api",
    "expression": "--cell-type-check, ideally with --expression-file",
    "lincs": "--lincs or --lincs-file",
    "genetics": "--genetics",
    "regulon": "--regulon-file",
    "upstream": "--upstream-network",
    "screen": "--screen",
}

# What you lose when a source is unavailable. Stated per source because
# "degraded" on its own does not tell anyone whether to rerun.
IMPACT = {
    "chip": "TF ranking rests on sequence alone; motif hits are not confirmed "
            "as bound. This is the single largest loss of evidence in the run.",
    "expression": "Factors are not checked for co-expression with your gene, "
                  "so a regulator absent from your tissue can still rank top.",
    "lincs": "No measured expression change anywhere in the report — every "
             "remaining route is inference.",
    "genetics": "No human causal evidence for the TF-to-gene link.",
    "regulon": "Transcriptional specificity falls back to coarse bundled "
               "tiers, so a global regulator may not be penalised enough.",
    "upstream": "TFs with no direct pharmacology return no route at all.",
    "screen": "No untested drug/receptor pairs are surfaced.",
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
        if not self.requested:
            return "not requested"
        if self.ok:
            return "ok"
        return "FAILED"

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

    def empirical_sources_ok(self) -> int:
        """How many sources carrying *measurement* rather than annotation."""
        return sum(1 for n in ("chip", "lincs", "genetics")
                   if self.usable(n))

    def to_dict(self) -> dict:
        return {
            "sources": {n: s.to_dict() for n, s in sorted(self.sources.items())},
            "n_failed": len(self.failures()),
            "failed": sorted(s.name for s in self.failures()),
            "empirical_sources_ok": self.empirical_sources_ok(),
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
    scan_verdict = (scan_quality or {}).get("verdict")
    n_emp = ledger.empirical_sources_ok()
    failed = [s.name for s in ledger.failures()]

    scan_bad = scan_verdict in ("no signal", "weak")
    reasons = []
    if scan_bad:
        extra = (f" — about {fdp:.0%} of the hit list is expected to be chance"
                 if isinstance(fdp, (int, float)) else "")
        reasons.append(f"the motif scan is '{scan_verdict}'{extra}")
    if n_emp == 0:
        reasons.append("no source of measured evidence is available "
                       "(ChIP, LINCS, genetics all absent or failed)")
    if failed:
        reasons.append("requested source(s) failed: " + ", ".join(failed))

    if scan_bad and n_emp == 0:
        level = "not interpretable"
        headline = (
            "This run cannot identify a target. The candidate list is a "
            "ranking of chance motif matches with nothing measured to test "
            "it against.")
    elif reasons:
        level = "degraded"
        headline = ("Read the ranking as a reading order, not as evidence. "
                    "At least one pillar of the analysis is missing.")
    else:
        level = "interpretable"
        headline = ("Scan is enriched and at least one measured source is "
                    "available. The ranking rests on more than annotation.")

    return {"level": level, "headline": headline, "reasons": reasons,
            "n_empirical_sources": n_emp,
            "scan_verdict": scan_verdict,
            "false_discovery_proportion": fdp}
