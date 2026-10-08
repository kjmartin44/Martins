"""Approved treatments for the condition a gene causes.

The rest of this package asks which drugs would change a gene's expression.
That is a different question from what is already prescribed for the
disease, and a reader looking at a report for HTT reasonably expects both
answers. Tetrabenazine, deutetrabenazine and valbenazine are approved for
the chorea of Huntington's disease. None of them acts on HTT: each blocks
the vesicular monoamine transporter (VMAT2, gene SLC18A2) and treats the
movement symptom. So they never appear in a ranking of drugs that lower
HTT, and no search that starts from the gene can find them.

This table lists them beside the ranking, as context. It is not part of any
score, and it is not a route to the gene.

It is hand-maintained and short, and it is not read from a live source on
purpose. Open Targets was checked for Huntington disease: it records
tetrabenazine and deutetrabenazine as approved, leaves out valbenazine, and
lists four ADHD drugs (atomoxetine, guanfacine, lisdexamfetamine and
methylphenidate) at the same stage. A table pulled from it would miss one of
the three drugs and add four that are not treatments for the disease.

A gene that is not listed has no entry here. That is not a statement that
the condition has no approved treatment.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class Treatment:
    generic: str
    brand: str
    approved: str               # year and month the indication was approved
    use: str                    # what it is approved to treat
    target_gene: str            # the gene whose product the drug acts on
    mechanism: str


@dataclass(frozen=True)
class Condition:
    name: str
    genes: tuple                # genes whose changes cause or drive it
    treatments: tuple
    source: str


# Approval dates are the ones given by the project owner. Confirm each against
# Drugs@FDA before quoting it: valbenazine was first approved in 2017 for
# tardive dyskinesia, and the date below is for Huntington's chorea.
CONDITIONS: tuple = (
    Condition(
        name="Huntington disease",
        genes=("HTT",),
        treatments=(
            Treatment("tetrabenazine", "Xenazine", "2008-08",
                      "chorea of Huntington's disease", "SLC18A2",
                      "VMAT2 inhibitor"),
            Treatment("deutetrabenazine", "Austedo", "2017-04",
                      "chorea of Huntington's disease", "SLC18A2",
                      "VMAT2 inhibitor"),
            Treatment("valbenazine", "Ingrezza", "2023-08",
                      "chorea of Huntington's disease", "SLC18A2",
                      "VMAT2 inhibitor"),
        ),
        source="FDA approvals, as supplied by the project owner",
    ),
)


def for_gene(symbol: str | None) -> dict | None:
    """The condition entry for a gene, or None when none is recorded.

    `acts_on_gene` is True for a treatment whose target is the gene itself.
    For every entry above it is False, and the note says so: these drugs
    treat the condition without touching the gene.
    """
    sym = (symbol or "").strip().upper()
    if not sym:
        return None
    for cond in CONDITIONS:
        if sym in cond.genes:
            rows = []
            for t in cond.treatments:
                d = asdict(t)
                d["acts_on_gene"] = t.target_gene.upper() == sym
                rows.append(d)
            direct = [r["generic"] for r in rows if r["acts_on_gene"]]
            verb = "acts" if len(direct) == 1 else "act"
            note = (f"{', '.join(direct)} {verb} on {sym} itself." if direct else
                    f"None of these acts on {sym}: they treat the condition "
                    f"through {', '.join(sorted({t.target_gene for t in cond.treatments}))}"
                    f", so none of them lowers {sym}.")
            return {"condition": cond.name, "gene": sym, "treatments": rows,
                    "note": note, "source": cond.source}
    return None
