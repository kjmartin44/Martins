"""Start from the factors known to regulate your gene, then look for drugs.

Why this exists
---------------
Every other factor route starts at the promoter: scan the sequence for motifs,
rank the factors whose motifs match, and ask about drugs for the top few. A
retrospective check against CTD found that most approved drugs recorded as
moving a gene never appeared, and the funnel showed why. Many act on a factor
the literature already lists as a regulator of the gene, and the motif scan did
not rank it, because the factor binds the gene through a site the scan cannot
see, through a partner, or at a distant enhancer.

This route reverses the search. It takes the factors recorded as regulating
the gene and asks what approved drugs act on them.

Where the factors come from
---------------------------
* CollecTRI: curated regulation, with a sign. An edge `FACTOR -> gene` with
  sign +1 says the factor raises the gene.
* Knockout data: experiments in which a factor was lost or raised and your gene
  was among the genes that changed (see `knockout.py`). A gene that falls when
  the factor is lost is one the factor activates. Only factors that are
  transcription factors in CollecTRI are read by default, because a knockout of
  any other gene moves hundreds of genes by indirect routes, and a drug at such
  a gene would be a guess. A factor must also respond in at least
  `MIN_RESPONSES` experiments in the same direction, since a single
  experiment lists the gene by chance about one time in eleven.

Factors the motif scan already ranked and found a drug for are skipped, since
they already have a route.

What it does not say
--------------------
It does not say the factor drives the gene in your cells, or that the drug
moves the factor in them. CollecTRI and CTD draw on the same literature, so a
drug found here and recorded in CTD is not an independent confirmation. The
signs are database annotations, and a factor can act both ways in different
cells. Direction comes from ChEMBL's recorded action at the target. Where that
is "unclear" the drug is listed and marked.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import knockout as knockout_mod

CATEGORY = "TF: known regulator"
SOURCE = "regulators"

COLLECTRI = "CollecTRI"
KNOCKOUT = "knockout data"

# Experiments in which the gene must move the same way before knockout data
# alone names a factor as a regulator.
MIN_RESPONSES = 2
# Drugs kept per factor. The list is already restricted to the direction the
# factor must move, and a factor with dozens of approved drugs (a steroid
# receptor) would otherwise fill the table.
MAX_DRUGS_PER_FACTOR = 10


@dataclass
class Regulator:
    symbol: str
    sign: int                       # +1 the factor raises the gene, -1 lowers it
    sources: list = field(default_factory=list)
    n_responses: int = 0            # knockout experiments that agree
    n_against: int = 0              # knockout experiments that disagree
    disagree: bool = False          # CollecTRI and knockout data give opposite signs
    via: str | None = None          # the recorded regulator this one acts through


def collectri_regulators(gene: str, signed_regulons: dict | None) -> dict:
    """`{factor: sign}` for the factors CollecTRI says regulate `gene`.

    An edge with no recorded mode (sign 0) is left out, because the direction a
    drug must move the factor depends on it.
    """
    gene = (gene or "").upper()
    out = {}
    for tf, targets in (signed_regulons or {}).items():
        sign = targets.get(gene)
        if sign in (1, -1):
            out[tf.upper()] = sign
    return out


def second_step(direct: dict, signed_regulons: dict | None, skip=()) -> list:
    """Factors that regulate a recorded regulator of the gene, in CollecTRI.

    `direct` is `{factor: sign}` for the gene's own regulators. A factor U with
    a signed edge to one of them, T, reaches the gene through T, with the sign
    of the two edges multiplied. Where two routes through different regulators
    disagree on the sign the factor is left out, because the direction a drug
    must move it depends on the sign. Factors already named directly, and those
    in `skip`, are not repeated.
    """
    skip = {s.upper() for s in skip}
    direct = {k.upper(): v for k, v in direct.items()}
    seen: dict = {}
    for u, targets in (signed_regulons or {}).items():
        u = u.upper()
        if u in direct or u in skip:
            continue
        for t, sign_t in direct.items():
            edge = targets.get(t)
            if edge in (1, -1):
                seen.setdefault(u, []).append((t, edge * sign_t))
    out = []
    for u, routes in seen.items():
        signs = {s for _, s in routes}
        if len(signs) != 1:
            continue
        sign = signs.pop()
        via = sorted({t for t, _ in routes})
        out.append(Regulator(u, sign, [COLLECTRI + ", one step removed"],
                             via=", ".join(via[:3])
                             + (f" and {len(via) - 3} more"
                                if len(via) > 3 else "")))
    return sorted(out, key=lambda r: r.symbol)


def candidates(gene: str, signed_regulons: dict | None, *,
               knockout_paths=None, tf_only: bool = True,
               min_responses: int = MIN_RESPONSES, skip=(),
               use_collectri: bool = True, use_knockout: bool = True,
               steps: int = 1) -> list:
    """The regulators of `gene` from both sources, as `Regulator` records.

    `skip` holds factors that already have a route and are left out. The two
    `use_` flags switch a source off, which `scripts/regulator_gain.py` uses to
    score each alone; the table is still read for the list of transcription
    factors when only knockout data is used.
    """
    skip = {s.upper() for s in skip}
    coll = (collectri_regulators(gene, signed_regulons) if use_collectri
            else {})
    found: dict = {}
    for tf, sign in coll.items():
        if tf not in skip:
            found[tf] = Regulator(tf, sign, [COLLECTRI])
    if use_knockout and knockout_mod.resolve_paths(knockout_paths):
        ko = knockout_mod.regulators_of(gene, knockout_paths)
        tfs = {t.upper() for t in (signed_regulons or {})}
        for tf, rec in ko.items():
            if tf in skip or (tf_only and tf not in tfs):
                continue
            n_for = rec["activates"] if rec["sign"] == 1 else rec["represses"]
            n_against = rec["represses"] if rec["sign"] == 1 else rec["activates"]
            if tf in found:
                # CollecTRI already names the factor. The knockout data adds
                # to it, or says the opposite, and either way is recorded.
                r = found[tf]
                r.sources.append(KNOCKOUT)
                r.disagree = bool(rec["sign"]) and rec["sign"] != r.sign
                agree = rec["activates"] if r.sign == 1 else rec["represses"]
                r.n_responses = agree
                r.n_against = (rec["represses"] if r.sign == 1
                               else rec["activates"])
            elif rec["sign"] and n_for >= min_responses:
                found[tf] = Regulator(tf, rec["sign"], [KNOCKOUT],
                                      n_responses=n_for, n_against=n_against)
    out = sorted(found.values(), key=lambda r: r.symbol)
    if steps >= 2 and use_collectri:
        # `coll` holds every direct regulator, including those skipped because
        # they already have a route: a skipped factor is still a way in.
        out += second_step({**coll, **{r.symbol: r.sign for r in out}},
                           signed_regulons, skip=skip)
    return out


def _wanted_drugs(r, goal: str, drug_index: dict, binders=None) -> tuple:
    """`(direction the factor must move, the drugs that move it that way)`."""
    need_up = (r.sign == 1) == (goal == "up")
    want = "up" if need_up else "down"
    drugs = [d for d in drug_index.get(r.symbol, [])
             if d.direction in (want, "unclear")]
    if not drugs and binders is not None:
        drugs = list(binders(r.symbol) or [])
    return want, drugs


def _record(r, want: str) -> dict:
    """The part of a regulator's row that does not depend on any drug."""
    return {"sign": r.sign, "need": want, "sources": list(r.sources),
            "n_responses": r.n_responses, "n_against": r.n_against,
            "disagree": r.disagree, "via": r.via}


def find_drugs(regs: list, goal: str, drug_index: dict, binders=None) -> dict:
    """`{factor: record}` for each regulator that has an approved drug.

    The factor must move in the direction `goal` asks for: up when the factor
    raises the gene and the goal is up, or lowers it and the goal is down.
    Drugs recorded with an unclear action are kept and marked.

    `binders` is an optional function `symbol -> [Drug]` (see
    `drugs.measured_binders`). It is asked only for a regulator the mechanism
    index gives nothing for, and what it returns has direction "unclear".

    The regulators this returns NOTHING for are not uninteresting; see
    `without_drugs`, which is the other half of the same walk.
    """
    out = {}
    for r in regs:
        want, drugs = _wanted_drugs(r, goal, drug_index, binders)
        if not drugs:
            continue
        out[r.symbol] = {
            **_record(r, want),
            "drugs": [d.to_dict() for d in drugs[:MAX_DRUGS_PER_FACTOR]],
            "n_drugs": len(drugs),
        }
    return out


def without_drugs(regs: list, goal: str, drug_index: dict,
                  binders=None) -> list:
    """The recorded regulators that no approved drug reaches.

    `find_drugs` drops these, and dropping them silently is the same
    mistake the evidence ledger exists to prevent: it leaves a report whose
    regulator section lists whatever happened to be druggable and never
    says which factors are actually recorded as running the gene.

    For BUB1B that was the whole finding. CollecTRI names FOXM1, MYC, TP53
    and ZNF143; none has an approved drug, so the section showed AHR, ESR2
    and three HDACs instead - factors that reach the gene far less
    directly and are there because somebody made a ligand for them.

    Returns rows in the same shape as `find_drugs`, with no `drugs` key, so
    a reader and a caller can treat the two together. No route is built
    from these: there is nothing to give anyone. They are the list of what
    to go and find a molecule for.
    """
    out = []
    for r in regs:
        want, drugs = _wanted_drugs(r, goal, drug_index, binders)
        if drugs:
            continue
        out.append({**_record(r, want), "symbol": r.symbol})
    return sorted(out, key=lambda d: d["symbol"])
