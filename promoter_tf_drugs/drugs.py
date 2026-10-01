"""TF -> approved drug lookup, filtered to pharmacology that moves the TF in
the direction you want.

Sources
  ChEMBL  /mechanism gives a curated action_type (AGONIST, ANTAGONIST,
          POSITIVE ALLOSTERIC MODULATOR, ...) plus max_phase on the molecule.
          max_phase == 4 is the approval flag.
  DGIdb   broader, noisier, aggregates DrugBank/ChEMBL/TTD/etc. Useful for
          recall; treat its interaction types as a hint, not a fact.

Direction logic
  To raise output from the target gene you want the *promoter-bound factor's
  activating function increased*. That means:
      activator TF  -> agonist / positive modulator / activator
      repressor TF  -> antagonist / inverse agonist / degrader
  Both cases are computed, and a drug is only reported with a direction that
  is actually consistent with the annotated role of that TF.

The honest expectation: of ~1600 human TFs, only the ~48 nuclear receptors and
a handful of others (AHR, HIF1A via PHD inhibitors, NFE2L2 via KEAP1
inhibitors, HSF1) have approved direct pharmacology. Empty results for a
bZIP or a zinc finger are the normal outcome, not a bug.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict

from .http import Client, HttpError

log = logging.getLogger(__name__)

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
DGIDB_GRAPHQL = "https://dgidb.org/api/graphql"

ACTIVATING_ACTIONS = {
    "AGONIST", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR",
    "POSITIVE MODULATOR", "ACTIVATOR", "OPENER", "STABILISER",
}
INHIBITING_ACTIONS = {
    "ANTAGONIST", "INVERSE AGONIST", "NEGATIVE ALLOSTERIC MODULATOR",
    "NEGATIVE MODULATOR", "INHIBITOR", "BLOCKER", "DEGRADER",
    "DISRUPTING AGENT", "SEQUESTERING AGENT",
}


@dataclass
class Drug:
    name: str
    chembl_id: str | None
    target_symbol: str
    action_type: str
    mechanism: str | None
    approved: bool
    max_phase: float | None
    first_approval: int | None
    withdrawn: bool
    oral: bool | None
    direction: str          # "up" | "down" | "unclear" - compared
                            # directly against desired_direction()
    source: str
    target_chembl_id: str | None = None   # needed to look up potency vs THIS target
    properties: dict | None = None        # molecule_properties, for tissue scoring
    interaction_score: float | None = None  # DGIdb confidence, if from DGIdb
    sources: str | None = None            # which databases claim this

    def to_dict(self) -> dict:
        return asdict(self)


def desired_direction(role: str, goal: str = "up") -> str:
    """Which way to push this TF to move the target gene in direction `goal`.

        goal="up"    activator -> agonise it,    repressor -> block it
        goal="down"  activator -> block it,      repressor -> agonise it

    Downregulation is by far the better-served direction pharmacologically:
    antagonists, inverse agonists and degraders vastly outnumber agonists, and
    they exist for receptors that have no useful agonist at all.
    """
    if goal not in ("up", "down"):
        raise ValueError(f"goal must be 'up' or 'down', got {goal!r}")
    if role in ("activator", "dual"):
        return "up" if goal == "up" else "down"
    if role == "repressor":
        return "down" if goal == "up" else "up"
    return "either"


def filter_drugs_by_goal(drug_map: dict, annotations: dict, goal: str = "up",
                         *, allow_unknown_direction: bool = False
                         ) -> tuple[dict, dict]:
    """Split drugs into (actionable, direction_unknown).

    A drug whose interaction has no recorded direction cannot be placed. The
    whole question here is which way it moves the transcription factor, so
    "unclear" is not weak evidence — it is unusable evidence, and mixing it
    into the recommendations makes the table look authoritative when every row
    says UNKNOWN. DGIdb carries many such interactions: a curated source
    recorded that the drug and gene are related, without saying how.

    They are kept and reported separately, because a known drug-gene
    relationship is still a lead worth chasing in the literature — just not a
    prediction.
    """
    out: dict[str, list[Drug]] = {}
    unknown: dict[str, list[Drug]] = {}

    for symbol, found in drug_map.items():
        role = (annotations.get(symbol) or {}).get("role", "unknown")
        want = desired_direction(role, goal)

        keep, unsure = [], []
        for d in found:
            if d.direction == "unclear":
                (keep if allow_unknown_direction else unsure).append(d)
            elif want == "either" or d.direction == want:
                keep.append(d)

        if keep:
            out[symbol] = sorted(keep, key=lambda d: (
                d.direction != want, d.withdrawn, -(d.max_phase or 0), d.name))
        if unsure:
            unknown[symbol] = sorted(
                unsure, key=lambda d: (-(d.interaction_score or 0), d.name))
    return out, unknown


# ------------------------------------------------------------------- ChEMBL
def _chembl_targets(client: Client, symbol: str, uniprot: str | None) -> list[str]:
    targets: list[dict] = []
    if uniprot:
        try:
            res = client.get(f"{CHEMBL}/target", {
                "target_components__accession": uniprot,
                "format": "json", "limit": 20,
            })
            targets = res.get("targets", [])
        except HttpError:
            pass
    if not targets:
        try:
            res = client.get(f"{CHEMBL}/target/search", {
                "q": symbol, "format": "json", "limit": 20,
            })
            targets = res.get("targets", [])
        except HttpError as exc:
            log.debug("ChEMBL target search failed for %s: %s", symbol, exc)
            return []

    ids = []
    for t in targets:
        if t.get("organism") != "Homo sapiens":
            continue
        if t.get("target_type") not in ("SINGLE PROTEIN", "PROTEIN COMPLEX"):
            continue
        comps = t.get("target_components") or []
        names = {symbol.upper()}
        matched = uniprot is not None and any(c.get("accession") == uniprot for c in comps)
        if not matched:
            for c in comps:
                for syn in c.get("target_component_synonyms") or []:
                    if str(syn.get("component_synonym", "")).upper() in names:
                        matched = True
        if matched:
            ids.append(t["target_chembl_id"])
    return ids


def _chembl_molecules(client: Client, ids: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        try:
            res = client.get(f"{CHEMBL}/molecule", {
                "molecule_chembl_id__in": ",".join(chunk),
                "format": "json", "limit": 50,
            })
        except HttpError:
            continue
        for m in res.get("molecules", []):
            out[m["molecule_chembl_id"]] = m
    return out


def chembl_drugs(client: Client, symbol: str, uniprot: str | None = None,
                 approved_only: bool = True) -> list[Drug]:
    drugs: list[Drug] = []
    for target_id in _chembl_targets(client, symbol, uniprot):
        try:
            res = client.get(f"{CHEMBL}/mechanism", {
                "target_chembl_id": target_id, "format": "json", "limit": 500,
            })
        except HttpError:
            continue
        mechs = res.get("mechanisms", [])
        if not mechs:
            continue
        mols = _chembl_molecules(client, [m["molecule_chembl_id"] for m in mechs
                                          if m.get("molecule_chembl_id")])
        for m in mechs:
            mol = mols.get(m.get("molecule_chembl_id"), {})
            max_phase = mol.get("max_phase")
            max_phase = float(max_phase) if max_phase is not None else None
            approved = max_phase is not None and max_phase >= 4
            if approved_only and not approved:
                continue
            action = (m.get("action_type") or "UNKNOWN").upper()
            if action in ACTIVATING_ACTIONS:
                direction = "up"
            elif action in INHIBITING_ACTIONS:
                direction = "down"
            else:
                direction = "unclear"
            props = mol.get("molecule_properties") or {}
            drugs.append(Drug(
                name=(mol.get("pref_name") or m.get("molecule_chembl_id") or "?").title(),
                chembl_id=m.get("molecule_chembl_id"),
                target_symbol=symbol,
                action_type=action,
                mechanism=m.get("mechanism_of_action"),
                approved=approved,
                max_phase=max_phase,
                first_approval=mol.get("first_approval"),
                withdrawn=bool(mol.get("withdrawn_flag")),
                oral=mol.get("oral"),
                direction=direction,
                source="ChEMBL",
                target_chembl_id=target_id,
                properties=props or None,
            ))
    return drugs


# -------------------------------------------------------------------- DGIdb
_DGIDB_QUERY = """
query($names: [String!]) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        drug { name conceptId approved }
        interactionTypes { type directionality }
        interactionScore
        sources { sourceDbName }
      }
    }
  }
}
"""


def dgidb_drugs(client: Client, symbols: list[str],
                approved_only: bool = True) -> list[Drug]:
    try:
        res = client.post(DGIDB_GRAPHQL, {
            "query": _DGIDB_QUERY,
            "variables": {"names": [s.upper() for s in symbols]},
        })
    except HttpError as exc:
        log.warning("DGIdb query failed: %s", exc)
        return []

    out: list[Drug] = []
    for node in (res.get("data", {}).get("genes", {}) or {}).get("nodes", []) or []:
        sym = node.get("name", "?")
        for inter in node.get("interactions") or []:
            drug = inter.get("drug") or {}
            approved = bool(drug.get("approved"))
            if approved_only and not approved:
                continue
            types = inter.get("interactionTypes") or []
            labels = [str(t.get("type", "")).upper() for t in types]
            dirs = {str(t.get("directionality", "")).upper() for t in types}
            if "ACTIVATING" in dirs or any(x in ACTIVATING_ACTIONS for x in labels):
                direction = "up"
            elif "INHIBITORY" in dirs or any(x in INHIBITING_ACTIONS for x in labels):
                direction = "down"
            else:
                direction = "unclear"
            srcs = sorted({(s0 or {}).get("sourceDbName") or ""
                           for s0 in (inter.get("sources") or []) if s0})
            try:
                iscore = float(inter.get("interactionScore"))
            except (TypeError, ValueError):
                iscore = None
            out.append(Drug(
                name=(drug.get("name") or "?").title(),
                chembl_id=drug.get("conceptId"),
                target_symbol=sym,
                action_type="/".join(labels) or "UNKNOWN",
                mechanism=None,
                approved=approved,
                max_phase=4.0 if approved else None,
                first_approval=None,
                withdrawn=False,
                oral=None,
                direction=direction,
                source="DGIdb",
                interaction_score=iscore,
                sources="; ".join(s for s in srcs if s) or None,
            ))
    return out


# ------------------------------------------------------------------ combine
def drugs_for_tfs(client: Client, annotations: dict, *, approved_only: bool = True,
                  use_dgidb: bool = True, goal: str = "up",
                  allow_unknown_direction: bool = False,
                  use_opentargets: bool = True):
    """Return {tf_symbol: [Drug]} keeping only pharmacologically useful direction."""
    collected: dict[str, list[Drug]] = {}

    dgidb_by_gene: dict[str, list[Drug]] = {}
    if use_dgidb and annotations:
        for d in dgidb_drugs(client, list(annotations), approved_only=approved_only):
            dgidb_by_gene.setdefault(d.target_symbol.upper(), []).append(d)

    for symbol, ann in annotations.items():
        found = chembl_drugs(client, symbol, ann.get("uniprot"), approved_only)
        found += dgidb_by_gene.get(symbol.upper(), [])
        # A third, independent curation. It finds compounds the other two
        # miss, and unlike DGIdb it carries an action type and a clinical
        # stage, which is what the direction filter and grading need.
        if use_opentargets:
            try:
                found += opentargets_drugs(client, symbol,
                                           approved_only=approved_only)
            except Exception as exc:                        # noqa: BLE001
                log.warning("Open Targets drugs unavailable for %s: %s",
                            symbol, str(exc)[:100])

        seen: dict[tuple, Drug] = {}
        for d in found:
            key = (d.name.lower(), d.target_symbol.upper())
            prev = seen.get(key)
            if prev is None:
                seen[key] = d
                continue
            # Same drug from both sources: keep the ChEMBL record (it has the
            # action type and approval year) but carry over DGIdb's score and
            # provenance instead of dropping them.
            if prev.source == "DGIdb" and d.source == "ChEMBL":
                d.interaction_score = d.interaction_score or prev.interaction_score
                d.sources = d.sources or prev.sources
                seen[key] = d
            elif prev.source == "ChEMBL" and d.source == "DGIdb":
                prev.interaction_score = prev.interaction_score or d.interaction_score
                prev.sources = prev.sources or d.sources
        if seen:
            collected[symbol] = list(seen.values())

    return filter_drugs_by_goal(collected, annotations, goal,
                                allow_unknown_direction=allow_unknown_direction)

# ====================================================================== #
# Open Targets as a third drug source.                                   #
#                                                                        #
# This API was already being called for genetics and tissue expression   #
# and never asked for drugs, which was a free source left on the table.  #
# It is a separate curation from ChEMBL and DGIdb, so it finds different #
# compounds: for VDR it returns seocalcitol, becocalcidiol and           #
# eldecalcitol, none of which appeared in a real CCR5 run's VDR list.    #
#                                                                        #
# It also carries the two things the grading needs and DGIdb lacks - an  #
# action type (AGONIST / ANTAGONIST / INHIBITOR) and a clinical stage -  #
# plus drug warnings, which feed the same withdrawn/black-box logic the  #
# rest of the package reads.                                             #
# ====================================================================== #
OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"

# The schema moves. `knownDrugs` was removed before API v26.9 and this is
# its replacement; if it disappears too, the lookup degrades to returning
# nothing rather than taking the run down.
_OT_DRUG_QUERY = """query($id: String!) {
  target(ensemblId: $id) {
    approvedSymbol
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          id name drugType maximumClinicalStage
          drugWarnings { warningType }
          mechanismsOfAction { rows { actionType mechanismOfAction } }
        }
      }
    }
  }
}"""

# Open Targets stages, mapped onto the max_phase the rest of the package uses.
_STAGE_PHASE = {
    "APPROVAL": 4.0, "PHASE_4": 4.0, "PHASE_3": 3.0,
    "PHASE_2": 2.0, "PHASE_1": 1.0, "PRECLINICAL": 0.0,
}


def opentargets_drugs(client: Client, symbol: str, *,
                      approved_only: bool = True) -> list[Drug]:
    """Drugs and clinical candidates for one target, from Open Targets."""
    from .genetics import _resolve_ensembl

    try:
        ensembl = _resolve_ensembl(client, symbol)
    except Exception:                                       # noqa: BLE001
        ensembl = None
    if not ensembl:
        return []

    try:
        res = client.post(OPENTARGETS, {"query": _OT_DRUG_QUERY,
                                        "variables": {"id": ensembl}})
    except Exception as exc:                                # noqa: BLE001
        log.warning("Open Targets drug lookup failed for %s: %s",
                    symbol, str(exc)[:120])
        return []
    if res.get("errors"):
        log.warning("Open Targets drug schema rejected the query for %s (%s) "
                    "- this source is skipped, the rest of the run is "
                    "unaffected", symbol, str(res["errors"])[:120])
        return []

    target = (res.get("data") or {}).get("target") or {}
    rows = ((target.get("drugAndClinicalCandidates") or {}).get("rows") or [])

    out: list[Drug] = []
    for row in rows:
        d = row.get("drug") or {}
        name = (d.get("name") or "").strip()
        if not name:
            continue
        stage = d.get("maximumClinicalStage") or row.get("maxClinicalStage")
        phase = _STAGE_PHASE.get(str(stage).upper())
        approved = phase == 4.0
        if approved_only and not approved:
            continue

        moa_rows = ((d.get("mechanismsOfAction") or {}).get("rows") or [])
        moa = moa_rows[0] if moa_rows else {}
        action = (moa.get("actionType") or "").strip().lower()

        # Open Targets records a withdrawal as a drug warning; the rest of
        # the package reads that flag, so it has to arrive the same way.
        warnings = d.get("drugWarnings") or []
        withdrawn = any("withdraw" in str(w.get("warningType", "")).lower()
                        for w in warnings)

        out.append(Drug(
            name=name.title() if name.isupper() else name,
            chembl_id=d.get("id"),
            target_symbol=symbol.upper(),
            action_type=action,
            mechanism=moa.get("mechanismOfAction"),
            approved=bool(approved),
            max_phase=phase,
            first_approval=None,
            withdrawn=withdrawn,
            oral=None,
            direction=_direction_for_action(action),
            source="OpenTargets",
        ))

    log.info("Open Targets: %d drug(s) for %s", len(out), symbol)
    return out


def _direction_for_action(action: str) -> str:
    """Map an Open Targets action type onto the direction the FILTER reads.

    The vocabulary matters more than it looks. `drugs_for_tfs` compares
    `d.direction == want` where want is "up", "down" or "either", and routes
    "unclear" to the unknown-direction bucket. The Drug dataclass comment
    says "increases TF activity", which is stale - following it meant every
    Open Targets drug failed the comparison and was dropped without even
    landing in the unknown bucket, for any factor with a role annotation.
    """
    a = (action or "").lower()
    # The DECREASE patterns are tested first, and that order is load-bearing:
    # "antagonist" contains "agonist" as a substring, so checking the
    # increase patterns first mapped every antagonist to "increases TF
    # activity" - the exact opposite of the truth, which would have put
    # blockers in the raise-the-gene list.
    if any(k in a for k in ("antagonist", "inverse agonist", "inhibitor",
                            "blocker", "degrader", "negative", "suppressor",
                            "antisense", "disrupt")):
        return "down"
    if any(k in a for k in ("agonist", "activator", "opener", "inducer",
                            "positive", "stabiliser", "stabilizer")):
        return "up"
    return "unclear"
