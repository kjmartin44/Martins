"""Does the drug reach the cells that express your gene?

Three things must overlap, and the pipeline otherwise checks none of them:

    where the drug goes   ∩   where the gene is expressed   ∩   where the TF is

Any empty intersection kills the route regardless of how good the motif, the
potency or the grade is. The third term is the one most often skipped: a
beautifully scored hepatocyte transcription factor is irrelevant if you are
treating neurons, and nothing upstream in this package notices.

CNS penetration: CNS MPO, not a rule of thumb
---------------------------------------------
The earlier MW/TPSA/logP check here was a rule of thumb. The published standard
is the CNS Multiparameter Optimisation score (Wager et al., ACS Chem Neurosci
2010): six physicochemical properties, each mapped through a desirability
function to 0-1 and summed to 0-6. A score >= 4 is the conventional threshold
for favourable CNS exposure, and it discriminates marketed CNS drugs from
peripheral ones considerably better than any single cutoff.

  ClogP     <= 3 ideal, >= 5 unfavourable
  ClogD7.4  <= 2 ideal, >= 4 unfavourable
  MW        <= 360 ideal, >= 500 unfavourable
  TPSA      40-90 ideal; both too polar and too apolar are penalised
  HBD       <= 0.5 ideal, >= 3.5 unfavourable
  pKa       <= 8 ideal, >= 10 unfavourable

ClogD and pKa need tooling beyond RDKit, so when they are unavailable the score
is computed over the four remaining parameters and **reported as such** rather
than silently rescaled — a 4-parameter MPO is not comparable to a 6-parameter
one and the output says which you have.

A worked limitation: atenolol scores 3.17/4 here, which reads as favourable,
yet it is a textbook example of poor brain penetration. The missing pKa term
(~9.6, strongly basic) is precisely what would have penalised it. Supply `pka`
and `logd` through `properties` when you have them; without them, treat a
borderline pass on a basic amine with suspicion.

The bigger reason CNS drugs fail
--------------------------------
Physicochemistry predicts passive diffusion. The commonest reason an otherwise
suitable molecule fails to reach brain is **active efflux by P-glycoprotein**
(ABCB1) at the endothelium. A compound can have an excellent MPO and still be
pumped straight back out. P-gp substrate prediction is genuinely poor, so this
module uses a curated list of known substrates plus a weak structural flag, and
says plainly that a negative result is not reassurance.

Other compartments
------------------
Be honest about the limits: **BBB is the only barrier with a decent predictive
model.** For liver, kidney, muscle or tumour, distribution depends on
transporters, protein binding and perfusion in ways physicochemistry does not
capture. Where a compartment is requested and no model exists, this module says
so instead of inventing a number.

Modality vs target compartment
------------------------------
A separate, absolute constraint that is easy to miss: **antibodies and other
biologics cannot reach intracellular targets.** Every transcription factor is
intracellular. An antibody against a TF is not a weak option, it is not an
option. Checked explicitly.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, asdict, field

log = logging.getLogger(__name__)

OPENTARGETS = "https://api.platform.opentargets.org/api/v4/graphql"

# Known or well-documented P-gp (ABCB1) substrates, lower-cased.
# An MPO needs at least this many of its six parameters to mean anything.
MIN_MPO_PARAMETERS = 3

PGP_SUBSTRATES = {
    "digoxin", "loperamide", "domperidone", "fexofenadine", "cyclosporine",
    "tacrolimus", "sirolimus", "everolimus", "paclitaxel", "docetaxel",
    "doxorubicin", "daunorubicin", "vinblastine", "vincristine", "etoposide",
    "imatinib", "dasatinib", "nilotinib", "lapatinib", "erlotinib", "gefitinib",
    "sunitinib", "sorafenib", "colchicine", "quinidine", "verapamil",
    "diltiazem", "rifampicin", "ritonavir", "saquinavir", "indinavir",
    "nelfinavir", "dabigatran", "rivaroxaban", "apixaban", "edoxaban",
    "ondansetron", "risperidone", "paliperidone", "amitriptyline",
    "nortriptyline", "morphine", "methadone", "phenytoin", "carbamazepine",
    "dexamethasone", "prednisolone", "hydrocortisone", "aldosterone",
}

BIOLOGIC_MODALITIES = {"monoclonal antibody", "antibody", "protein",
                       "enzyme replacement", "peptide"}


@dataclass
class DeliveryAssessment:
    drug: str
    compartment: str
    cns_mpo: float | None = None
    mpo_parameters: int = 0
    mpo_components: dict = field(default_factory=dict)
    likely_penetrant: bool | None = None
    pgp_substrate: bool | None = None
    compatible: bool = True
    multiplier: float = 1.0
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------- CNS MPO math
def _ramp(x: float, good: float, bad: float) -> float:
    """1.0 at or beyond `good`, 0.0 at or beyond `bad`, linear between."""
    if good < bad:
        if x <= good:
            return 1.0
        if x >= bad:
            return 0.0
        return (bad - x) / (bad - good)
    if x >= good:
        return 1.0
    if x <= bad:
        return 0.0
    return (x - bad) / (good - bad)


def _hump(x: float, lo_zero: float, lo_one: float,
          hi_one: float, hi_zero: float) -> float:
    """Rises to 1 across [lo_zero, lo_one], flat, falls to 0 across
    [hi_one, hi_zero]. TPSA is penalised at both ends."""
    if x <= lo_zero or x >= hi_zero:
        return 0.0
    if lo_one <= x <= hi_one:
        return 1.0
    if x < lo_one:
        return (x - lo_zero) / (lo_one - lo_zero)
    return (hi_zero - x) / (hi_zero - hi_one)


def cns_mpo(smiles: str | None = None, *, properties: dict | None = None
            ) -> tuple[float | None, int, dict]:
    """(score, n_parameters, components). Score is out of n_parameters."""
    props = dict(properties or {})

    if smiles:
        try:
            from rdkit import Chem, RDLogger
            from rdkit.Chem import Descriptors
            RDLogger.DisableLog("rdApp.*")
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                props.setdefault("mw", Descriptors.MolWt(mol))
                props.setdefault("logp", Descriptors.MolLogP(mol))
                props.setdefault("tpsa", Descriptors.TPSA(mol))
                props.setdefault("hbd", Descriptors.NumHDonors(mol))
        except ImportError:
            pass

    comp: dict[str, float] = {}
    if props.get("logp") is not None:
        comp["clogp"] = round(_ramp(float(props["logp"]), 3.0, 5.0), 3)
    if props.get("logd") is not None:
        comp["clogd"] = round(_ramp(float(props["logd"]), 2.0, 4.0), 3)
    if props.get("mw") is not None:
        comp["mw"] = round(_ramp(float(props["mw"]), 360.0, 500.0), 3)
    if props.get("tpsa") is not None:
        comp["tpsa"] = round(_hump(float(props["tpsa"]), 20, 40, 90, 120), 3)
    if props.get("hbd") is not None:
        comp["hbd"] = round(_ramp(float(props["hbd"]), 0.5, 3.5), 3)
    if props.get("pka") is not None:
        comp["pka"] = round(_ramp(float(props["pka"]), 8.0, 10.0), 3)

    if not comp:
        return None, 0, {}
    return round(sum(comp.values()), 2), len(comp), comp


def is_pgp_substrate(drug_name: str, properties: dict | None = None
                     ) -> tuple[bool | None, str | None]:
    """Known-substrate lookup first; structural guess only as a weak hint."""
    if drug_name and drug_name.strip().lower() in PGP_SUBSTRATES:
        return True, (f"{drug_name} is a known P-glycoprotein substrate: it will "
                      "be actively pumped out of brain endothelium regardless of "
                      "favourable physicochemistry")
    p = properties or {}
    try:
        mw = float(p.get("mw") or 0)
        tpsa = float(p.get("tpsa") or 0)
    except (TypeError, ValueError):
        return None, None
    if mw > 400 and tpsa > 75:
        return None, ("structurally consistent with P-gp efflux (MW>400, "
                      "TPSA>75), though substrate prediction is unreliable — "
                      "check experimentally before assuming brain exposure")
    return None, None


# ------------------------------------------------------------- the assessment
def assess_delivery(drug_name: str, *, compartment: str = "unknown",
                    smiles: str | None = None, properties: dict | None = None,
                    modality: str = "small molecule",
                    target_is_intracellular: bool = True) -> DeliveryAssessment:
    a = DeliveryAssessment(drug=drug_name, compartment=compartment)

    # --- absolute constraint: biologics cannot reach intracellular targets ---
    if modality.lower() in BIOLOGIC_MODALITIES and target_is_intracellular:
        a.compatible = False
        a.multiplier = 0.0
        a.notes.append(
            f"{modality} cannot reach an intracellular target. Every "
            "transcription factor is intracellular, so this is not a weak "
            "option — it is not an option.")
        return a

    if compartment != "cns":
        if compartment not in ("unknown", ""):
            a.notes.append(
                f"no validated physicochemical model exists for {compartment} "
                "exposure; distribution there depends on transporters, protein "
                "binding and perfusion. Use measured PK rather than a prediction.")
        return a

    score, n_params, comp = cns_mpo(smiles, properties=properties)
    a.cns_mpo, a.mpo_parameters, a.mpo_components = score, n_params, comp

    # Fewer than three parameters is not a multiparameter score. Reporting
    # "MPO 0.5/1" invites it to be read as a real assessment when it rests on
    # a single descriptor.
    if score is None or n_params < MIN_MPO_PARAMETERS:
        a.cns_mpo, a.mpo_components = None, {}
        have = ", ".join(comp) if comp else "none"
        a.notes.append(
            f"CNS exposure not assessed: only {n_params} of 6 physicochemical "
            f"parameters available ({have}). A score on fewer than "
            f"{MIN_MPO_PARAMETERS} is not meaningful — supply mw, tpsa, logp "
            f"(and ideally logd, pka) to get one.")
        a.multiplier = 0.7
        a.likely_penetrant = None
        return a

    # Threshold scales with how many parameters were available.
    threshold = 4.0 * (n_params / 6.0)
    a.likely_penetrant = score >= threshold
    a.multiplier = round(min(1.0, 0.25 + 0.75 * (score / max(n_params, 1))), 3)

    if n_params < 6:
        a.notes.append(
            f"CNS MPO computed over {n_params} of 6 parameters "
            f"(ClogD and pKa unavailable). Scored {score}/{n_params} against a "
            f"{threshold:.1f} threshold — not directly comparable to a published "
            "6-parameter MPO.")
    if not a.likely_penetrant:
        worst = sorted(comp.items(), key=lambda kv: kv[1])[:2]
        a.notes.append("unlikely to reach CNS at useful exposure; weakest terms: "
                       + ", ".join(f"{k} ({v})" for k, v in worst))

    pgp, note = is_pgp_substrate(drug_name, properties or
                                 {"mw": (comp and properties or {}).get("mw")})
    a.pgp_substrate = pgp
    if note:
        a.notes.append(note)
    if pgp:
        a.multiplier *= 0.3
        a.likely_penetrant = False
    elif a.likely_penetrant:
        a.notes.append("physicochemistry favours passive CNS entry, but this "
                       "does not exclude active efflux — a negative P-gp flag "
                       "here is not reassurance")
    return a


# --------------------------------------------------- tissue expression overlap
# Open Targets changes this schema between releases. Try the richest form
# first and fall back; a 400 means the shape is wrong, not that the gene is
# missing, so there is no point retrying the same query for every TF.
_EXPR_QUERIES = [
    """query($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    expressions { tissue { label } rna { level } protein { level } }
  }
}""",
    """query($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    expressions { tissue { label } rna { level } }
  }
}""",
]

# Set once the endpoint has proved unusable, so we warn once per run.
_EXPR_DISABLED = {"value": False, "reason": None}

# A locally-loaded expression table takes priority over the API entirely.
_LOCAL_EXPRESSION: dict[str, dict[str, float]] = {}


def expression_source_status() -> dict:
    """Is tissue-expression data actually available for this run?

    `tissue_expression` returns an empty set both when a gene is genuinely
    tissue-restricted and when the Open Targets schema changed and the lookup
    disabled itself. `expression_overlap` then reports "unknown" for every
    factor, which reads like a data gap in HPA rather than a broken endpoint.
    This says which it was.
    """
    if _LOCAL_EXPRESSION:
        return {"requested": True, "ok": True, "yielded": len(_LOCAL_EXPRESSION),
                "detail": (f"local expression table: {len(_LOCAL_EXPRESSION)} "
                           f"genes")}
    if _EXPR_DISABLED["value"]:
        return {"requested": True, "ok": False, "yielded": 0,
                "detail": (f"Open Targets expression lookup failed and was "
                           f"disabled for this run ({_EXPR_DISABLED['reason']}). "
                           f"Pass --expression-file with an HPA "
                           f"rna_tissue_consensus.tsv or a GTEx median-TPM "
                           f"matrix to remove the dependency entirely.")}
    return {"requested": True, "ok": True, "yielded": 0,
            "detail": "Open Targets expression lookup in use"}


def load_expression_file(path, *, min_value: float = 1.0) -> dict:
    """Load a tissue-expression table from disk, so this works with no network.

    Accepts the two formats people actually have:

      long  (Human Protein Atlas `rna_tissue_consensus.tsv`)
            columns: Gene, Gene name, Tissue, nTPM
      wide  (GTEx median TPM, or any gene-by-tissue matrix)
            first column gene symbol, remaining columns one tissue each

    HPA:  https://www.proteinatlas.org/about/download  (file 'RNA consensus
          tissue gene data')
    GTEx: https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
          ('gene median TPM')

    `min_value` is in the table's own units (nTPM for HPA, TPM for GTEx);
    1.0 is a conventional floor for "expressed".
    """
    import csv
    import gzip
    import io
    import zipfile
    from pathlib import Path

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found")

    # HPA ships a .zip and GTEx a .gz, and both are routinely saved under a
    # .tsv name. Sniff the magic bytes rather than trusting the extension.
    with open(path, "rb") as probe:
        magic = probe.read(4)

    if magic[:2] == b"PK":
        def opener():
            zf = zipfile.ZipFile(path)
            names = [n for n in zf.namelist()
                     if n.lower().endswith((".tsv", ".csv", ".txt"))]
            if not names:
                raise ValueError(
                    f"{path.name} is a zip archive with no .tsv/.csv inside "
                    f"(contains: {', '.join(zf.namelist()[:5])})")
            log.info("reading %s from inside %s", names[0], path.name)
            return io.TextIOWrapper(zf.open(names[0]), encoding="utf-8",
                                    errors="replace")
    elif magic[:2] == b"\x1f\x8b":
        def opener():
            return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8",
                                    errors="replace")
    elif magic[:1] == b"<":
        raise ValueError(
            f"{path.name} is an HTML page, not data — the download failed. "
            f"Fetch it again from https://www.proteinatlas.org/about/download")
    else:
        def opener():
            return open(path, encoding="utf-8", errors="replace")

    with opener() as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader, None)
        if not header:
            raise ValueError(f"{path.name} is empty")

        lower = [h.strip().strip('"').lower() for h in header]
        long_format = any(h in ("tissue", "tissue name") for h in lower)

        table: dict[str, dict[str, float]] = {}
        def _col(*names, default=None):
            """Find a column by NAME priority, not position.

            HPA's file has both 'Gene' (an Ensembl id) and 'Gene name' (the
            symbol). Scanning positionally grabs the id, which then matches
            nothing downstream.
            """
            for name in names:
                if name in lower:
                    return lower.index(name)
            return default

        if long_format:
            gi = _col("gene name", "gene_name", "symbol", "description",
                      "gene", default=0)
            ti = _col("tissue", "tissue name", default=1)
            vi = _col("ntpm", "tpm", "nx", "value", "expression",
                      default=len(lower) - 1)
            for row in reader:
                if len(row) <= max(gi, ti, vi):
                    continue
                try:
                    v = float(row[vi])
                except ValueError:
                    continue
                if v < min_value:
                    continue
                table.setdefault(row[gi].strip().upper(), {})[
                    row[ti].strip().lower()] = v
        else:
            # Wide: skip an Ensembl-id first column if a symbol column follows.
            gi = (1 if lower[0] in ("name", "gene_id", "ensembl", "id")
                  and len(lower) > 1 else 0)
            tissues = [h.strip().strip('"').lower()
                       for h in header[gi + 1:]]
            for row in reader:
                if len(row) <= gi + 1:
                    continue
                gene = row[gi].strip().upper()
                d = {}
                for t, cell in zip(tissues, row[gi + 1:]):
                    try:
                        v = float(cell)
                    except ValueError:
                        continue
                    if v >= min_value:
                        d[t] = v
                if d:
                    table[gene] = d

    # The wide-format branch accepts "a gene column plus numeric columns",
    # which any numeric matrix satisfies. Without this check a wrong file
    # loads as nonsense instead of erroring, which is worse than failing.
    if table:
        looks_like_genes = sum(1 for g in table if g[:1].isalpha())
        if looks_like_genes < 0.5 * len(table):
            raise ValueError(
                f"{path.name} does not look like an expression table: most row "
                f"labels are not gene symbols (saw {list(table)[:5]}). Expected "
                f"HPA rna_tissue_consensus.tsv or a GTEx median-TPM matrix.")
        if len(table) < 50:
            log.warning("only %d genes parsed from %s — if that seems low, "
                        "check this is the whole file and not a fragment",
                        len(table), path.name)

    if not table:
        raise ValueError(
            f"no expression rows parsed from {path.name}. Expected either a "
            f"long table with Tissue and nTPM/TPM columns (HPA "
            f"rna_tissue_consensus) or a wide gene-by-tissue matrix (GTEx "
            f"median TPM). Header seen: {', '.join(header[:6])}")
    log.info("loaded expression for %d genes from %s", len(table), path.name)
    _LOCAL_EXPRESSION.clear()
    _LOCAL_EXPRESSION.update(table)
    return table


def tissue_expression(client, symbol: str, *, min_level: int = 2) -> set:
    """Tissues where this gene is expressed above a threshold (HPA levels 0-3)."""
    from .genetics import _resolve_ensembl

    # A local table beats the API: faster, reproducible, and works offline.
    if _LOCAL_EXPRESSION:
        return set(_LOCAL_EXPRESSION.get(symbol.upper(), {}))

    if _EXPR_DISABLED["value"]:
        return set()

    ensembl = _resolve_ensembl(client, symbol)
    if not ensembl:
        return set()

    res, last = None, None
    for query in _EXPR_QUERIES:
        try:
            res = client.post(OPENTARGETS, {"query": query,
                                            "variables": {"ensemblId": ensembl}})
            if not (res.get("errors")):
                break
            last, res = res["errors"], None
        except Exception as exc:                          # noqa: BLE001
            last, res = exc, None

    if res is None:
        # A schema error will fail identically for every gene; warn once and
        # stop asking rather than emitting one line per transcription factor.
        _EXPR_DISABLED["value"] = True
        _EXPR_DISABLED["reason"] = str(last)[:160]
        log.warning(
            "Open Targets expression lookup unavailable (%s). Tissue "
            "co-expression checking is DISABLED for this run — TF/gene "
            "overlap will not be filtered. Everything else is unaffected.",
            _EXPR_DISABLED["reason"])
        return set()

    target = (res.get("data") or {}).get("target") or {}
    out = set()
    for e in target.get("expressions") or []:
        tissue = ((e.get("tissue") or {}).get("label") or "").lower()
        rna = (e.get("rna") or {}).get("level")
        prot = (e.get("protein") or {}).get("level")
        level = max([x for x in (rna, prot) if x is not None] or [0])
        if tissue and level >= min_level:
            out.add(tissue)
    return out


def expression_overlap(gene_tissues: set, tf_tissues: set) -> dict:
    """Do the gene and its regulator occur in the same tissues?

    An empty intersection means the predicted regulation cannot be happening
    where you care about it, whatever the motif score says.
    """
    if not gene_tissues or not tf_tissues:
        return {"known": False, "shared": [], "note":
                "tissue expression unavailable for one or both; cannot check "
                "whether the factor is present where the gene is expressed"}
    shared = sorted(gene_tissues & tf_tissues)
    result = {"known": True, "shared": shared[:15], "n_shared": len(shared),
              "n_gene": len(gene_tissues), "n_tf": len(tf_tissues)}
    if not shared:
        result["note"] = (
            "NO shared tissue: this factor and your gene are not expressed in "
            "the same places, so the predicted regulation cannot operate where "
            "you need it. Motif quality is irrelevant here.")
        result["compatible"] = False
    else:
        result["compatible"] = True
        result["note"] = f"co-expressed in {len(shared)} tissue(s)"
    return result
