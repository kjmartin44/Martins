"""Human peak plasma concentrations read from FDA drug labels.

Why this exists
---------------
The dose check compares the concentration a person reaches with the one the
evidence was seen at. ChEMBL is the only source of the first, and its Cmax
records are mostly mouse and rat studies from medicinal chemistry papers: of 26
FDA approved candidates for one gene, 13 had no ChEMBL Cmax record at all. The
label of an approved drug states a human value, in people dosed as prescribed.

What it does
------------
For a drug name it asks openFDA for the single ingredient labels, takes the
pharmacokinetics section (12.3), and looks for sentences that state a peak
concentration: "Cmax", "peak plasma concentration", "maximum concentration".
In those sentences it reads each number that has a concentration unit, with
the dose, whether the figure is single dose or steady state, and the sentence
it came from.

What it refuses to read, and why
--------------------------------
A label states many peaks, and the one a reader wants is the one for the
approved dose in an adult. A value is not used when its sentence or heading is
about a special population (renal or hepatic impairment, elderly, children, sex,
race), a drug interaction, food, or a metabolite. The labels are flattened to
running text with the table structure lost, so these checks look at the heading
the text sits under and at the sentence itself. This is what the earlier
openFDA attempt did not do, and a Cmax from an interaction study read as a blood
level is how calcitriol came out thousands of times too high.

A drug is written to the table only when its accepted values agree within
thirty fold (the rule the ChEMBL reader uses). Anything else goes to a review
file with the sentence, so a person can read it. Nothing is guessed: a drug
whose label has no peak sentence, no usable unit, or no molecular weight for a
mass unit gets no row.

What it cannot promise
----------------------
It is a text reader, not a pharmacologist. A peak for a dose that is not the
usual dose, a figure for a single formulation, or a sentence phrased in a way
the patterns do not know can all get through or be missed. Every row carries
the label, the dose and the sentence so it can be checked. The value is a total
plasma concentration unless the label gives a protein binding, in which case
the unbound fraction is written beside it.
"""

from __future__ import annotations

import logging
import re
import statistics
from collections import OrderedDict

log = logging.getLogger(__name__)

OPENFDA_LABEL = "https://api.fda.gov/drug/label.json"
MAX_SPREAD = 30.0               # the same rule exposure applies to ChEMBL
MAX_LABELS = 25

_PEAK = re.compile(
    r"\b(?:c\s?max|peak\s+(?:[\w-]+\s+){0,3}?(?:concentrations?|levels?)|"
    r"(?:maximum|maximal|highest)\s+(?:[\w-]+\s+){0,3}?"
    r"(?:concentrations?|levels?))", re.I)

# Headings the flattened text still carries. A segment under one of the
# first group is about a population or a condition of use and is skipped.
_SKIP_HEADINGS = (
    "Specific Populations", "Special Populations", "Renal Impairment",
    "Hepatic Impairment", "Pediatric", "Paediatric", "Geriatric", "Elderly",
    "Gender", "Race", "Drug Interaction", "Drug-Drug", "Interactions",
    "Effect of Food", "Food Effect", "Food Effects", "Pharmacogenomics",
    "Clinical Studies", "In Vitro Studies", "Cytochrome P450", "Pregnancy",
    "Lactation", "Body Weight", "Smoking")
_ALL_HEADINGS = _SKIP_HEADINGS + (
    "Absorption", "Distribution", "Metabolism", "Elimination", "Excretion",
    "Pharmacodynamics", "Mechanism of Action", "Steady State", "Bioavailability")
_HEADING = re.compile(r"(?<![\w-])(" + "|".join(
    re.escape(h) for h in sorted(_ALL_HEADINGS, key=len, reverse=True))
    + r")(?=\s+[A-Z(])")

_SKIP_WORDS = re.compile(
    r"\b(?:impair\w*|renal|hepatic|cirrhosis|elderly|geriatric|pediatric|"
    r"paediatric|children|child|adolescents?|infants?|neonat\w+|"
    r"females?|males?|women|men|sex|race|japanese|chinese|asian|"
    r"interaction\w*|co-?administ\w+|concomitant\w*|ketoconazole|ritonavir|"
    r"itraconazole|fluconazole|rifampin|rifampicin|cimetidine|"
    r"high[- ]fat|fed|fasted\s+state|with\s+food|food|meal|"
    r"metabolites?|ratio|geometric\s+mean\s+ratio|bioequivalen\w+|"
    r"dialysis|hemodialysis|obes\w+|crushed|sprinkled|suspension)\b", re.I)

_NUM = r"(?:\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?)"
# The units a concentration is stated in. The micro sign is lost to a
# replacement character in some labels, and the Greek mu is not the micro sign.
_MU = "[µμu�]"
_UNIT = (r"(?:(?:ng|mcg|pg|mg|" + _MU + r"g)\s?/\s?(?:mL|ml|L)\b"
         r"|(?:n|" + _MU + r"|m)(?:M\b|mol\s?/\s?L\b)|(?:n|" + _MU
         + r")mol\s?/\s?mL\b)")
_SD = r"(?:\s*(?:±|�|\+/-|\+-)\s*(?:" + _NUM + r"))?"
_PAREN = r"(?:\s*\([^)]{0,30}\))?"
_VALUE = re.compile(r"(?P<v>" + _NUM + r")" + _SD + _PAREN + r"\s*(?P<u>"
                    + _UNIT + r")", re.I)
_BARE = re.compile(r"(?P<v>" + _NUM + r")(?:\s*(?:±|�|\+/-|\+-)\s*"
                   + _NUM + r")(?!\s*(?:" + _UNIT + r"))", re.I)
_UNIT_FIRST = re.compile(
    r"(?:c\s?max|peak[^()]{0,30})\s*\(\s*(?P<u>" + _UNIT + r")\s*\)\s*[:|=-]?\s*"
    r"(?P<v>" + _NUM + r")", re.I)
_NOT_A_VALUE = re.compile(
    r"^\s*(?:%|mg|mcg|g\b|kg|h\b|hr|hours?|min|minutes?|days?|weeks?|months?|"
    r"years?|m\s?2|times|-fold|fold)", re.I)
# A dose is a mass not followed by "h/mL" (an AUC unit, as in "2.9 mcg h/mL") or
# "/mL". The number is a whole token, so "2,000 mg" is not read as "000 mg".
_DOSE = re.compile(r"(?<![\d,.])(" + _NUM + r")\s*(mg|mcg|g|" + _MU
                   + r"g)\b(?!\s?/\s?m[lL])(?!\s*[·.*]?\s*h\s*/)", re.I)
_DOSE_LIST = re.compile(r"(?:for|at|of)\s+((?:" + _NUM + r")(?:\s*,\s*(?:"
                        + _NUM + r"))*(?:,?\s*and\s+(?:" + _NUM
                        + r"))?)\s*(mg|mcg|g)\b", re.I)
# "0.6, 1.1, 1.4 and": figures that share the unit written after the last one.
_LEAD_LIST = re.compile(r"((?:" + _NUM + r")(?:\s*,\s*(?:" + _NUM
                        + r"))*,?\s*(?:and\s+)?)$", re.I)
_RANGE_END = re.compile(r"(?:range[^()]{0,25}?\bto|range[^()]{0,25}?[-–])\s*$",
                        re.I)


def clean(text: str) -> str:
    """Label text with the replacement characters made readable."""
    return re.sub(r"\s+", " ", text)


def _unit(raw: str) -> str:
    """The unit spelled the way `exposure.to_nm` reads it."""
    u = re.sub(r"\s+", "", raw)
    u = re.sub(r"^[µμ�]", "u", u)
    u = u.replace("mcg", "ug")
    low = u.lower()
    molar = {"um": "uM", "nm": "nM", "mm": "mM", "nmol/l": "nM",
             "umol/l": "uM", "mmol/l": "mM", "nmol/ml": "uM"}
    return molar.get(low, u)


def _to_float(token: str) -> float:
    return float(token.replace(",", ""))


def segments(text: str) -> list:
    """`[(heading, text)]`: the label text cut at the headings it still has."""
    parts, last, head = [], 0, ""
    for m in _HEADING.finditer(text):
        if m.start() > last:
            parts.append((head, text[last:m.start()]))
        head, last = m.group(1), m.end()
    parts.append((head, text[last:]))
    return [(h, t.strip()) for h, t in parts if t.strip()]


def _sentences(text: str) -> list:
    return [s for s in re.split(r"(?<=[.;])\s+(?=[A-Z0-9(])", text) if s]


def candidates(text: str) -> tuple:
    """`(accepted, rejected)` peak concentration statements in a label section.

    Each is `{"value", "unit", "dose", "kind", "sentence", "heading",
    "inherited"}`; a rejected one also has `why`.
    """
    accepted, rejected = [], []
    for heading, seg in segments(clean(text)):
        skipped_heading = any(heading.lower() == h.lower()
                              for h in _SKIP_HEADINGS)
        for sent in _sentences(seg):
            if not _PEAK.search(sent):
                continue
            found = _values(sent)
            if not found:
                continue
            why = None
            if skipped_heading:
                why = f"under the heading {heading!r}"
            else:
                m = _SKIP_WORDS.search(sent)
                if m:
                    why = f"sentence mentions {m.group(0)!r}"
            dose = _DOSE.search(sent)
            kind = ("steady state" if re.search(
                r"steady[- ]state|multiple[- ]dos|after\s+\d+\s+(?:days?|weeks?)|"
                r"repeated|continuous", sent, re.I)
                else "single dose" if re.search(
                    r"single|first dose|initial dose|one dose", sent, re.I)
                else "unspecified")
            dl = _DOSE_LIST.search(sent)
            doses = re.findall(_NUM, dl.group(1)) if dl else []
            listed = len(doses) == len(found) and len(found) > 1
            for k, (value, unit, inherited) in enumerate(found):
                one = (f"{doses[k]} {dl.group(2)}" if listed else
                       f"{dose.group(1)} {dose.group(2)}"
                       if dose and len(found) == 1 else None)
                rec = {"value": value, "unit": unit,
                       "dose": one,
                       "kind": kind, "sentence": sent[:300],
                       "heading": heading, "inherited": inherited}
                if why:
                    rejected.append(dict(rec, why=why))
                else:
                    accepted.append(rec)
    return accepted, rejected


def _values(sentence: str) -> list:
    """`[(value, unit, inherited)]` for each concentration in a sentence."""
    out, spans = [], []
    for m in _VALUE.finditer(sentence):
        if _NOT_A_VALUE.match(sentence[m.end():m.end() + 12]):
            continue
        before = sentence[:m.start("v")]
        if _RANGE_END.search(before):
            continue                    # the top of "(range 35 to 45 ng/mL)"
        unit = _unit(m.group("u"))
        lead = _LEAD_LIST.search(before)
        earlier = ([_to_float(t) for t in re.findall(_NUM, lead.group(1))]
                   if lead and ("," in lead.group(1)
                                or "and" in lead.group(1)) else [])
        for v in earlier:
            out.append((v, unit, True))
        out.append((_to_float(m.group("v")), unit, False))
        spans.append((m.start(), m.end()))
    m = _UNIT_FIRST.search(sentence)
    if m and not out:
        out.append((_to_float(m.group("v")), _unit(m.group("u")), False))
    units = {u for _, u, _ in out}
    if len(units) == 1:
        # "394 +/- 89 after the first dose and 138 +/- 139 ng/mL after a week":
        # the first figure has no unit of its own and takes the later one's.
        first_end = min(s for s, _ in spans) if spans else len(sentence)
        for b in _BARE.finditer(sentence):
            if b.start() < first_end and not _NOT_A_VALUE.match(
                    sentence[b.end():b.end() + 12]):
                out.insert(0, (_to_float(b.group("v")), next(iter(units)), True))
    return out


_BOUND = re.compile(
    r"(?:(?P<q>approximately|about|>|greater than|more than|over|~)\s*)?"
    r"(?P<p>\d+(?:\.\d+)?)\s*%\s*(?:of\s+(?:the\s+)?(?:drug|\w+)\s+)?"
    r"(?:is\s+|are\s+|was\s+)?(?:protein[- ]|plasma[- ]protein[- ])?bound", re.I)
_BINDING = re.compile(
    r"protein[- ]binding[^.]{0,60}?(?:(?P<q>approximately|about|>|greater than|"
    r"more than|over|~)\s*)?(?P<p>\d+(?:\.\d+)?)\s*%", re.I)


def unbound_fraction(text: str) -> tuple:
    """`(fraction or None, note)` from a protein binding statement.

    Used only when every statement agrees on one percentage, and only for
    statements that are not about a metabolite. A ">" or "greater than" gives
    the limit, so the fraction is an upper bound and the note says so.
    """
    found = set()
    for sent in _sentences(clean(text)):
        if re.search(r"metabolite|in vitro studies? (?:with|using) human "
                     r"serum", sent, re.I):
            continue
        for rx in (_BOUND, _BINDING):
            for m in rx.finditer(sent):
                # "6% to 12% bound" is a range, and either end alone would be
                # a claim the label did not make. Count it as a disagreement.
                if re.search(r"\d\s*%?\s*(?:to|-|–|and)\s*$",
                             sent[:m.start("p")], re.I) or re.match(
                        r"\s*%?\s*(?:to|-|–)\s*\d", sent[m.end("p"):], re.I):
                    found.add((-1.0, False))
                    continue
                found.add((float(m.group("p")),
                           bool(m.group("q")) and m.group("q") in (
                               ">", "greater than", "more than", "over")))
    pcts = {p for p, _ in found}
    if len(pcts) != 1:
        return None, ("no protein binding statement" if not pcts
                      else f"protein binding statements disagree: "
                           f"{sorted(pcts)}")
    p = pcts.pop()
    if not 0 < p < 100:
        return None, f"protein binding {p}% cannot give a fraction"
    upper = any(g for pp, g in found if pp == p)
    fu = round((100.0 - p) / 100.0, 4)
    return fu, (f"unbound fraction at most {fu} (label says >{p:g}% bound)"
                if upper else f"unbound fraction {fu} (label says {p:g}% bound)")


# ------------------------------------------------------------------ openFDA
def _generic_names(label: dict) -> list:
    return [str(x) for x in (label.get("openfda") or {}).get("generic_name")
            or []]


def fetch_labels(client, name: str, limit: int = MAX_LABELS) -> list:
    """Single ingredient labels for a drug name, newest first.

    Combination products are left out: their sections describe two drugs and
    their peaks cannot be told apart from the text.
    """
    from . import fda
    base = fda.base_name(name)
    if not base:
        return []
    got = []
    for field in ("generic_name", "substance_name"):
        q = f'openfda.{field}:"{name}"'
        try:
            res = client.get(OPENFDA_LABEL, {"search": q, "limit": limit,
                                             "sort": "effective_time:desc"})
        except Exception as exc:                            # noqa: BLE001
            # 404 is how openFDA says "no match"; anything else is worth a line.
            if "404" not in str(exc):
                log.debug("openFDA label search failed for %s: %s", name,
                          str(exc)[:100])
            continue
        got = (res or {}).get("results") or []
        if got:
            break
    out = []
    for lab in got:
        names = {fda.base_name(x.lower()) for x in _generic_names(lab)}
        subs = {fda.base_name(str(x).lower())
                for x in (lab.get("openfda") or {}).get("substance_name") or []}
        if base in names and len(subs or names) == 1:
            out.append(lab)
    return out


def pk_text(label: dict) -> str:
    """The pharmacokinetics section, or clinical pharmacology without one."""
    for key in ("pharmacokinetics", "clinical_pharmacology"):
        v = label.get(key)
        if v:
            return " ".join(v) if isinstance(v, list) else str(v)
    return ""


def describe(label: dict) -> str:
    o = label.get("openfda") or {}
    brand = (o.get("brand_name") or ["?"])[0]
    route = ", ".join(o.get("route") or []).lower()
    when = str(label.get("effective_time") or "")
    when = f"{when[:4]}-{when[4:6]}-{when[6:8]}" if len(when) >= 8 else when
    return f"FDA label {brand} ({route}), {when}, openFDA set {label.get('set_id', '?')}"


def read_drug(client, name: str, mw: float | None = None) -> dict:
    """One drug's label evidence: a table row, or the reason there is none.

    Returns `{"name", "status", "row", "reason", "accepted", "rejected",
    "labels"}` with status "row", "review" or "none". "row" carries the table
    columns. Labels that state the same figures are counted once.
    """
    from . import exposure
    out = {"name": name, "status": "none", "row": None, "reason": None,
           "accepted": [], "rejected": [], "labels": 0}
    labels = fetch_labels(client, name)
    out["labels"] = len(labels)
    if not labels:
        out["reason"] = "no single ingredient label found in openFDA"
        return out
    acc, rej, src, fu_note, fus = OrderedDict(), [], {}, None, []
    for lab in labels:
        text = pk_text(lab)
        if not text:
            continue
        a, r = candidates(text)
        rej += [dict(x, label=describe(lab)) for x in r]
        for c in a:
            key = (c["value"], c["unit"], c["dose"])
            acc.setdefault(key, dict(c, label=describe(lab), n_labels=0))
            acc[key]["n_labels"] += 1
        fu, note = unbound_fraction(text)
        if fu is not None:
            fus.append((fu, note))
    out["accepted"], out["rejected"] = list(acc.values()), rej
    if not acc:
        out["status"] = "review" if rej else "none"
        out["reason"] = ("every peak statement was set aside: "
                         + "; ".join(sorted({x["why"] for x in rej})[:3])
                         if rej else "no peak concentration sentence in the "
                         "pharmacokinetics section")
        return out
    nm = []
    for c in acc.values():
        v = exposure.to_nm(c["value"], c["unit"], mw)
        if v is not None:
            nm.append((v, c))
    if not nm:
        out["status"] = "review"
        out["reason"] = ("a peak is stated in a mass unit and no molecular "
                         "weight is known to convert it")
        return out
    vals = [v for v, _ in nm]
    if max(vals) / min(vals) > MAX_SPREAD:
        out["status"] = "review"
        out["reason"] = (f"{len(vals)} peak statements range from "
                         f"{min(vals):.3g} to {max(vals):.3g} nM, more than "
                         f"{MAX_SPREAD:g}-fold, so none is chosen")
        return out
    med = statistics.median(vals)
    near = min(nm, key=lambda t: abs(t[0] - med))[1]
    fu = None
    fu_note = None
    if fus and len({f for f, _ in fus}) == 1:
        fu, fu_note = fus[0]
    doses = sorted({c["dose"] for _, c in nm if c["dose"]})
    row = {"name": name, "cmax_nM": round(med, 4),
           "fu": fu if fu is not None else "",
           "dose": ", ".join(doses[:4]),
           "source": near["label"],
           "note": (f"median of {len(vals)} figure(s) in {len(labels)} label(s), "
                    f"range {min(vals):.3g} to {max(vals):.3g} nM; "
                    f"{near['kind']}; "
                    + (fu_note + "; " if fu_note else "")
                    + f"e.g. \"{near['sentence'][:140]}\"")}
    out.update({"status": "row", "row": row})
    return out
