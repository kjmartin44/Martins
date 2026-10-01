"""Collapse salt forms so the tables show distinct chemistry.

The problem this fixes
----------------------
ChEMBL lists every marketed form as its own entry. A CCR5 run returned, as
five separate "candidates" scoring identically:

    Dextrothyroxine · Dextrothyroxine Sodium · Liothyronine ·
    Liothyronine Sodium · Levothyroxine Sodium

That is one pharmacology occupying five rows. The cost is not cosmetic. Rows
are the scarce resource — tables print the top 25 — so duplicates push distinct
chemistry off the bottom of every table, and they inflate the candidate count
that the Pareto accounting reports, which makes a narrow result look broad.

Two kinds of relatedness, deliberately treated differently
----------------------------------------------------------
**Salts** are a counterion around an identical active moiety. Liothyronine
Sodium *is* liothyronine once it dissolves. These are collapsed: one row, with
the alternates named.

**Esters and prodrugs** are not. Fluticasone Propionate and Fluticasone Furoate
are different molecules with different potency, duration and formulation, and
the furoate was developed precisely because it is not the propionate. They are
*grouped and flagged* so you can see the family, never merged — merging them
would be a pharmacological claim this module has no business making.

Nothing is dropped. The collapsed names travel with the surviving row, because
the specific salt is what you order from a catalogue.
"""

from __future__ import annotations

import re

# Counterions and hydration states: the active moiety is unchanged.
SALTS = {
    "sodium", "potassium", "calcium", "magnesium", "zinc", "aluminum",
    "aluminium", "lithium", "ammonium",
    "lysine", "meglumine", "olamine", "diolamine", "tromethamine", "choline",
    "arginine", "benzathine", "procaine", "erbumine", "trometamol",
    "hydrochloride", "dihydrochloride", "hydrobromide", "hydroiodide",
    "hydrate", "dihydrate", "trihydrate", "monohydrate", "anhydrous",
    "sulfate", "sulphate", "bisulfate", "hemisulfate", "disulfate",
    "phosphate", "diphosphate", "biphosphate",
    "nitrate", "bromide", "chloride", "iodide", "fluoride",
    "mesylate", "besylate", "tosylate", "esylate", "camsylate", "edisylate",
    "napadisylate", "xinafoate", "isethionate", "estolate",
    "citrate", "dicitrate", "tartrate", "bitartrate", "maleate", "fumarate",
    "succinate", "malate", "oxalate", "gluconate", "lactate", "mucate",
    "pamoate", "embonate", "salicylate", "benzoate", "carbonate",
    "bicarbonate", "trifluoroacetate", "edetate", "aspartate", "glutamate",
    "saccharate", "adipate", "malonate", "sebacate",
}

# Esters, ethers and prodrug promoieties: a DIFFERENT molecule. Grouped for
# display, never merged.
ESTERS = {
    "acetate", "diacetate", "propionate", "dipropionate", "butyrate",
    "valerate", "caproate", "enanthate", "decanoate", "undecylenate",
    "undecanoate", "palmitate", "stearate", "laurate", "cypionate",
    "pivalate", "furoate", "aceponate", "acetonide", "hexacetonide",
    "dipivoxil", "axetil", "proxetil", "mofetil", "fosil", "medoxomil",
    "olamine", "sodium succinate", "phosphate sodium",
}

# Words that merely qualify the preparation.
NOISE = {"usp", "anhydrous", "micronized", "micronised", "purified"}

_SPLIT = re.compile(r"[\s,]+")


def _tokens(name: str) -> list[str]:
    return [t for t in _SPLIT.split(name.strip().lower()) if t]


def salt_free_name(name: str) -> str:
    """Strip trailing counterion and hydration words.

    Two rules, and the second is not obvious. Only *trailing* tokens are
    stripped, which preserves the head noun — the active moiety in every naming
    convention used here.

    But that alone is not enough. "Sodium Fluoride" and "Potassium Chloride"
    are drugs in their own right whose *every* token is a counterion word, so
    trailing-strip alone reduces them to "sodium" and "potassium" — silently
    merging two unrelated salts into one row. So a strip is refused when what
    survives is nothing but counterion words: a name made only of ions is the
    ion pair itself, not a salt of something else.
    """
    toks = _tokens(name)
    while len(toks) > 1 and (toks[-1] in SALTS or toks[-1] in NOISE):
        head = toks[:-1]
        if all(t in SALTS or t in NOISE for t in head):
            break          # the remainder is itself only ions — this IS the drug
        toks = head
    return " ".join(toks) if toks else name.strip().lower()


def ester_free_name(name: str) -> str:
    """Strip trailing ester/prodrug promoieties as well as salts.

    Used only to *group* a family for display. Never used to merge rows.
    """
    toks = _tokens(salt_free_name(name))
    while len(toks) > 1 and (toks[-1] in ESTERS or toks[-1] in SALTS):
        toks.pop()
    return " ".join(toks) if toks else name.strip().lower()


def is_salt_form(name: str) -> bool:
    return salt_free_name(name) != " ".join(_tokens(name))


def _display(names: list[str]) -> str:
    """Prefer the shortest name: the parent, not the marketed salt."""
    return min(names, key=lambda n: (len(n), n))


def collapse_grades(grades: list, *, key=None, name_of=None,
                    set_name=None) -> tuple[list, dict]:
    """Collapse salt forms of the same drug at the same target.

    Grouped on (salt-free name, target), so the same compound acting at two
    different TFs stays as two rows — those are genuinely different routes.
    The highest-scoring member survives and carries the alternates.

    Returns (kept, info) where info maps the surviving row's name to the list
    of collapsed alternates.
    """
    name_of = name_of or (lambda g: g.drug_name)
    key = key or (lambda g: (salt_free_name(name_of(g)),
                             getattr(g, "tf_symbol", None)))

    groups: dict = {}
    for g in grades:
        groups.setdefault(key(g), []).append(g)

    kept, info = [], {}
    for _, members in groups.items():
        if len(members) == 1:
            kept.append(members[0])
            continue
        # Order is already best-first from grade_all; keep that winner rather
        # than re-deriving a score here, so this cannot disagree with the
        # ranking that produced it.
        winner = members[0]
        others = [name_of(m) for m in members[1:]]
        kept.append(winner)
        info[name_of(winner)] = others
    return kept, info


def family_of(names) -> dict:
    """Group names into ester/prodrug families, for display only."""
    fams: dict = {}
    for n in names:
        fams.setdefault(ester_free_name(n), []).append(n)
    return {k: sorted(v) for k, v in fams.items() if len(v) > 1}
