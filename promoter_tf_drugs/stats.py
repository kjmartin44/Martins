"""Statistics shared by more than one analysis stage.

Benjamini-Hochberg correction is used on the LINCS signature scan
(`tfactivity`) and on the cell-line sensitivity correlations
(`sensitivity`). It lived in both as a copy; one copy is easier to keep
right.
"""

from __future__ import annotations

# A q-value that rounds to 0.0 is the most significant result possible, and
# `q or 1.0` would read it as the least - `0.0` is falsy. Every comparison
# goes through a helper that treats None as missing, and the correction
# floors its output here so a reported number never claims to be exactly
# zero either.
Q_FLOOR = 1e-6


def benjamini_hochberg(pvals: list) -> list:
    """Benjamini-Hochberg q-values, order preserved, floored at Q_FLOOR.

    None entries are skipped and stay None in the output; they do not count
    toward the number of tests.
    """
    vals = [(i, p) for i, p in enumerate(pvals) if p is not None]
    if not vals:
        return [None] * len(pvals)
    vals.sort(key=lambda t: t[1])
    m = len(vals)
    out = [None] * len(pvals)
    running = 1.0
    for rank in range(m, 0, -1):
        i, p = vals[rank - 1]
        running = min(running, p * m / rank)
        out[i] = max(Q_FLOOR, round(min(1.0, running), 6))
    return out
