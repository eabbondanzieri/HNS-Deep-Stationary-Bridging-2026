"""Exact rotation (circular-shift) test, shared by the analysis scripts.

A set of loci or pairs on a circular chromosome of n bins is shifted as a block by every offset s (all anchors by the
same s), so the null keeps pair distances, shared anchors and clustering and only moves the set to other genome
positions. Offsets within EXCLUDE bins of zero are left out (the shifted set would still overlap the real one), so the
null has exactly n - 1 - 2 * EXCLUDE members (452 on the 461-bin E. coli genome, 413 on B. subtilis' 422 bins).
Every offset is used once -- no random draws -- so the smallest attainable P is 1 / (members + 1) one-sided and
2 / (members + 1) two-sided (~0.0044 for E. coli): do not report P below that floor.
"""
import numpy as np

EXCLUDE = 4


def shifts(n, exclude=EXCLUDE):
    """all offsets used by the null, in order"""
    return np.arange(exclude + 1, n - exclude)


def p_value(obs, null, two_sided=True):
    """exact rotation P: (1 + #null at least as extreme) / (1 + #null); two-sided = 2 x the smaller tail, capped at 1"""
    null = np.asarray(null, float)
    k = len(null)
    hi = (1 + np.sum(null >= obs)) / (k + 1)
    if not two_sided:
        return hi
    lo = (1 + np.sum(null <= obs)) / (k + 1)
    return min(1.0, 2 * min(hi, lo))


def p_floor(n, two_sided=True, exclude=EXCLUDE):
    k = len(shifts(n, exclude))
    return (2 if two_sided else 1) / (k + 1)


def stars(p):
    """significance stars; *** needs P < 0.001, which a rotation test on a bacterial genome cannot reach"""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."
