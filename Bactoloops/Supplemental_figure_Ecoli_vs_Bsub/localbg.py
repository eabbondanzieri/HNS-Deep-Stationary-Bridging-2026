"""Local-background-corrected log2(obs/exp) maps, shared by the supplemental-figure scripts.

Distance-only log2(obs/exp) keeps large-scale structure that isn't a function of
genomic distance (in B. subtilis: the SMC arm-alignment X and the oriC/terminus
domains). The corrected map removes it the way it was done by hand in FIJI:
Process > Subtract Background, rolling ball radius 50 px (500 kb), default options,
applied to the macro's log2(obs/exp) map (HiC_map_log2OE: missing pixels -> 0,
3x3 average). Since 2026-10-06 this is the macro's rotation-invariant version
(ij_ops.periodic_background in ../python: ImageJ's rolling ball on a periodic copy
of the map, averaged over all 16 shrink-grid phases and symmetrized), computed in
Python, so the result does not depend on where the circular map starts.

    maps = corrected_maps({"name": ice_matrix_tif, ...})
    maps["name"] -> dict(raw=distance-only map, bs=corrected map, bg=raw - bs, E=expected per pixel)

Corrected maps are cached in corrected_maps/ next to this file (as float32 TIFFs)
and rebuilt only if missing.
"""
import os
import sys
import numpy as np
import tifffile

HERE = os.path.dirname(os.path.abspath(__file__)) + "/"
CACHE = HERE + "corrected_maps/"
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import ij_ops  # noqa: E402
PAD = 100
RADIUS = 50


def expected(m):
    n = len(m)
    i, j = np.indices(m.shape)
    d = np.minimum(abs(i - j), n - abs(i - j))
    ok = m > 0
    e = np.bincount(d[ok], m[ok], minlength=n // 2 + 1) / np.maximum(np.bincount(d[ok], minlength=n // 2 + 1), 1)
    return e[d], ok


def log2oe_macro(m):
    """the macro's HiC_map_log2OE: log2(obs/exp), missing -> 0, 3x3 average (circular)"""
    E, ok = expected(m)
    lo = np.where(ok & (E > 0), np.log2(np.where(ok, m, 1) / np.where(E > 0, E, 1)), 0.0)
    return sum(np.roll(np.roll(lo, a, 0), b, 1) for a in (-1, 0, 1) for b in (-1, 0, 1)) / 9


def corrected_maps(sources, radius=RADIUS):
    os.makedirs(CACHE, exist_ok=True)
    out = {}
    for name, tif in sources.items():
        m = tifffile.imread(tif).astype(float)
        raw = log2oe_macro(m)
        E, _ = expected(m)
        out[name] = dict(m=m, raw=raw, E=E)
        f = CACHE + f"{name}_log2oe_rb{radius}_rotinv.tif"
        if not os.path.exists(f):
            n = len(raw)
            padded = np.pad(raw.astype(np.float32), PAD, mode="wrap")
            bg = ij_ops.periodic_background(padded, n, PAD, radius)
            tifffile.imwrite(f, (padded - bg)[PAD:PAD + n, PAD:PAD + n].astype(np.float32))
        out[name]["bs"] = tifffile.imread(f).astype(float)
        out[name]["bg"] = raw - out[name]["bs"]   # the large-scale background the rolling ball removed (log2 units)
    return out


def distance_percentile(A, a, b):
    """fraction of pixels at the same circular distance with a lower value than A[a, b]"""
    n = len(A)
    d = abs(a - b) % n
    d = min(d, n - d)
    r = np.arange(n)
    v = A[r, (r + d) % n]
    return np.mean(v < A[a, b]), np.median(v)


def score_loop_table(loops_csv, A, out_csv, high=0.99, medium=0.95):
    """copy of a macro loop table with corrected-map columns added:
    log2_obs_exp_local, oe_percentile_local, confidence_local"""
    import pandas as pd
    L = pd.read_csv(loops_csv)
    vals, pcts = [], []
    for a, b in zip(L.bin1.round().astype(int), L.bin2.round().astype(int)):
        vals.append(A[a, b])
        pcts.append(distance_percentile(A, a, b)[0])
    L["log2_obs_exp_local"] = np.round(vals, 4)
    L["oe_percentile_local"] = np.round(pcts, 4)
    L["confidence_local"] = np.where(L.oe_percentile_local >= high, "high",
                                     np.where(L.oe_percentile_local >= medium, "medium", "low"))
    L.to_csv(out_csv, index=False)
    return L


def excess_fraction_map(d):
    """excess contacts over the LOCAL background, as a fraction of an anchor's non-local contacts:
    (observed - expected at that distance x the large-scale factor the rolling ball removed) / row total.
    d is one entry of corrected_maps()."""
    m, E, bg = d["m"], d["E"], d["bg"]
    n = len(m)
    i, j = np.indices(m.shape)
    dist = np.minimum(abs(i - j), n - abs(i - j))
    nonlocal_sum = np.median(np.where(dist >= 2, m, 0).sum(1)[m.sum(1) > 0])
    return np.where(m > 0, m - E * 2.0 ** bg, 0.0) / nonlocal_sum


def integrated_excess(X, a, b, half=2):
    """sum of X over a (2*half+1)^2 window around (a, b), wrapping circularly; x100 = % of anchor contacts"""
    n = len(X)
    rr = np.arange(a - half, a + half + 1) % n
    cc = np.arange(b - half, b + half + 1) % n
    return 100 * X[np.ix_(rr, cc)].sum()


def add_integrated_excess(table_csv, d, half=2):
    """add integrated_excess_local (% of anchor contacts, 50 x 50 kb window at 10 kb bins) to a corrected loop table"""
    import pandas as pd
    L = pd.read_csv(table_csv)
    X = excess_fraction_map(d)
    L["integrated_excess_local"] = [round(integrated_excess(X, a, b, half), 4)
                                    for a, b in zip(L.bin1.round().astype(int), L.bin2.round().astype(int))]
    L.to_csv(table_csv, index=False)
    return L
