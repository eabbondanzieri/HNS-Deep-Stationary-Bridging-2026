"""Integrated excess corrected for each map's baseline, for comparing loop strengths between maps.

The macro's integrated_excess (step 7c) sums observed - local background over a 5 x 5-bin window, as % of a typical
bin's non-local contacts. The local background comes from a rolling ball, which sits under the data (a lower
envelope), so even random bin pairs have a positive excess, and more so in noisier / shallower maps (e.g. 0.28% for
WT 96 h vs 0.8% for delta-ihfAB at ~1 Mb). Within one map that offset is nearly constant, but it inflates
comparisons between maps. Corrected excess = integrated excess - the map's mean integrated excess over ALL bin pairs
at the same circular distance (the expectation for a pair with no loop). It is ~0 for unlooped pairs in any map.

    B = baseline(tif, cache_name)          # tif: the ICE TIFF the macro ran on (recentred, TIFF bin order)
    B.W[a, b]                              # integrated excess (%) of the pair (a, b), = the macro's value
    B.base[d]                              # baseline at circular distance d (bins)
    B.corrected(a, b)                      # W - base, vectorised over arrays of bins
    corrected_loops(loop_table, B)         # corrected value for every row of a loop_calls_*.csv table
Maps are background-corrected with localbg.py (cached in Supplemental_figure_Ecoli_vs_Bsub/corrected_maps/).
"""
import os as _hns_os
def _hns_root():
    """data/code root: $HNS_ROOT, else the repository folder (marked by .hns_root)"""
    r = _hns_os.environ.get("HNS_ROOT")
    if r:
        return _hns_os.path.join(r, "")
    d = _hns_os.path.dirname(_hns_os.path.abspath(__file__))
    while not _hns_os.path.exists(_hns_os.path.join(d, ".hns_root")) and _hns_os.path.dirname(d) != d:
        d = _hns_os.path.dirname(d)
    return _hns_os.path.join(d, "")
HNS_ROOT = _hns_root()
import sys
import numpy as np
from scipy.ndimage import uniform_filter

sys.path.insert(0, HNS_ROOT + "Bactoloops/Supplemental_figure_Ecoli_vs_Bsub")
import localbg                      # noqa: E402

HALF = 2                            # the macro's integratedExcessHalf: 5 x 5-bin window


class baseline:
    def __init__(self, tif, cache):
        d = localbg.corrected_maps({cache: tif})[cache]
        m, E, bg = d["m"], d["E"], d["bg"]
        n = self.n = len(m)
        i, j = np.indices(m.shape)
        dist = np.minimum(abs(i - j), n - abs(i - j))
        nonlocal_med = np.median(np.where(dist >= 2, m, 0).sum(1)[m.sum(1) > 0])
        X = 100 * np.where(m > 0, m - E * 2.0 ** bg, 0.0) / nonlocal_med
        self.W = uniform_filter(X, size=2 * HALF + 1, mode="wrap") * (2 * HALF + 1) ** 2
        idx = np.arange(n)
        self.base = np.array([self.W[idx, (idx + k) % n].mean() for k in range(n // 2 + 1)])

    def dist(self, a, b):
        a, b = np.asarray(a) % self.n, np.asarray(b) % self.n
        return np.minimum((a - b) % self.n, (b - a) % self.n)

    def corrected(self, a, b):
        a, b = np.asarray(a) % self.n, np.asarray(b) % self.n
        return self.W[a, b] - self.base[self.dist(a, b)]


def corrected_loops(L, B):
    """baseline-corrected integrated excess for each loop (macro value minus the map's distance baseline)"""
    a, b = L.bin1.round().astype(int).values, L.bin2.round().astype(int).values
    return L.integrated_excess.values - B.base[B.dist(a, b)]
