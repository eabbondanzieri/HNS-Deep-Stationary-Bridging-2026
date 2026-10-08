"""Bin-size robustness of Bactoloops (WT 96 h rep1+rep2 average at 5, 10 and 20 kb; runs from tools/make_runs_multires.py
in FIJI_output/extra_runs/multires/). For each resolution: loops, mean H-NS ChIP z of the two anchors (ChIP binned at
the same resolution), ILD-ChIP r, and agreement with the 10-kb calls (anchors within one bin of the coarser map).
Output: FIJI_output/extra_runs/multires/multires_stats.csv"""
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
import numpy as np, pandas as pd
ROOT = HNS_ROOT + "Bactoloops/"
R = ROOT + "FIJI_output/extra_runs/multires/"; MR = ROOT + "iced_matrix_conversion/multires/"
G = 4607088
L = {kb: pd.read_csv(R + f"loop_calls_wt96_repavg_{kb}kb.csv") for kb in (5, 10, 20)}
bp = lambda T: np.c_[T.genome_bp1.values, T.genome_bp2.values].astype(float)
def match(A, B, tol):
    d = lambda a, b: np.minimum(abs(a - b) % G, G - abs(a - b) % G)
    return np.array([(((d(B[:, 0], a) <= tol) & (d(B[:, 1], b) <= tol)) | ((d(B[:, 0], b) <= tol) & (d(B[:, 1], a) <= tol))).any() for a, b in A])
rows = []
P10 = bp(L[10]); q = L[10].integrated_excess >= L[10].integrated_excess.quantile(0.75); hi = (L[10].confidence == "high").values
for kb in (5, 10, 20):
    T = L[kb]; z = pd.read_csv(MR + f"zscore_hns96_{kb}kb.csv").zscore.values
    gb = np.c_[T.genome_bin1.round().astype(int), T.genome_bin2.round().astype(int)] % len(z)
    ild = pd.read_csv(R + f"integrated_loop_density_wt96_repavg_{kb}kb.csv")
    ib = ((ild.genome_position_bp.values // (kb * 1000)).astype(int)) % len(z)
    zz = z[ib]; ok = np.isfinite(zz)
    row = dict(kb=kb, loops=len(T), high=int((T.confidence == "high").sum()), anchor_z=np.nanmean(np.nanmean(z[gb], axis=1)),
               ild_chip_r=np.corrcoef(ild.ILD_au.values[ok], zz[ok])[0, 1])
    if kb != 10:
        tol = max(kb, 10) * 1000
        P = bp(T)
        row.update(pct_match_10kb=100 * match(P, P10, tol).mean(), pct_10kb_recovered=100 * match(P10, P, tol).mean(),
                   pct_top_quartile_recovered=100 * match(P10[q.values], P, tol).mean(), pct_high_recovered=100 * match(P10[hi], P, tol).mean())
    rows.append(row)
D = pd.DataFrame(rows); D.to_csv(R + "multires_stats.csv", index=False); print(D.round(3).to_string(index=False))
