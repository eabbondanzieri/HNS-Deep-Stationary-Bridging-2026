"""H-NS dependence of looping across growth phases (for a rebuilt Supplementary Fig. 7).

Maps (E. coli, 10 kb, HiC-Pro ICE, recentred TIFFs): WT exp, delta-hns exp, WT 24 h, WT 96 h (rep1/rep2
average), delta-hns 96 h, delta-stpA 96 h. The anti-H-NS ChIP-seq z-score at 96 h stands in for every time point
(H-NS binding is very similar across growth phases and datasets).

Every map is scored on the same scale the macro uses for loop confidence: log2(obs/exp) (missing -> 0, 3x3
average) minus a 50-bin rolling-ball local background, padded circularly. A pixel's "enrichment" is
(value - median at its distance) / (95th percentile at its distance - median), so 0 = the distance-matched
median and 1 = the 95th percentile, comparable across maps.

1. Called loops (the macro, FIJI_output/): number, mean ChIP z of the two anchors, and a rotation null (the
   whole loop set shifted around the chromosome by every offset; rotation_test.py).
2. ChIP-defined candidate loops: every pair of H-NS ChIP peaks (local maxima within +/-2 bins, z >= Z_PEAK)
   at least 100 kb apart, and distance-matched control pairs between low-H-NS bins (z < 0), 20 per peak pair.
   For each map: median enrichment and % above the 95th percentile, peak pairs vs controls, plus a rotation
   test (the peak set shifted around the chromosome, keeping distances).
3. Integrated loop density (the macro's ILD, as a z-score) at the H-NS peak bins vs all other bins.
Output: hns_timecourse_numbers.csv (+ printed table).
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
import sys as _sys
_sys.path.insert(0, HNS_ROOT + "Bactoloops/")
import rotation_test as RT   # exact rotation test: every offset once
import sys
import numpy as np
import pandas as pd
import tifffile

ROOT = HNS_ROOT + "Bactoloops/"
IC = ROOT + "iced_matrix_conversion/"
sys.path.insert(0, ROOT + "python")
from ij_ops import subtract_background, circular_subtract_background          # noqa: E402  (ImageJ's rolling ball, ported exactly)

N, RC, PAD, R_BG, Z_PEAK, MIN_SEP, N_ROT = 461, 368, 100, 50, 2.0, 10, 2000
MAPS = [("WT exp", "cWX2200_exp_Azure_10percent_Sau3A_rep1_10000_iced.tif", "wt_Exp"),
        ("Δhns exp", "cWX2870_exp_Azure_10percent_rep1_10000_iced.tif", "HNSm_Exp"),
        ("WT 24 h", "cWX2200_24H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", "wt_24h"),
        ("WT 96 h", "cWX2200_96h_rep1rep2_average_10000_iced.tif", "wt_96h_repavg"),
        ("Δhns 96 h", "cWX2870_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", "HNSm_96h"),
        ("ΔstpA 96 h", "cWX3338_96h_Sau3AI_rep1_10000_iced.tif", "stpA_96h"),
        ("Lioy WT exp", "lioy_wt_exp_10000_iced.tif", "lioy_wt_exp"),
        ("Lioy WT stat", "lioy_wt_stat_10000_iced.tif", "lioy_wt_stat"),
        ("Lioy WT stat rep1", "lioy_wt_stat_rep1_10000_iced.tif", "lioy_wt_stat_rep1"),
        ("Lioy WT stat avg", "lioy_wt_stat_repavg_10000_iced.tif", "lioy_wt_stat_repavg")]
rng = np.random.default_rng(0)
r = np.arange(N)
D = np.minimum(np.abs(r[:, None] - r[None, :]), N - np.abs(r[:, None] - r[None, :]))
tiff_of = lambda g: (g - (RC - 1)) % N


def corrected(path):
    m = tifffile.imread(path).astype(np.float64)
    E = np.array([m[r, (r + d) % N][m[r, (r + d) % N] > 0].mean() for d in range(N // 2 + 1)])
    oe = np.where(m > 0, m / E[D], 1.0)
    lo = np.log2(oe)
    lo = sum(np.roll(np.roll(lo, a, 0), b, 1) for a in (-1, 0, 1) for b in (-1, 0, 1)) / 9
    lo = lo.astype(np.float32)
    return circular_subtract_background(lo, R_BG).astype(float)


def enrichment(A):
    """(value - distance median) / (distance 95th pct - distance median), per pixel"""
    med = np.array([np.median(A[r, (r + d) % N]) for d in range(N // 2 + 1)])
    p95 = np.array([np.percentile(A[r, (r + d) % N], 95) for d in range(N // 2 + 1)])
    return (A - med[D]) / np.maximum(p95[D] - med[D], 1e-9)


# ChIP (genome order) -> TIFF order
c = pd.read_csv(ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv")
zg = np.full(N, np.nan); zg[(c.start_bp // 10_000).astype(int) % N] = c.zscore.values
z = np.empty(N); z[tiff_of(np.arange(N))] = zg
zz = np.nan_to_num(z, nan=-9)
peaks = np.array([i for i in range(N) if zz[i] >= Z_PEAK and zz[i] == max(zz[(i + k) % N] for k in range(-2, 3))])
PP = np.array([(a, b) for k, a in enumerate(peaks) for b in peaks[k + 1:] if D[a, b] >= MIN_SEP])
low = np.where(zz < 0)[0]; lowset = set(low)
CP = []
for a, b in PP:
    d = D[a, b]
    cand = [(x, (x + d) % N) for x in low if (x + d) % N in lowset]
    CP += [cand[i] for i in rng.choice(len(cand), 20)]
CP = np.array(CP)
print(f"{len(peaks)} H-NS peaks (z >= {Z_PEAK}), {len(PP)} peak pairs >= {MIN_SEP * 10} kb, {len(CP)} control pairs")

rows = []
for name, tif, run in MAPS:
    A = corrected(IC + tif)
    Q = enrichment(A)
    pe, ce = Q[PP[:, 0], PP[:, 1]], Q[CP[:, 0], CP[:, 1]]
    rot = [np.median(Q[(PP[:, 0] + s) % N, (PP[:, 1] + s) % N]) for s in RT.shifts(N)]
    p_rot = RT.p_value(np.median(pe), rot)
    try:
        L = pd.read_csv(ROOT + f"FIJI_output/loop_calls_{run}.csv")
    except pd.errors.EmptyDataError:
        L = pd.DataFrame(columns=["bin1", "bin2"])
    n_loops = len(L)
    if n_loops:
        P = np.c_[L.bin1.round().astype(int), L.bin2.round().astype(int)] % N
        lz = np.nanmean(z[P], axis=1)
        lnull = [np.nanmean(np.nanmean(z[(P + s) % N], axis=1)) for s in RT.shifts(N)]
        loop_z, loop_null95 = np.nanmean(lz), np.percentile(lnull, 95)
        loop_p = RT.p_value(loop_z, lnull)
    else:
        loop_z = loop_null95 = loop_p = np.nan
    ild = pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{run}.csv").ILD_au.values
    ildz = (ild - ild.mean()) / ild.std()
    others = np.setdiff1d(np.arange(N), peaks)
    rows.append(dict(map=name, loops=n_loops, loop_mean_anchor_z=loop_z, loop_rot95=loop_null95, loop_p=loop_p,
                     peakpair_median_enrichment=np.median(pe), peakpair_pct_above95=100 * np.mean(pe > 1),
                     control_median_enrichment=np.median(ce), control_pct_above95=100 * np.mean(ce > 1),
                     peakpair_rot95=np.percentile(rot, 95), peakpair_p=p_rot,
                     ILDz_at_peaks=ildz[peaks].mean(), ILDz_elsewhere=ildz[others].mean(),
                     ILD_vs_ChIP_r=np.corrcoef(ildz[np.isfinite(z)], z[np.isfinite(z)])[0, 1]))
    np.save(ROOT + f"updated figures/_enrichment_{run}.npy", Q.astype(np.float32))
T = pd.DataFrame(rows)
T.to_csv(ROOT + "updated figures/hns_timecourse_numbers.csv", index=False)
pd.set_option("display.width", 250)
print(T.round(3).to_string(index=False))
np.savez(ROOT + "updated figures/_hns_timecourse_pairs.npz", peaks=peaks, PP=PP, CP=CP, z_tiff=z)
