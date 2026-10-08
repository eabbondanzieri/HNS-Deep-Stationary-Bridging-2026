"""WT 96 h (rep1, rep2) ICE maps at 5, 10 and 20 kb from HiC-Pro allValidPairs, averaged over the two replicates like
the 10 kb map the macro uses (cWX2200_96h_rep1rep2_average_10000_iced.tif), recentred so the TIFF starts at oriC
(genome 3,670,000 bp; 3,660,000 at 20 kb, the nearest bin edge). HiC-Pro binning: (pos - 1) // bin.
ICE: merge_and_ice.ice with no bin filtering (HiC-Pro's ICE filters none on these maps). Checked at 10 kb against
HiC-Pro's own ICE output. Also writes the H-NS ChIP z-score (anti-H-NS cWX2200 96 h rep2 / WGS input, per-base
coverage, as ChIP-seq_data/chip_zscore.py) on the same bins, in genome order.
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
import pandas as pd
import tifffile

ROOT = HNS_ROOT + "Bactoloops/"
sys.path.insert(0, ROOT + "iced_matrix_conversion/Dugar_HiC")
from merge_and_ice import ice          # noqa: E402

GLEN, ORI_BP = 4_607_088, 3_670_000
REPS = ["cWX2200_96h_Azure_10percent_Sau3AI_rep1", "cWX2200_96h_Azure_10percent_Sau3AI_rep2"]
pairs = []
for s in REPS:
    d = pd.read_csv(ROOT + f"Hi-C raw count files/allValidPairs/{s}.allValidPairs", sep="\t", header=None, usecols=[2, 5])
    pairs.append((d[2].values - 1, d[5].values - 1))
chip = pd.read_csv(ROOT + "ChIP-seq_data/anti-HNS_cWX2200_96hr_rep2_R1_clean.csv").iloc[:, 1].values.astype(float)
inp = pd.read_csv(ROOT + "ChIP-seq_data/WGS_cWX2200_96hr_rep2_R1_clean.csv").iloc[:, 1].values.astype(float)

for binbp in (5000, 10000, 20000):
    n = -(-GLEN // binbp)
    maps = []
    for a, b in pairs:
        i, j = a // binbp, b // binbp
        M = np.zeros((n, n)); np.add.at(M, (np.minimum(i, j), np.maximum(i, j)), 1)
        M = M + np.triu(M, 1).T
        maps.append(ice(M, filter_fraction=0.0))
    avg = (maps[0] + maps[1]) / 2
    start = ORI_BP // binbp
    r = (start + np.arange(n)) % n
    tifffile.imwrite(ROOT + f"iced_matrix_conversion/multires/wt96_repavg_{binbp // 1000}kb.tif", avg[np.ix_(r, r)].astype(np.float32))
    idx = np.minimum(np.arange(len(chip)) // binbp, n - 1)
    c = np.bincount(idx, weights=chip, minlength=n); w = np.bincount(idx, weights=inp, minlength=n)
    ratio = (c / c.sum()) / (w / w.sum()); z = (ratio - ratio.mean()) / ratio.std()
    pd.DataFrame({"bin": np.arange(n), "start_bp": np.arange(n) * binbp, "zscore": z}).to_csv(
        ROOT + f"iced_matrix_conversion/multires/zscore_hns96_{binbp // 1000}kb.csv", index=False)
    print(f"{binbp // 1000} kb: {n} bins, recenterBin {start + 1}")
    if binbp == 10000:
        ref = tifffile.imread(ROOT + "iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif")
        mine = avg[np.ix_(r, r)]
        print(f"  check vs HiC-Pro 10 kb average: max rel diff {np.max(np.abs(mine - ref) / np.maximum(ref, 1e-9)):.1e}")
        old = pd.read_csv(ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv").zscore.values
        print(f"  ChIP z vs zscore_wt_96h_10kb.csv: r = {np.corrcoef(z, old)[0, 1]:.4f}")
