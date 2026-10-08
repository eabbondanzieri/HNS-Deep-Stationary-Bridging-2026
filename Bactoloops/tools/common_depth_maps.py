"""Like-for-like maps for the hns-PAmCherry vs WT looping comparison (2026-10-07).
Every map is built from HiC-Pro valid pairs with the same steps (validpairs_to_iced_matrix.py: cis, deduplicated, >= 5 kb,
10-kb bins, ICE) and then subsampled to a common depth, the smallest library (hns-PAmCherry rep3, 5,732,683 pairs):
  full maps:      WT rep1, WT rep2, PAmCherry rep2, PAmCherry rep3
  common depth:   WT rep1, WT rep2, PAmCherry rep2 subsampled 3 times each to 5,732,683 pairs; PAmCherry rep3 as is
Output: iced_matrix_conversion/common_depth/<name>_10000_iced.{matrix,tif} (recentred at bin 368)"""
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
import subprocess, sys, os, numpy as np, pandas as pd
ROOT = HNS_ROOT + "Bactoloops/"
VP = ROOT + "Nudler CHIN paper/FullValidPairs/"; OUT = ROOT + "iced_matrix_conversion/common_depth/"; os.makedirs(OUT, exist_ok=True)
LIBS = {"wt96_rep1": "cWX2200_96h_Azure_10percent_Sau3A_rep1", "wt96_rep2": "cWX2200_96_H_Azure10percent_Sau3AI_rep2",
        "pam96_rep2": "cWX2882_96hr_rep2", "pam96_rep3": "cWX2882_96hr_rep3"}
BS, G = 10000, 4_607_088; NB = -(-G // BS)
# ICE as in validpairs_to_iced_matrix.py (copied from depth_matched_wt_maps.py)
def ice(i, j):
    M = np.zeros((NB, NB)); np.add.at(M, (i, j), 1); np.add.at(M, (j, i), 1)
    M[np.arange(NB), np.arange(NB)] /= 2
    X = M.copy(); total = np.triu(X).sum(); nz = (X != 0).sum(0); bad = nz == 0
    bad |= nz < np.percentile(nz[~bad], 2)
    X[bad, :] = 0; X[:, bad] = 0
    for it in range(100):
        s = X.sum(0); m = s[s != 0].mean(); s = s / m; s[s == 0] = 1
        X = X / s[:, None] / s[None, :]
        if np.abs(s - 1).max() < 0.1: break
    X *= total / np.triu(X).sum()
    return X




def write(name, i, j):
    X = ice(i, j); r, c = np.triu_indices(NB); v = X[r, c]; k = v != 0
    pd.DataFrame({"a": r[k] + 1, "b": c[k] + 1, "v": v[k]}).to_csv(OUT + name + "_10000_iced.matrix", sep="\t", header=False, index=False, float_format="%.6f")
    subprocess.run([sys.executable, ROOT + "iced_matrix_conversion/hic_circular_recenter_to_tiff.py", "--matrix", OUT + name + "_10000_iced.matrix",
                    "--bin-size", "10000", "--out", OUT + name + "_10000_iced.tif", "--recenter-bin", "368"], check=True, capture_output=True)


pairs = {}
for k, lib in LIBS.items():
    P = pd.read_csv(VP + lib + "_EcW3110del.bwt2pairs.validPairs", sep="\t", header=None, usecols=[1, 2, 3, 4, 5, 6],
                    names=["c1", "p1", "s1", "c2", "p2", "s2"])
    P = P[P.c1 == P.c2].drop_duplicates(["p1", "s1", "p2", "s2"]); P = P[(P.p2 - P.p1).abs() >= 5000]
    pairs[k] = (P.p1.values // BS, P.p2.values // BS); print(f"{k}: {len(P):,} pairs", flush=True)
TARGET = min(len(v[0]) for v in pairs.values()); print(f"common depth: {TARGET:,}")
for k, (i, j) in pairs.items():
    write(f"{k}_full", i, j)
    for d in range(3):
        if len(i) == TARGET:
            if d == 0: write(f"{k}_common0", i, j)
            continue
        s = np.random.default_rng(1000 * d + 7 * len(k) + int(k[-1])).choice(len(i), TARGET, replace=False)
        write(f"{k}_common{d}", i[s], j[s])
    print("  wrote", k, flush=True)
