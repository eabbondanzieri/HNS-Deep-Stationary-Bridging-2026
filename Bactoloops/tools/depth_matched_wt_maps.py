"""Depth-matched WT 96 h maps for the hns-PAmCherry comparison (2026-10-07).
The hns-PAmCherry 96 h map (cWX2882, one library) sums to 2.85 M; lab HiC-Pro maps sum to exactly 2x the deduplicated >= 5 kb
pair count (checked on the dhns, dstpA and dfis 96 h libraries: ratio 2.00), so it holds ~1.43 M pairs, ~14% of a WT 96 h replicate. Each WT
replicate's validPairs (cis, deduplicated, >= 5 kb, as in validpairs_to_iced_matrix.py) is randomly subsampled to the
same number of contacts (3 draws per replicate), ICE-normalized and written as a recentred 10-kb TIFF.
Output: iced_matrix_conversion/depth_matched/wt96_rep{1,2}_ds{0,1,2}_10000_iced.{matrix,tif}"""
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
import subprocess, sys, numpy as np, pandas as pd   # optional: argv[1] = target pairs, argv[2] = name tag (e.g. 5824556 m58)
ROOT = HNS_ROOT + "Bactoloops/"
VP = ROOT + "Nudler CHIN paper/FullValidPairs/"
OUT = ROOT + "iced_matrix_conversion/depth_matched/"
REPS = {"rep1": "cWX2200_96h_Azure_10percent_Sau3A_rep1", "rep2": "cWX2200_96_H_Azure10percent_Sau3AI_rep2"}
TARGET = int(sys.argv[1]) if len(sys.argv) > 1 else int(round(pd.read_csv(ROOT + "iced_matrix_conversion/cWX2882_96h_Sau3AI_rep1_10000_iced.matrix", sep="\t", header=None)[2].sum() / 2))   # HiC-Pro map total = 2 x pairs
TAG = sys.argv[2] if len(sys.argv) > 2 else ""
BS, G = 10000, 4_607_088; NB = -(-G // BS)


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


print(f"target pairs (hns-PAmCherry 96 h): {TARGET:,}")
for rep, lib in REPS.items():
    P = pd.read_csv(VP + lib + "_EcW3110del.bwt2pairs.validPairs", sep="\t", header=None, usecols=[1, 2, 3, 4, 5, 6],
                    names=["c1", "p1", "s1", "c2", "p2", "s2"])
    P = P[P.c1 == P.c2].drop_duplicates(["p1", "s1", "p2", "s2"]); P = P[(P.p2 - P.p1).abs() >= 5000]
    i_all, j_all = P.p1.values // BS, P.p2.values // BS
    print(f"{rep}: {len(P):,} pairs", flush=True)
    for d in range(3):
        k = np.random.default_rng(100 * d + int(rep[-1])).choice(len(P), TARGET, replace=False)
        X = ice(i_all[k], j_all[k])
        name = f"wt96_{rep}_ds{TAG}{d}_10000_iced"
        r, c = np.triu_indices(NB); v = X[r, c]; nzv = v != 0
        pd.DataFrame({"a": r[nzv] + 1, "b": c[nzv] + 1, "v": v[nzv]}).to_csv(OUT + name + ".matrix", sep="\t", header=False, index=False, float_format="%.6f")
        subprocess.run([sys.executable, ROOT + "iced_matrix_conversion/hic_circular_recenter_to_tiff.py", "--matrix", OUT + name + ".matrix",
                        "--bin-size", "10000", "--out", OUT + name + ".tif", "--recenter-bin", "368"], check=True, capture_output=True)
        print("  wrote", name, flush=True)
