"""H-NS at Bactoloops loop anchors (Results: 'H-NS is highly enriched at DNA-loop anchors').
10-kb bins, genome order (461 bins). ILD from the macro (TIFF bin 0 = genome bin 367).
1. Pearson r / Spearman rho between ILD and ChIP z-scores, exact circular-shift P (every shift >= 20 bins).
2. Overlap of bins > 1 SD: ILD z > 1 vs ChIP z > 1 (same bin, and within +/-1 bin).
3. Anchors: loop ends merged with gap 1 (28 for WT 96 h); ChIP z at anchors; strong ChIP peaks (z >= 2, adjacent bins merged)
   that are anchors.
Output: anchor_enrichment_summary.csv"""
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
import numpy as np, pandas as pd
from scipy.stats import spearmanr
B = HNS_ROOT + "Bactoloops/"
sys.path.insert(0, B); import rotation_test as RT
N = 461
def chip(f):
    c = pd.read_csv(B + "ChIP-seq_data/" + f); v = np.full(N, np.nan); v[(c.start_bp // 10000).astype(int) % N] = c.zscore.values; return v
def ild(run):
    t = pd.read_csv(B + f"FIJI_output/integrated_loop_density_{run}.csv").ILD_au.values
    g = np.empty(N); g[(367 + np.arange(N)) % N] = t; return (g - g.mean()) / g.std()
C = {"anti-H-NS 96 h": chip("zscore_wt_96h_10kb.csv"), "anti-mCherry 96 h": chip("zscore_mCherry_96h_10kb.csv"),
     "anti-mCherry exp": chip("zscore_mCherry_Exp_10kb.csv")}
I = {"WT 96 h": ild("wt_96h_repavg"), "WT exp": ild("wt_Exp")}
def shift(a, b):
    ok = np.isfinite(a) & np.isfinite(b); r0 = np.corrcoef(a[ok], b[ok])[0, 1]
    null = []
    for k in RT.shifts(N, exclude=20):
        bb = np.roll(b, k); o = np.isfinite(a) & np.isfinite(bb); null.append(np.corrcoef(a[o], bb[o])[0, 1])
    return r0, spearmanr(a[ok], b[ok])[0], (1 + np.sum(np.abs(null) >= abs(r0))) / (len(null) + 1)
rows = []
for x, y in [("WT 96 h", "anti-H-NS 96 h"), ("WT 96 h", "anti-mCherry 96 h"), ("WT exp", "anti-mCherry exp"),
             ("WT 96 h", "anti-mCherry exp"), ("WT exp", "anti-H-NS 96 h")]:
    r, rho, p = shift(I[x], C[y]); rows.append(dict(a=x + " ILD", b=y, r=r, rho=rho, p=p))
for x, y in [("anti-H-NS 96 h", "anti-mCherry 96 h"), ("anti-mCherry exp", "anti-mCherry 96 h")]:
    r, rho, p = shift(C[x], C[y]); rows.append(dict(a=x, b=y, r=r, rho=rho, p=p))
T = pd.DataFrame(rows); print(T.round(3).to_string(index=False))
cd = lambda a, b: np.minimum((a - b) % N, (b - a) % N)
def near(mask, other, w):
    idx = np.where(other)[0]; return np.array([np.any(cd(i, idx) <= w) for i in np.where(mask)[0]])
for il, ch in [("WT 96 h", "anti-H-NS 96 h"), ("WT 96 h", "anti-mCherry 96 h"), ("WT exp", "anti-mCherry exp")]:
    a, b = I[il] > 1, np.nan_to_num(C[ch]) > 1
    for w in (0, 1):
        print(f"z>1 overlap ({il} ILD / {ch}, +/-{w} bins): {a.sum()} ILD bins, {100*near(a,b,w).mean():.0f}% at ChIP; "
              f"{b.sum()} ChIP bins, {100*near(b,a,w).mean():.0f}% at ILD")
# anchors
L = pd.read_csv(B + "FIJI_output/loop_calls_wt_96h_repavg.csv")
ends = np.unique(np.r_[L.genome_bin1.round(), L.genome_bin2.round()].astype(int) % N)
def groups(bins, gap):
    bins = sorted(bins); out = [[bins[0]]]
    for b in bins[1:]:
        (out[-1].append(b) if b - out[-1][-1] <= gap + 1 else out.append([b]))
    if len(out) > 1 and (out[0][0] + N) - out[-1][-1] <= gap + 1: out[0] = out.pop() + out[0]
    return out
A = groups(ends, 0); print("anchors:", len(A))
z = C["anti-H-NS 96 h"]
az = [np.nanmax([z[(b + k) % N] for b in g for k in (-1, 0, 1)]) for g in A]
print("anchors with max ChIP z (+/-1 bin) >= 1:", sum(v >= 1 for v in az), ">= 2:", sum(v >= 2 for v in az), "of", len(A))
P = groups(np.where(np.nan_to_num(z) >= 2)[0], 0)
inA = [any(cd(b, e) <= 1 for b in g for e in ends) for g in P]
print("strong H-NS peaks (z>=2, merged):", len(P), "at anchors:", sum(inA), "not:", [g for g, f in zip(P, inA) if not f])
zm = C["anti-mCherry 96 h"]; azm = [np.nanmax([zm[(b + k) % N] for b in g for k in (-1, 0, 1)]) for g in A]
Pm = groups(np.where(np.nan_to_num(zm) >= 2)[0], 0); inAm = [any(cd(b, e) <= 1 for b in g for e in ends) for g in Pm]
print("mCherry 96 h: anchors z>=1", sum(v >= 1 for v in azm), "z>=2", sum(v >= 2 for v in azm), "; strong peaks", len(Pm), "at anchors", sum(inAm))
ze = C["anti-mCherry exp"]; aze = [np.nanmax([ze[(b + k) % N] for b in g for k in (-1, 0, 1)]) for g in A]
Pe = groups(np.where(np.nan_to_num(ze) >= 2)[0], 0); inAe = [any(cd(b, e) <= 1 for b in g for e in ends) for g in Pe]
print("mCherry exp: 96 h anchors z>=1", sum(v >= 1 for v in aze), "z>=2", sum(v >= 2 for v in aze), "; strong exp peaks", len(Pe), "at 96 h anchors", sum(inAe))
T.to_csv(B + "anchor_enrichment_summary.csv", index=False)
