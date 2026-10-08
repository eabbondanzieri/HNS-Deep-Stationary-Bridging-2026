"""Hairpin ("stem") test per element, from full valid pairs (balanced 500-bp band O/E as chin_fullpairs_balanced.py).
A hairpin folded at the element's centre brings together positions mirror-symmetric about the centre: pixels (c - d/2, c + d/2).
For each element and each separation d (2 bins = 1 kb up to the element length, in 1-bin steps):
  stem      O/E at the mirror-symmetric pixel (midpoint = element centre)
  offcentre O/E at the same separation with the midpoint moved by +/-2 bins (1 kb), kept only if both ends stay inside the element
            (separates a hairpin from general compaction of the element; needs elements >= d + 2 kb)
  slid      O/E at the same separation with the midpoint moved by +/-3..5 bins (1.5-2.5 kb), as in the aggregate maps
Per element: means over d of each class (offcentre over the d where it exists, stem restricted to the same d for that comparison).
Statistics over elements: ratio of means, bootstrap 95% CI (2,000), paired Wilcoxon P.
Elements: CHINs, CHIDs, OPCIDs (OPCIDs = transcription domains, negative control for a hairpin geometry). Output: stem_test_fullpairs.csv"""
import glob, os, numpy as np, pandas as pd
from scipy.stats import wilcoxon
exec(open("chin_fullpairs_balanced.py").read().split("rows = []")[0])
rng = np.random.default_rng(2)
def at(O, left, d):
    return O[left % NB, d] if 2 <= d < ND else np.nan
def element(O, s, e):
    c2 = (s // BS) + (e // BS)                                  # twice the centre bin
    n = e // BS - s // BS
    st, sl, st_oc, oc = [], [], [], []
    for d in range(2, min(n + 1, ND)):
        left = (c2 - d) // 2
        v = at(O, left, d); st.append(v)
        sl.append(np.nanmean([at(O, left + k, d) for k in (-5, -4, -3, 3, 4, 5)]))
        offs = [at(O, left + k, d) for k in (-2, 2) if left + k >= s // BS and left + k + d <= e // BS]
        if offs:
            st_oc.append(v); oc.append(np.nanmean(offs))
    m = lambda x: np.nanmean(x) if len(x) and np.isfinite(x).any() else np.nan
    return m(st), m(sl), m(st_oc), m(oc)
def compare(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5: return np.nan, np.nan, np.nan, np.nan, int(ok.sum())
    a, b = a[ok], b[ok]
    boot = [np.mean(a[i]) / np.mean(b[i]) for i in (rng.integers(0, len(a), len(a)) for _ in range(2000))]
    return np.mean(a) / np.mean(b), np.percentile(boot, 2.5), np.percentile(boot, 97.5), wilcoxon(a, b).pvalue, len(a)
LAB = {"cWX2200_exp_Azure_10percent_Sau3A_rep1": "WT exp r1", "cWX2200_exp_Azure10percent_Sau3AI_rep2": "WT exp r2", "cWX2870_exp_Azure_10percent_rep1": "dhns exp", "cWX2200_24_H_Azure10percent_Sau3AI_rep1": "WT 24 h",
       "cWX2200_96h_Azure_10percent_Sau3A_rep1": "WT 96 h r1", "cWX2200_96_H_Azure10percent_Sau3AI_rep2": "WT 96 h r2", "cWX2870_96HAzure10percent_Sau3AI_rep1": "dhns 96 h",
       "cWX3338_96h_Sau3AI_rep1": "dstpA 96 h r1", "cWX3338_96hr_rep2": "dstpA 96 h r2"}
rows = []
for f in sorted(glob.glob("FullValidPairs/*.validPairs")):
    O = oe(ice(band(f))); name = LAB[os.path.basename(f).split("_EcW3110del")[0]]
    for t in ("CHIN", "CHID", "OPCID"):
        C = X[(X.type == t) & ((X.e - X.s) >= 1500)]
        R = np.array([element(O, s, e) for s, e in zip(C.s.values, C.e.values)])
        r1 = compare(R[:, 0], R[:, 1]); r2 = compare(R[:, 2], R[:, 3])
        rows.append(dict(sample=name, set=t, stem_vs_slid=r1[0], slid_lo=r1[1], slid_hi=r1[2], slid_p=r1[3], n_slid=r1[4],
                         stem_vs_offcentre=r2[0], oc_lo=r2[1], oc_hi=r2[2], oc_p=r2[3], n_oc=r2[4]))
    print(name, flush=True)
T = pd.DataFrame(rows); T.to_csv("stem_test_fullpairs.csv", index=False)
order = list(LAB.values())
for t in ("CHIN", "CHID", "OPCID"):
    S = T[T.set == t].set_index("sample").loc[order]
    print(f"\n{t}: stem vs slid (95% CI, P, n) | stem vs off-centre inside the element (95% CI, P, n)")
    for s, r in S.iterrows():
        print(f"  {s:14s} {r.stem_vs_slid:.2f} ({r.slid_lo:.2f}-{r.slid_hi:.2f}, P {r.slid_p:.2g}, {r.n_slid}) | "
              f"{r.stem_vs_offcentre:.2f} ({r.oc_lo:.2f}-{r.oc_hi:.2f}, P {r.oc_p:.2g}, {r.n_oc})")
