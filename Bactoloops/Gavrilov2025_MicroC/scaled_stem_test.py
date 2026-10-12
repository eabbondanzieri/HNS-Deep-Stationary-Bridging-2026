"""Hairpin test on per-CHIN rescaled maps (pre-specified regions, fixed before running; see below).
Per CHIN (>= 2 kb; 141): 41 x 41 map of O/E sampled at t = (position - centre) / half-length in [-3, 3] (as scaled_pileups.py), from
the pooled balanced 500-bp maps; four-fold average (both diagonal reflections, 180-degree rotation) and a 2-sample Gaussian blur
(NaN-aware), as in Supplementary Fig. 9j. Regions, with a < b the two positions (t units), separation s = b - a, midpoint m = (a + b) / 2:
  stem   0.5 <= s <= 2 and |m| <= 0.25          (anti-diagonal band through the centre: a hairpin folded at the CHIN centre)
  flank  0.5 <= s <= 2 and 0.5 <= |m| <= 1.0    (same separations, midpoint shifted off-centre)
Score per CHIN = log2(mean O/E in stem / mean O/E in flank). Summary: median, bootstrap 95% CI (2,000), Wilcoxon signed-rank P vs 0.
Paired genotype comparisons on the same CHINs: WT exp vs dhns exp, WT 96 h vs dhns 96 h, WT 96 h vs WT exp, dstpA 96 h vs WT 96 h
(median paired difference of scores, bootstrap CI, Wilcoxon P). Output: scaled_stem_test.csv"""
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.stats import wilcoxon
exec(open("stem_test_fullpairs.py").read().split("LAB = {")[0])
FV = "FullValidPairs/"; SUF = "_EcW3110del.bwt2pairs.validPairs"
COND = {"WT exp": ["cWX2200_exp_Azure_10percent_Sau3A_rep1", "cWX2200_exp_Azure10percent_Sau3AI_rep2"],
        "dhns exp": ["cWX2870_exp_Azure_10percent_rep1"],
        "WT 24 h": ["cWX2200_24_H_Azure10percent_Sau3AI_rep1"],
        "WT 96 h": ["cWX2200_96h_Azure_10percent_Sau3A_rep1", "cWX2200_96_H_Azure10percent_Sau3AI_rep2"],
        "dhns 96 h": ["cWX2870_96HAzure10percent_Sau3AI_rep1"],
        "dstpA 96 h": ["cWX3338_96h_Sau3AI_rep1", "cWX3338_96hr_rep2"]}
C = X[(X.type == "CHIN") & ((X.e - X.s) >= 2000)]
T = np.linspace(-3, 3, 41)
A, B = np.meshgrid(T, T, indexing="ij"); SEP = np.abs(B - A); MID = (A + B) / 2
STEM = (SEP >= 0.5) & (SEP <= 2) & (np.abs(MID) <= 0.25)
FLANK = (SEP >= 0.5) & (SEP <= 2) & (np.abs(MID) >= 0.5) & (np.abs(MID) <= 1.0)
def blur(M, sg=2):
    w = np.isfinite(M).astype(float); v = np.where(w > 0, M, 0)
    num, den = gaussian_filter(v, sg, mode="constant"), gaussian_filter(w, sg, mode="constant")
    out = num / np.where(den > 0, den, np.nan); out[w == 0] = np.nan; return out
def fourfold(M):
    with np.errstate(all="ignore"): return np.nanmean(np.stack([M, M.T, M[::-1, ::-1].T, M[::-1, ::-1]]), axis=0)
rng = np.random.default_rng(4)
def summ(v):
    v = v[np.isfinite(v)]; b = [np.median(rng.choice(v, len(v))) for _ in range(2000)]
    return np.median(v), np.percentile(b, 2.5), np.percentile(b, 97.5), wilcoxon(v).pvalue, len(v)
scores = {}
for lab, libs in COND.items():
    O = oe(ice(sum(band(FV + l + SUF) for l in libs))); sc = []
    for s, e in zip(C.s.values, C.e.values):
        c, h = (s + e) / 2, (e - s) / 2
        bb = np.floor((c + T * h) / BS).astype(int)
        M = np.full((len(T), len(T)), np.nan)
        for i in range(len(T)):
            for j in range(len(T)):
                d = bb[j] - bb[i]; lo = bb[i] if d >= 0 else bb[j]; d = abs(d)
                if 2 <= d < ND: M[i, j] = O[lo % NB, d]
        M = blur(fourfold(M))
        st, fl = np.nanmean(M[STEM]), np.nanmean(M[FLANK])
        sc.append(np.log2(st / fl) if st > 0 and fl > 0 else np.nan)
    scores[lab] = np.array(sc); print(lab, flush=True)
rows = []
print("\nScaled stem test, per CHIN (stem / flank):")
for lab, v in scores.items():
    m, lo, hi, p, n = summ(v); rows.append(dict(test=lab, median_log2=m, lo=lo, hi=hi, p=p, n=n))
    print(f"  {lab:11s} {2**m:.2f}x ({2**lo:.2f}-{2**hi:.2f}), P {p:.2g}, n {n}")
print("\nPaired comparisons on the same CHINs (ratio of scores):")
for a, b in [("WT exp", "dhns exp"), ("WT 96 h", "dhns 96 h"), ("WT 96 h", "WT exp"), ("WT 24 h", "WT exp"), ("dstpA 96 h", "WT 96 h")]:
    m, lo, hi, p, n = summ(scores[a] - scores[b]); rows.append(dict(test=f"{a} vs {b}", median_log2=m, lo=lo, hi=hi, p=p, n=n))
    print(f"  {a} vs {b}: {2**m:.2f}x ({2**lo:.2f}-{2**hi:.2f}), P {p:.2g}, n {n}")
pd.DataFrame(rows).to_csv("scaled_stem_test.csv", index=False)
