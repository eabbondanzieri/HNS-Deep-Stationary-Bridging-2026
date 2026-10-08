"""Within-domain contact enrichment for OPCIDs (transcription-associated; Gavrilov, Shamovsky et al. 2025), CHIDs and CHINs, from full
valid pairs. Same balanced 500-bp band O/E as chin_fullpairs_balanced.py (inward + same-strand pairs, duplicates removed, ICE).
Domain score per element: mean O/E over all pixels with both ends inside the element and separation 1 kb to the element length.
Local control: the same window (same size) slid by +/- 1 and 2 element lengths; ratio = mean domain / mean local, bootstrap 95% CI
over elements, paired Wilcoxon P. Genome control: the element set shifted by every 10-kb offset; ratio = domain / null mean.
Output: domain_scores_fullpairs.csv"""
import glob, os, numpy as np, pandas as pd
from scipy.stats import wilcoxon
exec(open("chin_fullpairs_balanced.py").read().split("rows = []")[0])
rng = np.random.default_rng(1)
def dom(O, s, L):
    i0 = s // BS; n = max(3, L // BS); vals = []
    for x in range(n):
        for d in range(2, min(n - x, ND)):
            vals.append(O[(i0 + x) % NB, d])
    return np.nanmean(vals) if vals else np.nan
rows = []
for f in sorted(glob.glob("FullValidPairs/*.validPairs")):
    O = oe(ice(band(f))); name = os.path.basename(f).split("_EcW3110del")[0]; r = dict(sample=name)
    for t in ("OPCID", "CHID", "CHIN"):
        C = X[(X.type == t) & ((X.e - X.s) >= 1500) & ((X.e - X.s) < (ND - 2) * BS)]
        S, Lg = C.s.values, (C.e - C.s).values
        dm = np.array([dom(O, s, L) for s, L in zip(S, Lg)])
        lc = np.array([np.nanmean([dom(O, (s + k * L) % G, L) for k in (-2, -1, 1, 2)]) for s, L in zip(S, Lg)])
        ok = np.isfinite(dm) & np.isfinite(lc)
        boot = [np.mean(dm[ok][b]) / np.mean(lc[ok][b]) for b in (rng.integers(0, ok.sum(), ok.sum()) for _ in range(2000))]
        nul = np.array([np.nanmean([dom(O, (s + sh) % G, L) for s, L in zip(S, Lg)]) for sh in range(50_000, G - 50_000, 10_000)])
        r.update({f"{t}_n": int(ok.sum()), f"{t}_dom": dm[ok].mean(), f"{t}_vs_local": dm[ok].mean() / lc[ok].mean(),
                  f"{t}_lo": np.percentile(boot, 2.5), f"{t}_hi": np.percentile(boot, 97.5), f"{t}_p_local": wilcoxon(dm[ok], lc[ok]).pvalue,
                  f"{t}_vs_genome": dm[ok].mean() / nul.mean(), f"{t}_p_genome": (1 + np.sum(nul >= dm[ok].mean())) / (len(nul) + 1)})
    rows.append(r); print(name, flush=True)
T = pd.DataFrame(rows); T.to_csv("domain_scores_fullpairs.csv", index=False)
LAB = {"cWX2200_exp_Azure_10percent_Sau3A_rep1": "WT exp r1", "cWX2200_exp_Azure10percent_Sau3AI_rep2": "WT exp r2", "cWX2870_exp_Azure_10percent_rep1": "dhns exp", "cWX2200_24_H_Azure10percent_Sau3AI_rep1": "WT 24 h",
       "cWX2200_96h_Azure_10percent_Sau3A_rep1": "WT 96 h r1", "cWX2200_96_H_Azure10percent_Sau3AI_rep2": "WT 96 h r2", "cWX2870_96HAzure10percent_Sau3AI_rep1": "dhns 96 h",
       "cWX3338_96h_Sau3AI_rep1": "dstpA 96 h r1", "cWX3338_96hr_rep2": "dstpA 96 h r2"}
T["sample"] = T["sample"].map(LAB); T = T.set_index("sample").loc[list(LAB.values())]
for t in ("OPCID", "CHID", "CHIN"):
    print(f"\n{t} (n = {int(T[t+'_n'].iloc[0])}): domain O/E | vs local (95% CI, P) | vs genome-shifted (P)")
    for s, r in T.iterrows():
        print(f"  {s:14s} {r[t+'_dom']:.2f} | {r[t+'_vs_local']:.2f} ({r[t+'_lo']:.2f}-{r[t+'_hi']:.2f}, P {r[t+'_p_local']:.3g}) | {r[t+'_vs_genome']:.2f} (P {r[t+'_p_genome']:.3f})")
