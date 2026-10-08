"""Values for Supplementary Fig. 9 j-l (full valid pairs, no 5-kb cutoff); replicates pooled (each library deduplicated and binned,
raw counts summed, then ICE and O/E per distance; as pooled_stats.py).
  CHIN hairpin test: per CHIN, log2(mean O/E at mirror-symmetric pixels / mean O/E at off-centre pixels with the same separation
    inside the CHIN) (element() in stem_test_fullpairs.py); summary = median over CHINs, bootstrap 95% CI of the median (2,000),
    paired Wilcoxon signed-rank P.
  OPCID control: per OPCID, log2(within-domain O/E / same-size neighbouring windows) (dom() in domain_scores_fullpairs.py); same summary.
  CHID domain: per CHID, log2(within-domain O/E / same-size neighbouring windows); same summary.
  Pileups: mean O/E over CHINs in a 41 x 41-pixel window (+/-10 kb) around the pixel joining the two CHIN ends (pile() as
    chin_local_control.py), from the pooled maps.
Output: s9_fullpairs_stats.csv, s9_pileups.npz"""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
exec(open("stem_test_fullpairs.py").read().split("LAB = {")[0])
src = open("domain_scores_fullpairs.py").read(); exec(src[src.index("def dom("):src.index("rows = []", src.index("def dom("))])
src = open("chin_local_control.py").read(); exec(src[src.index("def pile("):src.index("rows = []", src.index("def pile("))])
FV = "FullValidPairs/"; SUF = "_EcW3110del.bwt2pairs.validPairs"
COND = {"WT exp": ["cWX2200_exp_Azure_10percent_Sau3A_rep1", "cWX2200_exp_Azure10percent_Sau3AI_rep2"],
        "dhns exp": ["cWX2870_exp_Azure_10percent_rep1"],
        "WT 24 h": ["cWX2200_24_H_Azure10percent_Sau3AI_rep1"],
        "WT 96 h": ["cWX2200_96h_Azure_10percent_Sau3A_rep1", "cWX2200_96_H_Azure10percent_Sau3AI_rep2"],
        "dhns 96 h": ["cWX2870_96HAzure10percent_Sau3AI_rep1"],
        "dstpA 96 h": ["cWX3338_96h_Sau3AI_rep1", "cWX3338_96hr_rep2"]}
rng = np.random.default_rng(3)
def summ(x, y):
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0); l = np.log2(x[ok] / y[ok])
    boot = [np.median(rng.choice(l, len(l))) for _ in range(2000)]
    return np.median(l), np.percentile(boot, 2.5), np.percentile(boot, 97.5), wilcoxon(x[ok], y[ok]).pvalue, int(ok.sum())
rows = []; piles = {}
for lab, libs in COND.items():
    O = oe(ice(sum(band(FV + l + SUF) for l in libs)))
    C = X[(X.type == "CHIN") & ((X.e - X.s) >= 1500)]
    R = np.array([element(O, s, e) for s, e in zip(C.s.values, C.e.values)])
    h = summ(R[:, 2], R[:, 3])
    P = X[(X.type == "OPCID") & ((X.e - X.s) >= 1500) & ((X.e - X.s) < (ND - 2) * BS)]
    dm = np.array([dom(O, s, L) for s, L in zip(P.s.values, (P.e - P.s).values)])
    lc = np.array([np.nanmean([dom(O, (s + k * L) % G, L) for k in (-2, -1, 1, 2)]) for s, L in zip(P.s.values, (P.e - P.s).values)])
    o = summ(dm, lc)
    Q = X[(X.type == "CHID") & ((X.e - X.s) < (ND - 2) * BS)]
    dq = np.array([dom(O, s, L) for s, L in zip(Q.s.values, (Q.e - Q.s).values)])
    lq = np.array([np.nanmean([dom(O, (s + k * L) % G, L) for k in (-2, -1, 1, 2)]) for s, L in zip(Q.s.values, (Q.e - Q.s).values)])
    cd = summ(dq, lq)
    Cp = X[(X.type == "CHIN") & ((X.e - X.s) >= 1000) & ((X.e - X.s) < (ND - 22) * BS)]
    piles[lab] = np.nanmean([pile(O, s // BS, e // BS - s // BS) for s, e in zip(Cp.s.values, Cp.e.values)], axis=0)
    rows.append(dict(condition=lab, libraries=len(libs), chin_log2=h[0], chin_lo=h[1], chin_hi=h[2], chin_p=h[3], chin_n=h[4],
                     opcid_log2=o[0], opcid_lo=o[1], opcid_hi=o[2], opcid_p=o[3], opcid_n=o[4],
                     chid_log2=cd[0], chid_lo=cd[1], chid_hi=cd[2], chid_p=cd[3], chid_n=cd[4]))
    print(f"{lab:11s} ({len(libs)} lib)  CHIN hairpin {2**h[0]:.2f}x ({2**h[1]:.2f}-{2**h[2]:.2f}, P {h[3]:.2g}, n {h[4]})  |  "
          f"OPCID domain {2**o[0]:.2f}x ({2**o[1]:.2f}-{2**o[2]:.2f}, P {o[3]:.2g}, n {o[4]})  |  CHID domain {2**cd[0]:.2f}x ({2**cd[1]:.2f}-{2**cd[2]:.2f}, P {cd[3]:.2g}, n {cd[4]})", flush=True)
pd.DataFrame(rows).to_csv("s9_fullpairs_stats.csv", index=False); np.savez("s9_pileups.npz", **piles)
