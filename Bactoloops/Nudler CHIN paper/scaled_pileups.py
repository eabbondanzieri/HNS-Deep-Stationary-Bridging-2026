"""CHIN aggregate maps centred on the CHIN centre and rescaled by CHIN length: axis t = (position - centre) / (half-length), so the two
CHIN ends sit at t = -1 and +1 for every CHIN (t from -3 to 3, 41 samples). Each sample takes the O/E of the 500-bp pixel containing it
(pooled balanced maps as s9_fullpairs_stats.py); pixels closer than 1 kb to the main diagonal are blank. CHINs >= 2 kb.
A hairpin folded at the centre gives a stripe along the anti-diagonal ending at (-1, +1). Output: scaled_pileups.npz"""
import numpy as np
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
out = {}
for lab, libs in COND.items():
    O = oe(ice(sum(band(FV + l + SUF) for l in libs)))
    acc = []
    for s, e in zip(C.s.values, C.e.values):
        c, h = (s + e) / 2, (e - s) / 2
        b = np.floor((c + T * h) / BS).astype(int)
        M = np.full((len(T), len(T)), np.nan)
        for i in range(len(T)):
            for j in range(len(T)):
                d = b[j] - b[i]
                if d < 0: d = -d; lo = b[j]
                else: lo = b[i]
                if 2 <= d < ND: M[i, j] = O[lo % NB, d]
        acc.append(M)
    out[lab] = np.nanmean(acc, axis=0); print(lab, len(acc), flush=True)
np.savez("scaled_pileups.npz", **out)
