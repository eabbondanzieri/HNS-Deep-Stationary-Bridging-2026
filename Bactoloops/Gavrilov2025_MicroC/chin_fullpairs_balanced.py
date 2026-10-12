"""CHIN / CHID arm contacts from full valid pairs, on a balanced 500-bp band matrix.
Pairs: cis, duplicates removed, inward and same-strand orientations only (outward pairs carry a library-specific re-ligation excess
at short range). Band matrix: 500-bp bins (circular), separations 2..120 bins (1-60 kb). Iterative correction (ICE) on the band
(bins with < 20% of median coverage masked), then O/E = balanced count / mean of its diagonal.
Score per element: mean O/E in the 3 x 3 pixels around (bin of start, bin of end); set score = mean over elements.
Null: the element set shifted by every 2-bin (1-kb) offset from 20 kb (separations kept); ratio = observed / null mean; one-sided P.
Elements: CHINs (246) and CHIDs (20 with ends >= 4 kb apart). Output: chin_fullpairs_balanced.csv"""
import glob, os, numpy as np, pandas as pd
G, BS, ND = 4_607_089, 500, 121
NB = -(-G // BS)
X = pd.read_excel("opcidsRenamed.xlsx", 0)
X.loc[X.OPCID_ID == "OPCID_40", "start_new"] = 2_522_632; X.loc[X.OPCID_ID == "CHID_20", "end_new"] = 3_950_261
X = X.dropna(subset=["start_new", "end_new"]); X["type"] = X.OPCID_ID.str.split("_").str[0]
X["s"], X["e"] = X[["start_new", "end_new"]].min(1).astype(int), X[["start_new", "end_new"]].max(1).astype(int)
def band(f):
    P = pd.read_csv(f, sep="\t", header=None, usecols=[1, 2, 3, 4, 5, 6], names=["c1", "p1", "s1", "c2", "p2", "s2"])
    P = P[P.c1 == P.c2]
    sw = P.p1.values > P.p2.values
    a = np.where(sw, P.p2.values, P.p1.values).astype(np.int64); b = np.where(sw, P.p1.values, P.p2.values).astype(np.int64)
    s1 = np.where(sw, P.s2.values, P.s1.values) == "+"; s2 = np.where(sw, P.s1.values, P.s2.values) == "+"
    _, ui = np.unique(a * 8 * G + b * 4 + s1 * 2 + s2, return_index=True)
    a, b, s1, s2 = a[ui], b[ui], s1[ui], s2[ui]
    keep = ~((~s1) & s2)
    a, b = a[keep], b[keep]
    i, j = a // BS, b // BS
    dd = j - i; wrap = dd > NB // 2
    i = np.where(wrap, j, i); dd = np.where(wrap, NB - dd, dd)              # circular: start from the left end of the short arc
    m = (dd >= 2) & (dd < ND)
    M = np.zeros((NB, ND)); np.add.at(M, (i[m] % NB, dd[m]), 1)
    return M
def ice(M, it=60):
    B = M.copy(); bias = np.ones(NB)
    for _ in range(it):
        cov = B.sum(1).copy()
        for k in range(2, ND): cov += np.roll(B[:, k], k)                     # contacts where the bin is the right end
        med = np.median(cov[cov > 0]); bad = cov < 0.2 * med
        s = np.where(bad, 0, cov / med); s[s == 0] = 1
        B = B / s[:, None]
        for k in range(2, ND): B[:, k] /= np.roll(s, -k)
        bias *= s
        if np.abs(s[~bad] - 1).max() < 1e-3: break
    B[bad] = np.nan
    for k in range(2, ND): B[np.roll(bad, -k), k] = np.nan
    return B
def oe(B):
    E = np.nanmean(B, 0); E[E == 0] = np.nan; return B / E
def score(O, S, E):
    i0, d0 = S // BS, (E // BS - S // BS) % NB                                # separations survive the circular shift
    v = [np.nanmean([O[(i + di) % NB, d + dd - di] for di in (-1, 0, 1) for dd in (-1, 0, 1) if 2 <= d + dd - di < ND])
         for i, d in zip(i0, d0)]
    return np.nanmean(v)
rows = []
for f in sorted(glob.glob("FullValidPairs/*.validPairs")):
    O = oe(ice(band(f))); name = os.path.basename(f).split("_EcW3110del")[0]; r = dict(sample=name)
    for t, mind in (("CHIN", 1000), ("CHID", 4000)):
        C = X[(X.type == t) & ((X.e - X.s) >= mind) & ((X.e - X.s) < (ND - 2) * BS)]
        S, E = C.s.values, C.e.values
        o = score(O, S, E); nul = np.array([score(O, (S + k * BS) % G, (E + k * BS) % G) for k in range(40, NB - 40, 2)])
        r.update({f"{t}_n": len(C), f"{t}_OE": o, f"{t}_ratio": o / nul.mean(), f"{t}_p": (1 + np.sum(nul >= o)) / (len(nul) + 1)})
    rows.append(r); print(pd.DataFrame([r]).round(3).to_string(index=False, header=len(rows) == 1), flush=True)
pd.DataFrame(rows).to_csv("chin_fullpairs_balanced.csv", index=False)
