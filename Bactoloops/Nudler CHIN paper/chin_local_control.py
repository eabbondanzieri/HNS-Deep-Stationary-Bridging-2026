"""Is the CHIN/CHID arm signal hairpin-specific? Same balanced 500-bp band O/E as chin_fullpairs_balanced.py.
For each element: arm pixel = mean O/E in the 3 x 3 pixels at (start, end); local control = the same 3 x 3 score with both ends slid
together by +/-2, 3 and 4 bins (1-2 kb; same separation, overlapping the element's region). Ratio = mean arm / mean local control
over elements; bootstrap 95% CI over elements; paired Wilcoxon P (arm vs local control per element).
Also saves a 41 x 41-pixel aggregate (pileup) of O/E around the arm pixel for each library: chin_pileups_fullpairs.npz."""
import glob, os, numpy as np, pandas as pd
from scipy.stats import wilcoxon
exec(open("chin_fullpairs_balanced.py").read().split("rows = []")[0])
rng = np.random.default_rng(0)
def px(O, i, d):
    v = [O[(i + di) % NB, d + dd - di] for di in (-1, 0, 1) for dd in (-1, 0, 1) if 2 <= d + dd - di < ND]
    return np.nanmean(v) if v else np.nan
def pile(O, i, d, h=20):
    M = np.full((2 * h + 1, 2 * h + 1), np.nan)
    for x in range(-h, h + 1):           # x: shift of the left end, y: shift of the right end
        for y in range(-h, h + 1):
            dd = d + y - x
            if 2 <= dd < ND: M[x + h, y + h] = O[(i + x) % NB, dd]
    return M
rows = []; piles = {}
for f in sorted(glob.glob("FullValidPairs/*.validPairs")):
    O = oe(ice(band(f))); name = os.path.basename(f).split("_EcW3110del")[0]; r = dict(sample=name)
    for t, mind in (("CHIN", 1000), ("CHID", 4000)):
        C = X[(X.type == t) & ((X.e - X.s) >= mind) & ((X.e - X.s) < (ND - 22) * BS)]
        i0, d0 = C.s.values // BS, (C.e.values // BS - C.s.values // BS)
        arm = np.array([px(O, i, d) for i, d in zip(i0, d0)])
        loc = np.array([np.nanmean([px(O, i + k, d) for k in (-4, -3, -2, 2, 3, 4)]) for i, d in zip(i0, d0)])
        ok = np.isfinite(arm) & np.isfinite(loc)
        boot = [np.mean(arm[ok][s]) / np.mean(loc[ok][s]) for s in (rng.integers(0, ok.sum(), ok.sum()) for _ in range(2000))]
        r.update({f"{t}_arm_OE": np.mean(arm[ok]), f"{t}_local_OE": np.mean(loc[ok]), f"{t}_arm/local": np.mean(arm[ok]) / np.mean(loc[ok]),
                  f"{t}_lo": np.percentile(boot, 2.5), f"{t}_hi": np.percentile(boot, 97.5), f"{t}_p": wilcoxon(arm[ok], loc[ok]).pvalue})
        piles[f"{name}|{t}"] = np.nanmean([pile(O, i, d) for i, d in zip(i0, d0)], axis=0)
    rows.append(r); print(pd.DataFrame([r]).round(3).to_string(index=False, header=len(rows) == 1), flush=True)
pd.DataFrame(rows).to_csv("chin_local_control.csv", index=False)
np.savez("chin_pileups_fullpairs.npz", **piles)
