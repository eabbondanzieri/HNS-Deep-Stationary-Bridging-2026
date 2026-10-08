"""Merge Dugar et al. B. subtilis Hi-C replicates and ICE-balance them.

Raw .cool files (integer counts, 10 kb bins, NC_000964.3) from GEO:
    GSM4288936  Hi-C-wt_st_1    GSM4818084  Hi-C-wt_st_2
    GSM4288937  Hi-C-rok_st_1   GSM4818085  Hi-C-rok_st_2

Merging = summing the raw counts of the two runs, then balancing the sum.
Balancing reproduces the GEO ICE matrix for GSM4288936 (r = 0.99997, the same
11 filtered bins): drop the lowest-coverage 2.5% of bins, then iterative
correction until every remaining row sum is equal.

Writes, for each condition, a float32 TIFF (not recentered; recenter bin 0)
for the ImageJ macro, next to this script:
    Dugar_WT_st_merged_iced.tif, Dugar_rok_st_merged_iced.tif
and the summed raw counts as a .cool in "cool and bed files/".
Run with --qc to print replicate agreement instead.
"""
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
import numpy as np
import pandas as pd
import cooler
import tifffile

ROOT = HNS_ROOT + "Bactoloops/"
COOLS = ROOT + "cool and bed files/"
OUT = ROOT + "iced_matrix_conversion/Dugar_HiC/"
RUNS = {"WT": ["GSM4288936_WT_st_HiC.cool", "GSM4818084_WT_st_HiC_run2.cool"],
        "rok": ["GSM4288937_Rok_st_HiC.cool", "GSM4818085_Rok_st_HiC_run2.cool"]}
FILTER_FRACTION = 0.025


def raw(f):
    return cooler.Cooler(COOLS + f).matrix(balance=False)[:].astype(float)


def ice(m, filter_fraction=FILTER_FRACTION, max_iter=500, eps=1e-6):
    m = m.copy()
    cov = m.sum(1)
    bad = cov < np.percentile(cov, 100 * filter_fraction)
    m[bad, :] = 0
    m[:, bad] = 0
    for _ in range(max_iter):
        s = m.sum(1)
        good = s > 0
        b = np.ones(len(m))
        b[good] = s[good] / s[good].mean()
        m = m / b[:, None] / b[None, :]
        if np.abs(b[good] - 1).max() < eps:
            break
    return m


def obs_exp_log2(m):
    n = len(m)
    i, j = np.indices(m.shape)
    d = np.minimum(abs(i - j), n - abs(i - j))
    ok = m > 0
    e = np.bincount(d[ok], m[ok], minlength=n // 2 + 1) / np.maximum(np.bincount(d[ok], minlength=n // 2 + 1), 1)
    lo = np.where(ok, np.log2(np.where(ok, m, 1) / np.where(e[d] > 0, e[d], 1)), 0.0)
    return sum(np.roll(np.roll(lo, a, 0), b, 1) for a in (-1, 0, 1) for b in (-1, 0, 1)) / 9


def scc(a, b, dmax=200):
    sm = lambda m: sum(np.roll(np.roll(m, u, 0), v, 1) for u in (-1, 0, 1) for v in (-1, 0, 1)) / 9
    a, b = sm(a), sm(b)
    n = len(a)
    r = np.arange(n)
    num = den = 0.0
    for d in range(2, dmax + 1):
        x, y = a[r, (r + d) % n], b[r, (r + d) % n]
        ok = (x > 0) & (y > 0)
        w = ok.sum() * x[ok].std() * y[ok].std()
        num += w * np.corrcoef(x[ok], y[ok])[0, 1]
        den += w
    return num / den


if __name__ == "__main__" and "--qc" in sys.argv:
    RC = {"RC1": 54, "RC2": 206, "RC3": 226, "RC6": 366, "RC7": 412, "RC8": 417}   # Dugar Supp. Table 3
    pairs = [("RC1", "RC6"), ("RC1", "RC7"), ("RC6", "RC8"), ("RC2", "RC3")]
    maps = {f"{c} run{k + 1}": ice(raw(f)) for c, fs in RUNS.items() for k, f in enumerate(fs)}
    print("contacts:", {f"{c} run{k + 1}": f"{raw(f).sum() / 1e6:.1f} M" for c, fs in RUNS.items() for k, f in enumerate(fs)})
    names = list(maps)
    print("\nstratum-adjusted correlation (20 kb - 2 Mb):")
    print(pd.DataFrame([[scc(maps[a], maps[b]) if a != b else 1.0 for b in names] for a in names],
                       index=names, columns=names).round(3).to_string())
    print("\nlog2(obs/exp), 3x3, at Dugar's Rok cluster pairs:")
    oe = {k: obs_exp_log2(m) for k, m in maps.items()}
    print(pd.DataFrame({k: {f"{a}-{b}": round(o[RC[a], RC[b]], 2) for a, b in pairs} for k, o in oe.items()}).to_string())
    sys.exit()

if __name__ == "__main__":
    for cond, files in RUNS.items():
        summed = sum(raw(f) for f in files)
        tifffile.imwrite(OUT + f"Dugar_{cond}_st_merged_iced.tif", ice(summed).astype(np.float32))
        bins = cooler.Cooler(COOLS + files[0]).bins()[:][["chrom", "start", "end"]]
        iu = np.triu_indices(len(summed))
        keep = summed[iu] > 0
        pixels = pd.DataFrame({"bin1_id": iu[0][keep], "bin2_id": iu[1][keep], "count": summed[iu][keep].astype(int)})
        cooler.create_cooler(COOLS + f"Dugar_{cond}_st_merged_raw.cool", bins, pixels)
        print(f"{cond}: merged {len(files)} runs, {summed.sum() / 1e6:.1f} M contacts -> Dugar_{cond}_st_merged_iced.tif")
