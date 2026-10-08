"""Compare the ImageJ macro's loop calls with published loop callers.

Inputs: build_inputs.py (contact maps), run_callers.sh (tool runs), and the macro runs
(Comparisons/runs/bactoloops/ for merged WT 96 h, FIJI_output/ for the others).

For every method and dataset, loops are read as pairs of 0-based genome bins (10 kb) and scored
three ways, all on the same maps:
  1. APA: mean observed/expected in a 9 x 9-bin window around each loop, on the dataset's
     balanced .cool (circular distance-expected; filtered bins ignored). Summaries: log2 of the
     centre pixel, and P2LL (centre / mean of the 3 x 3 lower-left corner, the bins nearest the
     diagonal) as in Juicer's APA.
  2. Loop counts on WT 96 h (merged, and each replicate alone), the delta-stpA, delta-fis and delta-ihf 96 h maps (each rep1, with
     less sequencing depth, i.e. more sampling noise), WT exponential and delta-hns 96 h.
  3. Mean ChIP-seq z-score of the two anchors (H-NS 96 h for E. coli -- the WT 96 h track also
     stands in for the mutant 96 h maps -- and Rok for B. subtilis), with a bootstrap 95% CI
     (2,000 resamples of the loops), against a distance-matched null: every loop set rotated
     around the circular chromosome by every offset (both anchors shifted together; rotation_test.py).
The macro is scored as three call sets: all calls; main pass only (no ridge-search loops);
medium + high confidence.

Because mean z falls as a method calls more (weaker) loops, each method's loops are also
ranked by the method's own strength score and the mean z of its top k loops computed for every
k (matched-number comparison): the macro by integrated_excess, cooltools and HiCCUPS by
observed / donut expected, Chromosight by its correlation score, Mustache by FDR, SIP by its
APA score, Peakachu by its random-forest probability.

Writes Comparisons/method_comparison.csv, method_comparison_topk.csv and
method_comparison_apa.npz.
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
import sys as _sys
_sys.path.insert(0, HNS_ROOT + "Bactoloops/")
import rotation_test as RT   # exact rotation test: every offset once
import glob
import numpy as np
import pandas as pd
import cooler

ROOT = HNS_ROOT + "Bactoloops/"
CMP = ROOT + "Comparisons/"
RUNS = CMP + "runs/"
BIN = 10_000
HALF = 4                                     # 9 x 9 APA window, as in the macro
DATASETS = ["wt96", "wt96r1", "wt96r2", "stpA96", "fis96", "ihf96", "wtexp", "hns96", "bsub"]
MACRO = {"wt96": RUNS + "bactoloops/loop_calls_wt_96h_merged.csv",
         "wt96r1": ROOT + "FIJI_output/loop_calls_wt_96h_rep1.csv",
         "wt96r2": ROOT + "FIJI_output/loop_calls_wt_96h_rep2.csv",
         "stpA96": ROOT + "FIJI_output/loop_calls_stpA_96h.csv",
         "fis96": ROOT + "FIJI_output/loop_calls_fis_96h.csv",
         "ihf96": ROOT + "FIJI_output/loop_calls_IHF_96h.csv",
         "wtexp": ROOT + "FIJI_output/loop_calls_wt_Exp.csv",
         "hns96": ROOT + "FIJI_output/loop_calls_HNSm_96h.csv",
         "bsub": ROOT + "FIJI_output/loop_calls_bsub_WT_st_merged.csv"}
# The mutant 96 h maps have no ChIP of their own; the WT 96 h H-NS track stands in for them.
HNS96 = ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv"
CHIP = {"wt96": HNS96, "wt96r1": HNS96, "wt96r2": HNS96, "stpA96": HNS96, "fis96": HNS96, "ihf96": HNS96,
        "bsub": ROOT + "Dugar et al data/rok_zscore_10kb.csv"}
METHODS = ["Bactoloops (all)", "Bactoloops (main pass)", "Bactoloops (medium+)",
           "cooltools", "Chromosight", "HiCCUPS", "Mustache", "SIP", "Peakachu"]
rng = np.random.default_rng(0)


def pairs(df, c1, c2, scale=BIN):
    if df is None or len(df) == 0:
        return np.zeros((0, 2), int)
    a = np.round(df[c1].astype(float).values / scale).astype(int)
    b = np.round(df[c2].astype(float).values / scale).astype(int)
    return np.sort(np.c_[a, b], axis=1)


def read_tsv(path, **kw):
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None


def load(method, d):
    """(pairs of genome bins, strength score -- higher = stronger) for one method and dataset"""
    none = (np.zeros((0, 2), int), np.array([]))
    if method.startswith("Bactoloops"):
        try:
            L = pd.read_csv(MACRO[d])
        except (FileNotFoundError, pd.errors.EmptyDataError):
            return none
        if len(L) == 0:
            return none
        if "main pass" in method:
            L = L[~L.note.astype(str).str.startswith(("anchor-guided", "ridge region"))]
        if "medium+" in method:
            L = L[L.confidence.isin(["high", "medium"])]
        return pairs(L, "genome_bin1", "genome_bin2", scale=1), L.integrated_excess.values
    if method == "cooltools":
        T = read_tsv(RUNS + f"cooltools/{d}.tsv")
        return (pairs(T, "start1", "start2"), (T["count"] / T["la_exp.donut.value"]).values) if T is not None else none
    if method == "Chromosight":
        T = read_tsv(RUNS + f"chromosight/{d}.tsv")
        return (pairs(T, "start1", "start2"), T.score.values) if T is not None else none
    if method == "Mustache":
        T = read_tsv(RUNS + f"mustache/{d}.tsv")
        return (pairs(T, "BIN1_START", "BIN2_START"), -T.FDR.values) if T is not None and len(T) else none
    if method == "Peakachu":
        T = read_tsv(RUNS + f"peakachu/{d}.tsv", header=None)
        return (pairs(T, 1, 4), T[6].values) if T is not None else none
    if method == "HiCCUPS":
        T = read_tsv(RUNS + f"hiccups/{d}/merged_loops.bedpe", comment="#", header=None)
        return (pairs(T, 1, 4), (T[11] / T[13]).values) if T is not None else none
    if method == "SIP":
        f = glob.glob(RUNS + f"sip/{d}/*Loops*.txt")
        T = read_tsv(f[0]) if f else None
        return (pairs(T, "x1", "y1"), T.APScoreAvg.values) if T is not None and len(T) else none


def oe_map(d):
    m = cooler.Cooler(CMP + f"inputs/{d}.cool").matrix(balance=True)[:]
    n = len(m)
    i, j = np.indices(m.shape)
    dist = np.minimum(abs(i - j), n - abs(i - j))
    ok = np.isfinite(m)
    exp = np.bincount(dist[ok], m[ok], minlength=n // 2 + 1) / np.maximum(np.bincount(dist[ok], minlength=n // 2 + 1), 1)
    return m / np.where(exp[dist] > 0, exp[dist], np.nan)


def apa(OE, P):
    n = len(OE)
    stack = []
    for a, b in P:
        r = np.arange(a - HALF, a + HALF + 1) % n
        c = np.arange(b - HALF, b + HALF + 1) % n
        stack.append(OE[np.ix_(r, c)])
    if not stack:
        return None, np.nan, np.nan
    A = np.nanmean(np.array(stack), axis=0)
    centre = A[HALF, HALF]
    ll = np.nanmean(A[-3:, :3])                  # rows nearest the diagonal (a < b), left columns
    return A, np.log2(centre), centre / ll


def chip_track(d, n):
    ch = pd.read_csv(CHIP[d])
    z = np.full(n, np.nan)
    z[(ch.start_bp // BIN).astype(int).values % n] = ch.zscore.values
    return z


rows, apas, topk = [], {}, []
for d in DATASETS:
    OE = oe_map(d)
    n = len(OE)
    z = chip_track(d, n) if d in CHIP else None
    for meth in METHODS:
        P, score = load(meth, d)
        dist = np.minimum(P[:, 1] - P[:, 0], n - (P[:, 1] - P[:, 0])) if len(P) else np.array([])
        r = dict(dataset=d, method=meth, n_loops=len(P),
                 median_distance_kb=float(np.median(dist) * 10) if len(P) else np.nan)
        A, c, p2ll = apa(OE, P)
        r.update(apa_centre_log2=c, apa_p2ll=p2ll)
        if A is not None:
            apas[f"{d}|{meth}"] = A
        if z is not None and len(P):
            s = np.nanmean(z[P], axis=1)
            r["chip_mean_z"] = np.nanmean(s)
            r["chip_median_z"] = np.nanmedian(s)
            r["chip_pct_both_z_gt1"] = 100 * np.mean((z[P[:, 0]] > 1) & (z[P[:, 1]] > 1))
            null = [np.nanmean(np.nanmean(z[(P + k) % n], axis=1)) for k in RT.shifts(n)]
            r["chip_null_mean_z"] = np.mean(null)
            r["chip_null_95th"] = np.percentile(null, 95)
            r["chip_p_rotation"] = RT.p_value(r["chip_mean_z"], null, two_sided=False)
            sf = s[np.isfinite(s)]
            boot = [np.mean(rng.choice(sf, len(sf))) for _ in range(2000)]
            r["chip_mean_z_ci_lo"], r["chip_mean_z_ci_hi"] = np.percentile(boot, [2.5, 97.5])
            ranked = s[np.argsort(-score, kind="stable")]
            for k in range(1, len(ranked) + 1):
                topk.append(dict(dataset=d, method=meth, k=k, mean_z=np.nanmean(ranked[:k])))
        rows.append(r)

T = pd.DataFrame(rows)
T.to_csv(CMP + "method_comparison.csv", index=False)
K = pd.DataFrame(topk)
K.to_csv(CMP + "method_comparison_topk.csv", index=False)
np.savez(CMP + "method_comparison_apa.npz", **apas)
pd.set_option("display.width", 250)
print(T.round(3).to_string(index=False))

print("\nMean ChIP z of each method's top k loops (by its own score):")
for d in CHIP:
    W = K[K.dataset == d].pivot(index="k", columns="method", values="mean_z")
    ks = [k for k in [5, 8, 10, 12, 15, 17, 25, 35, 41, 59, 79, 124, 148] if k in W.index]
    print(d); print(W.loc[ks, [m for m in METHODS if m in W.columns]].round(2).to_string())
