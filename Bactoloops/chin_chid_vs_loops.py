"""Exponential-phase CHINs and CHIDs (Gavrilov, Shamovsky et al. 2025; Micro-C, MG1655 LB) vs deep-stationary loops.
Elements: Nudler CHIN paper/opcidsRenamed.xlsx, lab-genome coordinates (typos fixed as in chin_analysis.py).
Density: bp covered by each element type in every 10-kb bin. Pearson r with the WT 96 h ILD and anti-H-NS ChIP z (96 h) and
the anti-mCherry exp ChIP z, exact circular-shift P (shifts >= 20 bins). Anchors: WT 96 h loop ends merged across adjacent
bins; overlap = element within one bin of an anchor, rotation null over every offset (rotation_test.shifts).
Output: chin_chid_vs_loops.csv"""
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
B = HNS_ROOT + "Bactoloops/"
sys.path.insert(0, B); import rotation_test as RT
N, BIN = 461, 10000
X = pd.read_excel(B + "Nudler CHIN paper/opcidsRenamed.xlsx", 0)
X.loc[X.OPCID_ID == "OPCID_40", "start_new"] = 2_522_632; X.loc[X.OPCID_ID == "CHID_20", "end_new"] = 3_950_261
X = X.dropna(subset=["start_new", "end_new"]); X["type"] = X.OPCID_ID.str.split("_").str[0]
X["s"], X["e"] = X[["start_new", "end_new"]].min(1).astype(int), X[["start_new", "end_new"]].max(1).astype(int)
def density(df):
    v = np.zeros(N)
    for s, e in zip(df.s, df.e):
        for b in range(s // BIN, e // BIN + 1):
            v[b % N] += max(0, min(e, (b + 1) * BIN) - max(s, b * BIN))
    return v
D = {t: density(X[X.type == t]) for t in ("CHIN", "CHID")}
def track(f):
    c = pd.read_csv(B + "ChIP-seq_data/" + f); v = np.zeros(N); v[(c.start_bp // BIN).astype(int) % N] = np.nan_to_num(c.zscore.values); return v
t = pd.read_csv(B + "FIJI_output/integrated_loop_density_wt_96h_repavg.csv").ILD_au.values; ild = np.empty(N); ild[(367 + np.arange(N)) % N] = t
TR = {"ILD WT 96 h": ild, "anti-H-NS 96 h": track("zscore_wt_96h_10kb.csv"), "anti-mCherry exp": track("zscore_mCherry_Exp_10kb.csv")}
def shift(a, b):
    r0 = np.corrcoef(a, b)[0, 1]; null = [np.corrcoef(a, np.roll(b, k))[0, 1] for k in RT.shifts(N, exclude=20)]
    return r0, (1 + np.sum(np.abs(null) >= abs(r0))) / (len(null) + 1)
rows = []
for t_, d in D.items():
    for k, v in TR.items():
        r, p = shift(d, v); rows.append(dict(element=t_, n=int((X.type == t_).sum()), track=k, r=r, p=p))
L = pd.read_csv(B + "FIJI_output/loop_calls_wt_96h_repavg.csv")
ends = np.unique(np.r_[L.genome_bin1.round(), L.genome_bin2.round()].astype(int) % N)
near = np.zeros(N, bool)
for e in ends:
    for k in (-1, 0, 1): near[(e + k) % N] = True
for t_ in ("CHIN", "CHID"):
    df = X[X.type == t_]; mid = ((df.s + df.e) // 2 // BIN).values % N
    obs = near[mid].mean(); null = np.array([near[(mid + k) % N].mean() for k in RT.shifts(N)])
    rows.append(dict(element=t_, n=len(df), track="fraction within 10 kb of a 96 h anchor", r=obs, p=RT.p_value(obs, null), null_mean=null.mean()))
    # anchors (merged) containing an element
    anc = []; cur = [ends[0]]
    for b in ends[1:]:
        (cur.append(b) if b - cur[-1] <= 1 else (anc.append(cur), cur := [b]))
    anc.append(cur)
    if len(anc) > 1 and anc[0][0] + N - anc[-1][-1] <= 1: anc[0] = anc.pop() + anc[0]
    has = np.mean([any(cd <= 1 for cd in [min((m - a) % N, (a - m) % N) for m in mid for a in g]) for g in anc])
    rows.append(dict(element=t_, n=len(anc), track="fraction of 96 h anchors with an element within 10 kb", r=has, p=np.nan))
T = pd.DataFrame(rows); T.to_csv(B + "chin_chid_vs_loops.csv", index=False); print(T.round(3).to_string(index=False))
