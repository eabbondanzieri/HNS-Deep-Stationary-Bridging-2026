"""Loop-call numbers quoted in the manuscript text (Results, Methods), from the Bactoloops outputs.
Run on the current FIJI_output (default) or on an archived run folder (argument) to check that the
definitions reproduce earlier numbers. Loops match when both ends lie within 2 bins (circularly, either
orientation); anchors are loop ends merged across adjacent bins; H-NS track = anti-H-NS 96 h z-score (10 kb).
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
import sys, numpy as np, pandas as pd
sys.path.insert(0, HNS_ROOT + "Bactoloops/")
import rotation_test as RT
ROOT = HNS_ROOT + "Bactoloops/"
D = sys.argv[1] if len(sys.argv) > 1 else ROOT + "FIJI_output/"
N = 461
z = np.full(N, np.nan); Z = pd.read_csv(ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv"); z[(Z.start_bp // 10000).values] = Z.zscore.values
def rd(name):
    try: L = pd.read_csv(D + f"loop_calls_{name}.csv")
    except (FileNotFoundError, pd.errors.EmptyDataError): return np.zeros((0, 2), int)
    return np.c_[L.genome_bin1.round().astype(int), L.genome_bin2.round().astype(int)] % N
cd = lambda a, b: np.minimum(np.abs(a - b) % N, N - np.abs(a - b) % N)
def matched(A, B, tol=2):
    if len(B) == 0: return np.zeros(len(A), bool)
    return np.array([(((cd(B[:, 0], a) <= tol) & (cd(B[:, 1], b) <= tol)) | ((cd(B[:, 0], b) <= tol) & (cd(B[:, 1], a) <= tol))).any() for a, b in A])
def anchors(P):
    ends = np.zeros(N, bool); ends[P.ravel()] = True
    regs, i0 = [], np.where(~ends)[0][0]
    cur = []
    for k in range(1, N + 1):
        b = (i0 + k) % N
        if ends[b]: cur.append(b)
        elif cur: regs.append(cur); cur = []
    if cur: regs.append(cur)
    return regs
def peaks(z):
    zz = np.nan_to_num(z, nan=-9)
    return np.array([i for i in range(N) if zz[i] >= 2 and zz[i] == max(zz[(i + k) % N] for k in range(-2, 3))])
W = rd("wt_96h_repavg"); A = anchors(W)
reg_z = [np.nanmax(z[[(b + k) % N for b in r for k in (-1, 0, 1)]]) for r in A]
pk = peaks(z); abins = np.array(sorted(b for r in A for b in r))
pk_dist = np.array([cd(abins, p).min() for p in pk])
print(f"WT 96 h: {len(W)} loops, {len(A)} anchor regions; anchors with z > 1 within 10 kb: {sum(v > 1 for v in reg_z)}/{len(A)}")
print(f"  strong H-NS peaks (z >= 2, local max +-2 bins): {len(pk)}; at an anchor (within 1 bin): {(pk_dist <= 1).sum()}; others at {sorted(10 * pk_dist[pk_dist > 1])} kb")
span = np.median(cd(W[:, 0], W[:, 1])) * 10
print(f"  median loop span {span / 1000:.2f} Mb")
for t, name in (("24 h", "wt_24h"), ("exp", "wt_Exp")):
    P = rd(name)
    if not len(P): print(t, "no loops"); continue
    m = matched(P, W).sum()
    null = [matched((P + s) % N, W).sum() for s in RT.shifts(N)]
    ends_at = sum(cd(abins, e).min() <= 1 for e in P.ravel())
    print(f"{t}: {len(P)} loops; matched a 96 h loop: {m} (rotation mean {np.mean(null):.2f}, P = {RT.p_value(m, null, two_sided=False):.3f} one-sided); "
          f"ends within 1 bin of a 96 h anchor: {ends_at}/{2 * len(P)}")
r1, r2 = rd("wt_96h_rep1"), rd("wt_96h_rep2")
print(f"replicates: {len(r1)} and {len(r2)} loops; rep1 loops matched in rep2: {matched(r1, r2).sum()}, rep2 in rep1: {matched(r2, r1).sum()}; "
      f"matched averaged map: {matched(r1, W).sum()} and {matched(r2, W).sum()}")
# strongly H-NS-bound regions as contiguous runs of bins with z >= 2 (circular)
hi = np.nan_to_num(z, nan=-9) >= 2
runs, cur, s0 = [], [], np.where(~hi)[0][0]
for k in range(1, N + 1):
    b = (s0 + k) % N
    if hi[b]: cur.append(b)
    elif cur: runs.append(cur); cur = []
if cur: runs.append(cur)
rd_ = [min(cd(abins, b).min() for b in r) for r in runs]
print(f"strong regions (runs of z >= 2): {len(runs)}; overlapping/adjacent to an anchor (<= 1 bin): {sum(d <= 1 for d in rd_)}; "
      f"exact overlap: {sum(d == 0 for d in rd_)}; others at {sorted(10 * d for d in rd_ if d > 1)} kb; distances of non-overlapping: {sorted(10 * d for d in rd_ if d > 0)} kb")
