"""Compare loop-network architecture in E. coli (this study) vs B. subtilis (Dugar et al.).

Loop strength is measured the same way in both species, as log2(observed/expected)
from the ICE-balanced matrix, where expected is the mean contact at the same
circular genomic distance -- so it doesn't depend on sequencing depth or on the
macro's intensity scaling.

Anchors are the loop-call end points, merged when within ANCHOR_GAP bins of each
other. Every pair of anchors (called as a loop or not) is then scored, which is
what lets us ask whether sub-threshold anchor pairs still interact.
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
from itertools import combinations
import numpy as np
import pandas as pd
import tifffile
from scipy import stats

ROOT = HNS_ROOT + "Bactoloops/"
DATASETS = {
    # name: (loop calls, ICE matrix in the SAME (tiff) coordinates as bin1/bin2,
    #        oriC genome bin, recenterBin, ChIP z-score track for the anchor protein)
    "E. coli": (ROOT + "FIJI_output/loop_calls_wt_96h_rep1.csv",
                ROOT + "iced_matrix_conversion/cWX2200_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif",
                367, 368,    # EcW3110del oriC at 3,671,577 bp (between mioC and gidA)
                ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv"),       # anti-H-NS
    "B. subtilis": (ROOT + "FIJI_output/loop_calls_rok_st.csv",
                    ROOT + "iced_matrix_conversion/Dugar_HiC/GSM4288936_WT_st_HiC_iced.tif",
                    0, 0,    # B. subtilis 168 oriC at 0
                    ROOT + "Dugar et al data/rok_zscore_10kb.csv"),  # anti-Rok
}
ANCHOR_GAP = 2      # bins; loop ends this close are the same anchor
MIN_DIST = 10       # bins; ignore anchor pairs closer than 100 kb (diagonal/ridge)
WIN = 1             # half-width of the window averaged around each pixel (3x3)
N_NULL = 200        # distance-matched random pixels per anchor pair
MATCH_TOL = 5       # bins; a weak loop end "lands on" a strong-loop anchor if within this
N_PERM = 2000
rng = np.random.default_rng(0)


def obs_exp(m):
    n = m.shape[0]
    i, j = np.indices(m.shape)
    d = np.minimum(abs(i - j), n - abs(i - j))
    ok = m > 0
    e = np.bincount(d[ok], m[ok], minlength=n // 2 + 1) / np.maximum(
        np.bincount(d[ok], minlength=n // 2 + 1), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ok, m / e[d], np.nan)


def score(oe, r, c):
    """mean log2(O/E) over a (2*WIN+1)^2 window, wrapping circularly"""
    n = oe.shape[0]
    rr = np.arange(r - WIN, r + WIN + 1) % n
    cc = np.arange(c - WIN, c + WIN + 1) % n
    return np.nanmean(np.log2(oe[np.ix_(rr, cc)]))


def cdist(a, b, n):
    d = abs(a - b) % n
    return min(d, n - d)


def gini(x):
    x = np.sort(np.asarray(x, float))
    x = x - min(0, x.min())
    n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum())


def cluster_anchors(bins, n):
    """single-linkage merge of loop-end bins on a circle"""
    u = np.sort(np.unique(bins))
    groups, cur = [], [u[0]]
    for b in u[1:]:
        if b - cur[-1] <= ANCHOR_GAP:
            cur.append(b)
        else:
            groups.append(cur); cur = [b]
    groups.append(cur)
    if len(groups) > 1 and (groups[0][0] + n) - groups[-1][-1] <= ANCHOR_GAP:
        groups[0] = groups[-1] + groups[0]; groups.pop()
    return groups


def strength_tiers(x):
    """split loops at the largest ratio between consecutive sorted strengths
    (searched within the top half, so the tail's noise can't win)"""
    s = np.sort(x)[::-1]
    ratio = s[:-1] / s[1:]
    k = int(np.argmax(ratio[: max(1, len(s) // 2)]))
    return s[k], ratio[k], k + 1        # cutoff value, gap size, n strong


def analyze(name, loops_csv, tif, ori_genome, recenter, chip_csv):
    L = pd.read_csv(loops_csv)
    m = tifffile.imread(tif).astype(float)
    n = m.shape[0]
    oe = obs_exp(m)
    b1 = L.bin1.round().astype(int).values
    b2 = L.bin2.round().astype(int).values
    L["oe"] = [score(oe, r, c) for r, c in zip(b1, b2)]

    groups = cluster_anchors(np.r_[b1, b2], n)
    bin2anchor = {b: k for k, g in enumerate(groups) for b in g}
    # representative bin of each anchor: the loop end used most often
    ends = pd.Series(np.r_[b1, b2])
    rep = [int(ends[ends.isin(g)].mode().iloc[0]) for g in groups]
    L["a1"] = [bin2anchor[b] for b in b1]
    L["a2"] = [bin2anchor[b] for b in b2]
    called = {tuple(sorted(p)) for p in zip(L.a1, L.a2)}
    A = len(groups)

    # every anchor pair, called or not, plus a distance-matched null for each
    allb = np.concatenate(groups)
    near_anchor = np.zeros(n, bool)
    for b in allb:
        near_anchor[np.arange(b - 3, b + 4) % n] = True
    rows = []
    for i, j in combinations(range(A), 2):
        d = cdist(rep[i], rep[j], n)
        if d < MIN_DIST:
            continue
        s = score(oe, rep[i], rep[j])
        null = []
        while len(null) < N_NULL:
            r = rng.integers(n)
            c = (r + d) % n
            if near_anchor[r] or near_anchor[c]:
                continue
            null.append(score(oe, r, c))
        null = np.array(null)
        rows.append(dict(i=i, j=j, dist=d, oe=s, called=(i, j) in called,
                         null_med=np.nanmedian(null),
                         pct=np.nanmean(null < s)))
    P = pd.DataFrame(rows)

    # additive "anchor stickiness" model: oe_ij ~ mu + a_i + a_j
    X = np.zeros((len(P), A))
    X[np.arange(len(P)), P.i] = 1
    X[np.arange(len(P)), P.j] = 1
    X = np.c_[np.ones(len(P)), X]
    beta, *_ = np.linalg.lstsq(X, P.oe.values, rcond=None)
    resid = P.oe.values - X @ beta
    r2 = 1 - resid.var() / P.oe.var()
    k = X.shape[1] - 1
    adj = lambda r: 1 - (1 - r) * (len(P) - 1) / (len(P) - k - 1)
    r2_perm = []
    for _ in range(N_PERM // 4):
        y = rng.permutation(P.oe.values)
        bb, *_ = np.linalg.lstsq(X, y, rcond=None)
        r2_perm.append(1 - (y - X @ bb).var() / y.var())
    p_r2 = (1 + np.sum(np.array(r2_perm) >= r2)) / (1 + len(r2_perm))

    # graph metrics on called loops
    edges = np.array(sorted(called))
    deg = np.bincount(edges.ravel(), minlength=A)
    density = len(called) / (A * (A - 1) / 2)
    centralization = ((deg.max() - deg).sum()) / ((A - 1) * (A - 2)) if A > 2 else np.nan

    # position relative to oriC (genome coordinates)
    gpos = lambda b: (b if recenter == 0 else (recenter - 1 + b) % n)
    ori_d = lambda b: cdist(gpos(b), ori_genome, n) / n   # fraction of genome, 0..0.5
    L["ori_dist"] = [(ori_d(r) + ori_d(c)) / 2 for r, c in zip(b1, b2)]
    ex = np.clip(2 ** L.oe - 1, 0, None)                # excess contact over expected
    near = L.ori_dist < 0.15                             # both-ish anchors within ~15% of genome of oriC
    frac_ori = ex[near].sum() / ex.sum()
    # null: rotate the whole loop set around the chromosome
    rot = []
    for s in range(n):
        od = [(cdist(gpos(r) + s, ori_genome, n) + cdist(gpos(c) + s, ori_genome, n)) / (2 * n)
              for r, c in zip(b1, b2)]
        rot.append(ex[np.array(od) < 0.15].sum() / ex.sum())
    p_ori = np.mean(np.array(rot) >= frac_ori)

    # strong vs weak tier, from the macro's own loop_strength. Anchor-guided
    # loops (n_detections == 0) have no particle, so no loop_strength: they
    # are left out of the tiers, and of the weak-loop tests below.
    main = (L.n_detections > 0).values
    cut, gap, n_strong = strength_tiers(L.loop_strength.values[main])
    strong = (L.loop_strength >= cut).values & main
    L["tier"] = np.where(strong, "strong", np.where(main, "weak", "anchor-guided"))
    weak = main & ~strong
    sbins = np.r_[b1[strong], b2[strong]]
    on_strong = lambda b: np.any([cdist(b, x, n) <= MATCH_TOL for x in sbins])
    weak_between_strong = np.mean([on_strong(r) and on_strong(c)
                                   for r, c in zip(b1[weak], b2[weak])])
    # null: same weak loops, rotated to a random position around the chromosome
    wb = [(r, c) for r, c in zip(b1[weak], b2[weak])]
    rot_ws = [np.mean([on_strong((r + sft) % n) and on_strong((c + sft) % n) for r, c in wb])
              for sft in range(1, n)]
    p_ws = (1 + np.sum(np.array(rot_ws) >= weak_between_strong)) / n
    weak_touch_strong = np.mean([on_strong(r) or on_strong(c) for r, c in wb])
    rot_wt = [np.mean([on_strong((r + sft) % n) or on_strong((c + sft) % n) for r, c in wb])
              for sft in range(1, n)]
    p_wt = (1 + np.sum(np.array(rot_wt) >= weak_touch_strong)) / n

    # anchor-protein ChIP enrichment at loop anchors (genome coords; track bin 1 = genome bin 0)
    # indexed by genome bin (start_bp // 10 kb), not row number: a ChIP track can
    # skip a bin (the Rok track has no 480-490 kb row), leaving NaN there
    chip = pd.read_csv(chip_csv)
    z = np.full(n, np.nan)
    z[(chip.start_bp // 10_000).astype(int).values % n] = chip.zscore.values
    zs = lambda bins: np.nanmean([z[gpos(b) % len(z)] for b in bins])
    z_strong = zs(np.unique(sbins))
    z_all = zs(np.unique(np.r_[b1, b2]))
    z_pct = np.mean([np.nanmean(rng.choice(z[np.isfinite(z)], len(np.unique(np.r_[b1, b2])), replace=False)) < z_all
                     for _ in range(N_PERM)])

    unc = P[~P.called]
    res = {
        "loops called": len(L),
        "anchors": A,
        "anchor pairs scored (>=100 kb)": len(P),
        "strong-tier loops (largest gap in loop_strength)": n_strong,
        "  gap between tiers (fold)": gap,
        "called loop log2(O/E), median": L.oe.median(),
        "called loop log2(O/E), IQR": f"{L.oe.quantile(.25):.2f}-{L.oe.quantile(.75):.2f}",
        "CV of called-loop excess (O/E-1)": ex.std() / ex.mean(),
        "Gini of called-loop excess": gini(ex),
        "share of excess in top 5 loops": np.sort(ex)[::-1][:5].sum() / ex.sum(),
        "Gini of macro loop_strength": gini(L.loop_strength.values[main]),
        "network density (called / possible anchor pairs)": density,
        "degree centralization (1 = star)": centralization,
        "max anchor degree": int(deg.max()),
        "uncalled anchor pairs, median log2(O/E)": unc.oe.median(),
        "uncalled anchor pairs, median null log2(O/E)": unc.null_med.median(),
        "uncalled anchor pairs above null 95th pct": np.mean(unc.pct > 0.95),
        "uncalled vs null, Wilcoxon p": stats.wilcoxon(unc.oe - unc.null_med).pvalue,
        "additive anchor model R^2 (all pairs)": r2,
        "  adjusted R^2": adj(r2),
        "  p (permuted pair scores)": p_r2,
        "weak loops with >=1 end on a strong-loop anchor": weak_touch_strong,
        "  p (random rotation of weak loops) ": p_wt,
        "weak loops joining two strong-loop anchors": weak_between_strong,
        "  p (random rotation of weak loops)": p_ws,
        "anchor ChIP z-score, all anchors": z_all,
        "  percentile vs random bins": z_pct,
        "anchor ChIP z-score, strong-tier anchors": z_strong,
        "excess within 15% of genome of oriC": frac_ori,
        "  p (random rotation of loop set)": p_ori,
    }
    return res, L, P, groups, rep, oe


def anchor_metrics(P, keep):
    """adjusted R^2 of the additive model and fraction of uncalled pairs above
    the null 95th pct, restricted to pairs among the anchors in `keep`"""
    Q = P[P.i.isin(keep) & P.j.isin(keep)]
    ids = {a: k for k, a in enumerate(sorted(keep))}
    X = np.zeros((len(Q), len(ids)))
    X[np.arange(len(Q)), Q.i.map(ids)] = 1
    X[np.arange(len(Q)), Q.j.map(ids)] = 1
    X = np.c_[np.ones(len(Q)), X]
    b, *_ = np.linalg.lstsq(X, Q.oe.values, rcond=None)
    r2 = 1 - (Q.oe.values - X @ b).var() / Q.oe.var()
    adj = 1 - (1 - r2) * (len(Q) - 1) / (len(Q) - len(ids) - 1)
    unc = Q[~Q.called]
    return adj, np.mean(unc.pct > 0.95)


def matched_comparison(big, small, n_iter=1000):
    """subsample the larger network (anchors, and loops for Gini) down to the
    smaller one's size, so size differences can't drive the comparison"""
    (_, Lb, Pb, gb, _, _), (_, Ls, Ps, gs, _, _) = big, small
    k = len(gs)
    sub = np.array([anchor_metrics(Pb, set(rng.choice(len(gb), k, replace=False)))
                    for _ in range(n_iter)])
    ex = lambda L: np.clip(2 ** L.oe.values - 1, 0, None)
    g_sub = [gini(rng.choice(ex(Lb), len(Ls), replace=False)) for _ in range(n_iter)]
    small_adj, small_frac = anchor_metrics(Ps, set(range(len(gs))))
    return {
        "anchors matched to": k,
        "adj R^2, subsampled big network (2.5-97.5%)": np.percentile(sub[:, 0], [2.5, 97.5]).round(2).tolist(),
        "adj R^2, small network": round(small_adj, 2),
        "  fraction of subsamples <= small": np.mean(sub[:, 0] <= small_adj),
        "uncalled above null, subsampled big (2.5-97.5%)": np.percentile(sub[:, 1], [2.5, 97.5]).round(2).tolist(),
        "uncalled above null, small network": round(small_frac, 2),
        "  fraction of subsamples <= small ": np.mean(sub[:, 1] <= small_frac),
        "loops matched to": len(Ls),
        "Gini, subsampled big (2.5-97.5%)": np.percentile(g_sub, [2.5, 97.5]).round(2).tolist(),
        "Gini, small network": round(gini(ex(Ls)), 2),
        "  fraction of subsamples >= small": np.mean(np.array(g_sub) >= gini(ex(Ls))),
    }


if __name__ == "__main__":
    out, full = {}, {}
    for name, args in DATASETS.items():
        full[name] = analyze(name, *args)
        res, L, P, groups, rep, oe = full[name]
        out[name] = res
        tag = name.replace(". ", "_").replace(" ", "_")
        L.to_csv(f"{tag}_loops_scored.csv", index=False)
        P.to_csv(f"{tag}_anchor_pairs.csv", index=False)
        print(name, "anchors (bins):", [f"{g[0]}-{g[-1]}" if len(g) > 1 else str(g[0]) for g in groups])
    T = pd.DataFrame(out)
    pd.set_option("display.width", 200)
    print(T.to_string())
    T.to_csv("loop_network_summary.csv")
    M = pd.Series(matched_comparison(full["E. coli"], full["B. subtilis"]))
    print("\nSize-matched comparison (E. coli subsampled to B. subtilis size):")
    print(M.to_string())
    M.to_csv("size_matched_comparison.csv")
