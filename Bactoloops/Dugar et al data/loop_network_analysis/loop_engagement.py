"""Loop strength without O/E normalization, and tests of "one hub" vs "stochastic pairs".

O/E divides by the expected contact frequency at each distance. The E. coli map
decays far more gently with distance than the B. subtilis map (at 2 Mb it is
still ~13% of its 10 kb contact level vs ~5%), so the same number of extra
contacts gives a smaller O/E in E. coli. Here loop strength is instead the
EXCESS contact (observed - expected), summed over a window around the loop and
expressed as a fraction of the anchor's non-local contacts (all contacts >= 2
bins away; ICE makes these near-identical for every bin). That is directly
comparable between datasets and adds up across partners.

1. Scalogram: each loop's excess fraction summed over windows of half-width
   r = 0..6 bins -- shows whether loops differ in peak height, width, or both.
2. Anchor engagement: for each anchor, the excess fraction summed over ALL
   other anchors (called or not), i.e. how much of that anchor's contact
   budget goes to the anchor network.
3. Budget test: if each anchor has a fixed contact budget shared among
   partners (stochastic/competitive pairing), per-partner strength falls as
   the number of partners rises; if anchors meet as one complex, per-partner
   strength doesn't depend on how many partners there are.
4. Distance test: in a single co-localized complex, anchor-anchor contact
   frequency shouldn't fall with genomic distance the way background contacts
   do; for encounter-driven pairing it should fall at roughly the background
   rate.
5. Contact level: the genomic separation at which background contact
   frequency equals the observed anchor-anchor contact ("equivalent
   distance"). A complex present in every cell should look like near
   neighbors (tens of kb).
"""
from itertools import combinations
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
import loop_network_stats as S

R_MAX = 6          # scalogram half-widths 0..R_MAX bins
R_SUM = 2          # window half-width (5x5) used for per-loop / per-anchor totals
MIN_D = 10         # bins; only pairs >= 100 kb apart (off the diagonal ridge)


def expected_and_excess(m):
    n = len(m)
    i, j = np.indices(m.shape)
    d = np.minimum(abs(i - j), n - abs(i - j))
    ok = m > 0
    e = np.bincount(d[ok], m[ok], minlength=n // 2 + 1) / np.maximum(np.bincount(d[ok], minlength=n // 2 + 1), 1)
    E = e[d]
    X = np.where(ok, m - E, 0.0)                      # excess contacts; missing bins contribute 0
    nonlocal_sum = np.median(np.where(d >= 2, m, 0).sum(1)[m.sum(1) > 0])
    # expected, excess and observed, all as fractions of an anchor's non-local contacts
    return e / nonlocal_sum, X / nonlocal_sum, m / nonlocal_sum


def window_sum(A, r, c, h):
    n = len(A)
    rr = np.arange(r - h, r + h + 1) % n
    cc = np.arange(c - h, c + h + 1) % n
    return A[np.ix_(rr, cc)].sum()


def equivalent_distance(obs_frac, e_frac):
    """smallest separation (bins) at which the background falls to obs_frac"""
    below = np.where(e_frac[1:] <= obs_frac)[0]       # from 1 bin: the diagonal itself isn't a separation
    return below[0] + 1 if len(below) else np.nan


def main():
    rows_loop, rows_anchor, summary = [], [], {}
    for name, (loops_csv, tif, *_rest) in S.DATASETS.items():
        m = tifffile.imread(tif).astype(float)
        n = len(m)
        e_frac, Xf, Of = expected_and_excess(m)

        L = pd.read_csv(loops_csv)
        b1 = L.bin1.round().astype(int).values
        b2 = L.bin2.round().astype(int).values
        dist = np.array([S.cdist(a, b, n) for a, b in zip(b1, b2)])
        keep = dist >= MIN_D
        for a, b, dd in zip(b1[keep], b2[keep], dist[keep]):
            scal = [window_sum(Xf, a, b, h) for h in range(R_MAX + 1)]
            obs = Of[a, b]
            rows_loop.append(dict(species=name, bin1=a, bin2=b, dist_kb=dd * 10, peak_frac=scal[0],
                                  **{f"excess_r{h}": v for h, v in enumerate(scal)},
                                  obs_frac=obs, equiv_dist_kb=equivalent_distance(obs, e_frac) * 10))

        # anchors and every anchor pair (called or not)
        groups = S.cluster_anchors(np.r_[b1, b2], n)
        ends = pd.Series(np.r_[b1, b2])
        rep = [int(ends[ends.isin(g)].mode().iloc[0]) for g in groups]
        bin2a = {bb: k for k, g in enumerate(groups) for bb in g}
        called = {tuple(sorted((bin2a[a], bin2a[b]))) for a, b in zip(b1, b2)}
        A = len(rep)
        pair = []
        for i, j in combinations(range(A), 2):
            dd = S.cdist(rep[i], rep[j], n)
            if dd < MIN_D:
                continue
            pair.append(dict(i=i, j=j, dist=dd, called=(i, j) in called,
                             excess=window_sum(Xf, rep[i], rep[j], R_SUM),
                             obs=np.mean([Of[(rep[i] + u) % n, (rep[j] + v) % n] for u in (-1, 0, 1) for v in (-1, 0, 1)]),
                             exp=e_frac[dd]))
        P = pd.DataFrame(pair)
        deg = np.zeros(A)
        for i, j in called:
            deg[i] += 1
            deg[j] += 1
        for k in range(A):
            mine = P[(P.i == k) | (P.j == k)]
            rows_anchor.append(dict(species=name, anchor=rep[k], degree=int(deg[k]),
                                    total_excess_all_anchors=mine.excess.sum(),
                                    total_excess_called=mine[mine.called].excess.sum(),
                                    per_partner_called=mine[mine.called].excess.mean() if mine.called.any() else np.nan))

        # distance test: slope of log(observed) vs log(distance) for anchor pairs vs for all pixels
        ok = P.obs > 0
        slope_pairs = stats.linregress(np.log10(P.dist[ok] * 10), np.log10(P.obs[ok]))
        dd_all = np.arange(MIN_D, n // 2 + 1)
        slope_bg = stats.linregress(np.log10(dd_all * 10), np.log10(e_frac[dd_all]))
        slope_called = stats.linregress(np.log10(P.dist[P.called] * 10), np.log10(P.obs[P.called]))
        # additive (a fixed extra contact wherever the anchors are -> excess flat with distance)
        # vs multiplicative (a fixed fold-increase over background -> excess falls like background)
        cx = P.called & (P.excess > 0)
        slope_excess = stats.linregress(np.log10(P.dist[cx] * 10), np.log10(P.excess[cx]))
        Aan = pd.DataFrame([r for r in rows_anchor if r["species"] == name])
        Aan = Aan[Aan.degree > 0]
        rho = stats.spearmanr(Aan.degree, Aan.per_partner_called)
        summary[name] = {
            "loops scored (>=100 kb)": int(keep.sum()),
            "loop excess, 5x5 window, median (% of anchor contacts)": 100 * np.median([r["excess_r2"] for r in rows_loop if r["species"] == name]),
            "loop peak excess, 1 pixel, median (%)": 100 * np.median([r["excess_r0"] for r in rows_loop if r["species"] == name]),
            "anchor total excess to ALL anchors, median (%)": 100 * Aan.total_excess_all_anchors.median(),
            "anchor total excess to called partners, median (%)": 100 * Aan.total_excess_called.median(),
            "partners per anchor (called), median": Aan.degree.median(),
            "budget test: Spearman(degree, per-partner excess)": f"rho={rho.correlation:.2f}, p={rho.pvalue:.2g}",
            "distance slope, background (log-log)": round(slope_bg.slope, 2),
            "distance slope, all anchor pairs": f"{slope_pairs.slope:.2f} (SE {slope_pairs.stderr:.2f})",
            "distance slope, called loops": f"{slope_called.slope:.2f} (SE {slope_called.stderr:.2f})",
            "distance slope, called loops' EXCESS": f"{slope_excess.slope:.2f} (SE {slope_excess.stderr:.2f})",
            "equivalent distance of called loops, median (kb)": np.nanmedian([r["equiv_dist_kb"] for r in rows_loop if r["species"] == name]),
            "loop genomic distance, median (kb)": np.median([r["dist_kb"] for r in rows_loop if r["species"] == name]),
        }

    pd.DataFrame(rows_loop).to_csv("loop_engagement_loops.csv", index=False)
    pd.DataFrame(rows_anchor).to_csv("loop_engagement_anchors.csv", index=False)
    T = pd.DataFrame(summary)
    T.to_csv("loop_engagement_summary.csv")
    pd.set_option("display.width", 200)
    print(T.to_string())


if __name__ == "__main__":
    main()
