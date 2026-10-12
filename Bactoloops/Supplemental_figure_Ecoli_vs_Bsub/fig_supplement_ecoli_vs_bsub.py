"""Draft supplemental figure: E. coli H-NS loop network vs B. subtilis Rok loop network.

All panels use the merged datasets:
  E. coli WT 96 h   -- average of rep1 + rep2 ICE matrices
                       (iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif,
                        loops FIJI_output/loop_calls_wt_96h_repavg.csv), H-NS ChIP (wt 96 h)
  B. subtilis WT st -- Dugar et al. GEO runs 1 + 2, raw counts summed then ICE-balanced
                       (iced_matrix_conversion/Dugar_HiC/Dugar_WT_st_merged_iced.tif,
                        loops FIJI_output/loop_calls_bsub_WT_st_merged.csv), Rok ChIP

Loop strength is measured against the local background (localbg.py): the macro's
log2(obs/exp) map with ImageJ's rolling ball (radius 50 px = 500 kb) subtracted, on a
circularly padded map. This removes large-scale structure that distance-only expected leaves
in (in B. subtilis the SMC arm-alignment X and the oriC/terminus domains). "Enriched" = above
the 95th percentile of all pixels at the same genomic distance in the corrected map.
Integrated excess = excess contacts over the local background summed over a 50 x 50 kb window,
as % of the anchor's non-local contacts (approximately how often a loop forms), baseline-corrected (minus the map's
mean at the same distance; excess_baseline.py) so that the two maps are comparable.

Panels
  A  circular loop networks, arcs weighted by integrated excess relative to the dataset mean
     (CircularLoopArcs/circular_ecoli_vs_bsub.py --metric integrated_excess --scale relative
     --min-confidence medium --panel), embedded as an image
  B  strength tiers: medium+ confidence loops ranked by integrated excess relative to the dataset mean
  C  every pair of loop anchors, split into called loops and uncalled pairs: corrected
     log2(obs/exp) above distance-matched background
  D  mean integrated excess per loop and per anchor (all partners), relative to E. coli,
     with bootstrap 95% CIs
  E  every pair of ChIP-seq peaks (independent of loop calls): enrichment above background,
     and the fraction enriched across peak thresholds
  G  silencing (expression higher in the deletion: log2 mutant/WT) at silencer-bound regions (ChIP z >= 2, 10-kb bins) in exponential vs stationary phase:
     mean log2(mutant/WT) expression at bound bins minus the genome mean, as fold, with bootstrap 95% CI and
     exact rotation P; above, the change from exponential to stationary phase. B. subtilis RNA-seq: Dugar et
     al. (GSE144475, one replicate per condition); E. coli: this study (CFU-normalized). From
     Dugar et al data/rok_vs_hns_silencing.py
  H  the local effect reported by Dugar et al.: expression in delta-rok within +/-50 kb of each of their eight Rok
     clusters, relative to the genome mean, exponential vs stationary (22 h)

Run from anywhere:  python fig_supplement_ecoli_vs_bsub.py   -> FigS_Ecoli_vs_Bsub_draft.*
Variants (panels A, B and the Gini only):
    --strength corrected_oe    fold over local background at the loop pixel -> *_localbg.*
    --strength loop_strength   the macro's loop_strength (main-pass loops only) -> *_loopstrength.*
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
import argparse
import sys
from itertools import combinations
import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.path
import matplotlib.patches
from matplotlib.colors import LinearSegmentedColormap

ROOT = HNS_ROOT + "Bactoloops/"
sys.path.insert(0, ROOT + "Dugar et al data/loop_network_analysis")
import loop_network_stats as S                                  # noqa: E402
from loop_engagement import window_sum                          # noqa: E402
sys.path.insert(0, ROOT + "Supplemental_figure_Ecoli_vs_Bsub")
import localbg                                                  # noqa: E402

OUT = ROOT + "Supplemental_figure_Ecoli_vs_Bsub/"
S.ANCHOR_GAP = 1                     # same anchor definition as the macro's anchorMergeGap
S.DATASETS = {
    "E. coli": (ROOT + "FIJI_output/loop_calls_wt_96h_repavg.csv",
                ROOT + "iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif",
                367, 368, ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv"),
    "B. subtilis": (ROOT + "FIJI_output/loop_calls_bsub_WT_st_merged.csv",
                    ROOT + "iced_matrix_conversion/Dugar_HiC/Dugar_WT_st_merged_iced.tif",
                    0, 0, ROOT + "Dugar et al data/rok_zscore_10kb.csv"),
}
# Integrated excess is baseline-corrected throughout (excess_baseline.py): the macro's value minus the map's mean
# integrated excess at the same distance, so the two organisms' maps (different noise, different rolling-ball offset)
# are compared on the same footing; unlooped pairs average ~0 in both.
import excess_baseline as EB                                    # noqa: E402
_CACHE = {"E. coli": "ecoli_wt96_avg", "B. subtilis": "bsub_wt_st_merged"}
BASE = {k: EB.baseline(v[1], _CACHE[k]) for k, v in S.DATASETS.items()}


def read_local(k):
    """the macro's loop table for dataset k, with integrated_excess replaced by the baseline-corrected value"""
    T_ = pd.read_csv(S.DATASETS[k][0])
    T_["integrated_excess_raw"] = T_.integrated_excess
    T_["integrated_excess"] = EB.corrected_loops(T_, BASE[k])
    return T_


ap = argparse.ArgumentParser()
ap.add_argument("--strength", choices=["corrected_oe", "loop_strength", "integrated_excess"],
                default="integrated_excess")
STRENGTH = ap.parse_args().strength
LOOPSTRENGTH = STRENGTH == "loop_strength"
INTEGRATED = STRENGTH == "integrated_excess"
SUFFIX = {"corrected_oe": "localbg", "loop_strength": "loopstrength", "integrated_excess": "integratedexcess"}[STRENGTH]
OUT_NAME = "FigS_Ecoli_vs_Bsub_draft" + ("" if INTEGRATED else "_" + SUFFIX)
CIRCLE_PNG = ROOT + f"CircularLoopArcs/circular_ecoli_vs_bsub_repavg_bsubmerged_relative_mediumconf_{SUFFIX}_panel.png"
# the macro's loop tables: log2_obs_exp, oe_percentile, confidence and integrated_excess are all
# measured on the local-background-corrected map (macro step 7c)
LOCAL_TABLES = {k: v[0] for k, v in S.DATASETS.items()}
CACHE_NAME = {"E. coli": "ecoli_wt96_avg", "B. subtilis": "bsub_wt_st_merged"}
ROK_DEP = ROOT + "Dugar et al data/loop_network_analysis/rok_dependence_loops.csv"
RC = {"RC1": 54, "RC2": 206, "RC3": 226, "RC4": 325, "RC5": 346, "RC6": 366, "RC7": 412, "RC8": 417}

COL = {"E. coli": "#2a78d6", "B. subtilis": "#eb6834"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
DIV = LinearSegmentedColormap.from_list("div", ["#1f5fa8", "#f0efec", "#c0392b"])
_hs_sys = __import__("sys"); _hs_sys.path.insert(0, ROOT + "updated figures/")
import house_style as HS    # noqa: E402
HS.apply()
from matplotlib.backends.backend_pdf import PdfPages
rng = np.random.default_rng(0)

# ---------------- analysis on the merged data ----------------
# S.analyze supplies the network structure (anchors, called pairs); every strength
# and percentile is then re-read from the local-background-corrected map.
full = {k: S.analyze(k, *v) for k, v in S.DATASETS.items()}
LB = localbg.corrected_maps({CACHE_NAME[k]: v[1] for k, v in S.DATASETS.items()})
corr = {k: LB[CACHE_NAME[k]]["bs"] for k in S.DATASETS}
for k in full:
    res, L, P, groups, rep, oe = full[k]
    A = corr[k]
    stats_ = [localbg.distance_percentile(A, rep[i], rep[j]) for i, j in zip(P.i, P.j)]
    P["oe"] = [A[rep[i], rep[j]] for i, j in zip(P.i, P.j)]
    P["pct"] = [q for q, _ in stats_]
    P["null_med"] = [med for _, med in stats_]
    P["p95"] = [np.percentile(A[np.arange(len(A)), (np.arange(len(A)) + S.cdist(rep[i], rep[j], len(A))) % len(A)], 95)
                for i, j in zip(P.i, P.j)]            # distance-matched 95th percentile
    L["oe"] = [A[a, b] for a, b in zip(L.bin1.round().astype(int), L.bin2.round().astype(int))]
    full[k] = (res, L, P, groups, rep, A)

eng = {}
for k, (loops_csv, tif, *_r) in S.DATASETS.items():
    d = LB[CACHE_NAME[k]]
    m, E, bg = d["m"], d["E"], d["bg"]
    n = len(m)
    i_, j_ = np.indices(m.shape)
    dist = np.minimum(abs(i_ - j_), n - abs(i_ - j_))
    nonlocal_sum = np.median(np.where(dist >= 2, m, 0).sum(1)[m.sum(1) > 0])
    # excess over the LOCAL background: expected at that distance x the large-scale factor the rolling ball removed
    Xf = np.where(m > 0, m - E * 2.0 ** bg, 0.0) / nonlocal_sum
    res, L, P, groups, rep, oe = full[k]
    b1 = L.bin1.round().astype(int).values
    b2 = L.bin2.round().astype(int).values
    B_ = BASE[k]                                                 # baseline-corrected integrated excess (%)
    far = np.array([S.cdist(a, b, n) >= 10 for a, b in zip(b1, b2)])
    per_loop = B_.corrected(b1[far], b2[far])
    per_anchor = []
    for i in range(len(rep)):
        js = [j for j in range(len(rep)) if j != i and S.cdist(rep[i], rep[j], n) >= 10]
        per_anchor.append(B_.corrected(np.full(len(js), rep[i]), np.array([rep[j] for j in js])).sum())
    eng[k] = (np.asarray(per_loop), np.array(per_anchor))

ec, bs = full["E. coli"], full["B. subtilis"]
k_anchor, k_loops = len(bs[3]), len(bs[1])
sub = np.array([S.anchor_metrics(ec[2], set(rng.choice(len(ec[3]), k_anchor, replace=False)))
                for _ in range(1000)])
bs_adj, bs_frac = S.anchor_metrics(bs[2], set(range(k_anchor)))
exf = lambda L: np.clip(2 ** np.asarray(L.oe.values, float) - 1, 0, None)
if LOOPSTRENGTH:   # Gini of the macro's loop_strength, main-pass loops only
    exf = lambda L: L.loop_strength.dropna().values
    k_loops = len(exf(bs[1]))
if INTEGRATED:     # Gini of integrated excess over the local background, all loops
    exf = lambda L: np.clip(L.integrated_excess.values, 0, None)
    for kk in full:   # the corrected tables are row-for-row copies of the macro tables
        full[kk][1]["integrated_excess"] = read_local(kk).integrated_excess.values
    k_loops = len(exf(bs[1]))
g_sub = np.array([S.gini(rng.choice(exf(ec[1]), k_loops, replace=False)) for _ in range(1000)])
g_bs = S.gini(exf(bs[1]))

# ChIP-seq peak pairs, independent of the loop calls. Peaks: local maxima of the
# 10 kb z-score track (within +/-3 bins) at or above a threshold; every pair of
# peaks >= 100 kb apart is scored as corrected log2(obs/exp) minus the median of
# all pixels at the same genomic distance, and called enriched above that
# distribution's 95th percentile. A peak pair "matches a called loop" when both
# ends are within PEAK_MATCH bins of the two ends of a called loop.
PEAK_THRESHOLDS = [1.5, 2, 2.5, 3]
PEAK_MAIN = 2
PEAK_MATCH = 2                       # bins (20 kb), the half-width of the integrated-excess window


def chip_peaks(z, thr, sep=3):
    n = len(z)
    return [i for i in range(n) if np.isfinite(z[i]) and z[i] >= thr
            and z[i] == np.nanmax(z[[(i + k) % n for k in range(-sep, sep + 1)]])]


peak_pairs = {}
for k, (loops_csv, tif, ori, recenter, chip_csv) in S.DATASETS.items():
    oe = full[k][5]
    n = len(oe)
    ch = pd.read_csv(chip_csv)
    zg = np.full(n, np.nan)
    zg[(ch.start_bp // 10_000).astype(int).values % n] = ch.zscore.values
    zt = np.array([zg[(recenter - 1 + t) % n] if recenter else zg[t] for t in range(n)])   # TIFF order
    ends = list(zip(full[k][1].bin1.round().astype(int), full[k][1].bin2.round().astype(int)))
    near = lambda a, b: any((S.cdist(a, x, n) <= PEAK_MATCH and S.cdist(b, y, n) <= PEAK_MATCH) or
                            (S.cdist(a, y, n) <= PEAK_MATCH and S.cdist(b, x, n) <= PEAK_MATCH) for x, y in ends)
    for thr in PEAK_THRESHOLDS:
        pk = chip_peaks(zt, thr)
        rows = []
        for a, b in combinations(pk, 2):
            if S.cdist(a, b, n) < 10:
                continue
            q, med = localbg.distance_percentile(oe, a, b)       # oe is the corrected map here
            d_ = min((a - b) % n, (b - a) % n)
            p95 = np.percentile(oe[np.arange(n), (np.arange(n) + d_) % n], 95)
            rows.append((oe[a, b] - med, q, near(a, b), p95 - med))
        peak_pairs[(k, thr)] = (len(pk), np.array(rows))

# ---------------- layout ----------------
def trimmed(path, pad=10):
    """image with its white border cropped off"""
    img = plt.imread(path)
    ink = np.where((img[..., :3] < 0.97).any(axis=2))
    r0, r1 = max(ink[0].min() - pad, 0), min(ink[0].max() + pad, img.shape[0])
    c0, c1 = max(ink[1].min() - pad, 0), min(ink[1].max() + pad, img.shape[1])
    return img[r0:r1, c0:c1]


circle_img = trimmed(CIRCLE_PNG)
FIG_W, LEFT, RIGHT = 500 / 72, 0.1, 0.98                  # house style: 500 pt wide, two pages (a-b, c-h)
row_a_in = FIG_W * (RIGHT - LEFT) * circle_img.shape[0] / circle_img.shape[1]   # height that fits the image width
ROW_MAPS_IN = 3.6                                   # two square loop maps side by side, plus titles
fig = fig1 = plt.figure(figsize=(FIG_W, ROW_MAPS_IN + row_a_in + 0.6))
outer = fig.add_gridspec(2, 1, height_ratios=[ROW_MAPS_IN, row_a_in], hspace=0.12,
                         left=LEFT, right=RIGHT, top=0.96, bottom=0.02)
fig2 = plt.figure(figsize=HS.size(500, 600))
lower = fig2.add_gridspec(3, 3, hspace=0.75, wspace=0.75, left=0.11, right=0.98, top=0.95, bottom=0.08)


def letter(ax, s, dx=-0.02, dy=1.0):
    HS.letter(ax, s.lower(), dx, dy)


# A. the macro's loop maps for B. subtilis WT and delta-rok: Loop_density_smoothed with every loop
#    circled (lighter, dashed-style colour = low confidence), ridge-search loops as green crosses,
#    and the ridge-mask region outlined in cyan
gMaps = outer[0].subgridspec(1, 2, wspace=0.12)
for c, (fname, title) in enumerate([("bsub_WT_st_merged", "$\\it{B. subtilis}$ WT"),
                                    ("bsub_rok_st_merged", "$\\it{B. subtilis}$ Δ$\\it{rok}$")]):
    ax = fig.add_subplot(gMaps[0, c])
    img = tifffile.imread(ROOT + f"FIJI_output/Loop_density_smoothed_loops_{fname}.tif")
    nb = len(pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{fname}.csv"))
    L_ = pd.read_csv(ROOT + f"FIJI_output/loop_calls_{fname}.csv")
    mb = nb * 10 / 1000
    ax.imshow(img, extent=[0, mb, mb, 0], interpolation="lanczos")
    ax.set_xticks([0, 1, 2, 3, 4]); ax.set_yticks([0, 1, 2, 3, 4])
    ax.tick_params(labelsize=HS.BASE, length=2)
    ax.set_xlabel("genome position (Mb)", fontsize=HS.BASE)
    if c == 0:
        ax.set_ylabel("genome position (Mb)", fontsize=HS.BASE)
    nr = int(L_.note.astype(str).str.startswith("ridge").sum())
    ax.set_title(f"{title}\n{len(L_)} loops ({nr} in the ridge region)", fontsize=HS.TITLE, pad=3)
    for sp_ in ax.spines.values():
        sp_.set_visible(True); sp_.set_color("#9a9a9a"); sp_.set_linewidth(0.5)
    if c == 0:
        letter(ax, "A", dx=-0.18, dy=1.0)

# B. circular networks
ax = fig.add_subplot(outer[1])
ax.imshow(circle_img, interpolation="lanczos")
ax.axis("off")
letter(ax, "B", dx=0.0, dy=0.97)

# B. distribution of loop strengths: kernel density estimate of log10(strength / dataset mean)
from scipy.stats import gaussian_kde                             # noqa: E402
fig = fig2
ax = fig.add_subplot(lower[0, 0])
grid = np.linspace(np.log10(0.08), np.log10(12), 400)
for k, (res, L, *_r) in full.items():
    Lc = read_local(k)                      # corrected strength and confidence
    Lm = Lc[Lc.confidence.isin(["high", "medium"])]
    v = (Lm.loop_strength.dropna().values if LOOPSTRENGTH else
         Lm.integrated_excess.values if INTEGRATED else 2 ** Lm.log2_obs_exp.values)
    v = v[v > 0]
    lr = np.log10(v / v.mean())
    dens = gaussian_kde(lr)(grid)                          # Scott's rule bandwidth
    ax.fill_between(10 ** grid, dens, color=COL[k], alpha=0.2, lw=0)
    ax.plot(10 ** grid, dens, color=COL[k], lw=1.5)
    ax.plot(10 ** lr, np.full(len(lr), -0.12 if k == "E. coli" else -0.3), "|", color=COL[k], ms=4, mew=0.6)
    peak = grid[np.argmax(dens)]
    ax.text(10 ** peak * (1.7 if k == "E. coli" else 4.5), dens.max() * (0.9 if k == "E. coli" else 0.6), f"$\\it{{{k}}}$\n(n = {len(lr)})",
            color=COL[k], fontsize=HS.BASE, va="center", ha="left")
ax.axvline(1, color=MUTED, lw=0.6, ls=":")
ax.axhline(0, color=MUTED, lw=0.6)
ax.set_xscale("log")
ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5, 10])
ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5", "10"])
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlim(0.08, 12)
ax.set_ylim(bottom=-0.45)
ax.set_ylabel("density (per log10 unit)")
ax.set_xlabel("loop_strength (area x intensity), relative to mean loop" if LOOPSTRENGTH else
              "integrated excess over local background,\nrelative to mean loop" if INTEGRATED else
              "fold over local background,\nrelative to mean loop")
ax.grid(axis="x", color=GRID, lw=0.4)
letter(ax, "C", dx=-0.22)

# C. every anchor pair, split into called loops and uncalled pairs: corrected log2(obs/exp) above the
#    distance-matched median, scaled to the distance-matched 95th percentile (0 = median, 1 = 95th pct).
#    C and E share this scale and y range.
ax = fig.add_subplot(lower[0, 1:3])
ce_top = max(max(np.max(r[:, 0] / r[:, 3]) for (kk, t), (_n, r) in peak_pairs.items() if t == PEAK_MAIN),
             max(np.max((P.oe - P.null_med) / (P.p95 - P.null_med)) for P in [full[k][2] for k in full]))
ce_bot = min(min(np.min(r[:, 0] / r[:, 3]) for (kk, t), (_n, r) in peak_pairs.items() if t == PEAK_MAIN),
             min(np.min((P.oe - P.null_med) / (P.p95 - P.null_med)) for P in [full[k][2] for k in full])) - 0.7
C_POS = {("E. coli", True): 0, ("E. coli", False): 1, ("B. subtilis", True): 2.4, ("B. subtilis", False): 3.4}
# Rotation test for each group: the median of the same anchor pairs after shifting the whole anchor set
# around the circular chromosome (keeps pair distances and shared anchors; distance-matched median and
# 95th percentile are unchanged by the shift). One-sided: is the observed median above the rotations?
N_ROT = None                 # every rotation offset is used once (rotation_test.shifts)
_rr = np.random.default_rng(11)
rot_test = {}
for (k, is_called) in C_POS:
    res_, L_, P_, groups_, rep_, A_ = full[k]
    n_ = len(A_)
    sel_ = P_[P_.called == is_called]
    a_ = np.array([rep_[i] for i in sel_.i]); b_ = np.array([rep_[j] for j in sel_.j])
    m_, q_ = sel_.null_med.values, sel_.p95.values
    obs_ = np.median((sel_.oe.values - m_) / (q_ - m_))
    null_ = np.array([np.median((A_[(a_ + s) % n_, (b_ + s) % n_] - m_) / (q_ - m_)) for s in RT.shifts(n_)])
    rot_test[(k, is_called)] = (obs_, RT.p_value(obs_, null_, two_sided=False), np.percentile(null_, 95), null_.max())
stars = RT.stars
for (k, is_called), x0 in C_POS.items():
    P = full[k][2]
    sel = P[P.called == is_called]
    y = ((sel.oe - sel.null_med) / (sel.p95 - sel.null_med)).values
    ax.text(x0, ce_top + 0.2, stars(rot_test[(k, is_called)][1]), ha="center", va="bottom", fontsize=HS.BASE, color=INK)
    vp = ax.violinplot(y, [x0], widths=0.75, showextrema=False)
    for body in vp["bodies"]:
        body.set_facecolor(COL[k])
        body.set_alpha(0.3 if is_called else 0.15)
        body.set_edgecolor("none")
    jit = rng.uniform(-0.18, 0.18, len(y))
    if is_called:
        ax.scatter(x0 + jit, y, s=5, facecolors=COL[k], edgecolors="none", alpha=0.9)
    else:
        ax.scatter(x0 + jit, y, s=5, facecolors="none", edgecolors=COL[k], lw=0.4, alpha=0.8)
    ax.hlines(np.median(y), x0 - 0.3, x0 + 0.3, color=INK, lw=1.2)
    ax.text(x0, ce_top + 1.3, f"{100 * np.mean(sel.pct > 0.5):.0f}%\n{100 * np.mean(sel.pct > 0.95):.0f}%",
            ha="center", va="bottom", fontsize=HS.BASE, color=COL[k], fontweight="bold")
    ax.text(x0, ce_bot + 0.1, f"n={len(y)}", ha="center", va="bottom", fontsize=HS.BASE, color=MUTED)
ax.axhline(0, color=MUTED, lw=0.6)
ax.axhline(1, color=MUTED, lw=0.6, ls=":")
ax.set_xticks(list(C_POS.values()))
ax.set_xticklabels(["called\nloops", "uncalled\npairs"] * 2, fontsize=HS.BASE)
for k, xc in [("E. coli", 0.5), ("B. subtilis", 2.9)]:
    ax.text(xc, -0.22, f"$\\it{{{k}}}$ ({'H-NS' if k == 'E. coli' else 'Rok'} loop anchors)", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=HS.BASE, color=COL[k])
ax.set_xlim(-0.6, 4.0)
ax.set_ylim(ce_bot, ce_top + 4.2)
ax.set_ylabel("enrichment over background\n(1 = 95th percentile)", fontsize=HS.BASE)
ax.tick_params(labelsize=HS.BASE)
letter(ax, "D", dx=-0.09)

# D. mean excess over local background (integrated excess), relative to E. coli, with bootstrap 95% CIs
ax = fig.add_subplot(lower[1, 0])
rel_d = {}
for g, (lab, idx) in enumerate([("per loop", 0), ("per anchor\n(all partners)", 1)]):
    ref = eng["E. coli"][idx]
    boot_ref = np.array([rng.choice(ref, len(ref)).mean() for _ in range(2000)])
    for j, k in enumerate(["E. coli", "B. subtilis"]):
        vals = eng[k][idx]
        boot = np.array([rng.choice(vals, len(vals)).mean() for _ in range(2000)])
        ratio = boot / (boot_ref if k == "B. subtilis" else ref.mean())
        m_ = vals.mean() / ref.mean()
        lo, hi = np.percentile(ratio, [2.5, 97.5])
        rel_d[(k, lab.split("\n")[0])] = (m_, lo, hi)
        x0 = g * 2.6 + j
        ax.bar(x0, m_, width=0.8, color=COL[k])
        ax.errorbar(x0, m_, yerr=[[m_ - lo], [hi - m_]], color=INK, lw=0.8, capsize=2)
        ax.text(x0, hi + 0.06, f"{m_:.1f}×", ha="center", va="bottom", fontsize=HS.BASE, color=INK)
ax.axhline(1, color=MUTED, lw=0.6, ls=":")
ax.set_xticks([0.5, 3.1])
ax.set_xticklabels(["per loop", "per anchor\n(all partners)"], fontsize=HS.BASE)
ax.set_xlim(-0.7, 4.3)
ax.set_ylim(0, max(v[2] for v in rel_d.values()) + 2.2)
ax.set_ylabel("mean integrated excess,\nrelative to $\\it{E. coli}$", fontsize=HS.BASE)
ax.tick_params(labelsize=HS.BASE)
ax.grid(axis="y", color=GRID, lw=0.4)
ax.set_axisbelow(True)
ax.legend(handles=[matplotlib.patches.Patch(color=COL[k], label=k) for k in ["E. coli", "B. subtilis"]],
          frameon=False, loc="upper right", handlelength=1, prop={"style": "italic", "size": HS.BASE})
letter(ax, "E", dx=-0.22)

# F. integrated loop density against the anchor protein's ChIP-seq z-score, bin by bin
from scipy.stats import pearsonr, spearmanr                     # noqa: E402
gF = lower[1, 1:3].subgridspec(1, 2, wspace=0.35)
ild_corr = {}
for c, (k, fname, prot) in enumerate([("E. coli", "wt_96h_repavg", "H-NS"), ("B. subtilis", "bsub_WT_st_merged", "Rok")]):
    ax = fig.add_subplot(gF[0, c])
    ild = pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{fname}.csv")
    nb = len(ild)
    ch = pd.read_csv(S.DATASETS[k][4])
    zg = np.full(nb, np.nan)
    zg[(ch.start_bp // 10_000).astype(int).values % nb] = ch.zscore.values
    gbin = (ild.genome_position_bp.values // 10_000).astype(int) % nb   # genome bins (positions may be signed)
    zz, yy = zg[gbin], ild.ILD_au.values
    ok = np.isfinite(zz) & np.isfinite(yy)
    r, pr = pearsonr(zz[ok], yy[ok]); rho, _ = spearmanr(zz[ok], yy[ok])
    ild_corr[k] = (r, rho, int(ok.sum()))
    ax.scatter(zz[ok], yy[ok], s=3, color=COL[k], alpha=0.5, lw=0)
    fit = np.polyfit(zz[ok], yy[ok], 1)
    xs = np.linspace(np.nanmin(zz[ok]), np.nanmax(zz[ok]), 50)
    ax.plot(xs, np.polyval(fit, xs), color=INK, lw=0.8)
    ax.text(0.04, 0.96, f"r = {r:.2f}\nρ = {rho:.2f}\nn = {ok.sum()} bins", transform=ax.transAxes,
            ha="left", va="top", fontsize=HS.BASE, color=INK)
    ax.set_xlabel(f"{prot} ChIP-seq z-score", fontsize=HS.BASE)
    if c == 0:
        ax.set_ylabel("integrated loop density (a.u.)", fontsize=HS.BASE)
    ax.set_title(f"$\\it{{{k}}}$", fontsize=HS.BASE, pad=3, color=COL[k])
    ax.tick_params(labelsize=HS.BASE)
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    if c == 0:
        letter(ax, "F", dx=-0.3, dy=1.0)

# G. silencing at silencer-bound regions in the two growth phases (Dugar et al. RNA-seq vs this study; same metric,
# 10-kb bins): mean log2(mutant/WT) at bins with ChIP z >= 2 minus the genome-wide mean. From
# Dugar et al data/rok_vs_hns_silencing.py (run that first; it writes rok_vs_hns_silencing.csv).
SIL = pd.read_csv(ROOT + "Dugar et al data/rok_vs_hns_silencing.csv")
ax = fig.add_subplot(lower[2, 0:2])
groups = [("E. coli", "E. coli (H-NS)", "H-NS ChIP z>=2", "H-NS", ["exponential", "deep stationary 96 h"], ["exp", "96 h"]),
          ("B. subtilis", "B. subtilis (Rok)", "Rok ChIP z>=2", "Rok", ["exponential", "stationary 22 h"], ["exp", "22 h"])]
for g, (k, org, setname, prot, phases, labs) in enumerate(groups):
    silg = SIL[(SIL.organism == org) & (SIL.set == setname)].set_index("phase")
    xs = [g * 2.6, g * 2.6 + 1]
    ys = [2 ** silg.loc[ph, "offset"] for ph in phases]
    ax.plot(xs, ys, color=COL[k], lw=0.9, zorder=2)
    for x, ph, lab in zip(xs, phases, labs):
        r = silg.loc[ph]
        ax.errorbar(x, 2 ** r.offset, yerr=[[2 ** r.offset - 2 ** r.lo], [2 ** r.hi - 2 ** r.offset]], fmt="o", ms=5,
                    color=COL[k], mec=INK, mew=0.4, lw=0.8, capsize=2, zorder=3)
        ax.text(x, 2 ** r.hi * 1.06, stars(r.p), ha="center", va="bottom", fontsize=HS.BASE)
        ax.text(x, 0.61, lab, ha="center", va="top", fontsize=HS.BASE)
    c = silg.loc["change st - ex"]
    top = max(2 ** silg.loc[phases[1], "hi"], 2 ** silg.loc[phases[0], "hi"])
    ax.text(np.mean(xs), top * 1.5, f"change {2 ** c.offset:.1f}× ({2 ** c.lo:.1f}–{2 ** c.hi:.1f}) {stars(c.p)}",
            ha="center", va="bottom", fontsize=HS.BASE, color=INK)
    ax.text(np.mean(xs), 0.47, f"$\\it{{{k}}}$ ({prot})", ha="center", va="top", fontsize=HS.BASE,
            color=COL[k], transform=ax.transData)
ax.axhline(1, color=MUTED, lw=0.6, ls=":")
ax.set_yscale("log", base=2)
ax.set_yticks([0.5, 1, 2, 4, 8]); ax.set_yticklabels(["0.5", "1", "2", "4", "8"])
ax.set_ylim(0.64, 16); ax.set_xlim(-0.7, 4.3); ax.set_xticks([])
ax.set_ylabel("H-NS or Rok silencing at bound\nregions, relative to genome (fold)", fontsize=HS.BASE)
ax.tick_params(labelsize=HS.BASE)
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "G", dx=-0.13)

# H. the local effect Dugar et al. reported: mean log2(delta-rok/WT) within +/-50 kb of each Rok cluster, minus the
# genome mean, exponential vs stationary phase (same bigWigs as G)
_LOC = pd.read_csv(ROOT + "Dugar et al data/rok_cluster_local.csv")   # from rok_vs_hns_silencing.py
ax = fig.add_subplot(lower[2, 2])
_RC = {"1": 54, "2": 206, "3": 226, "4": 325, "5": 346, "6": 366, "7": 412, "8": 417}
for i_, (nm, c_) in enumerate(_RC.items()):
    for j_, (ph, colr) in enumerate([("ex", "#c9c6bf"), ("st", COL["B. subtilis"])]):
        v = 2 ** float(_LOC.set_index("cluster").loc[f"RC{nm}", ph])
        ax.bar(i_ + (j_ - 0.5) * 0.38, v - 1, bottom=1, width=0.36, color=colr, edgecolor=INK, lw=0.3)
ax.axhline(1, color=MUTED, lw=0.6)
ax.set_xticks(range(8)); ax.set_xticklabels([f"RC{k}" for k in _RC], fontsize=HS.BASE, rotation=45, ha="right",
                                            rotation_mode="anchor")
ax.set_ylabel("Rok silencing within ±50 kb of the\ncluster, relative to genome (fold)", fontsize=HS.BASE)
ax.set_ylim(0.7, 1.6)
ax.legend(handles=[matplotlib.patches.Patch(color="#c9c6bf", label="exp"),
                   matplotlib.patches.Patch(color=COL["B. subtilis"], label="22 h")], fontsize=HS.BASE, frameon=False,
          loc="upper right", handlelength=1)
ax.tick_params(labelsize=HS.BASE)
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "H", dx=-0.3)

for f_ in (fig1, fig2):
    HS.check(f_)
with PdfPages(OUT + OUT_NAME + ".pdf") as pp:              # two pages: a-b, c-h
    pp.savefig(fig1, dpi=300); pp.savefig(fig2, dpi=300)
fig1.savefig(OUT + OUT_NAME + ".png", dpi=300)
fig2.savefig(OUT + OUT_NAME + "_p2.png", dpi=300)

summary = {"E. coli": {}, "B. subtilis": {}}
for k in summary:
    res, L, P, groups, rep, _A = full[k]
    unc = P[~P.called]
    Lc = read_local(k)
    summary[k]["map"] = "local-background corrected (rolling ball 50)"
    summary[k]["loops called"] = len(L)
    summary[k]["loops medium+ confidence (corrected)"] = int(Lc.confidence.isin(["high", "medium"]).sum())
    summary[k]["anchors"] = len(rep)
    summary[k]["Gini of called-loop excess"] = round(S.gini(exf(L)), 3)
    summary[k]["uncalled anchor pairs > median (%)"] = round(100 * np.mean(unc.pct > 0.5))
    summary[k]["uncalled anchor pairs > 95th pct (%)"] = round(100 * np.mean(unc.pct > 0.95))
    summary[k]["anchor stickiness adjusted R^2"] = round(S.anchor_metrics(P, set(range(len(rep))))[0], 3)
    summary[k]["median excess per loop (%)"] = np.median(eng[k][0])
    summary[k]["median excess per anchor (%)"] = np.median(eng[k][1])
summary["E. coli"]["subsampled adj R^2 2.5-97.5%"] = np.percentile(sub[:, 0], [2.5, 97.5]).round(2).tolist()
summary["E. coli"]["subsampled uncalled-above-null 2.5-97.5%"] = np.percentile(sub[:, 1], [2.5, 97.5]).round(2).tolist()
summary["E. coli"]["subsampled Gini 2.5-97.5%"] = np.percentile(g_sub, [2.5, 97.5]).round(2).tolist()
for (k, lab), (m_, lo, hi) in rel_d.items():
    summary[k][f"mean excess {lab}, relative to E. coli (95% bootstrap CI)"] = f"{m_:.2f} ({lo:.2f}-{hi:.2f})"
    summary[k][f"mean excess {lab} (%)"] = round(float(eng[k][0 if lab == "per loop" else 1].mean()), 2)
for k in summary:
    P = full[k][2]
    for is_called, lab in [(True, "called loops"), (False, "uncalled pairs")]:
        sel = P[P.called == is_called]
        summary[k][f"{lab}: n / % > median / % > 95th pct"] = (
            f"{len(sel)} / {100 * np.mean(sel.pct > 0.5):.0f}% / {100 * np.mean(sel.pct > 0.95):.0f}%")
from scipy.stats import skew                                    # noqa: E402
_rs = np.random.default_rng(5)
for k in summary:                    # skewness of the panel C distributions (medium+ loops, relative to mean)
    Lc = read_local(k); Lc = Lc[Lc.confidence.isin(["high", "medium"])]
    v = Lc.integrated_excess.values; v = v[v > 0]; v = v / v.mean()
    summary[k]["Gini of integrated excess, medium+ loops (panel C)"] = round(S.gini(v), 3)
    if k == "E. coli":
        nb_ = int((read_local("B. subtilis").confidence.isin(["high", "medium"])).sum())
        g_ = [S.gini(_rs.choice(v, nb_, replace=False)) for _ in range(5000)]
        summary[k]["Gini, medium+ loops, subsampled to the B. subtilis count 2.5-97.5%"] = np.percentile(g_, [2.5, 97.5]).round(2).tolist()
    for lab, x in [("linear", v), ("log10", np.log10(v))]:
        bt = np.percentile([skew(_rs.choice(x, len(x))) for _ in range(5000)], [2.5, 97.5])
        summary[k][f"skewness of integrated excess ({lab}; bootstrap 95% CI)"] = f"{skew(x):.2f} ({bt[0]:.2f}-{bt[1]:.2f})"
for (k, is_called), (obs_, p_, q95_, mx_) in rot_test.items():
    lab_ = "called loops" if is_called else "uncalled pairs"
    summary[k][f"panel D rotation test, {lab_}: median / rotation 95th pct / rotation max / p (all {len(RT.shifts(461))}/{len(RT.shifts(422))} rotations, one-sided)"] = (
        f"{obs_:.2f} / {q95_:.2f} / {mx_:.2f} / {p_:.4f} {stars(p_)}")
for k, (r, rho, nbin) in ild_corr.items():
    summary[k]["ILD vs ChIP z: Pearson r / Spearman rho / bins"] = f"{r:.3f} / {rho:.3f} / {nbin}"
for k in summary:
    for t in PEAK_THRESHOLDS:
        npk, rows = peak_pairs[(k, t)]
        summary[k][f"ChIP peaks z>={t}: peaks / pairs / % enriched / median above background"] = (
            f"{npk} / {len(rows)} / {100 * np.mean(rows[:, 1] > 0.95):.0f}% / {np.median(rows[:, 0]):+.2f}")
        summary[k][f"ChIP peaks z>={t}: % pairs > median"] = round(100 * np.mean(rows[:, 1] > 0.5))
        for m_, lab in [(True, "matching a called loop"), (False, "not matching a called loop")]:
            r = rows[rows[:, 2] == m_]
            summary[k][f"ChIP peaks z>={t}, pairs {lab}: n / % > median / % > 95th pct"] = (
                f"{len(r)} / {100 * np.mean(r[:, 1] > 0.5):.0f}% / {100 * np.mean(r[:, 1] > 0.95):.0f}%")
pd.DataFrame(summary).to_csv(OUT + OUT_NAME.replace("_draft", "_numbers") + ".csv")
print(pd.DataFrame(summary).to_string())
print("saved")
