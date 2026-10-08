"""Supplementary figure: H-NS association of the loops in every 96 h map (WT, NAP mutants, long-cell mutants).
Panels d-j of the previous NAP-mutant figure (figureS3_nap_mutants.py), extended to the long-cell mutants.

  a-g  Loops of WT 96 h (rep1 + rep2 average, as Fig. 2f), delta-stpA, delta-fis, delta-ihfAB, delta-zapA, delta-ftsP and
       delta-envC on the circular chromosome (oriC at top, ter at dif): arc width and opacity scale with baseline-corrected
       integrated excess (excess_baseline.py) on one absolute scale for all panels; grey ring, anti-H-NS ChIP-seq z-score
       (WT 96 h).
  h    Pearson r between ILD and anti-H-NS ChIP-seq (WT 96 h), bin by bin; exact circular-shift P (shifts >= 200 kb).
  i    Mean ChIP z-score of the two anchors of each called loop (box: median, IQR; whiskers: full range; points: loops);
       grey dash, 95th percentile of the loop set rotated by every offset; exact rotation P. delta-hns has no loops.
  j    Mean baseline-corrected integrated excess at every pair of H-NS ChIP-seq peaks (z >= 2 local maxima within +/-3
       bins, pairs >= 100 kb apart), independent of loop calls; 95% CI, bootstrap over peaks; grey dash, 95th percentile
       of the rotated peak set; exact rotation P.
Output: FigureS_HNS_comparison.pdf / .png, s_hns_comparison_stats.csv
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Wedge, PathPatch
from matplotlib.path import Path

ROOT = HNS_ROOT + "Bactoloops/"
OUT = ROOT + "updated figures/"
IC = ROOT + "iced_matrix_conversion/"
sys.path.insert(0, ROOT)
sys.path.insert(0, ROOT + "Supplemental_figure_Ecoli_vs_Bsub")
import rotation_test as RT                      # noqa: E402
import excess_baseline as EB                    # noqa: E402

N, RC, ORI, TER = 461, 368, 367, 155            # bins; TIFF recentre bin; oriC and dif genome bins
tiff_of = lambda g: (g - (RC - 1)) % N
genome_of = lambda t: (RC - 1 + t) % N
TICKS_KB = [3670, 0, 1000, 2000, 3000]
TICK_POS = [tiff_of(k // 10) for k in TICKS_KB]
# label, macro run, ICE TIFF, colour, localbg cache name
CASES = [("Δ$\\it{hns}$", "HNSm_96h", "cWX2870_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", "#1b6f6f", "hns_96h"),
         ("WT", "wt_96h_repavg", "cWX2200_96h_rep1rep2_average_10000_iced.tif", "#000000", "ecoli_wt96_avg"),
         ("Δ$\\it{stpA}$", "stpA_96h", "cWX3338_96h_Sau3AI_rep1_10000_iced.tif", "#a31f1f", "stpA_96h"),
         ("Δ$\\it{fis}$", "fis_96h", "cWX2478_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", "#2e6b1f", "fis_96h"),
         ("Δ$\\it{ihfAB}$", "IHF_96h", "cWX2943_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", "#1f2fa3", "ihf_96h"),
         ("Δ$\\it{zapA}$", "zapA_96h", "cWX3541_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif", "#c07a2c", "zapA_96h"),
         ("Δ$\\it{ftsP}$", "ftsP_96h", "cWX3543_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif", "#6a51a3", "ftsP_96h"),
         ("Δ$\\it{envC}$", "envC_96h", "cWX3550_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif", "#b0306a", "envC_96h")]
MUTANTS = CASES[2:]
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
GREY = LinearSegmentedColormap.from_list("grey", ["#f4f4f4", "#101010"])
EX_MIN, EX_MAX = 0.0, 3.0                       # arc scale (baseline-corrected excess, %); < 1% of loops above 3
import sys as _hs_sys
_hs_sys.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()
rng = np.random.default_rng(0)


def letter(ax, s, dx=-0.1, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=HS.LETTER, fontweight="bold", va="bottom", ha="right")


def loops_of(run):
    try:
        return pd.read_csv(ROOT + f"FIJI_output/loop_calls_{run}.csv")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["bin1", "bin2", "genome_bin1", "genome_bin2", "integrated_excess"])


c = pd.read_csv(ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv")
z_gen = np.full(N, np.nan); z_gen[(c.start_bp // 10_000).astype(int) % N] = c.zscore.values
z_tif = z_gen[[genome_of(t) for t in range(N)]]
ild = {run: pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{run}.csv").ILD_au.values for _, run, *_ in CASES}
L = {run: loops_of(run) for _, run, *_ in CASES}
cdist = lambda a, b: min((a - b) % N, (b - a) % N)

# ---------------- statistics ----------------
rows = []
# i: baseline-corrected integrated excess at H-NS peak pairs
BASE = {run: EB.baseline(IC + tif, cache) for _, run, tif, _, cache in CASES}
for _, run, *_ in CASES:
    if len(L[run]):
        L[run]["excess_corr"] = EB.corrected_loops(L[run], BASE[run])
peaks = [i for i in range(N) if np.isfinite(z_tif[i]) and z_tif[i] >= 2
         and z_tif[i] == np.nanmax(z_tif[[(i + k) % N for k in range(-3, 4)]])]


def pair_excess(B, pk):
    pr = np.array([(a, b) for a, b in combinations(pk, 2) if cdist(a, b) >= 10])
    return B.corrected(pr[:, 0], pr[:, 1])


for lab, run, tif, col, cache in CASES:
    B = BASE[run]
    obs = pair_excess(B, peaks).mean()
    boot = [pair_excess(B, sorted(set(rng.choice(peaks, len(peaks))))).mean() for _ in range(2000)]
    null = np.array([pair_excess(B, [(p + s_) % N for p in peaks]).mean() for s_ in RT.shifts(N)])
    # g: ILD vs ChIP
    iz = (ild[run] - ild[run].mean()) / ild[run].std(); ok = np.isfinite(z_tif)
    r0 = np.corrcoef(iz[ok], z_tif[ok])[0, 1]
    rnull = [np.corrcoef(np.roll(iz, k)[ok], z_tif[ok])[0, 1] for k in RT.shifts(N, exclude=20)]
    # h: mean anchor z per loop
    Lr = L[run]
    row = dict(case=lab, run=run, loops=len(Lr), ild_chip_r=r0,
               ild_chip_p=(1 + np.sum(np.abs(rnull) >= abs(r0))) / (len(rnull) + 1),
               peak_pair_excess=obs, excess_lo=np.percentile(boot, 2.5), excess_hi=np.percentile(boot, 97.5),
               excess_rot95=np.percentile(null, 95), excess_p=RT.p_value(obs, null))
    if len(Lr):
        P = np.c_[Lr.bin1.round().astype(int), Lr.bin2.round().astype(int)] % N
        f = lambda Q: np.nanmean(np.nanmean(z_tif[Q], axis=1))
        rot = [f((P + s) % N) for s in RT.shifts(N)]
        row.update(anchor_z=f(P), anchor_rot95=np.percentile(rot, 95), anchor_p=RT.p_value(f(P), rot))
        row["_per_loop"] = np.nanmean(z_tif[P], axis=1)
    rows.append(row)
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k != "_per_loop"}, flush=True)
S = pd.DataFrame(rows)
S.drop(columns="_per_loop").to_csv(OUT + "s_hns_comparison_stats.csv", index=False)
print(f"{len(peaks)} H-NS peaks (z >= 2), {len(pair_excess(BASE['wt_96h_repavg'], peaks))} pairs >= 100 kb apart")

# ---------------- figure ----------------
fig = plt.figure(figsize=HS.size(500, 600))
outer = fig.add_gridspec(2, 1, height_ratios=[1.45, 1], hspace=0.22, left=0.1, right=0.985, top=0.94, bottom=0.1)
mid = outer[0].subgridspec(2, 4, wspace=0.04, hspace=0.3)
R_IN, R_OUT = 1.0, 1.1
th = lambda g: np.pi / 2 - 2 * np.pi * ((g - ORI) % N) / N
circ = []
for k, (lab, run, tif, col, cache) in enumerate([CASES[1]] + MUTANTS):
    ax = fig.add_subplot(mid[k // 4, k % 4]); circ.append(ax)
    for g in range(N):
        v = z_gen[g]
        ax.add_patch(Wedge((0, 0), R_OUT, np.degrees(th(g + 1)), np.degrees(th(g)), width=R_OUT - R_IN,
                           facecolor="white" if not np.isfinite(v) else GREY(np.clip((v + 2) / 6, 0, 1)),
                           edgecolor="none", zorder=3))
    for _, r in L[run].sort_values("excess_corr").iterrows():
        a, b = int(round(r.genome_bin1)), int(round(r.genome_bin2))
        f = np.clip((r.excess_corr - EX_MIN) / (EX_MAX - EX_MIN), 0, 1)
        p1 = np.array([np.cos(th(a + 0.5)), np.sin(th(a + 0.5))]) * R_IN
        p2 = np.array([np.cos(th(b + 0.5)), np.sin(th(b + 0.5))]) * R_IN
        sep = min(abs(b - a), N - abs(b - a)) / N
        mpt = (p1 + p2) / 2
        kk = 0.9 * (1 - 2 * sep) / max(np.hypot(*mpt), 1e-6)
        ax.add_patch(PathPatch(Path([p1, mpt * kk, p2], [Path.MOVETO, Path.CURVE3, Path.CURVE3]), fc="none",
                               ec="#c8102e", lw=0.3 + 1.9 * f, alpha=0.12 + 0.88 * f, capstyle="round", zorder=2))
    ax.plot([0], [R_OUT + 0.07], marker="v", color="#f2c200", ms=4, mec=INK, mew=0.3)
    ax.text(0, R_OUT + 0.16, "ori", ha="center", va="bottom", fontsize=HS.BASE, style="italic")
    tx, ty = np.cos(th(TER)) * (R_OUT + 0.07), np.sin(th(TER)) * (R_OUT + 0.07)
    ax.plot([tx], [ty], marker="^", color="#2a78d6", ms=4, mec=INK, mew=0.3)
    ax.text(tx, ty - 0.1, "ter", ha="center", va="top", fontsize=HS.BASE, style="italic")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.45, 1.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{lab} 96 h\n{len(L[run])} loops", fontsize=HS.TITLE, pad=0, linespacing=1.1)
    letter(ax, "abcdefg"[k], dx=0.12, dy=0.98)
# arc key in the eighth slot
kax = fig.add_subplot(mid[1, 3]); kax.axis("off")
kax.set_xlim(0, 1); kax.set_ylim(0, 1)
kax.text(0.5, 0.8, "integrated excess\nabove baseline", ha="center", va="center", fontsize=HS.BASE)
for i, v in enumerate([0.0, 0.75, 1.5, 2.25, 3.0]):
    f = (v - EX_MIN) / (EX_MAX - EX_MIN); y = 0.62 - i * 0.1
    kax.plot([0.22, 0.6], [y, y], color="#c8102e", lw=0.3 + 1.9 * f, alpha=0.12 + 0.88 * f, solid_capstyle="round")
    kax.text(0.66, y, {0.0: "≤0", 3.0: "≥3%"}.get(v, f"{v:g}"), ha="left", va="center", fontsize=HS.BASE)
kax.text(0.5, 0.02, "grey ring:\nH-NS ChIP-seq z", ha="center", va="center", fontsize=HS.BASE, color=MUTED)

bot = outer[1].subgridspec(1, 3, wspace=0.55, width_ratios=[0.95, 1.2, 0.95])
xs = np.arange(len(CASES)); xl = [c_[0] for c_ in CASES]; cols = [c_[3] for c_ in CASES]
ax = fig.add_subplot(bot[0, 0])
ax.bar(xs, S.ild_chip_r, color=cols, width=0.7)
for x, (r_, p_) in enumerate(zip(S.ild_chip_r, S.ild_chip_p)):
    ax.text(x, max(r_, 0) + 0.02, RT.stars(p_), ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color=MUTED, lw=0.5)
ax.set_ylim(min(-0.05, S.ild_chip_r.min() - 0.05), 1.0)
ax.set_ylabel("Pearson r,\nILD vs ChIP-seq")
letter(ax, "h", dx=-0.32)
ax = fig.add_subplot(bot[0, 1])
for x, row in S.iterrows():
    if "_per_loop" not in row or not isinstance(row["_per_loop"], np.ndarray):
        ax.text(x, 1.0, "no\nloops", ha="center", va="bottom", fontsize=HS.BASE, color=MUTED)
        continue
    v = row["_per_loop"]
    bp = ax.boxplot([v], positions=[x], widths=0.55, patch_artist=True, showfliers=False, whis=(0, 100),
                    medianprops=dict(color=INK, lw=0.9), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    bp["boxes"][0].set_facecolor(cols[x]); bp["boxes"][0].set_alpha(0.45); bp["boxes"][0].set_edgecolor(cols[x])
    ax.scatter(x + rng.uniform(-0.16, 0.16, len(v)), v, s=1.5, color=INK, alpha=0.35, lw=0, zorder=3)
    ax.hlines(row.anchor_rot95, x - 0.3, x + 0.3, color="#9a9a9a", lw=0.9, ls=(0, (2, 1)), zorder=4)
    ax.text(x, np.nanmax(v) + 0.15, RT.stars(row.anchor_p), ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xlim(-0.6, len(CASES) - 0.4); ax.set_ylim(-1.2, 4.9)
ax.set_ylabel("anchor H-NS\nChIP-seq z per loop")
letter(ax, "i", dx=-0.26)
ax = fig.add_subplot(bot[0, 2])
ax.bar(xs, S.peak_pair_excess, color=cols, width=0.7)
ax.errorbar(xs, S.peak_pair_excess, yerr=[S.peak_pair_excess - S.excess_lo, S.excess_hi - S.peak_pair_excess],
            fmt="none", ecolor=INK, lw=0.7, capsize=1.5)
for x, row in S.iterrows():
    ax.hlines(row.excess_rot95, x - 0.35, x + 0.35, color="#9a9a9a", lw=0.9, ls=(0, (2, 1)), zorder=4)
    ax.text(x, max(row.excess_hi, row.excess_rot95) + 0.02 * S.excess_hi.max(), RT.stars(row.excess_p),
            ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color=MUTED, lw=0.5); ax.set_ylim(0, S.excess_hi.max() * 1.15)
ax.set_ylabel("excess at H-NS\npeak pairs (%)")
letter(ax, "j", dx=-0.32)
for ax in fig.axes[-3:]:
    ax.set_xticks(xs); ax.set_xticklabels(xl, rotation=45, ha="right", rotation_mode="anchor", fontsize=HS.BASE)
    ax.axvline(1.5, color=MUTED, lw=0.4, ls=":"); ax.axvline(4.5, color=MUTED, lw=0.4, ls=":")
    ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
HS.check(fig)
fig.savefig(OUT + "FigureS_HNS_comparison.pdf", dpi=600)
fig.savefig(OUT + "FigureS_HNS_comparison.png", dpi=300)
print("saved")
