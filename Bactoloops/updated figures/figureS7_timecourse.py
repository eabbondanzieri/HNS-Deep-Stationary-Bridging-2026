"""Supplementary Figure 7 (rebuilt): H-NS-dependent looping emerges during stationary phase.

Reads hns_timecourse_analysis.py's outputs (run that first) and the macro's outputs (FIJI_output/).
  a-e  loop maps (Loop_density_smoothed_loops) for WT exp, delta-hns exp, WT 24 h, WT 96 h, delta-hns 96 h,
       with the integrated loop density (ILD, median-subtracted, one shared scale) below
  f    mean H-NS ChIP z of the called loops' anchors; grey bar, 95th percentile of rotated loop sets
  g    mean baseline-corrected integrated excess at the H-NS ChIP peak pairs (z >= 2, >= 100 kb apart; excess_baseline.py:
       the macro's integrated excess minus the map's mean at the same distance, % of a typical bin's non-local
       contacts, ~0 for unlooped pairs in any map); 95% CI, bootstrap over peaks; grey, rotation null (mean and central
       95% of the rotated peak sets, every offset); stars, exact two-sided rotation test
  h    ILD z-score at the H-NS peaks vs all other bins (mean, bootstrap 95% CI)
(Panels i-m, the CHIN/CHID analysis, moved to main Figure 3 on 2026-10-06: figure3_hairpins.py.)
The 96 h anti-H-NS ChIP-seq stands in for all time points. Output: FigureS7_timecourse.pdf / .png
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
import excess_baseline as EB  # baseline-corrected integrated excess
ICE = {"wt_Exp": ("cWX2200_exp_Azure_10percent_Sau3A_rep1_10000_iced.tif", "wt_exp"),
       "HNSm_Exp": ("cWX2870_exp_Azure_10percent_rep1_10000_iced.tif", "hns_exp"),
       "wt_24h": ("cWX2200_24H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", "wt_24h"),
       "wt_96h_repavg": ("cWX2200_96h_rep1rep2_average_10000_iced.tif", "ecoli_wt96_avg"),
       "HNSm_96h": ("cWX2870_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", "hns_96h")}
import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = HNS_ROOT + "Bactoloops/"
OUT = ROOT + "updated figures/"
N, RC = 461, 368
tiff_of = lambda g: (g - (RC - 1)) % N
TICK_POS = [tiff_of(k // 10) for k in [3670, 0, 1000, 2000, 3000]]
MAPS = [("WT exp", "wt_Exp", "#9a7b4f"), ("Δ$\\it{hns}$ exp", "HNSm_Exp", "#7fb8b8"),
        ("WT 24 h", "wt_24h", "#5a4a32"), ("WT 96 h", "wt_96h_repavg", "#000000"),
        ("Δ$\\it{hns}$ 96 h", "HNSm_96h", "#1b6f6f")]
# panels f-h only: published 3C maps of Lioy et al. 2018 (GSE107301; MG1655, minimal medium, 30 C; converted to the
# lab genome, Lioy et al./lioy_to_lab_tiff.py), drawn as open symbols to the right of a divider
LIOY = []   # Lioy et al. 2018 maps removed from this figure (their H-NS-site signal is mostly H-NS-independent; see figure_lioy_novo_temp.py)
MAPS_B = MAPS + LIOY
OPEN = {run for _, run, _ in LIOY}
KEY = {"lioy_wt_stat_repavg": "Lioy WT stat avg", "lioy_wt_exp": "Lioy WT exp", "lioy_wt_stat": "Lioy WT stat", "lioy_wt_stat_rep1": "Lioy WT stat rep1",
       "wt_Exp": "WT exp", "HNSm_Exp": "Δhns exp", "wt_24h": "WT 24 h", "wt_96h_repavg": "WT 96 h", "HNSm_96h": "Δhns 96 h"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
import sys as _hs_sys
_hs_sys.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()
T = pd.read_csv(OUT + "hns_timecourse_numbers.csv").set_index("map")
PR = np.load(OUT + "_hns_timecourse_pairs.npz")
peaks, PP, CP, z = PR["peaks"], PR["PP"], PR["CP"], PR["z_tiff"]
rng = np.random.default_rng(0)
stars = RT.stars


def letter(ax, s, dx=-0.08, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=HS.LETTER, fontweight="bold", va="bottom", ha="right")


fig = plt.figure(figsize=HS.size(500, 480))
outer = fig.add_gridspec(2, 1, height_ratios=[0.75, 1], hspace=0.38, left=0.1, right=0.99, top=0.94, bottom=0.12)

# ---------------- a-e ----------------
top = outer[0].subgridspec(2, 5, height_ratios=[1, 0.36], hspace=0.12, wspace=0.1)
ild = {run: pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{run}.csv").ILD_au.values for _, run, _ in MAPS}
ild = {k: v - np.median(v) for k, v in ild.items()}
ymin, ymax = min(v.min() for v in ild.values()) * 1.05, max(v.max() for v in ild.values()) * 1.05
pairs = []
for c, (title, run, col) in enumerate(MAPS):
    ax = fig.add_subplot(top[0, c])
    ax.imshow(tifffile.imread(ROOT + f"FIJI_output/Loop_density_smoothed_loops_{run}.tif"), extent=[0, N, 0, N],
              origin="lower", interpolation="lanczos")
    ax.set_xticks(TICK_POS); ax.set_yticks(TICK_POS); ax.set_xticklabels([])
    ax.set_yticklabels([3670, 0, 1000, 2000, 3000] if c == 0 else [], fontsize=HS.BASE)
    n_loops = int(T.loc[KEY[run], "loops"])
    ax.set_title(f"{title}\n{n_loops} loops", fontsize=HS.TITLE, pad=3)
    letter(ax, "abcde"[c], dx=-0.32 if c == 0 else -0.04)
    ax2 = fig.add_subplot(top[1, c])
    ax2.plot(np.arange(N) + 0.5, ild[run], color=col, lw=0.6)
    ax2.set_xlim(0, N); ax2.set_ylim(ymin, ymax)
    ax2.set_xticks(TICK_POS); ax2.set_xticklabels([3670, 0, 1000, 2000, 3000], fontsize=HS.BASE, rotation=45, ha="right", rotation_mode="anchor")
    if c == 0:
        ax2.set_ylabel("ILD (a.u.)", fontsize=HS.BASE)
        ax.set_ylabel("Genome position (kb)", fontsize=HS.BASE)
    else:
        ax2.set_yticklabels([])
    pairs.append((ax, ax2))

# ---------------- f, g, h ----------------
bot = outer[1].subgridspec(1, 3, wspace=0.45, width_ratios=[1, 1.45, 1])
xs = np.arange(len(MAPS_B))
labels = [t.replace("\n", " ") for t, _, _ in MAPS_B]
fc = lambda run, col: "white" if run in OPEN else col


def divider(ax, ytext):
    if not LIOY:
        return
    ax.axvline(len(MAPS) - 0.5, color=MUTED, lw=0.6, ls=(0, (3, 2)))
    ax.text(len(MAPS) + 0.5 * (len(LIOY) - 1), ytext, "Lioy et al.\n2018", ha="center", va="top", fontsize=HS.BASE, color=MUTED)

# f: called-loop anchor z
ax = fig.add_subplot(bot[0, 0])
for i, (t, run, col) in enumerate(MAPS_B):
    rw = T.loc[KEY[run]]
    if rw.loops == 0 or not np.isfinite(rw.loop_mean_anchor_z):
        ax.text(i, 0.25, "no\nloops", ha="center", va="bottom", fontsize=HS.BASE, color=MUTED)
        continue
    ax.bar(i, rw.loop_rot95, 0.6, color="#e4e2dd", lw=0, zorder=1)
    ax.plot(i, rw.loop_mean_anchor_z, "o", color=fc(run, col), ms=5, mec=col if run in OPEN else INK, mew=0.4 if run not in OPEN else 1.0, zorder=3)
    ax.text(i, rw.loop_mean_anchor_z + 0.15, stars(rw.loop_p), ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=HS.BASE)
ax.set_ylabel("anchor H-NS\nChIP z (called loops)")
ax.set_ylim(-0.4, 3.6); ax.set_xlim(-0.6, len(MAPS_B) - 0.4)
divider(ax, 3.5)
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "f", dx=-0.3)

# g: mean baseline-corrected integrated excess at the H-NS peak pairs (95% CI, bootstrap over peaks -- pairs sharing a
# peak are not independent) vs the rotation null (the whole peak-pair set shifted around the chromosome; grey, mean and
# central 95% of the rotated-set means; every offset, rotation_test.py). Stars: exact two-sided rotation test of the mean
# (floor P = 0.0044, so at most **).
N_ROT_G, N_BOOT_G = 2000, 2000
ax = fig.add_subplot(bot[0, 1])
pk_idx = {p_: k for k, p_ in enumerate(peaks)}
ia, ib = np.array([pk_idx[a] for a in PP[:, 0]]), np.array([pk_idx[b] for b in PP[:, 1]])
g_rows = []
for i, (t, run, col) in enumerate(MAPS_B):
    B = EB.baseline(ROOT + "iced_matrix_conversion/" + ICE[run][0], ICE[run][1])
    pe = B.corrected(PP[:, 0], PP[:, 1])
    rr = np.random.default_rng(1)
    rot = np.array([B.corrected(PP[:, 0] + s_, PP[:, 1] + s_).mean() for s_ in RT.shifts(N)])
    boot = []
    for _ in range(N_BOOT_G):
        w = np.bincount(rr.integers(0, len(peaks), len(peaks)), minlength=len(peaks))
        wt = w[ia] * w[ib]
        if wt.sum():
            boot.append(np.average(pe, weights=wt))
    lo, hi = np.percentile(boot, [2.5, 97.5])
    p_mean = RT.p_value(pe.mean(), rot)                                                # exact, two-sided
    rlo, rhi = np.percentile(rot, [2.5, 97.5])
    ax.errorbar(i - 0.15, rot.mean(), yerr=[[rot.mean() - rlo], [rhi - rot.mean()]], fmt="s", color="#a9a7a1", ms=4,
                lw=1.0, capsize=2)
    ax.errorbar(i + 0.15, pe.mean(), yerr=[[pe.mean() - lo], [hi - pe.mean()]], fmt="o", color=fc(run, col), mec=col if run in OPEN else INK,
                mew=0.4 if run not in OPEN else 1.0, ms=5.5, lw=1.0, capsize=2, ecolor=INK)
    ax.text(i + 0.15, hi + 0.02, stars(p_mean), ha="center", va="bottom", fontsize=HS.BASE)
    g_rows.append((t, pe.mean(), lo, hi, rot.mean(), rlo, rhi, p_mean))
    print(f"g {t}: corrected excess {pe.mean():.3f}% ({lo:.3f}-{hi:.3f}), rotation mean {rot.mean():.3f} ({rlo:.3f}-{rhi:.3f}), P {p_mean:.4f}")
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=HS.BASE)
ax.set_xlim(-0.55, len(MAPS_B) - 0.45); ax.set_ylim(-0.12, 0.75)
divider(ax, 2.5)
ax.set_ylabel("excess at H-NS\npeak pairs (%)")
ax.set_title("peak pairs vs rotated (grey)", fontsize=HS.BASE, pad=2)
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "g", dx=-0.2)
print(pd.DataFrame(g_rows, columns=["map", "mean", "ci_lo", "ci_hi", "rot_mean", "rot_lo", "rot_hi", "p"]).round(3).to_string(index=False))

# h: ILD at H-NS peaks
ax = fig.add_subplot(bot[0, 2])
for i, (t, run, col) in enumerate(MAPS_B):
    v = pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{run}.csv").ILD_au.values
    vz = (v - v.mean()) / v.std()
    at = vz[peaks]
    boot = [rng.choice(at, len(at)).mean() for _ in range(2000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    ax.errorbar(i, at.mean(), yerr=[[at.mean() - lo], [hi - at.mean()]], fmt="o", color=fc(run, col), mec=col if run in OPEN else INK,
                mew=0.4 if run not in OPEN else 1.0, ms=5, lw=0.8, capsize=1.5, ecolor=col)
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xticks(xs); ax.set_xticklabels(labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=HS.BASE)
ax.set_ylabel("ILD z at H-NS peaks")
ax.set_xlim(-0.6, len(MAPS_B) - 0.4); ax.set_ylim(-0.7, 4.3)
divider(ax, 4.2)
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "h", dx=-0.3)

fig.canvas.draw()
for a_, b_ in pairs:                       # ILD track directly under its (square) map, bin for bin
    p, q = a_.get_position(), b_.get_position()
    b_.set_position([p.x0, p.y0 - 0.012 - q.height, p.width, q.height])
HS.check(fig)
fig.savefig(OUT + "FigureS7_timecourse.pdf", dpi=600)
fig.savefig(OUT + "FigureS7_timecourse.png", dpi=300)
print("saved")
