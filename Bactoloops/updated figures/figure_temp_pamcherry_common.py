"""TEMPORARY figure (2026-10-07): hns-PAmCherry (cWX2882) vs WT looping at 96 h, like for like.
All maps built from HiC-Pro valid pairs with the same steps and subsampled to a common depth (5,732,683 pairs, the smallest
library; tools/common_depth_maps.py): WT rep1 and rep2 and PAmCherry rep2 (3 subsamples each), PAmCherry rep3 (as is), and the
two-replicate averages; full-depth maps for WT rep1, rep2 and PAmCherry rep2 (all ~10.4 M pairs). Bactoloops runs in
FIJI_output/extra_runs/pamcherry_common/. (PAmCherry rep1 has no valid-pairs file and is not included.)
  a-d  loops at common depth: WT and PAmCherry two-replicate averages, one subsample of each PAmCherry replicate
  e    loops called per map (points: subsamples)
  f    call-independent loop strength: mean baseline-corrected integrated excess at the 128 WT loop positions
  g    per-loop comparison, PAmCherry average vs WT average at common depth; slope through the origin = looping retained
Output: FigureTemp_pamcherry_common.pdf / .png, temp_pamcherry_common_stats.csv, temp_pamcherry_common_retention.csv"""
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Wedge, PathPatch
from matplotlib.path import Path
from scipy.stats import wilcoxon

ROOT = HNS_ROOT + "Bactoloops/"
OUT = ROOT + "updated figures/"
CD = ROOT + "iced_matrix_conversion/common_depth/"
RUN = ROOT + "FIJI_output/extra_runs/pamcherry_common/"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "Supplemental_figure_Ecoli_vs_Bsub"); sys.path.insert(0, OUT)
import rotation_test as RT                      # noqa: E402
import excess_baseline as EB                    # noqa: E402
import house_style as HS                        # noqa: E402
HS.apply()
N, RC, ORI, TER = 461, 368, 367, 155
genome_of = lambda t: (RC - 1 + t) % N
PAM_C, WT_C = "#d1495b", "#000000"
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
GREY = LinearSegmentedColormap.from_list("grey", ["#f4f4f4", "#101010"])
EX_MIN, EX_MAX = 0.0, 3.0
rng = np.random.default_rng(0)

NAMES = ([f"{s}96_rep{r}_common{d}" for s, r in [("wt", 1), ("wt", 2), ("pam", 2)] for d in range(3)] + ["pam96_rep3_common0"] +
         ["wt96_avg2_common", "pam96_avg2_common", "wt96_rep1_full", "wt96_rep2_full", "pam96_rep2_full"])
c = pd.read_csv(ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv")
z_gen = np.full(N, np.nan); z_gen[(c.start_bp // 10_000).astype(int) % N] = c.zscore.values
z_tif = z_gen[[genome_of(t) for t in range(N)]]
cdist = lambda a, b: np.minimum((a - b) % N, (b - a) % N)
W = pd.read_csv(ROOT + "FIJI_output/loop_calls_wt_96h_repavg.csv")
WTL = np.c_[W.bin1.round().astype(int), W.bin2.round().astype(int)] % N       # the 128 WT loops (TIFF bins)
L, WL, rows = {}, {}, []
for n in NAMES:
    L[n] = pd.read_csv(RUN + f"loop_calls_{n}.csv")
    B = EB.baseline(CD + f"{n}_10000_iced.tif", "cd_" + n)
    if len(L[n]): L[n]["excess_corr"] = EB.corrected_loops(L[n], B)
    WL[n] = B.corrected(WTL[:, 0], WTL[:, 1])
    P = np.c_[L[n].bin1.round().astype(int), L[n].bin2.round().astype(int)] % N
    f = lambda Q: np.nanmean(np.nanmean(z_tif[Q], axis=1))
    rot = [f((P + s) % N) for s in RT.shifts(N)]
    ok = np.isfinite(z_tif); ild = pd.read_csv(RUN + f"integrated_loop_density_{n}.csv").ILD_au.values
    rows.append(dict(map=n, loops=len(L[n]), anchor_z=f(P), anchor_p=RT.p_value(f(P), rot),
                     ild_chip_r=np.corrcoef(ild[ok], z_tif[ok])[0, 1], wtloop_excess=WL[n].mean()))
S = pd.DataFrame(rows).set_index("map")
lib = lambda s, r: [f"{s}96_rep{r}_common{d}" for d in range(3)] if (s, r) != ("pam", 3) else ["pam96_rep3_common0"]
REPS = {("wt", 1): lib("wt", 1), ("wt", 2): lib("wt", 2), ("pam", 2): lib("pam", 2), ("pam", 3): lib("pam", 3)}
rep_wl = {k: np.mean([WL[n] for n in v], axis=0) for k, v in REPS.items()}            # per replicate, mean over subsamples
wt_ref = (rep_wl[("wt", 1)] + rep_wl[("wt", 2)]) / 2


def slope(v, ref):
    s0 = np.sum(v * ref) / np.sum(ref ** 2)
    bs = [np.sum(v[i] * ref[i]) / np.sum(ref[i] ** 2) for i in (rng.integers(0, len(v), len(v)) for _ in range(2000))]
    return s0, np.percentile(bs, 2.5), np.percentile(bs, 97.5)


ret = []
for lab, v, ref in [("PAmCherry rep2 (common depth)", rep_wl[("pam", 2)], wt_ref), ("PAmCherry rep3 (common depth)", rep_wl[("pam", 3)], wt_ref),
                    ("PAmCherry average vs WT average (common depth)", WL["pam96_avg2_common"], WL["wt96_avg2_common"]),
                    ("PAmCherry rep2 vs WT reps (full depth, ~10.4 M)", WL["pam96_rep2_full"], (WL["wt96_rep1_full"] + WL["wt96_rep2_full"]) / 2)]:
    s0, lo, hi = slope(v, ref)
    ret.append(dict(comparison=lab, slope=s0, lo=lo, hi=hi, mean_ratio=v.mean() / ref.mean(), wilcoxon_p=wilcoxon(v, ref).pvalue))
R = pd.DataFrame(ret)
rep_mean = {k: v.mean() for k, v in rep_wl.items()}
wt_m, pam_m = [rep_mean[("wt", 1)], rep_mean[("wt", 2)]], [rep_mean[("pam", 2)], rep_mean[("pam", 3)]]
print(f"per-replicate excess at WT loops (common depth): WT {wt_m}, PAmCherry {pam_m}; ratio of means {np.mean(pam_m) / np.mean(wt_m):.3f}; "
      f"range {min(pam_m) / max(wt_m):.2f}-{max(pam_m) / min(wt_m):.2f}")
pd.set_option("display.width", 200); print(S.round(3).to_string()); print(R.round(4).to_string())
S.to_csv(OUT + "temp_pamcherry_common_stats.csv"); R.to_csv(OUT + "temp_pamcherry_common_retention.csv", index=False)

# ---------------- figure ----------------
fig = plt.figure(figsize=HS.size(500, 520))
outer = fig.add_gridspec(2, 1, height_ratios=[0.8, 1], hspace=0.3, left=0.1, right=0.985, top=0.94, bottom=0.12)
top = outer[0].subgridspec(1, 4, wspace=0.05)
R_IN, R_OUT = 1.0, 1.1
th = lambda g: np.pi / 2 - 2 * np.pi * ((g - ORI) % N) / N
for n_, (n, ttl) in enumerate([("wt96_avg2_common", "WT, 2 reps"), ("pam96_avg2_common", "hns-PAmCherry,\n2 reps"),
                               ("pam96_rep2_common0", "hns-PAmCherry\nrep2"), ("pam96_rep3_common0", "hns-PAmCherry\nrep3")]):
    ax = fig.add_subplot(top[0, n_])
    for g in range(N):
        v = z_gen[g]
        ax.add_patch(Wedge((0, 0), R_OUT, np.degrees(th(g + 1)), np.degrees(th(g)), width=R_OUT - R_IN,
                           facecolor="white" if not np.isfinite(v) else GREY(np.clip((v + 2) / 6, 0, 1)), edgecolor="none", zorder=3))
    for _, r in L[n].sort_values("excess_corr").iterrows():
        a, b = int(round(r.genome_bin1)), int(round(r.genome_bin2))
        f = np.clip((r.excess_corr - EX_MIN) / (EX_MAX - EX_MIN), 0, 1)
        p1 = np.array([np.cos(th(a + 0.5)), np.sin(th(a + 0.5))]) * R_IN
        p2 = np.array([np.cos(th(b + 0.5)), np.sin(th(b + 0.5))]) * R_IN
        sep = min(abs(b - a), N - abs(b - a)) / N; mpt = (p1 + p2) / 2
        kk = 0.9 * (1 - 2 * sep) / max(np.hypot(*mpt), 1e-6)
        ax.add_patch(PathPatch(Path([p1, mpt * kk, p2], [Path.MOVETO, Path.CURVE3, Path.CURVE3]), fc="none",
                               ec="#c8102e", lw=0.3 + 1.9 * f, alpha=0.12 + 0.88 * f, capstyle="round", zorder=2))
    ax.plot([0], [R_OUT + 0.07], marker="v", color="#f2c200", ms=4, mec=INK, mew=0.3)
    ax.text(0, R_OUT + 0.16, "ori", ha="center", va="bottom", style="italic")
    tx, ty = np.cos(th(TER)) * (R_OUT + 0.07), np.sin(th(TER)) * (R_OUT + 0.07)
    ax.plot([tx], [ty], marker="^", color="#2a78d6", ms=4, mec=INK, mew=0.3)
    ax.text(tx, ty - 0.1, "ter", ha="center", va="top", style="italic")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.45, 1.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"{ttl}\n{len(L[n])} loops", fontsize=HS.BASE, pad=0)
    HS.letter(ax, "abcd"[n_], dx=0.12, dy=0.98)

bot = outer[1].subgridspec(1, 3, wspace=0.6, width_ratios=[1, 1, 1.15])
GR = [("WT\nrep1", REPS[("wt", 1)], WT_C), ("WT\nrep2", REPS[("wt", 2)], WT_C), ("PAm\nrep2", REPS[("pam", 2)], PAM_C), ("PAm\nrep3", REPS[("pam", 3)], PAM_C)]
for j, (col, ylab) in enumerate([("loops", "loops called\n(5.7 M pairs)"), ("wtloop_excess", "excess at the 128\nWT loops (%)")]):
    ax = fig.add_subplot(bot[0, j])
    for x, (lab, mem, colr) in enumerate(GR):
        v = S.loc[mem, col].values
        ax.bar(x, v.mean(), color=colr, alpha=0.8, width=0.65)
        ax.scatter(x + np.linspace(-0.12, 0.12, len(v)), v, s=9, color="white", edgecolor=INK, lw=0.4, zorder=3)
    ax.set_xticks(range(len(GR))); ax.set_xticklabels([g[0] for g in GR])
    ax.set_ylabel(ylab); ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
    HS.letter(ax, "ef"[j], dx=-0.32)
ax = fig.add_subplot(bot[0, 2])
x, y = WL["wt96_avg2_common"], WL["pam96_avg2_common"]
r0 = R.iloc[2]
ax.scatter(x, y, s=6, color=PAM_C, edgecolor=INK, lw=0.2, alpha=0.8)
lim = [min(x.min(), y.min(), 0) - 0.1, max(x.max(), y.max()) + 0.2]
ax.plot(lim, lim, color=MUTED, lw=0.6, ls=":"); ax.plot(lim, [r0.slope * v for v in lim], color=PAM_C, lw=1)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_box_aspect(1)
ax.set_xlabel("WT, 2 reps (%)"); ax.set_ylabel("hns-PAmCherry, 2 reps (%)")
ax.text(0.04, 0.96, f"slope {r0.slope:.2f}\n({r0.lo:.2f}–{r0.hi:.2f})", transform=ax.transAxes, va="top")
ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
HS.letter(ax, "g", dx=-0.3)
HS.check(fig)
fig.savefig(OUT + "FigureTemp_pamcherry_common.pdf"); fig.savefig(OUT + "FigureTemp_pamcherry_common.png", dpi=300)
print("saved")
