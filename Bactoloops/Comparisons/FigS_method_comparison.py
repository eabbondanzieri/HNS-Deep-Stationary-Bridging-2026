"""Supplemental figure: the ImageJ macro ("Bactoloop") vs published loop callers.

Same layout as the earlier "SuppFig Method comparison.pdf", with every method and the merged
datasets. Reads compare_methods.py's outputs (run that first) and its load()/chip_track().

  a     loops called on the E. coli maps: WT, delta-stpA, delta-fis, delta-ihf (96 h), WT exp and
        delta-hns 96 h
  b, c  mean ChIP-seq z-score of each loop's two anchors (H-NS / Rok); boxes: median and
        interquartile range, whiskers: full range; diamond, mean with bootstrap 95% CI; dashed,
        95th percentile of the rotation null
  d     mean H-NS ChIP z (bootstrap 95% CI) on each E. coli 96 h map
  e     replicate agreement: WT 96 h rep1 vs rep2, fraction of each replicate's top k loops (by the
        method's own score) found within 1 bin among the other's top k
  f, g  rotation test (rotation_test/): % of the loops called on the unrotated map that are called
        again (within 1 bin) when the circular map starts k bins later, E. coli WT 96 h and B. subtilis WT
  h     E. coli loops drawn on the circular chromosome with the H-NS ChIP-seq track as a grey
        ring; arc width and opacity scale with the mean ChIP z-score of the loop's two anchors
(The APA panels of earlier versions were removed on 2026-10-06; previous_versions/FigS_method_comparison_with_APA_2026-10-06.py.)

Inputs since 2026-09-29: the published callers run on raw counts (HiC-Pro allValidPairs ->
cool and bed files/raw_counts/); the macro on HiC-Pro's ICE maps of the same libraries.

The macro is shown as its medium + high confidence calls, ridge-search loops included.
Output: FigS_method_comparison.pdf / .png
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
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Wedge, PathPatch
from matplotlib.path import Path

CMP = HNS_ROOT + "Bactoloops/Comparisons/"
ns = {}
exec(open(CMP + "compare_methods.py").read().split("rows, apas, topk")[0], ns)   # load(), chip_track()
T = pd.read_csv(CMP + "method_comparison.csv")
ROT = pd.read_csv(CMP + "rotation_test/rotation_detail.csv")

BL = "Bactoloops (medium+)"
METHODS = [BL, "cooltools", "Chromosight", "HiCCUPS", "Mustache", "SIP", "Peakachu"]
NAME = {BL: "Bactoloops", **{m: m for m in METHODS[1:]}}
COL = {BL: "#e0506e", "cooltools": "#1a9e8f", "Chromosight": "#8e4fa6", "HiCCUPS": "#2a78d6",
       "Mustache": "#e0892b", "SIP": "#4d8f2f", "Peakachu": "#8a6d3b"}
SPECIES = {"wt96": ("E. coli", "H-NS, 96 h", 461, 367), "bsub": ("B. subtilis", "Rok, stationary", 422, 0)}
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
DIV = LinearSegmentedColormap.from_list("div", ["#1f5fa8", "#f0efec", "#c0392b"])
GREY = LinearSegmentedColormap.from_list("grey", ["#f4f4f4", "#101010"])
import sys as _hs
_hs.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402
HS.apply()
from matplotlib.backends.backend_pdf import PdfPages
VMAX = 1.25


def row(d, m):
    return T[(T.dataset == d) & (T.method == m)].iloc[0]


SPN = lambda sp: f"$\\it{{{sp}}}$"   # italic species name


def letter(ax, s, dx=-0.1, dy=1.02):
    HS.letter(ax, s.lower(), dx, dy)


def species_label(ax, d, x=-0.3):
    sp, prot, *_ = SPECIES[d]
    ax.text(x, 0.5, f"{SPN(sp)}\n{prot}", transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=HS.TITLE)


# two pages (house style: 500 x 600 pt, 8 pt minimum): a-d on page 1, e-h on page 2
fig = fig1 = plt.figure(figsize=HS.size(500, 600))
outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.2, 1.0], hspace=0.75, left=0.11, right=0.985, top=0.93, bottom=0.05)

# a. loop counts (values in Source Data)
ax = fig.add_subplot(outer[0])
w = 0.135
SETS = [("wt96", "WT", 1.0), ("stpA96", "Δ$\\it{stpA}$", 0.78), ("fis96", "Δ$\\it{fis}$", 0.6),
        ("ihf96", "Δ$\\it{ihfAB}$", 0.45), ("wtexp", "WT exp", 0.3), ("hns96", "Δ$\\it{hns}$", 0.16)]
for i, (d, lab, a) in enumerate(SETS):
    vals = [int(row(d, m).n_loops) for m in METHODS]
    x = np.arange(len(METHODS)) + (i - 2.5) * w
    ax.bar(x, vals, width=w * 0.9, color=[COL[m] for m in METHODS], alpha=a, edgecolor="none")
ax.set_xticks(range(len(METHODS)))
ax.set_xticklabels([NAME[m] for m in METHODS])
ax.set_ylabel("loops called\n($\\it{E. coli}$)")
ax.set_ylim(0, 195)
ax.grid(axis="y", color=GRID, lw=0.4)
ax.set_axisbelow(True)
ax.legend(handles=[Patch(color="#555555", alpha=a, label=lab) for _, lab, a in SETS], frameon=False,
          loc="lower center", ncol=6, bbox_to_anchor=(0.5, 1.0), title="96 h unless noted", handlelength=1.2, columnspacing=1.0)
letter(ax, "a", dx=-0.07, dy=1.12)

# b, c. ChIP z boxplots
gDE = outer[1].subgridspec(1, 2, wspace=0.3)
for j, (d, L) in enumerate([("wt96", "b"), ("bsub", "c")]):
    ax = fig.add_subplot(gDE[0, j])
    sp, prot, n, _ = SPECIES[d]
    z = ns["chip_track"](d, n)
    data = []
    for m in METHODS:
        P, _ = ns["load"](m, d)
        data.append(np.nanmean(z[P], axis=1))
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, whis=(0, 100),   # whiskers: full range
                    medianprops=dict(color=INK, lw=1), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    for b, m in zip(bp["boxes"], METHODS):
        b.set_facecolor(COL[m]); b.set_alpha(0.75); b.set_edgecolor(COL[m])
    for i, (s, m) in enumerate(zip(data, METHODS), start=1):
        ax.scatter(i + np.random.default_rng(i).uniform(-0.18, 0.18, len(s)), s, s=2.5, color=INK, alpha=0.35, lw=0,
                   zorder=3)
        ax.hlines(row(d, m).chip_null_95th, i - 0.35, i + 0.35, color="#9a9a9a", lw=0.9, ls=(0, (2, 1)), zorder=4)
        rw = row(d, m)                                       # mean with bootstrap 95% CI, beside the box
        ax.errorbar(i + 0.4, rw.chip_mean_z, yerr=[[rw.chip_mean_z - rw.chip_mean_z_ci_lo], [rw.chip_mean_z_ci_hi - rw.chip_mean_z]],
                    fmt="D", ms=2.2, color=INK, lw=0.8, capsize=1.2, zorder=5)
    ax.axhline(0, color=MUTED, lw=0.5, ls=":")
    ax.set_xticks(range(1, len(METHODS) + 1))
    ax.set_xticklabels([f"{NAME[m]} ({len(s)})" for m, s in zip(METHODS, data)], rotation=45,
                       ha="right", rotation_mode="anchor")
    ax.set_ylabel("mean ChIP z-score\n(two anchors)")
    ax.set_title(f"{SPN(sp)} / {prot}", fontsize=HS.TITLE, pad=3)
    ax.grid(axis="y", color=GRID, lw=0.4)
    ax.set_axisbelow(True)
    letter(ax, L, dx=-0.2, dy=1.0)

# d. mean H-NS z (bootstrap 95% CI) on every E. coli 96 h map; the WT 96 h H-NS track stands in for the mutants
ax = fig.add_subplot(outer[2])
MAPS96 = [("wt96", "WT", "o"), ("stpA96", "Δ$\\it{stpA}$", "s"), ("fis96", "Δ$\\it{fis}$", "^"), ("ihf96", "Δ$\\it{ihfAB}$", "v")]
for k, (d, lab, mk) in enumerate(MAPS96):
    for i, m in enumerate(METHODS):
        rw = row(d, m)
        if not rw.n_loops or not np.isfinite(rw.chip_mean_z):
            continue
        x = i + (k - 1.5) * 0.17
        ax.errorbar(x, rw.chip_mean_z, yerr=[[rw.chip_mean_z - rw.chip_mean_z_ci_lo], [rw.chip_mean_z_ci_hi - rw.chip_mean_z]],
                    fmt=mk, ms=3, color=COL[m], mec=INK, mew=0.3, lw=0.8, capsize=1.2)
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xticks(range(len(METHODS)))
ax.set_xticklabels([NAME[m] for m in METHODS])
ax.set_ylim(-0.8, 5.9)
ax.set_ylabel("mean H-NS ChIP z\n(± 95% CI)")
ax.grid(axis="y", color=GRID, lw=0.4)
ax.set_axisbelow(True)
ax.legend(handles=[plt.Line2D([], [], marker=mk, ls="none", color="#777777", mec=INK, mew=0.3, ms=3.5, label=lab)
                   for _, lab, mk in MAPS96], frameon=False, ncol=5, loc="upper right", handletextpad=0.1, columnspacing=1.0,
          title="$\\it{E. coli}$ 96 h maps:", alignment="left")
lg_ = ax.get_legend(); lg_._legend_box.align = "left"
letter(ax, "d", dx=-0.07, dy=1.0)

# ---------------- page 2 ----------------
fig = fig2 = plt.figure(figsize=HS.size(500, 600))
outer2 = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 1.55], hspace=0.5, left=0.11, right=0.985, top=0.955, bottom=0.03)
gFG2 = outer2[0].subgridspec(1, 2, width_ratios=[0.55, 0.45], wspace=0.05)

# e. replicate agreement, WT 96 h rep1 vs rep2 (separate libraries): each replicate's top k loops by the
# method's own score; fraction found (within 1 bin) among the other replicate's top k, averaged both ways
def _near(p, Q, n=461, tol=1):
    d = lambda a, b: np.minimum(np.abs(a - b) % n, n - np.abs(a - b) % n)
    return (((d(Q[:, 0], p[0]) <= tol) & (d(Q[:, 1], p[1]) <= tol)) | ((d(Q[:, 0], p[1]) <= tol) & (d(Q[:, 1], p[0]) <= tol))).any()


def _wilson(p, n, z=1.96):
    """95% Wilson interval for a proportion p out of n (n = k loops compared)"""
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def _found(A, B):
    return np.mean([_near(p, B) for p in A])


ax = fig.add_subplot(gFG2[0, 0])
rng_rep = np.random.default_rng(0)
for m in METHODS:
    (A, sa), (B, sb) = ns["load"](m, "wt96r1"), ns["load"](m, "wt96r2")
    kmax = min(len(A), len(B))
    if kmax < 3:
        continue
    A, B = A[np.argsort(-sa)], B[np.argsort(-sb)]
    ks = np.arange(3, kmax + 1)
    ov = [(_found(A[:k], B[:k]) + _found(B[:k], A[:k])) / 2 for k in ks]
    ax.plot(ks, ov, color=COL[m], lw=1.1, label=f"{NAME[m]} ({len(A)}/{len(B)})")
    lo, hi = _wilson(ov[-1], ks[-1])                  # 95% Wilson interval at the method's largest k
    ax.errorbar(ks[-1], ov[-1], yerr=[[ov[-1] - lo], [hi - ov[-1]]], fmt="o", ms=2.5, color=COL[m], lw=0.8, capsize=1.5)
    if m == BL:                                      # chance: rep2's top k rotated around the chromosome
        null = [np.mean([_found(A[:k], (B[:k] + s) % 461) for s in rng_rep.integers(5, 456, 50)]) for k in ks[::5]]
        ax.plot(ks[::5], null, color=MUTED, lw=0.8, ls=(0, (2, 1)), label="chance (rotated)")
ax.set_xscale("log")
ax.set_xticks([3, 10, 30, 100]); ax.set_xticklabels(["3", "10", "30", "100"])
ax.set_ylim(-0.03, 1.03)
ax.set_xlabel("top k loops per replicate")
ax.set_ylabel("shared between replicates")
ax.set_title("$\\it{E. coli}$ WT 96 h, rep1 vs rep2", fontsize=HS.TITLE, pad=3)
ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.04, 1.05), handlelength=1.2, labelspacing=0.25,
          title="loops rep1/rep2", alignment="left")
letter(ax, "e", dx=-0.15, dy=1.0)

# f, g. rotation test: the same circular map, started k bins later; % of the unrotated map's loops called again
gRot = outer2[1].subgridspec(1, 3, width_ratios=[0.38, 0.38, 0.24], wspace=0.3)
for j, (d, L) in enumerate([("wt96", "f"), ("bsub", "g")]):
    ax = fig.add_subplot(gRot[0, j])
    sp, prot, n, _ = SPECIES[d]
    R = ROT[(ROT.dataset == d) & (~ROT.missing)]
    ks = sorted(R.k.unique())
    for i, m in enumerate(METHODS):
        r = R[R.tool == NAME[m]].set_index("k").loc[ks]
        x = np.arange(len(ks)) + (i - 3) * 0.06
        y = 100 * r.kept / r.n0
        top = m == BL
        ax.plot(x, y, color=COL[m], lw=1.6 if top else 0.9, marker="o", ms=2.6 if top else 2.2,
                mec=COL[m], mew=0.3, zorder=5 if top else 3,
                label=f"{NAME[m]} (n={int(r.n0.iloc[0])})")
    ax.axvline(2.5, color=MUTED, lw=0.5, ls=":")
    ax.set_xticks(range(len(ks)))
    ax.set_xticklabels([str(k) for k in ks], rotation=90)
    ax.set_xlabel("rotation (10-kb bins)")
    ax.set_ylim(-4, 104)
    ax.set_ylabel("loops kept (%)" if j == 0 else "")
    ax.set_title(f"{SPN(sp)}", fontsize=HS.TITLE, pad=3)
    ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
    letter(ax, L, dx=-0.2 if j == 0 else -0.1, dy=1.0)
    if j == 0:
        handles = ax.get_legend_handles_labels()
lax = fig.add_subplot(gRot[0, 2]); lax.axis("off")
lax.legend(handles=[plt.Line2D([], [], color=COL[m], lw=1.2, marker="o", ms=2.5) for m in METHODS],
           labels=[NAME[m] for m in METHODS], frameon=False, loc="center left", handlelength=1.4, labelspacing=0.25,
           title="loops of the\nunrotated map\ncalled again", alignment="left")

# h. circular loop diagrams (E. coli); arcs scaled by the mean ChIP z of the two anchors
Z_LO, Z_HI = 0.0, 4.0                       # z at or below Z_LO -> thinnest/faintest, at or above Z_HI -> boldest
arc_lw = lambda f: 0.1 + 1.5 * f ** 2      # squared, so mid-range z stays light and dense panels do not fill in
arc_alpha = lambda f: 0.05 + 0.85 * f ** 2
gFG = outer2[2].subgridspec(2, 4, wspace=0.05, hspace=0.12)
R_IN, R_OUT = 1.0, 1.12
for r, (d, L) in enumerate([("wt96", "h")]):
    sp, prot, n, ori = SPECIES[d]
    z = ns["chip_track"](d, n)
    th = lambda b: np.pi / 2 - 2 * np.pi * ((b - ori) % n) / n          # oriC at top, clockwise
    for c, m in enumerate(METHODS):
        ax = fig.add_subplot(gFG[c // 4, c % 4])
        for b in range(n):                                             # ChIP ring, one wedge per bin
            v = z[b]
            col = "white" if not np.isfinite(v) else GREY(np.clip((v + 3) / 6, 0, 1))
            a1, a0 = np.degrees(th(b)), np.degrees(th(b + 1))
            ax.add_patch(Wedge((0, 0), R_OUT, a0, a1, width=R_OUT - R_IN, facecolor=col, edgecolor="none"))
        P, _ = ns["load"](m, d)
        zf = np.clip((np.nanmean(z[P], axis=1) - Z_LO) / (Z_HI - Z_LO), 0, 1) if len(P) else np.array([])
        for (a, b), f in sorted(zip(P, zf), key=lambda t: t[1]):      # weakest first, strongest on top
            p1 = np.array([np.cos(th(a + 0.5)), np.sin(th(a + 0.5))]) * R_IN
            p2 = np.array([np.cos(th(b + 0.5)), np.sin(th(b + 0.5))]) * R_IN
            sep = min(abs(b - a), n - abs(b - a)) / n                       # 0..0.5 of the circle
            mid = (p1 + p2) / 2
            f = 0.9 * (1 - 2 * sep) / max(np.hypot(*mid), 1e-6)
            ax.add_patch(PathPatch(Path([p1, mid * f, p2], [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
                                   fc="none", ec=COL[m], lw=arc_lw(f), alpha=arc_alpha(f)))
        ax.plot([0], [R_OUT + 0.07], marker="v", color="#f2c200", ms=3.2, mec=INK, mew=0.3)
        ax.text(0, R_OUT + 0.16, "ori", ha="center", va="bottom", fontsize=HS.BASE, style="italic")
        ter = 155 if n == 461 else 202      # dif on EcW3110del (1,553 kb); B. subtilis terC (~2,017 kb)
        tx, ty = np.cos(th(ter)), np.sin(th(ter))
        ax.plot([tx * (R_OUT + 0.07)], [ty * (R_OUT + 0.07)], marker="^", color="#2a78d6", ms=3.2, mec=INK, mew=0.3)
        ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.55)
        ax.set_aspect("equal"); ax.axis("off")
        ax.text(0.5, 0.0, f"{NAME[m]} ({len(P)})", transform=ax.transAxes, ha="center", va="top", fontsize=HS.BASE)
        if c == 0:
            letter(ax, L, dx=0.02, dy=0.95)

kax = fig.add_subplot(gFG[1, 3])                    # key for the arc scale
kax.axis("off")
kax.set_xlim(0, 1); kax.set_ylim(0, 1)
kax.text(0.05, 0.92, "mean H-NS z\nof the anchors", fontsize=HS.BASE, va="top")
for i, zv in enumerate([0, 1, 2, 3, 4]):
    y = 0.62 - i * 0.12
    f = (zv - Z_LO) / (Z_HI - Z_LO)
    kax.plot([0.05, 0.45], [y, y], color=INK, lw=arc_lw(f), alpha=arc_alpha(f), solid_capstyle="butt")
    kax.text(0.52, y, ("≤0" if zv == 0 else "≥4" if zv == 4 else str(zv)), fontsize=HS.BASE, va="center")

for f_ in (fig1, fig2):
    HS.check(f_)
with PdfPages(CMP + "FigS_method_comparison.pdf") as pp:      # two pages: a-d, e-h
    pp.savefig(fig1, dpi=300); pp.savefig(fig2, dpi=300)
fig1.savefig(CMP + "FigS_method_comparison.png", dpi=300)
fig2.savefig(CMP + "FigS_method_comparison_p2.png", dpi=300)
print("saved")
