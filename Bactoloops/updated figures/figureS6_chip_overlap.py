"""Supplementary Figure 6 (new data): H-NS binding at Bactoloops loop anchors, WT 96 h.

  a  Bins more than 1 SD above the mean (10-kb bins, z > 1) in the integrated loop density (ILD; WT 96 h, average of the
     two replicates) and in anti-H-NS ChIP-seq (WT 96 h): purple, both; grey, ILD only; red, ChIP only. Genome axis from
     oriC, as Fig. 2.
  b  Area-proportional Venn diagram of the bins in (a).
  c  Strongest 10-kb anti-H-NS z-score of each loop anchor (28; loop ends merged across adjacent bins) and of each H-NS-bound
     region that does not anchor a loop (ChIP z > 1 bins farther than one bin from any loop end, adjacent bins merged; 20).
  d  1-kb anti-H-NS ChIP-seq z-score (anti-H-NS / input) at two loop anchors and one non-anchor region; shaded, the
     30-kb window centred on the region's strongest 10-kb bin; red bars, 1-kb peaks (z >= 2, runs joined across a single
     sub-threshold bin).
  e, f  Per region, in that 30-kb window: span of the 1-kb peaks (first peak start to last peak end) and width of the
     widest peak. Mann-Whitney U test, two-sided.
Data: Bactoloops/anchor_clusters_1kb.py (anchor_clusters_1kb_z2.csv), anchor_enrichment.py.
Default: panels a-c only (the 1-kb analysis restates binding strength; see anchor_clusters_1kb.py). "python figureS6_chip_overlap.py alt"
draws the full a-g version as FigureS6_chip_overlap_alt_1kb.pdf / .png.
Output: FigureS6_chip_overlap.pdf / .png, s6_stats.csv
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
FULL = "alt" in sys.argv
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.stats import mannwhitneyu
from scipy.optimize import brentq

B = HNS_ROOT + "Bactoloops/"
OUT = B + "updated figures/"
N, RC = 461, 368
tiff_of = lambda g: (g - (RC - 1)) % N
TICKS_KB = [3670, 0, 1000, 2000, 3000]
TICK_POS = [tiff_of(k // 10) for k in TICKS_KB]
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
RED, GREY, BOTH = "#c0392b", "#8a8a8a", "#6a2c91"
ANC, NON = "#1f1f1e", "#d98880"
import sys as _hs_sys
_hs_sys.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()
letter = lambda ax, s, dx=-0.08, dy=1.02: ax.text(dx, dy, s, transform=ax.transAxes, fontsize=HS.LETTER, fontweight="bold",
                                                  va="bottom", ha="right")

# ---- data
c = pd.read_csv(B + "ChIP-seq_data/zscore_wt_96h_10kb.csv"); chip = np.full(N, np.nan)
chip[(c.start_bp // 10000).astype(int) % N] = c.zscore.values
t = pd.read_csv(B + "FIJI_output/integrated_loop_density_wt_96h_repavg.csv").ILD_au.values
ild = np.empty(N); ild[(367 + np.arange(N)) % N] = t; ild = (ild - ild.mean()) / ild.std()
hi_i, hi_c = ild > 1, np.nan_to_num(chip) > 1
nb, ni, nc = int((hi_i & hi_c).sum()), int(hi_i.sum()), int(hi_c.sum())
T = pd.read_csv(B + "anchor_clusters_1kb_z2.csv")
rd = lambda f: pd.read_csv(B + "ChIP-seq_data/" + f).iloc[:, 1].values.astype(float)
ch, inp = rd("anti-HNS_cWX2200_96hr_rep2_R1_clean.csv"), rd("WGS_cWX2200_96hr_rep2_R1_clean.csv")
G = len(ch); edges = np.arange(0, G, 1000); w = np.diff(np.r_[edges, G])
r = (np.add.reduceat(ch, edges) / w) / np.where(np.add.reduceat(inp, edges) > 0, np.add.reduceat(inp, edges) / w, np.nan)
z1 = (r - np.nanmean(r)) / np.nanstd(r)
above = np.nan_to_num(z1) >= 2; above[1:-1] |= above[:-2] & above[2:]
st = np.where(above & ~np.r_[False, above[:-1]])[0]; en = np.where(above & ~np.r_[above[1:], False])[0] + 1

if FULL:
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(3, 6, height_ratios=[0.75, 1, 1], hspace=0.65, wspace=1.1, left=0.11, right=0.98, top=0.96, bottom=0.08)
else:
    fig = plt.figure(figsize=HS.size(500, 290))
    gs = fig.add_gridspec(2, 6, height_ratios=[0.75, 1.25], hspace=0.7, wspace=1.1, left=0.11, right=0.98, top=0.93, bottom=0.13)

# a: binary tracks
ax = fig.add_subplot(gs[0, 0:4] if FULL else gs[0, 0:6])
for row, (mask, other, lab) in enumerate([(hi_i, hi_c, "Hi-C ILD\nWT 96 h"), (hi_c, hi_i, "ChIP H-NS\nWT 96 h")]):
    y0 = 1 - row
    for g in np.where(mask)[0]:
        col = BOTH if other[g] else (GREY if row == 0 else RED)
        ax.add_patch(plt.Rectangle((tiff_of(g), y0 + 0.08), 1, 0.84, color=col, lw=0))
    ax.text(-8, y0 + 0.5, lab, ha="right", va="center", fontsize=HS.BASE)
ax.axhline(1, color=INK, lw=0.5)
ax.set_xlim(0, N); ax.set_ylim(0, 2); ax.set_yticks([])
ax.set_xticks(TICK_POS); ax.set_xticklabels(TICKS_KB); ax.set_xlabel("Genome position (kb)")
letter(ax, "a", dx=-0.13 if FULL else -0.09)

# b: area-proportional Venn
ax = fig.add_subplot(gs[0, 4:6] if FULL else gs[1, 0:3]); ax.set_aspect("equal")
ra, rb = np.sqrt(nc / np.pi), np.sqrt(ni / np.pi)
def lens(d):
    if d >= ra + rb: return 0.0
    if d <= abs(ra - rb): return np.pi * min(ra, rb) ** 2
    a1 = ra ** 2 * np.arccos((d * d + ra * ra - rb * rb) / (2 * d * ra)); a2 = rb ** 2 * np.arccos((d * d + rb * rb - ra * ra) / (2 * d * rb))
    return a1 + a2 - 0.5 * np.sqrt((-d + ra + rb) * (d + ra - rb) * (d - ra + rb) * (d + ra + rb))
d = brentq(lambda x: lens(x) - nb, abs(ra - rb) + 1e-9, ra + rb - 1e-9)
RED_T, GREY_T, BOTH_T = "#f3c9c4", "#dcdcdc", "#d2bfe0"           # tints of the panel-a colours
ca = Circle((0, 0), ra, facecolor=RED_T, edgecolor="none"); ax.add_patch(ca)
cb = Circle((d, 0), rb, facecolor=GREY_T, edgecolor="none"); ax.add_patch(cb)
lensp = Circle((d, 0), rb, facecolor=BOTH_T, edgecolor="none"); ax.add_patch(lensp); lensp.set_clip_path(ca)
ax.add_patch(Circle((0, 0), ra, facecolor="none", edgecolor=RED, lw=0.9))
ax.add_patch(Circle((d, 0), rb, facecolor="none", edgecolor="#555555", lw=0.9))
xl, xr = d - rb, ra                                                   # lens edges on the axis
ax.text((-ra + xl) / 2, 0, f"{nc - nb}", ha="center", va="center", fontsize=HS.LETTER, fontweight="bold", color=RED)
ax.text((xl + xr) / 2, 0, f"{nb}", ha="center", va="center", fontsize=HS.LETTER, fontweight="bold", color=BOTH)
xo = (xr + d + rb) / 2                                                # ILD-only crescent: label outside with a leader
ax.plot([xo, d + rb + 0.9], [0, 0], color="#555555", lw=0.5)
ax.text(d + rb + 1.1, 0, f"{ni - nb}", ha="left", va="center", fontsize=HS.LETTER, fontweight="bold", color="#555555")
ax.text(-ra * 0.55, ra * 0.98, f"H-NS ChIP-seq\n{nc} bins", ha="right", va="bottom", fontsize=HS.BASE, color=RED)
ax.text(d + rb * 0.55, rb * 1.02, f"Hi-C ILD\n{ni} bins", ha="left", va="bottom", fontsize=HS.BASE, color="#555555")
ax.set_xlim(-ra - 3.5, d + rb + 3.5); ax.set_ylim(-ra - 0.3, ra + 2.6); ax.axis("off")
letter(ax, "b", dx=0.02 if FULL else 0.12)

# c: binding strength, anchors vs non-anchors
def strip(ax, col, ylab, fmt):
    a_, n_ = T[T.kind == "anchor"][col].values, T[T.kind == "non-anchor"][col].values
    for i, (v, cc) in enumerate([(a_, ANC), (n_, NON)]):
        ax.boxplot([v], positions=[i], widths=0.55, showfliers=False, whis=(0, 100), patch_artist=True,
                   boxprops=dict(facecolor=cc, alpha=0.25, edgecolor=cc), medianprops=dict(color=cc, lw=1.2),
                   whiskerprops=dict(color=cc, lw=0.7), capprops=dict(color=cc, lw=0.7))
        ax.scatter(i + np.random.default_rng(i).uniform(-0.14, 0.14, len(v)), v, s=6, color=cc, lw=0, alpha=0.8, zorder=3)
    p = mannwhitneyu(a_, n_).pvalue
    top = max(a_.max(), n_.max()); ax.plot([0, 0, 1, 1], [top * 1.06, top * 1.1, top * 1.1, top * 1.06], color=INK, lw=0.6)
    ax.text(0.5, top * 1.12, f"P = {p:.2g}" if p >= 0.001 else "P < 0.001", ha="center", va="bottom", fontsize=HS.BASE)
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"anchors\n(n = {len(a_)})", f"bound,\nno loops\n(n = {len(n_)})"], fontsize=HS.BASE)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, top * 1.3); ax.set_ylabel(ylab)
    ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
    return a_.mean(), n_.mean(), p
S = {}
ax = fig.add_subplot(gs[1, 0:2] if FULL else gs[1, 4:6]); S["max_z10"] = strip(ax, "max_z10", "strongest 10-kb H-NS\nChIP z-score", "{:.1f}"); letter(ax, "c", dx=-0.3)

if not FULL:
    pd.DataFrame([dict(measure="max_z10", anchors=S["max_z10"][0], non_anchor=S["max_z10"][1], p=S["max_z10"][2]),
                  dict(measure="bins_ILD_z1", anchors=ni), dict(measure="bins_ChIP_z1", anchors=nc), dict(measure="bins_both", anchors=nb)]
                 ).to_csv(OUT + "s6_stats.csv", index=False)
    print(f"bins: ILD {ni}, ChIP {nc}, both {nb}"); print(S)
    HS.check(fig)
    fig.savefig(OUT + "FigureS6_chip_overlap.pdf", dpi=600); fig.savefig(OUT + "FigureS6_chip_overlap.png", dpi=300)
    print("saved a-c"); sys.exit()
# d: 1-kb examples
EX = [("1200-1230", "anchor, 1,210–1,230 kb"), ("2730-2760", "anchor, 2,740–2,750 kb"), ("4530-4560", "no loops, 4,540–4,550 kb")]
sub = gs[1, 2:6].subgridspec(1, 3, wspace=0.12)
for k, (win, lab) in enumerate(EX):
    ax = fig.add_subplot(sub[0, k])
    lo, hi = map(int, win.split("-")); x0, x1 = lo - 20, hi + 20
    xs = np.arange(x0, x1)
    ax.axvspan(lo, hi, color=GRID, lw=0)
    ax.plot(xs + 0.5, z1[x0:x1], color=RED, lw=0.7)
    for s_, e_ in zip(st, en):
        if e_ > x0 and s_ < x1: ax.plot([max(s_, x0), min(e_, x1)], [-0.6, -0.6], color=RED, lw=2.2, solid_capstyle="butt")
    ax.axhline(2, color=MUTED, lw=0.5, ls=":")
    ax.set_xlim(x0, x1); ax.set_ylim(-1, 6)
    ax.set_xticks([x0 + 10, (x0 + x1) // 2, x1 - 10]); ax.set_xticklabels([f"{v:,}" for v in (x0 + 10, (x0 + x1) // 2, x1 - 10)], fontsize=HS.BASE)
    ax.set_title(lab, fontsize=HS.BASE, pad=2, color=ANC if k < 2 else RED)
    if k == 0: ax.set_ylabel("1-kb H-NS ChIP z-score")
    else: ax.set_yticklabels([])
    if k == 1: ax.set_xlabel("Genome position (kb)")
    if k == 0: letter(ax, "d", dx=-0.25)

# e, f: span and widest peak
ax = fig.add_subplot(gs[2, 0:2]); S["span_kb"] = strip(ax, "span_kb", "span of 1-kb peaks (kb)", "{:.0f}"); letter(ax, "e", dx=-0.3)
ax = fig.add_subplot(gs[2, 2:4]); S["widest_kb"] = strip(ax, "widest_kb", "widest 1-kb peak (kb)", "{:.0f}"); letter(ax, "f", dx=-0.3)
ax = fig.add_subplot(gs[2, 4:6]); S["n_peaks"] = strip(ax, "n_peaks", "1-kb peaks per region", "{:.0f}"); letter(ax, "g", dx=-0.3)

pd.DataFrame([dict(measure=k, anchors=v[0], non_anchor=v[1], p=v[2]) for k, v in S.items()] +
             [dict(measure="bins_ILD_z1", anchors=ni), dict(measure="bins_ChIP_z1", anchors=nc), dict(measure="bins_both", anchors=nb)]
             ).to_csv(OUT + "s6_stats_alt_1kb.csv", index=False)
print(f"bins: ILD {ni}, ChIP {nc}, both {nb} ({100*nb/ni:.0f}% of ILD, {100*nb/nc:.0f}% of ChIP)"); print(S)
fig.savefig(OUT + "FigureS6_chip_overlap_alt_1kb.pdf", dpi=600); fig.savefig(OUT + "FigureS6_chip_overlap_alt_1kb.png", dpi=300)
print("saved alt a-g")
