"""Supplementary Figure 3 (new order): Bactoloops finds loops at the same places in every 96 h map that loops.

  a-f  Loop maps (Loop_density_smoothed_loops; medium+ loops circled) for the NAP mutants delta-stpA, delta-fis,
       delta-ihfAB (a-c) and the long-cell mutants delta-zapA, delta-ftsP, delta-envC (d-f) at 96 h, with the integrated
       loop density (ILD, median-subtracted, one shared scale) below.
  g    Loops called in the three Fig. 2 maps (WT exp, WT 96 h, delta-hns 96 h) and these six: bar = all loops, split into
       loops matching a WT 96 h loop (both ends within 2 bins = 20 kb, either orientation, circular distance; filled) and
       loops that do not (hatched); grey % = share of the 124 WT 96 h loops also called in that map.
  h    Pearson r between the ILD tracks of all nine maps (10-kb bins), ordered so the seven looping 96 h maps come
       first and WT exp and delta-hns 96 h last.
WT 96 h = average of the two replicates (as Fig. 2). Output: FigureS3_conditions.pdf / .png, s3_conditions_stats.csv
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
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = HNS_ROOT + "Bactoloops/"
OUT = ROOT + "updated figures/"
N, RC = 461, 368
tiff_of = lambda g: (g - (RC - 1)) % N
TICKS_KB = [3670, 0, 1000, 2000, 3000]
TICK_POS = [tiff_of(k // 10) for k in TICKS_KB]
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
# label, macro run, colour
PANELS = [("Δ$\\it{stpA}$", "stpA_96h", "#a31f1f"), ("Δ$\\it{fis}$", "fis_96h", "#2e6b1f"),
          ("Δ$\\it{ihfAB}$", "IHF_96h", "#1f2fa3"), ("Δ$\\it{zapA}$", "zapA_96h", "#c07a2c"),
          ("Δ$\\it{ftsP}$", "ftsP_96h", "#6a51a3"), ("Δ$\\it{envC}$", "envC_96h", "#b0306a")]
ALL = [("WT 96 h", "wt_96h_repavg", "#000000")] + PANELS + \
      [("WT exp", "wt_Exp", "#9a7b4f"), ("Δ$\\it{hns}$ 96 h", "HNSm_96h", "#1b6f6f")]
import sys as _hs_sys
_hs_sys.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()


def letter(ax, s, dx=-0.08, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=HS.LETTER, fontweight="bold", va="bottom", ha="right")


def loops(run):
    try:
        L = pd.read_csv(ROOT + f"FIJI_output/loop_calls_{run}.csv")
    except pd.errors.EmptyDataError:
        return np.zeros((0, 2), int)
    return np.c_[L.genome_bin1.round().astype(int), L.genome_bin2.round().astype(int)] % N


cd = lambda a, b: np.minimum((a - b) % N, (b - a) % N)


def matched(P, Q):
    """for each loop in P: is there a loop in Q with both ends within 2 bins (either orientation)?"""
    if len(P) == 0 or len(Q) == 0:
        return np.zeros(len(P), bool)
    return np.array([np.any(((cd(a, Q[:, 0]) <= 2) & (cd(b, Q[:, 1]) <= 2)) |
                            ((cd(a, Q[:, 1]) <= 2) & (cd(b, Q[:, 0]) <= 2))) for a, b in P])


LP = {run: loops(run) for _, run, _ in ALL}
ILD = {run: pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{run}.csv").ILD_au.values for _, run, _ in ALL}
WT = LP["wt_96h_repavg"]
rows = []
for lab, run, col in ALL:
    m = matched(LP[run], WT)
    rec = matched(WT, LP[run]).mean() if len(LP[run]) else 0.0
    rows.append(dict(map=lab.replace("$\\it{", "").replace("}$", ""), run=run, loops=len(LP[run]),
                     matching_WT=int(m.sum()), not_matching=int((~m).sum()), pct_WT_recovered=100 * rec))
S = pd.DataFrame(rows)
R = np.corrcoef(np.array([ILD[run] for _, run, _ in ALL]))
S.to_csv(OUT + "s3_conditions_stats.csv", index=False)
print(S.round(1).to_string(index=False))
print(pd.DataFrame(R, index=S["map"], columns=S["map"]).round(2).to_string())

# ---------------- figure ----------------
from matplotlib.backends.backend_pdf import PdfPages
FIG_H = 600                                                        # one page, as Xindan's layout (2026-10-11)
fig = fig1 = fig2 = plt.figure(figsize=HS.size(500, FIG_H))        # maps a-f, then g, h
ALLG = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.2], hspace=0.36, left=0.1, right=0.9, top=0.95, bottom=0.085)
outer = ALLG
ildm = {run: ILD[run] - np.median(ILD[run]) for _, run, _ in PANELS}
ymin = min(v.min() for v in ildm.values()) * 1.05
ymax = max(v.max() for v in ildm.values()) * 1.05
pairs = []
for k, (lab, run, col) in enumerate(PANELS):
    g = outer[k // 3].subgridspec(2, 3, height_ratios=[1, 0.32], hspace=0.05, wspace=0.12)
    c = k % 3
    ax = fig.add_subplot(g[0, c])
    ax.imshow(tifffile.imread(ROOT + f"FIJI_output/Loop_density_smoothed_loops_{run}.tif"), extent=[0, N, 0, N],
              origin="lower", interpolation="lanczos")
    ax.set_xticks(TICK_POS); ax.set_yticks(TICK_POS); ax.set_xticklabels([])
    ax.set_yticklabels(TICKS_KB if c == 0 else [], fontsize=HS.BASE)
    ax.set_title(f"{lab} 96 h\n{len(LP[run])} loops", fontsize=HS.TITLE, pad=3)
    if c == 0:
        ax.set_ylabel("Genome position (kb)", fontsize=HS.BASE)
    letter(ax, "abcdef"[k], dx=-0.2 if c == 0 else -0.05)
    ax2 = fig.add_subplot(g[1, c])
    ax2.plot(np.arange(N) + 0.5, ildm[run], color=col, lw=0.6)
    ax2.set_xlim(0, N); ax2.set_ylim(ymin, ymax)
    ax2.set_xticks(TICK_POS); ax2.set_xticklabels([t if t != 3670 else "" for t in TICKS_KB], fontsize=HS.BASE)
    ax2.set_xlabel("Genome position (kb)", fontsize=HS.BASE)
    if c == 0:
        ax2.set_ylabel("ILD (a.u.)", fontsize=HS.BASE)
    else:
        ax2.set_yticklabels([])
    pairs.append((ax, ax2))

bot = ALLG[2].subgridspec(1, 2, wspace=0.38, width_ratios=[0.8, 1.2])
# g: loops called, matching WT 96 h or not
ax = fig.add_subplot(bot[0, 0])
order_g = [("WT exp", "wt_Exp"), ("WT 96 h", "wt_96h_repavg"), ("Δ$\\it{hns}$ 96 h", "HNSm_96h")] + \
          [(l, r) for l, r, _ in PANELS]
colors = {run: col for _, run, col in ALL}
for i, (lab, run) in enumerate(order_g):
    r = S.set_index("run").loc[run]
    ax.bar(i, r.matching_WT, color=colors[run], width=0.7, edgecolor=INK, lw=0.3)
    ax.bar(i, r.not_matching, bottom=r.matching_WT, color="white", width=0.7, edgecolor=colors[run], lw=0.8, hatch="////")
    ax.text(i, r.loops + 3, f"{r.loops}", ha="center", va="bottom", fontsize=HS.BASE)
    if run != "wt_96h_repavg" and r.loops:
        pass
ax.axvline(2.5, color=MUTED, lw=0.5, ls=":")
ax.set_xticks(range(len(order_g)))
ax.set_xticklabels([l + ("" if "96 h" in l or "exp" in l else " 96 h") for l, _ in order_g], rotation=45, ha="right",
                   rotation_mode="anchor", fontsize=HS.BASE)
ax.set_xlim(-0.6, len(order_g) - 0.4); ax.set_ylim(0, S.loops.max() * 1.55)
ax.set_ylabel("Loops called")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#8a8a8a", edgecolor=INK, lw=0.3, label="at a WT 96 h loop"),
                   Patch(facecolor="white", edgecolor="#8a8a8a", hatch="////", lw=0.8, label="elsewhere")],
          fontsize=HS.BASE, frameon=False, loc="upper left", handlelength=1.2)
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "g", dx=-0.22)

# h: ILD correlation matrix
ax = fig.add_subplot(bot[0, 1])
cmap = LinearSegmentedColormap.from_list("r", ["#f7f6f3", "#c9d9ec", "#2a78d6", "#0b2f66"])
im = ax.imshow(R, cmap=cmap, vmin=0, vmax=1)
labs = [l for l, _, _ in ALL]
ax.set_xticks(range(len(ALL))); ax.set_yticks(range(len(ALL)))
ax.set_xticklabels(labs, rotation=45, ha="right", rotation_mode="anchor", fontsize=HS.BASE)
ax.set_yticklabels(labs, fontsize=HS.BASE)
for i in range(len(ALL)):
    for j in range(len(ALL)):
        if (i < 7 and j < 7) or i == j:
            continue                                   # 96 h block: colour only (as Xindan's layout); values in Source Data
        ax.text(j, i, f"{R[i, j]:.2f}".replace("0.", ".").replace("-.00", ".00").replace("-", "\u2212"), ha="center", va="center", fontsize=HS.BASE,
                color="white" if R[i, j] > 0.6 else INK)
ax.axhline(6.5, color=INK, lw=0.8); ax.axvline(6.5, color=INK, lw=0.8)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03, ticks=[0, 0.5, 1]); cb.set_label("Pearson r of ILD", fontsize=HS.BASE)
cb.ax.tick_params(labelsize=HS.BASE)
letter(ax, "h", dx=-0.3, dy=1.06)

fig1.canvas.draw()
for a_, b_ in pairs:                               # ILD strip directly under its map
    p, q = a_.get_position(), b_.get_position()
    b_.set_position([p.x0, p.y0 - 0.008 - q.height, p.width, q.height])
HS.check(fig1)
if True:                                                   # one page: a-h
    fig1.savefig(OUT + "FigureS3_conditions.pdf", dpi=300)          # maps resampled at 300 dpi (above their native resolution); 600 dpi gave a 12.7 MB PDF
fig1.savefig(OUT + "FigureS3_conditions.png", dpi=300)
print("saved")
