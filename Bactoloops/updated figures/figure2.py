"""Figure 2 from the macro's current outputs (FIJI_output/) and the ChIP-seq z-score tracks.

a-c  Loop_density_smoothed_loops for E. coli WT exponential, WT 96 h and delta-hns 96 h (loops of
     medium or high confidence circled; crosses, loops found in the ridge region; cyan outline, the
     ridge-mask region), with the integrated loop density (ILD) below, median-subtracted, on one
     shared scale. Genome position in kb in TIFF order (maps centred on ter; start at 3670 kb).
d    Genome tracks, all as z-scores: ILD (WT exp, WT 96 h) and ChIP-seq (anti-H-NS in WT 96 h;
     anti-mCherry in cWX2882 at 96 h and in exponential phase).
e, f Loops of WT exponential and WT 96 h on the circular chromosome (oriC at top); arc width and
     opacity scale with baseline-corrected integrated excess (excess_baseline.py: the macro's integrated
     excess minus the map's mean at the same distance) on one absolute scale shared by both panels and
     Fig. S3 (EX_MIN..EX_MAX, % of a typical bin's non-local contacts, clipped); grey ring, H-NS ChIP-seq z-score for the same
     growth phase (anti-mCherry, cWX2882, for exponential phase; anti-H-NS for 96 h).
g, h ILD z-score against ChIP-seq z-score, bin by bin: WT exponential against anti-mCherry exponential
     (g) and WT 96 h against anti-H-NS 96 h (h), on shared axes; Pearson r in the legend (Spearman rho printed, for the caption).
Output: Figure2.pdf / .png (in this folder).
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
from matplotlib.patches import Wedge, PathPatch
from matplotlib.path import Path
from scipy.stats import spearmanr

ROOT = HNS_ROOT + "Bactoloops/"
OUT = ROOT + "updated figures/"
CHIP = ROOT + "ChIP-seq_data/"
N, BIN_KB, RECENTER, ORI = 461, 10, 368, 367            # bins; kb per bin; TIFF recentre bin; oriC genome bin (3,671.6 kb, between mioC and gidA on the lab genome)
TER_DIF = 155                                           # dif, 1,553,220 bp on the lab genome (W3110 1,592,464 minus the 39,244 bp deletion)
TER = TER_DIF
tiff_of = lambda g: (g - (RECENTER - 1)) % N            # genome bin -> TIFF position
genome_of = lambda t: (RECENTER - 1 + t) % N            # TIFF position -> genome bin
TICKS_KB = [3670, 0, 1000, 2000, 3000]
TICK_POS = [tiff_of(kb // BIN_KB) for kb in TICKS_KB]
MAPS = [("a", "WT exp", "wt_Exp", "#9a7b4f"), ("b", "WT 96 h", "wt_96h_repavg", "#000000"),
        ("c", "Δ$\\it{hns}$ 96 h", "HNSm_96h", "#1b6f6f")]
COL_HNS, COL_MCH96, COL_MCHEXP = "#c0392b", "#2e8b3a", "#1f4fbf"
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
GREY = LinearSegmentedColormap.from_list("grey", ["#f4f4f4", "#101010"])
import sys as _s0
_s0.path.insert(0, OUT)
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()


def chip_genome(fname):
    c = pd.read_csv(CHIP + fname)
    z = np.full(N, np.nan)
    z[(c.start_bp // 10_000).astype(int).values % N] = c.zscore.values
    return z


def ild_raw(name):
    return pd.read_csv(ROOT + f"FIJI_output/integrated_loop_density_{name}.csv").ILD_au.values   # TIFF order


zscore = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
to_tiff = lambda g: g[[genome_of(t) for t in range(N)]]          # genome-order track -> TIFF order
chip = {"hns96": chip_genome("zscore_wt_96h_10kb.csv"), "mch96": chip_genome("zscore_mCherry_96h_10kb.csv"),
        "mchexp": chip_genome("zscore_mCherry_Exp_10kb.csv")}
ild_z = {name: zscore(ild_raw(name)) for _, _, name, _ in MAPS}   # TIFF order
ild_z_genome = {k: v[[tiff_of(g) for g in range(N)]] for k, v in ild_z.items()}

fig = plt.figure(figsize=HS.size(500, 600))
outer = fig.add_gridspec(4, 1, height_ratios=[1.22, 0.7, 0.85, 0.5], hspace=0.3,
                         left=0.1, right=0.985, top=0.96, bottom=0.03)


def letter(ax, s, dx=-0.04, dy=1.02):
    HS.letter(ax, s, dx, dy)


# ---------------- a-c ----------------
gA = outer[0].subgridspec(2, 3, height_ratios=[1, 0.3], hspace=0.28, wspace=0.12)
ild_med = {name: ild_raw(name) - np.median(ild_raw(name)) for _, _, name, _ in MAPS}
ymax = max(v.max() for v in ild_med.values()) * 1.05
ymin = min(v.min() for v in ild_med.values()) * 1.05
pairs = []
for c, (L, title, name, colour) in enumerate(MAPS):
    ax = fig.add_subplot(gA[0, c])
    ax.imshow(tifffile.imread(ROOT + f"FIJI_output/Loop_density_smoothed_loops_{name}.tif"),
              extent=[0, N, 0, N], origin="lower", interpolation="lanczos")
    ax.set_xticks(TICK_POS); ax.set_yticks(TICK_POS)
    ax.set_xticklabels([]); ax.set_yticklabels(TICKS_KB if c == 0 else [])
    if c == 0:
        ax.set_ylabel("Genome position (kb)")
        for t in ax.get_yticklabels():
            t.set_rotation(90); t.set_va("center")
    ax.set_title(title, fontsize=HS.TITLE, pad=4)
    letter(ax, L, dx=-0.14 if c == 0 else -0.04)
    box = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85)
    lbl = dict(style="italic", fontsize=HS.BASE, bbox=box)
    ax.text(8, N - 8, "ori", va="top", ha="left", **lbl)
    ax.text(8, tiff_of(TER), "ter", va="center", ha="left", **lbl)
    ax.text(N - 8, 8, "ori", va="bottom", ha="right", **lbl)
    ax.text(tiff_of(TER), 8, "ter", va="bottom", ha="center", **lbl)
    ax2 = fig.add_subplot(gA[1, c])
    ax2.plot(np.arange(N) + 0.5, ild_med[name], color=colour, lw=0.7)
    ax2.axhline(0, color="#9a9a9a", lw=0.4)
    ax2.set_xlim(0, N); ax2.set_ylim(ymin, ymax)
    ax2.set_xticks(TICK_POS); ax2.set_xticklabels(TICKS_KB)
    ax2.set_xlabel("Genome position (kb)")
    if c == 0:
        ax2.set_ylabel("ILD (a.u.)")
    else:
        ax2.set_yticklabels([])
    pairs.append((ax, ax2))

# ---------------- d: z-scored tracks ----------------
gD = outer[1].subgridspec(1, 3, width_ratios=[0.17, 0.77, 0.06], wspace=0.0)   # right column: scale bar
ax = fig.add_subplot(gD[0, 1])
TRACKS = [("ILD WT exp", ild_z["wt_Exp"], MAPS[0][3]), ("ILD WT 96 h", ild_z["wt_96h_repavg"], MAPS[1][3]),
          None,
          ("anti-H-NS 96 h", to_tiff(chip["hns96"]), COL_HNS),
          ("anti-mCherry 96 h", to_tiff(chip["mch96"]), COL_MCH96),
          ("anti-mCherry exp", to_tiff(chip["mchexp"]), COL_MCHEXP)]
STEP = 6.0
y0 = 0.0
yticks, ylabs = [], []
for tr in TRACKS:
    if tr is None:
        ax.axhline(y0 + STEP * 0.62, color=INK, lw=0.6)
        y0 -= STEP * 0.35
        continue
    lab, v, colour = tr
    ax.plot(np.arange(N) + 0.5, np.clip(v, -2, 8) + y0, color=colour, lw=0.7)
    ax.text(-4, y0 + 0.3, lab.replace("\n", " "), color=colour, fontsize=HS.BASE, ha="right", va="center")
    yticks.append(y0); ylabs.append(lab)
    y0 -= STEP
ax.set_yticks([])
ax.set_xlim(0, N)
ax.set_ylim(y0 + STEP - 3.0, 7.0)
ax.set_xticks(TICK_POS); ax.set_xticklabels(TICKS_KB)
ax.set_xlabel("Genome position (kb)")
ax.plot([N + 8, N + 8], [yticks[-1], yticks[-1] + 4], color=INK, lw=0.9, clip_on=False)   # scale bar, outside the box
ax.text(N + 12, yticks[-1] + 2, "4 z", fontsize=HS.BASE, va="center", ha="left", clip_on=False)
letter(ax, "d", dx=-0.22, dy=0.98)

# ---------------- e, f: circular loop diagrams; g: ILD vs ChIP ----------------
gB = outer[2].subgridspec(1, 4, wspace=0.25, width_ratios=[1, 1, 0.85, 0.85])
circ_axes = []
ROW_EH = []                          # panels e-h: letters and titles aligned on one row (after layout)
R_IN, R_OUT = 1.0, 1.1
loops = {name: pd.read_csv(ROOT + f"FIJI_output/loop_calls_{name}.csv") for name in ["wt_Exp", "wt_96h_repavg"]}
import sys as _sys
_sys.path.insert(0, ROOT)
import excess_baseline as EB        # noqa: E402  integrated excess minus the map's mean at the same distance
_ICE = {"wt_Exp": ("cWX2200_exp_Azure_10percent_Sau3A_rep1_10000_iced.tif", "wt_exp"),
        "wt_96h_repavg": ("cWX2200_96h_rep1rep2_average_10000_iced.tif", "ecoli_wt96_avg")}
for name, (tif, cache) in _ICE.items():
    loops[name]["excess_corr"] = EB.corrected_loops(loops[name], EB.baseline(ROOT + "iced_matrix_conversion/" + tif, cache))
EX_MIN, EX_MAX = 0.0, 2.0           # baseline-corrected excess (%) -> thinnest/faintest .. thickest/darkest; shared with Fig. S3
th = lambda g: np.pi / 2 - 2 * np.pi * ((g - ORI) % N) / N          # oriC at top, clockwise
for c, (L_, title, name, ring, ring_lab, colour) in enumerate(
        [("e", "WT exp", "wt_Exp", chip["mchexp"], "anti-mCherry ChIP, exp", MAPS[0][3]),
         ("f", "WT 96 h", "wt_96h_repavg", chip["hns96"], "anti-H-NS ChIP, 96 h", MAPS[1][3])]):
    ax = fig.add_subplot(gB[0, c])
    circ_axes.append(ax)
    for g in range(N):
        v = ring[g]
        col = "white" if not np.isfinite(v) else GREY(np.clip((v + 2) / 6, 0, 1))
        ax.add_patch(Wedge((0, 0), R_OUT, np.degrees(th(g + 1)), np.degrees(th(g)), width=R_OUT - R_IN,
                           facecolor=col, edgecolor="none", zorder=3))   # ring drawn over the arc ends
    Lp = loops[name]
    for _, r in Lp.sort_values("excess_corr").iterrows():
        a, b = int(round(r.genome_bin1)), int(round(r.genome_bin2))
        f = np.clip((r.excess_corr - EX_MIN) / (EX_MAX - EX_MIN), 0, 1)
        p1 = np.array([np.cos(th(a + 0.5)), np.sin(th(a + 0.5))]) * R_IN
        p2 = np.array([np.cos(th(b + 0.5)), np.sin(th(b + 0.5))]) * R_IN
        sep = min(abs(b - a), N - abs(b - a)) / N
        mid = (p1 + p2) / 2
        k = 0.9 * (1 - 2 * sep) / max(np.hypot(*mid), 1e-6)
        ax.add_patch(PathPatch(Path([p1, mid * k, p2], [Path.MOVETO, Path.CURVE3, Path.CURVE3]), fc="none",
                               ec="#c8102e", lw=0.3 + 1.9 * f, alpha=0.12 + 0.88 * f,
                               capstyle="round", joinstyle="round", zorder=2))
    anchors = np.r_[Lp.genome_bin1.round().astype(int), Lp.genome_bin2.round().astype(int)]
    zmean = np.nanmean(np.nanmean(ring[np.c_[Lp.genome_bin1.round().astype(int), Lp.genome_bin2.round().astype(int)]], axis=1))
    ax.plot([0], [R_OUT + 0.07], marker="v", color="#f2c200", ms=4, mec=INK, mew=0.3)
    ax.text(0, R_OUT + 0.16, "ori", ha="center", va="bottom", fontsize=HS.BASE, style="italic")
    tx, ty = np.cos(th(TER)) * (R_OUT + 0.07), np.sin(th(TER)) * (R_OUT + 0.07)
    ax.plot([tx], [ty], marker="^", color="#2a78d6", ms=4, mec=INK, mew=0.3)
    ax.text(tx, ty - 0.1, "ter", ha="center", va="top", fontsize=HS.BASE, style="italic")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.42, 1.6)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=HS.TITLE, pad=0)
    print(f"{title}: {len(Lp)} loops; ring {ring_lab}; mean anchor z {zmean:.2f}")
    letter(ax, L_, dx=0.08, dy=1.0); ROW_EH.append(ax)

SC = [("g", "WT exp", chip["mchexp"], ild_z_genome["wt_Exp"], MAPS[0][3], "ChIP-seq z\n(anti-mCherry)"),
      ("h", "WT 96 h", chip["hns96"], ild_z_genome["wt_96h_repavg"], MAPS[1][3], "ChIP-seq z\n(anti-H-NS)")]
allx = np.concatenate([x[np.isfinite(x)] for _, _, x, _, _, _ in SC]); ally = np.concatenate([y[np.isfinite(y)] for _, _, _, y, _, _ in SC])
xl = (allx.min() - 0.3, allx.max() + 0.3); yl = (ally.min() - 0.3, ally.max() + 0.3)     # shared axes, g and h
for c, (L, lab, x, y, colour, xlab) in enumerate(SC):
    ax = fig.add_subplot(gB[0, 2 + c])
    ok = np.isfinite(x) & np.isfinite(y)
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    ax.scatter(x[ok], y[ok], s=4, color=colour, alpha=0.55, lw=0)
    fit = np.polyfit(x[ok], y[ok], 1)
    xs = np.linspace(np.nanmin(x[ok]), np.nanmax(x[ok]), 20)
    ax.plot(xs, np.polyval(fit, xs), color=colour, lw=0.9)
    ax.text(0.05, 0.95, f"$r$ = {r:.2f}", transform=ax.transAxes, fontsize=HS.BASE, ha="left", va="top")
    ax.set_xlim(*xl); ax.set_ylim(*yl)
    ax.set_title(lab, fontsize=HS.TITLE, pad=2)
    ax.set_xlabel(xlab, fontsize=HS.BASE)
    if c == 0:
        ax.set_ylabel("ILD z-score")
    else:
        ax.set_yticklabels([])
    ax.set_box_aspect(1)          # boxed: all four spines kept
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    letter(ax, L, dx=-0.3 if c == 0 else -0.08, dy=1.0); ROW_EH.append(ax)

# ---------------- i: cartoon (moved from Fig. 1g; vector drawing Fig2i_cartoon.pdf placed into the PDF below) ----------------
gI = outer[3].subgridspec(1, 3, width_ratios=[0.95, 1, 1], wspace=0.08)
iax = fig.add_subplot(gI[0, 1:3]); iax.axis("off")
CARTOON = OUT + "Fig2i_cartoon.pdf"
letter(iax, "i", dx=0.0, dy=0.92)

# arc key for e and f, directly under them (left of the cartoon)
kax = fig.add_subplot(gI[0, 0])
kax.text(0.5, 0.97, "integrated excess\nabove baseline (e, f)", ha="center", va="top", fontsize=HS.BASE)
for i, v in enumerate([0.0, 0.5, 1.0, 1.5, 2.0]):
    f = (v - EX_MIN) / (EX_MAX - EX_MIN); y = 0.6 - i * 0.14
    kax.plot([0.22, 0.58], [y, y], color="#c8102e", lw=0.3 + 1.9 * f, alpha=0.12 + 0.88 * f, solid_capstyle="round")
    kax.text(0.64, y, {0.0: "≤0", 2.0: "≥2%"}.get(v, f"{v:g}"), ha="left", va="center", fontsize=HS.BASE)
kax.set_xlim(0, 1); kax.set_ylim(0, 1); kax.axis("off")

fig.canvas.draw()                        # align each ILD track with its (square) map, bin for bin
for ax_map, ax_ild in pairs:
    p, q = ax_map.get_position(), ax_ild.get_position()
    ax_ild.set_position([p.x0, q.y0, p.width, q.height])
fig.canvas.draw()                        # e-h: panel letters on one row (titles stay with their panels)
_inv = fig.transFigure.inverted()
_lets = [[t for t in ax.texts if t.get_fontweight() == "bold" and t.get_text() in "efgh"][0] for ax in ROW_EH]
_fy = lambda ax, art: _inv.transform(ax.transAxes.transform(art.get_position()))
y_let = max(_fy(ax, t)[1] for ax, t in zip(ROW_EH, _lets))
for ax, t in zip(ROW_EH, _lets):
    x_ = _fy(ax, t)[0]; t.set_transform(fig.transFigure); t.set_position((x_, y_let))
import pymupdf                           # noqa: E402  places the vector cartoon (panel i) into the figure PDF
cart = pymupdf.open(CARTOON); cw, ch = cart[0].rect.width, cart[0].rect.height
W, H = fig.get_size_inches() * 72
p = iax.get_position()
TITLE_PT = 14                                                     # room above the cells for the 10-pt titles
x0, x1, y_top, y_bot = p.x0 * W, p.x1 * W, (1 - p.y1) * H + TITLE_PT, (1 - p.y0) * H
s_ = min((x1 - x0) / cw, (y_bot - y_top) / ch)                    # scale, aspect kept
w_, h_ = cw * s_, ch * s_; cx_ = (x0 + x1) / 2
RECT = pymupdf.Rect(cx_ - w_ / 2, y_top, cx_ + w_ / 2, y_top + h_)
for k, ttl in enumerate(["WT 96 h", "Δ$\\it{hns}$ 96 h"]):         # one title over each cell (cells fill the left and right halves)
    fig.text((RECT.x0 + w_ * (0.25 + 0.5 * k)) / W, 1 - (RECT.y0 - 2) / H, ttl, ha="center", va="bottom", fontsize=HS.TITLE)
HS.check(fig)
fig.savefig(OUT + "_Figure2_without_cartoon.pdf")
doc = pymupdf.open(OUT + "_Figure2_without_cartoon.pdf")
doc[0].show_pdf_page(RECT, cart, 0, keep_proportion=True)
doc[0].insert_text((RECT.x0 + w_ * 0.075, RECT.y0 + h_ * 0.27), "H-NS", fontsize=HS.BASE, color=(0.18, 0.49, 0.2),
                   fontname="ArialHNS", fontfile="/System/Library/Fonts/Supplemental/Arial.ttf")   # over the cell, as in the original
doc.save(OUT + "Figure2.pdf", garbage=3, deflate=True)
doc[0].get_pixmap(dpi=300).save(OUT + "Figure2.png")
print(f"panel i: cartoon {cw:.0f} x {ch:.0f} pt placed at {w_:.0f} x {h_:.0f} pt (scale {s_:.2f})")
for lab, x, y in [("96h vs anti-H-NS 96h", chip["hns96"], ild_z_genome["wt_96h_repavg"]),
                  ("96h vs anti-mCherry 96h", chip["mch96"], ild_z_genome["wt_96h_repavg"]),
                  ("exp vs anti-mCherry exp", chip["mchexp"], ild_z_genome["wt_Exp"])]:
    ok = np.isfinite(x) & np.isfinite(y)
    print(f"ILD {lab}: Spearman rho {spearmanr(x[ok], y[ok])[0]:.2f}, Pearson r {np.corrcoef(x[ok], y[ok])[0, 1]:.2f} ({ok.sum()} bins)")
print("saved")
