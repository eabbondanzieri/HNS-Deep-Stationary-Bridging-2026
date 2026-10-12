"""
figure7.py -- Figure 7 (summary model), following Xindan's layout of 2026-10-11 (Paper Update/Main Figures/From_Xindan/
Figure_7_schematic_model3.pdf).

Top: circular loop diagrams, exponential phase -> deep stationary phase, as figure7_summary_mock.py:
  inner ring  H-NS binding (anti-mCherry ChIP-seq z, same antibody in both phases)
  outer ring  H-NS silencing relative to the genome (RNA-seq, CFU), one scale for both phases, coloured from 2-fold
  arcs        Bactoloops loops (exp gold, 96 h black), darkness/width = baseline-corrected integrated excess
Bottom: Xindan's redrawn cells (Cell cartoons/Fig7_Exp.pdf, Fig7_96h.pdf) placed as vector drawings, with the legend of
the cell components (H-NS, AT-rich DNA, RNAP, ribosome) taken from her figure.

Output: updated figures/Figure7.pdf / .png
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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.path import Path

ROOT = HNS_ROOT
B = ROOT + "Bactoloops/"
sys.path.insert(0, B)
sys.path.insert(0, B + "updated figures/")
import excess_baseline as EB    # noqa: E402
import house_style as HS        # noqa: E402

HS.apply()
OUT = B + "updated figures/Figure7"
GENOME, ORIC, TER, BIN, N = 4_607_088, 3_671_577, 1_553_220, 10_000, 461
CHIP_COL, SIL_COL, ARC_COL, NODATA = "#c0392b", "#009999", "#1a1a1a", "#d9d9d9"
EXP_COL = "#9a7b4f"              # WT exp colour of Fig. 2
ZMAX = 3.0                       # ChIP z-score at full colour
SMIN, SMAX = 1.0, 3.0            # silencing, log2 relative to genome: white up to 2-fold, full colour at 8-fold
MIN_SPAN = 0                    # draw loops spanning >= MIN_SPAN (circular distance); 0 = all loops
EX_MIN, EX_MAX = 0.0, 2.0        # loop excess (%), as Fig. 2
LW, AL = (0.15, 0.9), (0.10, 0.95)

zscore = lambda v: (v - np.nanmean(v)) / np.nanstd(v)


def chip_z(fname):
    c = pd.read_csv(B + "ChIP-seq_data/" + fname)
    z = np.full(N, np.nan)
    z[(c.start_bp // BIN).astype(int) % N] = c.zscore.values
    return zscore(z)


def silencing(tp):
    r = pd.read_csv(ROOT + f"Spike-In/genome_positions/DESeq2_apeglm_CFU_{tp}_with_positions.txt", sep="\t")
    r = r[(r.GeneID != "hns") & r.Gene_start_lab.notna() & r.log2FoldChange.notna()]
    a = r.groupby(((r.Gene_start_lab - 1) // BIN).astype(int)).log2FoldChange.agg(["mean", "size"])
    s = -a["mean"].where(a["size"] >= 3).reindex(range(N)).values      # log2(delta-hns/WT)
    return s - np.nanmean(s)


def loops(name, tif, cache):
    L = pd.read_csv(B + f"FIJI_output/loop_calls_{name}.csv")
    L["excess"] = EB.corrected_loops(L, EB.baseline(B + "iced_matrix_conversion/" + tif, cache))
    d = (L.genome_bp1 - L.genome_bp2).abs() % GENOME
    L["span"] = np.minimum(d, GENOME - d)
    L.attrs["n_all"] = len(L)
    return L[L.span >= MIN_SPAN].sort_values("excess")


PHASES = [dict(title="Exponential phase", arc=EXP_COL, chip=chip_z("zscore_mCherry_Exp_10kb.csv"), sil=silencing("Exp"),
               L=loops("wt_Exp", "cWX2200_exp_Azure_10percent_Sau3A_rep1_10000_iced.tif", "wt_exp")),
          dict(title="Deep stationary phase", arc=ARC_COL, chip=chip_z("zscore_mCherry_96h_10kb.csv"), sil=silencing("96h"),
               L=loops("wt_96h_repavg", "cWX2200_96h_rep1rep2_average_10000_iced.tif", "ecoli_wt96_avg"))]

R_IN, R_CHIP, R_SIL, R_MARK, R_LAB = 1.0, (1.0, 1.12), (1.15, 1.27), 1.36, 1.52
deg = lambda bp: 90 - 360 * (np.asarray(bp) - ORIC) / GENOME


def xy(bp, r):
    t = np.deg2rad(deg(bp))
    return r * np.cos(t), r * np.sin(t)


def arc(bp1, bp2):
    sep = abs(bp1 - bp2) % GENOME
    sep = min(sep, GENOME - sep) / GENOME
    mid = (bp1 + bp2) / 2 if abs(bp1 - bp2) <= GENOME / 2 else (bp1 + bp2 + GENOME) / 2
    return Path([xy(bp1, R_IN), xy(mid, 0.9 * R_IN * (1 - 2 * sep)), xy(bp2, R_IN)],
                [Path.MOVETO, Path.CURVE3, Path.CURVE3])


cmap = lambda col: LinearSegmentedColormap.from_list("w_" + col, ["#ffffff", col])
CM_CHIP, CM_SIL = cmap(CHIP_COL), cmap(SIL_COL)
NZ, NS = plt.Normalize(0, ZMAX), plt.Normalize(SMIN, SMAX)
fr = lambda v: float(np.clip((v - EX_MIN) / (EX_MAX - EX_MIN), 0, 1))


def draw(ax, P):
    starts = np.arange(N) * BIN
    ends = np.minimum(starts + BIN, GENOME)
    for vals, (r0, r1), cm, nm in ((P["chip"], R_CHIP, CM_CHIP, NZ), (P["sil"], R_SIL, CM_SIL, NS)):
        cols = [to_rgba(NODATA) if not np.isfinite(v) else cm(nm(max(v, nm.vmin))) for v in vals]
        w = [patches.Wedge((0, 0), r1, deg(e), deg(s), width=r1 - r0) for s, e in zip(starts, ends)]
        ax.add_collection(PatchCollection(w, facecolors=cols, edgecolors=cols, linewidths=0.15, zorder=5))
        for r in (r0, r1):
            ax.add_patch(patches.Circle((0, 0), r, fill=False, lw=0.25, ec="#8c8c8c", zorder=6))
    for row in P["L"].itertuples():
        f = fr(row.excess)
        ax.add_patch(patches.PathPatch(arc(row.genome_bp1 + BIN / 2, row.genome_bp2 + BIN / 2), fc="none",
                                       ec=to_rgba(P["arc"], AL[0] + f * (AL[1] - AL[0])),
                                       lw=LW[0] + f * (LW[1] - LW[0]), zorder=3))
    for name, bp, mcol in (("ori", ORIC, "gold"), ("ter", TER, "dodgerblue")):
        ax.plot(*xy(bp, R_MARK), marker=(3, 0, deg(bp) + 90), ms=4, color=mcol, mec="black", mew=0.4, zorder=10)
        ax.text(*xy(bp, R_LAB), name, ha="center", va="center", fontsize=HS.BASE, style="italic", zorder=10)
    lim = R_LAB + 0.12
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect("equal"); ax.axis("off")



# ---------------- layout (pt), after Xindan's figure ----------------
import pymupdf                                                 # noqa: E402  places the vector cell drawings
CART = B + "Cell cartoons/"
# legend icons: Cell cartoons/Fig7_legend_icons.pdf, clipped from Xindan's Paper Update/Main Figures/From_Xindan/Figure_7_schematic_model3.pdf
WPT, HPT = 350, 205
fig = plt.figure(figsize=(WPT / 72, HPT / 72))
A = lambda x, y, w, h: fig.add_axes([x / WPT, 1 - (y + h) / HPT, w / WPT, h / HPT])   # pt from the top left
CIRC = 104                                                     # circle panel size (pt)
for i, P in enumerate(PHASES):
    x = 10 + i * 130
    draw(A(x, 22, CIRC, CIRC), P)
    fig.text((x + CIRC / 2) / WPT, 1 - 16 / HPT, P["title"], ha="center", va="bottom", fontsize=HS.LETTER)
fig.text(127 / WPT, 1 - 74 / HPT, "→", ha="center", va="center", fontsize=16)
fig.text(126.5 / WPT, 1 - 162 / HPT, "→", ha="center", va="center", fontsize=16)
for k, (lab, col) in enumerate((("ChIP-seq", CHIP_COL), ("RNA-seq", SIL_COL))):        # ring key
    y = 62 + k * 14
    fig.add_artist(plt.Line2D([262 / WPT] * 2, [1 - (y - 5) / HPT, 1 - (y + 5) / HPT], color=col, lw=2.2))
    fig.text(269 / WPT, 1 - y / HPT, lab, ha="left", va="center", fontsize=HS.TITLE, color=col)
HS.check(fig)
fig.savefig(OUT + "_without_cells.pdf", dpi=300)

doc = pymupdf.open(OUT + "_without_cells.pdf"); pg = doc[0]
c0 = pymupdf.open(CART + "Fig7_Exp.pdf"); cw, ch = c0[0].rect.width, c0[0].rect.height
CELL_W = 114; CELL_H = CELL_W * ch / cw
for i, f in enumerate(("Fig7_Exp.pdf", "Fig7_96h.pdf")):
    x = 3 + i * 134
    pg.show_pdf_page(pymupdf.Rect(x, 136, x + CELL_W, 136 + CELL_H), pymupdf.open(CART + f), 0, keep_proportion=True)
# legend of the cell components: the icons from Xindan's figure, labels set here
ICONS = pymupdf.open(CART + "Fig7_legend_icons.pdf")           # the icons column of her figure (clipped from XW)
ICON = pymupdf.Rect(368.5, 237.0, 388.3, 297.5)                # where that column sat on her page (pt), for the label offsets
pg.show_pdf_page(pymupdf.Rect(256, 133, 256 + ICON.width, 133 + ICON.height), ICONS, 0)
ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
for k, lab in enumerate(("H-NS", "AT-rich DNA", "RNAP", "Ribosome")):
    pg.insert_text((256 + 388.5 - ICON.x0, 133 + (236.3 + 9.1 - ICON.y0) + k * 14.6), lab,   # her label offsets
                   fontsize=HS.TITLE, fontname="ArialF7", fontfile=ARIAL)
doc.save(OUT + ".pdf", garbage=3, deflate=True)
pg.get_pixmap(dpi=600).save(OUT + ".png")
print("saved", OUT + ".pdf")
