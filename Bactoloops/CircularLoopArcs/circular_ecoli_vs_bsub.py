"""
circular_ecoli_vs_bsub.py

Circular loop diagrams of our E. coli (WT, 96 h) loop calls next to the
B. subtilis (Dugar et al. 2022, WT stationary phase) loop calls, both from
HiC_PeakFindingAlgorithm_with_symmetrization.ijm, in the same style as
circular_diagram_grayscale_ring.py: anchor-protein ChIP z-score as a grayscale
ring, called loops as arcs.

Arc darkness and width are set by each loop's log2(observed/expected) over the
local background (the macro's log2_obs_exp column) on ONE shared scale for both panels, so
arcs can be compared across species: B. subtilis shows a few very strong
loops over many weak ones, E. coli many loops of similar, moderate strength.
log2_obs_exp is used instead of loop_strength because it is defined for
every loop (ridge-region loops have no loop_strength) and doesn't depend on
the macro's per-dataset intensity scaling. Loops found by the ridge search (inside
the ridge mask) are drawn like all other loops.

Arcs bend toward the center in proportion to how far apart their anchors
are, so short-range loops stay near the ring instead of every arc passing
through the center. Each panel is rotated so oriC is at the top (ter at the
bottom); tick labels are genome position in Mb.

Layout, from the center out: arcs end exactly at the ring's inner edge
(R_IN), at the center of each loop anchor's 10 kb bin; the ChIP ring is
drawn as one annular segment per ChIP bin between R_IN and R_OUT, so each
bin covers exactly its own stretch of genome; tick marks, tick labels and
the oriC/ter markers sit entirely outside the ring so they don't cover any
ChIP data.

Usage:
    python circular_ecoli_vs_bsub.py                 # E. coli = average of WT 96 h rep1 + rep2 (default)
    python circular_ecoli_vs_bsub.py --scale linear  # arcs scaled by obs/exp (fold enrichment), not log2
    python circular_ecoli_vs_bsub.py --scale relative  # obs/exp divided by that dataset's mean loop obs/exp
    python circular_ecoli_vs_bsub.py --min-confidence medium  # drop the macro's "low" confidence loops
    python circular_ecoli_vs_bsub.py --bsub run1     # B. subtilis = GEO WT run 1 only (default: runs 1 + 2 merged)
    python circular_ecoli_vs_bsub.py --metric integrated_excess  # arcs weighted by the macro's integrated_excess
    python circular_ecoli_vs_bsub.py --metric loop_strength      # arcs weighted by the macro's loop_strength
All strengths and confidences are the macro's own columns, which are measured on the
local-background-corrected log2(obs/exp) map (rolling ball 50; see the macro README).
Options combine; each non-default choice adds a suffix to the output file name.
The average (iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif)
is the mean of the two replicates' ICE matrices, which have similar depth; its
loop calls are FIJI_output/loop_calls_wt_96h_repavg.csv.
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

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.colors import to_rgba

ROOT = HNS_ROOT + "Bactoloops/"
BINSIZE = 10_000

import sys as _sys
_sys.path.insert(0, ROOT)
import excess_baseline as EB    # noqa: E402
ICE_FOR = {"FIJI_output/loop_calls_wt_96h_repavg.csv": ("cWX2200_96h_rep1rep2_average_10000_iced.tif", "ecoli_wt96_avg"),
           "FIJI_output/loop_calls_bsub_WT_st_merged.csv": ("Dugar_HiC/Dugar_WT_st_merged_iced.tif", "bsub_wt_st_merged"),
           "FIJI_output/loop_calls_rok_st.csv": ("Dugar_HiC/GSM4288936_WT_st_HiC_iced.tif", "bsub_wt_run1")}

ECOLI = {
    "average": ("E. coli WT 96 h, rep1 + rep2 average (this study)", "FIJI_output/loop_calls_wt_96h_repavg.csv", "_repavg"),
}
ap = argparse.ArgumentParser()
ap.add_argument("--ecoli", choices=list(ECOLI), default="average")
BSUB = {
    "merged": ("B. subtilis WT stationary, runs 1 + 2 merged (Dugar et al. 2022)",
               "FIJI_output/loop_calls_bsub_WT_st_merged.csv", "_bsubmerged"),
    "run1": ("B. subtilis WT stationary, run 1 (Dugar et al. 2022)", "FIJI_output/loop_calls_rok_st.csv", ""),
}
ap.add_argument("--bsub", choices=list(BSUB), default="merged")
ap.add_argument("--metric", choices=["local", "loop_strength", "integrated_excess"], default="local",
                help="arc weight: the macro's log2_obs_exp (fold over the local background; default), its "
                     "loop_strength (area x intensity; ridge-region loops, which have no loop_strength, are "
                     "left out), or its integrated_excess (excess contacts over the local background summed "
                     "over a 50 x 50 kb window, % of a typical bin's non-local contacts). Confidence is always "
                     "the macro's confidence column")
ap.add_argument("--panel", action="store_true",
                help="omit the overall title (for use as a panel in a multi-panel figure); adds _panel to the name")
ap.add_argument("--scale", choices=["log", "linear", "relative"], default="log",
                help="arc darkness/width from log2(obs/exp) (default), from obs/exp itself, or from obs/exp "
                     "relative to each dataset's own mean loop obs/exp")
ap.add_argument("--min-confidence", choices=["low", "medium", "high"], default="low",
                help="keep only loops at or above this confidence (the macro's confidence column); "
                     "'low' keeps everything")
args = ap.parse_args()
CONF_RANK = {"low": 0, "medium": 1, "high": 2}
ecoli_title, ecoli_loops, OUT_SUFFIX = ECOLI[args.ecoli]
bsub_title, bsub_loops, bsub_suffix = BSUB[args.bsub]
OUT_SUFFIX += bsub_suffix
LINEAR = args.scale == "linear"
RELATIVE = args.scale == "relative"
if args.scale != "log":
    OUT_SUFFIX += "_" + args.scale
if args.min_confidence != "low":
    OUT_SUFFIX += f"_{args.min_confidence}conf"
LOCAL = True                         # the macro's log2_obs_exp is measured over the local background
LOOPSTRENGTH = args.metric == "loop_strength"
INTEGRATED = args.metric == "integrated_excess"
OUT_SUFFIX += {"local": "_localbg", "loop_strength": "_loopstrength",
               "integrated_excess": "_integratedexcess"}[args.metric]
if args.panel:
    OUT_SUFFIX += "_panel"

PANELS = [
    dict(title=ecoli_title,
         loops=ROOT + ecoli_loops,
         chip=ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv", chip_label="H-NS ChIP-seq z-score",
         genome_len=4_607_088, ori_bp=3_671_577, ter_bp=1_553_220),       # EcW3110del: oriC between mioC and gidA; dif
    dict(title=bsub_title,
         loops=ROOT + bsub_loops,
         chip=ROOT + "Dugar et al data/rok_zscore_10kb.csv", chip_label="Rok ChIP-seq z-score",
         genome_len=4_215_606, ori_bp=0, ter_bp=2_017_000,               # 168, NC_000964.3; terC approximate
         # Rok clusters, Dugar et al. 2022 Supplementary Table 3 (RC7 and RC8 are 45 kb apart: one label),
         # plus the Rok-bound site near yonX recruited to the terminus supercluster (their Extended Data Fig. 5)
         # (label, positions, where to put the label if not at the positions' mean -- nudged away from oriC/ter)
         features=[("RC1", [548_481]), ("RC2", [2_067_800], 2_140_000), ("RC3", [2_269_071], 2_360_000), ("RC4", [3_256_123]),
                   ("RC5", [3_465_450]), ("RC6", [3_669_973]), ("RC7/8", [4_129_166, 4_174_683], 4_010_000)],
         minor_features=[("yonX site", 2_205_000)]),
]

ARC_COLOR = "crimson"
# shared scale for arc darkness/width: log2(obs/exp) 0..2, or obs/exp 1..5 (fold enrichment) with --scale linear
# or obs/exp / dataset mean 0.5..2.5 with --scale relative (1 = that dataset's average loop)
OE_MIN, OE_MAX = (1.0, 5.0) if LINEAR else (0.5, 2.5) if RELATIVE else (0.0, 2.0)
ALPHA_MIN, ALPHA_MAX = 0.08, 1.0
LW_MIN, LW_MAX = 0.8, 4.0

# ---------------- Publication styling (as circular_diagram_grayscale_ring.py) ----------------
# House style (2026-10-07): drawn at 20 in wide but printed at ~5.66 in as Supplementary Fig. 14b, so every font is
# scaled by FS_K = 20 / 5.66 to print at >= 8 pt (titles 10 pt); the tight bounding box only makes the print larger.
_sys.path.insert(0, ROOT + "updated figures/")
import house_style as HS    # noqa: E402
HS.apply()
FS_K = 20 / 5.66
SUPTITLE_FS = 10 * FS_K
PANEL_TITLE_FS = 10 * FS_K - 4          # the panel title adds 4
TICK_LABEL_FS = 8 * FS_K + 2            # some labels use TICK_LABEL_FS - 2
LANDMARK_LABEL_FS = 8 * FS_K
CBAR_LABEL_FS = 8 * FS_K
CBAR_TICK_FS = 8 * FS_K
TICK_LW = 1.4
LANDMARK_MARKER_SIZE = 12

# ---------------- Radial layout (data units; the ring's inner edge is where arcs end) ----------------
R_IN, R_OUT = 1.0, 1.09          # ChIP ring
R_TICK = (1.11, 1.15)            # tick marks, just outside the ring
R_TICK_LABEL = 1.21
R_MARKER = 1.17                  # oriC/ter markers (triangles pointing at the ring)
R_LANDMARK_LABEL = 1.31
R_FEATURE_DOT = 1.115            # Rok-cluster dots, just outside the ring
R_FEATURE_LABEL = 1.2
FEATURE_COLOR = "#d35400"        # orange, as Dugar et al. mark Rok clusters


ORI_BP = 0   # set per panel below: the genome position drawn at 12 o'clock


def bp_to_xy(bp, genome_len, radius=1.0):
    theta = np.pi / 2 - 2 * np.pi * (np.asarray(bp) - ORI_BP) / genome_len
    return radius * np.cos(theta), radius * np.sin(theta)


def arc_path(bp1, bp2, genome_len):
    """quadratic Bezier between two anchors; the control point sits on the
    bisecting radius, closer to the center the further apart the anchors are"""
    sep = abs(bp1 - bp2) % genome_len
    sep = min(sep, genome_len - sep) / genome_len          # 0..0.5 of the circle
    mid = (bp1 + bp2) / 2 if abs(bp1 - bp2) <= genome_len / 2 else (bp1 + bp2 + genome_len) / 2
    xc, yc = bp_to_xy(mid, genome_len, radius=0.9 * R_IN * (1 - 2 * sep))
    x1, y1 = bp_to_xy(bp1, genome_len, R_IN)
    x2, y2 = bp_to_xy(bp2, genome_len, R_IN)
    return Path([(x1, y1), (xc, yc), (x2, y2)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])


def arc_value(log2_oe):
    if RELATIVE:
        oe = 2.0 ** log2_oe
        return oe / oe.mean()
    return 2.0 ** log2_oe if LINEAR else log2_oe


def scale(v, lo, hi):
    return lo + (hi - lo) * np.clip((v - OE_MIN) / (OE_MAX - OE_MIN), 0, 1)


fig, axes = plt.subplots(1, 2, figsize=(20, 10.5))
chip_axes = [fig.add_axes([0.075, 0.2, 0.012, 0.5]), fig.add_axes([0.905, 0.2, 0.012, 0.5])]

for ax, cax, P in zip(axes, chip_axes, PANELS):
    G = P["genome_len"]
    ORI_BP = P["ori_bp"]

    # ChIP ring: one annular segment per genome bin, from R_IN to R_OUT. Bins
    # the ChIP track has no row for (the Rok track skips 480-490 kb), or whose
    # z-score is NaN, are drawn hatched as "no data" rather than left as a gap.
    chip = pd.read_csv(P["chip"])
    n_bins = int(np.ceil(G / BINSIZE))
    z = np.full(n_bins, np.nan)
    z[(chip["start_bp"] // BINSIZE).astype(int).values] = chip["zscore"].values
    vmax = np.nanpercentile(np.abs(z), 98)
    norm = plt.Normalize(-vmax, vmax)
    deg = lambda bp: 90 - 360 * (bp - ORI_BP) / G       # same mapping as bp_to_xy, in degrees
    starts = np.arange(n_bins) * BINSIZE
    ends = np.minimum(starts + BINSIZE, G)
    have = np.isfinite(z)
    wedges = [patches.Wedge((0, 0), R_OUT, deg(e), deg(st), width=R_OUT - R_IN)
              for st, e in zip(starts[have], ends[have])]
    colors = plt.cm.Greys(norm(z[have]))
    ax.add_collection(PatchCollection(wedges, facecolors=colors, edgecolors=colors, linewidths=0.3, zorder=5))
    for st, e in zip(starts[~have], ends[~have]):
        ax.add_patch(patches.Wedge((0, 0), R_OUT, deg(e), deg(st), width=R_OUT - R_IN, facecolor="white",
                                   edgecolor="#9a9a9a", hatch="//////", linewidth=0.4, zorder=5))
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Greys), cax=cax)
    cb.set_label(P["chip_label"], fontsize=CBAR_LABEL_FS)
    cb.ax.tick_params(labelsize=CBAR_TICK_FS)
    if cax is chip_axes[0]:
        cb.ax.yaxis.set_ticks_position("left")
        cb.ax.yaxis.set_label_position("left")

    for tick_bp in range(0, G, 1_000_000):
        gap = min(min(abs(tick_bp - lm) % G, G - abs(tick_bp - lm) % G) for lm in (P["ori_bp"], P["ter_bp"]))
        if gap >= 0.03 * G:                    # tick mark would sit under an oriC/ter marker
            xi, yi = bp_to_xy(tick_bp, G, R_TICK[0])
            xo, yo = bp_to_xy(tick_bp, G, R_TICK[1])
            ax.plot([xi, xo], [yi, yo], color="black", linewidth=TICK_LW, zorder=6)
        near_feature = any(min(abs(tick_bp - f) % G, G - abs(tick_bp - f) % G) < 0.04 * G
                           for f in [bp for feat in P.get("features", []) for bp in feat[1]])
        if gap >= 0.05 * G and not near_feature:   # tick label would collide with oriC/ter or RC labels
            xl, yl = bp_to_xy(tick_bp, G, R_TICK_LABEL)
            ax.text(xl, yl, f"{tick_bp / 1e6:.0f}", ha="center", va="center", fontsize=TICK_LABEL_FS, zorder=6)

    L = pd.read_csv(P["loops"])
    if LOOPSTRENGTH:   # arcs weighted by loop_strength; 2 ** log2_obs_exp is used below as the arc value
        L = L[L.loop_strength.notna()].copy()
        L["log2_obs_exp"] = np.log2(L["loop_strength"])
    if INTEGRATED:     # arcs weighted by baseline-corrected integrated excess (excess_baseline.py: the macro's value
        # minus the map's mean integrated excess at the same distance, so maps with different noise are comparable)
        tif, cache = ICE_FOR["FIJI_output/" + P["loops"].split("FIJI_output/")[-1]]
        L["integrated_excess"] = EB.corrected_loops(L, EB.baseline(ROOT + "iced_matrix_conversion/" + tif, cache))
        L["log2_obs_exp"] = np.log2(L["integrated_excess"].clip(lower=1e-6))
    n_all = len(L)
    L = L[L.confidence.map(CONF_RANK) >= CONF_RANK[args.min_confidence]]
    L = L.sort_values("log2_obs_exp")      # weakest first, strongest on top; the relative scale's mean is over these
    L["arc_value"] = arc_value(L["log2_obs_exp"])
    for r in L.itertuples():
        anchor_guided = str(r.note).startswith(("anchor-guided", "ridge region"))   # ridge-region loops
        rgba = to_rgba(ARC_COLOR, scale(r.arc_value, ALPHA_MIN, ALPHA_MAX))
        # + BINSIZE / 2: genome_bp is the bin's start; the arc lands on the bin's center
        ax.add_patch(patches.PathPatch(arc_path(r.genome_bp1 + BINSIZE / 2, r.genome_bp2 + BINSIZE / 2, G),
                                       facecolor="none",
                                       edgecolor=rgba, linewidth=scale(r.arc_value, LW_MIN, LW_MAX),
                                       linestyle="-", zorder=3))

    for label, bp, mcolor in [("oriC", P["ori_bp"], "gold"), ("ter", P["ter_bp"], "dodgerblue")]:
        xm, ym = bp_to_xy(bp, G, R_MARKER)
        pointing_in = deg(bp) - 90 + 180       # marker rotation so the triangle points at the ring
        ax.plot(xm, ym, marker=(3, 0, pointing_in), markersize=LANDMARK_MARKER_SIZE, color=mcolor,
                markeredgecolor="black", markeredgewidth=1.2, zorder=10)
        xl, yl = bp_to_xy(bp, G, R_LANDMARK_LABEL)
        ax.text(xl, yl, "ori" if label == "oriC" else label, ha="center", va="center", fontsize=LANDMARK_LABEL_FS, style="italic", zorder=10)

    if (~have).any():
        xl, yl = bp_to_xy(starts[~have][0] + BINSIZE / 2, G, R_TICK_LABEL + 0.02)
        ax.text(xl, yl, "no data", ha="center", va="center", fontsize=TICK_LABEL_FS - 2, color="#6b6a66",
                style="italic", zorder=6)
    for feat in P.get("features", []):
        label, bps = feat[0], feat[1]
        label_bp = feat[2] if len(feat) > 2 else np.mean(bps)
        for bp in bps:
            xd, yd = bp_to_xy(bp, G, R_FEATURE_DOT)
            ax.plot(xd, yd, "o", ms=6, color=FEATURE_COLOR, mec="white", mew=0.6, zorder=9)
        xl, yl = bp_to_xy(label_bp, G, R_FEATURE_LABEL)
        ax.text(xl, yl, label, ha="center", va="center", fontsize=TICK_LABEL_FS, color=FEATURE_COLOR,
                fontweight="bold", zorder=9)
    for label, bp in P.get("minor_features", []):
        xd, yd = bp_to_xy(bp, G, R_FEATURE_DOT)
        ax.plot(xd, yd, "o", ms=4.5, color="#8a8a8a", mec="white", mew=0.5, zorder=9)
        # label left to the legend (house style: no room at >= 8 pt next to the terminus labels)
    n_ag = L.note.astype(str).str.startswith(("anchor-guided", "ridge region")).sum()
    what = "loop_strength" if LOOPSTRENGTH else "fold over local background" if LOCAL else "obs/exp"
    if INTEGRATED:
        mean_note = (f"\nmean loop excess {np.mean(2.0 ** L.log2_obs_exp):.2f}% of anchor contacts"
                     if RELATIVE else "")
    else:
        mean_note = ("" if LOOPSTRENGTH else f"\nmean loop {what} {np.mean(2.0 ** L.log2_obs_exp):.2f}") if RELATIVE else ""
    conf_note = "" if args.min_confidence == "low" else f" of {n_all} ({args.min_confidence}+ confidence)"
    if args.panel:   # figure panel: species name only; counts and means go in the caption
        species = " ".join(P["title"].split()[:2])
        protein = P["chip_label"].split()[0]
        ax.set_title(f"$\\it{{{species}}}$ ({protein})", fontsize=PANEL_TITLE_FS + 4, pad=10)
    else:
        ax.set_title(f"{P['title']}\n{len(L)} loops{conf_note}, {n_ag} ridge-region{mean_note}",
                     fontsize=PANEL_TITLE_FS,
                     fontweight="bold", pad=14)
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")

# shared arc-strength key, centered below the panels
kax = fig.add_axes([0.38, 0.06, 0.24, 0.07])
pad = 0.05 * (OE_MAX - OE_MIN)
kax.set_xlim(OE_MIN - pad, OE_MAX + pad)
kax.set_ylim(0, 1)
half = 0.045 * (OE_MAX - OE_MIN)
for v in np.linspace(OE_MIN + (OE_MAX - OE_MIN) / 8, OE_MAX, 8):
    kax.plot([v - half, v + half], [0.5, 0.5], color=to_rgba(ARC_COLOR, scale(v, ALPHA_MIN, ALPHA_MAX)),
             linewidth=scale(v, LW_MIN, LW_MAX), solid_capstyle="butt")
kax.set_yticks([])
for side in ("left", "right", "top"):
    kax.spines[side].set_visible(False)
if RELATIVE:
    kax.set_xticks([0.5, 1, 1.5, 2, 2.5])
    kax.set_xticklabels(["\u22640.5", "1\n(dataset mean)", "1.5", "2", "\u22652.5"], fontsize=CBAR_TICK_FS)
    kax.set_xlabel(("loop_strength (area x intensity)" if LOOPSTRENGTH else
                    "loop integrated excess over local background" if INTEGRATED else
                    "loop fold over local background" if LOCAL else "loop observed/expected")
                   + " relative to that dataset's mean loop", fontsize=CBAR_LABEL_FS)
elif LINEAR:
    kax.set_xticks([1, 2, 3, 4, 5])
    kax.set_xticklabels(["1", "2", "3", "4", "≥5"], fontsize=CBAR_TICK_FS)
    kax.set_xlabel("loop observed/expected (fold enrichment), same scale in both panels", fontsize=CBAR_LABEL_FS)
else:
    kax.set_xticks([0, 0.5, 1, 1.5, 2])
    kax.set_xticklabels(["0", "0.5", "1", "1.5", "≥2"], fontsize=CBAR_TICK_FS)
    kax.set_xlabel("loop log2(observed/expected), same scale in both panels", fontsize=CBAR_LABEL_FS)

if not args.panel:
    fig.suptitle("Loop networks: E. coli vs B. subtilis (oriC at top)\n"
                 "arc darkness and width = " + ("observed/expected relative to each dataset's mean loop" if RELATIVE
                                                else "observed/expected" if LINEAR else "log2(observed/expected)") +
                 "  |  ring = anchor-protein ChIP z-score",
                 fontsize=SUPTITLE_FS, y=1.0)
fig.savefig(ROOT + f"CircularLoopArcs/circular_ecoli_vs_bsub{OUT_SUFFIX}.png", dpi=300, bbox_inches="tight")
fig.savefig(ROOT + f"CircularLoopArcs/circular_ecoli_vs_bsub{OUT_SUFFIX}.pdf", bbox_inches="tight")
print("saved png and pdf")
