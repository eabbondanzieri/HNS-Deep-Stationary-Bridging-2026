"""Figure 3: H-NS hairpins persist into deep stationary phase, and the hairpin domains at loop anchors compact as loops form.
(Panels i-m of the previous Supplementary Fig. 9, moved to the main text on 2026-10-06.)

  a  genome tracks (z-scores, 10-kb bins): WT 96 h ILD, anti-H-NS ChIP-seq (96 h), and the density (bp per bin) of the
     exponential-phase CHINs, CHIDs and OPCIDs of Gavrilov, Shamovsky et al. 2025 (Gavrilov2025_MicroC/opcidsRenamed.xlsx, lab
     genome); correlations in ../chin_chid_vs_loops.py
  b  aggregate contact maps (log2 O/E, balanced 500-bp maps from full valid pairs without the 5-kb cutoff) around CHINs >= 2 kb,
     centred on the CHIN centre and rescaled by CHIN length (ends at -1, +1; Gavrilov2025_MicroC/scaled_pileups.py), averaged over
     the four orientations, 2-pixel Gaussian blur, each map relative to its own mean; dashed, anti-diagonal (hairpin stem);
     WT exp, WT 96 h and dstpA 96 h pool two libraries
  c  CHIN hairpin test on per-CHIN rescaled maps (Gavrilov2025_MicroC/scaled_stem_test.py); median over CHINs, bootstrap 95% CI,
     Wilcoxon P; brackets, paired WT vs dhns comparisons on the same CHINs; circles, two libraries pooled; diamonds, one library
  d  CHIDs: within-domain O/E vs same-size neighbouring windows (median over CHIDs; s9_fullpairs_stats.py)
  e  control: OPCIDs (transcription domains), same measure
Output: Figure3_hairpins.pdf / .png
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
_sys.path.insert(0, OUT)
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()
T = pd.read_csv(OUT + "hns_timecourse_numbers.csv").set_index("map")
PR = np.load(OUT + "_hns_timecourse_pairs.npz")
peaks, PP, CP, z = PR["peaks"], PR["PP"], PR["CP"], PR["z_tiff"]
rng = np.random.default_rng(0)
stars = RT.stars


def letter(ax, s, dx=-0.08, dy=1.02):
    HS.letter(ax, s, dx, dy)




def letter(ax, s, dx=-0.08, dy=1.02):
    HS.letter(ax, s, dx, dy)


FIG_H = 600                                   # pt; key above a (2026-10-11) and panel f, cells (Xindan, 2026-10-12)
F_PT = 64                                     # height of panel f at the bottom
fig = plt.figure(figsize=HS.size(500, FIG_H))
outer = fig.add_gridspec(3, 1, height_ratios=[0.72, 0.5, 1.0], hspace=0.42, left=0.1, right=0.985, top=1 - 43 / FIG_H, bottom=(60 + F_PT) / FIG_H)
_TR, _PG, _ROW = 0, 1, 2
# a: exponential-phase CHIDs / CHINs vs deep-stationary loops and H-NS
X = pd.read_excel(ROOT + "Gavrilov2025_MicroC/opcidsRenamed.xlsx", 0)
X.loc[X.OPCID_ID == "OPCID_40", "start_new"] = 2_522_632; X.loc[X.OPCID_ID == "CHID_20", "end_new"] = 3_950_261   # sheet typos
X = X.dropna(subset=["start_new", "end_new"]); X["type"] = X.OPCID_ID.str.split("_").str[0]
X["s"], X["e"] = X[["start_new", "end_new"]].min(1).astype(int), X[["start_new", "end_new"]].max(1).astype(int)
def _dens(df):
    v = np.zeros(N)
    for s_, e_ in zip(df.s, df.e):
        for b in range(s_ // 10000, e_ // 10000 + 1):
            v[b % N] += max(0, min(e_, (b + 1) * 10000) - max(s_, b * 10000))
    return v
_z = lambda v: (v - v.mean()) / v.std()
_ild = pd.read_csv(ROOT + "FIJI_output/integrated_loop_density_wt_96h_repavg.csv").ILD_au.values       # TIFF order
_c = pd.read_csv(ROOT + "ChIP-seq_data/zscore_wt_96h_10kb.csv"); _chip = np.zeros(N)
_chip[(_c.start_bp // 10000).astype(int) % N] = np.nan_to_num(_c.zscore.values)
to_tiff = lambda g_arr: g_arr[(RC - 1 + np.arange(N)) % N]                                           # genome order -> TIFF order
TRK = [("Hi-C ILD, WT 96 h", _z(_ild), "#000000"), ("H-NS ChIP-seq, 96 h", to_tiff(_z(_chip)), "#c0392b"),
       ("CHIN density, exp", to_tiff(_z(_dens(X[X.type == "CHIN"]))), "#2e7d32"),
       ("CHID density, exp", to_tiff(_z(_dens(X[X.type == "CHID"]))), "#6a2c91"),
       ("OPCID density, exp", to_tiff(_z(_dens(X[X.type == "OPCID"]))), "#d4780f")]      # in the order the text discusses them
ax = fig.add_subplot(outer[_TR])
STEP = 7
for k, (lab, v, col) in enumerate(TRK):
    ax.plot(np.arange(N) + 0.5, np.clip(v, -2, 9) - k * STEP, color=col, lw=0.6)
    ax.text(-6, -k * STEP + 1, lab, color=col, fontsize=HS.BASE, ha="right", va="center")
NT = len(TRK)
ax.set_xlim(0, N); ax.set_ylim(-(NT - 1) * STEP - 3, 9.5); ax.set_yticks([])
ax.set_xticks(TICK_POS); ax.set_xticklabels([3670, 0, 1000, 2000, 3000]); ax.set_xlabel("Genome position (kb)", fontsize=HS.BASE)
ax.plot([N + 4, N + 4], [-(NT - 1) * STEP, -(NT - 1) * STEP + 5], color="#1f1f1e", lw=0.9, clip_on=False)
ax.text(N + 7, -(NT - 1) * STEP + 2.5, "5 z", fontsize=HS.BASE, va="center", ha="left")
fig.canvas.draw()
p_ = ax.get_position(); ax.set_position([0.25, p_.y0, 0.69, p_.height])                          # room for the track labels
letter(ax, "a", dx=-0.3)
# key above panel a: a small sketch of each element type, its name and median length (Xindan, 2026-10-11)
fig.canvas.draw()
_pa = ax.get_position(); _y = _pa.y1 + 6 / FIG_H                    # baseline of the key, just above the tracks
_xs = np.linspace(-0.5, 0.5, 200); _g = lambda c, w: np.exp(-((_xs - c) / w) ** 2)
_t = np.linspace(0, 4 * np.pi, 400)
ICONS = {"CHIN": (_xs, _g(0, 0.09)),                                                    # one hairpin
         "CHID": (_xs, 0.55 * _g(-0.3, 0.07) + 1.0 * _g(0, 0.07) + 0.7 * _g(0.3, 0.07)),  # a cluster of hairpins
         "OPCID": (np.r_[np.linspace(-0.55, -0.2, 20), 0.2 * np.cos(_t) + 0.09 * np.cos(4.3 * _t) + 0.02 * _t - 0.12,
                         np.linspace(0.25, 0.55, 20)],                                     # a tangle with two tails
                   np.r_[np.full(20, 0.1), 0.38 * np.sin(1.3 * _t) + 0.16 * np.sin(3.7 * _t) + 0.5, np.full(20, 0.1)])}
KEY = [("CHIN", "CHIN (hairpin)", "#2e7d32", 0.0), ("CHID", "CHID (hairpin domain)", "#6a2c91", 0.33),
       ("OPCID", "OPCID (highly transcribed operons)", "#d4780f", 0.66)]
_x0 = 0.1                                                         # left edge of the track labels
KEY_SLOTS = []
for typ, name, col, fx in KEY:
    med = (X[X.type == typ].e - X[X.type == typ].s).median() / 1000
    xt = _x0 + fx * (_pa.x1 - _x0)
    t_ = fig.text(xt, _y, f"{name}\nmedian {med:.1f} kb", ha="left", va="bottom", fontsize=HS.BASE, color=col, linespacing=1.15)
    fig.canvas.draw()
    bb = t_.get_window_extent().transformed(fig.transFigure.inverted())
    iw = min(0.075, 0.985 - (bb.x1 + 0.008))                      # keep the sketch inside the figure
    iax = fig.add_axes([bb.x1 + 0.008, _y - 2 / FIG_H, iw, 30 / FIG_H]); iax.axis("off")
    KEY_SLOTS.append((typ, bb.x1 + 0.008, _y - 2 / FIG_H, iw, 30 / FIG_H))     # Xindan's vector drawing goes here (below)
_la = [t for t in ax.texts if t.get_text() == "a" and t.get_fontweight() == "bold"][0]   # letter a: up to the key row
_lx = fig.transFigure.inverted().transform(ax.transAxes.transform(_la.get_position()))[0]   # keep its x
_la.set_transform(fig.transFigure); _la.set_position((_lx, _y + 24 / FIG_H)); _la.set_va("top")

# b-e: short-range H-NS structure from full valid pairs (no 5-kb cutoff); replicates pooled (Gavrilov2025_MicroC/s9_fullpairs_stats.py)
NUD = ROOT + "Gavrilov2025_MicroC/"
PZ = np.load(NUD + "scaled_pileups.npz"); SS = pd.read_csv(NUD + "s9_fullpairs_stats.csv").set_index("condition")
from scipy.ndimage import gaussian_filter
def _blur(M, sg=2):                                    # 2-pixel Gaussian, NaN-aware (blurred values / blurred weights)
    w = np.isfinite(M).astype(float); v = np.where(w > 0, M, 0)
    num, den = gaussian_filter(v, sg, mode="constant"), gaussian_filter(w, sg, mode="constant")
    out = num / np.where(den > 0, den, np.nan); out[w == 0] = np.nan; return out
def _fourfold(M):                                      # average over both diagonal reflections and the 180-degree rotation
    with np.errstate(all="ignore"):
        return np.nanmean(np.stack([M, M.T, M[::-1, ::-1].T, M[::-1, ::-1]]), axis=0)
CND = [("WT exp", "WT exp", "#9a7b4f"), ("dhns exp", "Δ$\\it{hns}$ exp", "#7fb8b8"), ("WT 24 h", "WT 24 h", "#5a4a32"),
       ("WT 96 h", "WT 96 h", "#000000"), ("dhns 96 h", "Δ$\\it{hns}$ 96 h", "#1b6f6f"), ("dstpA 96 h", "Δ$\\it{stpA}$ 96 h", "#a31f1f")]
PILES = CND
SC = pd.read_csv(NUD + "scaled_stem_test.csv").set_index("test")      # rescaled per-CHIN stem / flank test (scaled_stem_test.py)
for k_ in SS.index:
    SS.loc[k_, "chin_log2"], SS.loc[k_, "chin_lo"], SS.loc[k_, "chin_hi"], SS.loc[k_, "chin_p"] = SC.loc[k_, ["median_log2", "lo", "hi", "p"]].values
pg = outer[_PG].subgridspec(1, 9, wspace=0.12, width_ratios=[0.3, 1, 1, 1, 1, 1, 1, 0.09, 0.3])
row = outer[_ROW].subgridspec(1, 3, wspace=0.62)
paxes = []
for c, (key, lab, col) in enumerate(PILES):
    ax = fig.add_subplot(pg[0, c + 1])
    with np.errstate(divide="ignore"): M = _blur(_fourfold(np.log2(PZ[key])))
    M = M - np.nanmean(M)                                   # each map relative to its own mean (local contrast)
    # Hi-C orientation as in the other figures: position increases upwards, main diagonal from bottom left to top right
    # (the map is four-fold symmetric, so this vertical flip equals the requested 90-degree rotation); the dashed
    # hairpin stem (anti-diagonal) then runs from top left to bottom right
    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.4, vmax=0.4, extent=[-3.075, 3.075, -3.075, 3.075], origin="lower", interpolation="nearest")
    ax.plot([-3, 3], [3, -3], color="#1f1f1e", lw=0.6, ls=(0, (2, 1.5)))
    for v_ in (-1, 1): ax.axvline(v_, color="#1f1f1e", lw=0.35, ls=":"); ax.axhline(v_, color="#1f1f1e", lw=0.35, ls=":")
    ax.set_title(lab, fontsize=HS.TITLE, pad=3, color=col)
    ax.set_xticks([-2, 0, 2]); ax.set_yticks([-2, 0, 2]); ax.tick_params(pad=1); ax.tick_params(labelsize=HS.BASE); paxes.append(ax)
    if c == 0: ax.set_ylabel("position (CHIN\nhalf-lengths)", fontsize=HS.BASE); letter(ax, "b", dx=-0.6)
    else: ax.set_yticklabels([])
cax = fig.add_subplot(pg[0, 7]); cb = fig.colorbar(im, cax=cax, ticks=[-0.4, 0, 0.4]); cb.ax.tick_params(labelsize=HS.BASE)
cax.set_title("log$_2$\nO/E", fontsize=HS.BASE, pad=3)
CBFIX = cax
paxes[2].text(1.08, -0.3, "position relative to the CHIN center (CHIN half-lengths)", transform=paxes[2].transAxes, ha="center", va="top", fontsize=HS.BASE)
def dots(ax, pre, ylab, ylim, yt):
    for g, (key, lab, col) in enumerate(CND):
        r = SS.loc[key]; m, lo, hi = 2 ** r[pre + "_log2"], 2 ** r[pre + "_lo"], 2 ** r[pre + "_hi"]
        ax.errorbar(g, m, yerr=[[m - lo], [hi - m]], fmt="o" if r.libraries > 1 else "D", ms=3.8 if r.libraries > 1 else 3.2,
                    color=col, mec="#1f1f1e", mew=0.4, lw=0.8, capsize=1.2, zorder=3)
        p = r[pre + "_p"]; star = "**" if p < 0.01 else ("*" if p < 0.05 else "n.s.")
        ax.text(g, hi * 1.03, star, ha="center", va="bottom", fontsize=HS.BASE)
    ax.axhline(1, color="#6b6a66", lw=0.5, ls=":")
    ax.set_yscale("log", base=2); ax.set_yticks(yt); ax.set_yticklabels([f"{v:g}" for v in yt]); ax.minorticks_off()
    ax.set_xticks(range(len(CND))); ax.set_xticklabels([c[1] for c in CND], rotation=45, ha="right", rotation_mode="anchor", fontsize=HS.BASE)
    ax.set_xlim(-0.6, len(CND) - 0.4); ax.set_ylim(*ylim); ax.set_ylabel(ylab, fontsize=HS.BASE); ax.tick_params(labelsize=HS.BASE)
    ax.grid(axis="y", color="#e4e2dd", lw=0.4); ax.set_axisbelow(True)
ax = fig.add_subplot(row[0, 0])
dots(ax, "chin", "stem / flank\ncontacts", (0.78, 1.95), [0.8, 1, 1.2, 1.4, 1.6])
for (g1, g2, key) in [(0, 1, "WT exp vs dhns exp"), (3, 4, "WT 96 h vs dhns 96 h")]:    # paired comparisons on the same CHINs
    p_ = SC.loc[key, "p"]; y_ = 1.7 if g1 == 0 else 1.62
    ax.plot([g1, g1, g2, g2], [y_ * 0.985, y_, y_, y_ * 0.985], color="#1f1f1e", lw=0.6)
    ax.text((g1 + g2) / 2, y_ * 1.01, "**" if p_ < 0.01 else ("*" if p_ < 0.05 else "n.s."), ha="center", va="bottom", fontsize=HS.BASE)
ax.set_title("CHINs", fontsize=HS.BASE, pad=3)
letter(ax, "c", dx=-0.36)
ax = fig.add_subplot(row[0, 1])
dots(ax, "chid", "within / neighboring\nwindows", (0.75, 1.6), [0.8, 1, 1.2, 1.4])
ax.set_title("CHIDs", fontsize=HS.BASE, pad=3)
letter(ax, "d", dx=-0.42)
ax = fig.add_subplot(row[0, 2])
dots(ax, "opcid", "within / neighboring\nwindows", (0.85, 1.3), [0.9, 1, 1.1, 1.2])
ax.set_title("OPCIDs", fontsize=HS.BASE, pad=3)
letter(ax, "e", dx=-0.42)
# f: hairpins persist (cells drawn by Xindan; placed below as vector art). Positions in pt from the top-left of the figure,
# taken from her layout (Paper Update/Figure_3_HNS_hairpins4.pdf, shifted by her page offset)
F_TOP = FIG_H - 58                                 # top of the cells
CELLS = [("WT exp", "Fig3_exp.pdf", 61), ("WT 96 h", "Fig3_96h.pdf", 243)]
for ttl, _, x in CELLS:
    fig.text((x + 51) / 500, 1 - (F_TOP - 3) / FIG_H, ttl, ha="center", va="bottom", fontsize=HS.TITLE)
fig.text(200 / 500, 1 - (F_TOP + 7) / FIG_H, "hairpins persist", ha="center", va="top", fontsize=HS.TITLE)
fig.text(200 / 500, 1 - (F_TOP + 31) / FIG_H, "→", ha="center", va="center", fontsize=16)
fig.text(10 / 500, 1 - (F_TOP - 16) / FIG_H, "f", fontsize=HS.LETTER, fontweight="bold", ha="left", va="top")
fig.canvas.draw()
_p = paxes[-1].get_position(); _c = CBFIX.get_position()
CBFIX.set_position([_c.x0, _p.y0, _c.width, _p.height])
HS.check(fig)
fig.savefig(OUT + "_Figure3_without_drawings.pdf", dpi=600)
# key drawings: Xindan's vector art (Cell cartoons/<type>_drawing.pdf), fitted in each slot, bottom-left aligned
import pymupdf                                     # noqa: E402
W_, H_ = fig.get_size_inches() * 72
doc = pymupdf.open(OUT + "_Figure3_without_drawings.pdf"); pg = doc[0]
DRAW_H = {"CHIN": 24, "CHID": 26, "OPCID": 16}     # pt; heights relative to each other as in her key
for typ, x, y, w, h in KEY_SLOTS:
    art = pymupdf.open(ROOT + f"Cell cartoons/{typ}_drawing.pdf"); aw, ah = art[0].rect.width, art[0].rect.height
    s_ = min(DRAW_H[typ] / ah, w * W_ / aw)
    x0, yb = x * W_, (1 - y) * H_ - 1                # bottom-left of the slot (PDF y runs down)
    pg.show_pdf_page(pymupdf.Rect(x0, yb - ah * s_, x0 + aw * s_, yb), art, 0)
for _, fn, x in CELLS:                              # panel f cells
    art = pymupdf.open(ROOT + "Cell cartoons/" + fn); aw, ah = art[0].rect.width, art[0].rect.height
    pg.show_pdf_page(pymupdf.Rect(x, F_TOP, x + 102, F_TOP + 102 * ah / aw), art, 0)
ICONS3 = pymupdf.open(ROOT + "Cell cartoons/Fig3_legend_icons.pdf")       # legend icons from her figure (H-NS, AT-rich DNA)
ICON3 = pymupdf.Rect(411.0, 604.0, 425.5, 629.5)
pg.show_pdf_page(pymupdf.Rect(359.5, F_TOP + 10, 359.5 + ICON3.width, F_TOP + 10 + ICON3.height), ICONS3, 0)
for k, lab in enumerate(("H-NS", "AT-rich DNA")):
    pg.insert_text((374, F_TOP + 17 + k * 15), lab, fontsize=HS.TITLE, fontname="ArialF3",
                   fontfile="/System/Library/Fonts/Supplemental/Arial.ttf")
doc.save(OUT + "Figure3_hairpins.pdf", garbage=3, deflate=True)
pg.get_pixmap(dpi=300).save(OUT + "Figure3_hairpins.png")
print("saved")
