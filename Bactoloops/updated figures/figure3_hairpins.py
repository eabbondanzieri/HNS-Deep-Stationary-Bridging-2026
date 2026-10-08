"""Figure 3: H-NS hairpins persist into deep stationary phase, and the hairpin domains at loop anchors compact as loops form.
(Panels i-m of the previous Supplementary Fig. 9, moved to the main text on 2026-10-06.)

  a  genome tracks (z-scores, 10-kb bins): WT 96 h ILD, anti-H-NS ChIP-seq (96 h), and the density (bp per bin) of the
     exponential-phase CHIDs and CHINs of Gavrilov, Shamovsky et al. 2025 (Nudler CHIN paper/opcidsRenamed.xlsx, lab
     genome); correlations in ../chin_chid_vs_loops.py
  b  aggregate contact maps (log2 O/E, balanced 500-bp maps from full valid pairs without the 5-kb cutoff) around CHINs >= 2 kb,
     centred on the CHIN centre and rescaled by CHIN length (ends at -1, +1; Nudler CHIN paper/scaled_pileups.py), averaged over
     the four orientations, 2-pixel Gaussian blur, each map relative to its own mean; dashed, anti-diagonal (hairpin stem);
     WT exp, WT 96 h and dstpA 96 h pool two libraries
  c  CHIN hairpin test on per-CHIN rescaled maps (Nudler CHIN paper/scaled_stem_test.py); median over CHINs, bootstrap 95% CI,
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


fig = plt.figure(figsize=HS.size(500, 520))
outer = fig.add_gridspec(3, 1, height_ratios=[0.6, 0.5, 1.0], hspace=0.42, left=0.1, right=0.985, top=0.965, bottom=0.115)
_TR, _PG, _ROW = 0, 1, 2
# a: exponential-phase CHIDs / CHINs vs deep-stationary loops and H-NS
X = pd.read_excel(ROOT + "Nudler CHIN paper/opcidsRenamed.xlsx", 0)
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
       ("CHID density, exp", to_tiff(_z(_dens(X[X.type == "CHID"]))), "#6a2c91"),
       ("CHIN density, exp", to_tiff(_z(_dens(X[X.type == "CHIN"]))), "#2e7d32")]
ax = fig.add_subplot(outer[_TR])
STEP = 7
for k, (lab, v, col) in enumerate(TRK):
    ax.plot(np.arange(N) + 0.5, np.clip(v, -2, 9) - k * STEP, color=col, lw=0.6)
    ax.text(-6, -k * STEP + 1, lab, color=col, fontsize=HS.BASE, ha="right", va="center")
ax.set_xlim(0, N); ax.set_ylim(-3 * STEP - 3, 9.5); ax.set_yticks([])
ax.set_xticks(TICK_POS); ax.set_xticklabels([3670, 0, 1000, 2000, 3000]); ax.set_xlabel("Genome position (kb)", fontsize=HS.BASE)
ax.plot([N + 4, N + 4], [-3 * STEP, -3 * STEP + 5], color="#1f1f1e", lw=0.9, clip_on=False)
ax.text(N + 7, -3 * STEP + 2.5, "5 z", fontsize=HS.BASE, va="center", ha="left")
fig.canvas.draw()
p_ = ax.get_position(); ax.set_position([0.25, p_.y0, 0.69, p_.height])                          # room for the track labels
letter(ax, "a", dx=-0.3)

# b-e: short-range H-NS structure from full valid pairs (no 5-kb cutoff); replicates pooled (Nudler CHIN paper/s9_fullpairs_stats.py)
NUD = ROOT + "Nudler CHIN paper/"
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
ax.set_title("hairpins (CHINs)", fontsize=HS.BASE, pad=3)
letter(ax, "c", dx=-0.42)
ax = fig.add_subplot(row[0, 1])
dots(ax, "chid", "within / neighboring\nwindows", (0.75, 1.6), [0.8, 1, 1.2, 1.4])
ax.set_title("hairpin domains\n(CHIDs)", fontsize=HS.BASE, pad=3)
letter(ax, "d", dx=-0.42)
ax = fig.add_subplot(row[0, 2])
dots(ax, "opcid", "within / neighboring\nwindows", (0.85, 1.3), [0.9, 1, 1.1, 1.2])
ax.set_title("operon domains\n(OPCIDs)", fontsize=HS.BASE, pad=3)
letter(ax, "e", dx=-0.42)
fig.canvas.draw()
_p = paxes[-1].get_position(); _c = CBFIX.get_position()
CBFIX.set_position([_c.x0, _p.y0, _c.width, _p.height])
HS.check(fig)
fig.savefig(OUT + "Figure3_hairpins.pdf", dpi=600)
fig.savefig(OUT + "Figure3_hairpins.png", dpi=300)
print("saved")
