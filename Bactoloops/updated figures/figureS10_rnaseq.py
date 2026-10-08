"""Supplementary Figure 10 (reworked): RNA-seq, WT vs delta-hns, CFU (spike-in per viable cell) normalization.
Old panel c (Venn diagram) is dropped; old d and e are merged into one comparison.

  a, b  log2 fold change (WT / delta-hns, apeglm-shrunk) against mean normalized expression (DESeq2 baseMean with
        CFU size factors, i.e. reads per CFU up to a constant), exponential phase (a) and 96 h (b). Coloured: padj < 0.05
        and |log2 FC| >= 1 (as the Fig. 4 volcano plots, so the counts match); hns circled.
  c     mean log2 FC per 10-kb bin (genes whose start falls in the bin, hns left out, >= 3 genes; the RNA-seq statistic
        of Fig. 4d) at 96 h against exponential phase, coloured by H-NS ChIP z (96 h). Pearson r, exact circular-shift P.
  d     Exponential phase in 10-kb bins, as z-scores (compare Fig. 4d at 96 h): Hi-C integrated loop density (WT exp),
        H-NS-mCherry ChIP-seq (anti-mCherry, exp) and H-NS silencing = -log2 FC (exp; high = higher in delta-hns). Orange stars, as Fig. 4c:
        flagellar/chemotaxis gene clusters (flg; flh/che + fli), strongly higher in WT. Pearson r between
        tracks with exact circular-shift P (every shift >= 200 kb).
  e-g   Robustness of Fig. 4 to the DESeq2 normalization: figure4.py rerun with CFU (spike-in per viable cell; Fig. 4),
        OD (spike-in per OD unit) and Count (DESeq2 median-of-ratios on read counts) size factors (fig4_summary_<norm>.csv).
        e: per-bin Pearson r of RNA-seq (-log2 FC, 96 h) with H-NS ChIP and with loop density (Fig. 4e, f); exact
        circular-shift P. f: mean RNA-seq z at the anchors of the 124
        WT 96 h loops (Fig. 4g) with peak-bootstrap 95% CI; grey dash, 95th percentile of the rotated loop set; exact
        rotation P. g: genes with padj < 0.05 and |log2 FC| >= 1 (Fig. 4a, b), higher or lower in WT. h: H-NS silencing at
        H-NS-bound regions relative to the genome (Fig. 4h): mean -log2 FC of the genes in bins with 96 h H-NS ChIP z >= 2
        minus the genome-wide mean, exp and 96 h, with bootstrap 95% CI (over bound bins); stars: exact rotation P of the
        change.
Output: FigureS10.pdf / .png
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = HNS_ROOT
B = ROOT + "Bactoloops/"
OUT = B + "updated figures/"
sys.path.insert(0, B)
import rotation_test as RT                     # noqa: E402

PADJ, LFC = 0.05, 1.0
COL_UP, COL_DOWN, COL_NS = "#e8524a", "#4a55d6", "#9a9a9a"
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#e4e2dd"
N, BIN, GENOME_LEN, ORIC = 461, 10_000, 4_607_088, 3_671_577
import sys as _hs_sys
_hs_sys.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()


def letter(ax, s, dx=-0.13, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=HS.LETTER, fontweight="bold", va="bottom", ha="right")


def load(tp):
    r = pd.read_csv(ROOT + f"Spike-In/genome_positions/DESeq2_apeglm_CFU_{tp}_with_positions.txt", sep="\t")
    r = r.dropna(subset=["padj", "log2FoldChange"])
    r["cls"] = np.where((r.padj < PADJ) & (r.log2FoldChange >= LFC), "up",
                        np.where((r.padj < PADJ) & (r.log2FoldChange <= -LFC), "down", "ns"))
    return r


D = {tp: load(tp) for tp in ["Exp", "96h"]}
from matplotlib.backends.backend_pdf import PdfPages
fig = fig1 = plt.figure(figsize=HS.size(500, 450))             # page 1: a-d
outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.05], hspace=0.36, left=0.09, right=0.94, top=0.96, bottom=0.08)
top = outer[0].subgridspec(1, 3, wspace=0.42, width_ratios=[1, 1, 1.15])

# ---- a, b ----
xlo = min(np.log10(d.baseMean[d.baseMean > 0]).min() for d in D.values()) - 0.2
xhi = max(np.log10(d.baseMean.max()) for d in D.values()) + 0.2
ylim = max(d.log2FoldChange.abs().max() for d in D.values()) * 1.08
for c, (L, tp, title) in enumerate([("a", "Exp", "exp"), ("b", "96h", "96 h")]):
    ax = fig.add_subplot(top[0, c])
    d = D[tp]
    x = np.log10(d.baseMean.clip(lower=1e-3))
    for cls, col in [("ns", COL_NS), ("down", COL_DOWN), ("up", COL_UP)]:
        m = d.cls == cls
        ax.scatter(x[m], d.log2FoldChange[m], s=3, color=col, alpha=0.55, lw=0, rasterized=True)
    h = d[d.GeneID == "hns"].iloc[0]
    ax.scatter([np.log10(h.baseMean)], [h.log2FoldChange], s=28, facecolor="none", edgecolor=INK, lw=0.6)
    ax.text(np.log10(h.baseMean) - 0.35, h.log2FoldChange, r"$\it{hns}$", ha="right", va="center", fontsize=HS.BASE)
    ax.text(0.04, 0.05, f"n={(d.cls == 'down').sum():,}", transform=ax.transAxes, color=COL_DOWN, fontsize=HS.BASE)
    ax.text(0.04, 0.93, f"n={(d.cls == 'up').sum():,}", transform=ax.transAxes, color=COL_UP, fontsize=HS.BASE, va="top")
    ax.axhline(0, color=MUTED, lw=0.5, ls=":")
    ax.set_xlim(xlo, xhi); ax.set_ylim(-ylim, ylim)
    ax.set_xlabel("log$_{10}$ mean expression (per CFU)"); ax.set_ylabel("log$_2$ FC (WT / Δ$\\it{hns}$)")
    ax.set_title(title, fontsize=HS.TITLE, pad=3)
    letter(ax, L)

# ---- c: mean log2 FC per 10-kb bin, exp vs 96 h (the RNA-seq statistic of Fig. 4d) ----
chip = pd.read_csv(B + "ChIP-seq_data/zscore_wt_96h_10kb.csv")
zbin = np.full(N, np.nan); zbin[(chip.start_bp // BIN).astype(int) % N] = chip.zscore.values


def bin_mean(d):
    d = d[(d.GeneID != "hns") & d.Gene_start_lab.notna()]
    a = d.groupby(((d.Gene_start_lab - 1) // BIN).astype(int)).log2FoldChange.agg(["mean", "size"])
    return a["mean"].where(a["size"] >= 3).reindex(range(N)).values


be, b9 = bin_mean(D["Exp"]), bin_mean(D["96h"])
ok = np.isfinite(be) & np.isfinite(b9)
ax = fig.add_subplot(top[0, 2])
sc = ax.scatter(be[ok], b9[ok], c=np.clip(zbin[ok], -1, 4), cmap="Reds", vmin=-1, vmax=4, s=7, lw=0.2, edgecolor="#8a8a8a")
r_bin = np.corrcoef(be[ok], b9[ok])[0, 1]
null = []
for k in RT.shifts(N, exclude=20):
    bb = np.roll(b9, k); o = np.isfinite(be) & np.isfinite(bb); null.append(np.corrcoef(be[o], bb[o])[0, 1])
p_bin = (1 + np.sum(np.abs(null) >= abs(r_bin))) / (len(null) + 1)
lo_, hi_ = np.nanmin([be[ok].min(), b9[ok].min()]) - 0.3, np.nanmax([be[ok].max(), b9[ok].max()]) + 0.3
ax.plot([lo_, hi_], [lo_, hi_], color=MUTED, lw=0.5, ls=(0, (3, 2)))
ax.axhline(0, color=MUTED, lw=0.5, ls=":"); ax.axvline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xlim(lo_, hi_); ax.set_ylim(lo_, hi_); ax.set_box_aspect(1)
ax.set_xlabel("mean log$_2$ FC per 10-kb bin, exp"); ax.set_ylabel("mean log$_2$ FC per 10-kb bin, 96 h")
ax.text(0.04, 0.96, f"$r$ = {r_bin:.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=HS.BASE)
cb = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.03); cb.set_label("H-NS ChIP z (96 h)", fontsize=HS.BASE); cb.ax.tick_params(labelsize=HS.BASE)
ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "c", dx=-0.3)
bound = zbin >= 2
print(f"c: {ok.sum()} bins; r {r_bin:.2f} (circular-shift P {p_bin:.4f}), rho {spearmanr(be[ok], b9[ok])[0]:.2f}; "
      f"H-NS-bound bins (z>=2): mean log2FC exp {np.nanmean(be[bound]):.2f}, 96 h {np.nanmean(b9[bound]):.2f}; "
      f"others exp {np.nanmean(be[~bound]):.2f}, 96 h {np.nanmean(b9[~bound]):.2f}")

# ---- d: exponential-phase tracks ----
zscore = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
rna = D["Exp"]; rna = rna[(rna.GeneID != "hns") & rna.Gene_start_lab.notna()]
agg = rna.groupby(((rna.Gene_start_lab - 1) // BIN).astype(int)).log2FoldChange.agg(["mean", "size"])
rna_bin = agg["mean"].where(agg["size"] >= 3).reindex(range(N)).values
mch = pd.read_csv(B + "ChIP-seq_data/zscore_mCherry_Exp_10kb.csv")
mch_bin = np.full(N, np.nan); mch_bin[(mch.start_bp // BIN).astype(int) % N] = mch.zscore.values
ild_t = pd.read_csv(B + "FIJI_output/integrated_loop_density_wt_Exp.csv").ILD_au.values
ild_bin = np.empty(N); ild_bin[(367 + np.arange(N)) % N] = ild_t
T = [("Hi-C loop density", zscore(ild_bin), "#9a7b4f"), ("H-NS-mCherry ChIP-seq", zscore(mch_bin), "#c0392b"),
     ("H-NS silencing (RNA-seq)", -zscore(rna_bin), "#999900")]
centre = (np.arange(N) * BIN + BIN / 2) / 1000
x_kb = np.where(centre * 1000 > ORIC, centre - GENOME_LEN / 1000, centre); order = np.argsort(x_kb)
x0 = (ORIC - GENOME_LEN) / 1000
bot = outer[1].subgridspec(1, 3, width_ratios=[0.15, 0.67, 0.18], wspace=0.04)   # label column | tracks | correlations
ax = fig.add_subplot(bot[0, 1])
STEP = 6.5
for i, (lab, v, col) in enumerate(T):
    ax.plot(x_kb[order], np.clip(v[order], -3.5, 6) - i * STEP, color=col, lw=0.7)
    ax.text(x0 - 60, -i * STEP + 0.4, lab, color=col, fontsize=HS.BASE, ha="right", va="center")
# orange stars (as Fig. 4c): flagellar/chemotaxis clusters, strongly higher in WT -> downward spikes in -log2 FC
yr = np.clip(T[2][1], -3.5, 6) - 2 * STEP
for lo, hi in [(1131, 1142), (1925, 1986)]:
    m = (x_kb >= lo - 5) & (x_kb <= hi + 5)
    ax.plot([x_kb[m].mean()], [yr[m].min() - 1.0], marker="*", ms=7, color="#f08c1a", mec=INK, mew=0.4, clip_on=False, zorder=5)
ax.text(x0 - 60, -2 * STEP - 2.2, "exp, z-scores", color=INK, fontsize=HS.BASE, ha="right", va="center", style="italic")
ax.set_xlim(x0, ORIC / 1000); ax.set_ylim(-2 * STEP - 5.2, 7)
ax.set_xticks([x0, 0, 1000, 2000, 3000]); ax.set_xticklabels([3670, 0, 1000, 2000, 3000])
ax.set_xlabel("Genome position (kb)"); ax.set_yticks([])
xr = ORIC / 1000 + 60
ax.plot([xr, xr], [-2 * STEP, -2 * STEP + 4], color=INK, lw=0.9, clip_on=False)
ax.text(xr + 40, -2 * STEP + 2, "4 z", fontsize=HS.BASE, va="center", ha="left")
letter(ax, "d", dx=-0.23, dy=1.0)


def shift_r(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    r0 = np.corrcoef(a[ok], b[ok])[0, 1]
    null = []
    for k in RT.shifts(N, exclude=20):
        bb = np.roll(b, k); o = np.isfinite(a) & np.isfinite(bb)
        null.append(np.corrcoef(a[o], bb[o])[0, 1])
    return r0, (1 + np.sum(np.abs(null) >= abs(r0))) / (len(null) + 1)


axt = fig.add_subplot(bot[0, 2]); axt.axis("off")
lines = [("loop density vs ChIP", T[0][1], T[1][1]), ("loop density vs RNA-seq", T[0][1], T[2][1]),
         ("RNA-seq vs ChIP", T[2][1], T[1][1])]
axt.text(0, 0.95, "Pearson r (exp)", fontsize=HS.BASE, fontweight="bold", va="top", transform=axt.transAxes)
for j, (lab, a, b_) in enumerate(lines):
    r, p = shift_r(a, b_)
    axt.text(0, 0.8 - j * 0.16, f"{lab}\n$r$ = {r:.2f} (P = {p:.3f})", fontsize=HS.BASE, va="top", transform=axt.transAxes)
    print(f"d: {lab}: r {r:.2f}, circular-shift P {p:.4f}")

# ---- e-g: Fig. 4 under the three DESeq2 normalizations ----
NORMS = [("CFU", "CFU (Fig. 5)", INK), ("OD", "OD", "#7b6fb0"), ("Count", "Count", "#c9a227")]
SM = {k: pd.read_csv(OUT + f"fig4_summary_{k}.csv").iloc[0] for k, _, _ in NORMS}
fig = fig2 = plt.figure(figsize=HS.size(500, 470))             # page 2: e-h
gEG = fig.add_gridspec(2, 2, wspace=0.45, hspace=0.42, left=0.12, right=0.98, top=0.88, bottom=0.1)
dx = {"CFU": -0.22, "OD": 0, "Count": 0.22}

ax = fig.add_subplot(gEG[0, 0])
MEAS = [("r_e", "p_e", "ChIP"), ("r_f", "p_f", "loops")]
for i, (key, pk, _) in enumerate(MEAS):
    for k, lab, col in NORMS:
        v = SM[k][key]
        ax.plot(i + dx[k], v, "o", ms=4.5, color=col, mec=INK, mew=0.4, label=lab if i == 0 else None, zorder=3)
    if pk is not None:
        ax.text(i, max(SM[k][key] for k, _, _ in NORMS) + 0.06, RT.stars(max(SM[k][pk] for k, _, _ in NORMS)),
                ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xticks(range(len(MEAS))); ax.set_xticklabels([m[2] for m in MEAS], fontsize=HS.BASE)
ax.set_xlim(-0.55, len(MEAS) - 0.45); ax.set_ylim(0, 1.0)
ax.set_ylabel("Pearson r with H-NS\nsilencing (10-kb bins, 96 h)")
fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.55, 1.0),
           title="DESeq2 normalization")
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "e", dx=-0.2)

ax = fig.add_subplot(gEG[0, 1])
for j, (k, lab, col) in enumerate(NORMS):
    m, lo, hi = SM[k].g_rna_mean, SM[k].g_rna_lo, SM[k].g_rna_hi
    ax.errorbar(j, m, yerr=[[m - lo], [hi - m]], fmt="o", ms=4.5, color=col, mec=INK, mew=0.4, lw=0.8, capsize=1.5, zorder=3)
    ax.hlines(SM[k].g_rna_rot95, j - 0.3, j + 0.3, color="#9a9a9a", lw=0.9, ls=(0, (2, 1)))
    ax.text(j, hi + 0.08, RT.stars(SM[k].g_rna_p), ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xticks(range(3)); ax.set_xticklabels([n[0] for n in NORMS], fontsize=HS.BASE)
ax.set_xlim(-0.6, 2.6); ax.set_ylim(-0.3, 2.7)
ax.set_ylabel("mean H-NS silencing z\nat loop anchors (96 h)")
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "f", dx=-0.2)

ax = fig.add_subplot(gEG[1, 0])
GRP = [("up_Exp", "exp\nup"), ("down_Exp", "exp\ndown"), ("up_96h", "96 h\nup"), ("down_96h", "96 h\ndown")]
for i, (key, _) in enumerate(GRP):
    for k, lab, col in NORMS:
        v = SM[k][key]
        ax.bar(i + dx[k], v, width=0.2, color=col, edgecolor=INK, lw=0.3)
        ax.text(i + dx[k], v + 40, f"{v:,}", ha="center", va="bottom", fontsize=HS.BASE, rotation=90)
ax.set_xticks(range(len(GRP))); ax.set_xticklabels([g[1] for g in GRP], fontsize=HS.BASE)
ax.set_xlim(-0.55, len(GRP) - 0.45); ax.set_ylim(0, 3700)
ax.set_ylabel("DE genes")
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "g", dx=-0.2)

ax = fig.add_subplot(gEG[1, 1])
for k, lab, col in NORMS:
    xs = np.array([0, 1]) + dx[k]
    ys = [SM[k].h_Exp, SM[k]["h_96h"]]
    ax.plot(xs, ys, color=col, lw=0.6, zorder=2)
    for j, tp in enumerate(["Exp", "96h"]):
        ax.errorbar(xs[j], ys[j], yerr=[[ys[j] - SM[k][f"h_{tp}_lo"]], [SM[k][f"h_{tp}_hi"] - ys[j]]], fmt="o", ms=4,
                    color=col, mec=INK, mew=0.4, lw=0.8, capsize=1.2, zorder=3)
    ax.text(1 + dx[k], SM[k]["h_96h_hi"] + 0.12, RT.stars(SM[k].h_diff_p), ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color=MUTED, lw=0.5, ls=":")
ax.set_xticks([0, 1]); ax.set_xticklabels(["exp", "96 h"])
ax.set_xlim(-0.5, 1.5); ax.set_ylim(-0.4, 4.1)
ax.set_ylabel("silencing at bound genes\n(log$_2$, rel. to genome)")
ax.grid(axis="y", color=GRID, lw=0.4); ax.set_axisbelow(True)
letter(ax, "h", dx=-0.2)

for f_ in (fig1, fig2):
    HS.check(f_)
with PdfPages(OUT + "FigureS10.pdf") as pp:                  # two pages: a-d, e-h
    pp.savefig(fig1, dpi=600); pp.savefig(fig2, dpi=600)
fig1.savefig(OUT + "FigureS10.png", dpi=300)
fig2.savefig(OUT + "FigureS10_p2.png", dpi=300)
print("saved")
