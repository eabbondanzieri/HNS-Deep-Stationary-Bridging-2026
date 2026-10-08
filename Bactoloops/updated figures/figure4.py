"""Figure 4 (in progress), from the old Figure 4 layout. Data come from the "Volcano plots from DESeq2
output" project (../Spike-In/).

a, b  Volcano plots of RNA-seq, WT vs delta-hns, exponential phase (a) and 96 h (b): DESeq2 with apeglm
      shrinkage and CFU (spike-in per cell) normalization (Spike-In/DESeq2_apeglm_CFU_<Exp|96h>.txt).
      log2 FC > 0 = higher in WT. Genes with padj < 0.05 and |log2 FC| >= 1 are coloured (same cutoffs
      as Spike-In/Volcano_and_Counts.R). Both panels share axes; the y axis is capped (CFU: -log10 padj = 25) and the
      few genes above it (CFU exp: gadA/B/C, hdeA/B, hns) are drawn as triangles on the top edge and named.
c     log2 FC along the genome, one stick per gene (all genes tested, hns left out; "sigonly" restricts to
      padj < 0.05), at the gene's start, exp track above the 96 h track; left, each track's log2 FC histogram
      (fraction of that track's genes per 1-log2 bin). Same data as Spike-In/Genome_Position_Plot.R
      (genome_position_output/LFC_vs_genome_CFU_genes.txt). Genome axis runs oriC -> oriC on the lab
      genome (W3110 minus 1,398,953-1,438,196; oriC 3,671,577 bp), labelled as in Figure 2. Orange
      stars: clusters of flagellar/chemotaxis genes higher in WT in exponential phase.
d     Genome tracks at 96 h, each as a z-score over the 461 10-kb bins: Hi-C integrated loop density
      (WT 96 h, the macro's current output, as in Figure 2), H-NS ChIP-seq (anti-H-NS, cWX2200) and
      H-NS silencing = -(mean log2 FC) (CFU; log2(delta-hns/WT)) of all genes whose start falls in the bin (hns left
      out; every bin has >= 3 genes); positive = higher in delta-hns. The RNA-seq statistic is the one in Spike-In/Bin_Correlation.R
      (reproduced exactly; that script used an older ILD file, so its correlations differ slightly).
e-g   RNA-seq vs H-NS ChIP and loop density per bin (e, f); RNA-seq and ChIP at the anchors of the WT 96 h loops (g).
h     H-NS silencing at H-NS-bound regions relative to the genome: mean -log2 FC (= log2 delta-hns/WT) of the genes in
      the bins with 96 h H-NS ChIP z >= 2 minus the genome-wide mean, exp and 96 h (removes the global per-cell shift at
      96 h). 95% CI: bootstrap over bound bins; grey dash: 95th percentile of the rotated bound-bin set; bracket: exact
      rotation P of the change. "Silencing" (positive = higher without H-NS) is used for d-h; a-c keep log2 FC = WT/delta-hns.
Normalization: CFU by default; "python figure4.py OD" or "... Count" redraws every panel with that DESeq2 normalization
(for the robustness supplement), written to Figure4_<norm>.pdf / .png. Panel c plots every gene; adding "sigonly" (e.g.
"python figure4.py sigonly") restricts it to padj < 0.05, written to Figure4[_<norm>]_sigonly. Panel c is built from
Spike-In/genome_positions/DESeq2_apeglm_<norm>_<tp>_with_positions.txt (padj < 0.05, hns left out; identical to
LFC_vs_genome_CFU_genes.txt for CFU).
Output: Figure4.pdf / .png (CFU) or Figure4_<norm>.pdf / .png (in this folder).
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
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = HNS_ROOT
OUT = ROOT + "Bactoloops/updated figures/"
ALLGENES = "sigonly" not in _sys.argv[1:]                      # panel c: every gene (default); "sigonly" = padj < 0.05
_args = [a for a in _sys.argv[1:] if a != "sigonly"]
NORM = _args[0] if _args else "CFU"                             # CFU (default), OD or Count
assert NORM in ("CFU", "OD", "Count"), NORM
DE = ROOT + "Spike-In/DESeq2_apeglm_" + NORM + "_{}.txt"
POSF = ROOT + "Spike-In/genome_positions/DESeq2_apeglm_" + NORM + "_{}_with_positions.txt"
OUTNAME = ("Figure4" if NORM == "CFU" else f"Figure4_{NORM}") + ("" if ALLGENES else "_sigonly")
PADJ, LFC = 0.05, 1.0
COL_UP, COL_DOWN, COL_NS = "#e8524a", "#4a55d6", "#9a9a9a"
COL_TRACK = {"Exp": "#999900", "96h": "#009999"}          # Spike-In/Genome_Position_Plot.R colours, swapped: teal (96 h) reads better on white and next to red
OFFSET = {"Exp": 20, "96h": 0}
GENOME_LEN, ORIC = 4_607_088, 3_671_577
INK, GUIDE = "#1f1f1e", "#6666ff"
STATS = dict(norm=NORM)                                    # key numbers, saved to fig4_summary_<norm>.csv for Fig. S10
import sys as _hs_sys
_hs_sys.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()


def letter(ax, s, dx=-0.22, dy=1.02):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=HS.LETTER, fontweight="bold", va="bottom", ha="right")


def load(tp):
    r = pd.read_csv(DE.format(tp), sep="\t").dropna(subset=["padj", "log2FoldChange"])
    r["y"] = -np.log10(r.padj.clip(lower=1e-300))
    r["cls"] = np.where((r.padj < PADJ) & (r.log2FoldChange >= LFC), "up",
                        np.where((r.padj < PADJ) & (r.log2FoldChange <= -LFC), "down", "ns"))
    return r


fig = plt.figure(figsize=HS.size(500, 600))
outer = fig.add_gridspec(4, 1, height_ratios=[0.85, 1.15, 0.8, 1.05], hspace=0.42, left=0.08, right=0.985, top=0.965, bottom=0.055)
TRACK_W = [0.16, 0.8, 0.04]                                    # c and d share these columns, so their genome axes line up

# ---------------- a, b: volcano plots ----------------
data = {tp: load(tp) for tp in ["Exp", "96h"]}
xlim = max(d.log2FoldChange.abs().max() for d in data.values()) * 1.08
# -log10 padj ceiling: a few genes (CFU exp: gadA/B/C, hdeA/B, hns) would otherwise compress the rest. The cap is the
# smallest multiple of 5 that leaves <= 6 genes above it in either panel; those are drawn as triangles on the top edge.
YCAP = 5 * np.ceil(max(np.sort(d.y.values)[-7] for d in data.values()) / 5)
ymax = YCAP * 1.06
gAB = outer[0].subgridspec(1, 2, wspace=0.24)
for c, (L, tp, title) in enumerate([("a", "Exp", "exp"), ("b", "96h", "96 h")]):
    ax = fig.add_subplot(gAB[0, c])
    d = data[tp]
    over = d.y > YCAP
    for cls, col, a in [("ns", COL_NS, 0.5), ("down", COL_DOWN, 0.55), ("up", COL_UP, 0.55)]:
        s = d[(d.cls == cls) & ~over]
        ax.scatter(s.log2FoldChange, s.y, s=3, color=col, alpha=a, lw=0, rasterized=True)
        s = d[(d.cls == cls) & over]
        ax.scatter(s.log2FoldChange, np.full(len(s), YCAP), s=14, marker="^", color=col, edgecolor=INK, lw=0.3,
                   clip_on=False, zorder=4)
    n_up, n_down = (d.cls == "up").sum(), (d.cls == "down").sum()
    ax.text(0.05, 0.62, f"n={n_down:,}", transform=ax.transAxes, color=COL_DOWN, fontsize=HS.BASE, ha="left")
    ax.text(0.95, 0.62, f"n={n_up:,}", transform=ax.transAxes, color=COL_UP, fontsize=HS.BASE, ha="right")
    h = d[d.GeneID == "hns"].iloc[0]
    hy = min(h.y, YCAP)
    ax.scatter([h.log2FoldChange], [hy], s=28, facecolor="none", edgecolor=INK, lw=0.6, clip_on=False, zorder=5)
    ax.text(h.log2FoldChange - 1.3, hy - (1.6 if h.y > YCAP else 0), r"$\it{hns}$", ha="right", va="center", fontsize=HS.BASE)
    cap_genes = d[over & (d.GeneID != "hns")]
    if len(cap_genes):                                         # name the off-scale genes beside their triangles
        _nm = [r"$\it{" + g.replace("_", r"\_") + "}$" for g in sorted(cap_genes.GeneID)]
        names = (", ".join(_nm[:2]) + ",\n" + ", ".join(_nm[2:]) if len(_nm) > 3 else ", ".join(_nm)) if len(_nm) <= 6 \
            else f"{len(cap_genes)} genes above"
        if cap_genes.log2FoldChange.mean() < 0:                # label on the side of the triangles facing the centre
            ax.text(cap_genes.log2FoldChange.max() + 0.8, YCAP * 0.985, names, fontsize=HS.BASE, ha="left", va="top")
        else:
            ax.text(cap_genes.log2FoldChange.min() - 0.8, YCAP * 0.985, names, fontsize=HS.BASE, ha="right", va="top")
    ax.axhline(YCAP, color="#b5b3ad", lw=0.4, ls=(0, (2, 2)), zorder=0)
    step = next(v for v in (5, 10, 20, 25, 50, 100) if YCAP / v <= 6)
    yt = np.arange(0, YCAP + 1, step)
    yt = yt[yt < YCAP * 0.9]
    ax.set_yticks(list(yt) + [YCAP]); ax.set_yticklabels([f"{v:.0f}" for v in yt] + [f"≥{YCAP:.0f}"])
    ax.set_xlim(-xlim, xlim); ax.set_ylim(-ymax * 0.03, ymax)
    if len(over):
        print(f"{tp}: capped at -log10 padj {YCAP}: " + ", ".join(f"{g} ({y:.0f})" for g, y in zip(d[over].GeneID, d[over].y)))
    ax.set_xlabel("log$_2$ (FC)")
    ax.set_ylabel("−log$_{10}$ ($p_{adj}$)")
    ax.set_title(title, fontsize=HS.TITLE, pad=3)
    letter(ax, L, dx=-0.1, dy=1.0)
    STATS[f"up_{tp}"], STATS[f"down_{tp}"] = n_up, n_down
    print(f"{tp}: {len(d)} genes tested, up {n_up}, down {n_down}; hns log2FC {h.log2FoldChange:.1f}, padj {h.padj:.1e}")

# ---------------- c: log2 FC along the genome, with histograms ----------------
_g = []
for _tp in ["96h", "Exp"]:
    _r = pd.read_csv(POSF.format(_tp), sep="\t")
    _r = _r[((_r.padj < PADJ) | ALLGENES) & (_r.GeneID != "hns") & _r.Gene_start_lab.notna() & _r.log2FoldChange.notna()]
    _g.append(pd.DataFrame({"Track": _tp, "GeneID": _r.GeneID, "log2FoldChange": _r.log2FoldChange,
                            "Position_from_oriC_axis": np.where(_r.Gene_start_lab > ORIC, _r.Gene_start_lab - GENOME_LEN, _r.Gene_start_lab)}))
G = pd.concat(_g, ignore_index=True)
YLIM = (-11, 29)
gC = outer[1].subgridspec(1, 3, width_ratios=TRACK_W, wspace=0.03)
axh = fig.add_subplot(gC[0, 0])
ax = fig.add_subplot(gC[0, 1], sharey=axh)
x0_kb = (ORIC - GENOME_LEN) / 1000
for tr in ["96h", "Exp"]:
    g = G[G.Track == tr]
    x = g.Position_from_oriC_axis.values / 1000
    y = g.log2FoldChange.values
    ax.vlines(x, OFFSET[tr], np.clip(OFFSET[tr] + y, *YLIM), color=COL_TRACK[tr], lw=0.45, alpha=0.45)
    edges = np.arange(np.floor(y.min()), np.ceil(y.max()) + 1)
    frac = np.histogram(y, edges)[0] / len(y)
    axh.barh(edges[:-1] + OFFSET[tr] + 0.5, frac, height=1.0, color=COL_TRACK[tr], edgecolor=INK, lw=0.3)
    axh.text(0.04, OFFSET[tr] + (7.5 if tr == "Exp" else 9), "exp" if tr == "Exp" else "96 h",
             transform=axh.get_yaxis_transform(), color=COL_TRACK[tr], fontsize=HS.TITLE, fontweight="bold",
             ha="left", va="center")                            # top-left of each track, clear of the bars
    print(f"c {tr}: {len(g)} genes{'' if ALLGENES else f' with padj < {PADJ}'} ({(y > 0).sum()} higher in WT, {(y < 0).sum()} lower)")
# stars: clusters of genes higher in WT in exponential phase (flagellar and chemotaxis genes; see
# the Spike-In session's cluster table). The flh/che (1,925-1,940 kb) and fli (1,964-1,986 kb) clusters
# are ~40 kb apart, too close to separate at this width, so they share one star.
STAR_WINDOWS = [("flg", 1131, 1142), ("flh/che + fli", 1925, 1986)]
ge = G[(G.Track == "Exp") & (G.log2FoldChange > 0)]
for name, lo, hi in STAR_WINDOWS:
    w = ge[(ge.Position_from_oriC_axis / 1000 >= lo) & (ge.Position_from_oriC_axis / 1000 <= hi)]
    xs = w.Position_from_oriC_axis.mean() / 1000
    ax.plot([xs], [OFFSET["Exp"] + w.log2FoldChange.max() + 1.6], marker="*", ms=7, color="#f08c1a",
            mec=INK, mew=0.4, clip_on=False, zorder=5)
    STATS[f"star_{name}_n"], STATS[f"star_{name}_lfc_median"] = len(w), w.log2FoldChange.median()
    print(f"star {name}: {len(w)} genes, {lo}-{hi} kb, log2FC {w.log2FoldChange.min():.1f} to {w.log2FoldChange.max():.1f}")
for a_ in (ax, axh):
    for yg in (10, 20):
        a_.axhline(yg, color=GUIDE, lw=0.4, ls=(0, (4, 3)))
ax.set_xlim(x0_kb, ORIC / 1000)
ax.set_xticks([x0_kb, 0, 1000, 2000, 3000]); ax.set_xticklabels([3670, 0, 1000, 2000, 3000])
ax.set_xlabel("Genome position (kb)")
ax.tick_params(labelleft=False)
axh.set_ylim(*YLIM)
yt = [-10, -5, 0, 5, 15, 20, 25]                               # per-track ticks: 96 h around 0, exp around 20
axh.set_yticks(yt); axh.set_yticklabels([f"{v - (20 if v > 10 else 0)}".replace("-", "−") for v in yt])
axh.set_ylabel("log$_2$ (FC)")
axh.set_xlim(axh.get_xlim()[1] * 1.08, 0)                     # bars grow leftward, toward the axis label
axh.set_xticks([])
axh.set_xlabel("Norm. gene\ncounts")
letter(axh, "c", dx=-0.35, dy=1.0)

# ---------------- d: 96 h tracks as z-scores ----------------
N, BIN = 461, 10_000
zscore = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
rna = pd.read_csv(POSF.format("96h"), sep="\t")
rna = rna[(rna.GeneID != "hns") & rna.Gene_start_lab.notna() & rna.log2FoldChange.notna()]
agg = rna.groupby(((rna.Gene_start_lab - 1) // BIN).astype(int)).log2FoldChange.agg(["mean", "size"])
rna_bin = agg["mean"].where(agg["size"] >= 3).reindex(range(N)).values          # genome order
chip = pd.read_csv(ROOT + "Bactoloops/ChIP-seq_data/zscore_wt_96h_10kb.csv")
chip_bin = np.full(N, np.nan); chip_bin[(chip.start_bp // BIN).astype(int) % N] = chip.zscore.values
ild_tiff = pd.read_csv(ROOT + "Bactoloops/FIJI_output/integrated_loop_density_wt_96h_repavg.csv").ILD_au.values
ild_bin = np.empty(N); ild_bin[(367 + np.arange(N)) % N] = ild_tiff                 # TIFF bin 0 = genome bin 367 (oriC)
centre_kb = (np.arange(N) * BIN + BIN / 2) / 1000
x_kb = np.where(centre_kb * 1000 > ORIC, centre_kb - GENOME_LEN / 1000, centre_kb)  # same oriC-rotated axis as c
order = np.argsort(x_kb)
rna_z = -zscore(rna_bin)          # sign flipped: high = lower in WT than in delta-hns (repressed by H-NS)
TRACKS_D = [("Hi-C loop density", zscore(ild_bin), "#000000"),
            ("H-NS ChIP-seq", zscore(chip_bin), "#c0392b"),
            ("H-NS silencing (RNA-seq)", rna_z, COL_TRACK["96h"])]
gD = outer[2].subgridspec(1, 3, width_ratios=TRACK_W, wspace=0.03)
ax = fig.add_subplot(gD[0, 1])
STEP = 6.5
for i, (lab, v, col) in enumerate(TRACKS_D):
    y0 = -i * STEP
    ax.plot(x_kb[order], np.clip(v[order], -3.5, 6) + y0, color=col, lw=0.7)
    ax.text(x0_kb - 60, y0 + 0.4, lab, color=col, fontsize=HS.BASE, ha="right", va="center")
ax.text(x0_kb - 60, -2 * STEP - 2.2, "96 h, z-scores", color=INK, fontsize=HS.BASE, ha="right", va="center", style="italic")
ax.set_xlim(x0_kb, ORIC / 1000); ax.set_ylim(-2 * STEP - 4, 7)
ax.set_xticks([x0_kb, 0, 1000, 2000, 3000]); ax.set_xticklabels([3670, 0, 1000, 2000, 3000])
ax.set_xlabel("Genome position (kb)"); ax.set_yticks([])
xr = ORIC / 1000 + 60
ax.plot([xr, xr], [-2 * STEP, -2 * STEP + 4], color=INK, lw=0.9, clip_on=False)
ax.text(xr + 40, -2 * STEP + 2, "4 z", fontsize=HS.BASE, va="center", ha="left")
letter(ax, "d", dx=-0.2, dy=1.0)
# ---------------- e, f: RNA-seq against ChIP and loop density, per bin; g: loop anchors ----------------
chip_z, ild_z = zscore(chip_bin), zscore(ild_bin)
rng = np.random.default_rng(0)


def shift_p(a, b):
    """two-sided circular-shift p for Pearson r: b rotated around the chromosome, keeping each track's smoothness"""
    r0 = np.corrcoef(a, b)[0, 1]
    null = [np.corrcoef(a, np.roll(b, k))[0, 1] for k in RT.shifts(N, exclude=20)]   # every shift >= 200 kb, exact
    return (1 + np.sum(np.abs(null) >= abs(r0))) / (len(null) + 1)


def partial(y, x, z):
    res = lambda v: v - np.polyval(np.polyfit(z, v, 1), z)
    return np.corrcoef(res(y), res(x))[0, 1]


gEG = outer[3].subgridspec(1, 4, wspace=0.55, width_ratios=[1, 1, 1.05, 0.8])
for c, (L, x, xlab, other, oname, pcol) in enumerate([("e", chip_z, "ChIP-seq z-score", ild_z, "loop density", "#c0392b"),
                                                     ("f", ild_z, "loop density z-score", chip_z, "ChIP", "#000000")]):   # colours as the tracks in d
    ax = fig.add_subplot(gEG[0, c])
    r = np.corrcoef(x, rna_z)[0, 1]
    ax.scatter(x, rna_z, s=4, color=pcol, alpha=0.5, lw=0)
    fit = np.polyfit(x, rna_z, 1); xs = np.linspace(x.min(), x.max(), 20)
    ax.plot(xs, np.polyval(fit, xs), color=INK, lw=0.9)
    ax.set_box_aspect(1)
    ax.set_xlabel(xlab); ax.set_ylabel("H-NS silencing z-score")
    ax.text(0.04, 0.96, f"$r$ = {r:.2f}", transform=ax.transAxes, fontsize=HS.BASE, ha="left", va="top")
    ax.grid(color="#e4e2dd", lw=0.4); ax.set_axisbelow(True)
    letter(ax, L, dx=-0.32, dy=1.02)
    STATS[f"r_{L}"], STATS[f"p_{L}"], STATS[f"partial_{L}"] = r, shift_p(rna_z, x), partial(rna_z, x, other)
    print(f"{L}: RNA vs {xlab}: Pearson r {r:.2f}, Spearman rho {spearmanr(x, rna_z)[0]:.2f}, "
          f"circular-shift p {shift_p(rna_z, x):.4f}; controlling for {oname}: partial r {partial(rna_z, x, other):.2f}")

# g: mean z of each loop's two anchors (124 WT 96 h loops), as in the method comparison figure
Lp = pd.read_csv(ROOT + "Bactoloops/FIJI_output/loop_calls_wt_96h_repavg.csv")
P = np.c_[Lp.genome_bin1.round().astype(int), Lp.genome_bin2.round().astype(int)] % N
ax = fig.add_subplot(gEG[0, 2])
GROUPS = [("H-NS\nChIP-seq", chip_z, "#c0392b"), ("H-NS\nsilencing", rna_z, COL_TRACK["96h"])]
data = [np.mean(v[P], axis=1) for _, v, _ in GROUPS]
bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False, whis=(0, 100),     # whiskers: full range
                medianprops=dict(color=INK, lw=1), whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
for i, ((lab, v, col), b, d_) in enumerate(zip(GROUPS, bp["boxes"], data), start=1):
    b.set_facecolor(col); b.set_alpha(0.55); b.set_edgecolor(col)
    ax.scatter(i + np.random.default_rng(i).uniform(-0.16, 0.16, len(d_)), d_, s=2.5, color=INK, alpha=0.4, lw=0, zorder=3)
    null = [np.mean(v[(P + k) % N]) for k in RT.shifts(N)]              # loop set rotated around the chromosome
    ax.hlines(np.percentile(null, 95), i - 0.33, i + 0.33, color="#9a9a9a", lw=0.9, ls=(0, (2, 1)), zorder=4)
    boot = [np.mean(rng.choice(d_, len(d_))) for _ in range(2000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    ax.errorbar(i + 0.4, d_.mean(), yerr=[[d_.mean() - lo], [hi - d_.mean()]], fmt="D", ms=2.4, color=INK, lw=0.8,
                capsize=1.2, zorder=5)
    p_rot = RT.p_value(d_.mean(), null)
    k_ = "chip" if "ChIP" in lab else "rna"
    STATS.update({f"g_{k_}_mean": d_.mean(), f"g_{k_}_lo": lo, f"g_{k_}_hi": hi, f"g_{k_}_rot95": np.percentile(null, 95),
                  f"g_{k_}_p": p_rot})
    print(f"g: {lab.replace(chr(10), ' ')}: mean {d_.mean():.2f} (95% CI {lo:.2f}-{hi:.2f}), median {np.median(d_):.2f}, "
          f"rotation null mean {np.mean(null):.2f}, 95th {np.percentile(null, 95):.2f}, p {p_rot:.4f}")
ax.axhline(0, color="#6b6a66", lw=0.5, ls=":")
ax.set_xticks([1, 2]); ax.set_xticklabels([g[0] for g in GROUPS])
ax.set_xlim(0.45, 2.75)
ax.set_ylabel("Mean anchor z-score")
ax.set_title(f"WT 96 h loops (n = {len(Lp)})", fontsize=HS.BASE, pad=3)
ax.grid(axis="y", color="#e4e2dd", lw=0.4); ax.set_axisbelow(True)
letter(ax, "g", dx=-0.3, dy=1.0)

# h: repression at H-NS-bound regions relative to the whole genome, exp and 96 h. Genes (tested at both time points,
# hns left out) in the 27 bins with 96 h H-NS ChIP z >= 2 (same footprint in exp: r = 0.93 with the exp mCherry ChIP);
# offset = their mean log2 FC minus the genome-wide mean log2 FC, which removes the global per-cell shift at 96 h.
# CIs: bootstrap over bound bins (genes in a bin move together); P: exact rotation of the bound-bin set.
GE = {tp: pd.read_csv(POSF.format(tp), sep="\t") for tp in ("Exp", "96h")}
GE = {tp: g[(g.GeneID != "hns") & g.Gene_start_lab.notna() & g.log2FoldChange.notna()].set_index("GeneID") for tp, g in GE.items()}
common = GE["Exp"].index.intersection(GE["96h"].index)
gbin = ((GE["Exp"].loc[common].Gene_start_lab.values - 1) // BIN).astype(int)
FC = {tp: -GE[tp].loc[common].log2FoldChange.values for tp in GE}   # H-NS silencing: -log2 FC = log2(delta-hns/WT)
bound = np.where(chip_bin >= 2)[0]
offset = lambda f, bs: f[np.isin(gbin, bs)].mean() - f.mean()
members = {b: np.where(gbin == b)[0] for b in bound}
boot = []
for _ in range(2000):
    w = np.concatenate([members[b] for b in rng.choice(bound, len(bound))])
    boot.append([FC[tp][w].mean() - FC[tp].mean() for tp in ("Exp", "96h")])
boot = np.array(boot)
rot = np.array([[offset(FC[tp], (bound + k) % N) for tp in ("Exp", "96h")] for k in RT.shifts(N)])
ax = fig.add_subplot(gEG[0, 3])
obs = [offset(FC[tp], bound) for tp in ("Exp", "96h")]
for j, (tp, lab) in enumerate([("Exp", "exp"), ("96h", "96 h")]):
    lo, hi = np.percentile(boot[:, j], [2.5, 97.5])
    ax.errorbar(j, obs[j], yerr=[[obs[j] - lo], [hi - obs[j]]], fmt="o", ms=5, color=COL_TRACK[tp], mec=INK, mew=0.4,
                lw=0.9, capsize=1.5, zorder=4)
    ax.hlines(np.percentile(rot[:, j], 95), j - 0.25, j + 0.25, color="#9a9a9a", lw=0.9, ls=(0, (2, 1)))
    ax.text(j + 0.13, obs[j], f"{2 ** obs[j]:.1f}×", fontsize=HS.BASE, ha="left", va="center")
    p_ = RT.p_value(obs[j], rot[:, j])
    STATS.update({f"h_{tp}": obs[j], f"h_{tp}_lo": lo, f"h_{tp}_hi": hi, f"h_{tp}_p": p_})
d_obs = obs[1] - obs[0]; d_lo, d_hi = np.percentile(boot[:, 1] - boot[:, 0], [2.5, 97.5])
d_p = RT.p_value(d_obs, rot[:, 1] - rot[:, 0])
STATS.update(h_diff=d_obs, h_diff_lo=d_lo, h_diff_hi=d_hi, h_diff_p=d_p, h_bound_bins=len(bound), h_genes=int(np.isin(gbin, bound).sum()))
ax.plot([0, 1], obs, color="#9a9a9a", lw=0.6, zorder=2)
yb = 3.6
ax.plot([0, 0, 1, 1], [yb - 0.12, yb, yb, yb - 0.12], color=INK, lw=0.6)
ax.text(0.5, yb + 0.04, RT.stars(d_p), ha="center", va="bottom", fontsize=HS.BASE)
ax.axhline(0, color="#6b6a66", lw=0.5, ls=":")
ax.set_xticks([0, 1]); ax.set_xticklabels(["exp", "96 h"])
ax.set_xlim(-0.5, 1.9); ax.set_ylim(-0.9, 4.0)
ax.set_ylabel("silencing at bound genes\n(log$_2$, rel. to genome)")
ax.set_title(f"{STATS['h_genes']} genes, {len(bound)} bins", fontsize=HS.BASE, pad=3)
ax.grid(axis="y", color="#e4e2dd", lw=0.4); ax.set_axisbelow(True)
letter(ax, "h", dx=-0.42, dy=1.0)
print(f"h: offset exp {obs[0]:.2f} ({STATS['h_Exp_lo']:.2f} to {STATS['h_Exp_hi']:.2f}, P {STATS['h_Exp_p']:.4f}), "
      f"96 h {obs[1]:.2f} ({STATS['h_96h_lo']:.2f} to {STATS['h_96h_hi']:.2f}, P {STATS['h_96h_p']:.4f}); "
      f"change {d_obs:.2f} ({d_lo:.2f} to {d_hi:.2f}), rotation P {d_p:.4f}")

if NORM != "CFU":
    fig.text(0.985, 0.995, f"DESeq2 normalization: {NORM}", ha="right", va="top", fontsize=HS.BASE, color="#6b6a66")
if not ALLGENES:
    pd.DataFrame([STATS]).to_csv(OUT + f"fig4_summary_{NORM}.csv", index=False)
HS.check(fig)
fig.savefig(OUT + OUTNAME + ".pdf", dpi=600)
fig.savefig(OUT + OUTNAME + ".png", dpi=300)
print("saved")
