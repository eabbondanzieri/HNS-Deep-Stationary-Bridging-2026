"""Supplementary Figure 2 (rebuilt): how Bactoloops calls loops, step by step, on E. coli WT 96 h.

Intermediate images come from the Python implementation run with --dump-intermediates (python/bactoloops.py;
identical to the ImageJ macro's steps); the final calls (panel h) are the macro's (FIJI_output/). Zoomed
panels show genome 2,000-3,000 kb on both axes (the region enlarged in Fig. 1).

  a  input map (row sums scaled, diagonal smoothed) padded periodically (circular chromosome)
  b  zoom: input;  c  rolling-ball background (100 kb);  d  background-subtracted, ridge mask outlined
  e  loop density (ridge masked, noise floor removed, smoothed), thresholded; particles (watershed-split)
     and the candidate loops they give
  f  distance-matched confidence -- option A: one loop against every pixel at its distance;
     option B: every candidate against the 95th/99th percentile at each distance
  g  near-diagonal (ridge-region) candidates: the HiCCUPS-style donut test (kernels; one E. coli candidate
     that fails, B. subtilis RC2-RC3 that passes)
  h  final calls and the integrated loop density (ILD)
Usage: python figureS2_method.py [A|B]   (default: both) -> FigureS2_method_option{A,B}.pdf/.png
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
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

ROOT = HNS_ROOT + "Bactoloops/"
DUMP = HNS_ROOT + "Bactoloops/FIJI_output/python_runs/"
OUT = ROOT + "updated figures/"
sys.path.insert(0, ROOT + "python")
from bactoloops import refine_pick                                       # noqa: E402

E = np.load(DUMP + "intermediates_wt_96h_repavg.npz")
B = np.load(DUMP + "intermediates_bsub_WT_st_merged.npz")
N, PH, UP = int(E["n"]), int(E["pad_half"]), int(E["upsample"])
RC = 368
tiff_of = lambda g: (g - (RC - 1)) % N
T0, T1 = tiff_of(200), tiff_of(300) + 1                                    # genome 2,000-3,000 kb, in TIFF bins
Z = slice(T0 + PH, T1 + PH)                                              # the same window in the padded images
INK, MUTED = "#1f1f1e", "#6b6a66"
RIDGE_C, CAND_C = "#00c8e0", "#1a2fb0"                              # ridge outline as in the macro
CONF_C = {"high": "#1a2fb0", "medium": "#5b9bd5", "low": "#bdbdbd"}   # blues read on the Fire LUT (the macro circles loops in blue)
# maps before background subtraction: the collaborators' Hi-C scheme (linear white-blue-purple-black ramp,
# map / row sum clipped to 0-0.01, diagonal and +/-1 off-diagonals set to the maximum)
COLLAB = LinearSegmentedColormap.from_list("collab", np.array([
    (255, 255, 255), (162, 192, 222), (140, 137, 187), (140, 87, 167), (140, 45, 143),
    (120, 20, 120), (90, 15, 90), (60, 10, 60), (30, 5, 30), (0, 0, 0)]) / 255, N=256)
ZUPPER = 0.01
# maps after background subtraction: ImageJ's Fire LUT, inverted, 0-displayMax -- as the macro displays them
_fr = [0, 0, 1, 25, 49, 73, 98, 122, 146, 162, 173, 184, 195, 207, 217, 229, 240, 252, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255]
_fg = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 14, 35, 57, 79, 101, 117, 133, 147, 161, 175, 190, 205, 219, 234, 248, 255, 255, 255, 255]
_fb = [0, 61, 96, 130, 165, 192, 220, 227, 210, 181, 151, 122, 93, 64, 35, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 35, 98, 160, 223, 255, 255, 255]
_x = np.linspace(0, 255, 32)
FIRE = np.stack([np.interp(np.arange(256), _x, c) for c in (_fr, _fg, _fb)], 1) / 255
FIRE_INV = ListedColormap(FIRE[::-1], name="fire_inverted")
DISPLAY_MAX = 3000                                                       # the macro's displayMax
import sys as _hs_sys
_hs_sys.path.insert(0, HNS_ROOT + "Bactoloops/updated figures/")
import house_style as HS    # noqa: E402  Arial, 8 pt minimum, 12 pt bold letters, 10 pt titles
HS.apply()


def letter(ax, s, dx=-0.08, dy=1.03):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=HS.LETTER, fontweight="bold", va="bottom", ha="right")


def zoom_axes(ax, title):
    ticks = [0, 50, 100]
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    ax.set_xticklabels(["2.0", "2.5", "3.0"]); ax.set_yticklabels([])
    ax.set_title(title, fontsize=HS.BASE, pad=2)


def collab_norm(img, diag_fill=True):
    """the collaborators' display: padded-image values / row sum, diagonal and +/-1 set to the cut-off"""
    core = E["scaled_padded"][PH:PH + N, PH:PH + N]
    m = img / np.median(core.sum(1))
    if diag_fill:
        pi, pj = np.indices(m.shape)
        d = (pi - pj) % N
        m = np.where((d == 0) | (d == 1) | (d == N - 1), ZUPPER, m)
    return np.clip(m, 0, ZUPPER)


def show_collab(ax, img, extent=None):
    ax.imshow(img, cmap=COLLAB, vmin=0, vmax=ZUPPER,
              origin="lower", interpolation="nearest", extent=extent)


# ---------------- confidence null (every pixel at a distance, refined the same way as the loop) ----------------
refA, avgA = E["refine_map"], E["log2oe_corrected"]


def null_values(d):
    ks = np.arange(N)
    out = np.empty(N)
    for k in ks:
        kc = (k + d) % N
        s1, s2 = refine_pick(k, kc, refA, N, 1, 3)
        out[k] = avgA[(k + s1) % N, (kc + s2) % N]
    return out


circ = lambda a, b: min(abs(a - b) % N, N - abs(a - b) % N)
scored = pd.read_csv(DUMP + "scored_macro_wt_96h_repavg.csv")   # the macro run with keepLowConfidence = true: all 145 candidates
final = pd.read_csv(ROOT + "FIJI_output/loop_calls_wt_96h_repavg.csv")
cands = pd.read_csv(DUMP + "candidates_wt_96h_repavg.csv")
inz = lambda L: L[L.bin1.between(T0, T1 - 1) & L.bin2.between(T0, T1 - 1)]

# ---------------- donut test ----------------
def donut_ratios(D, i, j, p=2, w=5):
    n, ph = int(D["n"]), int(D["pad_half"])
    core = D["scaled_padded"][ph:ph + n, ph:ph + n].astype(float)
    expd = D["expected"]
    rows, cols = np.arange(i - w, i + w + 1), np.arange(j - w, j + w + 1)
    di, dj = np.meshgrid(rows - i, cols - j, indexing="ij")
    O = core[np.ix_(rows % n, cols % n)]
    dd = np.minimum(abs(rows[:, None] % n - cols[None, :] % n), n - abs(rows[:, None] % n - cols[None, :] % n))
    Ed = expd[dd]
    ok = (O > 0) & (dd >= 2)
    peak = (abs(di) <= p) & (abs(dj) <= p)
    kern = {"donut": ~peak & (di != 0) & (dj != 0), "lower-left": (di > 0) & (dj < 0) & ~peak,
            "horizontal": (abs(di) <= 1) & (abs(dj) > p), "vertical": (abs(dj) <= 1) & (abs(di) > p)}
    d0 = circ(i, j)
    return {k: core[i % n, j % n] / (O[m & ok].sum() / Ed[m & ok].sum() * expd[d0]) for k, m in kern.items()}, kern


def ecoli_ridge_fail():
    """a near-diagonal local maximum in the zoom window that reaches the 95th percentile but fails the donut test"""
    ridge_in = E["ridge_keep_mask"][PH:PH + N, PH:PH + N] < 0.5
    best = None
    for i in range(T0, T1):
        for d in range(3, 25):
            j = i + d
            if j >= T1 or not ridge_in[i, j]:
                continue
            v = refA[i, j]
            if v < refA[i - 1:i + 2, j - 1:j + 2].max():
                continue
            pct = np.mean(null_values(d) < avgA[i, j])
            if pct < 0.95:
                continue
            r, _ = donut_ratios(E, i, j)
            ok = all(x >= 1.5 for x in r.values()) and max(r["donut"], r["lower-left"]) >= 1.75
            if not ok and (best is None or pct > best[2]):
                best = (i, j, pct, r)
    return best


def build(option):
    fig = plt.figure(figsize=HS.size(500, 600))
    g = fig.add_gridspec(4, 1, height_ratios=[1, 1.15, 0.85, 1.0], hspace=0.5, left=0.07, right=0.985, top=0.955, bottom=0.06)
    r1 = g[0].subgridspec(1, 4, wspace=0.22)
    sp = E["scaled_padded"]
    vmin, vmax = np.percentile(sp[sp > 0], 5), np.percentile(sp[sp > 0], 99.7)

    # a: padded input
    ax = fig.add_subplot(r1[0, 0])
    SPN = collab_norm(sp)
    show_collab(ax, SPN)
    ax.add_patch(Rectangle((PH - 0.5, PH - 0.5), N, N, fill=False, ec=INK, lw=0.8, ls=(0, (3, 2))))
    ax.add_patch(Rectangle((T0 + PH - 0.5, T0 + PH - 0.5), T1 - T0, T1 - T0, fill=False, ec="#1f5fa8", lw=1.0))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("padded map", fontsize=HS.BASE, pad=2)
    letter(ax, "a", dx=-0.04)
    # b, c, d: zooms
    ax = fig.add_subplot(r1[0, 1]); show_collab(ax, SPN[Z, Z]); zoom_axes(ax, "input (zoom)"); letter(ax, "b", dx=-0.04)
    ax.set_xlabel("Genome position (Mb)")
    ax.set_ylabel("")
    ax = fig.add_subplot(r1[0, 2]); show_collab(ax, collab_norm(E["rolling_ball_bg"], diag_fill=False)[Z, Z])
    zoom_axes(ax, "background"); letter(ax, "c", dx=-0.04)
    ax = fig.add_subplot(r1[0, 3])
    rb = E["bg_subtracted"][Z, Z] - (float(E["shot_noise_median"]) + 0.5 * float(E["shot_noise_sd"]))   # noise floor, as the display map
    ax.imshow(np.clip(rb, 0, None), cmap=FIRE_INV, vmin=0, vmax=DISPLAY_MAX, origin="lower", interpolation="nearest")
    ax.contour(E["ridge_keep_mask"][Z, Z] < 0.5, levels=[0.5], colors=[RIDGE_C], linewidths=0.8)
    zoom_axes(ax, "background removed"); letter(ax, "d", dx=-0.04)

    # e: particles and candidates
    r2 = g[1].subgridspec(1, 2, wspace=0.45, width_ratios=[1, 1.6])
    ax = fig.add_subplot(r2[0, 0])
    ld = E["loop_density_smoothed_padded"][Z, Z]
    ax.imshow(ld, cmap=FIRE_INV, vmin=0, vmax=DISPLAY_MAX, origin="lower",
              interpolation="nearest", extent=[-0.5, T1 - T0 - 0.5, -0.5, T1 - T0 - 0.5])
    pb = E["particles_binary"][(T0 + PH) * UP:(T1 + PH) * UP, (T0 + PH) * UP:(T1 + PH) * UP]
    ax.contour(np.arange(pb.shape[1]) / UP - 0.5 + 0.5 / UP, np.arange(pb.shape[0]) / UP - 0.5 + 0.5 / UP, pb.astype(float),
               levels=[0.5], colors=[CAND_C], linewidths=0.4)
    c = inz(cands)
    ax.plot(c.bin2 - T0, c.bin1 - T0, "+", color="#00c8e0", ms=3, mew=0.6)
    ax.plot(c.bin1 - T0, c.bin2 - T0, "+", color="#00c8e0", ms=3, mew=0.6)
    zoom_axes(ax, "candidates"); ax.set_xlabel("Genome position (Mb)"); letter(ax, "e", dx=-0.04)
    ax.set_xlim(-0.5, T1 - T0 - 0.5); ax.set_ylim(-0.5, T1 - T0 - 0.5)

    # f: confidence
    ax = fig.add_subplot(r2[0, 1])
    if option == "A":
        z = inz(final).sort_values("integrated_excess", ascending=False).iloc[0]
        b1, b2 = int(round(z.bin1)), int(round(z.bin2))
        d = circ(b1, b2)
        nv = null_values(d)
        val = avgA[b1, b2]
        ax.hist(nv, bins=40, color="#c9c7c1", edgecolor="white", lw=0.3)
        for q, ls, lab in [(95, (0, (3, 2)), "95th"), (99, "-", "99th")]:
            x = np.percentile(nv, q); ax.axvline(x, color=MUTED, lw=0.7, ls=ls)
            ax.text(x, ax.get_ylim()[1] * 0.97, f" {lab}", fontsize=HS.BASE, color=MUTED, va="top")
        ax.annotate("loop", xy=(val, 0), xytext=(val, ax.get_ylim()[1] * 0.45), ha="center", fontsize=HS.BASE, color=CONF_C["high"],
                    arrowprops=dict(arrowstyle="-|>", color=CONF_C["high"], lw=0.8))
        ax.set_xlabel("log$_2$ (obs/exp), background removed")
        ax.set_ylabel("pixels at the\nsame distance")
        ax.set_title(f"{d * 10}-kb loop vs all pixels at {d * 10} kb: {100 * np.mean(nv < val):.1f}th percentile", fontsize=HS.BASE, pad=2)
    else:
        s = scored.copy()
        s["d"] = [circ(int(round(a)), int(round(b))) for a, b in zip(s.bin1, s.bin2)]
        cache = {}
        for d in s.d.unique():
            nv = null_values(d); cache[d] = (np.median(nv), np.percentile(nv, 95))
        s["rel"] = [(v - cache[d][0]) / (cache[d][1] - cache[d][0]) for v, d in zip(s.log2_obs_exp, s.d)]
        ax.axhspan(-1, 1, color="#ecebe7", lw=0, zorder=0)
        ax.axhline(1, color=MUTED, lw=0.7, ls=(0, (3, 2))); ax.axhline(0, color=MUTED, lw=0.5, ls=":")
        ax.text(2330, 1, "95th pct", fontsize=HS.BASE, color=MUTED, va="bottom", ha="right")
        ax.text(2330, 0, "median", fontsize=HS.BASE, color=MUTED, va="bottom", ha="right")
        for conf in ["low", "medium", "high"]:
            k = s[s.confidence == conf]
            ax.scatter(k.d * 10, k.rel, s=7 if conf != "low" else 6, color=CONF_C[conf], lw=0, label=f"{conf} ({len(k)})", zorder=3)
        ax.set_ylim(-0.3, s.rel.max() * 1.08); ax.set_xlim(0, 2340)
        ax.set_xlabel("loop distance (kb)")
        ax.set_ylabel("enrichment over the distance-matched\nbackground (1 = 95th percentile)")
        ax.legend(frameon=False, fontsize=HS.BASE, loc="upper right", handletextpad=0.2, title="candidates", title_fontsize=HS.BASE)
        ax.set_title(f"all {len(s)} candidates vs pixels at the same distance", fontsize=HS.BASE, pad=2)
    letter(ax, "f", dx=-0.17)

    # g: donut test
    gg = g[2].subgridspec(1, 2, wspace=0.15, width_ratios=[1.0, 1.1])
    gk = gg[0].subgridspec(1, 4, wspace=0.25)
    _, kern = donut_ratios(B, 206, 227)
    for c_, (k, col) in enumerate(zip(["donut", "lower-left", "horizontal", "vertical"], ["#6baed6", "#fd8d3c", "#74c476", "#9e9ac8"])):
        ax = fig.add_subplot(gk[0, c_])
        img = np.ones((11, 11, 3))
        img[kern[k]] = matplotlib.colors.to_rgb(col)
        img[5, 5] = matplotlib.colors.to_rgb("#b2182b")
        ax.imshow(img, origin="upper", interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title({"lower-left": "lower\nleft", "horizontal": "horiz.", "vertical": "vert."}.get(k, k), fontsize=HS.BASE, pad=1.5)
        if c_ == 0:
            letter(ax, "g", dx=-0.25, dy=1.0)
    fail = ecoli_ridge_fail()
    passr, _ = donut_ratios(B, 206, 227)
    gz = gg[1].subgridspec(1, 2, wspace=0.6)
    DIVC = LinearSegmentedColormap.from_list("div", ["#1f5fa8", "#f7f7f7", "#b2182b"])
    examples = [(E, fail[0], fail[1], fail[3], f"$\\it{{E. coli}}$, {circ(fail[0], fail[1]) * 10} kb:\nfails"),
                (B, 206, 227, passr, "$\\it{B. subtilis}$ RC2–RC3:\npasses")]
    for c_, (D, i0, j0, r, title) in enumerate(examples):
        n, ph = int(D["n"]), int(D["pad_half"])
        core = D["scaled_padded"][ph:ph + n, ph:ph + n].astype(float)
        rows, cols = np.arange(i0 - 5, i0 + 6) % n, np.arange(j0 - 5, j0 + 6) % n
        dd = np.minimum(abs(rows[:, None] - cols[None, :]), n - abs(rows[:, None] - cols[None, :]))
        w = np.log2(np.where(core[np.ix_(rows, cols)] > 0, core[np.ix_(rows, cols)], np.nan) / D["expected"][dd])
        w = np.ma.masked_where(dd < 2, w)                       # the diagonal is not used
        ax = fig.add_subplot(gz[0, c_])
        cm = DIVC.copy(); cm.set_bad("#d9d9d9")
        ax.imshow(w, cmap=cm, vmin=-2, vmax=2, origin="upper", interpolation="nearest")
        ax.add_patch(Rectangle((2.5, 2.5), 5, 5, fill=False, ec=INK, lw=0.6, ls=(0, (2, 1))))
        ax.add_patch(Rectangle((4.5, 4.5), 1, 1, fill=False, ec=INK, lw=0.9))
        ax.set_xticks([]); ax.set_yticks([])
        ok = all(v >= 1.5 for v in r.values()) and max(r["donut"], r["lower-left"]) >= 1.75
        ax.set_title(title, fontsize=HS.BASE, pad=1.5, color=INK if ok else MUTED)
        ax.text(0.5, -0.05, "donut {:.2f}, LL {:.2f}\nH {:.2f}, V {:.2f}".format(
            r["donut"], r["lower-left"], r["horizontal"], r["vertical"]), transform=ax.transAxes, ha="center", va="top",
            fontsize=HS.BASE, color="#b2182b" if ok else MUTED)
    # pass criteria and the colour scale are given in the legend

    # h: final calls + ILD
    r3 = g[3].subgridspec(1, 2, wspace=0.25, width_ratios=[1, 2.3])
    ax = fig.add_subplot(r3[0, 0])
    disp = E["display"][T0:T1, T0:T1]
    ax.imshow(disp, cmap=FIRE_INV, vmin=0, vmax=DISPLAY_MAX, origin="lower", interpolation="nearest")
    ax.contour(E["ridge_keep_mask"][Z, Z] < 0.5, levels=[0.5], colors=[RIDGE_C], linewidths=0.6)
    f = inz(final)
    for conf in ["medium", "high"]:
        k = f[f.confidence == conf]
        for a_, b_ in [(k.bin2, k.bin1), (k.bin1, k.bin2)]:
            ax.scatter(a_ - T0, b_ - T0, s=70, facecolor="none", edgecolor=CONF_C[conf], lw=1.1)
    zoom_axes(ax, f"final calls ({len(f)} here)"); ax.set_xlabel("Genome position (Mb)"); letter(ax, "h", dx=-0.04)
    ax = fig.add_subplot(r3[0, 1])
    ild = pd.read_csv(ROOT + "FIJI_output/integrated_loop_density_wt_96h_repavg.csv").ILD_au.values
    ax.plot(np.arange(N), ild, color=INK, lw=0.7)
    ax.axvspan(T0, T1, color="#1f5fa8", alpha=0.08, lw=0)
    ax.set_xlim(0, N)
    tk = [tiff_of(k // 10) for k in [3670, 0, 1000, 2000, 3000]]
    ax.set_xticks(tk); ax.set_xticklabels([3670, 0, 1000, 2000, 3000])
    ax.set_xlabel("Genome position (kb)"); ax.set_ylabel("ILD (a.u.)")
    ax.set_title(f"integrated loop density (ILD); {len(final)} loops genome-wide", fontsize=HS.BASE, pad=2)
    HS.check(fig)
    fig.savefig(OUT + f"FigureS2_method_option{option}.pdf", dpi=600)
    fig.savefig(OUT + f"FigureS2_method_option{option}.png", dpi=300)
    print("saved option", option, "| donut fail example:", None if not fail else (fail[0], fail[1], round(fail[2], 3), {k: round(v, 2) for k, v in fail[3].items()}))


for opt in (sys.argv[1:] or ["A"]):
    build(opt)
