"""House style for the H-NS paper figures (Xindan's print-size rules, 2026-10-07).

    import house_style as HS
    HS.apply()                                   # rcParams: Arial, 8 pt base, boxed axes, gene names in Arial italic
    fig = plt.figure(figsize=HS.size(500, 600))  # canvas in points (max 500 wide x 600 tall)
    HS.letter(ax, "a")                           # 12 pt bold panel letter
    ax.set_title("WT 96 h", fontsize=HS.TITLE)   # 10 pt second-level heading
    HS.check(fig)                                # assert canvas and font limits before saving

Rules: Arial throughout; gene names italic ($\\it{hns}$ renders in Arial italic); panel letters bold 12 pt; second-level headings
10 pt; nothing below 8 pt; ori and ter italic; canvas at final print size, at most 500 x 600 pt. A figure that does not fit is
split across pages (same figure number, panel letters continue).
"""
import matplotlib.pyplot as plt

LETTER, TITLE, BASE = 12, 10, 8           # pt
MAX_W, MAX_H = 500, 600                   # pt
FONT = "Arial"


def apply():
    plt.rcParams.update({
        "font.family": FONT, "font.size": BASE, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic", "mathtext.bf": f"{FONT}:bold",
        "mathtext.sf": FONT, "mathtext.default": "regular",
        "axes.titlesize": BASE, "axes.labelsize": BASE, "xtick.labelsize": BASE, "ytick.labelsize": BASE,
        "legend.fontsize": BASE, "figure.titlesize": TITLE,
        "axes.linewidth": 0.6, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.pad": 2, "ytick.major.pad": 2,
        "axes.spines.top": True, "axes.spines.right": True, "legend.frameon": False,
    })


def size(w_pt, h_pt):
    """figsize in inches from a canvas in points"""
    assert w_pt <= MAX_W and h_pt <= MAX_H, (w_pt, h_pt)
    return w_pt / 72, h_pt / 72


def letter(ax, s, dx=-0.1, dy=1.02, **kw):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=LETTER, fontweight="bold", va="bottom", ha="right", **kw)


def fig_letter(fig, s, x, y, **kw):
    fig.text(x, y, s, fontsize=LETTER, fontweight="bold", va="top", ha="left", **kw)


def check(fig, min_pt=BASE):
    """canvas within 500 x 600 pt and every visible text >= 8 pt (letters may be 12)"""
    w, h = fig.get_size_inches() * 72
    assert w <= MAX_W + 0.5 and h <= MAX_H + 0.5, f"canvas {w:.0f} x {h:.0f} pt"
    texts = list(fig.findobj(plt.Text))
    for ax in fig.axes:
        lg = ax.get_legend()
        if lg is not None:
            texts += list(lg.get_texts()) + [lg.get_title()]
    small = sorted({round(t.get_fontsize(), 2) for t in texts
                    if t.get_visible() and t.get_text().strip() and t.get_fontsize() < min_pt - 1e-6})
    assert not small, f"text below {min_pt} pt: {small}"
    return w, h
