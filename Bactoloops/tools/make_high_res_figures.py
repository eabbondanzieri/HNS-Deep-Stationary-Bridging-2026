"""High-resolution copies of the figures that show Loop_density_smoothed_loops images (Fig. 2, Supplementary Figs. 4, 9
and 15), for production. The regular figures are kept small for email (raster panels resampled at 100-300 dpi).

Each figure script is run with two changes: large images (>= 400 px on a side, i.e. the loop maps and the circle panels)
are embedded at their own pixel resolution (interpolation "none": no resampling), and PDFs are saved at 600 dpi.
The copies are labelled like the collected figures and written to
  Paper Update/Main Figures/high_res/          Figure_2_..._high_res.pdf
  Paper Update/Supplementary Figures/high_res/ Supplementary_Figure_NN_..._high_res.pdf
after which each script is run again unchanged, so the regular figure files are as before."""
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
import os, shutil, subprocess, sys

B = HNS_ROOT + "Bactoloops/"
PU = B + "Paper Update/"
PY = sys.executable
TOOLS = B + "tools/"
# script, source PDF written by it, destination (folder, name), label
FIGS = [("updated figures/figure2.py", "updated figures/Figure2.pdf",
         ("Main Figures", "Figure_2_loop_anchors_and_HNS"), "Figure 2"),
        ("updated figures/figureS3_loop_conditions.py", "updated figures/FigureS3_conditions.pdf",
         ("Supplementary Figures", "Supplementary_Figure_04_loops_in_96h_mutants"), "Supplementary Figure 4"),
        ("updated figures/figureS7_timecourse.py", "updated figures/FigureS7_timecourse.pdf",
         ("Supplementary Figures", "Supplementary_Figure_09_looping_time_course"), "Supplementary Figure 9"),
        ("Supplemental_figure_Ecoli_vs_Bsub/fig_supplement_ecoli_vs_bsub.py", "Supplemental_figure_Ecoli_vs_Bsub/FigS_Ecoli_vs_Bsub_draft.pdf",
         ("Supplementary Figures", "Supplementary_Figure_15_Ecoli_vs_Bsubtilis"), "Supplementary Figure 15")]   # renumbered 2026-10-12

PATCH = r'''
import sys, runpy, numpy as np, matplotlib
matplotlib.use("Agg")
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
_imshow, _save, _pp = Axes.imshow, Figure.savefig, PdfPages.savefig
def imshow(self, X, *a, **k):
    if min(np.shape(X)[:2]) >= 400:
        k["interpolation"] = "none"                 # embed the image's own pixels, no resampling
    return _imshow(self, X, *a, **k)
def savefig(self, fname, *a, **k):
    if str(fname).lower().endswith(".pdf"):
        k["dpi"] = 600
    return _save(self, fname, *a, **k)
def ppsave(self, figure=None, **k):
    k["dpi"] = 600
    return _pp(self, figure, **k)
Axes.imshow, Figure.savefig, PdfPages.savefig = imshow, savefig, ppsave
sys.argv = [sys.argv[1]]
runpy.run_path(sys.argv[0], run_name="__main__")
'''


def run(script, patched):
    d, f = os.path.split(B + script)
    py = PY
    cmd = [py, "-c", PATCH, f] if patched else [py, f]
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"{script} failed:\n{p.stderr[-2000:]}")


for script, src, (folder, name), label in FIGS:
    run(script, patched=True)
    os.makedirs(PU + folder + "/high_res", exist_ok=True)
    dst = PU + folder + f"/high_res/{name}_high_res.pdf"
    subprocess.run(["swift", TOOLS + "add_label_pdf.swift", B + src, dst, label], check=True, capture_output=True)   # Arial Bold label
    print(f"{label}: {os.path.getsize(dst) / 1e6:.1f} MB -> {dst}")
    run(script, patched=False)                       # restore the regular (email-size) figure files
