"""Inputs of Supplementary Fig. 2 (the method figure), from the Python implementation of Bactoloops:
  - every intermediate image and the candidate table for E. coli WT 96 h and B. subtilis WT (--dump-intermediates);
  - all 145 scored candidates of WT 96 h, low-confidence ones kept (keepLowConfidence = True), as scored_macro_wt_96h_repavg.csv.
Output: FIJI_output/python_runs/"""
import os, shutil, subprocess, sys, tempfile
import os as _hns_os
def _hns_root():
    r = _hns_os.environ.get("HNS_ROOT")
    if r:
        return _hns_os.path.join(r, "")
    d = _hns_os.path.dirname(_hns_os.path.abspath(__file__))
    while not _hns_os.path.exists(_hns_os.path.join(d, ".hns_root")) and _hns_os.path.dirname(d) != d:
        d = _hns_os.path.dirname(d)
    return _hns_os.path.join(d, "")
B = _hns_root() + "Bactoloops/"; OUT = B + "FIJI_output/python_runs/"; os.makedirs(OUT, exist_ok=True)
BL = B + "python/bactoloops.py"
for tif, rc, rid in [("iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif", 368, "wt_96h_repavg"),
                     ("iced_matrix_conversion/Dugar_HiC/Dugar_WT_st_merged_iced.tif", 0, "bsub_WT_st_merged")]:
    subprocess.run([sys.executable, BL, B + tif, "--recenter-bin", str(rc), "--out", OUT, "--id", rid, "--dump-intermediates"], check=True)
with tempfile.TemporaryDirectory() as t:     # keep low-confidence candidates: run with keepLowConfidence = True
    code = open(BL).read().replace("keepLowConfidence = False", "keepLowConfidence = True", 1)
    keep = os.path.join(os.path.dirname(BL), "_bactoloops_keeplow.py")
    open(keep, "w").write(code)
    try:
        subprocess.run([sys.executable, keep, B + "iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif",
                        "--recenter-bin", "368", "--out", t, "--id", "keep"], check=True)
    finally:
        os.remove(keep)
    shutil.copy(os.path.join(t, "loop_calls_keep.csv"), OUT + "scored_macro_wt_96h_repavg.csv")
print("Supplementary Fig. 2 inputs in", OUT)
