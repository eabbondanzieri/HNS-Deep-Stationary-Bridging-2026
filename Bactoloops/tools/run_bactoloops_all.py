"""Run Bactoloops on every contact map used in the paper (outputs in Bactoloops/FIJI_output/, named as the paper's scripts expect).
Engine: the Python implementation (default; Bactoloops/python/bactoloops.py, reproduces the macro) or the ImageJ macro
(--engine fiji; needs Fiji, set FIJI=/path/to/fiji). Maps are 10-kb ICE TIFFs (iced_matrix_conversion/), recentred at bin 368 for
E. coli and not recentred (0) for B. subtilis.
Usage: python run_bactoloops_all.py [--engine python|fiji] [run_id ...]"""
import argparse, os, subprocess, sys
import os as _hns_os
def _hns_root():
    r = _hns_os.environ.get("HNS_ROOT")
    if r:
        return _hns_os.path.join(r, "")
    d = _hns_os.path.dirname(_hns_os.path.abspath(__file__))
    while not _hns_os.path.exists(_hns_os.path.join(d, ".hns_root")) and _hns_os.path.dirname(d) != d:
        d = _hns_os.path.dirname(d)
    return _hns_os.path.join(d, "")
HNS_ROOT = _hns_root()
B = HNS_ROOT + "Bactoloops/"
IC = "iced_matrix_conversion/"
RUNS = {  # run id: (map, recenter bin)
    "wt_Exp": (IC + "cWX2200_exp_Azure_10percent_Sau3A_rep1_10000_iced.tif", 368),
    "wt_24h": (IC + "cWX2200_24H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", 368),
    "wt_96h_rep1": (IC + "cWX2200_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif", 368),
    "wt_96h_rep2": (IC + "cWX2200_96H_Azure_10percent_Sau3AI_rep2_10000_iced.tif", 368),
    "wt_96h_repavg": (IC + "cWX2200_96h_rep1rep2_average_10000_iced.tif", 368),
    "HNSm_Exp": (IC + "cWX2870_exp_Azure_10percent_rep1_10000_iced.tif", 368),
    "HNSm_96h": (IC + "cWX2870_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", 368),
    "stpA_96h": (IC + "cWX3338_96h_Sau3AI_rep1_10000_iced.tif", 368),
    "stpA_96h_rep2": (IC + "cWX3338_96h_Sau3AI_rep2_10000_iced.tif", 368),
    "fis_96h": (IC + "cWX2478_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", 368),
    "IHF_96h": (IC + "cWX2943_96H_Azure_10percent_Sau3AI_rep1_10000_iced.tif", 368),
    "zapA_96h": (IC + "cWX3541_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif", 368),
    "ftsP_96h": (IC + "cWX3543_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif", 368),
    "envC_96h": (IC + "cWX3550_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif", 368),
    "hnsPAmCherry_96h": (IC + "cWX2882_96h_Sau3AI_rep1_10000_iced.tif", 368),
    "bsub_WT_st_merged": (IC + "Dugar_HiC/Dugar_WT_st_merged_iced.tif", 0),
    "bsub_rok_st_merged": (IC + "Dugar_HiC/Dugar_rok_st_merged_iced.tif", 0),
    "rok_st": (IC + "Dugar_HiC/GSM4288936_WT_st_HiC_iced.tif", 0),
}
# hns-PAmCherry vs WT at a common depth (Supplementary Note 2): maps from tools/common_depth_maps.py
for s, r in [("wt", 1), ("wt", 2), ("pam", 2)]:
    for d in range(3):
        RUNS[f"{s}96_rep{r}_common{d}"] = (IC + f"common_depth/{s}96_rep{r}_common{d}_10000_iced.tif", 368)
for n in ["pam96_rep3_common0", "wt96_avg2_common", "pam96_avg2_common", "wt96_rep1_full", "wt96_rep2_full", "pam96_rep2_full"]:
    RUNS[n] = (IC + f"common_depth/{n}_10000_iced.tif", 368)
OUTDIR = {k: ("FIJI_output/extra_runs/pamcherry_common/" if ("_common" in k or k.endswith("_full")) else "FIJI_output/") for k in RUNS}

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0]); ap.add_argument("--engine", choices=["python", "fiji"], default="python")
ap.add_argument("runs", nargs="*"); a = ap.parse_args()
todo = a.runs or list(RUNS)
macro = open(B + "HiC_PeakFindingAlgorithm_with_symmetrization.ijm").read()
i0 = macro.index('Dialog.create("HiC Loop Finding");'); i1 = macro.index("skipSaving      = Dialog.getCheckbox();") + len("skipSaving      = Dialog.getCheckbox();")
for k in todo:
    tif, rc = RUNS[k]; out = B + OUTDIR[k]; os.makedirs(out, exist_ok=True)
    if not os.path.exists(B + tif):
        print(f"skip {k}: {tif} not found"); continue
    if a.engine == "python":
        subprocess.run([sys.executable, B + "python/bactoloops.py", B + tif, "--recenter-bin", str(rc), "--bin-size", "10000", "--out", out, "--id", k], check=True)
    else:   # the macro's dialog replaced by fixed settings; one Fiji process per map
        body = macro[:i0] + f'sourceTitle = getTitle();\nsaveDir = "{out}";\nimageIdentifier = "{k}";\nrecenterBin = {rc};\nbinSizeBp = 10000;\nskipSaving = false;' + macro[i1:]
        ijm = out + f"run_{k}.ijm"
        open(ijm, "w").write(f'open("{B + tif}");\n' + body + f'\nFile.saveString(getInfo("log"), "{out}log_{k}.txt");\neval("script", "System.exit(0);");\n')
        subprocess.run([os.environ.get("FIJI", "fiji"), "--allow-multiple", "--no-splash", "--", "--headless", "-macro", ijm], check=True)
    print("done", k, flush=True)
