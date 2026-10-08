"""Build the contact maps every loop caller in the comparison runs on.

Datasets (10 kb bins):
  wt96   E. coli WT 96 h, rep1 + rep2 merged (this study)
  wtexp  E. coli WT exponential, rep1
  hns96  E. coli delta-hns 96 h, rep1 (cWX2870)
  stpA96 E. coli delta-stpA 96 h, rep1 (cWX3338)
  fis96  E. coli delta-fis 96 h, rep1 (cWX2478)
  ihf96  E. coli delta-ihf 96 h, rep1 (cWX2943)
  bsub   B. subtilis WT stationary, GEO runs 1 + 2 merged (Dugar et al. 2022)

For each: inputs/<name>.cool with raw counts and cooler-balanced weights (what cooltools,
Chromosight, Mustache and Peakachu read), and inputs/<name>.short.txt.gz, the same counts in
juicer "short with score" format for building a .hic (HiCCUPS, SIP).

E. coli counts are raw integers since 2026-09-29: cool and bed files/raw_counts/, built from
HiC-Pro's allValidPairs by build_raw_cools.py (identical to HiC-Pro's raw matrices). The merged
WT 96 h map is the sum of the two replicates' raw counts. B. subtilis counts are raw integers.
(Before that, the E. coli inputs held HiC-Pro's ICE values in the count column; those runs are in
previous_runs_2026-09-29_ice_as_counts/.)

With --write-tiff, also writes the merged WT 96 h ICE matrix as a TIFF for the ImageJ macro,
recentered on bin 368 like the other E. coli TIFFs (iced_matrix_conversion/
cWX2200_96h_rep1rep2_merged_10000_iced.tif). Off by default, so rebuilding the caller inputs never
changes the macro's input map.
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
import gzip
import os
import sys
import numpy as np
import pandas as pd
import cooler
import tifffile

ROOT = HNS_ROOT + "Bactoloops/"
C = ROOT + "cool and bed files/"
OUT = ROOT + "Comparisons/inputs/"
sys.path.insert(0, ROOT + "iced_matrix_conversion/Dugar_HiC")
from merge_and_ice import ice          # noqa: E402

RAW = C + "raw_counts/"
SOURCES = {
    "wt96": [RAW + "wt96_rep1_raw_10000.cool", RAW + "wt96_rep2_raw_10000.cool"],
    "wtexp": [RAW + "wtexp_rep2_raw_10000.cool"],
    "wt96r1": [RAW + "wt96_rep1_raw_10000.cool"],   # single replicates, for checking the merged map
    "wt96r2": [RAW + "wt96_rep2_raw_10000.cool"],
    "hns96": [RAW + "hns96_rep1_raw_10000.cool"],
    "bsub": [C + "Dugar_WT_st_merged_raw.cool"],
    "stpA96": [RAW + "stpA96_rep1_raw_10000.cool"],
    "fis96": [RAW + "fis96_rep1_raw_10000.cool"],
    "ihf96": [RAW + "ihf96_rep1_raw_10000.cool"],
}
ECOLI_BINS = RAW + "wt96_rep1_raw_10000.cool"      # bin table for maps read from a .matrix


def read_map(path):
    if not path.endswith(".matrix"):
        return cooler.Cooler(path).matrix(balance=False)[:].astype(float)
    t = pd.read_csv(path, sep="\t", header=None)
    n = cooler.Cooler(ECOLI_BINS).info["nbins"]
    m = np.zeros((n, n))
    m[t[0] - 1, t[1] - 1] = t[2]
    m[t[1] - 1, t[0] - 1] = t[2]
    return m
RECENTER = 368


def write(name, srcs):
    c0 = cooler.Cooler(ECOLI_BINS if srcs[0].endswith(".matrix") else srcs[0])
    bins = c0.bins()[:][["chrom", "start", "end"]]
    m = sum(read_map(s) for s in srcs)
    iu = np.triu_indices(len(m))
    keep = m[iu] > 0
    pix = pd.DataFrame({"bin1_id": iu[0][keep], "bin2_id": iu[1][keep], "count": m[iu][keep]})
    integer = np.allclose(pix["count"], np.round(pix["count"]))
    if integer:
        pix["count"] = np.round(pix["count"]).astype(int)
    uri = OUT + name + ".cool"
    cooler.create_cooler(uri, bins, pix, dtypes={"count": int if integer else float}, ordered=True)
    clr = cooler.Cooler(uri)
    cooler.balance_cooler(clr, store=True)          # cooler defaults: mad-max 5, min-nnz 10, ignore-diags 2
    chrom = bins.chrom.iloc[0]
    with gzip.open(OUT + name + ".short.txt.gz", "wt") as fh:   # str chr pos frag str chr pos frag score
        for a, b, v in zip(pix.bin1_id, pix.bin2_id, pix["count"]):
            fh.write(f"0 {chrom} {bins.start.iloc[a]} 0 0 {chrom} {bins.start.iloc[b]} 1 {v}\n")
    with open(OUT + name + ".chrom.sizes", "w") as fh:
        fh.write(f"{chrom}\t{clr.chromsizes[chrom]}\n")
    print(f"{name}: {len(m)} bins, {m[iu].sum():,.0f} contacts, integer counts: {integer}, chrom {chrom}")
    return m


os.makedirs(OUT, exist_ok=True)
WRITE_TIFF = "--write-tiff" in sys.argv
only = [a for a in sys.argv[1:] if not a.startswith("--")]   # optional: build just these datasets
for name, srcs in SOURCES.items():
    if only and name not in only:
        continue
    m = write(name, srcs)
    if name == "wt96" and WRITE_TIFF:
        iced = ice(m, filter_fraction=0.0)              # HiC-Pro ICE on these maps filters no bins
        tif = np.roll(np.roll(iced, -(RECENTER - 1), 0), -(RECENTER - 1), 1).astype(np.float32)
        tifffile.imwrite(ROOT + "iced_matrix_conversion/cWX2200_96h_rep1rep2_merged_10000_iced.tif", tif)
        avg = tifffile.imread(ROOT + "iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif")
        ok = (tif > 0) & (avg > 0)
        print(f"  merged ICE TIFF vs the rep1/rep2 ICE average: r = {np.corrcoef(tif[ok], avg[ok])[0, 1]:.6f}")
