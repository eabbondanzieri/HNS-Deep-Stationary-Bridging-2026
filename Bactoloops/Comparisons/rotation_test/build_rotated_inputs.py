"""Rotated copies of the comparison inputs, for the rotation test of every loop caller.

A circular chromosome has no natural first bin, so a caller's loops should not depend on where
the contact map starts. For each dataset and offset k, the raw count matrix is rotated by k bins
(bin k of the genome becomes bin 0) and written exactly as build_inputs.py writes the originals:
a .cool with raw counts and cooler's balancing weights (cooltools, Chromosight, Mustache,
Peakachu), the juicer "short" file for the .hic (HiCCUPS, SIP), and chrom.sizes; plus, for the
ImageJ macro, the ICE TIFF rotated the same way. Calls on a rotated map are mapped back to genome
bins by adding k (mod n).

Writes rotation_test/inputs/<dataset>_r<k>.{cool,short.txt.gz,chrom.sizes,tif}.
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
import numpy as np
import pandas as pd
import cooler
import tifffile

ROOT = HNS_ROOT + "Bactoloops/"
IN = ROOT + "Comparisons/inputs/"
OUT = ROOT + "Comparisons/rotation_test/inputs/"
SHIFTS = {"wt96": [0, 1, 2, 3, 57, 118, 230, 333, 400], "bsub": [0, 1, 2, 3, 97, 150, 211, 300]}
TIF = {"wt96": (ROOT + "iced_matrix_conversion/cWX2200_96h_rep1rep2_merged_10000_iced.tif", 368),
       "bsub": (ROOT + "iced_matrix_conversion/Dugar_HiC/Dugar_WT_st_merged_iced.tif", 0)}

os.makedirs(OUT, exist_ok=True)
for d, ks in SHIFTS.items():
    c0 = cooler.Cooler(IN + d + ".cool")
    bins = c0.bins()[:][["chrom", "start", "end"]]
    m0 = c0.matrix(balance=False)[:]
    n = len(m0)
    chrom = bins.chrom.iloc[0]
    tif0 = tifffile.imread(TIF[d][0])          # genome bin g sits at TIFF bin (g - (rc - 1)) mod n
    rc0 = TIF[d][1]
    for k in ks:
        name = f"{d}_r{k}"
        m = np.roll(m0, -k, axis=(0, 1))       # rotated bin i = genome bin (i + k) mod n
        iu = np.triu_indices(n)
        keep = m[iu] > 0
        pix = pd.DataFrame({"bin1_id": iu[0][keep], "bin2_id": iu[1][keep], "count": m[iu][keep].astype(int)})
        uri = OUT + name + ".cool"
        cooler.create_cooler(uri, bins, pix, dtypes={"count": int}, ordered=True)
        cooler.balance_cooler(cooler.Cooler(uri), store=True)
        with gzip.open(OUT + name + ".short.txt.gz", "wt") as fh:
            for a, b, v in zip(pix.bin1_id, pix.bin2_id, pix["count"]):
                fh.write(f"0 {chrom} {bins.start.iloc[a]} 0 0 {chrom} {bins.start.iloc[b]} 1 {v}\n")
        with open(OUT + name + ".chrom.sizes", "w") as fh:
            fh.write(f"{chrom}\t{c0.chromsizes[chrom]}\n")
        # TIFF for the macro: TIFF bin t = genome bin (t + k) mod n, i.e. recenter bin k + 1
        g0 = (rc0 - 1) if rc0 else 0
        tifffile.imwrite(OUT + name + ".tif", np.roll(tif0, -(k - g0), axis=(0, 1)).astype(np.float32))
        print(name, "done")
