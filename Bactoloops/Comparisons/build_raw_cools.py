"""Raw-count 10 kb .cool files for the E. coli Hi-C samples, from HiC-Pro's allValidPairs.

Input:  Hi-C raw count files/ (from the collaborators): allValidPairs, 10 kb BED, ICE biases, and
        ICED matrices for four samples.
Output: cool and bed files/raw_counts/<name>_raw_10000.cool, one per sample, plus the WT 96 h
        rep1 + rep2 sum. Each holds integer counts; bins carry
          weight        -- cooler's ICE balance (cooler.balance_cooler defaults), what cooltools,
                           mustache etc. read as "balanced";
          hicpro_weight -- 1 / HiC-Pro's ICE bias, so raw * w_i * w_j reproduces HiC-Pro's ICED matrix.

Binning: HiC-Pro assigns a read at 1-based position p to bin (p - 1) // 10000. With that rule the
matrix built here equals HiC-Pro's raw matrix exactly (checked below as ICED * b_i * b_j, which is
integer, for every sample with an ICED matrix; for the others the check is that raw / (b_i b_j) has
constant row sums, which is what ICE converges to).
Genome: EcW3110del = W3110 (NC_007779) minus 1,398,953-1,438,196; 4,607,088 bp, 461 bins.
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
import os
import numpy as np
import pandas as pd
import cooler

ROOT = HNS_ROOT + "Bactoloops/"
SRC = ROOT + "Hi-C raw count files/"
OUT = ROOT + "cool and bed files/raw_counts/"
N, BIN, CHROM, GLEN = 461, 10_000, "EcW3110del", 4_607_088
SAMPLES = {   # short name: HiC-Pro sample
    "wt96_rep1": "cWX2200_96h_Azure_10percent_Sau3AI_rep1",
    "wt96_rep2": "cWX2200_96h_Azure_10percent_Sau3AI_rep2",
    "wtexp_rep2": "cWX2200_exp_Azure_10percent_Sau3AI_rep2",
    "hns96_rep1": "cWX2870_96h_Azure_10percent_Sau3AI_rep1",
    "stpA96_rep1": "cWX3338_96h_Sau3AI_rep1",
    "fis96_rep1": "cWX2478_96HAzure10percent_Sau3AI_rep1",
    "ihf96_rep1": "cWX2943_96HAzure10percent_Sau3AI_rep1",
}
os.makedirs(OUT, exist_ok=True)
bins = pd.DataFrame({"chrom": CHROM, "start": np.arange(N) * BIN, "end": np.minimum(np.arange(1, N + 1) * BIN, GLEN)})


def raw_matrix(sample):
    M = np.zeros((N, N), dtype=np.int64)
    n = 0
    for ch in pd.read_csv(SRC + f"allValidPairs/{sample}.allValidPairs", sep="\t", header=None,
                          usecols=[1, 2, 4, 5], chunksize=5_000_000):
        assert set(ch[1]) | set(ch[4]) == {CHROM}
        a, b = (ch[2].values - 1) // BIN, (ch[5].values - 1) // BIN
        np.add.at(M, (np.minimum(a, b), np.maximum(a, b)), 1)
        n += len(ch)
    return M, n                                          # upper triangle incl. diagonal


def write_cool(name, U, hicpro_bias=None):
    i, j = np.nonzero(U)
    uri = OUT + f"{name}_raw_10000.cool"
    cooler.create_cooler(uri, bins, pd.DataFrame({"bin1_id": i, "bin2_id": j, "count": U[i, j]}),
                         dtypes={"count": np.int32}, ordered=True, symmetric_upper=True)
    clr = cooler.Cooler(uri)
    w, info = cooler.balance_cooler(clr, store=True)     # cooler's ICE -> "weight"
    if hicpro_bias is not None:
        with clr.open("r+") as f:
            f["bins"].create_dataset("hicpro_weight", data=1.0 / hicpro_bias)
    return uri, info


summary = []
for name, sample in SAMPLES.items():
    U, n = raw_matrix(sample)
    S = U + np.triu(U, 1).T
    b = np.loadtxt(SRC + f"biases/{sample}_10000_iced.matrix.biases")
    ice_ours = S / np.outer(b, b)
    rows = ice_ours.sum(1)
    check = f"ICE row sums (ours / HiC-Pro biases) spread {rows.std() / rows.mean():.1e}"
    fice = SRC + f"ICED/{sample}_10000_iced.matrix"
    if os.path.exists(fice):
        I = pd.read_csv(fice, sep="\t", header=None)
        ice = np.zeros((N, N)); ice[I[0] - 1, I[1] - 1] = I[2]
        diff = int(np.abs(np.round(ice * np.outer(b, b)) - S).sum())
        check += f"; vs HiC-Pro ICED: {diff} counts differ"
    uri, info = write_cool(name, U, b)
    print(f"{name}: {n:,} pairs, {U.sum():,} in matrix; {check}; cooler balance converged={info['converged']}")
    summary.append((name, sample, n))
    if name == "wt96_rep1":
        U1 = U
    if name == "wt96_rep2":
        uri, info = write_cool("wt96_merged", U1 + U)        # rep1 + rep2 summed; cooler weights only
        print(f"wt96_merged: {(U1 + U).sum():,} pairs; cooler balance converged={info['converged']}")

pd.DataFrame(summary, columns=["name", "hicpro_sample", "valid_pairs"]).to_csv(OUT + "samples.csv", index=False)
print("written to", OUT)
