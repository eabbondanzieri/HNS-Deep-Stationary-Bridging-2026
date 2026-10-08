"""1-kb ChIP/input tracks of Supplementary Fig. 5 (old Fig. S4): per-base coverage of ChIP and input each divided by its
total, binned to 1 kb (mean), ChIP / input. Writes ChIP-seq_data/s5_tracks_1kb.csv (genome order, lab genome)."""
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
import numpy as np, pandas as pd, sys
D = HNS_ROOT + "Bactoloops/ChIP-seq_data/"
def cov(f): return pd.read_csv(D + f, usecols=[1]).iloc[:, 0].values.astype(float)
def track(chip, inp, bs=1000):
    c, i = cov(chip), cov(inp); c, i = c / c.sum(), i / i.sum()
    n = -(-len(c) // bs); pad = n * bs - len(c)
    cb = np.r_[c, np.full(pad, np.nan)].reshape(n, bs); ib = np.r_[i, np.full(pad, np.nan)].reshape(n, bs)
    return np.nanmean(cb, 1) / np.nanmean(ib, 1)
PAIRS = {k: v for k, v in [l.split("=") for l in sys.argv[1:]]} if len(sys.argv) > 1 else None

TRACKS = [("anti-H-NS, dhns, exp", "anti-HNS_cWX2870_Exp_rep2_R1_clean.csv", "WGS_cWX2870_Exp_rep2_R1_clean.csv"),
          ("anti-H-NS, WT, 96 h", "anti-HNS_cWX2200_96hr_rep2_R1_clean.csv", "WGS_cWX2200_96hr_rep2_R1_clean.csv"),
          ("anti-H-NS, dhns, 96 h", "anti-HNS_cWX2870_96hr_rep2_R1_clean.csv", "WGS_cWX2870_96hr_rep2_R1_clean.csv"),
          ("anti-mCherry, hns-PAmCherry, exp", "anti-mCherry_cWX2882_Exp_rep1_R1_clean.csv", "WGS_cWX2882_Exp_rep1_R1_clean.csv"),
          ("anti-mCherry, WT (no mCherry), exp", "anti_mCherry_cWX2200_Exp_rep1_R1_clean.csv", "WGS_cWX2200_Exp_rep2_R1_clean.csv"),
          ("anti-mCherry, hns-PAmCherry, 96 h", "anti-mCherry_cWX2882_96hr_rep1_R1_clean.csv", "WGS_cWX2882_96hr_rep1_R1_clean.csv"),
          ("anti-mCherry, WT (no mCherry), 96 h", "anti_mCherry_cWX2200_96hr_rep1_R1_clean.csv", "WGS_cWX2200_96hr_rep2_R1_clean.csv")]
if __name__ == "__main__":
    cols = {lab: track(c, i) for lab, c, i in TRACKS}
    n = len(next(iter(cols.values())))
    out = pd.DataFrame({"bin_start_kb": np.arange(n), **{k: np.round(v, 4) for k, v in cols.items()}})
    out.to_csv(D + "s5_tracks_1kb.csv", index=False); print(out.shape)
