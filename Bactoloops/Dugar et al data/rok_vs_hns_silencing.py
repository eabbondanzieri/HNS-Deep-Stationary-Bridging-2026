"""Does silencing at silencer-bound regions strengthen when loops form -- Rok (B. subtilis; Dugar et al. 2022) vs H-NS
(E. coli; this study)? Same metric for both, in 10-kb bins (the ChIP bins):

  per bin   d = log2((mutant + 1) / (WT + 1)) of the bin's expression (positive = higher without the silencer)
  offset    mean d over silencer-bound bins minus the genome-wide mean d  (relative silencing; normalization cancels)
  change    offset(stationary) - offset(exponential)
  stats     95% CI: bootstrap over bound bins (2,000); P: exact rotation test (bound-bin set shifted by every offset >= 5
            bins; rotation_test.py), two-sided

B. subtilis: GSE144475 RNA-seq coverage bigWigs (one replicate per condition; exponential 5 h, stationary 22 h, LB),
summed per 10-kb bin; bound bins = Rok ChIP z >= 2 (rok_zscore_10kb.csv) or Dugar's eight Rok clusters (+/-1 bin).
E. coli: our RNA-seq (raw counts / CFU spike-in size factors, mean of 3 replicates), genes summed per 10-kb bin by start
position; bound bins = anti-H-NS ChIP z >= 2 (96 h). Run with the hicexplorer env python (pyBigWig).
Output: rok_vs_hns_silencing.csv, rok_cluster_local.csv (log2 local effect +/-50 kb around each Rok cluster)
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
from io import StringIO
import numpy as np
import pandas as pd
import pyBigWig

ROOT = HNS_ROOT
B = ROOT + "Bactoloops/"
sys.path.insert(0, B)
import rotation_test as RT                     # noqa: E402

rng = np.random.default_rng(0)
BIN = 10_000


def stats(d, bound, label):
    n = len(d)
    off = lambda bs: np.nanmean(d[bs]) - np.nanmean(d)
    obs = off(bound)
    boot = [off(rng.choice(bound, len(bound))) for _ in range(2000)]
    null = np.array([off((bound + s) % n) for s in RT.shifts(n)])
    return dict(set=label, bound_bins=len(bound), offset=obs, lo=np.percentile(boot, 2.5), hi=np.percentile(boot, 97.5),
                p=RT.p_value(obs, null), boot=np.array(boot), null=null)


def change(dA, dB, bound):
    n = len(dA)
    off = lambda d, bs: np.nanmean(d[bs]) - np.nanmean(d)
    obs = off(dB, bound) - off(dA, bound)
    boot = []
    for _ in range(2000):
        bs = rng.choice(bound, len(bound)); boot.append(off(dB, bs) - off(dA, bs))
    null = np.array([off(dB, (bound + s) % n) - off(dA, (bound + s) % n) for s in RT.shifts(n)])
    return obs, np.percentile(boot, 2.5), np.percentile(boot, 97.5), RT.p_value(obs, null)


rows = []
# ---------------- B. subtilis / Rok ----------------
G = B + "Dugar et al data/GSE144475_RAW/"
L = 4_215_606
nb = -(-L // BIN)


def binsum(f):
    bw = pyBigWig.open(G + f)
    v = np.nan_to_num(np.array(bw.values("NC_000964.3", 0, L), float))
    return np.add.reduceat(v, np.arange(0, L, BIN))[:nb] / 1000.0


E = {k: binsum(f) for k, f in [("WT_ex", "GSM4288928_WT_ex_RNAseq.bigwig"), ("mut_ex", "GSM4288929_Rok_ex_RNAseq.bigwig"),
                                ("WT_st", "GSM4288930_WT_st_RNAseq.bigwig"), ("mut_st", "GSM4288931_Rok_st_RNAseq.bigwig")]}
rok = pd.read_csv(B + "Dugar et al data/rok_zscore_10kb.csv")
zr = np.full(nb, np.nan); zr[(rok.start_bp // BIN).astype(int)] = rok.zscore.values
RC = {"RC1": 54, "RC2": 206, "RC3": 226, "RC4": 325, "RC5": 346, "RC6": 366, "RC7": 412, "RC8": 417}
sets_b = {"Rok ChIP z>=2": np.where(zr >= 2)[0],
          "Dugar Rok clusters RC1-8 (+/-1 bin)": np.unique([(c + k) % nb for c in RC.values() for k in (-1, 0, 1)])}
dB = {ph: np.log2((E[f"mut_{ph}"] + 1) / (E[f"WT_{ph}"] + 1)) for ph in ("ex", "st")}
for lab, bs in sets_b.items():
    for ph in ("ex", "st"):
        r = stats(dB[ph], bs, lab); r.update(organism="B. subtilis (Rok)", phase={"ex": "exponential", "st": "stationary 22 h"}[ph])
        rows.append(r)
    c = change(dB["ex"], dB["st"], bs)
    rows.append(dict(organism="B. subtilis (Rok)", set=lab, phase="change st - ex", bound_bins=len(bs), offset=c[0],
                     lo=c[1], hi=c[2], p=c[3]))

# local effect around each Rok cluster (the measure Dugar et al. reported): +/-50 kb, relative to the genome mean
loc = [dict(cluster=k, bin=c, **{ph: float(np.mean(dB[ph][[(c + j) % nb for j in range(-5, 6)]]) - np.mean(dB[ph]))
                                for ph in ("ex", "st")}) for k, c in RC.items()]
pd.DataFrame(loc).to_csv(B + "Dugar et al data/rok_cluster_local.csv", index=False)

# ---------------- E. coli / H-NS ----------------
rd = lambda f: pd.read_csv(StringIO(open(f, newline="").read().replace("\r\n", "\n").replace("\r", "\n")), sep="\t", index_col=0)
RCt = rd(ROOT + "Spike-In/Raw_Counts_Data.txt")
SF = rd(ROOT + "Spike-In/BestSizeFactors.txt").SizeFactor
Nrm = RCt / SF[RCt.columns]
pos = pd.read_csv(ROOT + "Spike-In/genome_positions/gene_positions.txt", sep="\t").set_index("GeneID").Gene_start_lab
Nrm = Nrm.loc[Nrm.index.intersection(pos.index).drop("hns", errors="ignore")]
ne = 461
gb = ((pos.loc[Nrm.index] - 1) // BIN).astype(int).values % ne
binE = lambda cols: np.bincount(gb, Nrm[cols].mean(1).values, minlength=ne)
EE = {"WT_ex": binE([f"WT_Exp_{i}" for i in (1, 2, 3)]), "mut_ex": binE([f"dHNS_Exp_{i}" for i in (1, 2, 3)]),
      "WT_st": binE([f"WT_96h_{i}" for i in (1, 2, 3)]), "mut_st": binE([f"dHNS_96h_{i}" for i in (1, 2, 3)])}
ch = pd.read_csv(B + "ChIP-seq_data/zscore_wt_96h_10kb.csv")
ze = np.full(ne, np.nan); ze[(ch.start_bp // BIN).astype(int) % ne] = ch.zscore.values
# scale each condition to a common total so the +1 pseudocount means the same in every sample
tot = np.mean([v.sum() for v in EE.values()])
EE = {k: v / v.sum() * tot * (E["WT_ex"].sum() / E["WT_ex"].sum()) for k, v in EE.items()}
scale = np.mean([v.mean() for v in E.values()]) / np.mean([v.mean() for v in EE.values()])
EE = {k: v * scale for k, v in EE.items()}                    # same mean bin value as the B. subtilis tracks
dE = {ph: np.log2((EE[f"mut_{ph}"] + 1) / (EE[f"WT_{ph}"] + 1)) for ph in ("ex", "st")}
bs = np.where(ze >= 2)[0]
for ph in ("ex", "st"):
    r = stats(dE[ph], bs, "H-NS ChIP z>=2"); r.update(organism="E. coli (H-NS)", phase={"ex": "exponential", "st": "deep stationary 96 h"}[ph])
    rows.append(r)
c = change(dE["ex"], dE["st"], bs)
rows.append(dict(organism="E. coli (H-NS)", set="H-NS ChIP z>=2", phase="change st - ex", bound_bins=len(bs), offset=c[0],
                 lo=c[1], hi=c[2], p=c[3]))

T = pd.DataFrame(rows)[["organism", "set", "phase", "bound_bins", "offset", "lo", "hi", "p"]]
T["fold"] = 2 ** T.offset
T.to_csv(B + "Dugar et al data/rok_vs_hns_silencing.csv", index=False)
pd.set_option("display.width", 220)
print(T.round(3).to_string(index=False))
