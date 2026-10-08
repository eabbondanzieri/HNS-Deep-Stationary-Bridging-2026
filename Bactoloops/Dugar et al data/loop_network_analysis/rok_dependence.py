"""Which B. subtilis WT loops depend on Rok? Uses the delta-rok maps as a control.

For every loop called in the merged WT map (FIJI_output/loop_calls_bsub_WT_st_merged.csv):
- log2(obs/exp), 3x3, local-background corrected (50-bin rolling ball), in merged WT and merged
  delta-rok -- the same map the macro's log2_obs_exp is read from;
- Rok dependence = WT - delta-rok. Its significance is judged against the same
  difference at every pixel of the same genomic distance (distance-matched null);
- the same difference within each sequencing batch (WT run 1 - rok run 1, WT run 2 -
  rok run 2), since the two runs differ globally more than the genotypes do;
- whether the loop lies on the SMC arm-alignment stripe near oriC (i + j close to 0
  mod N and both anchors within 600 kb of oriC), and the anchors' Rok ChIP z-score.
"""
import sys
import numpy as np
import pandas as pd
import loop_network_stats as S

sys.path.insert(0, S.ROOT + "iced_matrix_conversion/Dugar_HiC")
from merge_and_ice import raw, ice, obs_exp_log2, RUNS  # noqa: E402

N = 422
LOOPS = S.ROOT + "FIJI_output/loop_calls_bsub_WT_st_merged.csv"
CHIP = S.ROOT + "Dugar et al data/rok_zscore_10kb.csv"
RC = {"RC1": 54, "RC2": 206, "RC3": 226, "RC4": 325, "RC5": 346, "RC6": 366, "RC7": 412, "RC8": 417, "yonX": 220}

maps = {
    "WT": ice(sum(raw(f) for f in RUNS["WT"])),
    "rok": ice(sum(raw(f) for f in RUNS["rok"])),
    "WT1": ice(raw(RUNS["WT"][0])), "WT2": ice(raw(RUNS["WT"][1])),
    "rok1": ice(raw(RUNS["rok"][0])), "rok2": ice(raw(RUNS["rok"][1])),
}
sys.path.insert(0, S.ROOT + "python")
from ij_ops import subtract_background, circular_subtract_background  # noqa: E402   (ImageJ's rolling ball, ported exactly)


def corrected(m, pad=100, radius=50):
    """log2(obs/exp), 3x3-averaged, minus a 50-bin rolling-ball local background -- the macro's
    HiC_map_log2OE (padded periodically so the circular chromosome has no edges)."""
    lo = obs_exp_log2(m).astype(np.float32)
    return circular_subtract_background(lo, radius).astype(float)


oe = {k: corrected(m) for k, m in maps.items()}
delta = oe["WT"] - oe["rok"]
r = np.arange(N)


def dist(a, b):
    d = abs(a - b) % N
    return min(d, N - d)


def delta_pct(a, b):
    d = dist(a, b)
    v = delta[r, (r + d) % N]
    return np.mean(v < delta[a, b])


chip = pd.read_csv(CHIP)
z = np.full(N, np.nan)
z[(chip.start_bp // 10_000).astype(int)] = chip.zscore.values
zmax = lambda b: np.nanmax(z[[(b + k) % N for k in (-1, 0, 1)]])
ori_d = lambda b: min(b, N - b)

L = pd.read_csv(LOOPS)
rows = []
for x in L.itertuples():
    a, b = int(round(x.bin1)), int(round(x.bin2))
    name = next((f"{p}-{q}" for p in RC for q in RC if p != q and abs(RC[p] - a) <= 3 and abs(RC[q] - b) <= 3), "")
    s = (a + b) % N
    rows.append(dict(bin1=a, bin2=b, dist_kb=dist(a, b) * 10, rok_clusters=name,
                     note="ridge region" if x.note.startswith(("anchor", "ridge")) else "main", wt_confidence=x.confidence,
                     oe_WT=oe["WT"][a, b], oe_rok=oe["rok"][a, b], rok_dependence=delta[a, b],
                     dependence_pct=delta_pct(a, b),
                     dep_batch1=oe["WT1"][a, b] - oe["rok1"][a, b], dep_batch2=oe["WT2"][a, b] - oe["rok2"][a, b],
                     arm_stripe=min(s, N - s) <= 25 and max(ori_d(a), ori_d(b)) <= 60,
                     anchor_rok_z=min(zmax(a), zmax(b))))
T = pd.DataFrame(rows)
T["rok_dependent"] = (T.dependence_pct >= 0.95) & (T.dep_batch1 > 0) & (T.dep_batch2 > 0)
T.sort_values("rok_dependence", ascending=False).round(3).to_csv("rok_dependence_loops.csv", index=False)
# check: the merged WT map matches the corrected map the figure scripts use (Fiji's rolling ball)
sys.path.insert(0, S.ROOT + "Supplemental_figure_Ecoli_vs_Bsub")
import localbg  # noqa: E402
_ref = localbg.corrected_maps({"bsub_wt_st_merged": S.ROOT + "iced_matrix_conversion/Dugar_HiC/Dugar_WT_st_merged_iced.tif"})["bsub_wt_st_merged"]["bs"]
print(f"corrected WT map vs the figure's (Fiji) corrected map: max |diff| {np.abs(oe['WT'] - _ref).max():.2e}")

pd.set_option("display.width", 220)
print(T.sort_values("rok_dependence", ascending=False)[
    ["bin1", "bin2", "dist_kb", "rok_clusters", "note", "wt_confidence", "oe_WT", "oe_rok", "rok_dependence",
     "dependence_pct", "dep_batch1", "dep_batch2", "arm_stripe", "anchor_rok_z", "rok_dependent"]].round(2).to_string(index=False))
print("\nRok-dependent (top 5% of distance-matched WT-rok difference AND positive in both batches):")
for grp, g in [("all", T), ("high", T[T.wt_confidence == "high"]), ("medium", T[T.wt_confidence == "medium"]),
               ("low", T[T.wt_confidence == "low"]), ("on arm stripe", T[T.arm_stripe]), ("off arm stripe", T[~T.arm_stripe])]:
    print(f"  {grp:15s} {g.rok_dependent.sum():2d}/{len(g):2d}   median WT {g.oe_WT.median():.2f}, delta-rok {g.oe_rok.median():.2f}")
print("batch agreement (Spearman, dependence in batch 1 vs batch 2):",
      round(pd.Series(T.dep_batch1).corr(pd.Series(T.dep_batch2), method="spearman"), 2))
print("Rok ChIP z at anchors vs dependence (Spearman):",
      round(T.anchor_rok_z.corr(T.rok_dependence, method="spearman"), 2))
