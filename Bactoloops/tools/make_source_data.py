"""Source Data workbook: the existing sheets of Way Source Data.xlsx, renamed to the final supplementary-figure numbers,
plus one sheet per panel (or panel group) of the new and rebuilt figures, in the same simple format (a header row, then
columns of values).

The plotted values are taken from the figure scripts themselves: each script is run with saving switched off and the
arrays it draws are read from its namespace (tools/_ns_probe.run), so the numbers here are exactly those in the figures.
Images of contact maps (Fig. 2a-c, Supplementary Figs. 2, 4a-f, 9a-e, 14a, upper parts) are not repeated: the Hi-C data are in GEO;
their integrated loop density tracks are included.

Output: Paper Update/Way Source Data_EA.xlsx (the original is not modified)
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
import sys
import numpy as np
import pandas as pd
import openpyxl

B = HNS_ROOT + "Bactoloops/"
sys.path.insert(0, B + "tools")
from _ns_probe import run   # noqa: E402

SRC = B + "Paper Update/Way Source Data.xlsx"
OUT = B + "Paper Update/Way Source Data_EA.xlsx"
N, BIN, ORI = 461, 10000, 367
RENAME = {"Fig 3b": "Fig 4b", "Fig 3c": "Fig 4c", "Fig 5b-WT": "Fig 6b-WT", "Fig 5b-dHNS": "Fig 6b-dHNS", "Fig 5c": "Fig 6c",
          "Supplementary Fig 8a": "Supplementary Fig 10a", "Supplementary Fig 8b": "Supplementary Fig 10b",
          "Supplementary Fig 8c": "Supplementary Fig 10c", "Supplementary Fig 9b-c": "Supplementary Fig 11b-c",
          "Supplementary Fig 11a - WT": "Supplementary Fig 14a - WT", "Supplementary Fig 11a - delHNS": "Supplementary Fig 14a - delHNS",
          "Supplementary Fig 11b": "Supplementary Fig 14b"}
SHEETS = {}                                       # name -> DataFrame, in figure order


def add(name, df):
    assert len(name) <= 31, name
    SHEETS[name] = df.reset_index(drop=True)


def clean(s):
    return s.replace("$\\it{", "").replace("}$", "").replace("Δ", "Δ")


gkb = lambda gb: (np.asarray(gb) % N) * BIN // 1000            # genome bin -> bin start (kb, lab genome)
tiff_kb = gkb(ORI + np.arange(N))                              # TIFF/plot order (map starts at oriC side, 3670 kb)


def loops_table(L, extra=None):
    g1, g2 = L.genome_bin1.round().astype(int).values % N, L.genome_bin2.round().astype(int).values % N
    a, b = np.minimum(g1, g2), np.maximum(g1, g2)
    d = pd.DataFrame({"anchor1_kb": gkb(a), "anchor2_kb": gkb(b), "confidence": L.confidence.values,
                      "integrated_excess": L.integrated_excess.round(4).values})
    for k, v in (extra or {}).items():
        d[k] = v
    return d


# ---------------------------------------------------------------- Figure 2
ns = run(B + "updated figures/figure2.py")
lab2 = {"wt_Exp": "WT exp", "wt_96h_repavg": "WT 96 h", "HNSm_96h": "dhns 96 h"}
add("Fig 2a-c", pd.DataFrame({"genome_position_kb": tiff_kb,
                              **{f"ILD {lab2[k]} (median-subtracted)": np.round(v, 3) for k, v in ns["ild_med"].items()}}))
add("Fig 2d,g,h", pd.DataFrame({"genome_position_kb": gkb(np.arange(N)),
                                "ILD z WT exp": np.round(ns["ild_z_genome"]["wt_Exp"], 4),
                                "ILD z WT 96 h": np.round(ns["ild_z_genome"]["wt_96h_repavg"], 4),
                                "anti-H-NS ChIP z WT 96 h": np.round(ns["chip"]["hns96"], 4),
                                "anti-mCherry ChIP z hns-PAmCherry 96 h": np.round(ns["chip"]["mch96"], 4),
                                "anti-mCherry ChIP z hns-PAmCherry exp": np.round(ns["chip"]["mchexp"], 4)}))
for k, panel in (("wt_Exp", "Fig 2e"), ("wt_96h_repavg", "Fig 2f")):
    L = ns["loops"][k]
    add(panel, loops_table(L, {"baseline-corrected integrated excess (%)": L.excess_corr.round(4).values}))

# ---------------------------------------------------------------- Figure 3 (hairpins; Supplementary Fig. 9i-m before 2026-10-06)
ns = run(B + "updated figures/figure3_hairpins.py")
add("Fig 3a", pd.DataFrame({"genome_position_kb": gkb(np.arange(N)), **{lab: np.round(v, 4) for lab, v, _ in ns["TRK"]}}))
t = np.round(np.linspace(-3, 3, 41), 2); rows = []
for key, lab, _ in ns["PILES"]:
    with np.errstate(divide="ignore"):
        M = ns["_blur"](ns["_fourfold"](np.log2(ns["PZ"][key])))
    M = M - np.nanmean(M)
    rows.append(pd.DataFrame({"condition": clean(lab), "position_1": np.repeat(t, 41), "position_2": np.tile(t, 41), "log2_OE_rel_to_mean": np.round(M.ravel(), 4)}))
add("Fig 3b", pd.concat(rows))
add("Fig 3c", ns["SC"].reset_index().assign(stem_flank_ratio=lambda x: np.round(2 ** x.median_log2, 4),
                                                          ci_lo_ratio=lambda x: np.round(2 ** x.lo, 4), ci_hi_ratio=lambda x: np.round(2 ** x.hi, 4))
    [["test", "stem_flank_ratio", "ci_lo_ratio", "ci_hi_ratio", "p", "n"]])
SS = ns["SS"].reset_index()
for col, panel in (("chid", "Fig 3d"), ("opcid", "Fig 3e")):
    add(panel, pd.DataFrame({"condition": SS.condition, "libraries": SS.libraries, "within_domain_enrichment": np.round(2 ** SS[f"{col}_log2"], 4),
                             "ci_lo": np.round(2 ** SS[f"{col}_lo"], 4), "ci_hi": np.round(2 ** SS[f"{col}_hi"], 4), "p": SS[f"{col}_p"], "n": SS[f"{col}_n"]}))


# ---------------------------------------------------------------- Figure 5 (RNA-seq; Figure 4 before 2026-10-06)
ns = run(B + "updated figures/figure4.py")
vol = []
for tp in ("Exp", "96h"):
    d = pd.read_csv(ns["DE"].format(tp), sep="\t")
    vol.append(d[["GeneID", "log2FoldChange", "padj"]].rename(columns={"log2FoldChange": f"log2FC WT/dhns {tp}", "padj": f"padj {tp}"}))
add("Fig 5a-b", vol[0].merge(vol[1], on="GeneID", how="outer"))
G = ns["G"].rename(columns={"log2FoldChange": "log2FC WT/dhns", "Position_from_oriC_axis": "position_kb (axis from oriC)"})
add("Fig 5c", G)
add("Fig 5d-f", pd.DataFrame({"genome_position_kb": gkb(np.arange(N)), "ILD z WT 96 h": np.round(ns["ild_z"], 4),
                              "anti-H-NS ChIP z 96 h": np.round(ns["chip_z"], 4), "H-NS silencing z 96 h": np.round(ns["rna_z"], 4)}))
P = ns["P"]
add("Fig 5g", pd.DataFrame({"anchor1_kb": gkb(P.min(1)), "anchor2_kb": gkb(P.max(1)),
                            "mean anti-H-NS ChIP z of anchors": np.round(ns["chip_z"][P].mean(1), 4),
                            "mean H-NS silencing z of anchors": np.round(ns["rna_z"][P].mean(1), 4)}))
fs = pd.read_csv(B + "updated figures/fig4_summary_CFU.csv").iloc[0]
add("Fig 5h", pd.DataFrame([dict(phase=p, silencing_offset_log2=round(fs[f"h_{k}"], 4), ci_lo=round(fs[f"h_{k}_lo"], 4),
                                 ci_hi=round(fs[f"h_{k}_hi"], 4), P_rotation=round(fs[f"h_{k}_p"], 4))
                            for p, k in (("exp", "Exp"), ("96 h", "96h"), ("change 96 h - exp", "diff"))]))

# ---------------------------------------------------------------- Supplementary Fig. 2f (distance-matched confidence)
import matplotlib.axes
_cap = {}
_hist, _ann = matplotlib.axes.Axes.hist, matplotlib.axes.Axes.annotate
def _h(self, x, *a, **k):
    _cap.setdefault("hist", np.asarray(x, float)); return _hist(self, x, *a, **k)
def _a(self, text, *a, **k):
    if text == "loop": _cap["loop"] = k.get("xy", a[0] if a else None)[0]
    return _ann(self, text, *a, **k)
matplotlib.axes.Axes.hist, matplotlib.axes.Axes.annotate = _h, _a
try:
    run(B + "updated figures/figureS2_method.py", ["A"])
finally:
    matplotlib.axes.Axes.hist, matplotlib.axes.Axes.annotate = _hist, _ann
nv = _cap["hist"]
add("Supplementary Fig 2f", pd.DataFrame({"pixel": np.arange(1, len(nv) + 1),
                                          "log2_OE_local_bg_removed (all pixels at 290 kb)": np.round(nv, 4)}))
add("Supplementary Fig 2f loop", pd.DataFrame([dict(loop_log2_OE=round(float(_cap["loop"]), 4), percentile=round(100 * float(np.mean(nv < _cap["loop"])), 1),
                                                    pct95=round(float(np.percentile(nv, 95)), 4), pct99=round(float(np.percentile(nv, 99)), 4))]))

# ---------------------------------------------------------------- Supplementary Fig. 5 (ChIP-seq, 1-kb ChIP/input; tools/s5_tracks.py)
s5 = pd.read_csv(B + "ChIP-seq_data/s5_tracks_1kb.csv").rename(columns={"bin_start_kb": "genome_position_kb"})
add("Supplementary Fig 5a-d", s5)
REG = [(-60, 60), (2380, 2510), (690, 805)]                    # the three regions of e and f (kb from the start; -60 = 4547)
def region(cols):
    out = []
    for lo, hi in REG:
        kb = np.arange(lo, hi) % len(s5)
        out.append(pd.DataFrame({"region": f"{lo}-{hi} kb", "genome_position_kb": np.arange(lo, hi), **{c: s5[c].values[kb] for c in cols}}))
    return pd.concat(out)
add("Supplementary Fig 5e", region(["anti-H-NS, WT, 96 h", "anti-mCherry, hns-PAmCherry, 96 h"]))
add("Supplementary Fig 5f", region(["anti-mCherry, hns-PAmCherry, exp", "anti-mCherry, hns-PAmCherry, 96 h"]))

# ---------------------------------------------------------------- Supplementary Fig. 3 (loop callers)
ns = run(B + "Comparisons/FigS_method_comparison.py")
T, METHODS, NAME = ns["T"], ns["METHODS"], ns["NAME"]
SETS = [d for d, *_ in ns["SETS"]]
add("Supplementary Fig 3a", T[T.method.isin(METHODS) & T.dataset.isin(SETS)].assign(method=lambda x: x.method.map(NAME))
    .pivot(index="method", columns="dataset", values="n_loops").reindex([NAME[m] for m in METHODS])[SETS].reset_index())
nsC = ns["ns"]
for d, panel in (("wt96", "Supplementary Fig 3b"), ("bsub", "Supplementary Fig 3c")):
    n = 461 if d == "wt96" else 422
    z = nsC["chip_track"](d, n); rows = []
    for m in METHODS:
        Pm, _ = nsC["load"](m, d)
        for a, b in Pm:
            rows.append(dict(method=NAME[m], anchor1_kb=int(a) * 10, anchor2_kb=int(b) * 10, mean_anchor_ChIP_z=round(float(np.nanmean(z[[a, b]])), 4)))
    add(panel, pd.DataFrame(rows))
D3 = T[T.method.isin(METHODS) & T.dataset.isin(["wt96", "stpA96", "fis96", "ihf96"])]
add("Supplementary Fig 3d", D3.assign(method=D3.method.map(NAME))[["method", "dataset", "n_loops", "chip_mean_z", "chip_mean_z_ci_lo", "chip_mean_z_ci_hi"]].round(4))
rows = []
for m in METHODS:
    (A_, sa), (B_, sb) = nsC["load"](m, "wt96r1"), nsC["load"](m, "wt96r2")
    kmax = min(len(A_), len(B_))
    if kmax < 3:
        continue
    A_, B_ = A_[np.argsort(-sa)], B_[np.argsort(-sb)]
    for k in range(3, kmax + 1):
        rows.append(dict(method=NAME[m], top_k=k, fraction_shared=round((ns["_found"](A_[:k], B_[:k]) + ns["_found"](B_[:k], A_[:k])) / 2, 4)))
add("Supplementary Fig 3e", pd.DataFrame(rows))
R = pd.read_csv(B + "Comparisons/rotation_test/rotation_detail.csv")
R = R[~R.missing]
for d, panel in (("wt96", "Supplementary Fig 3f"), ("bsub", "Supplementary Fig 3g")):
    r = R[R.dataset == d]
    add(panel, pd.DataFrame({"method": r.tool, "rotation_bins": r.k, "loops_unrotated": r.n0.astype(int), "loops_rotated": r.nk.astype(int),
                             "kept": r.kept.astype(int), "pct_kept": (100 * r.kept / r.n0).round(2), "lost": r.lost.astype(int), "gained": r.gained.astype(int)}))
add("Supplementary Fig 3h", SHEETS["Supplementary Fig 3b"].copy())

# ---------------------------------------------------------------- Supplementary Fig. 4 (loops in every 96 h map)
ns = run(B + "updated figures/figureS3_loop_conditions.py")
add("Supplementary Fig 4a-f", pd.DataFrame({"genome_position_kb": tiff_kb,
                                            **{f"ILD {run_} (median-subtracted)": np.round(v, 3) for run_, v in ns["ildm"].items()}}))
add("Supplementary Fig 4g", ns["S"].rename(columns={"pct_WT_recovered": "pct_of_WT_96h_loops_also_called"}))
labs = [clean(l).replace("\n", " ") for l in ns["labs"]]
add("Supplementary Fig 4h", pd.DataFrame(np.round(ns["R"], 4), columns=labs).assign(map=labs)[["map"] + labs])

# ---------------------------------------------------------------- Supplementary Fig. 7 (H-NS at anchors, every map)
ns = run(B + "updated figures/figureS_hns_comparison.py")
rows, zrows = [], []
for lab, run_, *_ in ns["CASES"]:
    L = ns["L"][run_]
    if not len(L):
        continue
    t = loops_table(L, {"baseline-corrected integrated excess (%)": L.excess_corr.round(4).values})
    t.insert(0, "map", clean(lab)); rows.append(t)
    Pm = np.c_[L.genome_bin1.round().astype(int), L.genome_bin2.round().astype(int)] % N
    zrows.append(pd.DataFrame({"map": clean(lab), "anchor1_kb": gkb(Pm.min(1)), "anchor2_kb": gkb(Pm.max(1)),
                               "mean_anchor_H-NS_ChIP_z": np.round(np.nanmean(ns["z_gen"][Pm], axis=1), 4)}))
add("Supplementary Fig 7a-g", pd.concat(rows))
S7 = ns["S"].assign(case=lambda x: x.case.map(clean))
add("Supplementary Fig 7h", S7[["case", "loops", "ild_chip_r", "ild_chip_p"]].round(4))
add("Supplementary Fig 7i", pd.concat(zrows))
add("Supplementary Fig 7j", S7[["case", "peak_pair_excess", "excess_lo", "excess_hi", "excess_rot95", "excess_p"]].round(4))

# ---------------------------------------------------------------- Supplementary Fig. 8 (H-NS binding at anchors)
z = np.full(N, np.nan); Zc = pd.read_csv(B + "ChIP-seq_data/zscore_wt_96h_10kb.csv"); z[(Zc.start_bp // BIN).values % N] = Zc.zscore.values
ild = np.empty(N); ild[(ORI + np.arange(N)) % N] = pd.read_csv(B + "FIJI_output/integrated_loop_density_wt_96h_repavg.csv").ILD_au.values
ildz = (ild - ild.mean()) / ild.std()
add("Supplementary Fig 8a-b", pd.DataFrame({"genome_position_kb": gkb(np.arange(N)), "ILD z WT 96 h": np.round(ildz, 4),
                                            "anti-H-NS ChIP z 96 h": np.round(z, 4), "ILD z > 1": ildz > 1, "ChIP z > 1": np.nan_to_num(z) > 1}))
T8 = pd.read_csv(B + "anchor_clusters_1kb_z2.csv")
add("Supplementary Fig 8c", T8[["kind", "region_kb", "max_z10"]].rename(columns={"max_z10": "strongest_10kb_H-NS_ChIP_z"}))

# ---------------------------------------------------------------- Supplementary Fig. 9 (time course, CHINs)
ns = run(B + "updated figures/figureS7_timecourse.py")
lab9 = {run_: clean(lab) for lab, run_, _ in ns["MAPS"]}
add("Supplementary Fig 9a-e", pd.DataFrame({"genome_position_kb": tiff_kb,
                                            **{f"ILD {lab9[k]} (median-subtracted)": np.round(v - np.median(v), 3) for k, v in ns["ild"].items()}}))
T9 = ns["T"].reset_index()
rows = []
for lab, run_, _ in ns["MAPS"]:
    try:
        L = pd.read_csv(B + f"FIJI_output/loop_calls_{run_}.csv")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        continue
    Pm = np.c_[L.genome_bin1.round().astype(int), L.genome_bin2.round().astype(int)] % N
    rows.append(pd.DataFrame({"map": clean(lab), "anchor1_kb": gkb(Pm.min(1)), "anchor2_kb": gkb(Pm.max(1)),
                              "mean_anchor_H-NS_ChIP_z": np.round(np.nanmean(z[Pm], axis=1), 4)}))
add("Supplementary Fig 9f", pd.concat(rows))
add("Supplementary Fig 9f summary", T9[T9.columns[:5]].round(4))
add("Supplementary Fig 9g", pd.DataFrame(ns["g_rows"], columns=["map", "mean_excess_pct", "ci_lo", "ci_hi", "rotation_mean", "rotation_2.5pct", "rotation_97.5pct", "P_rotation"])
    .assign(map=lambda x: x.map.map(clean)).round(4))
add("Supplementary Fig 9h", T9[[T9.columns[0], "ILDz_at_peaks", "ILDz_elsewhere"]].round(4))
# ---------------------------------------------------------------- Supplementary Fig. 12 (RNA-seq)
ns = run(B + "updated figures/figureS10_rnaseq.py")
ma = []
for tp in ("Exp", "96h"):
    d = ns["D"][tp]
    ma.append(d[["GeneID", "baseMean", "log2FoldChange", "padj"]].rename(columns={"baseMean": f"baseMean {tp}", "log2FoldChange": f"log2FC WT/dhns {tp}", "padj": f"padj {tp}"}))
add("Supplementary Fig 12a-b", ma[0].merge(ma[1], on="GeneID", how="outer"))
add("Supplementary Fig 12c", pd.DataFrame({"genome_position_kb": gkb(np.arange(N)), "mean log2FC exp": np.round(ns["be"], 4),
                                           "mean log2FC 96 h": np.round(ns["b9"], 4), "anti-H-NS ChIP z 96 h": np.round(ns["zbin"], 4)}))
add("Supplementary Fig 12d", pd.DataFrame({"genome_position_kb": gkb(np.arange(N)), **{lab: np.round(v, 4) for lab, v, _ in ns["T"]}}))
F = pd.concat([pd.read_csv(B + f"updated figures/fig4_summary_{n}.csv") for n in ("CFU", "OD", "Count")])
add("Supplementary Fig 13a-d", F[["norm", "r_e", "p_e", "r_f", "p_f", "g_rna_mean", "g_rna_lo", "g_rna_hi", "g_rna_rot95", "g_rna_p",
                                  "up_Exp", "down_Exp", "up_96h", "down_96h", "h_Exp", "h_Exp_lo", "h_Exp_hi", "h_96h", "h_96h_lo", "h_96h_hi",
                                  "h_diff", "h_diff_lo", "h_diff_hi", "h_diff_p"]].round(4))

# ---------------------------------------------------------------- Supplementary Fig. 15 (E. coli vs B. subtilis; was 14)
ns = run(B + "Supplemental_figure_Ecoli_vs_Bsub/fig_supplement_ecoli_vs_bsub.py")
rows = []
for sp in ("E. coli", "B. subtilis"):
    L = ns["full"][sp][1]; n = 461 if sp == "E. coli" else 422
    Pm = np.c_[L.genome_bin1.round().astype(int), L.genome_bin2.round().astype(int)] % n
    ex = L["excess_corr"].values if "excess_corr" in L.columns else L.integrated_excess.values
    rows.append(pd.DataFrame({"organism": sp, "anchor1_kb": Pm.min(1) * 10, "anchor2_kb": Pm.max(1) * 10,
                              "baseline-corrected integrated excess (%)": np.round(ex, 4),
                              "relative to dataset mean": np.round(ex / np.mean(ex), 4)}))
add("Supplementary Fig 15b-c", pd.concat(rows))
rows = []
for sp in ("E. coli", "B. subtilis"):
    Pp = ns["full"][sp][2]
    rows.append(Pp.assign(organism=sp, anchor1_kb=Pp.i * 10, anchor2_kb=Pp.j * 10, distance_kb=Pp.dist * 10)
                [["organism", "anchor1_kb", "anchor2_kb", "distance_kb", "called", "oe", "null_med", "p95", "pct"]]
                .rename(columns={"oe": "corrected_log2_OE", "null_med": "distance_matched_median", "p95": "distance_matched_95th", "pct": "percentile"}).round(4))
add("Supplementary Fig 15d", pd.concat(rows))
rows = []
for sp in ("E. coli", "B. subtilis"):
    per_loop, per_anchor = ns["eng"][sp]
    rows += [dict(organism=sp, unit="loop", value_pct=round(float(v), 4)) for v in per_loop]
    rows += [dict(organism=sp, unit="anchor", value_pct=round(float(v), 4)) for v in per_anchor]
add("Supplementary Fig 15e", pd.DataFrame(rows))
add("Supplementary Fig 15e summary", pd.DataFrame([dict(organism=sp, unit=u, relative_to_E_coli=round(v[0], 4), ci_lo=round(v[1], 4), ci_hi=round(v[2], 4))
                                                   for (sp, u), v in ns["rel_d"].items()]))
rows = []
for sp, chipf, n in (("E. coli", B + "ChIP-seq_data/zscore_wt_96h_10kb.csv", 461), ("B. subtilis", B + "Dugar et al data/rok_zscore_10kb.csv", 422)):
    run_ = "wt_96h_repavg" if sp == "E. coli" else "bsub_WT_st_merged"
    il = pd.read_csv(B + f"FIJI_output/integrated_loop_density_{run_}.csv")
    gb = (il.genome_position_bp.values // BIN).astype(int) % n
    zz = np.full(n, np.nan); c = pd.read_csv(chipf); zz[(c.start_bp // BIN).values % n] = c.zscore.values
    rows.append(pd.DataFrame({"organism": sp, "genome_position_kb": gb * 10, "ILD": il.ILD_au.values, "ChIP_z": np.round(zz[gb], 4)}))
add("Supplementary Fig 15f", pd.concat(rows))
add("Supplementary Fig 15g", ns["SIL"].round(4))
add("Supplementary Fig 15h", pd.read_csv(B + "Dugar et al data/rok_cluster_local.csv").rename(
    columns={"bin": "genome_bin", "ex": "mean log2(drok/WT) +-50 kb, exp, minus genome mean", "st": "mean log2(drok/WT) +-50 kb, 22 h, minus genome mean"}).round(4))

# ---------------------------------------------------------------- Figure 7 (top: loops, H-NS binding, silencing; figure7.py)
ns = run(B + "updated figures/figure7.py")
rows, lrows = [], []
for P in ns["PHASES"]:
    ph = "exp" if P["title"].startswith("Exp") else "96 h"
    rows.append(pd.DataFrame({"growth_phase": ph, "genome_position_kb": gkb(np.arange(N)),
                              "H-NS ChIP z (anti-mCherry), inner ring": np.round(P["chip"], 4),
                              "H-NS silencing, log2(dhns/WT) relative to genome mean, outer ring": np.round(P["sil"], 4)}))
    L = P["L"]; g1, g2 = L.genome_bin1.round().astype(int).values % N, L.genome_bin2.round().astype(int).values % N
    lrows.append(pd.DataFrame({"growth_phase": ph, "anchor1_kb": gkb(np.minimum(g1, g2)), "anchor2_kb": gkb(np.maximum(g1, g2)),
                               "baseline-corrected integrated excess (%)": np.round(L["excess"].values, 4)}))
add("Fig 7 rings", pd.concat(rows))
add("Fig 7 loops", pd.concat(lrows))

# ---------------------------------------------------------------- write
wb = openpyxl.load_workbook(SRC)
for old, new in RENAME.items():
    wb[old].title = new
for name, df in SHEETS.items():
    ws = wb.create_sheet(name)
    ws.append([str(c) for c in df.columns])
    for row in df.itertuples(index=False):
        ws.append([None if (isinstance(v, float) and np.isnan(v)) else (v.item() if hasattr(v, "item") else v) for v in row])
order = ["Fig 2", "Fig 3", "Fig 4", "Fig 5", "Fig 6", "Fig 7"] + [f"Supplementary Fig {i}" for i in range(1, 16)]
def key(t):
    for i, p in enumerate(order):
        if t.startswith(p) and (not p[-1].isdigit() or not t[len(p):len(p) + 1].isdigit()):
            return (i, t)
    return (99, t)
wb._sheets.sort(key=lambda ws: key(ws.title))
note = wb.create_sheet("README", 0)
for r in [["Source Data for Way et al."],
          ["One sheet per figure panel (or panel group): a header row, then the plotted values."],
          ["Genome positions are 10-kb bin starts (kb) on the lab genome (W3110 NC_007779.1 minus 1,398,953-1,438,196; 4,607,088 bp); "
           "B. subtilis positions are on 168 (NC_000964.3)."],
          ["Contact-map images (Fig. 1, Fig. 2a-c maps, Supplementary Figs. 1c, 2, 4a-f, 6c, 9a-e, 15a) are not repeated here: "
           "the Hi-C data are deposited in GEO; the integrated loop density (ILD) tracks drawn under the maps are included."],
          ["Loop calls and anchors are also in Supplementary Dataset 1."],
          ["Supplementary Fig 5: 1-kb ChIP/input (each track divided by its total reads). The WT exponential-phase anti-H-NS track "
           "(Supplementary Fig 5a, left) was not available in the analysis folder and is not included."],
          ["Supplementary Fig 4h: the correlation matrix sheet gives every pair of maps (the figure prints values only for WT exp and dhns 96 h)."],
          ["Fig 3f and the lower part of Fig 7 are schematics (no data)."],
          ["Renumbered 2026-10-06: Supplementary Fig 10a-c, 11b-c, and 13a-b were Supplementary Fig 8, 9, and 11; Figs 4-6 were Figs 3-5 (Fig 3, H-NS hairpins, was Supplementary Fig 9i-m)."],
          ["Renumbered 2026-10-12: Supplementary Fig 12 was split into 12 (a-d) and 13 (the former 12e-h, now 13a-d); the former Supplementary Figs 13 and 14 are now 14 and 15."]]:
    note.append(r)
wb.save(OUT)
print(len(SHEETS), "new sheets;", len(wb.sheetnames), "sheets in all ->", OUT)
print(", ".join(wb.sheetnames))
