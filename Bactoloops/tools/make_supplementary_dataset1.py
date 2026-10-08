"""Supplementary Dataset 1: Bactoloops loop calls and loop anchors.

Sheets
  README               column definitions and coordinate systems
  WT_96h_loops         the 128 loops of E. coli WT 96 h (average of two biological replicates; Figs. 1, 2)
  WT_96h_anchors       the anchor regions those loops connect (loop ends merged across adjacent 10-kb bins), with H-NS
                       binding, CHIDs, genes and H-NS silencing at 96 h
  Strong_HNS_regions   runs of 10-kb bins with anti-H-NS 96 h ChIP-seq z >= 2, and whether each anchors a loop
  All_loop_calls       every Bactoloops call in every map analysed in the paper (E. coli and B. subtilis)

Inputs: FIJI_output/ (rotation-invariant Bactoloops, 2026-10-06), ChIP-seq_data/zscore_wt_96h_10kb.csv,
Nudler CHIN paper/opcidsRenamed.xlsx (CHIDs, lab genome), Spike-In/genome_positions (genes, lab genome; DESeq2 CFU).
Output: Paper Update/Supplementary Data/Supplementary_Dataset_1_loops_and_anchors.xlsx
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

B = HNS_ROOT + "Bactoloops/"
OUT = B + "Paper Update/Supplementary Data/Supplementary_Dataset_1_loops_and_anchors.xlsx"
N, BIN = 461, 10000
GENOME = 4_607_088                              # EcW3110del: W3110 (NC_007779.1) minus 1,398,953-1,438,196
DEL_START, DEL_LEN = 1_398_953, 39_244
cd = lambda a, b: np.minimum(np.abs(a - b) % N, N - np.abs(a - b) % N)


def w3110(x):
    """lab-genome position (1-based bp) -> W3110 NC_007779.1 position"""
    return x + DEL_LEN if x >= DEL_START else x


def bin_range(b, n=N, genome=GENOME):
    s = b * BIN + 1
    return s, min((b + 1) * BIN, genome)


def calls(name):
    try:
        return pd.read_csv(B + f"FIJI_output/loop_calls_{name}.csv")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def pairs(L, n=N):
    return np.c_[L.genome_bin1.round().astype(int), L.genome_bin2.round().astype(int)] % n


def matched(P, Q, tol=2, n=N):
    d = lambda a, b: np.minimum(np.abs(a - b) % n, n - np.abs(a - b) % n)
    return np.array([bool(len(Q)) and bool((((d(Q[:, 0], a) <= tol) & (d(Q[:, 1], b) <= tol)) |
                                            ((d(Q[:, 0], b) <= tol) & (d(Q[:, 1], a) <= tol))).any()) for a, b in P])


# ---- H-NS ChIP-seq (genome bins), CHIDs, genes ----
z = np.full(N, np.nan)
Z = pd.read_csv(B + "ChIP-seq_data/zscore_wt_96h_10kb.csv"); z[(Z.start_bp // BIN).values % N] = Z.zscore.values
X = pd.read_excel(B + "Nudler CHIN paper/opcidsRenamed.xlsx", 0)
X.loc[X.OPCID_ID == "OPCID_40", "start_new"] = 2_522_632; X.loc[X.OPCID_ID == "CHID_20", "end_new"] = 3_950_261
X = X.dropna(subset=["start_new", "end_new"]); X["type"] = X.OPCID_ID.str.split("_").str[0]
X["s"], X["e"] = X[["start_new", "end_new"]].min(1).astype(int), X[["start_new", "end_new"]].max(1).astype(int)
CHID = X[X.type == "CHID"].copy(); CHID["mid"] = ((CHID.s + CHID.e) // 2 // BIN) % N
G = pd.read_csv(B + "../Spike-In/genome_positions/DESeq2_apeglm_CFU_96h_with_positions.txt", sep="\t")
G["lo"] = G[["Gene_start_lab", "Gene_end_lab"]].min(1); G["hi"] = G[["Gene_start_lab", "Gene_end_lab"]].max(1)
G["silencing"] = -G.log2FoldChange                                # log2(delta-hns / WT), positive = silenced by H-NS

# ---- WT 96 h loops and anchors ----
L = calls("wt_96h_repavg"); P = pairs(L)
ends = np.zeros(N, bool); ends[P.ravel()] = True
regions, cur, s0 = [], [], int(np.where(~ends)[0][0])
for k in range(1, N + 1):
    b = (s0 + k) % N
    if ends[b]:
        cur.append(b)
    elif cur:
        regions.append(cur); cur = []
if cur:
    regions.append(cur)
regions.sort(key=lambda r: r[0])
anchor_of = {b: i + 1 for i, r in enumerate(regions) for b in r}

anc_rows = []
for i, r in enumerate(regions, start=1):
    s, _ = bin_range(r[0]); _, e = bin_range(r[-1])
    loops_here = [j for j, (a, b) in enumerate(P) if anchor_of[a] == i or anchor_of[b] == i]
    partners = sorted({anchor_of[b] if anchor_of[a] == i else anchor_of[a] for a, b in P[loops_here]} - {i})
    near = [(b + k) % N for b in r for k in (-1, 0, 1)]
    chids = CHID[[min(cd(m, np.array(r))) <= 1 for m in CHID.mid]].OPCID_ID.tolist()
    g = G[(G.hi >= s) & (G.lo <= e)].sort_values("lo")
    anc_rows.append(dict(anchor_id=f"A{i:02d}", start_bp=s, end_bp=e, start_bp_W3110=w3110(s), end_bp_W3110=w3110(e),
                         length_kb=(e - s + 1) / 1000, n_loops=len(loops_here), n_partner_anchors=len(partners),
                         partner_anchors=", ".join(f"A{p:02d}" for p in partners),
                         HNS_ChIP_z_max_within_10kb=round(float(np.nanmax(z[near])), 2),
                         HNS_bound=bool(np.nanmax(z[near]) > 1), CHIDs_within_10kb=", ".join(chids),
                         n_genes=len(g), genes=", ".join(g.GeneID.astype(str)),
                         mean_HNS_silencing_96h_log2=round(float(g.silencing.mean()), 2) if len(g) else np.nan))
A = pd.DataFrame(anc_rows)

loop_rows = []
for j, ((a, b), (_, r)) in enumerate(zip(P, L.iterrows()), start=1):
    a, b = sorted((a, b))
    (s1, e1), (s2, e2) = bin_range(a), bin_range(b)
    loop_rows.append(dict(loop_id=f"L{j:03d}", anchor1=f"A{anchor_of[a]:02d}", anchor2=f"A{anchor_of[b]:02d}",
                          end1_start_bp=s1, end1_end_bp=e1, end2_start_bp=s2, end2_end_bp=e2,
                          end1_start_bp_W3110=w3110(s1), end2_start_bp_W3110=w3110(s2),
                          span_kb=int(cd(a, b)) * 10, confidence=r.confidence, oe_percentile=round(r.oe_percentile, 3),
                          log2_obs_exp=round(r.log2_obs_exp, 3), integrated_excess=round(r.integrated_excess, 3),
                          HNS_ChIP_z_end1=round(z[a], 2), HNS_ChIP_z_end2=round(z[b], 2),
                          call="ridge region (donut test)" if str(r.note).startswith("ridge") else "main pass"))
LP = pd.DataFrame(loop_rows).sort_values(["end1_start_bp", "end2_start_bp"]).reset_index(drop=True)
LP["loop_id"] = [f"L{j:03d}" for j in range(1, len(LP) + 1)]

# ---- strongly H-NS-bound regions (runs of z >= 2) ----
hi = np.nan_to_num(z, nan=-9) >= 2
runs, cur, s0 = [], [], int(np.where(~hi)[0][0])
for k in range(1, N + 1):
    b = (s0 + k) % N
    if hi[b]:
        cur.append(b)
    elif cur:
        runs.append(cur); cur = []
if cur:
    runs.append(cur)
runs.sort(key=lambda r: r[0])
abins = np.array(sorted(anchor_of))
S = pd.DataFrame([dict(region_id=f"H{i:02d}", start_bp=bin_range(r[0])[0], end_bp=bin_range(r[-1])[1],
                       start_bp_W3110=w3110(bin_range(r[0])[0]), end_bp_W3110=w3110(bin_range(r[-1])[1]),
                       HNS_ChIP_z_max=round(float(np.nanmax(z[r])), 2),
                       distance_to_nearest_anchor_kb=int(min(cd(abins, b).min() for b in r)) * 10,
                       anchors=", ".join(sorted({f"A{anchor_of[b]:02d}" for b in r if b in anchor_of})))
                  for i, r in enumerate(runs, start=1)])
S["anchors_a_loop"] = S.distance_to_nearest_anchor_kb == 0

# ---- every call in every map ----
MAPS = [("wt_96h_repavg", "E. coli", "WT", "96 h", "average of rep1 + rep2 (Figs. 1, 2)"),
        ("wt_96h_rep1", "E. coli", "WT", "96 h", "replicate 1"), ("wt_96h_rep2", "E. coli", "WT", "96 h", "replicate 2"),
        ("wt_24h", "E. coli", "WT", "24 h", ""), ("wt_Exp", "E. coli", "WT", "exponential", ""),
        ("HNSm_Exp", "E. coli", "Δhns", "exponential", ""), ("HNSm_96h", "E. coli", "Δhns", "96 h", ""),
        ("stpA_96h", "E. coli", "ΔstpA", "96 h", "replicate 1"), ("stpA_96h_rep2", "E. coli", "ΔstpA", "96 h", "replicate 2"),
        ("fis_96h", "E. coli", "Δfis", "96 h", ""), ("IHF_96h", "E. coli", "ΔihfAB", "96 h", ""),
        ("zapA_96h", "E. coli", "ΔzapA", "96 h", ""), ("ftsP_96h", "E. coli", "ΔftsP", "96 h", ""),
        ("envC_96h", "E. coli", "ΔenvC", "96 h", ""), ("hnsPAmCherry_96h", "E. coli", "hns-PAmCherry", "96 h", ""),
        ("bsub_WT_st_merged", "B. subtilis", "WT", "stationary", "Dugar et al. 2022, runs 1 + 2"),
        ("bsub_rok_st_merged", "B. subtilis", "Δrok", "stationary", "Dugar et al. 2022, runs 1 + 2")]
rows, counts = [], []
for name, sp, strain, cond, note in MAPS:
    T = calls(name)
    counts.append(dict(map=name, organism=sp, strain=strain, condition=cond, note=note, loops=len(T)))
    if not len(T):
        continue
    n = 461 if sp == "E. coli" else 422
    Q = pairs(T, n)
    m = matched(Q, P) if sp == "E. coli" else [np.nan] * len(Q)
    for (a, b), (_, r), mm in zip(Q, T.iterrows(), m):
        a, b = sorted((a, b))
        rows.append(dict(map=name, organism=sp, strain=strain, condition=cond, end1_start_bp=a * BIN + 1, end2_start_bp=b * BIN + 1,
                         span_kb=int(min(abs(a - b), n - abs(a - b))) * 10, confidence=r.confidence,
                         oe_percentile=round(r.oe_percentile, 3), log2_obs_exp=round(r.log2_obs_exp, 3),
                         integrated_excess=round(r.integrated_excess, 3), matches_WT_96h_loop=mm))
ALL = pd.DataFrame(rows)

README = [
    ("Supplementary Dataset 1", "DNA loops called by Bactoloops and the loop anchors they connect."),
    ("Method", "Bactoloops (rotation-invariant version, 2026-10-06), 10-kb bins, medium- and high-confidence loops (above the 95th and 99th "
               "percentile of all pixels at the same circular genomic distance); see Methods and Supplementary Fig. 2."),
    ("E. coli coordinates", "1-based bp on the lab genome (W3110, NC_007779.1, with the deletion 1,398,953-1,438,196; 4,607,088 bp). "
                            "Columns ending _W3110 give the same positions on NC_007779.1 (+39,244 bp after the deletion). "
                            "Each loop end is a 10-kb bin; start/end give the bin."),
    ("B. subtilis coordinates", "1-based bp on B. subtilis 168 (NC_000964.3), as in Dugar et al. 2022 (Mol Cell 82:804)."),
    ("WT_96h_loops", f"The {len(LP)} loops of E. coli WT at 96 h (average of two biological replicates). confidence: medium (>= 95th "
                     "percentile) or high (>= 99th); oe_percentile: rank among all pixels at the same circular distance; log2_obs_exp: "
                     "3x3-averaged log2(observed/expected) after removal of the 500-kb local background; integrated_excess: contacts "
                     "above the local background summed over 50 x 50 kb, as % of a typical bin's non-local contacts; HNS_ChIP_z: "
                     "anti-H-NS ChIP-seq z-score (96 h) of each end bin; call: main pass or ridge-region search (donut test)."),
    ("WT_96h_anchors", f"The {len(A)} anchor regions: loop ends merged across adjacent 10-kb bins. HNS_ChIP_z_max_within_10kb: strongest "
                       "anti-H-NS 96 h ChIP-seq z-score within 10 kb; HNS_bound: that value > 1. CHIDs_within_10kb: chromosomal hairpin "
                       "domains (Gavrilov, Shamovsky et al. 2025; exponential phase) whose midpoint lies within 10 kb. genes: genes "
                       "overlapping the anchor; mean_HNS_silencing_96h_log2: their mean log2(Δhns/WT) RNA-seq fold change at 96 h "
                       "(DESeq2, CFU normalization; positive = silenced by H-NS)."),
    ("Strong_HNS_regions", f"The {len(S)} runs of 10-kb bins with anti-H-NS 96 h ChIP-seq z-score >= 2 and their distance to the "
                           "nearest loop anchor (0 = overlaps an anchor)."),
    ("All_loop_calls", "Every Bactoloops call in every map analysed in the paper. matches_WT_96h_loop (E. coli only): both ends within "
                       "20 kb of the ends of a WT 96 h loop, either orientation, circular distance. Loop counts per map below."),
]
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with pd.ExcelWriter(OUT, engine="openpyxl") as W:
    R = pd.DataFrame(README, columns=["item", "description"])
    R.to_excel(W, "README", index=False)
    pd.DataFrame(counts).to_excel(W, "README", index=False, startrow=len(R) + 2)
    LP.to_excel(W, "WT_96h_loops", index=False)
    A.to_excel(W, "WT_96h_anchors", index=False)
    S.to_excel(W, "Strong_HNS_regions", index=False)
    ALL.to_excel(W, "All_loop_calls", index=False)
    for ws in W.book.worksheets:
        ws.freeze_panes = "A2"
        for col in ws.columns:
            width = max(len(str(c.value)) if c.value is not None else 0 for c in col[:60])
            ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 60 if ws.title != "README" else 120)
print(f"{len(LP)} loops, {len(A)} anchors (bound: {A.HNS_bound.sum()}, with CHID: {(A.CHIDs_within_10kb != '').sum()}), "
      f"{len(S)} strong regions ({S.anchors_a_loop.sum()} anchored; others at {sorted(S.distance_to_nearest_anchor_kb[S.distance_to_nearest_anchor_kb > 0])} kb), "
      f"{len(ALL)} calls in {sum(c['loops'] > 0 for c in counts)} maps")
print("->", OUT)
