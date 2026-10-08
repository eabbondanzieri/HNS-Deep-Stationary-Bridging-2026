"""Supplementary Dataset 2: RNA-seq of WT and delta-hns cells in exponential phase and at 96 h, gene by gene.

Sheets
  README                    column definitions, normalizations and counts
  RNAseq_all_genes          every gene: position, H-NS binding, loop anchors, and the DESeq2 results (WT vs delta-hns) for both growth
                            phases and the three normalizations (CFU, used in the main text; OD; read count)
  Flagellar_chemotaxis      the clusters of flagellar and chemotaxis genes with higher RNA levels in WT cells in exponential phase
                            (Fig. 5c, orange stars)

Inputs: Spike-In/genome_positions/DESeq2_apeglm_<CFU|OD|Count>_<Exp|96h>_with_positions.txt (raw-count DESeq2 with apeglm
shrinkage, 2026-10-01; gene positions on the lab genome), ChIP-seq_data/zscore_wt_96h_10kb.csv, Supplementary Dataset 1 (anchors).
Output: Paper Update/Supplementary Data/Supplementary_Dataset_2_RNA-seq.xlsx
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
SP = os.path.dirname(B.rstrip("/")) + "/Spike-In/"
OUT = B + "Paper Update/Supplementary Data/Supplementary_Dataset_2_RNA-seq.xlsx"
D1 = B + "Paper Update/Supplementary Data/Supplementary_Dataset_1_loops_and_anchors.xlsx"
GENOME, ORIC, BIN, N = 4_607_088, 3_671_577, 10_000, 461
PADJ, LFC = 0.05, 1.0
NORMS = [("CFU", "spike-in reads x CFU collected (main text)"), ("OD", "spike-in reads x OD units collected"),
         ("Count", "DESeq2 median-of-ratios on gene counts, no spike-ins")]
PHASES = [("Exp", "exp"), ("96h", "96h")]
# Fig. 5c orange stars: clusters of genes higher in WT in exponential phase, on the figure's genome axis (kb from oriC,
# positions after oriC wrapped to negative values), as STAR_WINDOWS in figure4.py; the flh/che and fli clusters share one star
CLUSTERS = [("flg", 1131, 1142), ("flh/che", 1925, 1940), ("fli", 1964, 1986)]

res = {}
for nm, _ in NORMS:
    for tp, _ in PHASES:
        r = pd.read_csv(SP + f"genome_positions/DESeq2_apeglm_{nm}_{tp}_with_positions.txt", sep="\t").set_index("GeneID")
        chk = pd.read_csv(SP + f"DESeq2_apeglm_{nm}_{tp}.txt", sep="\t").set_index("GeneID")
        assert np.allclose(r.log2FoldChange, chk.loc[r.index].log2FoldChange, equal_nan=True), (nm, tp)   # current DESeq2 results
        res[(nm, tp)] = r
base = res[("CFU", "96h")]
G = pd.DataFrame(index=base.index)
G["strand"] = base.Strand
G["start_bp"], G["end_bp"] = base[["Gene_start_lab", "Gene_end_lab"]].min(1).astype("Int64"), base[["Gene_start_lab", "Gene_end_lab"]].max(1).astype("Int64")
G["start_bp_W3110"], G["end_bp_W3110"] = base[["Gene_start_W3110", "Gene_end_W3110"]].min(1).astype("Int64"), base[["Gene_start_W3110", "Gene_end_W3110"]].max(1).astype("Int64")
G["gene_start_5prime_bp"] = base.Gene_start_lab.astype("Int64")

# H-NS binding of the gene's 10-kb bin (5' end), and WT 96 h loop anchors (Supplementary Dataset 1)
c = pd.read_csv(B + "ChIP-seq_data/zscore_wt_96h_10kb.csv")
zb = np.full(N, np.nan); zb[(c.start_bp // BIN).astype(int) % N] = c.zscore.values
gb = ((G.gene_start_5prime_bp.astype(float) - 1) // BIN)
G["HNS_ChIP_z_96h_bin"] = [zb[int(b) % N] if np.isfinite(b) else np.nan for b in gb]
A = pd.read_excel(D1, sheet_name="WT_96h_anchors")
anc = []
for s, e in zip(G.start_bp.astype(float), G.end_bp.astype(float)):
    hit = A[(A.start_bp <= e) & (A.end_bp >= s)] if np.isfinite(s) else A.iloc[:0]
    anc.append(", ".join(hit.anchor_id.astype(str)))
G["WT_96h_loop_anchor"] = anc

# Fig. 5c clusters
ax = np.where(G.gene_start_5prime_bp.astype(float) > ORIC, G.gene_start_5prime_bp.astype(float) - GENOME, G.gene_start_5prime_bp.astype(float)) / 1000
G["Fig5c_cluster"] = ""
for name, lo, hi in CLUSTERS:
    G.loc[(ax >= lo) & (ax <= hi), "Fig5c_cluster"] = name

for tp, tl in PHASES:
    for nm, _ in NORMS:
        r = res[(nm, tp)].reindex(G.index)          # same gene order as G (the call column below is positional)
        G[f"{tl}_{nm}_baseMean"] = r.baseMean
        G[f"{tl}_{nm}_log2FC_WT_vs_dhns"] = r.log2FoldChange
        G[f"{tl}_{nm}_lfcSE"] = r.lfcSE
        G[f"{tl}_{nm}_padj"] = r.padj
        G[f"{tl}_{nm}_call"] = np.where((r.padj < PADJ) & (r.log2FoldChange >= LFC), "higher in WT",
                                        np.where((r.padj < PADJ) & (r.log2FoldChange <= -LFC), "higher in Δhns", ""))
G = G.reset_index().rename(columns={"GeneID": "gene"})
G = G.sort_values("start_bp", na_position="last").reset_index(drop=True)
F = G[G.Fig5c_cluster != ""][["gene", "Fig5c_cluster", "strand", "start_bp", "end_bp", "start_bp_W3110", "end_bp_W3110",
                             "exp_CFU_log2FC_WT_vs_dhns", "exp_CFU_padj", "exp_CFU_call", "96h_CFU_log2FC_WT_vs_dhns", "96h_CFU_padj", "96h_CFU_call"]]

counts = []
for tp, tl in PHASES:
    for nm, desc in NORMS:
        cl = G[f"{tl}_{nm}_call"]
        counts.append(dict(growth_phase="exponential" if tp == "Exp" else "96 h", normalization=nm,
                           higher_in_WT=int((cl == "higher in WT").sum()), higher_in_dhns=int((cl == "higher in Δhns").sum())))
C = pd.DataFrame(counts)
README = [
    ("Supplementary Dataset 2", "RNA-seq of E. coli WT (cWX2200) and Δhns (cWX2870) cells in exponential phase and at 96 h: DESeq2 results "
                                "for every gene."),
    ("Method", "Three biological replicates per strain and condition; ERCC RNA spike-ins. DESeq2 on raw gene counts with apeglm "
               "log2 fold-change shrinkage (see Methods). Fold changes are WT / Δhns (log2): positive = higher RNA level in WT, negative = "
               "higher in Δhns (silenced by H-NS; H-NS silencing in the paper is the negative of this value). Calls: padj < 0.05 and "
               "|log2 fold change| >= 1."),
    ("Normalizations", "; ".join(f"{nm}: {d}" for nm, d in NORMS) + ". The main text uses CFU; Supplementary Fig. 12 compares the three. "
                       "The number of genes called differs between normalizations at 96 h because the read-count normalization removes "
                       "the global difference in RNA per cell by design."),
    ("Coordinates", "1-based bp on the lab genome (W3110, NC_007779.1, with the deletion 1,398,953-1,438,196; 4,607,088 bp); columns "
                    "ending _W3110 give NC_007779.1 positions. gene_start_5prime_bp: the start-codon end (strand-aware)."),
    ("HNS_ChIP_z_96h_bin", "anti-H-NS ChIP-seq z-score (WT, 96 h) of the 10-kb bin containing the gene's 5' end."),
    ("WT_96h_loop_anchor", "WT 96 h loop anchor(s) overlapping the gene (Supplementary Dataset 1, sheet WT_96h_anchors)."),
    ("Fig5c_cluster", "flagellar and chemotaxis gene clusters with higher RNA levels in WT cells in exponential phase (Fig. 5c, orange "
                      "stars; the flh/che and fli clusters, ~40 kb apart, share one star): flg, flh/che (with tap, tar, mot), fli."),
    ("Columns per condition", "<phase>_<normalization>_baseMean, _log2FC_WT_vs_dhns, _lfcSE, _padj, _call (phase: exp or 96h)."),
    ("hns", "hns itself is deleted in Δhns; its fold change is not meaningful and it is left out of the genome-wide analyses."),
]
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with pd.ExcelWriter(OUT, engine="openpyxl") as W:
    R = pd.DataFrame(README, columns=["item", "description"])
    R.to_excel(W, sheet_name="README", index=False)
    C.to_excel(W, sheet_name="README", index=False, startrow=len(R) + 2)
    G.to_excel(W, sheet_name="RNAseq_all_genes", index=False)
    F.to_excel(W, sheet_name="Flagellar_chemotaxis", index=False)
    for ws in W.book.worksheets:
        ws.freeze_panes = "B2" if ws.title != "README" else "A2"
        for col in ws.columns:
            width = max(len(str(c_.value)) if c_.value is not None else 0 for c_ in col[:60])
            ws.column_dimensions[col[0].column_letter].width = min(max(10, width + 2), 40 if ws.title != "README" else 120)
print(f"{len(G)} genes; counts:\n{C.to_string(index=False)}")
print("Fig. 5c clusters:", {k: int((G.Fig5c_cluster == k).sum()) for k, _, _ in CLUSTERS},
      "| of these higher in WT (exp, CFU):", int(((G.Fig5c_cluster != "") & (G.exp_CFU_call == "higher in WT")).sum()))
print("genes at WT 96 h loop anchors:", int((G.WT_96h_loop_anchor != "").sum()))
print("->", OUT)
