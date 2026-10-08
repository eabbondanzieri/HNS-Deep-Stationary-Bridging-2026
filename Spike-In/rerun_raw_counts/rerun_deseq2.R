# Re-run of DEseq2_wSpikeIn_TotalAnalysis.R (same model: all 12 samples, design ~ Condition, custom size factors,
# apeglm shrinkage of WT vs dHNS at each time point) on two count tables:
#   raw    Raw_Counts_Data.txt  (CLC "Total gene reads", Mapped to EcW3110del) -- what the Excel size factors were made for
#   old    spike_normalized_named.tabular  (raw counts x a within-group median-of-ratios factor, rounded) -- the original input
# Normalizations: CFU = BestSizeFactors.txt, OD = ODSizeFactors.txt, Count = DESeq2 median-of-ratios (estimateSizeFactors).
# Output: rerun_raw_counts/DESeq2_apeglm_<norm>_<tp>_<counts>.txt (positive log2FC = higher in WT, as before)
suppressMessages({library(DESeq2); library(apeglm)})
setwd(file.path(Sys.getenv("HNS_ROOT", unset = "../.."), "Spike-In"))
rd <- function(f) { x <- readLines(f, warn = FALSE); x <- unlist(strsplit(x, "\r")); read.delim(text = x, row.names = 1, check.names = FALSE) }
meta <- rd("HNS_metadata.txt")
sfs <- list(CFU = rd("BestSizeFactors.txt"), OD = rd("ODSizeFactors.txt"), Count = NULL)
for (ct in c("raw", "old")) {
  counts <- as.matrix(rd(if (ct == "raw") "Raw_Counts_Data.txt" else "spike_normalized_named.tabular"))
  storage.mode(counts) <- "integer"
  m <- meta[colnames(counts), , drop = FALSE]
  for (nm in names(sfs)) for (tp in c("Exp", "96h")) {
    m$Condition <- relevel(factor(m$Condition), ref = paste0("dHNS_", tp))
    dds <- DESeqDataSetFromMatrix(counts, m, ~ Condition)
    if (is.null(sfs[[nm]])) dds <- estimateSizeFactors(dds) else sizeFactors(dds) <- sfs[[nm]][colnames(dds), "SizeFactor"]
    dds <- DESeq(dds, quiet = TRUE)
    res <- lfcShrink(dds, coef = paste0("Condition_WT_", tp, "_vs_dHNS_", tp), type = "apeglm", quiet = TRUE)
    write.table(as.data.frame(res), sprintf("rerun_raw_counts/DESeq2_apeglm_%s_%s_%s.txt", nm, tp, ct), sep = "\t", quote = FALSE)
    cat(ct, nm, tp, "up", sum(res$padj < 0.05 & res$log2FoldChange >= 1, na.rm = TRUE),
        "down", sum(res$padj < 0.05 & res$log2FoldChange <= -1, na.rm = TRUE), "\n")
  }
}
