"""Put the DESeq2 results of rerun_deseq2.R where the figure scripts read them: Spike-In/DESeq2_apeglm_<CFU|OD|Count>_<Exp|96h>.txt,
with the gene name as a GeneID column and genes sorted by adjusted P (genes without one last). Values are copied unchanged."""
import pandas as pd
import os as _hns_os
def _hns_root():
    r = _hns_os.environ.get("HNS_ROOT")
    if r:
        return _hns_os.path.join(r, "")
    d = _hns_os.path.dirname(_hns_os.path.abspath(__file__))
    while not _hns_os.path.exists(_hns_os.path.join(d, ".hns_root")) and _hns_os.path.dirname(d) != d:
        d = _hns_os.path.dirname(d)
    return _hns_os.path.join(d, "")
SP = _hns_root() + "Spike-In/"
for nm in ["CFU", "OD", "Count"]:
    for tp in ["Exp", "96h"]:
        r = pd.read_csv(SP + f"rerun_raw_counts/DESeq2_apeglm_{nm}_{tp}_raw.txt", sep="\t").reset_index().rename(columns={"index": "GeneID"})
        r.sort_values("padj", kind="mergesort", na_position="last").to_csv(SP + f"DESeq2_apeglm_{nm}_{tp}.txt", sep="\t", index=False)
print("installed 6 DESeq2 result tables")
