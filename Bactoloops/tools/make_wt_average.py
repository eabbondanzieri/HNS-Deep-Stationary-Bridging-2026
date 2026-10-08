"""The WT 96 h replicate-average map used throughout the paper: the mean of the two replicates' 10-kb ICE TIFFs.
Output: iced_matrix_conversion/cWX2200_96h_rep1rep2_average_10000_iced.tif"""
import tifffile
import os as _hns_os
def _hns_root():
    r = _hns_os.environ.get("HNS_ROOT")
    if r:
        return _hns_os.path.join(r, "")
    d = _hns_os.path.dirname(_hns_os.path.abspath(__file__))
    while not _hns_os.path.exists(_hns_os.path.join(d, ".hns_root")) and _hns_os.path.dirname(d) != d:
        d = _hns_os.path.dirname(d)
    return _hns_os.path.join(d, "")
IC = _hns_root() + "Bactoloops/iced_matrix_conversion/"
r1 = tifffile.imread(IC + "cWX2200_96h_Azure_10percent_Sau3A_rep1_10000_iced.tif")
r2 = tifffile.imread(IC + "cWX2200_96H_Azure_10percent_Sau3AI_rep2_10000_iced.tif")
tifffile.imwrite(IC + "cWX2200_96h_rep1rep2_average_10000_iced.tif", ((r1 + r2) / 2).astype("float32"))
print("wrote the WT 96 h replicate average")
