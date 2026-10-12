# Data

Data are not stored in this repository. Put each item at the path shown (relative to the repository folder, or to `HNS_ROOT`);
`python run_all.py` then rebuilds everything downstream. Processed files that are not on GEO (marked **Zenodo**) are archived
at Zenodo: DOI to be added.

## This study (GEO GSE306451)

| path | contents | source |
|---|---|---|
| `Bactoloops/iced_matrix_conversion/<sample>_10000_iced.matrix` | HiC-Pro 10-kb ICE matrices of every *E. coli* Hi-C library | GEO / **Zenodo** |
| `Bactoloops/Hi-C raw count files/allValidPairs/<sample>.allValidPairs` | HiC-Pro valid pairs (multi-resolution test, raw-count maps for the loop-caller comparison) | **Zenodo** |
| `Bactoloops/Gavrilov2025_MicroC/FullValidPairs/<library>_EcW3110del.bwt2pairs.validPairs` | HiC-Pro valid pairs per library, without the 5-kb minimum separation (hairpin analysis; hns-PAmCherry depth-matched comparison) | **Zenodo** |
| `Bactoloops/ChIP-seq_data/<antibody>_<sample>_clean.csv`, `WGS_<sample>_clean.csv` | per-base ChIP-seq and input coverage | GEO / **Zenodo** |
| `Spike-In/Raw_Counts_Data.txt`, `HNS_metadata.txt`, `BestSizeFactors.txt`, `ODSizeFactors.txt`, `NameAndRegion.txt` | RNA-seq gene counts, sample table, spike-in size factors (CFU, OD), gene annotation on the lab genome | GEO / **Zenodo** |

All sequencing data are mapped to the lab genome: *E. coli* K-12 W3110 (GenBank AP009048.1) with bases 1,398,953-1,438,196
deleted (4,607,088 bp).

## Published data

| path | contents | source |
|---|---|---|
| `Bactoloops/iced_matrix_conversion/Dugar_HiC/GSM4288936_WT_st_HiC.cool`, `GSM4818084_WT_st_HiC_run2.cool`, `GSM4288937_Rok_st_HiC.cool`, `GSM4818085_Rok_st_HiC_run2.cool` | *B. subtilis* WT and delta-rok stationary-phase Hi-C (Dugar et al. 2022) | GEO GSE144475 |
| `Bactoloops/Dugar et al data/GSE144475_RAW/GSM42889{28,29,30,31}_*_RNAseq.bigwig` | *B. subtilis* RNA-seq coverage (Dugar et al. 2022) | GEO GSE144475 |
| `Bactoloops/ChIP-seq_data/Bsub_rok/coverage_matrix.tab` | Rok ChIP-seq and input coverage per 10-kb bin (deepTools multiBigwigSummary of the Dugar et al. Rok ChIP-seq) | **Zenodo** |
| `Bactoloops/Gavrilov2025_MicroC/opcidsRenamed.xlsx` | CHIN, CHID and OPCID coordinates (Gavrilov et al. 2025) converted to the lab genome | **Zenodo** |
| `Bactoloops/Comparisons/tools/` | `juicer_tools_1.22.01.jar` (HiCCUPS), `SIP-1.6.5.jar`, `peakachu_models/high-confidence.*.10kb.w6.pkl` | the tools' own releases |

## Made by the scripts

Everything else, including the recentred TIFF maps, the Bactoloops outputs (`Bactoloops/FIJI_output/`), ChIP-seq z-score tracks,
DESeq2 results, figures, Source Data and Supplementary Dataset 1, is made by `run_all.py`.
