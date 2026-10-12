# HNS-Deep-Stationary-Bridging-2026

Analysis code for

> Way, L. E. *et al.* **H-NS organizes genome-wide DNA looping linked to gene silencing in starved *Escherichia coli*.**
> bioRxiv 2025.10.14.682352 (2025); doi:10.1101/2025.10.14.682352.

This repository holds the Hi-C, ChIP-seq and RNA-seq analyses and the figures made from them. The loop caller, **Bactoloops**,
is maintained separately at [eabbondanzieri/bactoloops](https://github.com/eabbondanzieri/bactoloops); the version used for the
paper is included here (`Bactoloops/HiC_PeakFindingAlgorithm_with_symmetrization.ijm`, `Bactoloops/python/`) so that the
analysis is reproducible as published. Other code for the paper: sequencing data processing,
[xindanwanglab/Way-HNS-2026](https://github.com/xindanwanglab/Way-HNS-2026); microscopy and single-particle tracking,
[BiteenMatlab/H-NS_DNA_bridging](https://github.com/BiteenMatlab/H-NS_DNA_bridging).

## Layout

Code and data share one folder tree (the folder names are those used during the project). Data are not in the repository:
download them into the same tree (see [DATA.md](DATA.md)) or point `HNS_ROOT` at a folder that has this layout.

| folder | contents |
|---|---|
| `Bactoloops/iced_matrix_conversion/` | contact maps: HiC-Pro ICE matrices -> recentred 10-kb TIFFs; WT replicate average; depth-matched maps; *B. subtilis* maps (`Dugar_HiC/`) |
| `Bactoloops/ChIP-seq_data/` | ChIP-seq enrichment z-score tracks (H-NS, H-NS-PAmCherry; Rok in `Bsub_rok/`) |
| `Spike-In/` | RNA-seq: DESeq2 with ERCC spike-ins (R), gene positions |
| `Bactoloops/tools/` | running Bactoloops on every map, numbers quoted in the text, Source Data, Supplementary Dataset 1 |
| `Bactoloops/updated figures/` | Figs. 2, 3, 5 and Supplementary Figs. 2, 4, 7, 8, 9, 12; shared figure style (`house_style.py`) |
| `Bactoloops/Comparisons/` | comparison with six published loop callers, including the rotation test (Supplementary Note 1, Supplementary Fig. 3) |
| `Bactoloops/Nudler CHIN paper/` | short-range analysis of chromosomal hairpins (CHINs), hairpin domains (CHIDs) and OPCIDs from full valid pairs (Fig. 3) |
| `Bactoloops/Dugar et al data/`, `Bactoloops/Supplemental_figure_Ecoli_vs_Bsub/`, `Bactoloops/CircularLoopArcs/` | reanalysis of *B. subtilis* Hi-C, Rok ChIP-seq, and RNA-seq (Dugar et al. 2022); Supplementary Fig. 15 |
| `Bactoloops/` (top level) | shared statistics: exact rotation and circular-shift tests (`rotation_test.py`), baseline-corrected integrated excess (`excess_baseline.py`), H-NS at loop anchors, CHIN/CHID vs loops |

## Figures and scripts

| figure | script (inputs made by) |
|---|---|
| Fig. 2 | `Bactoloops/updated figures/figure2.py` (cell drawings in `Bactoloops/Cell cartoons/`) |
| Fig. 3 | `Bactoloops/updated figures/figure3_hairpins.py` (drawings in `Cell cartoons/`; `Nudler CHIN paper/scaled_pileups.py`, `scaled_stem_test.py`, `s9_fullpairs_stats.py`; `updated figures/hns_timecourse_analysis.py`; `chin_chid_vs_loops.py`) |
| Fig. 5 | `Bactoloops/updated figures/figure4.py` (file name from an earlier numbering; `Spike-In/rerun_raw_counts/rerun_deseq2.R`, `Spike-In/make_gene_positions.py`) |
| Fig. 7 | `Bactoloops/updated figures/figure7.py` (cell drawings in `Bactoloops/Cell cartoons/`) |
| Supplementary Fig. 2 | `Bactoloops/updated figures/figureS2_method.py` (`tools/s2_inputs.py`) |
| Supplementary Fig. 3, Supplementary Note 1 | `Bactoloops/Comparisons/FigS_method_comparison.py` (`build_raw_cools.py`, `build_inputs.py`, `run_callers.sh`, `rotation_test/`, `compare_methods.py`) |
| Supplementary Fig. 4 | `Bactoloops/updated figures/figureS3_loop_conditions.py` |
| Supplementary Fig. 5 (tracks) | `Bactoloops/tools/s5_tracks.py` |
| Supplementary Fig. 7 | `Bactoloops/updated figures/figureS_hns_comparison.py` |
| Supplementary Fig. 8 | `Bactoloops/updated figures/figureS6_chip_overlap.py` |
| Supplementary Fig. 9 | `Bactoloops/updated figures/figureS7_timecourse.py` (`hns_timecourse_analysis.py`) |
| Supplementary Figs. 12 and 13 | `Bactoloops/updated figures/figureS10_rnaseq.py` (`figure4.py` run with `CFU`, `OD`, and `Count`) |
| Supplementary Fig. 15 | `Bactoloops/Supplemental_figure_Ecoli_vs_Bsub/fig_supplement_ecoli_vs_bsub.py` (`CircularLoopArcs/circular_ecoli_vs_bsub.py`; `Dugar et al data/`) |
| Supplementary Note 2 (hns-PAmCherry looping) | `Bactoloops/tools/common_depth_maps.py`, `Bactoloops/updated figures/figure_temp_pamcherry_common.py` |
| Supplementary Methods 1 (robustness) | `Bactoloops/iced_matrix_conversion/wt96_multires.py`, `tools/make_runs_multires.py`, `tools/multires_stats.py` |
| numbers in the text | `Bactoloops/tools/manuscript_numbers.py`, `anchor_enrichment.py`, `chin_chid_vs_loops.py` |
| Source Data; Supplementary Datasets 1 and 2; high-resolution figure copies | `Bactoloops/tools/make_source_data.py`; `Bactoloops/tools/make_supplementary_dataset1.py`, `make_supplementary_dataset2.py`; `make_high_res_figures.py` |

Figs. 1, 4, and 6 and Supplementary Figs. 1, 5 (other than the tracks), 6, 10, 11, and 14 show microscopy, Western blots,
and Hi-C maps prepared with the lab pipeline; their code is in the repositories linked above.

## Running

```bash
conda env create -f environment.yml && conda activate hns-paper
python run_all.py --list          # the 52 steps, in dependency order
python run_all.py                 # everything (hours; the loop-caller comparison takes longest)
python run_all.py --from "Fig. 2" # from a given step onwards
python run_all.py --only "Fig. 3" "Supplementary Fig. 9"
```

Each script also runs on its own from its folder. Scripts find the data through `HNS_ROOT` (an environment variable; by
default the folder holding this README, marked by `.hns_root`).

Bactoloops runs use the Python implementation by default (`tools/run_bactoloops_all.py`), which reproduces the ImageJ macro's
calls; `--engine fiji` runs the macro itself (set `FIJI=/path/to/fiji`). The multi-resolution test runs the macro in Fiji.

**Other software.** R with DESeq2 and apeglm (RNA-seq). The loop-caller comparison needs the published callers: cooltools,
Chromosight and Mustache (`LOOPCOMPARE_BIN`), HiCExplorer and cooler (`HICEXPLORER_BIN`, `HICEXPLORER_PYTHON`), Peakachu with its
pretrained 10-kb model (`PEAKACHU_BIN`), and the Juicer tools and SIP jars in `Bactoloops/Comparisons/tools/` (see DATA.md).

## Earlier version

`previous_submission/HiC_PeakFindingAlgorithm.ijm` is the macro posted with the first submission (integrated loop density
only, before loop calling and the circular treatment of the map). It is kept for reference; the analysis uses the current
version.

## License

MIT (see `LICENSE`).
