"""Run the paper's analysis in dependency order. Each step runs in its own folder with HNS_ROOT set to this repository (or to
$HNS_ROOT), so data downloaded into the same folder layout (see DATA.md) is found. Steps whose inputs are missing fail and are
listed at the end; later steps that need their outputs then fail too.
Usage: python run_all.py [--list] [--from STEP] [--only STEP ...]"""
import argparse, os, subprocess, sys, time
ROOT = os.environ.get("HNS_ROOT", os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
B, UF, CMP, NUD = "Bactoloops", "Bactoloops/updated figures", "Bactoloops/Comparisons", "Bactoloops/Gavrilov2025_MicroC"
DUG, DGL, SB = "Bactoloops/Dugar et al data", "Bactoloops/Dugar et al data/loop_network_analysis", "Bactoloops/Supplemental_figure_Ecoli_vs_Bsub"
STEPS = [  # (name, folder, command)
    # 1. inputs: contact maps, ChIP-seq and RNA-seq tracks
    ("maps: E. coli ICE matrices -> TIFFs", "Bactoloops/iced_matrix_conversion", [PY, "batch_hic_matrix_to_tiff.py", "--bin-size", "10000", "--pattern", "*_10000_iced.matrix", "--recenter-bin", "368"]),
    ("maps: WT 96 h replicate average", "Bactoloops/tools", [PY, "make_wt_average.py"]),
    ("maps: B. subtilis (Dugar et al.) merged + ICE", "Bactoloops/iced_matrix_conversion/Dugar_HiC", [PY, "merge_and_ice.py"]),
    ("maps: hns-PAmCherry vs WT at common depth", "Bactoloops/tools", [PY, "common_depth_maps.py"]),
    ("maps: WT 96 h at 5/10/20 kb", "Bactoloops/iced_matrix_conversion", [PY, "wt96_multires.py"]),
    ("ChIP-seq z-score track: zscore_wt_96h_10kb.csv", "Bactoloops/ChIP-seq_data", [PY, "chip_zscore.py", "--chip", "anti-HNS_cWX2200_96hr_rep2_R1_clean.csv", "--input", "WGS_cWX2200_96hr_rep2_R1_clean.csv", "--bin-size", "10000", "--out", "zscore_wt_96h_10kb.csv"]),
    ("ChIP-seq z-score track: zscore_mCherry_96h_10kb.csv", "Bactoloops/ChIP-seq_data", [PY, "chip_zscore.py", "--chip", "anti-mCherry_cWX2882_96hr_rep1_R1_clean.csv", "--input", "WGS_cWX2882_96hr_rep1_R1_clean.csv", "--bin-size", "10000", "--out", "zscore_mCherry_96h_10kb.csv"]),
    ("ChIP-seq z-score track: zscore_mCherry_Exp_10kb.csv", "Bactoloops/ChIP-seq_data", [PY, "chip_zscore.py", "--chip", "anti-mCherry_cWX2882_Exp_rep1_R1_clean.csv", "--input", "WGS_cWX2882_Exp_rep1_R1_clean.csv", "--bin-size", "10000", "--out", "zscore_mCherry_Exp_10kb.csv"]),
    ("ChIP-seq z-score track: Rok (B. subtilis)", "Bactoloops/ChIP-seq_data/Bsub_rok", [PY, "coverage_to_zscore.py", "coverage_matrix.tab", "rok_zscore_10kb.csv"]),
    ("ChIP-seq z-score track: Rok, copy for the B. subtilis scripts", "Bactoloops", ["cp", "ChIP-seq_data/Bsub_rok/rok_zscore_10kb.csv", "Dugar et al data/rok_zscore_10kb.csv"]),
    ("RNA-seq: DESeq2 (R)", "Spike-In/rerun_raw_counts", ["Rscript", "rerun_deseq2.R"]),
    ("RNA-seq: DESeq2 results into place", "Spike-In/rerun_raw_counts", [PY, "install_results.py"]),
    ("RNA-seq: gene positions", "Spike-In", [PY, "make_gene_positions.py"]),
    # 2. Bactoloops on every map
    ("Bactoloops runs", "Bactoloops/tools", [PY, "run_bactoloops_all.py"]),
    ("Bactoloops: Supplementary Fig. 2 intermediates", "Bactoloops/tools", [PY, "s2_inputs.py"]),
    ("Bactoloops runs: 5/10/20 kb (Fiji)", "Bactoloops/FIJI_output/extra_runs/multires", ["bash", "-c", f'{PY} "{ROOT}/Bactoloops/tools/make_runs_multires.py" "$PWD" && '
        'for m in run_*.ijm; do "${FIJI:-fiji}" --allow-multiple --no-splash -- --headless -macro "$PWD/$m"; done']),
    # 3. analyses
    ("time course (24 h, peak pairs)", UF, [PY, "hns_timecourse_analysis.py"]),
    ("CHIN/CHID vs loops", B, [PY, "chin_chid_vs_loops.py"]),
    ("H-NS at loop anchors", B, [PY, "anchor_enrichment.py"]),
    ("hairpins: per-CHIN rescaled maps", NUD, [PY, "scaled_pileups.py"]),
    ("hairpins: stem test", NUD, [PY, "scaled_stem_test.py"]),
    ("hairpins: CHID/OPCID domain scores", NUD, [PY, "s9_fullpairs_stats.py"]),
    ("B. subtilis: Rok dependence of loops", DGL, [PY, "rok_dependence.py"]),
    ("B. subtilis: silencing vs H-NS", DUG, [PY, "rok_vs_hns_silencing.py"]),
    ("loop-caller comparison: inputs", CMP, [PY, "build_raw_cools.py"]),
    ("loop-caller comparison: inputs (2)", CMP, [PY, "build_inputs.py"]),
    ("loop-caller comparison: run callers", CMP, ["bash", "-c", "for t in cooltools chromosight mustache peakachu hiccups sip; do bash run_callers.sh $t; done"]),
    ("loop-caller comparison: rotated maps", "Bactoloops/Comparisons/rotation_test", [PY, "build_rotated_inputs.py"]),
    ("loop-caller comparison: run on rotated maps", "Bactoloops/Comparisons/rotation_test", ["bash", "run_rotated.sh"]),
    ("loop-caller comparison: rotation summary", "Bactoloops/Comparisons/rotation_test", [PY, "analyze_rotation.py"]),
    ("loop-caller comparison: summary", CMP, [PY, "compare_methods.py"]),
    ("robustness: 5/10/20 kb", "Bactoloops/tools", [PY, "multires_stats.py"]),
    # 4. figures, numbers, source data
    ("Fig. 2", UF, [PY, "figure2.py"]),
    ("Fig. 3", UF, [PY, "figure3_hairpins.py"]),
    ("Fig. 5 (CFU; OD and Count for Supplementary Fig. 13)", UF, ["bash", "-c", f"for n in CFU OD Count; do {PY} figure4.py $n; done"]),
    ("Supplementary Fig. 2", UF, [PY, "figureS2_method.py", "A"]),
    ("Supplementary Fig. 3", CMP, [PY, "FigS_method_comparison.py"]),
    ("Supplementary Fig. 4", UF, [PY, "figureS3_loop_conditions.py"]),
    ("Supplementary Fig. 7", UF, [PY, "figureS_hns_comparison.py"]),
    ("Supplementary Fig. 8", UF, [PY, "figureS6_chip_overlap.py"]),
    ("Supplementary Fig. 9", UF, [PY, "figureS7_timecourse.py"]),
    ("Supplementary Figs. 12 and 13", UF, [PY, "figureS10_rnaseq.py"]),
    ("Supplementary Fig. 15b (circle panel)", "Bactoloops/CircularLoopArcs", [PY, "circular_ecoli_vs_bsub.py", "--metric", "integrated_excess", "--scale", "relative", "--min-confidence", "medium", "--panel"]),
    ("Supplementary Fig. 15", SB, [PY, "fig_supplement_ecoli_vs_bsub.py"]),
    ("Fig. 7", UF, [PY, "figure7.py"]),
    ("Supplementary Note 2 (hns-PAmCherry)", UF, [PY, "figure_temp_pamcherry_common.py"]),
    ("Supplementary Fig. 5 tracks", "Bactoloops/tools", [PY, "s5_tracks.py"]),
    ("numbers in the text", "Bactoloops/tools", [PY, "manuscript_numbers.py"]),
    ("Supplementary Dataset 1", "Bactoloops/tools", [PY, "make_supplementary_dataset1.py"]),
    ("Supplementary Dataset 2", "Bactoloops/tools", [PY, "make_supplementary_dataset2.py"]),
    ("Source Data", "Bactoloops/tools", [PY, "make_source_data.py"]),
    ("high-resolution figure copies", "Bactoloops/tools", [PY, "make_high_res_figures.py"]),
]
ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--list", action="store_true"); ap.add_argument("--from", dest="start"); ap.add_argument("--only", nargs="*")
a = ap.parse_args()
if a.list:
    for i, (n, d, c) in enumerate(STEPS, 1): print(f"{i:2d}  {n:60s} {d}")
    sys.exit()
env = dict(os.environ, HNS_ROOT=ROOT)
fails, on = [], a.start is None
for n, d, c in STEPS:
    on = on or (a.start is not None and a.start.lower() in n.lower())
    if not on or (a.only and not any(o.lower() in n.lower() for o in a.only)): continue
    t0 = time.time(); print(f"== {n}", flush=True)
    os.makedirs(os.path.join(ROOT, d), exist_ok=True)
    r = subprocess.run(c, cwd=os.path.join(ROOT, d), env=env)
    print(f"   {'ok' if r.returncode == 0 else 'FAILED'} ({time.time() - t0:.0f} s)", flush=True)
    if r.returncode: fails.append(n)
print("\nfailed steps:" if fails else "\nall steps ran", *fails, sep="\n  ")
