"""
coverage_to_zscore.py

Converts the deeptools multiBigwigSummary output (raw per-10kb-bin coverage
for Rok ChIP and input) into the same z-score format used throughout this
project for the H-NS ChIP-seq data:

    bin, start_bp, end_bp, chip_mean, input_mean, ratio, zscore

- chip_mean / input_mean: mean per-base coverage in that 10kb bin
- ratio: chip_mean / input_mean
- zscore: the ratio, z-scored (mean 0, std 1) across all bins genome-wide

Every bin in the coverage matrix is kept, so the output has one row per
genome bin and `bin` always equals start_bp // bin size + 1. Bins where the
input coverage is missing, zero, or below MIN_INPUT_FRACTION of the median
input (too few reads for a reliable ratio) get NaN for ratio and zscore, and
are left out of the mean/std used for z-scoring. This matches chip_zscore.py.

In the Rok data this masks two bins: 480-490 kb, where neither the ChIP nor
the input BAM has any reads (476-491 kb is empty in both), and 1,970-1,980 kb,
where the input has ~5% of its usual coverage.

Usage:
    python coverage_to_zscore.py coverage_matrix.tab rok_zscore_10kb.csv
"""

import sys
import numpy as np
import pandas as pd

MIN_INPUT_FRACTION = 0.2   # input below this fraction of the genome-wide median input -> NaN

def main():
    if len(sys.argv) != 3:
        print("Usage: python coverage_to_zscore.py <coverage_matrix.tab> <output.csv>")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    # multiBigwigSummary --outRawCounts format: a header line starting with
    # "#" listing chrom/start/end/<bigwig filenames>, then one row per bin.
    with open(in_path) as f:
        header_line = f.readline().lstrip("#").strip()
    col_names = header_line.split("\t")

    df = pd.read_csv(in_path, sep="\t", skiprows=1, names=col_names)

    # The last two columns are the two bigwig tracks, in the order given to
    # multiBigwigSummary (rok_chip.bw first, input.bw second per the command above)
    chip_col, input_col = col_names[-2], col_names[-1]

    df = df.rename(columns={col_names[0]: "chrom", col_names[1]: "start_bp", col_names[2]: "end_bp"})
    df["chip_mean"] = pd.to_numeric(df[chip_col], errors="coerce")
    df["input_mean"] = pd.to_numeric(df[input_col], errors="coerce")

    # Keep every bin; mask (NaN) the ones whose input is missing, zero, or too
    # low for a reliable ratio, rather than dropping them -- dropping a row
    # shifts every later row off its genome bin.
    df = df.sort_values("start_bp").reset_index(drop=True)
    min_input = MIN_INPUT_FRACTION * df["input_mean"].median()
    ok = df["chip_mean"].notna() & df["input_mean"].notna() & (df["input_mean"] >= min_input)
    for _, r in df[~ok].iterrows():
        print(f"Masked {int(r.start_bp)}-{int(r.end_bp)}: chip={r.chip_mean}, input={r.input_mean} "
              f"(input threshold {min_input:.0f})")

    df["ratio"] = np.where(ok, df["chip_mean"] / df["input_mean"], np.nan)
    mean_ratio = df.loc[ok, "ratio"].mean()
    std_ratio = df.loc[ok, "ratio"].std()
    df["zscore"] = (df["ratio"] - mean_ratio) / std_ratio

    bin_size = int(df["end_bp"].iloc[0] - df["start_bp"].iloc[0])
    df["bin"] = df["start_bp"] // bin_size + 1  # 1-indexed, matching the H-NS z-score file's convention

    out = df[["bin", "start_bp", "end_bp", "chip_mean", "input_mean", "ratio", "zscore"]]
    out.to_csv(out_path, index=False)

    print(f"Wrote {len(out)} bins ({(~ok).sum()} masked) to {out_path}")
    print(f"Mean ratio: {mean_ratio:.4f}, std: {std_ratio:.4f}")
    print(out.head())


if __name__ == "__main__":
    main()
