#!/usr/bin/env python3
"""
chip_zscore.py

Computes a ChIP-seq (or any IP) enrichment z-score track from two per-base
coverage files: a ChIP/IP track and a WGS/input track used to normalize it.
This is the same pipeline used to build the anti-HNS/WGS z-score track
compared against the Hi-C loop calls: bin both tracks, divide ChIP by input
bin-for-bin, then z-score the resulting ratio.

PIPELINE
--------
    1. Bin both tracks to --bin-size bp (mean coverage per bin; bin 1 =
       positions 1..bin_size, bin 2 = bin_size+1..2*bin_size, etc. -- the
       last bin may be partial if the genome length isn't an exact multiple
       of bin_size).
    2. ratio = ChIP_bin_mean / Input_bin_mean, bin by bin.
    3. z-score = (ratio - mean(ratio)) / std(ratio), over all finite bins
       (a bin with zero input coverage gives a non-finite ratio and is
       excluded from the mean/std, left as NaN in the output rather than
       silently dropped).

INPUT FORMAT
------------
A two-column CSV with a header -- "Position,Value" (or any two column
names; the first two columns are used positionally). Position is 1-based
and contiguous, one row per base pair -- the format FIJI/ImageJ's "Plot
Profile"/coverage export and similar tools produce.

USAGE
-----
    Single pair:
        python3 chip_zscore.py --chip anti-HNS_sample.csv --input WGS_sample.csv \\
            --bin-size 10000 --out sample_zscore_10kb.csv

    Batch -- process every ChIP/Input PAIR found in a folder. Pairing is done
    by filename: a ChIP file "<chip-prefix><suffix>" is paired with the
    Input file "<input-prefix><suffix>" in the same folder (same <suffix>) --
    e.g. "anti-HNS_cWX2200_96hr_rep2_R1_clean.csv" pairs with
    "WGS_cWX2200_96hr_rep2_R1_clean.csv" by default.
        python3 chip_zscore.py --batch-dir /path/to/folder \\
            --chip-prefix "anti-HNS_" --input-prefix "WGS_" \\
            --bin-size 10000 --out-dir /path/to/folder/zscores

Requires: numpy, pandas (pip install numpy pandas)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load_coverage(path):
    df = pd.read_csv(path)
    if df.shape[1] < 2:
        sys.exit(f"{path}: expected at least 2 columns (Position, Value), got {df.shape[1]}")
    df = df.iloc[:, :2].copy()
    df.columns = ["Position", "Value"]
    df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
    df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
    bad = df[["Position", "Value"]].isna().any(axis=1)
    if bad.any():
        first_bad = df.loc[bad].iloc[0]
        sys.exit(f"{path}: {int(bad.sum())} row(s) aren't numeric Position/Value pairs -- "
                  f"is this really a per-base coverage file? First bad row "
                  f"(0-indexed row {bad.idxmax()}): {first_bad.to_dict()}")
    return df


def bin_track(df, bin_size, label):
    df = df.copy()
    # 1-based bin number: bin 1 = positions 1..bin_size, bin 2 = bin_size+1..2*bin_size, etc.
    df["bin"] = ((df["Position"] - 1) // bin_size + 1).astype(np.int64)
    return df.groupby("bin")["Value"].mean().rename(f"{label}_mean")


def compute_zscore(chip_file, input_file, bin_size, out_file, ddof=1):
    chip = load_coverage(chip_file)
    inp = load_coverage(input_file)

    if (chip["Position"].min() != inp["Position"].min()
            or chip["Position"].max() != inp["Position"].max()):
        print(
            f"WARNING: {chip_file} and {input_file} don't cover the same Position "
            f"range ({chip['Position'].min():.0f}-{chip['Position'].max():.0f} vs "
            f"{inp['Position'].min():.0f}-{inp['Position'].max():.0f}) -- proceeding "
            f"anyway, but double check these are really the same sample/genome.",
            file=sys.stderr,
        )

    chip_binned = bin_track(chip, bin_size, "chip")
    input_binned = bin_track(inp, bin_size, "input")

    df = pd.concat([chip_binned, input_binned], axis=1).sort_index().reset_index()
    n_bins = len(df)
    genome_length = int(max(chip["Position"].max(), inp["Position"].max()))

    df["start_bp"] = (df["bin"] - 1) * bin_size
    df["end_bp"] = np.minimum(df["bin"] * bin_size, genome_length)

    df["ratio"] = df["chip_mean"] / df["input_mean"]

    finite = np.isfinite(df["ratio"])
    n_bad = int((~finite).sum())
    if n_bad:
        print(
            f"{Path(chip_file).name}: {n_bad} bin(s) have a non-finite ratio (zero "
            f"input coverage) -- excluded from the mean/std used for the z-score, "
            f"left as NaN in the output.",
            file=sys.stderr,
        )

    mu = df.loc[finite, "ratio"].mean()
    sigma = df.loc[finite, "ratio"].std(ddof=ddof)
    df["zscore"] = (df["ratio"] - mu) / sigma

    out_cols = ["bin", "start_bp", "end_bp", "chip_mean", "input_mean", "ratio", "zscore"]
    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    df[out_cols].to_csv(out_file, index=False)
    print(
        f"{Path(chip_file).name} / {Path(input_file).name} -> {out_file}  "
        f"({n_bins} bins, ratio mean={mu:.4f}, std={sigma:.4f})"
    )
    return df


def find_batch_pairs(batch_dir, chip_prefix, input_prefix):
    batch_dir = Path(batch_dir)
    chip_files = sorted(
        p for p in batch_dir.iterdir() if p.is_file() and p.name.startswith(chip_prefix)
    )

    pairs, unmatched = [], []
    for chip_path in chip_files:
        suffix = chip_path.name[len(chip_prefix):]
        input_path = batch_dir / f"{input_prefix}{suffix}"
        if input_path.exists():
            pairs.append((chip_path, input_path, suffix))
        else:
            unmatched.append(chip_path)

    return pairs, unmatched


def main():
    p = argparse.ArgumentParser(
        description="Compute a ChIP-seq (or similar IP) enrichment z-score track from a "
        "ChIP coverage file normalized against a WGS/input coverage file, binned "
        "to a given resolution. Runs on a single file pair, or in batch over "
        "every matching pair found in a folder."
    )
    p.add_argument("--chip", help="ChIP/IP per-base coverage CSV (Position,Value)")
    p.add_argument("--input", help="WGS/input per-base coverage CSV (Position,Value), normalizes --chip")
    p.add_argument("--out", help="Output CSV path (single-pair mode)")

    p.add_argument("--batch-dir", help="Folder to scan for ChIP/Input pairs (batch mode)")
    p.add_argument(
        "--chip-prefix", default="anti-HNS_",
        help="Filename prefix identifying a ChIP file in batch mode (default: 'anti-HNS_')",
    )
    p.add_argument(
        "--input-prefix", default="WGS_",
        help="Filename prefix identifying the matching Input/WGS file in batch mode (default: 'WGS_')",
    )
    p.add_argument("--out-dir", help="Output folder for batch mode (default: --batch-dir itself)")
    p.add_argument(
        "--out-suffix", default="_zscore",
        help="Inserted before the resolution and extension in batch-mode output filenames "
        "(default: '_zscore', e.g. anti-HNS_sample_zscore_10kb.csv)",
    )

    p.add_argument("--bin-size", type=int, required=True, help="Bin size in bp (e.g. 10000)")
    p.add_argument(
        "--ddof", type=int, default=1,
        help="Delta degrees of freedom for the std used in the z-score (default: 1, sample std)",
    )

    args = p.parse_args()

    single_mode = bool(args.chip or args.input or args.out)
    batch_mode = args.batch_dir is not None

    if single_mode and batch_mode:
        sys.exit("Use either --chip/--input/--out (single pair) or --batch-dir (batch mode), not both.")

    if batch_mode:
        pairs, unmatched = find_batch_pairs(args.batch_dir, args.chip_prefix, args.input_prefix)
        if not pairs:
            sys.exit(
                f"No ChIP/Input pairs found in {args.batch_dir} matching --chip-prefix "
                f"'{args.chip_prefix}' / --input-prefix '{args.input_prefix}'."
            )
        out_dir = Path(args.out_dir) if args.out_dir else Path(args.batch_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"Found {len(pairs)} pair(s); {len(unmatched)} unmatched ChIP file(s).")
        for chip_path, input_path, suffix in pairs:
            suffix_stem = Path(suffix).stem
            out_name = f"{args.chip_prefix}{suffix_stem}{args.out_suffix}_{args.bin_size}bp.csv"
            compute_zscore(chip_path, input_path, args.bin_size, out_dir / out_name, ddof=args.ddof)

        if unmatched:
            print("\nChIP files with no matching Input file (skipped):")
            for u in unmatched:
                print(f"  {u.name}")

    else:
        if not (args.chip and args.input and args.out):
            sys.exit("Single-pair mode requires --chip, --input, and --out (or use --batch-dir instead).")
        compute_zscore(Path(args.chip), Path(args.input), args.bin_size, Path(args.out), ddof=args.ddof)


if __name__ == "__main__":
    main()
