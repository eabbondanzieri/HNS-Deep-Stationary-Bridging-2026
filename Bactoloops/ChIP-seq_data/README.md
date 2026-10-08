# chip_zscore.py

Computes a ChIP-seq (or any IP) enrichment z-score track from two per-base
coverage files: a ChIP/IP track and a WGS/input track used to normalize it.

This is the pipeline validated against the Hi-C loop-calling macro
(`HiC_PeakFindingAlgorithm_with_symmetrization.ijm`): bin both tracks, divide
ChIP by input bin-for-bin, then z-score the resulting ratio. The script
reproduces that earlier calculation exactly and generalizes it into a
reusable, documented tool with single-file and batch modes.

## Requirements

- Python 3
- `numpy`, `pandas`

```
pip install numpy pandas
```

## Pipeline

1. **Bin** both tracks to `--bin-size` bp, taking the mean coverage per bin.
   Bin 1 covers positions 1..bin_size, bin 2 covers bin_size+1..2*bin_size,
   and so on. The last bin may be partial if the genome length isn't an
   exact multiple of the bin size.
2. **Ratio** = ChIP bin mean / Input bin mean, bin by bin.
3. **Z-score** = (ratio − mean(ratio)) / std(ratio), computed over all finite
   bins. A bin with zero input coverage gives a non-finite ratio; it's
   excluded from the mean/std and left as `NaN` in the output rather than
   silently dropped.

## Input format

A two-column CSV with a header — `Position,Value` (any two column names
work; the first two columns are used positionally). Position is 1-based and
contiguous, one row per base pair — the format FIJI/ImageJ's "Plot
Profile"/coverage export and similar tools produce.

## Usage

### Single pair

```
python3 chip_zscore.py --chip anti-HNS_sample.csv --input WGS_sample.csv \
    --bin-size 10000 --out sample_zscore_10kb.csv
```

Requires `--chip`, `--input`, `--out`, and `--bin-size`.

### Batch — every ChIP/Input pair in a folder

Pairing is done by filename: a ChIP file `<chip-prefix><suffix>` is paired
with the Input file `<input-prefix><suffix>` in the same folder (same
`<suffix>`) — e.g. `anti-HNS_cWX2200_96hr_rep2_R1_clean.csv` pairs with
`WGS_cWX2200_96hr_rep2_R1_clean.csv` by default.

```
python3 chip_zscore.py --batch-dir /path/to/folder \
    --chip-prefix "anti-HNS_" --input-prefix "WGS_" \
    --bin-size 10000 --out-dir /path/to/folder/zscores
```

ChIP files with no matching Input file are reported and skipped rather than
failing the run.

## Options

| Flag | Mode | Description |
|---|---|---|
| `--chip` | single | ChIP/IP per-base coverage CSV |
| `--input` | single | WGS/input per-base coverage CSV, normalizes `--chip` |
| `--out` | single | Output CSV path |
| `--batch-dir` | batch | Folder to scan for ChIP/Input pairs |
| `--chip-prefix` | batch | Filename prefix identifying a ChIP file (default `anti-HNS_`) |
| `--input-prefix` | batch | Filename prefix identifying the matching Input file (default `WGS_`) |
| `--out-dir` | batch | Output folder (default: `--batch-dir` itself) |
| `--out-suffix` | batch | Inserted before resolution/extension in output filenames (default `_zscore`) |
| `--bin-size` | both | Bin size in bp, e.g. `10000` (required) |
| `--ddof` | both | Delta degrees of freedom for the z-score std (default `1`, sample std) |

Single-pair mode is triggered by any of `--chip` / `--input` / `--out`;
batch mode is triggered by `--batch-dir`. The two modes are mutually
exclusive.

## Output columns

| Column | Meaning |
|---|---|
| `bin` | 1-based bin number |
| `start_bp` / `end_bp` | Bin boundaries in bp |
| `chip_mean` / `input_mean` | Mean coverage in the bin, for each track |
| `ratio` | `chip_mean / input_mean` |
| `zscore` | Z-scored ratio, `NaN` for non-finite-ratio bins |

## Notes

- If the ChIP and Input files don't cover the same Position range, a warning
  is printed to stderr but the run proceeds — useful to double-check the
  files really are the same sample/genome.
- Batch-mode output filenames follow
  `{chip_prefix}{suffix_stem}{out_suffix}_{bin_size}bp.csv`, e.g.
  `anti-HNS_cWX2200_96hr_rep2_R1_clean_zscore_10000bp.csv`.
