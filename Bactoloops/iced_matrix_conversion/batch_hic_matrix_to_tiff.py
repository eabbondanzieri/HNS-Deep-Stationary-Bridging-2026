#!/usr/bin/env python3
"""
batch_hic_matrix_to_tiff.py

Batch wrapper around hic_circular_recenter_to_tiff.py: converts every
matching sparse Hi-C matrix file in a directory (default: the current
directory) to a 32-bit float TIFF, using the same --bin-size/--id-mode/
--log1p/--recenter-bin settings for all of them.

Must live in the same folder as hic_circular_recenter_to_tiff.py -- it
imports build_and_convert() from it rather than duplicating the logic, so
any future fix to the single-file script is picked up here automatically.

USAGE
-----
    # Convert every *.matrix file in the current directory, same as before:
    python3 batch_hic_matrix_to_tiff.py --bin-size 10000

    # Different input/output locations, a custom pattern, and recentering
    # applied to every file:
    python3 batch_hic_matrix_to_tiff.py --bin-size 10000 \\
        --dir ./matrices --out-dir ./tiffs --pattern "*_iced.matrix" \\
        --recenter-bin 368

    Other optional flags: [--log1p] [--id-offset 0] [--id-mode auto|bp|index]
                          [--max-pixels 4e7] [--overwrite]

By default, existing output files are left alone (skipped) so re-running
the script after adding a few new .matrix files doesn't redo everything --
pass --overwrite to force reconversion. One file failing (e.g. a corrupt
or unrelated .matrix file) is reported but does not stop the rest of the
batch; a summary of successes/skips/failures is printed at the end.
"""

import argparse
import sys
import traceback
from pathlib import Path

try:
    from hic_circular_recenter_to_tiff import build_and_convert
except ImportError:
    sys.exit(
        "Could not import build_and_convert() from hic_circular_recenter_to_tiff.py.\n"
        "Make sure batch_hic_matrix_to_tiff.py is in the SAME folder as "
        "hic_circular_recenter_to_tiff.py, and that you're running this from "
        "that folder (or that folder is on your PYTHONPATH)."
    )


def main():
    p = argparse.ArgumentParser(
        description="Convert every matching Hi-C sparse matrix file in a "
                    "directory to a 32-bit float TIFF.")
    p.add_argument("--bin-size", type=int, required=True, help="Resolution in bp")
    p.add_argument("--dir", default=".", help="Directory to scan for matrix files (default: current directory)")
    p.add_argument("--out-dir", default=None, help="Directory to write TIFFs into (default: same as --dir)")
    p.add_argument("--pattern", default="*.matrix", help="Glob pattern for input files (default: *.matrix)")
    p.add_argument("--id-mode", choices=["auto", "bp", "index"], default="auto")
    p.add_argument("--id-offset", type=int, default=0)
    p.add_argument("--max-pixels", type=float, default=4e7)
    p.add_argument("--log1p", action="store_true")
    p.add_argument("--recenter-bin", type=int, default=None)
    p.add_argument("--overwrite", action="store_true",
                    help="Reconvert and overwrite files that already have a matching .tif output")
    args = p.parse_args()

    in_dir = Path(args.dir)
    out_dir = Path(args.out_dir) if args.out_dir else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_files = sorted(in_dir.glob(args.pattern))
    if not matrix_files:
        sys.exit(f"No files matching '{args.pattern}' found in {in_dir.resolve()}")

    print(f"Found {len(matrix_files)} file(s) matching '{args.pattern}' in {in_dir.resolve()}")

    converted, skipped, failed = [], [], []

    for i, matrix_file in enumerate(matrix_files, start=1):
        out_file = out_dir / (matrix_file.stem + ".tif")
        print(f"\n[{i}/{len(matrix_files)}] {matrix_file.name}")

        if out_file.exists() and not args.overwrite:
            print(f"  Skipping: {out_file.name} already exists (use --overwrite to redo it)")
            skipped.append(matrix_file.name)
            continue

        try:
            build_and_convert(
                matrix_file=str(matrix_file),
                bin_size=args.bin_size,
                out_file=str(out_file),
                id_mode=args.id_mode,
                log1p_transform=args.log1p,
                id_offset=args.id_offset,
                max_pixels=args.max_pixels,
                recenter_bin=args.recenter_bin,
            )
            converted.append(matrix_file.name)
        except SystemExit as e:
            # build_and_convert() uses sys.exit(...) on data problems (e.g. a
            # bad --recenter-bin, or bin ids that don't line up) -- catch that
            # here so one bad file doesn't kill the whole batch.
            print(f"  FAILED: {e}", file=sys.stderr)
            failed.append((matrix_file.name, str(e)))
        except Exception as e:
            print(f"  FAILED: {e}", file=sys.stderr)
            traceback.print_exc()
            failed.append((matrix_file.name, str(e)))

    print("\n" + "=" * 60)
    print(f"Done: {len(converted)} converted, {len(skipped)} skipped, {len(failed)} failed")
    if failed:
        print("\nFailed files:")
        for name, err in failed:
            print(f"  - {name}: {err}")
        sys.exit(1)


if __name__ == "__main__":
    main()
