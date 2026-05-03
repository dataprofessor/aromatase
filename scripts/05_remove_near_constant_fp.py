"""
Remove near-constant fingerprint columns using a standard deviation (SD) threshold.

For each fingerprint CSV in data/fingerprints/, computes the SD of every
numeric column and drops columns with SD < 0.1.  Writes filtered files
to data/fingerprints_filtered/ with the same base name.

Usage:
    python scripts/05_remove_near_constant_fp.py                 # all fingerprints
    python scripts/05_remove_near_constant_fp.py maccs            # only MACCS
    python scripts/05_remove_near_constant_fp.py maccs pubchem    # MACCS + PubChem
"""

import csv
import math
import os
import sys

_DIR = os.path.dirname(__file__)
INPUT_DIR = os.path.join(_DIR, "..", "data", "fingerprints")
OUTPUT_DIR = os.path.join(_DIR, "..", "data", "fingerprints_filtered")
SD_THRESHOLD = 0.1


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def detect_fp_columns(fieldnames):
    """Return fingerprint columns (everything except molecule_chembl_id)."""
    return [c for c in fieldnames if c != "molecule_chembl_id"]


def compute_sd(rows, col):
    """Compute population SD for a column, skipping empty values."""
    vals = []
    for r in rows:
        v = r[col]
        if v != "":
            vals.append(float(v))
    if len(vals) < 2:
        return 0.0
    mean = sum(vals) / len(vals)
    variance = sum((v - mean) ** 2 for v in vals) / len(vals)
    return math.sqrt(variance)


def filter_fingerprint(input_path, output_path):
    """Filter near-constant columns from a fingerprint CSV."""
    basename = os.path.basename(input_path)
    print(f"\n{'='*60}")
    print(f"Processing: {basename}")
    print(f"{'='*60}")

    rows = load_csv(input_path)
    with open(input_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

    fp_cols = detect_fp_columns(fieldnames)
    print(f"  Rows: {len(rows)}")
    print(f"  Fingerprint columns: {len(fp_cols)}")

    # Compute SD for each column
    col_sd = {}
    for c in fp_cols:
        col_sd[c] = compute_sd(rows, c)

    kept = [c for c in fp_cols if col_sd[c] >= SD_THRESHOLD]
    removed = [c for c in fp_cols if col_sd[c] < SD_THRESHOLD]

    print(f"\n  SD threshold: {SD_THRESHOLD}")
    print(f"  Removed (SD < {SD_THRESHOLD}): {len(removed)}")
    print(f"  Remaining: {len(kept)}")

    if removed:
        print(f"\n  Removed columns:")
        for c in removed:
            vals = [float(r[c]) for r in rows if r[c] != ""]
            n = len(vals)
            freq = (sum(vals) / n * 100) if n else 0
            print(f"    {c}: SD={col_sd[c]:.4f}, freq={freq:.1f}%")

    # Write filtered output
    out_fieldnames = ["molecule_chembl_id"] + kept
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  Wrote: {output_path}")
    print(f"  Dimensions: {len(rows)} rows x {len(out_fieldnames)} columns")

    return len(fp_cols), len(kept), len(removed)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Collect fingerprint files
    all_fps = sorted(
        f for f in os.listdir(INPUT_DIR) if f.endswith("_fp.csv")
    )

    # Filter by name if arguments given (e.g. "maccs", "pubchem")
    args = [a.lower() for a in sys.argv[1:]]
    if args:
        fps = [f for f in all_fps if any(a in f.lower() for a in args)]
    else:
        fps = all_fps

    if not fps:
        print(f"No fingerprint files found matching: {args or 'all'}")
        print(f"Available: {all_fps}")
        return

    print(f"Fingerprint files to process: {len(fps)}")
    for f in fps:
        print(f"  - {f}")

    # Process each file
    total_original = 0
    total_kept = 0
    total_removed = 0

    for f in fps:
        input_path = os.path.join(INPUT_DIR, f)
        output_path = os.path.join(OUTPUT_DIR, f.replace("_fp.csv", "_fp_filtered.csv"))
        n_orig, n_kept, n_removed = filter_fingerprint(input_path, output_path)
        total_original += n_orig
        total_kept += n_kept
        total_removed += n_removed

    # Grand summary
    print(f"\n{'='*60}")
    print(f"Grand Summary")
    print(f"{'='*60}")
    print(f"  Files processed: {len(fps)}")
    print(f"  Total original columns: {total_original}")
    print(f"  Total removed (SD < {SD_THRESHOLD}): {total_removed}")
    print(f"  Total remaining: {total_kept}")
    print(f"  Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
