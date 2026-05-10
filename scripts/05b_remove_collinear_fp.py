"""
Remove collinear fingerprint columns using Pearson correlation.

For each filtered fingerprint CSV in data/fingerprints_filtered/, computes
the pairwise Pearson correlation matrix and iteratively removes columns
involved in the most |r| >= 0.90 pairs (greedy heuristic, ties broken by
lower variance).

Writes reduced files to data/fingerprints_decorrelated/.

Usage:
    python scripts/05b_remove_collinear_fp.py                 # all fingerprints
    python scripts/05b_remove_collinear_fp.py maccs            # only MACCS
    python scripts/05b_remove_collinear_fp.py maccs pubchem    # MACCS + PubChem
"""

import csv
import os
import sys

import numpy as np

_DIR = os.path.dirname(__file__)
INPUT_DIR = os.path.join(_DIR, "..", "data", "fingerprints_filtered")
OUTPUT_DIR = os.path.join(_DIR, "..", "data", "fingerprints_decorrelated")
CORR_THRESHOLD = 0.90


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    return rows, fieldnames


def build_matrix(rows, fp_cols):
    """Build a numpy matrix from fingerprint columns."""
    n = len(rows)
    m = len(fp_cols)
    matrix = np.zeros((n, m), dtype=np.float64)
    for i, r in enumerate(rows):
        for j, c in enumerate(fp_cols):
            v = r[c]
            matrix[i, j] = float(v) if v != "" else 0.0
    return matrix


def find_collinear_pairs(corr_matrix, threshold):
    """Return set of (i, j) pairs where |r| >= threshold (i < j)."""
    n = corr_matrix.shape[0]
    abs_corr = np.abs(corr_matrix)
    # Zero out diagonal and lower triangle
    mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    pairs = set()
    rows_idx, cols_idx = np.where((abs_corr >= threshold) & mask)
    for i, j in zip(rows_idx, cols_idx):
        pairs.add((int(i), int(j)))
    return pairs


def greedy_remove(pairs, variances):
    """Greedily remove columns involved in the most correlated pairs.

    At each step, remove the column appearing in the most pairs.
    Ties broken by lower variance (remove the less variable column).
    """
    # Build a working copy of pairs
    active_pairs = set(pairs)
    removed = set()

    while active_pairs:
        # Count how many pairs each column is involved in
        counts = {}
        for i, j in active_pairs:
            counts[i] = counts.get(i, 0) + 1
            counts[j] = counts.get(j, 0) + 1

        # Find column with the most pairs; break ties by lower variance
        worst = max(counts, key=lambda c: (counts[c], -variances[c]))
        removed.add(worst)

        # Remove all pairs involving this column
        active_pairs = {(i, j) for i, j in active_pairs if i != worst and j != worst}

    return removed


def process_fingerprint(input_path, output_path):
    """Remove collinear columns from a single fingerprint CSV."""
    basename = os.path.basename(input_path)
    print(f"\n{'='*60}")
    print(f"Processing: {basename}")
    print(f"{'='*60}")

    rows, fieldnames = load_csv(input_path)
    fp_cols = [c for c in fieldnames if c != "molecule_chembl_id"]
    print(f"  Rows: {len(rows)}")
    print(f"  Fingerprint columns: {len(fp_cols)}")

    if len(fp_cols) < 2:
        print("  Skipping (fewer than 2 columns)")
        # Write as-is
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        return len(fp_cols), len(fp_cols), 0, 0

    # Build numeric matrix
    X = build_matrix(rows, fp_cols)

    # Compute Pearson correlation matrix
    # np.corrcoef returns NaN for zero-variance columns, but those were
    # already removed in step 05; handle gracefully anyway
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.corrcoef(X, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0)

    # Find collinear pairs
    pairs = find_collinear_pairs(corr, CORR_THRESHOLD)
    print(f"\n  Correlated pairs (|r| >= {CORR_THRESHOLD}): {len(pairs)}")

    if not pairs:
        print("  No collinear columns found — writing unchanged")
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        return len(fp_cols), len(fp_cols), 0, len(pairs)

    # Compute column variances for tie-breaking
    variances = np.var(X, axis=0)

    # Greedy removal
    removed_idx = greedy_remove(pairs, variances)
    removed_names = sorted([fp_cols[i] for i in removed_idx])
    kept_cols = [c for i, c in enumerate(fp_cols) if i not in removed_idx]

    print(f"  Columns removed: {len(removed_idx)}")
    print(f"  Columns remaining: {len(kept_cols)}")

    # Show a sample of removed columns (up to 20)
    if removed_names:
        show = removed_names[:20]
        print(f"\n  Removed columns (showing {len(show)}/{len(removed_names)}):")
        for c in show:
            idx = fp_cols.index(c)
            # Find its highest correlation with a kept column
            abs_row = np.abs(corr[idx])
            abs_row[idx] = 0.0  # exclude self
            max_partner = np.argmax(abs_row)
            max_r = abs_row[max_partner]
            print(f"    {c}: max |r|={max_r:.4f} with {fp_cols[max_partner]}")

    # Write filtered output
    out_fieldnames = ["molecule_chembl_id"] + kept_cols
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  Wrote: {output_path}")
    print(f"  Dimensions: {len(rows)} rows x {len(out_fieldnames)} columns")

    return len(fp_cols), len(kept_cols), len(removed_idx), len(pairs)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Collect filtered fingerprint files
    all_fps = sorted(
        f for f in os.listdir(INPUT_DIR) if f.endswith("_filtered.csv")
    )

    # Filter by name if arguments given (e.g. "maccs", "pubchem")
    args = [a.lower() for a in sys.argv[1:]]
    if args:
        fps = [f for f in all_fps if any(a in f.lower() for a in args)]
    else:
        fps = all_fps

    if not fps:
        print(f"No filtered fingerprint files found matching: {args or 'all'}")
        print(f"Available: {all_fps}")
        return

    print(f"Fingerprint files to process: {len(fps)}")
    print(f"Correlation threshold: |r| >= {CORR_THRESHOLD}")
    for f in fps:
        print(f"  - {f}")

    # Process each file
    total_original = 0
    total_kept = 0
    total_removed = 0
    total_pairs = 0

    results = []
    for f in fps:
        input_path = os.path.join(INPUT_DIR, f)
        output_name = f.replace("_filtered.csv", "_decorrelated.csv")
        output_path = os.path.join(OUTPUT_DIR, output_name)
        n_orig, n_kept, n_removed, n_pairs = process_fingerprint(input_path, output_path)
        total_original += n_orig
        total_kept += n_kept
        total_removed += n_removed
        total_pairs += n_pairs
        results.append({
            "name": f,
            "n_orig": n_orig,
            "n_kept": n_kept,
            "n_removed": n_removed,
            "n_pairs": n_pairs,
        })

    # Grand summary
    print(f"\n{'='*60}")
    print(f"Grand Summary")
    print(f"{'='*60}")
    print(f"  Files processed: {len(fps)}")
    print(f"  Correlation threshold: |r| >= {CORR_THRESHOLD}")
    print(f"  Total correlated pairs found: {total_pairs}")
    print(f"  Total original columns: {total_original}")
    print(f"  Total removed: {total_removed}")
    print(f"  Total remaining: {total_kept}")
    print(f"  Output directory: {OUTPUT_DIR}")

    print(f"\n  {'Fingerprint':<45} {'Original':>8} {'Removed':>8} {'Kept':>8} {'Pairs':>8}")
    print(f"  {'-'*45} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
    for r in results:
        name = r["name"].replace("aromatase_", "").replace("_fp_filtered.csv", "")
        print(f"  {name:<45} {r['n_orig']:>8} {r['n_removed']:>8} {r['n_kept']:>8} {r['n_pairs']:>8}")


if __name__ == "__main__":
    main()
