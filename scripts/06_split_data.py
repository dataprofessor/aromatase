"""
Split filtered fingerprint data into train/test sets (80/20).

Two splitting strategies, both operating at the molecule level to prevent
data leakage (all rows for a given molecule go to the same set):

  1. Random split — reproducible shuffle with fixed seed
  2. Kennard-Stone split — maximises chemical diversity in the training set
     by iteratively selecting the most distant molecule from the current
     training pool

Reads from data/fingerprints_filtered/ and writes to data/splits/.

Usage:
    python scripts/06_split_data.py                 # all fingerprints
    python scripts/06_split_data.py maccs            # only MACCS
    python scripts/06_split_data.py maccs pubchem    # MACCS + PubChem
"""

import csv
import os
import random
import sys

import numpy as np

_DIR = os.path.dirname(__file__)
INPUT_DIR = os.path.join(_DIR, "..", "data", "fingerprints_filtered")
OUTPUT_DIR = os.path.join(_DIR, "..", "data", "splits")
TRAIN_RATIO = 0.8
RANDOM_SEED = 42


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    return rows, fieldnames


def get_unique_molecules(rows):
    """Return ordered list of unique molecule_chembl_id values."""
    seen = set()
    ordered = []
    for r in rows:
        mid = r["molecule_chembl_id"]
        if mid not in seen:
            seen.add(mid)
            ordered.append(mid)
    return ordered


def build_fp_matrix(rows, fp_cols, mol_ids):
    """Build a numpy matrix of fingerprints for unique molecules.

    For molecules with multiple rows, uses the first row's fingerprint
    (fingerprints are identical across assay types for the same molecule).
    """
    mol_to_row = {}
    for r in rows:
        mid = r["molecule_chembl_id"]
        if mid not in mol_to_row:
            mol_to_row[mid] = r

    matrix = np.zeros((len(mol_ids), len(fp_cols)), dtype=np.float64)
    for i, mid in enumerate(mol_ids):
        row = mol_to_row[mid]
        for j, c in enumerate(fp_cols):
            v = row[c]
            matrix[i, j] = float(v) if v != "" else 0.0
    return matrix


def kennard_stone(X, n_train):
    """Kennard-Stone algorithm for selecting a diverse training set.

    Parameters
    ----------
    X : np.ndarray, shape (n_samples, n_features)
    n_train : int, number of samples to select for training

    Returns
    -------
    train_idx : list of int, indices selected for training
    test_idx : list of int, remaining indices for testing
    """
    n = X.shape[0]

    # Precompute pairwise Euclidean distances
    # Using (a-b)^2 = a^2 + b^2 - 2ab for efficiency
    sq_norms = np.sum(X ** 2, axis=1)
    dist_matrix = np.sqrt(
        np.maximum(sq_norms[:, None] + sq_norms[None, :] - 2.0 * X @ X.T, 0.0)
    )

    # Start with the two most distant samples
    i, j = np.unravel_index(np.argmax(dist_matrix), dist_matrix.shape)
    selected = [int(i), int(j)]
    remaining = set(range(n)) - {i, j}

    # min_dist_to_selected[k] = min distance from k to any selected sample
    min_dist = np.minimum(dist_matrix[i], dist_matrix[j])

    while len(selected) < n_train:
        # Among remaining, pick the one with the largest min-distance
        candidates = np.array(list(remaining))
        best_local = np.argmax(min_dist[candidates])
        best_idx = candidates[best_local]

        selected.append(int(best_idx))
        remaining.discard(int(best_idx))

        # Update min distances
        min_dist = np.minimum(min_dist, dist_matrix[best_idx])

    test_idx = sorted(remaining)
    return selected, test_idx


def random_split(mol_ids, n_train, seed):
    """Random split of molecule IDs."""
    ids = list(mol_ids)
    rng = random.Random(seed)
    rng.shuffle(ids)
    return ids[:n_train], ids[n_train:]


def write_split(rows, fieldnames, mol_set, output_path):
    """Write rows whose molecule_chembl_id is in mol_set."""
    mol_set = set(mol_set)
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            if r["molecule_chembl_id"] in mol_set:
                writer.writerow(r)


def process_fingerprint(input_path, output_dir):
    """Run both split strategies on a single fingerprint file."""
    basename = os.path.basename(input_path)
    stem = basename.replace("_filtered.csv", "")
    print(f"\n{'='*60}")
    print(f"Processing: {basename}")
    print(f"{'='*60}")

    rows, fieldnames = load_csv(input_path)
    fp_cols = [c for c in fieldnames if c != "molecule_chembl_id"]
    mol_ids = get_unique_molecules(rows)

    n_molecules = len(mol_ids)
    n_train = int(n_molecules * TRAIN_RATIO)
    n_test = n_molecules - n_train

    print(f"  Total rows: {len(rows)}")
    print(f"  Unique molecules: {n_molecules}")
    print(f"  FP columns: {len(fp_cols)}")
    print(f"  Train molecules: {n_train} ({n_train/n_molecules*100:.1f}%)")
    print(f"  Test molecules: {n_test} ({n_test/n_molecules*100:.1f}%)")

    # ── Random split ──
    print(f"\n  [1/2] Random split (seed={RANDOM_SEED})...")
    train_mols_rand, test_mols_rand = random_split(mol_ids, n_train, RANDOM_SEED)

    train_path = os.path.join(output_dir, f"{stem}_random_train.csv")
    test_path = os.path.join(output_dir, f"{stem}_random_test.csv")
    write_split(rows, fieldnames, train_mols_rand, train_path)
    write_split(rows, fieldnames, test_mols_rand, test_path)

    train_rows = sum(1 for r in rows if r["molecule_chembl_id"] in set(train_mols_rand))
    test_rows = sum(1 for r in rows if r["molecule_chembl_id"] in set(test_mols_rand))
    print(f"    Train: {train_rows} rows ({len(train_mols_rand)} molecules)")
    print(f"    Test:  {test_rows} rows ({len(test_mols_rand)} molecules)")

    # ── Kennard-Stone split ──
    print(f"\n  [2/2] Kennard-Stone split...")
    X = build_fp_matrix(rows, fp_cols, mol_ids)
    train_idx, test_idx = kennard_stone(X, n_train)

    train_mols_ks = [mol_ids[i] for i in train_idx]
    test_mols_ks = [mol_ids[i] for i in test_idx]

    train_path = os.path.join(output_dir, f"{stem}_kennard_stone_train.csv")
    test_path = os.path.join(output_dir, f"{stem}_kennard_stone_test.csv")
    write_split(rows, fieldnames, train_mols_ks, train_path)
    write_split(rows, fieldnames, test_mols_ks, test_path)

    train_rows_ks = sum(1 for r in rows if r["molecule_chembl_id"] in set(train_mols_ks))
    test_rows_ks = sum(1 for r in rows if r["molecule_chembl_id"] in set(test_mols_ks))
    print(f"    Train: {train_rows_ks} rows ({len(train_mols_ks)} molecules)")
    print(f"    Test:  {test_rows_ks} rows ({len(test_mols_ks)} molecules)")

    return {
        "name": basename,
        "n_rows": len(rows),
        "n_molecules": n_molecules,
        "n_train_mol": n_train,
        "n_test_mol": n_test,
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Collect filtered fingerprint files
    all_fps = sorted(
        f for f in os.listdir(INPUT_DIR) if f.endswith("_filtered.csv")
    )

    # Filter by name if arguments given
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
    print(f"Split ratio: {TRAIN_RATIO:.0%} train / {1-TRAIN_RATIO:.0%} test")
    print(f"Methods: random (seed={RANDOM_SEED}), Kennard-Stone")
    for f in fps:
        print(f"  - {f}")

    results = []
    for f in fps:
        input_path = os.path.join(INPUT_DIR, f)
        info = process_fingerprint(input_path, OUTPUT_DIR)
        results.append(info)

    # Grand summary
    print(f"\n{'='*60}")
    print(f"Grand Summary")
    print(f"{'='*60}")
    print(f"  Files processed: {len(results)}")
    print(f"  Output directory: {OUTPUT_DIR}")
    print(f"  Files written: {len(results) * 4} (2 methods x train/test)")
    print(f"\n  {'Fingerprint':<50} {'Molecules':>10} {'Train':>7} {'Test':>6}")
    print(f"  {'-'*50} {'-'*10} {'-'*7} {'-'*6}")
    for r in results:
        name = r["name"].replace("aromatase_", "").replace("_fp_filtered.csv", "")
        print(f"  {name:<50} {r['n_molecules']:>10} {r['n_train_mol']:>7} {r['n_test_mol']:>6}")


if __name__ == "__main__":
    main()
