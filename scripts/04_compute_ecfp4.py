"""
Compute ECFP4 (Morgan radius=2) fingerprints for the aromatase bioactivity dataset.

Not included in the final output (replaced by PubChem FP) but script is available.
"""

import csv
import os
from rdkit import Chem
from rdkit.Chem import AllChem

_DIR = os.path.dirname(__file__)
INPUT_FILE = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_curated.csv")
OUTPUT_FILE = os.path.join(_DIR, "..", "data", "fingerprints", "aromatase_ecfp4_fp.csv")
NBITS = 1024
ECFP_COLS = [f"ECFP4_{i}" for i in range(NBITS)]


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def compute_ecfp4(smiles):
    """Compute ECFP4 (Morgan radius=2) fingerprint. Returns dict or None."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=NBITS)
    return {f"ECFP4_{i}": int(fp[i]) for i in range(NBITS)}


def main():
    print(f"Loading {INPUT_FILE}...")
    rows = load_csv(INPUT_FILE)
    print(f"  Loaded {len(rows)} rows")

    # Compute fingerprints for unique SMILES
    unique_smiles = set(r["canonical_smiles"] for r in rows if r["canonical_smiles"])
    print(f"\nComputing ECFP4 (radius=2, {NBITS} bits) for {len(unique_smiles)} unique SMILES...")

    fp_cache = {}
    failed = 0
    for i, smi in enumerate(unique_smiles):
        fp = compute_ecfp4(smi)
        if fp is not None:
            fp_cache[smi] = fp
        else:
            failed += 1
        if (i + 1) % 500 == 0:
            print(f"  Processed {i + 1}/{len(unique_smiles)}")

    print(f"  Computed: {len(fp_cache)}, Failed: {failed}")

    # Merge fingerprints into rows
    no_fp = 0
    for row in rows:
        smi = row["canonical_smiles"]
        fp = fp_cache.get(smi)
        if fp:
            row.update(fp)
        else:
            no_fp += 1
            row.update({c: "" for c in ECFP_COLS})

    # Write output
    original_cols = [c for c in rows[0].keys() if not c.startswith("ECFP4_")]
    fieldnames = original_cols + ECFP_COLS

    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    print(f"\n{'='*50}")
    print(f"Results saved to: {OUTPUT_FILE}")
    print(f"{'='*50}")
    maccs_cols = [c for c in original_cols if c.startswith("MACCS_")]
    base_cols = [c for c in original_cols if not c.startswith("MACCS_")]
    print(f"  Rows: {len(rows)}")
    print(f"  Columns: {len(fieldnames)} ({len(base_cols)} base + {len(maccs_cols)} MACCS + {len(ECFP_COLS)} ECFP4)")
    print(f"  Rows with ECFP4: {len(rows) - no_fp}")
    print(f"  Rows without ECFP4: {no_fp}")

    # Bit density stats
    on_counts = []
    for row in rows:
        bits_on = sum(1 for c in ECFP_COLS if row.get(c) == "1")
        on_counts.append(bits_on)

    avg_on = sum(on_counts) / len(on_counts)
    min_on = min(on_counts)
    max_on = max(on_counts)
    print(f"\n  ECFP4 bit density:")
    print(f"    Avg bits ON: {avg_on:.1f}/{NBITS}")
    print(f"    Min bits ON: {min_on}")
    print(f"    Max bits ON: {max_on}")


if __name__ == "__main__":
    main()
