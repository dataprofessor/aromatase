"""
Compute MACCS fingerprints (166 bits) for the aromatase bioactivity dataset.

Reads aromatase_bioactivity_clean.csv, computes MACCS keys using RDKit
for each unique molecule, and writes:
  - data/processed/aromatase_bioactivity_curated.csv (base columns only)
  - data/fingerprints/aromatase_maccs_fp.csv (molecule_chembl_id + 166 MACCS bits)
"""

import csv
import os
from rdkit import Chem
from rdkit.Chem import MACCSkeys

_DIR = os.path.dirname(__file__)
INPUT_FILE = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_clean.csv")
OUTPUT_CURATED = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_curated.csv")
OUTPUT_MACCS = os.path.join(_DIR, "..", "data", "fingerprints", "aromatase_maccs_fp.csv")
MACCS_COLS = [f"MACCS_{i}" for i in range(1, 167)]  # bits 1-166 (bit 0 is always 0)


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def compute_maccs(smiles):
    """Compute MACCS keys for a SMILES string. Returns dict of bit values or None."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = MACCSkeys.GenMACCSKeys(mol)
    return {f"MACCS_{i}": int(fp[i]) for i in range(1, 167)}


def main():
    print(f"Loading {INPUT_FILE}...")
    rows = load_csv(INPUT_FILE)
    print(f"  Loaded {len(rows)} rows")

    # Compute fingerprints for unique SMILES
    unique_smiles = set(r["canonical_smiles"] for r in rows if r["canonical_smiles"])
    print(f"\nComputing MACCS keys for {len(unique_smiles)} unique SMILES...")

    fp_cache = {}
    failed = 0
    for i, smi in enumerate(unique_smiles):
        fp = compute_maccs(smi)
        if fp is not None:
            fp_cache[smi] = fp
        else:
            failed += 1
        if (i + 1) % 500 == 0:
            print(f"  Processed {i + 1}/{len(unique_smiles)}")

    print(f"  Computed: {len(fp_cache)}, Failed: {failed}")

    # Merge fingerprints into rows
    out_rows = []
    no_fp = 0
    for row in rows:
        smi = row["canonical_smiles"]
        fp = fp_cache.get(smi)
        if fp:
            row.update(fp)
        else:
            no_fp += 1
            row.update({c: "" for c in MACCS_COLS})
        out_rows.append(row)

    # Assign bioactivity class based on pchembl_value
    #   active:       pchembl_value > 7
    #   inactive:     pchembl_value < 6
    #   intermediate: 6 <= pchembl_value <= 7
    class_counts = {"active": 0, "inactive": 0, "intermediate": 0, "": 0}
    for row in out_rows:
        pval = row.get("pchembl_value", "")
        if pval == "":
            row["bioactivity_class"] = ""
            class_counts[""] += 1
        else:
            pval = float(pval)
            if pval > 7:
                row["bioactivity_class"] = "active"
                class_counts["active"] += 1
            elif pval < 6:
                row["bioactivity_class"] = "inactive"
                class_counts["inactive"] += 1
            else:
                row["bioactivity_class"] = "intermediate"
                class_counts["intermediate"] += 1

    print(f"\n  Bioactivity class distribution:")
    print(f"    active (pChEMBL > 7):       {class_counts['active']}")
    print(f"    intermediate (6-7):         {class_counts['intermediate']}")
    print(f"    inactive (pChEMBL < 6):     {class_counts['inactive']}")
    print(f"    unassigned (no pChEMBL):    {class_counts['']}")

    # Write curated bioactivity file (no fingerprints)
    original_cols = [c for c in rows[0].keys() if c not in MACCS_COLS]

    with open(OUTPUT_CURATED, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=original_cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"\n  Wrote {OUTPUT_CURATED}: {len(out_rows)} rows x {len(original_cols)} columns")

    # Write MACCS fingerprint file (molecule_chembl_id + MACCS bits)
    maccs_fieldnames = ["molecule_chembl_id"] + MACCS_COLS

    with open(OUTPUT_MACCS, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=maccs_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    # Summary
    print(f"  Wrote {OUTPUT_MACCS}: {len(out_rows)} rows x {len(maccs_fieldnames)} columns")
    print(f"\n{'='*50}")
    print(f"Summary")
    print(f"{'='*50}")
    print(f"  Rows: {len(out_rows)}")
    print(f"  Rows with fingerprints: {len(out_rows) - no_fp}")
    print(f"  Rows without fingerprints: {no_fp}")

    # Bit frequency summary
    bit_sums = {c: 0 for c in MACCS_COLS}
    counted = 0
    for row in out_rows:
        if row.get("MACCS_1") != "":
            counted += 1
            for c in MACCS_COLS:
                bit_sums[c] += int(row[c])

    print(f"\n  Top 10 most common MACCS bits:")
    sorted_bits = sorted(bit_sums.items(), key=lambda x: x[1], reverse=True)
    for name, count in sorted_bits[:10]:
        pct = count / counted * 100 if counted else 0
        print(f"    {name}: {count}/{counted} ({pct:.1f}%)")

    print(f"\n  Bottom 5 rarest MACCS bits (non-zero):")
    nonzero_bits = [(n, c) for n, c in sorted_bits if c > 0]
    for name, count in nonzero_bits[-5:]:
        pct = count / counted * 100 if counted else 0
        print(f"    {name}: {count}/{counted} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
