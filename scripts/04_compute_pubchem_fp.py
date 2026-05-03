"""
Compute PubChem fingerprints (881 bits) using PaDEL batch mode.

Reads aromatase_bioactivity_fp.csv (has MACCS + ECFP4), computes PubChem
fingerprints via PaDEL batch processing with SMILES input, replaces ECFP4
columns with 881 PubChem FP columns, and overwrites aromatase_bioactivity_fp.csv.

Note: SMILES (.smi) input is used instead of SDF to avoid aromaticity
perception differences that cause ring-count bit mismatches in Section 2
(bits 115-262) of the PubChem fingerprint.
"""

import csv
import os
import tempfile
import time

from padelpy import padeldescriptor

_DIR = os.path.dirname(__file__)
INPUT_FILE = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_curated.csv")
OUTPUT_FILE = os.path.join(_DIR, "..", "data", "fingerprints", "aromatase_pubchem_fp.csv")
PUBCHEM_COLS = [f"PubchemFP{i}" for i in range(881)]


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_smi(smiles_map, smi_path):
    """Write unique molecules to SMILES (.smi) file. smiles_map: {chembl_id: smiles}."""
    written = 0
    skipped = []
    with open(smi_path, "w", encoding="utf-8") as f:
        for cid, smi in smiles_map.items():
            if not smi or not smi.strip():
                skipped.append(cid)
                continue
            f.write(f"{smi}\t{cid}\n")
            written += 1
    return written, skipped


def run_padel_batch(smi_path, output_csv):
    """Run PaDEL fingerprint calculation in batch mode with SMILES input."""
    padeldescriptor(
        mol_dir=smi_path,
        d_file=output_csv,
        fingerprints=True,
        removesalt=True,
        retainorder=True,
        threads=1,
    )


def main():
    print(f"Loading {INPUT_FILE}...")
    rows = load_csv(INPUT_FILE)
    print(f"  Loaded {len(rows)} rows")

    # Collect unique molecules
    smiles_map = {}
    for r in rows:
        cid = r["molecule_chembl_id"]
        smi = r["canonical_smiles"]
        if cid not in smiles_map and smi:
            smiles_map[cid] = smi
    print(f"  Unique molecules: {len(smiles_map)}")

    # Create temp files
    tmpdir = tempfile.mkdtemp()
    smi_path = os.path.join(tmpdir, "molecules.smi")
    fp_csv_path = os.path.join(tmpdir, "fingerprints.csv")

    # Write SMILES file
    print(f"\nWriting SMILES file...")
    written, skipped = write_smi(smiles_map, smi_path)
    print(f"  Written: {written}, Skipped: {len(skipped)}")
    if skipped:
        print(f"  Skipped IDs: {skipped[:10]}")

    # Run PaDEL batch
    print(f"\nRunning PaDEL batch fingerprint calculation (SMILES input)...")
    start = time.time()
    run_padel_batch(smi_path, fp_csv_path)
    elapsed = time.time() - start
    print(f"  PaDEL completed in {elapsed:.1f}s ({elapsed/written*1000:.0f}ms per molecule)")

    # Load PaDEL results
    with open(fp_csv_path, newline="", encoding="utf-8") as f:
        fp_rows = list(csv.DictReader(f))
    print(f"  PaDEL output rows: {len(fp_rows)}")

    # Build lookup: chembl_id -> PubChem FP dict
    fp_map = {}
    for fr in fp_rows:
        name = fr.get("Name", "")
        fp_dict = {c: fr.get(c, "0") for c in PUBCHEM_COLS}
        fp_map[name] = fp_dict

    print(f"  Mapped {len(fp_map)} molecules to PubChem FPs")

    # Build output: molecule_chembl_id + PubChem FP columns only
    out_fieldnames = ["molecule_chembl_id"] + PUBCHEM_COLS

    missing_fp = 0
    for row in rows:
        cid = row["molecule_chembl_id"]
        fp = fp_map.get(cid)
        if fp:
            row.update(fp)
        else:
            missing_fp += 1
            row.update({c: "" for c in PUBCHEM_COLS})

    # Write output
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Cleanup temp files
    os.remove(smi_path)
    os.remove(fp_csv_path)
    os.rmdir(tmpdir)

    # Summary
    print(f"\n{'='*55}")
    print(f"Results saved to: {OUTPUT_FILE}")
    print(f"{'='*55}")
    print(f"  Rows: {len(rows)}")
    print(f"  Columns: {len(out_fieldnames)} (1 ID + {len(PUBCHEM_COLS)} PubChem)")
    print(f"  Rows with PubChem FP: {len(rows) - missing_fp}")
    print(f"  Rows without PubChem FP: {missing_fp}")

    # Bit density
    on_counts = []
    for row in rows:
        bits_on = sum(1 for c in PUBCHEM_COLS if row.get(c) == "1")
        on_counts.append(bits_on)

    avg_on = sum(on_counts) / len(on_counts) if on_counts else 0
    print(f"\n  PubChem FP bit density:")
    print(f"    Avg bits ON: {avg_on:.1f}/{len(PUBCHEM_COLS)}")
    print(f"    Min bits ON: {min(on_counts)}")
    print(f"    Max bits ON: {max(on_counts)}")

    # Top 10 most common bits
    bit_sums = {c: 0 for c in PUBCHEM_COLS}
    counted = len(rows) - missing_fp
    for row in rows:
        for c in PUBCHEM_COLS:
            if row.get(c) == "1":
                bit_sums[c] += 1

    print(f"\n  Top 10 most common PubChem bits:")
    sorted_bits = sorted(bit_sums.items(), key=lambda x: x[1], reverse=True)
    for name, count in sorted_bits[:10]:
        pct = count / counted * 100 if counted else 0
        print(f"    {name}: {count}/{counted} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
