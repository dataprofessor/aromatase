"""
Compute AtomPairs2D Count fingerprints (780 counts) using PaDEL batch mode with SMILES input.

Uses PaDEL-Descriptor's AtomPairs2DFingerprintCount. These are frequency counts (not binary).
Column names are descriptive (e.g., APC2D1_C_C, APC2D1_C_N) rather than numeric indices.
SMILES (.smi) input is used instead of SDF to ensure consistent aromaticity perception.
"""

import csv
import os
import tempfile
import time

from padelpy import padeldescriptor

_DIR = os.path.dirname(__file__)
INPUT_FILE = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_curated.csv")
OUTPUT_FILE = os.path.join(_DIR, "..", "data", "fingerprints", "aromatase_atompairs2d_count_fp.csv")

NBITS = 780


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_smi(smiles_map, smi_path):
    """Write unique molecules to SMILES (.smi) file."""
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


def make_descriptor_xml(output_path):
    """Create a custom descriptors.xml enabling only AtomPairs2DFingerprintCount."""
    import importlib
    padel_dir = os.path.dirname(importlib.import_module("padelpy").__file__)
    src_xml = os.path.join(padel_dir, "PaDEL-Descriptor", "descriptors.xml")
    with open(src_xml, encoding="utf-8") as f:
        xml = f.read()
    xml = xml.replace(
        '<Descriptor name="AtomPairs2DFingerprintCount" value="false"/>',
        '<Descriptor name="AtomPairs2DFingerprintCount" value="true"/>',
    )
    xml = xml.replace(
        '<Descriptor name="PubchemFingerprinter" value="true"/>',
        '<Descriptor name="PubchemFingerprinter" value="false"/>',
    )
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(xml)


def run_padel_batch(smi_path, output_csv, desc_xml):
    padeldescriptor(
        mol_dir=smi_path, d_file=output_csv, descriptortypes=desc_xml,
        fingerprints=True, removesalt=True, retainorder=True, threads=-1,
    )


def main():
    print(f"Loading {INPUT_FILE}...")
    rows = load_csv(INPUT_FILE)
    print(f"  Loaded {len(rows)} rows")

    smiles_map = {}
    for r in rows:
        cid = r["molecule_chembl_id"]
        smi = r["canonical_smiles"]
        if cid not in smiles_map and smi:
            smiles_map[cid] = smi
    print(f"  Unique molecules: {len(smiles_map)}")

    tmpdir = tempfile.mkdtemp()
    smi_path = os.path.join(tmpdir, "molecules.smi")
    fp_csv_path = os.path.join(tmpdir, "fingerprints.csv")
    desc_xml = os.path.join(tmpdir, "descriptors.xml")

    print(f"\nWriting SMILES file...")
    written, skipped = write_smi(smiles_map, smi_path)
    print(f"  Written: {written}, Skipped: {len(skipped)}")

    make_descriptor_xml(desc_xml)

    print(f"\nRunning PaDEL batch AtomPairs2D Count calculation (SMILES input)...")
    start = time.time()
    run_padel_batch(smi_path, fp_csv_path, desc_xml)
    elapsed = time.time() - start
    print(f"  PaDEL completed in {elapsed:.1f}s ({elapsed / written * 1000:.0f}ms per molecule)")

    with open(fp_csv_path, newline="", encoding="utf-8") as f:
        fp_rows = list(csv.DictReader(f))
    print(f"  PaDEL output rows: {len(fp_rows)}")

    # Discover column names dynamically (descriptive names like APC2D1_C_C)
    fp_cols = [c for c in fp_rows[0].keys() if c != "Name"]
    print(f"  Discovered {len(fp_cols)} APC2D count columns")

    fp_map = {}
    for fr in fp_rows:
        name = fr.get("Name", "")
        fp_dict = {c: fr.get(c, "0") for c in fp_cols}
        fp_map[name] = fp_dict
    print(f"  Mapped {len(fp_map)} molecules to APC2D counts")

    out_fieldnames = ["molecule_chembl_id"] + fp_cols
    missing_fp = 0
    for row in rows:
        cid = row["molecule_chembl_id"]
        fp = fp_map.get(cid)
        if fp:
            row.update(fp)
        else:
            missing_fp += 1
            row.update({c: "" for c in fp_cols})

    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    os.remove(smi_path)
    os.remove(fp_csv_path)
    os.remove(desc_xml)
    os.rmdir(tmpdir)

    print(f"\n{'=' * 55}")
    print(f"Results saved to: {OUTPUT_FILE}")
    print(f"{'=' * 55}")
    print(f"  Rows: {len(rows)}")
    print(f"  Columns: {len(out_fieldnames)} (1 ID + {len(fp_cols)} APC2D counts)")
    print(f"  Rows with APC2D counts: {len(rows) - missing_fp}")
    print(f"  Rows without APC2D counts: {missing_fp}")


if __name__ == "__main__":
    main()
