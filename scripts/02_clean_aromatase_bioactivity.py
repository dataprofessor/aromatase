"""
Deduplicate aromatase bioactivity data with SD filtering and add InChIKeys.

For duplicate groups (same molecule + standard_type):
  - Compute SD of pchembl_value across measurements
  - Remove groups with SD > 3 (inconsistent data)
  - Collapse groups with SD <= 3 into a single row with mean values
  - Keep singleton groups as-is

Fetches InChIKeys from the ChEMBL molecule API.
"""

import csv
import json
import math
import os
import time
import urllib.request
import urllib.error
from collections import defaultdict

_DIR = os.path.dirname(__file__)
INPUT_FILE = os.path.join(_DIR, "..", "data", "raw", "aromatase_bioactivity.csv")
OUTPUT_FILE = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_clean.csv")
SD_THRESHOLD = 3.0
BATCH_SIZE = 50
MOL_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

OUTPUT_COLUMNS = [
    "molecule_chembl_id",
    "molecule_pref_name",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "sd_pchembl",
    "n_measurements",
    "assay_chembl_id",
    "assay_description",
    "assay_type",
    "document_chembl_id",
    "document_journal",
    "document_year",
    "target_chembl_id",
    "target_pref_name",
    "inchi_key",
]


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def calc_sd(vals):
    """Sample standard deviation."""
    n = len(vals)
    if n < 2:
        return 0.0
    mean = sum(vals) / n
    return math.sqrt(sum((x - mean) ** 2 for x in vals) / (n - 1))


def aggregate_groups(rows):
    """Group by (molecule_chembl_id, standard_type), filter by SD, and average."""
    groups = defaultdict(list)
    for r in rows:
        key = (r["molecule_chembl_id"], r["standard_type"])
        groups[key].append(r)

    result = []
    removed_groups = 0
    removed_rows = 0
    averaged_groups = 0

    for key, group in groups.items():
        # Collect pchembl_values for SD calculation
        pchembl_vals = []
        for r in group:
            pv = r.get("pchembl_value", "")
            if pv:
                try:
                    pchembl_vals.append(float(pv))
                except ValueError:
                    pass

        # If multiple pchembl measurements, check SD
        if len(pchembl_vals) > 1:
            sd = calc_sd(pchembl_vals)
            if sd > SD_THRESHOLD:
                removed_groups += 1
                removed_rows += len(group)
                continue

            # SD <= threshold: collapse to mean
            averaged_groups += 1
            mean_pchembl = sum(pchembl_vals) / len(pchembl_vals)

            # Average the standard_value too
            std_vals = []
            for r in group:
                sv = r.get("standard_value", "")
                if sv:
                    try:
                        std_vals.append(float(sv))
                    except ValueError:
                        pass
            mean_std_val = sum(std_vals) / len(std_vals) if std_vals else ""

            # Use first row as template for metadata
            out = dict(group[0])
            out["standard_value"] = f"{mean_std_val:.4f}" if mean_std_val != "" else ""
            out["pchembl_value"] = f"{mean_pchembl:.2f}"
            out["sd_pchembl"] = f"{sd:.4f}"
            out["n_measurements"] = str(len(group))
            result.append(out)
        else:
            # Singleton or single pchembl: keep as-is
            out = dict(group[0])
            sd = 0.0 if len(pchembl_vals) == 1 else ""
            out["sd_pchembl"] = f"{sd:.4f}" if isinstance(sd, float) else ""
            out["n_measurements"] = str(len(group))
            result.append(out)

    return result, removed_groups, removed_rows, averaged_groups


def fetch_inchikeys(chembl_ids):
    """Fetch InChIKeys for a list of molecule ChEMBL IDs in batches."""
    inchikey_map = {}
    batches = [chembl_ids[i:i + BATCH_SIZE] for i in range(0, len(chembl_ids), BATCH_SIZE)]

    for i, batch in enumerate(batches):
        ids_param = ",".join(batch)
        url = (
            f"{MOL_URL}?molecule_chembl_id__in={ids_param}"
            f"&only=molecule_chembl_id,molecule_structures&limit={BATCH_SIZE}"
        )

        for attempt in range(3):
            try:
                req = urllib.request.Request(url)
                req.add_header("Accept", "application/json")
                with urllib.request.urlopen(req, timeout=60) as resp:
                    data = json.loads(resp.read().decode("utf-8"))
                break
            except (urllib.error.URLError, TimeoutError) as e:
                if attempt < 2:
                    print(f"  Retry {attempt + 1} for batch {i + 1}: {e}")
                    time.sleep(2 ** attempt)
                else:
                    raise

        for mol in data.get("molecules", []):
            mid = mol["molecule_chembl_id"]
            structs = mol.get("molecule_structures") or {}
            inchikey_map[mid] = structs.get("standard_inchi_key", "")

        fetched = min((i + 1) * BATCH_SIZE, len(chembl_ids))
        print(f"  InChIKeys: fetched {fetched}/{len(chembl_ids)} molecules")

        if i < len(batches) - 1:
            time.sleep(0.3)

    return inchikey_map


def main():
    print(f"Loading {INPUT_FILE}...")
    rows = load_csv(INPUT_FILE)
    print(f"  Loaded {len(rows)} rows")

    # Group, filter by SD, and average
    print(f"\nGrouping by (molecule, standard_type) and filtering SD(pchembl) > {SD_THRESHOLD}...")
    clean_rows, removed_groups, removed_rows, averaged_groups = aggregate_groups(rows)
    print(f"  Groups removed (SD > {SD_THRESHOLD}): {removed_groups} ({removed_rows} rows)")
    print(f"  Groups averaged (SD <= {SD_THRESHOLD}): {averaged_groups}")
    print(f"  Remaining rows: {len(clean_rows)}")

    # Fetch InChIKeys
    mol_ids = sorted(set(r["molecule_chembl_id"] for r in clean_rows))
    print(f"\nFetching InChIKeys for {len(mol_ids)} unique molecules...")
    inchikey_map = fetch_inchikeys(mol_ids)

    # Add InChIKey to each row
    missing_inchikey = 0
    for row in clean_rows:
        row["inchi_key"] = inchikey_map.get(row["molecule_chembl_id"], "")
        if not row["inchi_key"]:
            missing_inchikey += 1

    # Write output
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(clean_rows)

    # Summary
    from collections import Counter
    type_counts = Counter(r["standard_type"] for r in clean_rows)

    print(f"\n{'='*55}")
    print(f"Results saved to: {OUTPUT_FILE}")
    print(f"{'='*55}")
    print(f"\n{'Metric':<40} {'Value':>8}")
    print(f"{'-'*40} {'-'*8}")
    print(f"{'Original rows':<40} {len(rows):>8}")
    print(f"{'Groups removed (SD > 3)':<40} {removed_groups:>8}")
    print(f"{'Rows removed':<40} {removed_rows:>8}")
    print(f"{'Groups averaged':<40} {averaged_groups:>8}")
    print(f"{'Final rows':<40} {len(clean_rows):>8}")
    print(f"{'Unique molecules':<40} {len(mol_ids):>8}")
    print(f"{'With InChIKey':<40} {len(clean_rows) - missing_inchikey:>8}")
    print(f"{'Missing InChIKey':<40} {missing_inchikey:>8}")
    print(f"\n{'Type':<12} {'Count':>8}")
    print(f"{'-'*12} {'-'*8}")
    for stype in ["Ki", "IC50", "pIC50"]:
        print(f"{stype:<12} {type_counts.get(stype, 0):>8}")
    print(f"{'-'*12} {'-'*8}")
    print(f"{'TOTAL':<12} {len(clean_rows):>8}")


if __name__ == "__main__":
    main()
