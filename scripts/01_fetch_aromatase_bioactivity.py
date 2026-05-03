"""
Fetch bioactivity data (Ki, IC50, pIC50, pChEMBL) for Aromatase (CHEMBL1978)
from the ChEMBL REST API and save to CSV.
"""

import csv
import json
import os
import time
import urllib.request
import urllib.error

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TARGET_ID = "CHEMBL1978"  # Human Aromatase (CYP19A1)
STANDARD_TYPES = ["Ki", "IC50", "pIC50"]
PAGE_SIZE = 1000
OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "aromatase_bioactivity.csv")

COLUMNS = [
    "activity_id",
    "molecule_chembl_id",
    "molecule_pref_name",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "assay_chembl_id",
    "assay_description",
    "assay_type",
    "document_chembl_id",
    "document_journal",
    "document_year",
    "target_chembl_id",
    "target_pref_name",
]


def fetch_page(standard_type, offset=0):
    """Fetch one page of activity results from ChEMBL."""
    url = (
        f"{BASE_URL}?target_chembl_id={TARGET_ID}"
        f"&standard_type={standard_type}"
        f"&limit={PAGE_SIZE}&offset={offset}"
    )
    for attempt in range(3):
        try:
            req = urllib.request.Request(url)
            req.add_header("Accept", "application/json")
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt < 2:
                print(f"  Retry {attempt + 1} for {standard_type} offset={offset}: {e}")
                time.sleep(2 ** attempt)
            else:
                raise


def fetch_all(standard_type):
    """Fetch all pages for a given standard_type."""
    offset = 0
    all_activities = []
    while True:
        data = fetch_page(standard_type, offset)
        activities = data.get("activities", [])
        all_activities.extend(activities)
        total = data["page_meta"]["total_count"]
        print(f"  {standard_type}: fetched {len(all_activities)}/{total}")
        if data["page_meta"]["next"] is None:
            break
        offset += PAGE_SIZE
    return all_activities


def extract_row(activity):
    """Extract relevant fields from an activity record."""
    return {col: activity.get(col, "") for col in COLUMNS}


def main():
    all_rows = []
    type_counts = {}

    for stype in STANDARD_TYPES:
        print(f"Fetching {stype} data for Aromatase ({TARGET_ID})...")
        activities = fetch_all(stype)
        rows = [extract_row(a) for a in activities]
        all_rows.extend(rows)
        type_counts[stype] = len(rows)

    # Write CSV
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(all_rows)

    # Summary
    print(f"\n{'='*50}")
    print(f"Results saved to: {OUTPUT_FILE}")
    print(f"{'='*50}")
    print(f"{'Type':<12} {'Count':>8}")
    print(f"{'-'*12} {'-'*8}")
    for stype, count in type_counts.items():
        print(f"{stype:<12} {count:>8}")
    print(f"{'-'*12} {'-'*8}")
    print(f"{'TOTAL':<12} {sum(type_counts.values()):>8}")

    # pChEMBL coverage
    pchembl_count = sum(1 for r in all_rows if r.get("pchembl_value"))
    print(f"\nRecords with pChEMBL value: {pchembl_count}/{len(all_rows)}")


if __name__ == "__main__":
    main()
