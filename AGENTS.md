# Aromatase Bioactivity Dataset

## Overview

Bioactivity data for the human **Aromatase** (CYP19A1) protein target, retrieved from the ChEMBL database (target ID: **CHEMBL1978**). The dataset contains Ki, IC50, and pIC50 measurements with molecular fingerprints for QSAR/ML modeling.

## Project Structure

```
aromatase/
├── AGENTS.md
├── scripts/
│   ├── 01_fetch_aromatase_bioactivity.py  # Step 1: Fetch from ChEMBL API
│   ├── 02_clean_aromatase_bioactivity.py  # Step 2: Deduplicate + InChIKeys
│   ├── 03_eda_aromatase.py                # Step 3: Exploratory data analysis
│   ├── 04_compute_maccs_fp.py             # Step 4: MACCS keys (166 bits)
│   ├── 04_compute_pubchem_fp.py           # Step 4: PubChem FP (881 bits)
│   ├── 04_compute_substructure_fp.py      # Step 4: SubstructureFP (307 bits)
│   ├── 04_compute_klekota_roth_fp.py      # Step 4: Klekota-Roth FP (4860 bits)
│   ├── 04_compute_cdk_fp.py              # Step 4: CDK Fingerprinter (1024 bits)
│   ├── 04_compute_cdk_ext_fp.py          # Step 4: CDK Extended FP (1024 bits)
│   ├── 04_compute_cdk_graphonly_fp.py     # Step 4: CDK GraphOnly FP (1024 bits)
│   ├── 04_compute_estate_fp.py           # Step 4: E-State FP (79 bits)
│   ├── 04_compute_atompairs2d_fp.py      # Step 4: AtomPairs2D FP (780 bits)
│   ├── 04_compute_substructure_count_fp.py  # Step 4: Substructure Count (307 counts)
│   ├── 04_compute_klekota_roth_count_fp.py  # Step 4: Klekota-Roth Count (4860 counts)
│   ├── 04_compute_atompairs2d_count_fp.py   # Step 4: AtomPairs2D Count (780 counts)
│   ├── 04_compute_ecfp4.py               # Step 4: ECFP4 (1024 bits, optional)
│   └── data/
│       └── smarts_substructure.txt      # 307 SMARTS patterns (CDK/Inte:Ligand)
├── data/
│   ├── raw/
│   │   └── aromatase_bioactivity.csv            # 5,097 rows, 17 cols
│   ├── processed/
│   │   ├── aromatase_bioactivity_clean.csv      # 3,774 rows, 19 cols
│   │   └── aromatase_bioactivity_curated.csv    # 3,774 rows, 19 cols
│   └── fingerprints/
│       ├── aromatase_maccs_fp.csv               # 3,774 rows, 167 cols
│       ├── aromatase_pubchem_fp.csv             # 3,774 rows, 882 cols
│       ├── aromatase_substructure_fp.csv        # 3,774 rows, 308 cols
│       ├── aromatase_klekota_roth_fp.csv       # 3,774 rows, 4861 cols
│       ├── aromatase_cdk_fp.csv                # 3,774 rows, 1025 cols
│       ├── aromatase_cdk_ext_fp.csv            # 3,774 rows, 1025 cols
│       ├── aromatase_cdk_graphonly_fp.csv       # 3,774 rows, 1025 cols
│       ├── aromatase_estate_fp.csv             # 3,774 rows, 80 cols
│       ├── aromatase_atompairs2d_fp.csv        # 3,774 rows, 781 cols
│       ├── aromatase_substructure_count_fp.csv # 3,774 rows, 308 cols
│       ├── aromatase_klekota_roth_count_fp.csv # 3,774 rows, 4861 cols
│       └── aromatase_atompairs2d_count_fp.csv  # 3,774 rows, 781 cols
```

## Data Pipeline

### 1. Data Retrieval (`scripts/01_fetch_aromatase_bioactivity.py`)

- Queries the ChEMBL REST API (`/activity.json`) for aromatase (CHEMBL1978)
- Retrieves all records for `standard_type` in {Ki, IC50, pIC50} with pagination (1000/page)
- Includes `pchembl_value` (-log10 molar potency) from each record
- Output: `data/raw/aromatase_bioactivity.csv` (5,097 rows, 17 columns)

### 2. Deduplication and Cleaning (`scripts/02_clean_aromatase_bioactivity.py`)

- Groups records by `(molecule_chembl_id, standard_type)`
- For groups with multiple pchembl_value measurements:
  - Computes SD of pchembl_value
  - **Removes** groups with SD > 3 (inconsistent data) — 2 groups / 4 rows removed
  - **Averages** groups with SD <= 3 into single rows with mean values — 668 groups collapsed
- Singleton groups kept as-is
- Fetches **InChIKeys** from ChEMBL molecule API in batches of 50
- Adds columns: `sd_pchembl`, `n_measurements`, `inchi_key`
- Output: `data/processed/aromatase_bioactivity_clean.csv` (3,774 rows, 19 columns)

### 3. Exploratory Data Analysis (`scripts/03_eda_aromatase.py`)

- Exploratory data analysis of the cleaned bioactivity dataset

### 4a. MACCS Fingerprints (`scripts/04_compute_maccs_fp.py`)

- Computes 166-bit MACCS structural keys using RDKit
- Each bit represents a specific structural feature (e.g., aromatic ring, nitrogen in ring, ester group)
- Outputs:
  - `data/processed/aromatase_bioactivity_curated.csv` — base bioactivity data (19 columns)
  - `data/fingerprints/aromatase_maccs_fp.csv` — `molecule_chembl_id` + 166 MACCS bits

### 4b. PubChem Fingerprints (`scripts/04_compute_pubchem_fp.py`)

- Computes 881-bit PubChem substructure fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Writes a `.smi` file (SMILES + molecule ID) and passes it to PaDEL (~91ms/molecule)
- 881 bits across 7 sections: element counts, ring counts, atom pairs, atom environments, SMARTS substructures
- Cross-validated against per-molecule PaDEL: **exact match** on all 881 bits (5/5 test molecules)
- Output: `data/fingerprints/aromatase_pubchem_fp.csv` — `molecule_chembl_id` + 881 PubChem bits

> **Note**: An earlier version used SDF input, which caused 3-6 bit differences in Section 2 (ring count bits 115-262) due to aromaticity perception differences in RDKit's SDF writer vs CDK's SMILES parser. Switching to `.smi` input eliminated these discrepancies entirely.

### 4c. SubstructureFP (`scripts/04_compute_substructure_fp.py`)

- Computes 307-bit functional group fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Uses CDK's Inte:Ligand classification (307 SMARTS patterns covering amines, ketones, lactams, heterocycles, etc.)
- Generates a custom `descriptors.xml` at runtime to enable only SubstructureFingerprinter
- ~22 ms/mol (~76s for full dataset)
- Cross-validated against per-molecule PaDEL: **exact match** on all 307 bits (5/5 test molecules)
- Output: `data/fingerprints/aromatase_substructure_fp.csv` — `molecule_chembl_id` + 307 SubFP bits

### 4d. Klekota-Roth Fingerprints (`scripts/04_compute_klekota_roth_fp.py`)

- Computes 4860-bit Klekota-Roth substructure fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Uses CDK's KlekotaRothFingerprinter (4860 SMARTS patterns defined by Klekota & Roth, 2008)
- Generates a custom `descriptors.xml` at runtime to enable only KlekotaRothFingerprinter
- ~187 ms/mol (~637s for full dataset with multithreading, `threads=-1`)
- Cross-validated against per-molecule PaDEL: **exact match** on all 4860 bits (5/5 test molecules)
- Output: `data/fingerprints/aromatase_klekota_roth_fp.csv` — `molecule_chembl_id` + 4860 KRFP bits
- Bit density: avg 45.8/4860 bits ON (0.9%), range 1-198

### 4e. CDK Fingerprinter (`scripts/04_compute_cdk_fp.py`)

- Computes 1024-bit CDK hashed fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Uses CDK's default Fingerprinter (path-based hashed fingerprint, Daylight-type)
- ~5 ms/mol (~17s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 1024 bits (3/3 test molecules)
- Output: `data/fingerprints/aromatase_cdk_fp.csv` — `molecule_chembl_id` + 1024 FP bits
- Bit density: avg 251.0/1024 bits ON (24.5%)

### 4f. CDK Extended Fingerprinter (`scripts/04_compute_cdk_ext_fp.py`)

- Computes 1024-bit CDK Extended fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Extends the default CDK Fingerprinter with additional ring and atom-type features
- ~8 ms/mol (~28s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 1024 bits (3/3 test molecules)
- Output: `data/fingerprints/aromatase_cdk_ext_fp.csv` — `molecule_chembl_id` + 1024 ExtFP bits
- Bit density: avg 265.0/1024 bits ON (25.9%)

### 4g. CDK GraphOnly Fingerprinter (`scripts/04_compute_cdk_graphonly_fp.py`)

- Computes 1024-bit CDK GraphOnly fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Uses only graph topology (ignores element types), capturing pure connectivity patterns
- ~4 ms/mol (~14s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 1024 bits (3/3 test molecules)
- Output: `data/fingerprints/aromatase_cdk_graphonly_fp.csv` — `molecule_chembl_id` + 1024 GraphFP bits
- Bit density: avg 126.3/1024 bits ON (12.3%)

### 4h. E-State Fingerprinter (`scripts/04_compute_estate_fp.py`)

- Computes 79-bit E-State fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Based on electrotopological state atom types (Hall & Kier)
- ~6 ms/mol (~22s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 79 bits (3/3 test molecules)
- Output: `data/fingerprints/aromatase_estate_fp.csv` — `molecule_chembl_id` + 79 EStateFP bits
- Bit density: avg 8.6/79 bits ON (10.9%)

### 4i. AtomPairs2D Fingerprinter (`scripts/04_compute_atompairs2d_fp.py`)

- Computes 780-bit AtomPairs2D fingerprints using PaDEL-Descriptor in **batch SMILES mode**
- Encodes presence of atom-type pairs at various topological distances
- ~6 ms/mol (~21s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 780 bits (3/3 test molecules)
- Output: `data/fingerprints/aromatase_atompairs2d_fp.csv` — `molecule_chembl_id` + 780 AD2D bits
- Bit density: avg 36.5/780 bits ON (4.7%)

### 4j. Substructure Count (`scripts/04_compute_substructure_count_fp.py`)

- Computes 307 frequency counts using PaDEL-Descriptor in **batch SMILES mode**
- Same 307 SMARTS patterns as SubstructureFP but records occurrence count instead of binary presence
- ~14 ms/mol (~46s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 307 counts (3/3 test molecules)
- Output: `data/fingerprints/aromatase_substructure_count_fp.csv` — `molecule_chembl_id` + 307 SubFPC counts

### 4k. Klekota-Roth Count (`scripts/04_compute_klekota_roth_count_fp.py`)

- Computes 4860 frequency counts using PaDEL-Descriptor in **batch SMILES mode**
- Same 4860 SMARTS patterns as KRFP but records occurrence count instead of binary presence
- ~206 ms/mol (~701s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 4860 counts (3/3 test molecules)
- Output: `data/fingerprints/aromatase_klekota_roth_count_fp.csv` — `molecule_chembl_id` + 4860 KRFPC counts

### 4l. AtomPairs2D Count (`scripts/04_compute_atompairs2d_count_fp.py`)

- Computes 780 frequency counts using PaDEL-Descriptor in **batch SMILES mode**
- Same atom-pair types as AtomPairs2D but records frequency count instead of binary presence
- Column names are descriptive (e.g., `APC2D1_C_C`, `APC2D1_C_N`) rather than numeric indices
- ~5 ms/mol (~16s for full dataset with multithreading)
- Cross-validated against per-molecule PaDEL: **exact match** on all 780 counts (3/3 test molecules)
- Output: `data/fingerprints/aromatase_atompairs2d_count_fp.csv` — `molecule_chembl_id` + 780 APC2D counts

### 4m. ECFP4 Fingerprints (`scripts/04_compute_ecfp4.py`)

- Computes 1024-bit Morgan/ECFP4 circular fingerprints (radius=2) using RDKit
- Not included in final output (replaced by PubChem FP) but script is available
- Output: `data/fingerprints/aromatase_ecfp4_fp.csv` (if run)

## Output Files

All fingerprint CSVs use `molecule_chembl_id` as the join key to link back to bioactivity data.

| File | Rows | Columns | Description |
|---|---:|---:|---|
| `data/processed/aromatase_bioactivity_curated.csv` | 3,774 | 19 | Bioactivity measurements, molecule info, assay metadata |
| `data/fingerprints/aromatase_maccs_fp.csv` | 3,774 | 167 | molecule_chembl_id + 166 MACCS key bits |
| `data/fingerprints/aromatase_pubchem_fp.csv` | 3,774 | 882 | molecule_chembl_id + 881 PubChem FP bits |
| `data/fingerprints/aromatase_substructure_fp.csv` | 3,774 | 308 | molecule_chembl_id + 307 SubFP bits |
| `data/fingerprints/aromatase_klekota_roth_fp.csv` | 3,774 | 4,861 | molecule_chembl_id + 4860 KRFP bits |
| `data/fingerprints/aromatase_cdk_fp.csv` | 3,774 | 1,025 | molecule_chembl_id + 1024 CDK FP bits |
| `data/fingerprints/aromatase_cdk_ext_fp.csv` | 3,774 | 1,025 | molecule_chembl_id + 1024 CDK Extended bits |
| `data/fingerprints/aromatase_cdk_graphonly_fp.csv` | 3,774 | 1,025 | molecule_chembl_id + 1024 CDK GraphOnly bits |
| `data/fingerprints/aromatase_estate_fp.csv` | 3,774 | 80 | molecule_chembl_id + 79 E-State bits |
| `data/fingerprints/aromatase_atompairs2d_fp.csv` | 3,774 | 781 | molecule_chembl_id + 780 AtomPairs2D bits |
| `data/fingerprints/aromatase_substructure_count_fp.csv` | 3,774 | 308 | molecule_chembl_id + 307 SubFPC counts |
| `data/fingerprints/aromatase_klekota_roth_count_fp.csv` | 3,774 | 4,861 | molecule_chembl_id + 4860 KRFPC counts |
| `data/fingerprints/aromatase_atompairs2d_count_fp.csv` | 3,774 | 781 | molecule_chembl_id + 780 APC2D counts |

### Base Columns

| Column | Description |
|---|---|
| `molecule_chembl_id` | ChEMBL molecule identifier |
| `molecule_pref_name` | Preferred molecule name |
| `canonical_smiles` | Canonical SMILES string |
| `standard_type` | Measurement type: Ki, IC50, or pIC50 |
| `standard_relation` | Relation operator (=, >, <) |
| `standard_value` | Measured value (mean if averaged) |
| `standard_units` | Units (nM for Ki/IC50, uM for pIC50) |
| `pchembl_value` | -log10 molar potency (mean if averaged) |
| `sd_pchembl` | SD of pchembl_value across duplicate measurements |
| `n_measurements` | Number of original measurements averaged |
| `assay_chembl_id` | ChEMBL assay identifier |
| `assay_description` | Assay description text |
| `assay_type` | Assay type (B=binding, F=functional) |
| `document_chembl_id` | Source publication ChEMBL ID |
| `document_journal` | Journal name |
| `document_year` | Publication year |
| `target_chembl_id` | Target ID (CHEMBL1978) |
| `target_pref_name` | Target name (Aromatase) |
| `inchi_key` | InChIKey molecular identifier |

### Record Counts by Type

| Type | Count |
|---|---:|
| Ki | 516 |
| IC50 | 3,211 |
| pIC50 | 47 |
| **Total** | **3,774** |

## Dependencies

- **Python 3.9+**
- **RDKit** (`pip install rdkit-pypi`) — SMILES parsing, MACCS keys, ECFP4
- **PaDEL-Descriptor** (`pip install padelpy`) — PubChem FP, SubstructureFP, Klekota-Roth, CDK FP, CDK Extended, CDK GraphOnly, E-State, AtomPairs2D (requires Java runtime)
- Standard library: `csv`, `json`, `urllib`, `math`, `time`, `collections`

## Reproducibility

Run scripts from anywhere (paths are relative to script location):

```bash
python scripts/01_fetch_aromatase_bioactivity.py   # ~2 min (API calls)
python scripts/02_clean_aromatase_bioactivity.py   # ~4 min (API calls for InChIKeys)
python scripts/04_compute_maccs_fp.py              # ~5 sec
python scripts/04_compute_pubchem_fp.py            # ~5 min (PaDEL batch)
python scripts/04_compute_substructure_fp.py       # ~76 sec (PaDEL batch)
python scripts/04_compute_klekota_roth_fp.py       # ~11 min (PaDEL batch, multithreaded)
python scripts/04_compute_cdk_fp.py                # ~17 sec (PaDEL batch, multithreaded)
python scripts/04_compute_cdk_ext_fp.py            # ~28 sec (PaDEL batch, multithreaded)
python scripts/04_compute_cdk_graphonly_fp.py       # ~14 sec (PaDEL batch, multithreaded)
python scripts/04_compute_estate_fp.py             # ~22 sec (PaDEL batch, multithreaded)
python scripts/04_compute_atompairs2d_fp.py        # ~21 sec (PaDEL batch, multithreaded)
python scripts/04_compute_substructure_count_fp.py  # ~46 sec (PaDEL batch, multithreaded)
python scripts/04_compute_klekota_roth_count_fp.py  # ~12 min (PaDEL batch, multithreaded)
python scripts/04_compute_atompairs2d_count_fp.py   # ~16 sec (PaDEL batch, multithreaded)
```

## Notes for Future Native Python Implementation

Both PubChem FP and SubstructureFP currently use PaDEL (Java/CDK) for computation. A native Python/RDKit implementation was attempted for both and revealed a recurring issue: **RDKit and CDK perceive aromaticity differently**, which causes bit-level mismatches in SMARTS-based fingerprints.

Key findings for anyone porting these fingerprints to pure Python:

- **Always use SMILES input** (not SDF) when interfacing with PaDEL/CDK. SDF files generated by RDKit embed RDKit's aromaticity model, which CDK reinterprets differently. SMILES lets CDK apply its own perception from scratch.
- **PubChem FP**: Sections 1 (element counts) and 3-7 (SMARTS substructures) are straightforward to port. Section 2 (ESSSR ring counts, bits 115-262) requires CDK's envelope ring algorithm, which differs from RDKit's ring perception.
- **SubstructureFP**: 306/307 SMARTS patterns parse in RDKit (bit 298 "Salt" needs a custom formal-charge check). However, ~10 bits consistently disagree due to aromaticity: tautomerizable groups (bits 300-301), annelated/bridged rings (279-280), oxoarene (177), CH-acidic (304-305), and cis/trans double bonds (289-290).
- **SMARTS patterns source**: CDK's `SMARTS_InteLigand.txt` (307 patterns, LGPL, by Christian Laggner/Inte:Ligand). A copy is stored in `scripts/data/smarts_substructure.txt`.
- **Performance target**: Native RDKit SMARTS matching runs at ~5 ms/mol (vs ~22 ms/mol PaDEL), so a ~4x speedup is achievable if the aromaticity differences are acceptable or can be resolved.

## Known Issues

### cuML SVR (Linear) crash on high-dimensional data

**Error**: `klekota_roth_count/random/SVR (Linear)` — cuML's SVR with `kernel='linear'` crashes with:
```
exception occurred! file=/__w/cuml/cuml/cpp/src/svm/kernelcache.cuh line=487:
Working set has already been initialized!
```

**Context**: Occurs on multiple fingerprints with `cuSVR(kernel='linear')` during 10-fold cross-validation:
- `klekota_roth_count` (515 features, 3022 samples) — random split
- `klekota_roth_count` (515 features, 3022 samples) — kennard_stone split
- `atompairs2d_count` (226 features, 3022 samples) — random split
- `atompairs2d_count` (226 features, 3022 samples) — kennard_stone split

The error originates in cuML's C++ SMO solver (`SmoSolver::Solve` → `KernelCache::InitWorkingSet`). Likely a cuML bug with the working set initialization being called twice on linear SVR. Affects count-based fingerprints more frequently than binary fingerprints.

**Workaround**: The notebook's error handling (`try/except`) catches and skips this model/fingerprint combination. Results for other fingerprints and SVR (RBF) are unaffected.

**Future fix options**:
1. Fall back to sklearn's `SVR(kernel='linear')` for fingerprints with >400 features
2. Replace `cuSVR(kernel='linear')` with `cuml.linear_model.LinearSVR` (different solver, avoids the SMO path)
3. Wait for cuML upstream fix (file issue at github.com/rapidsai/cuml)

## To-Do

1. **Run regression/classification models with collinear-removed data** — Re-run both the regression (`07_build_models_colab.ipynb`) and classification (`08_build_classification_models_colab.ipynb`) notebooks using the fingerprint data after removing highly collinear features. Compare performance to the full-feature models.

2. **Feature importance analysis** — Identify the most predictive fingerprint bits/features for aromatase inhibition. Use permutation importance, SHAP values, or model-native feature importances (Random Forest, XGBoost) to rank descriptors across fingerprint types.

3. **Statistical comparison of bioactivity classes** — Perform statistical tests (e.g., Kruskal-Wallis, Mann-Whitney U with Bonferroni correction) to compare molecular descriptor distributions across the 3 bioactivity classes (active, intermediate, inactive). Identify which descriptors significantly differentiate the classes.
