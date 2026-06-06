# Aromatase Bioactivity Dataset

## Overview

Bioactivity data for the human **Aromatase** (CYP19A1) protein target, retrieved from the ChEMBL database (target ID: **CHEMBL1978**). The dataset contains Ki, IC50, and pIC50 measurements with molecular fingerprints for QSAR/ML modeling.

## Project Structure

```
aromatase/
├── AGENTS.md
├── scripts/
│   ├── 01_fetch_aromatase_bioactivity.py      # Step 1: Fetch from ChEMBL API
│   ├── 02_clean_aromatase_bioactivity.py      # Step 2: Deduplicate + InChIKeys
│   ├── 03_eda_aromatase.py                    # Step 3: Exploratory data analysis
│   ├── 04_compute_maccs_fp.py                 # Step 4: MACCS keys (166 bits)
│   ├── 04_compute_pubchem_fp.py               # Step 4: PubChem FP (881 bits)
│   ├── 04_compute_substructure_fp.py          # Step 4: SubstructureFP (307 bits)
│   ├── 04_compute_klekota_roth_fp.py          # Step 4: Klekota-Roth FP (4860 bits)
│   ├── 04_compute_cdk_fp.py                   # Step 4: CDK Fingerprinter (1024 bits)
│   ├── 04_compute_cdk_ext_fp.py               # Step 4: CDK Extended FP (1024 bits)
│   ├── 04_compute_cdk_graphonly_fp.py          # Step 4: CDK GraphOnly FP (1024 bits)
│   ├── 04_compute_estate_fp.py                # Step 4: E-State FP (79 bits)
│   ├── 04_compute_atompairs2d_fp.py           # Step 4: AtomPairs2D FP (780 bits)
│   ├── 04_compute_substructure_count_fp.py    # Step 4: Substructure Count (307 counts)
│   ├── 04_compute_klekota_roth_count_fp.py    # Step 4: Klekota-Roth Count (4860 counts)
│   ├── 04_compute_atompairs2d_count_fp.py     # Step 4: AtomPairs2D Count (780 counts)
│   ├── 04_compute_ecfp4.py                    # Step 4: ECFP4 (1024 bits, optional)
│   ├── 05_remove_near_constant_fp.py          # Step 5: Remove near-constant features (SD < 0.1)
│   ├── 05b_remove_collinear_fp.py             # Step 5b: Remove collinear features (|r| >= 0.90)
│   ├── 06_split_data.py                       # Step 6: Train/test split (random + Kennard-Stone)
│   ├── 07_build_models.py                     # Step 7: Regression models (CLI)
│   ├── 07b_build_models.py                    # Step 7b: Regression models on decorrelated FPs (CLI)
│   ├── 08c_build_classification_models.py     # Step 8c: Classification models (CLI, 3-tier)
│   ├── 08d_build_classification_models.py     # Step 8d: Classification on decorrelated FPs (CLI, 3-tier)
│   └── data/
│       └── smarts_substructure.txt            # 307 SMARTS patterns (CDK/Inte:Ligand)
├── notebooks/
│   ├── 03_eda_aromatase.ipynb                 # EDA notebook
│   ├── 06_split_visualization.ipynb           # Train/test split visualization
│   ├── 07_build_models_colab.ipynb            # Regression (Colab GPU, original FPs)
│   ├── 07b_build_models_decorrelated_colab.ipynb  # Regression (Colab GPU, decorrelated FPs)
│   ├── 08_build_classification_models_colab.ipynb  # Classification (Colab GPU, original FPs)
│   ├── 08b_build_classification_models_decorrelated_colab.ipynb  # Classification (decorrelated)
│   ├── 08c_build_classification_models_colab.ipynb  # Classification 3-tier (original)
│   └── 08d_build_classification_models_decorrelated_colab.ipynb  # Classification 3-tier (decorrelated)
├── streamlit_app/
│   ├── streamlit_app.py                       # Main Streamlit app
│   ├── data_utils.py                          # Data loading utilities
│   └── app_pages/                             # Page modules (overview, chemical_space, etc.)
├── data/
│   ├── raw/
│   │   └── aromatase_bioactivity.csv                    # 5,097 rows, 17 cols
│   ├── processed/
│   │   ├── aromatase_bioactivity_clean.csv              # 3,774 rows, 19 cols
│   │   └── aromatase_bioactivity_curated.csv            # 3,774 rows, 19 cols
│   ├── fingerprints/                                    # Raw fingerprints (12 files)
│   ├── fingerprints_filtered/                           # After near-constant removal (12 files)
│   ├── fingerprints_decorrelated/                       # After collinearity removal (12 files)
│   ├── splits/                                          # Train/test splits (48 files)
│   ├── models/                                          # 07 regression results
│   ├── models_decorrelated/                             # 07b regression results
│   ├── models_classification/                           # 08c classification (default mode)
│   ├── models_classification_balanced/                  # 08c classification (balanced mode)
│   ├── models_classification_zips/                      # 08c per-model zip files
│   ├── models_classification_decorrelated/              # 08d classification (default mode)
│   ├── models_classification_balanced_decorrelated/     # 08d classification (balanced mode)
│   └── models_classification_decorrelated_zips/         # 08d per-model zip files
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

### 4. Fingerprint Computation (12 fingerprint types)

See individual script descriptions below for details on each fingerprint type.

### 5. Near-Constant Feature Removal (`scripts/05_remove_near_constant_fp.py`)

- Removes columns with SD < 0.1 from each fingerprint CSV
- Input: `data/fingerprints/` → Output: `data/fingerprints_filtered/`
- Reduces total features from ~15,000 to 4,897

### 5b. Collinear Feature Removal (`scripts/05b_remove_collinear_fp.py`)

- Computes Pearson intercorrelation matrix for each filtered fingerprint
- Greedy removal of features with |r| >= 0.90 (removes the feature involved in most correlated pairs first)
- Input: `data/fingerprints_filtered/` → Output: `data/fingerprints_decorrelated/`
- **Total reduction: 4,897 → 4,213 features (684 removed, 14%)**

| Fingerprint | Before | After | Removed | % |
|---|---:|---:|---:|---:|
| pubchem | 439 | 245 | 194 | 44% |
| klekota_roth_count | 515 | 376 | 139 | 27% |
| cdk_graphonly | 775 | 667 | 108 | 14% |
| klekota_roth | 439 | 339 | 100 | 23% |
| atompairs2d_count | 226 | 183 | 43 | 19% |
| atompairs2d | 183 | 142 | 41 | 22% |
| maccs | 133 | 111 | 22 | 17% |
| cdk_ext | 1002 | 992 | 10 | 1% |
| substructure_count | 70 | 61 | 9 | 13% |
| cdk | 1020 | 1012 | 8 | <1% |
| substructure | 65 | 57 | 8 | 12% |
| estate | 30 | 28 | 2 | 7% |

### 6. Train/Test Splitting (`scripts/06_split_data.py`)

- Reads from `data/fingerprints_decorrelated/`
- Two split strategies (both at molecule level to prevent data leakage):
  1. **Random split** — 80/20 with seed=42
  2. **Kennard-Stone split** — maximizes chemical diversity in training set
- Output: `data/splits/` (48 files: 12 FPs × 2 strategies × train/test)

### 7. Regression Model Building

**07 — Original fingerprints** (`scripts/07_build_models.py`)
- 16 regression models × 12 FPs × 2 splits = 384 fits
- Results: `data/models/` (24 test CSV files, 14-16 models each)

**07b — Decorrelated fingerprints** (`scripts/07b_build_models.py`)
- Same 16 models on decorrelated data
- Results: `data/models_decorrelated/` (24 test CSV files, 16 models each)

**Regression algorithms:**
1. Ridge Regression
2. Lasso Regression
3. ElasticNet
4. k-Nearest Neighbors
5. SVR (RBF)
6. SVR (Linear)
7. Decision Tree
8. Random Forest (500 trees)
9. Extra Trees (500 trees)
10. Gradient Boosting (100 trees)
11. XGBoost (500 trees)
12. AdaBoost (500 trees)
13. MLP Regressor (128-64)
14. Gaussian Process
15. Bayesian Ridge
16. PLS Regression

**CLI usage:**
```bash
python scripts/07_build_models.py --fp maccs --split random        # single FP
python scripts/07_build_models.py --fp cdk --split kennard_stone   # KS split
python scripts/07b_build_models.py --fp pubchem --split random     # decorrelated
```

### 8. Classification Model Building

**Task:** Predict 3-class bioactivity (active: pChEMBL > 7, intermediate: 6–7, inactive: < 6)

**08c — Original fingerprints** (`scripts/08c_build_classification_models.py`)
- 14 classification models × 12 FPs × default + balanced modes
- Results: `data/models_classification/` + `data/models_classification_balanced/`
- Per-model zips: `data/models_classification_zips/` (252 files)

**08d — Decorrelated fingerprints** (`scripts/08d_build_classification_models.py`)
- Same models on decorrelated data
- Results: `data/models_classification_decorrelated/` + `data/models_classification_balanced_decorrelated/`
- Per-model zips: `data/models_classification_decorrelated_zips/` (252 files)

**Classification algorithms (14 models):**

*Tier 1 — Fast (support class_weight):*
1. Ridge Classifier
2. Logistic Regression (L1)
3. Logistic Regression (ElasticNet)
4. Decision Tree

*Tier 1 — Fast (no class_weight):*
5. k-Nearest Neighbors
6. Gaussian Naive Bayes
7. Linear Discriminant Analysis

*Tier 2 — Medium:*
8. Gradient Boosting (100 trees)
9. XGBoost (500 trees)
10. AdaBoost (500 trees)
11. MLP Classifier (128-64)

*Tier 3 — Slow (support class_weight):*
12. SVC (RBF)
13. SVC (Linear)
14. Random Forest (500 trees)

**Two training modes:**
- **Default**: All 14 models without class weighting
- **Balanced**: 7 models with `class_weight='balanced'` (Ridge, LogReg ×2, SVC ×2, DecTree, RF)

**CLI usage:**
```bash
python scripts/08c_build_classification_models.py                        # all
python scripts/08c_build_classification_models.py --fp maccs --tier 1    # fast only
python scripts/08c_build_classification_models.py --tier 2 --mode balanced
python scripts/08d_build_classification_models.py --fp pubchem --tier 1-2
```

## Colab Notebooks

All notebooks clone from `https://github.com/dataprofessor/aromatase.git` and use GPU acceleration via cuML + XGBoost CUDA.

| Notebook | Task | Data | Notes |
|---|---|---|---|
| `07_build_models_colab.ipynb` | Regression | Original FPs | 14 models, GPU |
| `07b_build_models_decorrelated_colab.ipynb` | Regression | Decorrelated FPs | 14 models, GPU |
| `08_build_classification_models_colab.ipynb` | Classification | Original FPs | 14 models, default+balanced |
| `08b_build_classification_models_decorrelated_colab.ipynb` | Classification | Decorrelated FPs | 14 models, default+balanced |
| `08c_build_classification_models_colab.ipynb` | Classification | Original FPs | 3-tier with incremental zip downloads |
| `08d_build_classification_models_decorrelated_colab.ipynb` | Classification | Decorrelated FPs | 3-tier with incremental zip downloads |

All notebooks include a final cell that zips results and triggers `google.colab.files.download()`.

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

### Classification Class Definitions

| Class | pChEMBL Range | IC50 Range | Count (train) |
|---|---|---|---:|
| Active | > 7 | < 100 nM | 893 (33.2%) |
| Intermediate | 6–7 | 100 nM – 1 μM | 716 (26.6%) |
| Inactive | < 6 | > 1 μM | 1,082 (40.2%) |

## Dependencies

- **Python 3.9+**
- **RDKit** (`pip install rdkit-pypi`) — SMILES parsing, MACCS keys, ECFP4
- **PaDEL-Descriptor** (`pip install padelpy`) — PubChem FP, SubstructureFP, Klekota-Roth, CDK FP, CDK Extended, CDK GraphOnly, E-State, AtomPairs2D (requires Java runtime)
- **scikit-learn** — All ML models (Ridge, Lasso, SVM, RF, etc.)
- **XGBoost** (`pip install xgboost`) — Gradient boosted trees
- **NumPy** — Array operations, correlation matrices
- Standard library: `csv`, `json`, `urllib`, `math`, `time`, `collections`, `zipfile`

### Colab-specific dependencies (installed at runtime)
- **cuML** (`cuml-cu12`) — GPU-accelerated sklearn-compatible models
- **cuDF** (`cudf-cu12`) — GPU DataFrame operations
- **tqdm** — Progress bars

## Reproducibility

Run scripts from anywhere (paths are relative to script location):

```bash
# Data pipeline
python scripts/01_fetch_aromatase_bioactivity.py         # ~2 min (API calls)
python scripts/02_clean_aromatase_bioactivity.py         # ~4 min (API calls for InChIKeys)
python scripts/04_compute_maccs_fp.py                    # ~5 sec
python scripts/04_compute_pubchem_fp.py                  # ~5 min (PaDEL batch)
python scripts/04_compute_substructure_fp.py             # ~76 sec (PaDEL batch)
python scripts/04_compute_klekota_roth_fp.py             # ~11 min (PaDEL batch, multithreaded)
python scripts/04_compute_cdk_fp.py                      # ~17 sec
python scripts/04_compute_cdk_ext_fp.py                  # ~28 sec
python scripts/04_compute_cdk_graphonly_fp.py            # ~14 sec
python scripts/04_compute_estate_fp.py                   # ~22 sec
python scripts/04_compute_atompairs2d_fp.py              # ~21 sec
python scripts/04_compute_substructure_count_fp.py       # ~46 sec
python scripts/04_compute_klekota_roth_count_fp.py       # ~12 min
python scripts/04_compute_atompairs2d_count_fp.py        # ~16 sec

# Feature selection
python scripts/05_remove_near_constant_fp.py             # ~10 sec
python scripts/05b_remove_collinear_fp.py                # ~30 sec

# Splitting
python scripts/06_split_data.py                          # ~5 sec

# Model building (local CLI)
python scripts/07_build_models.py --fp maccs --split random
python scripts/07b_build_models.py --fp maccs --split random
python scripts/08c_build_classification_models.py --tier 1
python scripts/08d_build_classification_models.py --tier 1
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

**Workaround**: The notebook's error handling (`try/except`) catches and skips this model/fingerprint combination. The CLI scripts (`07_build_models.py`) use sklearn's SVR which works without issue. Results for other fingerprints and SVR (RBF) are unaffected.

**Resolution**: The missing SVR (Linear) results for these 4 combinations were filled using sklearn locally via the CLI scripts.

### XGBoost + sklearn 1.6 tag incompatibility

**Error**: `'super' object has no attribute '__sklearn_tags__'` when using `cross_val_predict()` or `clone()` with XGBoost 1.6.x and sklearn >= 1.6.

**Context**: sklearn 1.6 introduced a new `__sklearn_tags__` protocol that XGBoost 1.x doesn't implement.

**Workaround**: The CLI scripts use manual CV loops for XGBoost (bypass `cross_val_predict`/`clone`):
```python
if name == "XGBoost":
    y_pred_cv = np.zeros_like(y_train)
    for tr_idx, val_idx in cv.split(X_train, y_train):
        m_cv = type(model)(**model.get_params())
        m_cv.fit(X_train[tr_idx], y_train[tr_idx])
        y_pred_cv[val_idx] = m_cv.predict(X_train[val_idx])
```

This issue does not affect Colab notebooks (which use cuML + compatible XGBoost versions).

### XGBoost libomp architecture mismatch (macOS)

**Error**: `XGBoost Library (libxgboost.dylib) could not be loaded` — needs x86_64 `libomp.dylib` but only arm64 available via Homebrew.

**Context**: When Python is x86_64 (miniconda under Rosetta) but brew-installed libomp is arm64.

**Workaround**: Use XGBoost 1.6.1 (`pip install xgboost==1.6.1`) which bundles its own compatible libomp. XGBoost 2.x requires a matching-architecture libomp. This is a local dev environment issue only — does not affect Colab.

## Model Results Summary

### Regression (07/07b) — Test R² (top performers)

Models trained on 3,774 molecules, 80/20 split, 10-fold CV.

Best models consistently across fingerprints:
- **Random Forest** and **XGBoost**: R² ≈ 0.45–0.55
- **Gradient Boosting**: R² ≈ 0.40–0.50
- **k-Nearest Neighbors**: R² ≈ 0.35–0.45

### Classification (08c/08d) — Test Balanced Accuracy (top performers)

3-class prediction (active/intermediate/inactive):
- **Random Forest** (balanced): BalAcc ≈ 0.63–0.67
- **XGBoost**: BalAcc ≈ 0.58–0.65
- **MLP Classifier**: BalAcc ≈ 0.63–0.66
- **SVC (RBF)**: BalAcc ≈ 0.59–0.65

Best fingerprints: MACCS, PubChem, CDK (hashed path-based fingerprints capture relevant structural features for aromatase inhibition).

## To-Do

1. **Feature importance analysis** — Identify the most predictive fingerprint bits/features for aromatase inhibition. Use permutation importance, SHAP values, or model-native feature importances (Random Forest, XGBoost) to rank descriptors across fingerprint types.

2. **Statistical comparison of bioactivity classes** — Perform statistical tests (e.g., Kruskal-Wallis, Mann-Whitney U with Bonferroni correction) to compare molecular descriptor distributions across the 3 bioactivity classes (active, intermediate, inactive). Identify which descriptors significantly differentiate the classes.

3. **Applicability domain analysis** — Define the chemical space coverage of the training set and flag test predictions that fall outside this domain.

4. **Consensus/ensemble models** — Combine predictions from the top 3-5 models per fingerprint type to potentially improve robustness.
