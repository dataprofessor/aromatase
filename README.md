# Aromatase (CYP19A1) QSAR Modeling

Bioactivity prediction for human **Aromatase** (CYP19A1) using molecular fingerprints and machine learning. Data sourced from ChEMBL (target: CHEMBL1978).

## Quick Start

```bash
git clone https://github.com/dataprofessor/aromatase.git
cd aromatase
pip install numpy pandas scipy scikit-learn xgboost rdkit-pypi padelpy matplotlib seaborn tqdm
```

**Java runtime** required for PaDEL-Descriptor (fingerprint computation only).

All data, fingerprints, splits, and model results are included in the repo — no need to re-run the pipeline from scratch.

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_fetch_aromatase_bioactivity.py` | Fetch bioactivity data from ChEMBL API |
| 2 | `02_clean_aromatase_bioactivity.py` | Deduplicate, average replicates, add InChIKeys |
| 3 | `03_eda_aromatase.py` | Exploratory data analysis |
| 4 | `04_compute_*.py` | Compute 12 molecular fingerprint types |
| 5 | `05_remove_near_constant_fp.py` | Remove low-variance features (SD < 0.1) |
| 5b | `05b_remove_collinear_fp.py` | Remove collinear features (\|r\| >= 0.90) |
| 6 | `06_split_data.py` | Train/test split (random + Kennard-Stone) |
| 7 | `07_build_models.py` | Regression models (16 algorithms, CLI) |
| 7b | `07b_build_models.py` | Regression on decorrelated fingerprints |
| 8c | `08c_build_classification_models.py` | Classification (14 algorithms, 3-tier, CLI) |
| 8d | `08d_build_classification_models.py` | Classification on decorrelated fingerprints |

## Fingerprints (12 types)

MACCS (166-bit), PubChem (881-bit), Substructure (307-bit), Klekota-Roth (4860-bit), CDK (1024-bit), CDK Extended (1024-bit), CDK GraphOnly (1024-bit), E-State (79-bit), AtomPairs2D (780-bit), plus count versions of Substructure, Klekota-Roth, and AtomPairs2D.

**Feature selection**: Near-constant removal (SD < 0.1) followed by collinearity removal (|r| >= 0.90) reduces features from ~15,000 → 4,897 → 4,213.

## Colab Notebooks

| Notebook | Task | Colab |
|---|---|---|
| `07_build_models_colab.ipynb` | Regression (original FPs) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dataprofessor/aromatase/blob/main/notebooks/07_build_models_colab.ipynb) |
| `07b_build_models_decorrelated_colab.ipynb` | Regression (decorrelated) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dataprofessor/aromatase/blob/main/notebooks/07b_build_models_decorrelated_colab.ipynb) |
| `08c_build_classification_models_colab.ipynb` | Classification 3-tier (original) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dataprofessor/aromatase/blob/main/notebooks/08c_build_classification_models_colab.ipynb) |
| `08d_build_classification_models_decorrelated_colab.ipynb` | Classification 3-tier (decorrelated) | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dataprofessor/aromatase/blob/main/notebooks/08d_build_classification_models_decorrelated_colab.ipynb) |

Uses cuML (GPU) for Ridge, Lasso, ElasticNet, KNN, SVR, Random Forest and XGBoost with CUDA. Each notebook auto-zips and downloads results upon completion.

## Dataset

- **3,774 records** (516 Ki, 3,211 IC50, 47 pIC50)
- **Target variable**: pChEMBL value (-log10 molar potency)
- **Classification**: Active (pChEMBL > 7), Intermediate (6–7), Inactive (< 6)
- **Train/test**: 80/20 split at molecule level (seed=42)

## CLI Usage

```bash
# Regression (16 models per fingerprint)
python scripts/07_build_models.py --fp maccs --split random
python scripts/07b_build_models.py --fp pubchem --split kennard_stone

# Classification (14 models, 3-tier training)
python scripts/08c_build_classification_models.py --fp maccs --tier 1        # fast models only
python scripts/08c_build_classification_models.py --tier 1-2 --mode balanced # fast+medium, balanced
python scripts/08d_build_classification_models.py                            # all (decorrelated)
```

## Requirements

```
numpy, pandas, scipy, scikit-learn, xgboost, rdkit-pypi, padelpy, matplotlib, seaborn, tqdm
```

PaDEL-Descriptor (Java) required for fingerprint computation scripts (step 04). The Colab notebooks additionally install `cuml-cu12` and `cudf-cu12` for GPU acceleration.

## Documentation

See [`AGENTS.md`](AGENTS.md) for detailed pipeline documentation, known issues, and implementation notes.
