# Aromatase (CYP19A1) QSAR Modeling

Bioactivity prediction for human **Aromatase** (CYP19A1) using molecular fingerprints and machine learning. Data sourced from ChEMBL (target: CHEMBL1978).

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_fetch_aromatase_bioactivity.py` | Fetch bioactivity data from ChEMBL API |
| 2 | `02_clean_aromatase_bioactivity.py` | Deduplicate, average replicates, add InChIKeys |
| 3 | `03_eda_aromatase.py` | Exploratory data analysis |
| 4 | `04_compute_*.py` | Compute 12 molecular fingerprint types |
| 5 | `05_remove_near_constant_fp.py` | Remove low-variance features (SD < 0.1) |
| 6 | `06_split_data.py` | Train/test split (random + Kennard-Stone) |
| 7 | `07_build_models.py` | Train 14 regression models per fingerprint |
| 7 | `07_build_models_colab.ipynb` | GPU-accelerated model building (Colab) |

## Fingerprints (12 types)

MACCS (166-bit), PubChem (881-bit), Substructure (307-bit), Klekota-Roth (4860-bit), CDK (1024-bit), CDK Extended (1024-bit), CDK GraphOnly (1024-bit), E-State (79-bit), AtomPairs2D (780-bit), plus count versions of Substructure, Klekota-Roth, and AtomPairs2D.

## Colab Notebook

GPU-accelerated model building notebook for Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dataprofessor/aromatase/blob/main/notebooks/07_build_models_colab.ipynb)

Uses cuML (GPU) for Ridge, Lasso, ElasticNet, KNN, SVR, Random Forest and XGBoost with CUDA. Runs all 12 fingerprints x 2 splits x 14 models in ~25 min on G4 GPU.

## Dataset

- **3,774 records** (516 Ki, 3,211 IC50, 47 pIC50)
- **Target variable**: pChEMBL value (-log10 molar potency)
- **Train/test**: 80/20 split at molecule level (seed=42)

## Requirements

```
numpy, pandas, scipy, scikit-learn, xgboost, rdkit
```

PaDEL-Descriptor (Java) required for fingerprint computation steps only.
