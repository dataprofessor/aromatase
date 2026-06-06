"""
Build and evaluate 14 classification models predicting bioactivity class from fingerprints.

3-tier training with per-model zip file generation:
  - Tier 1 (Fast): Ridge, LogReg, DecTree, KNN, NB, LDA
  - Tier 2 (Medium): Gradient Boosting, XGBoost, AdaBoost, MLP
  - Tier 3 (Slow): SVC (RBF/Linear), Random Forest

Each model's results are saved to a separate zip file immediately after completion.

Uses the random or kennard_stone 80/20 split (seed=42) and evaluates with
Accuracy, Balanced Accuracy, F1 (weighted), and MCC on three sets:
training, 10-fold stratified cross-validation, and external test.

Two modes:
  - Default: all 14 models without class weighting
  - Balanced: 7 models with class_weight='balanced'

Usage:
    python scripts/08c_build_classification_models.py                         # all FPs, random split
    python scripts/08c_build_classification_models.py --fp maccs              # single FP
    python scripts/08c_build_classification_models.py --split kennard_stone   # KS split
    python scripts/08c_build_classification_models.py --tier 1                # only fast models
    python scripts/08c_build_classification_models.py --tier 1-2              # fast + medium
    python scripts/08c_build_classification_models.py --mode balanced         # only balanced mode
"""

import argparse
import csv
import math
import os
import time
import zipfile

import numpy as np
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, matthews_corrcoef
from sklearn.base import clone

_DIR = os.path.dirname(os.path.abspath(__file__))
SPLITS_DIR = os.path.join(_DIR, "..", "data", "splits")
CURATED = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_curated.csv")
OUTPUT_DIR_DEFAULT = os.path.join(_DIR, "..", "data", "models_classification")
OUTPUT_DIR_BALANCED = os.path.join(_DIR, "..", "data", "models_classification_balanced")
ZIP_DIR = os.path.join(_DIR, "..", "data", "models_classification_zips")

RANDOM_STATE = 42
N_FOLDS = 10

CSV_FIELDS = ["model", "accuracy", "balanced_accuracy", "f1_weighted", "mcc", "time_s"]

# Class definitions
CLASS_LABELS = ["inactive", "intermediate", "active"]
CLASS_ENCODING = {"inactive": 0, "intermediate": 1, "active": 2}

# Bioactivity class thresholds
# Active: pChEMBL > 7, Intermediate: 6 <= pChEMBL <= 7, Inactive: pChEMBL < 6

FINGERPRINTS = [
    "maccs", "pubchem", "substructure", "klekota_roth",
    "cdk", "cdk_ext", "cdk_graphonly", "estate",
    "atompairs2d", "substructure_count", "klekota_roth_count", "atompairs2d_count",
]

# 3-tier model split
TIER1_FAST = [
    "Ridge Classifier", "Logistic Regression (L1)", "Logistic Regression (ElasticNet)",
    "Decision Tree", "k-Nearest Neighbors", "Gaussian Naive Bayes",
    "Linear Discriminant Analysis",
]
TIER2_MEDIUM = ["Gradient Boosting", "XGBoost", "AdaBoost", "MLP Classifier"]
TIER3_SLOW = ["SVC (RBF)", "SVC (Linear)", "Random Forest"]


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    return rows, fieldnames


def pchembl_to_class(pchembl):
    """Convert pChEMBL value to bioactivity class."""
    if pchembl > 7:
        return 2  # active
    elif pchembl >= 6:
        return 1  # intermediate
    else:
        return 0  # inactive


def prepare_data(fp_path, curated_rows):
    """Load fingerprint split and join with bioactivity class. Returns X, y arrays."""
    fp_rows, fp_fields = load_csv(fp_path)
    fp_cols = [c for c in fp_fields if c != "molecule_chembl_id"]

    # Build lookup: molecule_chembl_id -> class label
    target_lookup = {}
    for r in curated_rows:
        mid = r["molecule_chembl_id"]
        pval = r["pchembl_value"]
        if pval != "":
            if mid not in target_lookup:
                target_lookup[mid] = pchembl_to_class(float(pval))

    # Join and filter
    X_list = []
    y_list = []
    for r in fp_rows:
        mid = r["molecule_chembl_id"]
        if mid in target_lookup:
            features = [float(r[c]) if r[c] != "" else 0.0 for c in fp_cols]
            X_list.append(features)
            y_list.append(target_lookup[mid])

    return np.array(X_list, dtype=np.float32), np.array(y_list, dtype=np.int32), fp_cols


def compute_metrics(y_true, y_pred):
    """Compute classification metrics."""
    acc = accuracy_score(y_true, y_pred)
    bal_acc = balanced_accuracy_score(y_true, y_pred)
    f1_w = f1_score(y_true, y_pred, average="weighted")
    mcc = matthews_corrcoef(y_true, y_pred)
    return acc, bal_acc, f1_w, mcc


def get_models(balanced=False):
    """Return list of (name, model) tuples for all 14 classification algorithms."""
    from sklearn.linear_model import RidgeClassifier, LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.svm import SVC
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import (
        RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier,
    )
    from sklearn.neural_network import MLPClassifier
    from sklearn.naive_bayes import GaussianNB
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    try:
        from xgboost import XGBClassifier
        HAS_XGBOOST = True
    except Exception:
        HAS_XGBOOST = False

    cw = "balanced" if balanced else None

    # Models that support class_weight
    balanced_models = [
        ("Ridge Classifier", RidgeClassifier(alpha=1.0, class_weight=cw)),
        ("Logistic Regression (L1)",
         LogisticRegression(penalty="l1", solver="saga", C=1.0,
                            max_iter=1000, random_state=RANDOM_STATE, class_weight=cw)),
        ("Logistic Regression (ElasticNet)",
         LogisticRegression(penalty="elasticnet", solver="saga", C=1.0,
                            l1_ratio=0.5, max_iter=1000, random_state=RANDOM_STATE, class_weight=cw)),
        ("SVC (RBF)", SVC(kernel="rbf", C=1.0, random_state=RANDOM_STATE, class_weight=cw)),
        ("SVC (Linear)", SVC(kernel="linear", C=1.0, random_state=RANDOM_STATE, class_weight=cw)),
        ("Decision Tree", DecisionTreeClassifier(random_state=RANDOM_STATE, class_weight=cw)),
        ("Random Forest", RandomForestClassifier(
            n_estimators=500, random_state=RANDOM_STATE, class_weight=cw, n_jobs=-1)),
    ]

    # Models that do NOT support class_weight
    other_models = [
        ("k-Nearest Neighbors", KNeighborsClassifier(n_neighbors=5)),
        ("Gradient Boosting", GradientBoostingClassifier(
            n_estimators=100, random_state=RANDOM_STATE)),
    ]
    if HAS_XGBOOST:
        other_models.append(("XGBoost", XGBClassifier(
            n_estimators=500, random_state=RANDOM_STATE, n_jobs=-1,
            verbosity=0, objective="multi:softprob", num_class=3,
            eval_metric="mlogloss")))
    other_models += [
        ("AdaBoost", AdaBoostClassifier(n_estimators=500, random_state=RANDOM_STATE)),
        ("MLP Classifier", MLPClassifier(
            hidden_layer_sizes=(128, 64), max_iter=1000,
            random_state=RANDOM_STATE, early_stopping=True)),
        ("Gaussian Naive Bayes", GaussianNB()),
        ("Linear Discriminant Analysis", LinearDiscriminantAnalysis()),
    ]

    if balanced:
        return balanced_models
    else:
        return balanced_models + other_models


def append_row(path, row_dict):
    """Append a single row to a CSV file, creating it with header if needed."""
    file_exists = os.path.exists(path) and os.path.getsize(path) > 0
    with open(path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row_dict)


def load_completed_models(path):
    """Return set of model names already in a results file."""
    if not os.path.exists(path):
        return set()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return {r["model"] for r in reader}


def zip_model_results(model_name, fp_name, split_name, output_dir, mode_label):
    """Create a zip file for a single model's results across the current FP/split."""
    # Sanitize model name for filename
    safe_name = model_name.lower().replace(" ", "_").replace("(", "").replace(")", "")
    zip_name = f"08c_{fp_name}_{split_name}_{mode_label}_{safe_name}.zip"
    zip_path = os.path.join(ZIP_DIR, zip_name)

    prefix = f"{fp_name}_{split_name}_split"
    files_to_zip = [
        os.path.join(output_dir, f"{prefix}_train.csv"),
        os.path.join(output_dir, f"{prefix}_cv.csv"),
        os.path.join(output_dir, f"{prefix}_test.csv"),
    ]

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for fpath in files_to_zip:
            if os.path.exists(fpath):
                zf.write(fpath, os.path.basename(fpath))

    return zip_path


def parse_tier(tier_str):
    """Parse tier range like '1', '2', '1-3'. Returns list of tier numbers."""
    if "-" in tier_str:
        start, end = tier_str.split("-")
        return list(range(int(start), int(end) + 1))
    return [int(tier_str)]


def get_tier_models(tier_nums):
    """Return list of model names for specified tiers."""
    models = []
    if 1 in tier_nums:
        models.extend(TIER1_FAST)
    if 2 in tier_nums:
        models.extend(TIER2_MEDIUM)
    if 3 in tier_nums:
        models.extend(TIER3_SLOW)
    return models


def main():
    parser = argparse.ArgumentParser(
        description="Build classification models for aromatase (3-tier, per-model zip)."
    )
    parser.add_argument("--fp", default=None,
                        help="Fingerprint name (default: all 12)")
    parser.add_argument("--split", default="random", choices=["random", "kennard_stone"],
                        help="Split strategy (default: random)")
    parser.add_argument("--tier", default="1-3",
                        help="Tier range: '1', '2', '3', '1-2', '1-3' (default: 1-3)")
    parser.add_argument("--mode", default="both", choices=["default", "balanced", "both"],
                        help="Training mode (default: both)")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR_DEFAULT, exist_ok=True)
    os.makedirs(OUTPUT_DIR_BALANCED, exist_ok=True)
    os.makedirs(ZIP_DIR, exist_ok=True)

    # Determine fingerprints to process
    fps = [args.fp] if args.fp else FINGERPRINTS

    # Determine tiers
    tier_nums = parse_tier(args.tier)
    tier_model_names = get_tier_models(tier_nums)
    tier_label = f"Tier {args.tier}"

    # Determine modes
    modes = []
    if args.mode in ("default", "both"):
        modes.append(("default", OUTPUT_DIR_DEFAULT, False))
    if args.mode in ("balanced", "both"):
        modes.append(("balanced", OUTPUT_DIR_BALANCED, True))

    print("=" * 70)
    print(f"Classification Models (08c) — {tier_label}")
    print("=" * 70)
    print(f"  Fingerprints: {', '.join(fps)}")
    print(f"  Split: {args.split}")
    print(f"  Tiers: {tier_nums} ({len(tier_model_names)} models)")
    print(f"  Modes: {[m[0] for m in modes]}")
    print(f"  Zip output: {ZIP_DIR}")

    # Load curated data
    curated_rows, _ = load_csv(CURATED)
    print(f"\n  Loaded curated data: {len(curated_rows)} rows")

    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    total_zips = 0

    for fp_name in fps:
        fp_stem = f"aromatase_{fp_name}_fp"
        train_path = os.path.join(SPLITS_DIR, f"{fp_stem}_{args.split}_train.csv")
        test_path = os.path.join(SPLITS_DIR, f"{fp_stem}_{args.split}_test.csv")

        if not os.path.exists(train_path):
            print(f"\n  WARNING: {train_path} not found, skipping {fp_name}")
            continue

        # Load data
        print(f"\n{'='*70}")
        print(f"  Fingerprint: {fp_name} ({args.split} split)")
        print(f"{'='*70}")

        X_train, y_train, fp_cols = prepare_data(train_path, curated_rows)
        X_test, y_test, _ = prepare_data(test_path, curated_rows)
        print(f"    X_train: {X_train.shape}, X_test: {X_test.shape}")

        # Class distribution
        unique, counts = np.unique(y_train, return_counts=True)
        for u, c in zip(unique, counts):
            print(f"    Class {CLASS_LABELS[u]}: {c} ({c/len(y_train)*100:.1f}%)")

        for mode_label, output_dir, is_balanced in modes:
            all_models = get_models(balanced=is_balanced)
            # Filter to tier models
            models_to_use = [(n, m) for n, m in all_models if n in tier_model_names]

            if not models_to_use:
                continue

            # Check completed
            prefix = f"{fp_name}_{args.split}_split"
            out_train = os.path.join(output_dir, f"{prefix}_train.csv")
            out_cv = os.path.join(output_dir, f"{prefix}_cv.csv")
            out_test = os.path.join(output_dir, f"{prefix}_test.csv")

            completed = load_completed_models(out_train)
            models_to_run = [(n, m) for n, m in models_to_use if n not in completed]

            if not models_to_run:
                print(f"\n    [{mode_label}] All {len(models_to_use)} models already completed.")
                continue

            print(f"\n    [{mode_label}] Running {len(models_to_run)}/{len(models_to_use)} models...")

            for i, (name, model) in enumerate(models_to_run, 1):
                print(f"\n      [{i}/{len(models_to_run)}] {name}...", end=" ", flush=True)
                start = time.time()

                try:
                    # Fit
                    model_instance = clone(model) if name != "XGBoost" else model
                    model_instance.fit(X_train, y_train)

                    # Predict
                    y_pred_train = model_instance.predict(X_train)
                    y_pred_test = model_instance.predict(X_test)

                    # Cross-validation (manual for XGBoost due to sklearn tag incompatibility)
                    if name == "XGBoost":
                        y_pred_cv = np.zeros_like(y_train)
                        for tr_idx, val_idx in cv.split(X_train, y_train):
                            m_cv = type(model)(**model.get_params())
                            m_cv.fit(X_train[tr_idx], y_train[tr_idx])
                            y_pred_cv[val_idx] = m_cv.predict(X_train[val_idx])
                    else:
                        model_cv = clone(model)
                        y_pred_cv = cross_val_predict(model_cv, X_train, y_train, cv=cv, n_jobs=-1)

                    elapsed = time.time() - start

                    # Metrics
                    acc_tr, bal_tr, f1_tr, mcc_tr = compute_metrics(y_train, y_pred_train)
                    acc_cv, bal_cv, f1_cv, mcc_cv = compute_metrics(y_train, y_pred_cv)
                    acc_te, bal_te, f1_te, mcc_te = compute_metrics(y_test, y_pred_test)

                    # Save results
                    append_row(out_train, {"model": name, "accuracy": acc_tr,
                                           "balanced_accuracy": bal_tr, "f1_weighted": f1_tr,
                                           "mcc": mcc_tr, "time_s": elapsed})
                    append_row(out_cv, {"model": name, "accuracy": acc_cv,
                                        "balanced_accuracy": bal_cv, "f1_weighted": f1_cv,
                                        "mcc": mcc_cv, "time_s": elapsed})
                    append_row(out_test, {"model": name, "accuracy": acc_te,
                                          "balanced_accuracy": bal_te, "f1_weighted": f1_te,
                                          "mcc": mcc_te, "time_s": elapsed})

                    # Zip this model's results
                    zip_path = zip_model_results(name, fp_name, args.split, output_dir, mode_label)
                    total_zips += 1

                    print(f"done ({elapsed:.1f}s)")
                    print(f"        BalAcc={bal_te:.4f} F1={f1_te:.4f} MCC={mcc_te:.4f}")
                    print(f"        -> {os.path.basename(zip_path)}")

                except Exception as e:
                    print(f"ERROR: {e}")
                    continue

    # Final summary
    print(f"\n{'='*70}")
    print(f"COMPLETE")
    print(f"{'='*70}")
    print(f"  Total zip files generated: {total_zips}")
    print(f"  Zip directory: {ZIP_DIR}")
    print(f"  Results directories:")
    print(f"    Default:  {OUTPUT_DIR_DEFAULT}")
    print(f"    Balanced: {OUTPUT_DIR_BALANCED}")


if __name__ == "__main__":
    main()
