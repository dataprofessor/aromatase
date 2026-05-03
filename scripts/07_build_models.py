"""
Build and evaluate 16 regression models predicting pChEMBL values from fingerprints.

Uses the random 80/20 split (seed=42) and evaluates with R², RMSE, and MAE
on three sets: training, 10-fold cross-validation, and external test.

Results are saved incrementally into 3 separate CSV files per run:
  - {fp}_{split}_train.csv
  - {fp}_{split}_cv.csv
  - {fp}_{split}_test.csv

Supports batch mode to run a subset of models at a time (avoids crashes).

Algorithms (one or more representative per ML family):
  1. Ridge Regression           (Linear)
  2. Lasso Regression           (Linear)
  3. ElasticNet                 (Linear)
  4. k-Nearest Neighbors        (Instance-based)
  5. SVR (RBF kernel)           (Support Vector)
  6. SVR (linear kernel)        (Support Vector)
  7. Decision Tree              (Decision Tree)
  8. Random Forest              (Bagging)
  9. Extra Trees                (Bagging)
  10. Gradient Boosting          (Boosting)
  11. XGBoost                    (Boosting)
  12. AdaBoost                   (Boosting)
  13. MLP Regressor              (Neural Network)
  14. Gaussian Process           (Gaussian Process)
  15. Bayesian Ridge             (Bayesian)
  16. PLS Regression             (Partial Least Squares)

Usage:
    python scripts/07_build_models.py                    # all 16 models (resumable)
    python scripts/07_build_models.py --batch 1          # only model #1
    python scripts/07_build_models.py --batch 5          # only model #5
    python scripts/07_build_models.py --batch 1-5        # models 1-5
    python scripts/07_build_models.py --fp maccs --split random   # (default)
"""

import argparse
import csv
import math
import os
import time

import numpy as np
from sklearn.model_selection import cross_val_predict, KFold

_DIR = os.path.dirname(__file__)
SPLITS_DIR = os.path.join(_DIR, "..", "data", "splits")
CURATED = os.path.join(_DIR, "..", "data", "processed", "aromatase_bioactivity_curated.csv")
OUTPUT_DIR = os.path.join(_DIR, "..", "data", "models")

RANDOM_STATE = 42
N_FOLDS = 10

CSV_FIELDS = ["model", "r2", "rmse", "mae", "time_s"]


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    return rows, fieldnames


def prepare_data(fp_path, curated_rows):
    """Load fingerprint split and join with pchembl_value. Returns X, y arrays."""
    fp_rows, fp_fields = load_csv(fp_path)
    fp_cols = [c for c in fp_fields if c != "molecule_chembl_id"]

    # Build lookup: molecule_chembl_id -> pchembl_value
    target_lookup = {}
    for r in curated_rows:
        mid = r["molecule_chembl_id"]
        pval = r["pchembl_value"]
        if pval != "":
            if mid not in target_lookup:
                target_lookup[mid] = float(pval)

    # Join and filter
    X_list = []
    y_list = []
    for r in fp_rows:
        mid = r["molecule_chembl_id"]
        if mid in target_lookup:
            features = [float(r[c]) if r[c] != "" else 0.0 for c in fp_cols]
            X_list.append(features)
            y_list.append(target_lookup[mid])

    return np.array(X_list), np.array(y_list), fp_cols


def compute_metrics(y_true, y_pred):
    """Compute R², RMSE, MAE."""
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    rmse = math.sqrt(np.mean((y_true - y_pred) ** 2))
    mae = np.mean(np.abs(y_true - y_pred))
    return r2, rmse, mae


def get_models():
    """Return list of (name, model) tuples for all 16 algorithms."""
    from sklearn.linear_model import Ridge, Lasso, ElasticNet, BayesianRidge
    from sklearn.neighbors import KNeighborsRegressor
    from sklearn.svm import SVR
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.ensemble import (
        RandomForestRegressor,
        ExtraTreesRegressor,
        GradientBoostingRegressor,
        AdaBoostRegressor,
    )
    from sklearn.neural_network import MLPRegressor
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.cross_decomposition import PLSRegression
    from xgboost import XGBRegressor

    models = [
        ("Ridge Regression", Ridge(alpha=1.0, solver="svd")),
        ("Lasso Regression", Lasso(alpha=0.1, random_state=RANDOM_STATE)),
        ("ElasticNet", ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=RANDOM_STATE)),
        ("k-Nearest Neighbors", KNeighborsRegressor(n_neighbors=5)),
        ("SVR (RBF)", SVR(kernel="rbf", C=1.0, epsilon=0.1)),
        ("SVR (Linear)", SVR(kernel="linear", C=1.0, epsilon=0.1)),
        ("Decision Tree", DecisionTreeRegressor(random_state=RANDOM_STATE)),
        ("Random Forest", RandomForestRegressor(
            n_estimators=500, random_state=RANDOM_STATE, n_jobs=2
        )),
        ("Extra Trees", ExtraTreesRegressor(
            n_estimators=500, random_state=RANDOM_STATE, n_jobs=2
        )),
        ("Gradient Boosting", GradientBoostingRegressor(
            n_estimators=100, random_state=RANDOM_STATE
        )),
        ("XGBoost", XGBRegressor(
            n_estimators=500, random_state=RANDOM_STATE, n_jobs=2,
            verbosity=0
        )),
        ("AdaBoost", AdaBoostRegressor(
            n_estimators=500, random_state=RANDOM_STATE
        )),
        ("MLP Regressor", MLPRegressor(
            hidden_layer_sizes=(128, 64), max_iter=1000,
            random_state=RANDOM_STATE, early_stopping=True
        )),
        ("Gaussian Process", GaussianProcessRegressor(
            random_state=RANDOM_STATE, n_restarts_optimizer=2
        )),
        ("Bayesian Ridge", BayesianRidge()),
        ("PLS Regression", PLSRegression(n_components=10)),
    ]
    return models


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


def parse_batch(batch_str):
    """Parse batch range like '1-5' or '6-10'. Returns (start, end) 1-indexed."""
    parts = batch_str.split("-")
    if len(parts) == 2:
        return int(parts[0]), int(parts[1])
    elif len(parts) == 1:
        v = int(parts[0])
        return v, v
    else:
        raise ValueError(f"Invalid batch format: {batch_str}. Use e.g. '1-5' or '6-10'.")


def main():
    parser = argparse.ArgumentParser(description="Build ML regression models for aromatase.")
    parser.add_argument("--fp", default="maccs", help="Fingerprint name (default: maccs)")
    parser.add_argument("--split", default="random", choices=["random", "kennard_stone"],
                        help="Split strategy (default: random)")
    parser.add_argument("--batch", default=None,
                        help="Model range to run, e.g. '1-5', '6-10', '11-16' (default: all)")
    args = parser.parse_args()

    # Resolve paths
    fp_stem = f"aromatase_{args.fp}_fp"
    train_path = os.path.join(SPLITS_DIR, f"{fp_stem}_{args.split}_train.csv")
    test_path = os.path.join(SPLITS_DIR, f"{fp_stem}_{args.split}_test.csv")

    if not os.path.exists(train_path):
        print(f"ERROR: Train file not found: {train_path}")
        print(f"Available splits:")
        for f in sorted(os.listdir(SPLITS_DIR)):
            print(f"  {f}")
        return

    # Output file paths
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    prefix = f"{args.fp}_{args.split}_split"
    out_train = os.path.join(OUTPUT_DIR, f"{prefix}_train.csv")
    out_cv = os.path.join(OUTPUT_DIR, f"{prefix}_cv.csv")
    out_test = os.path.join(OUTPUT_DIR, f"{prefix}_test.csv")

    print("=" * 70)
    print(f"ML Regression Models: {args.fp.upper()} Fingerprint ({args.split} split, seed=42)")
    print("=" * 70)
    print(f"\n  Output files:")
    print(f"    Train: {out_train}")
    print(f"    CV:    {out_cv}")
    print(f"    Test:  {out_test}")

    # Load curated data for target values
    curated_rows, _ = load_csv(CURATED)
    print(f"\n  Loaded curated data: {len(curated_rows)} rows")

    # Prepare train/test
    print(f"\n  Preparing training data...")
    X_train, y_train, fp_cols = prepare_data(train_path, curated_rows)
    print(f"    X_train: {X_train.shape}")
    print(f"    y_train: {y_train.shape} (mean={y_train.mean():.2f}, std={y_train.std():.2f})")

    print(f"\n  Preparing test data...")
    X_test, y_test, _ = prepare_data(test_path, curated_rows)
    print(f"    X_test: {X_test.shape}")
    print(f"    y_test: {y_test.shape} (mean={y_test.mean():.2f}, std={y_test.std():.2f})")

    # Get models and apply batch filter
    all_models = get_models()
    if args.batch:
        start, end = parse_batch(args.batch)
        models = all_models[start - 1:end]
        print(f"\n  Batch mode: models {start}-{end} ({len(models)} models)")
    else:
        models = all_models
        print(f"\n  Running all {len(models)} models")

    # Check which models already completed (skip them)
    completed = load_completed_models(out_train)
    models_to_run = [(name, m) for name, m in models if name not in completed]
    if len(models_to_run) < len(models):
        skipped = len(models) - len(models_to_run)
        print(f"  Skipping {skipped} already-completed model(s)")

    if not models_to_run:
        print(f"\n  All models in this batch already completed. Nothing to do.")
        return

    # Train and evaluate — one model at a time (checkpoint after each)
    cv = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    total = len(models_to_run)

    print(f"\n{'='*70}")
    print(f"Training {total} models (with {N_FOLDS}-fold CV)...")
    print(f"Each model is checkpointed immediately after completion.")
    print(f"{'='*70}")

    for i, (name, model) in enumerate(models_to_run, 1):
        print(f"\n  [{i:2d}/{total}] {name}...")
        print(f"        Training...", end=" ", flush=True)
        start = time.time()

        # Fit on full training set
        model.fit(X_train, y_train)

        # Predict on train and test
        y_pred_train = model.predict(X_train)
        y_pred_test = model.predict(X_test)

        # PLS returns 2D array
        if y_pred_train.ndim > 1:
            y_pred_train = y_pred_train.ravel()
        if y_pred_test.ndim > 1:
            y_pred_test = y_pred_test.ravel()

        # 10-fold cross-validation on training set
        from sklearn.base import clone
        model_cv = clone(model)
        y_pred_cv = cross_val_predict(model_cv, X_train, y_train, cv=cv, n_jobs=1)
        if y_pred_cv.ndim > 1:
            y_pred_cv = y_pred_cv.ravel()

        elapsed = time.time() - start

        # Compute metrics
        r2_train, rmse_train, mae_train = compute_metrics(y_train, y_pred_train)
        r2_cv, rmse_cv, mae_cv = compute_metrics(y_train, y_pred_cv)
        r2_test, rmse_test, mae_test = compute_metrics(y_test, y_pred_test)

        # Checkpoint: save immediately to 3 separate files
        append_row(out_train, {"model": name, "r2": r2_train, "rmse": rmse_train, "mae": mae_train, "time_s": elapsed})
        append_row(out_cv, {"model": name, "r2": r2_cv, "rmse": rmse_cv, "mae": mae_cv, "time_s": elapsed})
        append_row(out_test, {"model": name, "r2": r2_test, "rmse": rmse_test, "mae": mae_test, "time_s": elapsed})

        print(f"done ({elapsed:.1f}s) [SAVED]")
        print(f"        Train: R²={r2_train:.4f}, RMSE={rmse_train:.4f}")
        print(f"        CV:    R²={r2_cv:.4f}, RMSE={rmse_cv:.4f}")
        print(f"        Test:  R²={r2_test:.4f}, RMSE={rmse_test:.4f}")
        print(f"        ✓ Checkpoint {len(completed)+i}/{len(completed)+total} saved")

    # Print summary of all results in the files (including previous batches)
    print(f"\n{'='*70}")
    print(f"Batch complete. Summary of all results so far:")
    print(f"{'='*70}")

    if os.path.exists(out_test):
        rows, _ = load_csv(out_test)
        rows.sort(key=lambda x: float(x["r2"]), reverse=True)
        print(f"\n  {'#':<3} {'Model':<25} {'R²(test)':<10} {'RMSE(test)':<12}")
        print(f"  {'-'*3} {'-'*25} {'-'*10} {'-'*12}")
        for i, r in enumerate(rows, 1):
            print(f"  {i:<3} {r['model']:<25} {float(r['r2']):<10.4f} {float(r['rmse']):<12.4f}")

    print(f"\n  Results saved to:")
    print(f"    {out_train}")
    print(f"    {out_cv}")
    print(f"    {out_test}")
    print(f"\n{'='*70}")


if __name__ == "__main__":
    main()
