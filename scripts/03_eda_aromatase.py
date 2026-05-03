#!/usr/bin/env python3
"""Exploratory Data Analysis of aromatase curated bioactivity dataset."""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# ── Paths ────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_FILE = os.path.join(PROJECT_DIR, "data", "processed", "aromatase_bioactivity_curated.csv")
OUT_DIR = os.path.join(PROJECT_DIR, "data", "eda")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load ─────────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_FILE)
print(f"Dataset: {df.shape[0]} rows × {df.shape[1]} columns")
print(f"Columns: {list(df.columns)}\n")

# ══════════════════════════════════════════════════════════════════════════
# 1. BASIC STATISTICS
# ══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("1. BASIC STATISTICS")
print("=" * 70)

print(f"\nDataset shape: {df.shape}")
print(f"Unique molecules: {df['molecule_chembl_id'].nunique()}")
print(f"Unique assays: {df['assay_chembl_id'].nunique()}")
print(f"Unique documents: {df['document_chembl_id'].nunique()}")

print("\n── Data types ──")
print(df.dtypes.to_string())

print("\n── Numeric summary ──")
numeric_cols = ["standard_value", "pchembl_value", "sd_pchembl", "n_measurements"]
print(df[numeric_cols].describe().round(3).to_string())

# ══════════════════════════════════════════════════════════════════════════
# 2. MISSING DATA
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. MISSING DATA ANALYSIS")
print("=" * 70)

missing = df.isnull().sum()
missing_pct = (missing / len(df) * 100).round(2)
missing_df = pd.DataFrame({"missing_count": missing, "missing_pct": missing_pct})
missing_df = missing_df[missing_df["missing_count"] > 0].sort_values("missing_pct", ascending=False)
if len(missing_df) > 0:
    print(missing_df.to_string())
else:
    print("No missing values in any column.")

# Check empty strings in key text columns
for col in ["molecule_pref_name", "canonical_smiles", "assay_description", "document_journal"]:
    empty = (df[col].fillna("") == "").sum()
    if empty > 0:
        print(f"  '{col}' has {empty} empty/null values ({empty/len(df)*100:.1f}%)")

# ══════════════════════════════════════════════════════════════════════════
# 3. MEASUREMENT TYPE DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. MEASUREMENT TYPE DISTRIBUTION")
print("=" * 70)

type_counts = df["standard_type"].value_counts()
print(type_counts.to_string())

print("\n── standard_relation distribution ──")
print(df["standard_relation"].value_counts().to_string())

print("\n── assay_type distribution ──")
print(df["assay_type"].value_counts().to_string())

# ══════════════════════════════════════════════════════════════════════════
# 4. pchembl_value DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. pChEMBL VALUE DISTRIBUTION")
print("=" * 70)

pvals = df["pchembl_value"].dropna()
print(f"Count: {len(pvals)}")
print(f"Mean:  {pvals.mean():.3f}")
print(f"Median:{pvals.median():.3f}")
print(f"Std:   {pvals.std():.3f}")
print(f"Min:   {pvals.min():.3f}")
print(f"Max:   {pvals.max():.3f}")
print(f"IQR:   {pvals.quantile(0.25):.3f} – {pvals.quantile(0.75):.3f}")

# Per standard_type
print("\n── pchembl_value by standard_type ──")
print(df.groupby("standard_type")["pchembl_value"].describe().round(3).to_string())

# ══════════════════════════════════════════════════════════════════════════
# 5. DEDUPLICATION STATISTICS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. DEDUPLICATION / AVERAGING STATISTICS")
print("=" * 70)

n_meas = df["n_measurements"]
print(f"Singleton records (n=1): {(n_meas == 1).sum()} ({(n_meas == 1).mean()*100:.1f}%)")
print(f"Averaged records (n>1):  {(n_meas > 1).sum()} ({(n_meas > 1).mean()*100:.1f}%)")
print(f"Max measurements averaged: {n_meas.max()}")
print(f"\nn_measurements distribution:")
print(n_meas.value_counts().sort_index().head(15).to_string())

sd_nonzero = df.loc[df["sd_pchembl"] > 0, "sd_pchembl"]
print(f"\nsd_pchembl (non-zero): n={len(sd_nonzero)}, mean={sd_nonzero.mean():.3f}, max={sd_nonzero.max():.3f}")

# ══════════════════════════════════════════════════════════════════════════
# 6. CHEMICAL SPACE (RDKit descriptors)
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("6. CHEMICAL SPACE ANALYSIS")
print("=" * 70)

smiles_list = df["canonical_smiles"].tolist()
mols = [Chem.MolFromSmiles(s) if pd.notna(s) and s != "" else None for s in smiles_list]
valid_mask = [m is not None for m in mols]
print(f"Valid SMILES: {sum(valid_mask)} / {len(mols)} ({sum(valid_mask)/len(mols)*100:.1f}%)")

# Compute descriptors
mw = [Descriptors.MolWt(m) if m else np.nan for m in mols]
logp = [Descriptors.MolLogP(m) if m else np.nan for m in mols]
tpsa = [Descriptors.TPSA(m) if m else np.nan for m in mols]
hba = [Descriptors.NumHAcceptors(m) if m else np.nan for m in mols]
hbd = [Descriptors.NumHDonors(m) if m else np.nan for m in mols]
rotb = [Descriptors.NumRotatableBonds(m) if m else np.nan for m in mols]
rings = [Descriptors.RingCount(m) if m else np.nan for m in mols]
aromatic_rings = [Descriptors.NumAromaticRings(m) if m else np.nan for m in mols]
heavy_atoms = [Descriptors.HeavyAtomCount(m) if m else np.nan for m in mols]
fsp3 = [Descriptors.FractionCSP3(m) if m else np.nan for m in mols]

df["MW"] = mw
df["LogP"] = logp
df["TPSA"] = tpsa
df["HBA"] = hba
df["HBD"] = hbd
df["RotBonds"] = rotb
df["RingCount"] = rings
df["AromaticRings"] = aromatic_rings
df["HeavyAtoms"] = heavy_atoms
df["FractionCSP3"] = fsp3

desc_cols = ["MW", "LogP", "TPSA", "HBA", "HBD", "RotBonds", "RingCount", "AromaticRings", "HeavyAtoms", "FractionCSP3"]
print("\n── Molecular descriptor statistics ──")
print(df[desc_cols].describe().round(2).to_string())

# Lipinski Rule-of-5
ro5_pass = sum(1 for i in range(len(df)) if valid_mask[i]
               and mw[i] <= 500 and logp[i] <= 5 and hba[i] <= 10 and hbd[i] <= 5)
ro5_total = sum(valid_mask)
print(f"\nLipinski Rule-of-5: {ro5_pass}/{ro5_total} pass ({ro5_pass/ro5_total*100:.1f}%)")

# Violations breakdown
v_mw = sum(1 for i in range(len(df)) if valid_mask[i] and mw[i] > 500)
v_logp = sum(1 for i in range(len(df)) if valid_mask[i] and logp[i] > 5)
v_hba = sum(1 for i in range(len(df)) if valid_mask[i] and hba[i] > 10)
v_hbd = sum(1 for i in range(len(df)) if valid_mask[i] and hbd[i] > 5)
print(f"  MW > 500:  {v_mw} ({v_mw/ro5_total*100:.1f}%)")
print(f"  LogP > 5:  {v_logp} ({v_logp/ro5_total*100:.1f}%)")
print(f"  HBA > 10:  {v_hba} ({v_hba/ro5_total*100:.1f}%)")
print(f"  HBD > 5:   {v_hbd} ({v_hbd/ro5_total*100:.1f}%)")

# ══════════════════════════════════════════════════════════════════════════
# 7. TEMPORAL ANALYSIS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("7. TEMPORAL ANALYSIS")
print("=" * 70)

year = df["document_year"].dropna().astype(int)
print(f"Publication years: {year.min()} – {year.max()}")
print(f"Peak year: {year.mode().iloc[0]} ({(year == year.mode().iloc[0]).sum()} records)")
print("\n── Records per decade ──")
decade = (year // 10 * 10)
print(decade.value_counts().sort_index().to_string())

# ══════════════════════════════════════════════════════════════════════════
# 8. POTENCY CLASSES
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("8. POTENCY CLASSIFICATION")
print("=" * 70)

def potency_class(pchembl):
    if pd.isna(pchembl):
        return "Unknown"
    if pchembl >= 8:
        return "Very potent (pChEMBL >= 8)"
    elif pchembl >= 6:
        return "Potent (6 <= pChEMBL < 8)"
    elif pchembl >= 5:
        return "Moderate (5 <= pChEMBL < 6)"
    else:
        return "Weak (pChEMBL < 5)"

df["potency_class"] = df["pchembl_value"].apply(potency_class)
print(df["potency_class"].value_counts().to_string())

# Potency by standard_type
print("\n── Potency class by standard_type ──")
ct = pd.crosstab(df["standard_type"], df["potency_class"])
print(ct.to_string())

# ══════════════════════════════════════════════════════════════════════════
# 9. TOP ASSAYS AND DOCUMENTS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("9. TOP ASSAYS AND DOCUMENTS")
print("=" * 70)

print("\n── Top 10 assays by record count ──")
top_assays = df["assay_chembl_id"].value_counts().head(10)
for aid, cnt in top_assays.items():
    desc = df.loc[df["assay_chembl_id"] == aid, "assay_description"].iloc[0]
    desc_short = desc[:80] + "..." if len(str(desc)) > 80 else desc
    print(f"  {aid}: {cnt:>4} records — {desc_short}")

print("\n── Top 10 journals ──")
print(df["document_journal"].value_counts().head(10).to_string())

# ══════════════════════════════════════════════════════════════════════════
# 10. CORRELATION: pChEMBL vs. descriptors
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("10. CORRELATION: pChEMBL vs. MOLECULAR DESCRIPTORS")
print("=" * 70)

for col in desc_cols:
    corr = df[["pchembl_value", col]].dropna().corr().iloc[0, 1]
    print(f"  pchembl_value vs {col:15s}: r = {corr:+.3f}")


# ══════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("GENERATING PLOTS...")
print("=" * 70)

# ── Figure 1: Overview dashboard ─────────────────────────────────────────
fig = plt.figure(figsize=(18, 14))
fig.suptitle("Aromatase (CYP19A1) Bioactivity Dataset — EDA Overview", fontsize=16, fontweight="bold", y=0.98)
gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

# 1a: pChEMBL distribution
ax1 = fig.add_subplot(gs[0, 0])
ax1.hist(df["pchembl_value"].dropna(), bins=50, color="#2196F3", edgecolor="white", alpha=0.85)
ax1.axvline(df["pchembl_value"].median(), color="red", linestyle="--", label=f'Median={df["pchembl_value"].median():.2f}')
ax1.set_xlabel("pChEMBL value")
ax1.set_ylabel("Count")
ax1.set_title("pChEMBL Distribution")
ax1.legend(fontsize=8)

# 1b: pChEMBL by standard_type (box-like via hist)
ax2 = fig.add_subplot(gs[0, 1])
for stype, color in [("IC50", "#2196F3"), ("Ki", "#FF9800"), ("pIC50", "#4CAF50")]:
    subset = df.loc[df["standard_type"] == stype, "pchembl_value"].dropna()
    ax2.hist(subset, bins=40, alpha=0.6, color=color, label=f"{stype} (n={len(subset)})", edgecolor="white")
ax2.set_xlabel("pChEMBL value")
ax2.set_ylabel("Count")
ax2.set_title("pChEMBL by Measurement Type")
ax2.legend(fontsize=8)

# 1c: Measurement type pie
ax3 = fig.add_subplot(gs[0, 2])
type_counts.plot.pie(ax=ax3, autopct="%1.1f%%", colors=["#2196F3", "#FF9800", "#4CAF50"],
                     startangle=90, textprops={"fontsize": 9})
ax3.set_ylabel("")
ax3.set_title("Measurement Types")

# 1d: MW distribution
ax4 = fig.add_subplot(gs[1, 0])
ax4.hist(df["MW"].dropna(), bins=50, color="#9C27B0", edgecolor="white", alpha=0.85)
ax4.axvline(500, color="red", linestyle="--", alpha=0.7, label="Ro5 cutoff (500)")
ax4.set_xlabel("Molecular Weight (Da)")
ax4.set_ylabel("Count")
ax4.set_title("Molecular Weight Distribution")
ax4.legend(fontsize=8)

# 1e: LogP distribution
ax5 = fig.add_subplot(gs[1, 1])
ax5.hist(df["LogP"].dropna(), bins=50, color="#009688", edgecolor="white", alpha=0.85)
ax5.axvline(5, color="red", linestyle="--", alpha=0.7, label="Ro5 cutoff (5)")
ax5.set_xlabel("LogP")
ax5.set_ylabel("Count")
ax5.set_title("LogP Distribution")
ax5.legend(fontsize=8)

# 1f: TPSA distribution
ax6 = fig.add_subplot(gs[1, 2])
ax6.hist(df["TPSA"].dropna(), bins=50, color="#FF5722", edgecolor="white", alpha=0.85)
ax6.axvline(140, color="red", linestyle="--", alpha=0.7, label="Veber cutoff (140)")
ax6.set_xlabel("TPSA (Å²)")
ax6.set_ylabel("Count")
ax6.set_title("TPSA Distribution")
ax6.legend(fontsize=8)

# 1g: MW vs LogP colored by potency
ax7 = fig.add_subplot(gs[2, 0])
pot_colors = {"Very potent (pChEMBL >= 8)": "#D32F2F",
              "Potent (6 <= pChEMBL < 8)": "#FF9800",
              "Moderate (5 <= pChEMBL < 6)": "#2196F3",
              "Weak (pChEMBL < 5)": "#9E9E9E"}
for label, color in pot_colors.items():
    sub = df[df["potency_class"] == label]
    ax7.scatter(sub["MW"], sub["LogP"], s=8, alpha=0.5, color=color, label=label)
ax7.set_xlabel("MW (Da)")
ax7.set_ylabel("LogP")
ax7.set_title("Chemical Space (MW vs LogP)")
ax7.legend(fontsize=7, markerscale=2)

# 1h: Records per year
ax8 = fig.add_subplot(gs[2, 1])
year_counts = df["document_year"].dropna().astype(int).value_counts().sort_index()
ax8.bar(year_counts.index, year_counts.values, color="#3F51B5", edgecolor="white", alpha=0.85)
ax8.set_xlabel("Publication Year")
ax8.set_ylabel("Records")
ax8.set_title("Records by Publication Year")
ax8.tick_params(axis="x", rotation=45)

# 1i: Potency class bar
ax9 = fig.add_subplot(gs[2, 2])
pot_order = ["Very potent (pChEMBL >= 8)", "Potent (6 <= pChEMBL < 8)",
             "Moderate (5 <= pChEMBL < 6)", "Weak (pChEMBL < 5)"]
pot_vals = [df["potency_class"].value_counts().get(p, 0) for p in pot_order]
bars = ax9.barh(pot_order, pot_vals, color=[pot_colors[p] for p in pot_order], edgecolor="white")
for bar, val in zip(bars, pot_vals):
    ax9.text(bar.get_width() + 20, bar.get_y() + bar.get_height()/2, str(val),
             va="center", fontsize=9)
ax9.set_xlabel("Count")
ax9.set_title("Potency Classification")
ax9.invert_yaxis()

fig.savefig(os.path.join(OUT_DIR, "eda_overview.png"), dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {os.path.join(OUT_DIR, 'eda_overview.png')}")

# ── Figure 2: Descriptor correlations with pChEMBL ──────────────────────
fig2, axes2 = plt.subplots(2, 3, figsize=(16, 10))
fig2.suptitle("Molecular Descriptors vs. pChEMBL Value", fontsize=14, fontweight="bold")
plot_descs = ["MW", "LogP", "TPSA", "HBA", "RotBonds", "AromaticRings"]
colors2 = ["#9C27B0", "#009688", "#FF5722", "#2196F3", "#FF9800", "#795548"]

for ax, desc, c in zip(axes2.ravel(), plot_descs, colors2):
    valid = df[["pchembl_value", desc]].dropna()
    ax.scatter(valid[desc], valid["pchembl_value"], s=6, alpha=0.3, color=c)
    corr = valid.corr().iloc[0, 1]
    ax.set_xlabel(desc)
    ax.set_ylabel("pChEMBL")
    ax.set_title(f"{desc} vs pChEMBL (r={corr:+.3f})")

fig2.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "eda_descriptors_vs_pchembl.png"), dpi=150, bbox_inches="tight")
plt.close(fig2)
print(f"  Saved: {os.path.join(OUT_DIR, 'eda_descriptors_vs_pchembl.png')}")

# ── Figure 3: Descriptor distributions ──────────────────────────────────
fig3, axes3 = plt.subplots(2, 5, figsize=(20, 8))
fig3.suptitle("Molecular Descriptor Distributions", fontsize=14, fontweight="bold")
all_desc = ["MW", "LogP", "TPSA", "HBA", "HBD", "RotBonds", "RingCount", "AromaticRings", "HeavyAtoms", "FractionCSP3"]
colors3 = ["#9C27B0", "#009688", "#FF5722", "#2196F3", "#FF9800", "#795548", "#E91E63", "#3F51B5", "#00BCD4", "#8BC34A"]

for ax, desc, c in zip(axes3.ravel(), all_desc, colors3):
    vals = df[desc].dropna()
    ax.hist(vals, bins=40, color=c, edgecolor="white", alpha=0.85)
    ax.set_title(f"{desc}\nμ={vals.mean():.1f}, σ={vals.std():.1f}", fontsize=10)
    ax.set_xlabel(desc)

fig3.tight_layout()
fig3.savefig(os.path.join(OUT_DIR, "eda_descriptor_distributions.png"), dpi=150, bbox_inches="tight")
plt.close(fig3)
print(f"  Saved: {os.path.join(OUT_DIR, 'eda_descriptor_distributions.png')}")

# ── Figure 4: n_measurements distribution ────────────────────────────────
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(12, 5))
fig4.suptitle("Measurement Averaging Statistics", fontsize=14, fontweight="bold")

n_dist = df["n_measurements"].value_counts().sort_index()
ax4a.bar(n_dist.index[:15], n_dist.values[:15], color="#3F51B5", edgecolor="white")
ax4a.set_xlabel("n_measurements")
ax4a.set_ylabel("Count")
ax4a.set_title("Distribution of Averaged Measurements")

# sd_pchembl for averaged records
sd_vals = df.loc[df["n_measurements"] > 1, "sd_pchembl"]
ax4b.hist(sd_vals, bins=40, color="#FF9800", edgecolor="white", alpha=0.85)
ax4b.set_xlabel("SD of pChEMBL")
ax4b.set_ylabel("Count")
ax4b.set_title(f"SD of pChEMBL (averaged records, n={len(sd_vals)})")

fig4.tight_layout()
fig4.savefig(os.path.join(OUT_DIR, "eda_averaging_stats.png"), dpi=150, bbox_inches="tight")
plt.close(fig4)
print(f"  Saved: {os.path.join(OUT_DIR, 'eda_averaging_stats.png')}")

print("\n✓ EDA complete. All plots saved to:", OUT_DIR)
