"""Shared data loading and descriptor computation for the Streamlit EDA app."""

import os
import pandas as pd
import numpy as np
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors


DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
CURATED_FILE = os.path.join(DATA_DIR, "processed", "aromatase_bioactivity_curated.csv")


@st.cache_data
def load_curated_data():
    """Load the curated bioactivity dataset."""
    return pd.read_csv(CURATED_FILE)


@st.cache_data
def compute_descriptors(df):
    """Compute RDKit molecular descriptors and potency classes."""
    smiles_list = df["canonical_smiles"].tolist()
    mols = [Chem.MolFromSmiles(s) if pd.notna(s) and s != "" else None for s in smiles_list]

    desc_funcs = {
        "MW": Descriptors.MolWt,
        "LogP": Descriptors.MolLogP,
        "TPSA": Descriptors.TPSA,
        "HBA": Descriptors.NumHAcceptors,
        "HBD": Descriptors.NumHDonors,
        "RotBonds": Descriptors.NumRotatableBonds,
        "RingCount": Descriptors.RingCount,
        "AromaticRings": Descriptors.NumAromaticRings,
        "HeavyAtoms": Descriptors.HeavyAtomCount,
        "FractionCSP3": Descriptors.FractionCSP3,
    }

    out = df.copy()
    for name, func in desc_funcs.items():
        out[name] = [func(m) if m else np.nan for m in mols]

    out["potency_class"] = out["pchembl_value"].apply(_potency_class)
    return out


DESC_COLS = [
    "MW", "LogP", "TPSA", "HBA", "HBD",
    "RotBonds", "RingCount", "AromaticRings", "HeavyAtoms", "FractionCSP3",
]

POTENCY_ORDER = [
    "Very potent (pChEMBL >= 8)",
    "Potent (6 <= pChEMBL < 8)",
    "Moderate (5 <= pChEMBL < 6)",
    "Weak (pChEMBL < 5)",
]

POTENCY_COLORS = {
    "Very potent (pChEMBL >= 8)": "#D32F2F",
    "Potent (6 <= pChEMBL < 8)": "#FF9800",
    "Moderate (5 <= pChEMBL < 6)": "#2196F3",
    "Weak (pChEMBL < 5)": "#9E9E9E",
}


def _potency_class(pchembl):
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
