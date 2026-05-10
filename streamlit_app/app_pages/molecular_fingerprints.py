"""Molecular Fingerprints page — intercorrelation heatmap of fingerprint features."""

import os
import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

alt.data_transformers.disable_max_rows()

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")
FILTERED_DIR = os.path.join(DATA_DIR, "fingerprints_filtered")

FINGERPRINTS = [
    "maccs", "pubchem", "substructure", "klekota_roth",
    "cdk", "cdk_ext", "cdk_graphonly", "estate",
    "atompairs2d", "substructure_count", "klekota_roth_count", "atompairs2d_count",
]

FP_LABELS = {
    "maccs": "MACCS Keys",
    "pubchem": "PubChem FP",
    "substructure": "SubstructureFP",
    "klekota_roth": "Klekota-Roth FP",
    "cdk": "CDK FP",
    "cdk_ext": "CDK Extended FP",
    "cdk_graphonly": "CDK GraphOnly FP",
    "estate": "E-State FP",
    "atompairs2d": "AtomPairs2D FP",
    "substructure_count": "Substructure Count",
    "klekota_roth_count": "Klekota-Roth Count",
    "atompairs2d_count": "AtomPairs2D Count",
}


@st.cache_data
def load_fingerprint(fp_name):
    """Load a filtered fingerprint CSV (before collinearity removal)."""
    path = os.path.join(FILTERED_DIR, f"aromatase_{fp_name}_fp_filtered.csv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path)


@st.cache_data
def compute_correlation(fp_name):
    """Compute Pearson correlation matrix and cluster-ordered indices."""
    df = load_fingerprint(fp_name)
    if df is None:
        return None, None, None

    fp_cols = [c for c in df.columns if c != "molecule_chembl_id"]
    X = df[fp_cols].values.astype(np.float64)

    # Pearson correlation
    corr = np.corrcoef(X, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0)

    # Hierarchical clustering for row/column ordering
    # Convert correlation to distance: d = 1 - |r|
    dist = 1.0 - np.abs(corr)
    np.fill_diagonal(dist, 0.0)
    # Ensure symmetry and no negatives (numerical noise)
    dist = np.clip((dist + dist.T) / 2, 0, 2)
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    order = leaves_list(Z)

    return corr, order, fp_cols


# ── Page content ──
st.header("Molecular Fingerprints")
st.caption("Intercorrelation heatmap of fingerprint feature columns (before collinearity removal).")

# ── Controls ──
col_fp, col_thresh = st.columns([2, 1])
with col_fp:
    fp_choice = st.selectbox(
        "Fingerprint",
        FINGERPRINTS,
        format_func=lambda x: f"{FP_LABELS[x]}",
    )
with col_thresh:
    threshold = st.slider("Highlight threshold |r| ≥", 0.3, 0.9, 0.7, 0.05)

# ── Compute correlation ──
with st.spinner(f"Computing correlation matrix for {FP_LABELS[fp_choice]}..."):
    corr, order, fp_cols = compute_correlation(fp_choice)

if corr is None:
    st.warning(f"Filtered fingerprint file not found for {fp_choice}.")
    st.stop()

n_features = len(fp_cols)

# ── Summary metrics ──
# Extract upper triangle values (excluding diagonal)
triu_idx = np.triu_indices(n_features, k=1)
pairwise_r = corr[triu_idx]
abs_pairwise = np.abs(pairwise_r)

n_pairs = len(pairwise_r)
mean_abs_r = float(abs_pairwise.mean())
max_abs_r = float(abs_pairwise.max())
pairs_above = int((abs_pairwise >= threshold).sum())

with st.container(horizontal=True):
    st.metric("Features", f"{n_features}", border=True)
    st.metric("Total pairs", f"{n_pairs:,}", border=True)
    st.metric("Mean |r|", f"{mean_abs_r:.4f}", border=True)
    st.metric("Max |r|", f"{max_abs_r:.4f}", border=True)
    st.metric(f"Pairs |r| ≥ {threshold}", f"{pairs_above:,}", border=True)

# ── Intercorrelation heatmap ──
with st.container(border=True):
    st.subheader("Intercorrelation Heatmap")
    st.caption(
        f"Pairwise Pearson correlation of {n_features} features in "
        f"{FP_LABELS[fp_choice]}. Rows and columns are reordered by hierarchical "
        f"clustering (average linkage) to group correlated features together."
    )

    # Reorder correlation matrix by clustering
    corr_ordered = corr[np.ix_(order, order)]
    cols_ordered = [fp_cols[i] for i in order]

    # Melt for Altair
    progress = st.progress(0, text="Preparing heatmap data...")
    melted = []
    total = n_features * n_features
    for i in range(n_features):
        for j in range(n_features):
            melted.append({
                "x": cols_ordered[j],
                "y": cols_ordered[i],
                "r": float(corr_ordered[i, j]),
            })
        # Update progress per row
        progress.progress((i + 1) / n_features, text=f"Preparing heatmap data... ({(i+1)*n_features:,}/{total:,} cells)")
    melted_df = pd.DataFrame(melted)
    progress.progress(1.0, text="Rendering heatmap...")

    heatmap = (
        alt.Chart(melted_df)
        .mark_rect()
        .encode(
            x=alt.X("x:N", sort=cols_ordered, title=None,
                    axis=alt.Axis(labels=False, ticks=False, domainColor="white")),
            y=alt.Y("y:N", sort=cols_ordered, title=None,
                    axis=alt.Axis(labels=False, ticks=False, domainColor="white")),
            color=alt.Color("r:Q",
                            scale=alt.Scale(scheme="redblue", domain=[-1, 1]),
                            title="Pearson r"),
            tooltip=[
                alt.Tooltip("x:N", title="Feature 1"),
                alt.Tooltip("y:N", title="Feature 2"),
                alt.Tooltip("r:Q", title="r", format=".4f"),
            ],
        )
        .properties(
            height=600,
        )
    )
    st.altair_chart(heatmap, use_container_width=True)
    progress.empty()

# ── Correlation distribution histogram ──
with st.container(border=True):
    st.subheader("Correlation Distribution")
    st.caption("Distribution of all pairwise |r| values (upper triangle only).")

    hist_df = pd.DataFrame({"abs_r": abs_pairwise})

    hist_chart = (
        alt.Chart(hist_df)
        .mark_bar(opacity=0.85, color="steelblue")
        .encode(
            x=alt.X("abs_r:Q", bin=alt.Bin(maxbins=50), title="|Pearson r|"),
            y=alt.Y("count()", title="Count"),
        )
        .properties(height=300)
    )

    # Add threshold rule
    rule = (
        alt.Chart(pd.DataFrame({"threshold": [threshold]}))
        .mark_rule(color="red", strokeDash=[4, 4], strokeWidth=2)
        .encode(x="threshold:Q")
    )

    st.altair_chart(hist_chart + rule, use_container_width=True)

# ── Top correlated pairs table ──
with st.container(border=True):
    st.subheader(f"Correlated Pairs (|r| ≥ {threshold})")

    if pairs_above == 0:
        st.info(f"No feature pairs with |r| ≥ {threshold} found.")
    else:
        # Build pairs table
        pair_rows = []
        for idx in range(len(pairwise_r)):
            if abs_pairwise[idx] >= threshold:
                i, j = triu_idx[0][idx], triu_idx[1][idx]
                pair_rows.append({
                    "Feature 1": fp_cols[i],
                    "Feature 2": fp_cols[j],
                    "r": float(pairwise_r[idx]),
                    "|r|": float(abs_pairwise[idx]),
                })

        pairs_df = pd.DataFrame(pair_rows).sort_values("|r|", ascending=False).reset_index(drop=True)
        st.dataframe(
            pairs_df.head(200),
            use_container_width=True,
            hide_index=True,
            column_config={
                "r": st.column_config.NumberColumn(format="%.4f"),
                "|r|": st.column_config.NumberColumn(format="%.4f"),
            },
        )
        if len(pairs_df) > 200:
            st.caption(f"Showing top 200 of {len(pairs_df)} pairs.")
