"""Chemical space page — molecular descriptors, drug-likeness, PCA & t-SNE."""

import os
import streamlit as st
import pandas as pd
import numpy as np
import altair as alt
from data_utils import (
    load_curated_data, compute_descriptors,
    DESC_COLS, POTENCY_ORDER, POTENCY_COLORS,
)

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")
SPLITS_DIR = os.path.join(DATA_DIR, "splits")
CURATED_PATH = os.path.join(DATA_DIR, "processed", "aromatase_bioactivity_curated.csv")

df = compute_descriptors(load_curated_data())

st.header("Chemical space analysis")

tab_lipinski, tab_pca = st.tabs(["Lipinski Descriptors", "PCA & t-SNE"])

# ═══════════════════════════════════════════════════════════════════════════
# TAB 1: Lipinski Descriptors
# ═══════════════════════════════════════════════════════════════════════════
with tab_lipinski:
    # ── Lipinski Rule-of-5 KPIs ──
    ro5_mask = (df["MW"] <= 500) & (df["LogP"] <= 5) & (df["HBA"] <= 10) & (df["HBD"] <= 5)

    with st.container(horizontal=True):
        st.metric("Ro5 compliant", f"{ro5_mask.sum():,} ({ro5_mask.mean()*100:.1f}%)", border=True)
        st.metric("MW > 500", f"{(df['MW'] > 500).sum()} ({(df['MW'] > 500).mean()*100:.1f}%)", border=True)
        st.metric("LogP > 5", f"{(df['LogP'] > 5).sum()} ({(df['LogP'] > 5).mean()*100:.1f}%)", border=True)
        st.metric("HBA > 10", f"{(df['HBA'] > 10).sum()} ({(df['HBA'] > 10).mean()*100:.1f}%)", border=True)
        st.metric("HBD > 5", f"{(df['HBD'] > 5).sum()} ({(df['HBD'] > 5).mean()*100:.1f}%)", border=True)

    # ── Chemical space scatter ──
    with st.container(border=True):
        st.subheader("Chemical space explorer")
        col_ctrl1, col_ctrl2, col_ctrl3 = st.columns(3)
        x_axis = col_ctrl1.selectbox("X axis", DESC_COLS, index=0)
        y_axis = col_ctrl2.selectbox("Y axis", DESC_COLS, index=1)
        color_by = col_ctrl3.selectbox("Color by", ["potency_class", "standard_type"])

        plot_df = df.dropna(subset=[x_axis, y_axis, "pchembl_value"]) if color_by == "potency_class" else df.dropna(subset=[x_axis, y_axis])

        if color_by == "potency_class":
            plot_df = plot_df[plot_df["potency_class"] != "Unknown"]
            color_enc = alt.Color("potency_class:N", scale=alt.Scale(
                domain=POTENCY_ORDER,
                range=[POTENCY_COLORS[p] for p in POTENCY_ORDER],
            ), title="Potency")
        else:
            color_enc = alt.Color("standard_type:N", scale=alt.Scale(
                domain=["IC50", "Ki", "pIC50"],
                range=["#2196F3", "#FF9800", "#4CAF50"],
            ), title="Type")

        chart = (
            alt.Chart(plot_df)
            .mark_circle(size=20, opacity=0.5)
            .encode(
                x=alt.X(f"{x_axis}:Q", title=x_axis),
                y=alt.Y(f"{y_axis}:Q", title=y_axis),
                color=color_enc,
                tooltip=["molecule_chembl_id", x_axis, y_axis, "pchembl_value", "standard_type"],
            )
            .interactive()
            .properties(height=500)
        )
        st.altair_chart(chart)

    # ── Descriptor distributions ──
    with st.container(border=True):
        st.subheader("Descriptor distributions")

        selected_descs = st.multiselect(
            "Select descriptors",
            DESC_COLS,
            default=["MW", "LogP", "TPSA", "HBA", "HBD", "RotBonds"],
        )

        if selected_descs:
            cols = st.columns(min(len(selected_descs), 3))
            for i, desc in enumerate(selected_descs):
                with cols[i % len(cols)]:
                    vals = df[desc].dropna()
                    hist_df = pd.DataFrame({desc: vals})
                    chart = (
                        alt.Chart(hist_df)
                        .mark_bar(opacity=0.85)
                        .encode(
                            x=alt.X(f"{desc}:Q", bin=alt.Bin(maxbins=40), title=desc),
                            y=alt.Y("count()", title="Count"),
                        )
                        .properties(height=200, title=f"{desc} (mean={vals.mean():.1f})")
                    )
                    st.altair_chart(chart)

    # ── Descriptor statistics table ──
    with st.container(border=True):
        st.subheader("Descriptor statistics")
        st.dataframe(
            df[DESC_COLS].describe().round(2),
            column_config={col: st.column_config.NumberColumn(format="%.2f") for col in DESC_COLS},
        )

# ═══════════════════════════════════════════════════════════════════════════
# TAB 2: PCA & t-SNE
# ═══════════════════════════════════════════════════════════════════════════
with tab_pca:

    FINGERPRINTS_AVAILABLE = [
        "maccs", "pubchem", "substructure", "klekota_roth",
        "cdk", "cdk_ext", "cdk_graphonly", "estate",
        "atompairs2d", "substructure_count", "klekota_roth_count", "atompairs2d_count",
    ]

    @st.cache_data
    def load_split_data(fp_name, split_name, subset):
        """Load a split file."""
        path = os.path.join(SPLITS_DIR, f"aromatase_{fp_name}_fp_{split_name}_{subset}.csv")
        if not os.path.exists(path):
            return None
        return pd.read_csv(path)

    @st.cache_data
    def get_target_map():
        """Load pchembl target values."""
        curated = pd.read_csv(CURATED_PATH, usecols=["molecule_chembl_id", "pchembl_value"])
        curated = curated.dropna(subset=["pchembl_value"]).drop_duplicates("molecule_chembl_id")
        return curated.set_index("molecule_chembl_id")["pchembl_value"].to_dict()

    @st.cache_data
    def compute_pca(fp_name):
        """Compute PCA on all molecules for a given fingerprint."""
        from sklearn.decomposition import PCA

        train_rand = load_split_data(fp_name, "random", "train")
        test_rand = load_split_data(fp_name, "random", "test")
        if train_rand is None or test_rand is None:
            return None

        all_mols = pd.concat([train_rand, test_rand], ignore_index=True).drop_duplicates("molecule_chembl_id")
        fp_cols = [c for c in all_mols.columns if c != "molecule_chembl_id"]

        X = all_mols[fp_cols].values.astype(np.float32)
        pca = PCA(n_components=2, random_state=42)
        coords = pca.fit_transform(X)

        all_mols["PC1"] = coords[:, 0]
        all_mols["PC2"] = coords[:, 1]
        var_explained = pca.explained_variance_ratio_

        return all_mols[["molecule_chembl_id", "PC1", "PC2"]], var_explained

    @st.cache_data
    def compute_tsne(fp_name):
        """Compute t-SNE on all molecules for a given fingerprint."""
        from sklearn.manifold import TSNE

        train_rand = load_split_data(fp_name, "random", "train")
        test_rand = load_split_data(fp_name, "random", "test")
        if train_rand is None or test_rand is None:
            return None

        all_mols = pd.concat([train_rand, test_rand], ignore_index=True).drop_duplicates("molecule_chembl_id")
        fp_cols = [c for c in all_mols.columns if c != "molecule_chembl_id"]

        X = all_mols[fp_cols].values.astype(np.float32)
        tsne = TSNE(n_components=2, random_state=42, perplexity=30, n_iter=1000)
        coords = tsne.fit_transform(X)

        all_mols["tSNE1"] = coords[:, 0]
        all_mols["tSNE2"] = coords[:, 1]

        return all_mols[["molecule_chembl_id", "tSNE1", "tSNE2"]]

    # ── Controls ──
    col_fp, col_method = st.columns(2)
    with col_fp:
        fp_choice = st.selectbox("Fingerprint", FINGERPRINTS_AVAILABLE, index=0)
    with col_method:
        method = st.radio("Method", ["PCA", "t-SNE"], horizontal=True)

    target_map = get_target_map()

    # Load split molecule IDs
    train_rand_df = load_split_data(fp_choice, "random", "train")
    test_rand_df = load_split_data(fp_choice, "random", "test")
    train_ks_df = load_split_data(fp_choice, "kennard_stone", "train")
    test_ks_df = load_split_data(fp_choice, "kennard_stone", "test")

    if train_rand_df is None:
        st.warning(f"Split files not found for {fp_choice}")
        st.stop()

    train_rand_ids = set(train_rand_df["molecule_chembl_id"])
    test_rand_ids = set(test_rand_df["molecule_chembl_id"])
    train_ks_ids = set(train_ks_df["molecule_chembl_id"]) if train_ks_df is not None else set()
    test_ks_ids = set(test_ks_df["molecule_chembl_id"]) if test_ks_df is not None else set()

    # Compute embeddings
    if method == "PCA":
        result = compute_pca(fp_choice)
        if result is None:
            st.error("Could not compute PCA")
            st.stop()
        coords_df, var_explained = result
        x_col, y_col = "PC1", "PC2"
        x_title = f"PC1 ({var_explained[0]:.1%})"
        y_title = f"PC2 ({var_explained[1]:.1%})"
    else:
        with st.spinner("Computing t-SNE (may take ~30s)..."):
            coords_df = compute_tsne(fp_choice)
        if coords_df is None:
            st.error("Could not compute t-SNE")
            st.stop()
        x_col, y_col = "tSNE1", "tSNE2"
        x_title = "t-SNE 1"
        y_title = "t-SNE 2"

    # Assign split labels
    def assign_split(mol_id, train_ids, test_ids):
        if mol_id in train_ids:
            return "Train"
        elif mol_id in test_ids:
            return "Test"
        return "Unknown"

    # Build plot DataFrames
    plot_rand = coords_df.copy()
    plot_rand["split"] = plot_rand["molecule_chembl_id"].apply(
        lambda x: assign_split(x, train_rand_ids, test_rand_ids)
    )
    plot_rand["pchembl_value"] = plot_rand["molecule_chembl_id"].map(target_map)

    plot_ks = coords_df.copy()
    plot_ks["split"] = plot_ks["molecule_chembl_id"].apply(
        lambda x: assign_split(x, train_ks_ids, test_ks_ids)
    )
    plot_ks["pchembl_value"] = plot_ks["molecule_chembl_id"].map(target_map)

    # ── Side-by-side scatter plots ──
    st.subheader(f"{method} — Random vs Kennard-Stone ({fp_choice.upper()})")

    col1, col2 = st.columns(2)

    for col, plot_data, title in [
        (col1, plot_rand, "Random Split"),
        (col2, plot_ks, "Kennard-Stone Split"),
    ]:
        with col:
            st.markdown(f"**{title}**")
            chart = (
                alt.Chart(plot_data)
                .mark_circle(size=15, opacity=0.5)
                .encode(
                    x=alt.X(f"{x_col}:Q", title=x_title),
                    y=alt.Y(f"{y_col}:Q", title=y_title),
                    color=alt.Color("split:N",
                                   scale=alt.Scale(domain=["Train", "Test"],
                                                   range=["steelblue", "coral"]),
                                   title="Set"),
                    tooltip=["molecule_chembl_id", "split",
                            alt.Tooltip("pchembl_value:Q", format=".2f")],
                )
                .interactive()
                .properties(height=400)
            )
            st.altair_chart(chart, use_container_width=True)

    # ── pChEMBL distribution comparison ──
    st.subheader("pChEMBL Distribution: Train vs Test")

    col1, col2 = st.columns(2)

    for col, plot_data, title in [
        (col1, plot_rand, "Random Split"),
        (col2, plot_ks, "Kennard-Stone Split"),
    ]:
        with col:
            st.markdown(f"**{title}**")
            hist_data = plot_data.dropna(subset=["pchembl_value"])
            chart = (
                alt.Chart(hist_data)
                .mark_bar(opacity=0.6)
                .encode(
                    x=alt.X("pchembl_value:Q", bin=alt.Bin(maxbins=30), title="pChEMBL value"),
                    y=alt.Y("count()", title="Count", stack=None),
                    color=alt.Color("split:N",
                                   scale=alt.Scale(domain=["Train", "Test"],
                                                   range=["steelblue", "coral"]),
                                   title="Set"),
                )
                .properties(height=250)
            )
            st.altair_chart(chart, use_container_width=True)

    # ── Dispersity Metrics ──
    st.subheader("Data Dispersity")
    st.caption(
        "Measures how spread out the train and test points are in the embedding space. "
        "Higher values = more dispersed/spread out."
    )

    def compute_dispersity(coords_2d):
        """Compute dispersity metrics for a set of 2D points.

        Returns: mean_dist_to_centroid, convex_hull_area, std_spread
        """
        from scipy.spatial import ConvexHull

        if len(coords_2d) < 4:
            return 0.0, 0.0, 0.0

        centroid = coords_2d.mean(axis=0)
        dists_to_centroid = np.sqrt(np.sum((coords_2d - centroid) ** 2, axis=1))
        mean_dist = float(dists_to_centroid.mean())

        try:
            hull = ConvexHull(coords_2d)
            hull_area = float(hull.volume)  # In 2D, 'volume' is area
        except Exception:
            hull_area = 0.0

        std_spread = float(np.sqrt(coords_2d.var(axis=0).sum()))

        return mean_dist, hull_area, std_spread

    # Compute for each split/set combination
    dispersity_rows = []
    for plot_data, split_label in [(plot_rand, "Random"), (plot_ks, "Kennard-Stone")]:
        for set_name in ["Train", "Test"]:
            subset = plot_data[plot_data["split"] == set_name]
            if subset.empty:
                continue
            coords = subset[[x_col, y_col]].values
            mean_dist, hull_area, std_spread = compute_dispersity(coords)
            dispersity_rows.append({
                "Split Method": split_label,
                "Set": set_name,
                "N": len(subset),
                "Mean Dist to Centroid": mean_dist,
                "Convex Hull Area": hull_area,
                "Std Spread": std_spread,
            })

    if dispersity_rows:
        disp_df = pd.DataFrame(dispersity_rows)
        st.dataframe(
            disp_df.style.format({
                "Mean Dist to Centroid": "{:.2f}",
                "Convex Hull Area": "{:.1f}",
                "Std Spread": "{:.2f}",
            }),
            use_container_width=True,
            hide_index=True,
        )
