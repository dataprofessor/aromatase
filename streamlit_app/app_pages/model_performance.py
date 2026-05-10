"""Model Performance page — compare random vs Kennard-Stone split results."""

import os
import streamlit as st
import pandas as pd
import altair as alt

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")
MODELS_DIR = os.path.join(DATA_DIR, "models")

FINGERPRINTS = [
    "maccs", "pubchem", "substructure", "klekota_roth",
    "cdk", "cdk_ext", "cdk_graphonly", "estate",
    "atompairs2d", "substructure_count", "klekota_roth_count", "atompairs2d_count",
]

SPLITS = ["random", "kennard_stone"]

METRIC_LABELS = {"r2": "R²", "rmse": "RMSE", "mae": "MAE"}

EXCLUDED_MODELS = ["Extra Trees", "Gaussian Process"]


SPLITS_DIR = os.path.join(DATA_DIR, "splits")


@st.cache_data
def _get_descriptor_counts():
    """Count descriptors per fingerprint from the split train files."""
    counts = {}
    for fp in FINGERPRINTS:
        path = os.path.join(SPLITS_DIR, f"aromatase_{fp}_fp_random_train.csv")
        if os.path.exists(path):
            header = pd.read_csv(path, nrows=0)
            # Subtract 1 for the molecule_chembl_id column
            counts[fp] = len(header.columns) - 1
    return counts


@st.cache_data
def load_all_results():
    """Load all test result CSVs into a single DataFrame."""
    frames = []
    for fp in FINGERPRINTS:
        for split in SPLITS:
            path = os.path.join(MODELS_DIR, f"{fp}_{split}_split_test.csv")
            if os.path.exists(path):
                df = pd.read_csv(path)
                df["fingerprint"] = fp
                df["split"] = split
                frames.append(df)
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    # Remove models that were not part of the final run
    combined = combined[~combined["model"].isin(EXCLUDED_MODELS)]
    # Add descriptor count per fingerprint
    desc_counts = _get_descriptor_counts()
    combined["n_descriptors"] = combined["fingerprint"].map(desc_counts)
    return combined


results = load_all_results()

if results.empty:
    st.warning("No model results found. Run the model building pipeline first.")
    st.stop()

# ── Header ──
st.header("Model Performance")
st.caption(
    f"{len(results)} results: {results['fingerprint'].nunique()} fingerprints × "
    f"{results['split'].nunique()} splits × {results['model'].nunique()} models"
)

# ── Controls ──
selected_fps = FINGERPRINTS
filtered = results[results["fingerprint"].isin(selected_fps)]

# ── Side-by-side heatmaps: Random vs Kennard-Stone ──
# Initialize metric in session state for use in the subheader before the widget renders
if "metric" not in st.session_state:
    st.session_state["metric"] = "r2"

with st.container(border=True):
    st.subheader(f"Test {METRIC_LABELS[st.session_state['metric']]} — Random vs Kennard-Stone")
    st.caption("Heatmap comparing model performance across fingerprints for each data split. Rows are models, columns are fingerprints.")
    col_metric, col_sort_metric, _ = st.columns(3)
    with col_metric:
        metric = st.radio("Metric", options=["r2", "rmse", "mae"], horizontal=True, key="metric")
    with col_sort_metric:
        sort_order = st.radio("Sort axes by", options=["Performance", "Alphabetical"], horizontal=True, key="sort_metric")

    col1, col2 = st.columns(2)

    for col, split_name, title in [
        (col1, "random", "Random Split"),
        (col2, "kennard_stone", "Kennard-Stone Split"),
    ]:
        with col:
            st.markdown(f"**{title}**")
            split_df = filtered[filtered["split"] == split_name]
            if split_df.empty:
                st.info(f"No results for {title}")
                continue

            pivot = split_df.pivot_table(
                index="model", columns="fingerprint", values=metric, aggfunc="first"
            )
            if sort_order == "Alphabetical":
                cols_order = sorted([f for f in selected_fps if f in pivot.columns])
                pivot = pivot[cols_order]
                pivot = pivot.loc[sorted(pivot.index)]
            else:
                # Sort by this split's own performance data
                cols_available = [f for f in selected_fps if f in pivot.columns]
                col_means = pivot[cols_available].mean(axis=0)
                cols_order = list(col_means.sort_values(ascending=(metric != "r2")).index)
                pivot = pivot[cols_order]
                # Sort rows by mean performance
                pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=(metric != "r2")).index]

            # Melt for Altair
            melted = pivot.reset_index().melt(id_vars="model", var_name="fingerprint", value_name=metric)

            # Color scale
            if metric == "r2":
                color_scale = alt.Scale(scheme="redyellowgreen", domain=[0, 0.7])
            else:
                color_scale = alt.Scale(scheme="redyellowgreen", reverse=True)

            chart = (
                alt.Chart(melted)
                .mark_rect()
                .encode(
                    x=alt.X("fingerprint:N", sort=cols_order, title=None,
                            axis=alt.Axis(labelLimit=200, tickColor="white", domainColor="white",
                                          labelColor="white")),
                    y=alt.Y("model:N", sort=list(pivot.index), title=None,
                            axis=alt.Axis(labelLimit=200, tickColor="white", domainColor="white",
                                          labelColor="white")),
                    color=alt.Color(f"{metric}:Q", scale=color_scale, title=METRIC_LABELS[metric]),
                    tooltip=["model", "fingerprint", alt.Tooltip(f"{metric}:Q", format=".4f")],
                )
                .properties(height=max(400, len(pivot.index) * 35))
            )
            st.altair_chart(chart, use_container_width=True)

# ── Side-by-side heatmaps: Training Time ──
with st.container(border=True):
    st.subheader("Training Time (seconds) — Random vs Kennard-Stone")
    st.caption("Heatmap showing training time for each model-fingerprint combination. Green = fast, red = slow.")
    sort_order_time = st.radio("Sort axes by", options=["Performance", "Alphabetical"], horizontal=True, key="sort_time")

    col_t1, col_t2 = st.columns(2)

    for col, split_name, title in [
        (col_t1, "random", "Random Split"),
        (col_t2, "kennard_stone", "Kennard-Stone Split"),
    ]:
        with col:
            st.markdown(f"**{title}**")
            split_df = filtered[filtered["split"] == split_name]
            if split_df.empty:
                st.info(f"No results for {title}")
                continue

            pivot = split_df.pivot_table(
                index="model", columns="fingerprint", values="time_s", aggfunc="first"
            )
            if sort_order_time == "Alphabetical":
                cols_order = sorted([f for f in selected_fps if f in pivot.columns])
                pivot = pivot[cols_order]
                pivot = pivot.loc[sorted(pivot.index)]
            else:
                # Sort by this split's own time data (fastest first)
                cols_available = [f for f in selected_fps if f in pivot.columns]
                col_means = pivot[cols_available].mean(axis=0)
                cols_order = list(col_means.sort_values(ascending=True).index)
                pivot = pivot[cols_order]
                # Sort rows by mean time (fastest on top)
                pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=True).index]

            melted = pivot.reset_index().melt(id_vars="model", var_name="fingerprint", value_name="time_s")

            time_chart = (
                alt.Chart(melted)
                .mark_rect()
                .encode(
                    x=alt.X("fingerprint:N", sort=cols_order, title=None,
                            axis=alt.Axis(labelLimit=200, tickColor="white", domainColor="white",
                                          labelColor="white")),
                    y=alt.Y("model:N", sort=list(pivot.index), title=None,
                            axis=alt.Axis(labelLimit=200, tickColor="white", domainColor="white",
                                          labelColor="white")),
                    color=alt.Color("time_s:Q",
                                    scale=alt.Scale(scheme="redyellowgreen", reverse=True),
                                    title="Time (s)"),
                    tooltip=["model", "fingerprint", alt.Tooltip("time_s:Q", format=".1f")],
                )
                .properties(height=max(400, len(pivot.index) * 35))
            )
            st.altair_chart(time_chart, use_container_width=True)

# ── Best model per fingerprint comparison ──
with st.container(border=True):
    st.subheader("Best Model per Fingerprint")
    st.caption("Bar chart showing the best-performing model's score for each fingerprint, comparing both split methods side by side.")

    sort_split = st.radio("Sort by split", options=["random", "kennard_stone"], horizontal=True)

    best_rows = []
    for fp in selected_fps:
        for split in SPLITS:
            subset = filtered[(filtered["fingerprint"] == fp) & (filtered["split"] == split)]
            if subset.empty:
                continue
            if metric == "r2":
                best = subset.loc[subset["r2"].idxmax()]
            else:
                best = subset.loc[subset[metric].idxmin()]
            best_rows.append({
                "fingerprint": fp,
                "split": split,
                "best_model": best["model"],
                metric: best[metric],
            })

    if best_rows:
        best_df = pd.DataFrame(best_rows)

        # Bar chart comparing best metric per FP for both splits
        chart_data = best_df.copy()
        chart_data["split_label"] = chart_data["split"].map(
            {"random": "Random", "kennard_stone": "Kennard-Stone"}
        )

        if sort_order == "Alphabetical":
            bar_fp_order = sorted(selected_fps)
        else:
            # Sort by best metric value in the selected split
            sort_subset = best_df[best_df["split"] == sort_split]
            if not sort_subset.empty:
                bar_fp_order = list(
                    sort_subset.sort_values(metric, ascending=(metric != "r2"))["fingerprint"]
                )
            else:
                bar_fp_order = sorted(selected_fps)

        bar = (
            alt.Chart(chart_data)
            .mark_bar()
            .encode(
                x=alt.X("fingerprint:N", sort=bar_fp_order, title="Fingerprint"),
                y=alt.Y(f"{metric}:Q", title=METRIC_LABELS[metric]),
                color=alt.Color("split_label:N", title="Split",
                               scale=alt.Scale(range=["steelblue", "coral"])),
                xOffset="split_label:N",
                tooltip=["fingerprint", "split_label", "best_model",
                        alt.Tooltip(f"{metric}:Q", format=".4f")],
            )
            .properties(height=350)
        )
        st.altair_chart(bar, use_container_width=True)

# ── Detailed table ──
with st.container(border=True):
    st.subheader("Detailed Results")
    st.caption("Full table of all model results with metrics and training time. Filter by split or model to narrow down.")

    col_split, col_model = st.columns(2)
    with col_split:
        split_filter = st.selectbox("Split", ["All"] + SPLITS)
    with col_model:
        models_available = sorted(filtered["model"].unique())
        model_filter = st.selectbox("Model", ["All"] + models_available)

    table_df = filtered.copy()
    if split_filter != "All":
        table_df = table_df[table_df["split"] == split_filter]
    if model_filter != "All":
        table_df = table_df[table_df["model"] == model_filter]

    # Sort by chosen metric
    ascending = metric != "r2"
    table_df = table_df.sort_values(metric, ascending=ascending)

    st.dataframe(
        table_df[["fingerprint", "n_descriptors", "split", "model", "r2", "rmse", "mae", "time_s"]].round(4),
        use_container_width=True,
        hide_index=True,
        column_config={"n_descriptors": st.column_config.NumberColumn("# Descriptors")},
    )
