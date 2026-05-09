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
    return pd.concat(frames, ignore_index=True)


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
col_fp, col_metric = st.columns(2)

with col_fp:
    selected_fps = st.multiselect(
        "Fingerprints",
        options=FINGERPRINTS,
        default=FINGERPRINTS,
    )

with col_metric:
    metric = st.radio("Metric", options=["r2", "rmse", "mae"], horizontal=True)

filtered = results[results["fingerprint"].isin(selected_fps)]

# ── Side-by-side heatmaps: Random vs Kennard-Stone ──
st.subheader(f"Test {METRIC_LABELS[metric]} — Random vs Kennard-Stone")

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
        # Reorder columns by selected fingerprints order
        cols_order = [f for f in selected_fps if f in pivot.columns]
        pivot = pivot[cols_order]

        # Sort rows by mean metric
        ascending = metric != "r2"
        pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=ascending).index]

        # Color scale
        if metric == "r2":
            color_scale = alt.Scale(scheme="redyellowgreen", domain=[0, 0.7])
        else:
            color_scale = alt.Scale(scheme="redyellowgreen", reverse=True)

        # Melt for Altair
        melted = pivot.reset_index().melt(id_vars="model", var_name="fingerprint", value_name=metric)

        chart = (
            alt.Chart(melted)
            .mark_rect()
            .encode(
                x=alt.X("fingerprint:N", sort=cols_order, title=None),
                y=alt.Y("model:N", sort=list(pivot.index), title=None),
                color=alt.Color(f"{metric}:Q", scale=color_scale, title=METRIC_LABELS[metric]),
                tooltip=["model", "fingerprint", alt.Tooltip(f"{metric}:Q", format=".4f")],
            )
            .properties(height=400)
        )
        st.altair_chart(chart, use_container_width=True)

# ── Best model per fingerprint comparison ──
st.subheader("Best Model per Fingerprint")

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
    pivot_best = best_df.pivot_table(
        index="fingerprint", columns="split", values=metric, aggfunc="first"
    )
    # Reorder
    cols_present = [s for s in SPLITS if s in pivot_best.columns]
    pivot_best = pivot_best[cols_present]
    pivot_best = pivot_best.loc[[f for f in selected_fps if f in pivot_best.index]]

    # Bar chart comparing best R² per FP for both splits
    chart_data = best_df.copy()
    chart_data["split_label"] = chart_data["split"].map(
        {"random": "Random", "kennard_stone": "Kennard-Stone"}
    )

    bar = (
        alt.Chart(chart_data)
        .mark_bar()
        .encode(
            x=alt.X("fingerprint:N", sort=selected_fps, title="Fingerprint"),
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
st.subheader("Detailed Results")

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
    table_df[["fingerprint", "split", "model", "r2", "rmse", "mae", "time_s"]].round(4),
    use_container_width=True,
    hide_index=True,
)
