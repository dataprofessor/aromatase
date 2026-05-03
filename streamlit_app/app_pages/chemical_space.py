"""Chemical space page — molecular descriptors and drug-likeness."""

import streamlit as st
import pandas as pd
import altair as alt
from data_utils import (
    load_curated_data, compute_descriptors,
    DESC_COLS, POTENCY_ORDER, POTENCY_COLORS,
)

df = compute_descriptors(load_curated_data())

st.header("Chemical space analysis")

# ── Lipinski Rule-of-5 KPIs ──────────────────────────────────────────────
ro5_mask = (df["MW"] <= 500) & (df["LogP"] <= 5) & (df["HBA"] <= 10) & (df["HBD"] <= 5)
n_total = len(df)

with st.container(horizontal=True):
    st.metric("Ro5 compliant", f"{ro5_mask.sum():,} ({ro5_mask.mean()*100:.1f}%)", border=True)
    st.metric("MW > 500", f"{(df['MW'] > 500).sum()} ({(df['MW'] > 500).mean()*100:.1f}%)", border=True)
    st.metric("LogP > 5", f"{(df['LogP'] > 5).sum()} ({(df['LogP'] > 5).mean()*100:.1f}%)", border=True)
    st.metric("HBA > 10", f"{(df['HBA'] > 10).sum()} ({(df['HBA'] > 10).mean()*100:.1f}%)", border=True)
    st.metric("HBD > 5", f"{(df['HBD'] > 5).sum()} ({(df['HBD'] > 5).mean()*100:.1f}%)", border=True)

# ── Chemical space scatter ───────────────────────────────────────────────
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

# ── Descriptor distributions ─────────────────────────────────────────────
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

# ── Descriptor statistics table ──────────────────────────────────────────
with st.container(border=True):
    st.subheader("Descriptor statistics")
    st.dataframe(
        df[DESC_COLS].describe().round(2),
        column_config={col: st.column_config.NumberColumn(format="%.2f") for col in DESC_COLS},
    )
