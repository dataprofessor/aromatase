"""Overview page — KPI metrics and dataset summary."""

import streamlit as st
import pandas as pd
import altair as alt
from data_utils import load_curated_data, compute_descriptors, POTENCY_ORDER, POTENCY_COLORS

df_raw = load_curated_data()
df = compute_descriptors(df_raw)

# ── KPI row ──────────────────────────────────────────────────────────────
st.header("Dataset at a glance")

with st.container(horizontal=True):
    st.metric("Total records", f"{len(df):,}", border=True)
    st.metric("Unique molecules", f"{df['molecule_chembl_id'].nunique():,}", border=True)
    st.metric("Unique assays", f"{df['assay_chembl_id'].nunique()}", border=True)
    st.metric("Publications", f"{df['document_chembl_id'].nunique()}", border=True)

pvals = df["pchembl_value"].dropna()
with st.container(horizontal=True):
    st.metric("Median pChEMBL", f"{pvals.median():.2f}", border=True)
    st.metric("pChEMBL range", f"{pvals.min():.1f} - {pvals.max():.1f}", border=True)
    st.metric("Lipinski Ro5 pass", f"{((df['MW'] <= 500) & (df['LogP'] <= 5) & (df['HBA'] <= 10) & (df['HBD'] <= 5)).sum():,} / {len(df):,}", border=True)
    st.metric("Years span", f"{int(df['document_year'].min())} - {int(df['document_year'].max())}", border=True)

# ── Charts row ───────────────────────────────────────────────────────────
col1, col2 = st.columns(2)

with col1:
    with st.container(border=True):
        st.subheader("Measurement type breakdown")
        type_df = df["standard_type"].value_counts().reset_index()
        type_df.columns = ["Type", "Count"]
        chart = (
            alt.Chart(type_df)
            .mark_arc(innerRadius=50)
            .encode(
                theta=alt.Theta("Count:Q"),
                color=alt.Color("Type:N", scale=alt.Scale(
                    domain=["IC50", "Ki", "pIC50"],
                    range=["#2196F3", "#FF9800", "#4CAF50"],
                )),
                tooltip=["Type", "Count"],
            )
            .properties(height=300)
        )
        st.altair_chart(chart)

with col2:
    with st.container(border=True):
        st.subheader("Potency classification")
        pot_df = df[df["potency_class"] != "Unknown"]["potency_class"].value_counts()
        pot_data = pd.DataFrame({
            "Class": POTENCY_ORDER,
            "Count": [pot_df.get(p, 0) for p in POTENCY_ORDER],
        })
        chart = (
            alt.Chart(pot_data)
            .mark_bar(cornerRadiusTopRight=4, cornerRadiusBottomRight=4)
            .encode(
                x=alt.X("Count:Q", title="Count"),
                y=alt.Y("Class:N", sort=POTENCY_ORDER, title=""),
                color=alt.Color("Class:N", scale=alt.Scale(
                    domain=POTENCY_ORDER,
                    range=[POTENCY_COLORS[p] for p in POTENCY_ORDER],
                ), legend=None),
                tooltip=["Class", "Count"],
            )
            .properties(height=300)
        )
        st.altair_chart(chart)

# ── Missing data ─────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Missing data summary")
    missing = df_raw.isnull().sum()
    missing_pct = (missing / len(df_raw) * 100).round(2)
    missing_df = pd.DataFrame({"Column": missing.index, "Missing": missing.values, "Percent": missing_pct.values})
    missing_df = missing_df[missing_df["Missing"] > 0].sort_values("Percent", ascending=False).reset_index(drop=True)
    if len(missing_df) > 0:
        st.dataframe(
            missing_df,
            hide_index=True,
            column_config={
                "Percent": st.column_config.NumberColumn(format="%.2f%%"),
            },
        )
    else:
        st.success("No missing values in any column.")

# ── Data preview ─────────────────────────────────────────────────────────
with st.expander("Raw data preview", icon=":material/table_chart:"):
    st.dataframe(df_raw.head(100), hide_index=True)
