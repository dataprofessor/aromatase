"""Temporal & assay analysis page."""

import streamlit as st
import pandas as pd
import altair as alt
from data_utils import load_curated_data, compute_descriptors

df = compute_descriptors(load_curated_data())

st.header("Temporal & assay analysis")

# ── Publication timeline ─────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Records by publication year")
    year_df = (
        df.dropna(subset=["document_year"])
        .groupby("document_year")
        .size()
        .reset_index(name="Records")
    )
    year_df.columns = ["Year", "Records"]
    year_df["Year"] = year_df["Year"].astype(int)

    chart = (
        alt.Chart(year_df)
        .mark_bar(color="#3F51B5", opacity=0.85)
        .encode(
            x=alt.X("Year:O", title="Publication year"),
            y=alt.Y("Records:Q", title="Records"),
            tooltip=["Year", "Records"],
        )
        .properties(height=350)
    )
    st.altair_chart(chart)

    peak_year = year_df.loc[year_df["Records"].idxmax()]
    st.caption(
        f"Peak year: **{int(peak_year['Year'])}** with {int(peak_year['Records'])} records. "
        f"Span: {int(year_df['Year'].min())} - {int(year_df['Year'].max())}."
    )

# ── Records per decade ───────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Records per decade")
    decade_df = df.dropna(subset=["document_year"]).copy()
    decade_df["Decade"] = (decade_df["document_year"] // 10 * 10).astype(int).astype(str) + "s"
    decade_counts = decade_df["Decade"].value_counts().sort_index().reset_index()
    decade_counts.columns = ["Decade", "Records"]

    chart = (
        alt.Chart(decade_counts)
        .mark_bar(cornerRadiusTopRight=4, cornerRadiusTopLeft=4)
        .encode(
            x=alt.X("Decade:N", title="Decade"),
            y=alt.Y("Records:Q", title="Records"),
            color=alt.Color("Decade:N", legend=None),
            tooltip=["Decade", "Records"],
        )
        .properties(height=300)
    )
    st.altair_chart(chart)

# ── pChEMBL over time ────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Median pChEMBL by year")
    yearly_pchembl = (
        df.dropna(subset=["document_year", "pchembl_value"])
        .groupby("document_year")["pchembl_value"]
        .median()
        .reset_index()
    )
    yearly_pchembl.columns = ["Year", "Median pChEMBL"]
    yearly_pchembl["Year"] = yearly_pchembl["Year"].astype(int)

    chart = (
        alt.Chart(yearly_pchembl)
        .mark_line(point=True, color="#FF9800")
        .encode(
            x=alt.X("Year:O", title="Publication year"),
            y=alt.Y("Median pChEMBL:Q", title="Median pChEMBL", scale=alt.Scale(zero=False)),
            tooltip=["Year", "Median pChEMBL"],
        )
        .properties(height=300)
    )
    st.altair_chart(chart)

# ── Top assays ───────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Top assays by record count")

    n_assays = st.slider("Number of assays to show", 5, 20, 10)
    top_assays = df["assay_chembl_id"].value_counts().head(n_assays)
    assay_data = []
    for aid, cnt in top_assays.items():
        desc = df.loc[df["assay_chembl_id"] == aid, "assay_description"].iloc[0]
        desc_short = desc[:80] + "..." if len(str(desc)) > 80 else str(desc)
        assay_data.append({"Assay": aid, "Records": cnt, "Description": desc_short})
    assay_df = pd.DataFrame(assay_data)

    st.dataframe(assay_df, hide_index=True, column_config={
        "Assay": st.column_config.TextColumn(width="small"),
        "Records": st.column_config.NumberColumn(width="small"),
        "Description": st.column_config.TextColumn(width="large"),
    })

# ── Top journals ─────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Top journals")
    journal_df = (
        df["document_journal"]
        .dropna()
        .value_counts()
        .head(10)
        .reset_index()
    )
    journal_df.columns = ["Journal", "Records"]

    chart = (
        alt.Chart(journal_df)
        .mark_bar(cornerRadiusTopRight=4, cornerRadiusBottomRight=4)
        .encode(
            x=alt.X("Records:Q", title="Records"),
            y=alt.Y("Journal:N", sort="-x", title=""),
            color=alt.value("#3F51B5"),
            tooltip=["Journal", "Records"],
        )
        .properties(height=300)
    )
    st.altair_chart(chart)

# ── Assay type breakdown ─────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Assay type and relation breakdown")
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Assay type**")
        atype = df["assay_type"].value_counts().reset_index()
        atype.columns = ["Type", "Count"]
        atype["Label"] = atype["Type"].map({"B": "Binding", "A": "ADMET/Functional"})
        st.dataframe(atype[["Label", "Count"]], hide_index=True)

    with col2:
        st.markdown("**Standard relation**")
        rel = df["standard_relation"].value_counts().reset_index()
        rel.columns = ["Relation", "Count"]
        st.dataframe(rel, hide_index=True)
