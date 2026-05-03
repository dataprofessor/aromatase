"""Potency analysis page — pChEMBL distributions and measurement types."""

import streamlit as st
import pandas as pd
import altair as alt
from data_utils import load_curated_data, compute_descriptors, DESC_COLS

df = compute_descriptors(load_curated_data())

st.header("pChEMBL value analysis")

# ── Filters ──────────────────────────────────────────────────────────────
with st.container(horizontal=True):
    selected_types = st.multiselect(
        "Measurement type",
        options=df["standard_type"].unique().tolist(),
        default=df["standard_type"].unique().tolist(),
    )
    relation_filter = st.multiselect(
        "Relation",
        options=df["standard_relation"].dropna().unique().tolist(),
        default=["="],
    )

filtered = df[
    df["standard_type"].isin(selected_types)
    & df["standard_relation"].isin(relation_filter)
]
pvals = filtered["pchembl_value"].dropna()

st.caption(f"Showing {len(filtered):,} records ({len(pvals):,} with pChEMBL values)")

# ── Distribution ─────────────────────────────────────────────────────────
col1, col2 = st.columns(2)

with col1:
    with st.container(border=True):
        st.subheader("pChEMBL distribution")
        if len(pvals) > 0:
            hist_data = pd.DataFrame({"pChEMBL": pvals})
            chart = (
                alt.Chart(hist_data)
                .mark_bar(color="#2196F3", opacity=0.85)
                .encode(
                    x=alt.X("pChEMBL:Q", bin=alt.Bin(maxbins=50), title="pChEMBL value"),
                    y=alt.Y("count()", title="Count"),
                )
                .properties(height=350)
            )
            median_rule = (
                alt.Chart(pd.DataFrame({"x": [pvals.median()]}))
                .mark_rule(color="red", strokeDash=[4, 4], size=2)
                .encode(x="x:Q")
            )
            st.altair_chart(chart + median_rule)
            st.caption(f"Red line = median ({pvals.median():.2f})")
        else:
            st.info("No data matches the current filters.")

with col2:
    with st.container(border=True):
        st.subheader("pChEMBL by measurement type")
        ptype = filtered.dropna(subset=["pchembl_value"])
        if len(ptype) > 0:
            chart = (
                alt.Chart(ptype)
                .mark_bar(opacity=0.6)
                .encode(
                    x=alt.X("pchembl_value:Q", bin=alt.Bin(maxbins=40), title="pChEMBL value"),
                    y=alt.Y("count()", title="Count", stack=None),
                    color=alt.Color("standard_type:N", scale=alt.Scale(
                        domain=["IC50", "Ki", "pIC50"],
                        range=["#2196F3", "#FF9800", "#4CAF50"],
                    )),
                )
                .properties(height=350)
            )
            st.altair_chart(chart)

# ── Summary stats ────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Summary statistics by measurement type")
    if len(ptype) > 0:
        stats = ptype.groupby("standard_type")["pchembl_value"].describe().round(3)
        st.dataframe(stats)

# ── Deduplication ────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Deduplication statistics")
    n_meas = filtered["n_measurements"]
    col_a, col_b, col_c = st.columns(3)
    col_a.metric("Singleton (n=1)", f"{(n_meas == 1).sum():,}", border=True)
    col_b.metric("Averaged (n>1)", f"{(n_meas > 1).sum():,}", border=True)
    col_c.metric("Max averaged", f"{n_meas.max()}", border=True)

    averaged = filtered[filtered["n_measurements"] > 1]
    if len(averaged) > 0:
        sd_data = pd.DataFrame({"SD pChEMBL": averaged["sd_pchembl"]})
        chart = (
            alt.Chart(sd_data)
            .mark_bar(color="#FF9800", opacity=0.85)
            .encode(
                x=alt.X("SD pChEMBL:Q", bin=alt.Bin(maxbins=40), title="SD of pChEMBL"),
                y=alt.Y("count()", title="Count"),
            )
            .properties(height=250)
        )
        st.altair_chart(chart)
        st.caption(f"Mean SD = {averaged['sd_pchembl'].mean():.3f} across {len(averaged):,} averaged records")

# ── Correlations ─────────────────────────────────────────────────────────
with st.container(border=True):
    st.subheader("Descriptor correlations with pChEMBL")
    corrs = []
    for col in DESC_COLS:
        valid = filtered[["pchembl_value", col]].dropna()
        if len(valid) > 2:
            r = valid.corr().iloc[0, 1]
            corrs.append({"Descriptor": col, "Pearson r": round(r, 3)})
    if corrs:
        corr_df = pd.DataFrame(corrs).sort_values("Pearson r", key=abs, ascending=False)
        chart = (
            alt.Chart(corr_df)
            .mark_bar()
            .encode(
                x=alt.X("Pearson r:Q", scale=alt.Scale(domain=[-0.3, 0.3])),
                y=alt.Y("Descriptor:N", sort=corr_df["Descriptor"].tolist()),
                color=alt.condition(
                    alt.datum["Pearson r"] > 0,
                    alt.value("#4CAF50"),
                    alt.value("#F44336"),
                ),
                tooltip=["Descriptor", "Pearson r"],
            )
            .properties(height=300)
        )
        st.altair_chart(chart)
