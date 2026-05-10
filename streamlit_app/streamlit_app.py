"""Aromatase (CYP19A1) Bioactivity EDA — Multi-page Streamlit app."""

import streamlit as st

st.set_page_config(
    page_title="Aromatase EDA",
    page_icon=":material/science:",
    layout="wide",
)

page = st.navigation(
    [
        st.Page("app_pages/overview.py", title="Overview", icon=":material/dashboard:"),
        st.Page("app_pages/potency.py", title="Potency analysis", icon=":material/bar_chart:"),
        st.Page("app_pages/chemical_space.py", title="Chemical space", icon=":material/hexagon:"),
        st.Page("app_pages/molecular_fingerprints.py", title="Molecular Fingerprints", icon=":material/fingerprint:"),
        st.Page("app_pages/temporal.py", title="Temporal & assays", icon=":material/calendar_today:"),
        st.Page("app_pages/model_performance.py", title="Model performance", icon=":material/model_training:"),
    ],
    position="top",
)

page.run()
