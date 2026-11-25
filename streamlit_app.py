"""Interactive dashboard:  streamlit run streamlit_app.py"""
import numpy as np, matplotlib.pyplot as plt, streamlit as st
from fin_core import fin, fd_fin, h_forced, annular_eta

st.title("Engine Fin Convection Dashboard")
V = st.sidebar.slider("Air velocity (m/s)", 0.5, 35.0, 10.0)
L = st.sidebar.slider("Fin height (mm)", 10, 80, 50)/1000
t = st.sidebar.slider("Fin thickness (mm)", 1.0, 8.0, 3.0)/1000
k = st.sidebar.selectbox("Material k (W/mK)", [45, 52, 156, 180, 390], index=3)
Tw = st.sidebar.slider("Base temperature (°C)", 100, 350, 200)
h, Re = h_forced(V, Tw=Tw); Q, eta, m = fin(h, k, L, t, Tw)
c = st.columns(4); c[0].metric("h (W/m²K)", f"{h:.1f}"); c[1].metric("Re", f"{Re:,.0f}")
c[2].metric("Q (W)", f"{Q:.1f}"); c[3].metric("Efficiency", f"{eta:.3f}")
Qfd, th = fd_fin(h, k=k, L=L, t=t, Tw=Tw)
fig, ax = plt.subplots(figsize=(7, 4)); ax.plot(np.linspace(0, L*1000, len(th)), th+40, lw=2)
ax.set(xlabel="Fin length (mm)", ylabel="T (°C)", title="Temperature along the fin"); ax.grid(alpha=.3); st.pyplot(fig)
st.caption(f"Finite-difference Q = {Qfd:.1f} W; annular-fin efficiency (r1=40 mm) = {annular_eta(h, k, t, 0.04, L):.3f}")
