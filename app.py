import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from engineering_math import (
    vogel_ipr, fetkovich_ipr, calculate_vlp_curve, find_intersection
)

st.set_page_config(page_title="OptiFlow Nodal Analysis", layout="wide")

st.title("OptiFlow : Multiphase Nodal Analysis Software")
st.markdown("Optimize tubing size and analyze well performance using rigorous engineering models.")

# --- Sidebar Inputs ---
st.sidebar.header("Reservoir Parameters")
ipr_method = st.sidebar.selectbox("IPR Model", ["Vogel (Two-Phase)", "Fetkovich"])
pr = st.sidebar.number_input("Reservoir Pressure (Pr) [psia]", value=3000.0, step=100.0)
q_test = st.sidebar.number_input("Test Flow Rate (Q_test) [STB/d]", value=500.0, step=50.0)
pwf_test = st.sidebar.number_input("Test Bottom-hole Pressure (Pwf_test) [psia]", value=2000.0, step=100.0)

if ipr_method == "Vogel (Two-Phase)":
    pb = st.sidebar.number_input("Bubble Point Pressure (Pb) [psia]", value=2500.0, step=100.0)
else:
    n = st.sidebar.number_input("Fetkovich exponent (n)", value=1.0, min_value=0.5, max_value=1.0, step=0.05)

st.sidebar.header("Wellbore/Fluid Parameters")
pvt_method = st.sidebar.selectbox("PVT Correlation", ["Vasquez-Beggs", "Standing"])
depth = st.sidebar.number_input("Well Depth [ft]", value=8000.0, step=500.0)
p_wh = st.sidebar.number_input("Wellhead Pressure (Pwh) [psia]", value=300.0, step=50.0)
t_wh = st.sidebar.number_input("Wellhead Temperature [F]", value=100.0, step=10.0)
t_bh = st.sidebar.number_input("Bottomhole Temperature [F]", value=180.0, step=10.0)

water_cut = st.sidebar.slider("Water Cut [fraction]", 0.0, 1.0, 0.2, 0.05)
gor = st.sidebar.number_input("Gas-Oil Ratio (GOR) [scf/STB]", value=800.0, step=100.0)
api = st.sidebar.number_input("Oil API Gravity", value=35.0, step=1.0)
sg_gas = st.sidebar.number_input("Gas Specific Gravity (air=1)", value=0.7, step=0.05)
roughness = st.sidebar.number_input("Tubing Roughness [in]", value=0.0006, format="%f")

st.sidebar.header("Tubing Optimization")
# Common tubing IDs in inches
available_tubings = [1.995, 2.441, 2.992, 3.958, 4.892]
selected_tubings = st.sidebar.multiselect(
    "Select Tubing Internal Diameters (in)", 
    available_tubings, 
    default=[2.441, 2.992, 3.958]
)

if not selected_tubings:
    st.warning("Please select at least one tubing size to generate VLP curves.")
    st.stop()

# --- Calculations ---
with st.spinner("Calculating IPR and VLP..."):
    # 1. Calculate IPR
    if ipr_method == "Vogel (Two-Phase)":
        q_ipr, p_ipr, qmax = vogel_ipr(pr, pb, q_test, pwf_test)
    else:
        q_ipr, p_ipr, qmax = fetkovich_ipr(pr, q_test, pwf_test, n)
        
    # Set up VLP flow rates based on IPR max rate to cover the required range
    q_vlp_range = np.linspace(0, max(qmax*1.2, 1000), 20)
    
    # 2. Calculate VLP for each tubing size and find intersections
    vlp_results = {}
    metrics_data = []
    
    for d_in in selected_tubings:
        p_vlp = calculate_vlp_curve(p_wh, t_wh, t_bh, depth, q_vlp_range, water_cut, gor, api, sg_gas, d_in, roughness, pvt_method=pvt_method)
        vlp_results[d_in] = p_vlp
        
        q_op, p_op = find_intersection(q_ipr, p_ipr, q_vlp_range, p_vlp)
        
        if q_op is not None and p_op is not None:
            pi = q_op / (pr - p_op) if pr != p_op else 0
            metrics_data.append({
                "Tubing ID (in)": d_in,
                "AOFP [STB/d]": f"{qmax:.1f}",
                "Operating Rate [STB/d]": f"{q_op:.1f}",
                "Operating BHP [psia]": f"{p_op:.1f}",
                "Calculated PI [STB/d/psi]": f"{pi:.3f}",
                "Status": "Flowing"
            })
        else:
            metrics_data.append({
                "Tubing ID (in)": d_in,
                "AOFP [STB/d]": f"{qmax:.1f}",
                "Operating Rate [STB/d]": "0.0",
                "Operating BHP [psia]": "N/A",
                "Calculated PI [STB/d/psi]": "N/A",
                "Status": "Dead Well"
            })
            
# --- Visualization ---

# Check if any well is dead and display warning
dead_wells = [m["Tubing ID (in)"] for m in metrics_data if m["Status"] == "Dead Well"]
if dead_wells:
    st.error(f"⚠️ **Warning**: The well is dead (no intersection) for tubing size(s): {', '.join(map(str, dead_wells))} in. The VLP curve lies entirely above the IPR curve, indicating insufficient reservoir pressure or liquid loading.")

fig = go.Figure()

# Plot IPR
fig.add_trace(go.Scatter(
    x=q_ipr, y=p_ipr, 
    mode='lines', 
    name=f'IPR ({ipr_method})',
    line=dict(color='white', width=3)
))

# Plot VLPs
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
for i, d_in in enumerate(selected_tubings):
    p_vlp = vlp_results[d_in]
    fig.add_trace(go.Scatter(
        x=q_vlp_range, y=p_vlp, 
        mode='lines', 
        name=f'VLP ({d_in}")',
        line=dict(color=colors[i % len(colors)], width=2)
    ))
    
    # Plot Operating Point if flowing
    q_op_str = next((m["Operating Rate [STB/d]"] for m in metrics_data if m["Tubing ID (in)"] == d_in), "0.0")
    if q_op_str != "0.0":
        q_op = float(q_op_str)
        p_op = float(next((m["Operating BHP [psia]"] for m in metrics_data if m["Tubing ID (in)"] == d_in)))
        fig.add_trace(go.Scatter(
            x=[q_op], y=[p_op],
            mode='markers',
            name=f'Op. Pt ({d_in}")',
            marker=dict(size=10, color=colors[i % len(colors)], line=dict(width=2, color='white')),
            showlegend=False,
            hovertemplate=f"Tubing: {d_in} in<br>Rate: {q_op:.1f} STB/d<br>BHP: {p_op:.1f} psia"
        ))

fig.update_layout(
    title="Nodal Analysis (IPR vs VLP)",
    xaxis_title="Liquid Rate [STB/d]",
    yaxis_title="Bottom-hole Pressure [psia]",
    hovermode="x unified",
    template="plotly_white",
    height=600,
    xaxis=dict(range=[0, max(qmax*1.1, max(q_vlp_range))])
)

st.plotly_chart(fig, use_container_width=True)

# --- Output Metrics Table ---
st.subheader("Results Summary")
df_metrics = pd.DataFrame(metrics_data)
st.dataframe(df_metrics, use_container_width=True, hide_index=True)
