# OptiFlow : Multiphase Nodal Analysis Software

OptiFlow is an interactive web application designed for Petroleum Engineers to perform Nodal Analysis in vertical petroleum wells. The software calculates and visualizes Inflow Performance Relationship (IPR) and Tubing Performance Relationship (TPR/VLP) curves, identifies operating points, and allows users to optimize tubing sizes for multiphase flow.

## Benchmarking & Validation
- **PROSPER Comparison**: The custom software outputs (Vogel IPR and Beggs & Brill VLP) were strictly benchmarked against industry-standard PROSPER software. This validates the core correlation behaviors and operating point intersections while highlighting the analytical differences in commercial flow-regime transition smoothing and static fluid density PVT evaluations.
- **Rigorous Physics-Based IPR vs. Vogel**: Included in the repository is an Excel analysis comparing the empirical Vogel IPR method against a rigorous mathematical model. Rather than relying on standard empirical assumptions, this physics-based approach calculates the flow rate ($q$) for each flowing bottom-hole pressure ($P_{wf}$) by dynamically integrating the $k_{ro}/(\mu_o B_o)$ function with respect to pressure, utilizing granular PVT and relative permeability data across all pressure steps.

## Features
- **IPR Models**: Vogel's method for two-phase flow and the Fetkovich method.
- **Multiphase VLP/TPR**: Uses the rigorous Beggs & Brill multiphase flow correlation. Accounts for flow regime changes, friction, hydrostatic elevation, and density changes.
- **PVT Properties**: Dynamically calculates fluid properties (Rs, Bo, Bg, viscosities) up the wellbore using standard Black Oil empirical correlations (Standing and Vasquez-Beggs).
- **Interactive UI**: Powered by Streamlit with interactive Plotly visualizations. Allows tubing size optimization via a multi-select feature.
- **Edge Case Handling**: Identifies dead wells where the VLP curve lies entirely above the IPR curve, displaying a warning while still plotting the curves to allow visual analysis of the pressure gap.

## Installation and Execution

1. **Clone or Download** the repository to your local machine.
2. **Set up a Virtual Environment**:
   ```bash
   python -m venv venv
   # On Windows
   .\venv\Scripts\activate
   # On Mac/Linux
   source venv/bin/activate
   ```
3. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
4. **Run the Application**:
   ```bash
   streamlit run app.py
   ```
   This will start a local web server and open the app in your default web browser (typically at `http://localhost:8501`).

## Multiphase Assumptions and Limitations

- **Vertical Flow**: The current Beggs & Brill implementation assumes a vertical wellbore (90 degrees). Highly deviated or horizontal sections are not currently supported by this VLP calculation script.
- **Temperature Gradient**: A linear temperature gradient is assumed between the bottomhole temperature ($T_{bh}$) and the wellhead temperature ($T_{wh}$). Heat transfer dynamics to the surrounding formation are not rigorously modeled.
- **PVT Correlations**: Black oil correlations (Standing, Vasquez-Beggs) have limited accuracy for volatile oils, gas condensates, or highly heavy oils (API < 15).
- **Steady-State**: The model assumes steady-state flow and does not account for transient behaviors (like wellbore storage or temporary liquid loading un-loading phenomena).
