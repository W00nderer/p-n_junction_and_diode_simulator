import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output
# ----- Constants -------------------------------------------
q = 1.602e-19      # elementary charge [C]
k = 1.381e-23      # Boltzmann constant [J/K]
eps0 = 8.854e-14      # permittivity of free space [F/cm]
epsr = 11.7           # relative permittivity of silicon
epss = eps0 * epsr    # permittivity of silicon [F/cm]
ni = 1.5e10         # intrinsic carrier concentration [cm⁻³]
Eg = 1.12           # bandgap energy [eV]
T = 300            # temperature [K]
Vt = k * T / q     # thermal voltage [V]

# -------- Doping concentration -----------------------------
Na = 1e17             # acceptor concentration [cm⁻³]
Nd = 1e16             # donor concentration [cm⁻³]

# ------- Equilibrium parameters -----------------------------
V0 = Vt * np.log(Na * Nd / ni**2)               # built-in potential [V]
W0 = np.sqrt(2 * epss * V0 / q * (1/Na + 1/Nd)) # equilibrium depletion region width [cm]

# ------- App design constants -------------------------------
font_title = "Playfair Display, serif"
font_text = "'Times New Roman', Times, serif"
colors = {"blue": "#4d6bc0",
          "pink": "#ca5dca",
          "orange": "#e07b3a",
          "green": "#2ca02c",
          "gray": "#808080",
          "dark_gray": "#505050"}
fills = {
    "blue":   "rgba(77, 107, 192, 0.2)",
    "pink":   "rgba(202, 93, 202, 0.2)",
    "gray":   "rgba(128, 128, 128, 0.1)"
}
external_stylesheets = ["https://fonts.googleapis.com/css2?family=Playfair+Display:ital,wght@0,400..900;1,400..900&display=swap"]
# Making sure the calculations are correct
assert np.isclose(Na * W0 * Nd/(Na+Nd), Nd * W0 * Na/(Na+Nd)), "Charge neutrality violated"

# -------- Junction Physics -----------------------------------
# Computing the quasi fermi levels (specifically the split in the depletion region):
def _quasi_fermi_levels(x, xpo_um, xno_um, Efn_bulk, Efp_bulk):
    Efn = np.zeros_like(x) # initialize with the same amount of points as x
    Efp = np.zeros_like(x)

    decay_length_um = 0.6 # where the split should end and the levels come back together

    Efn[x >= -xpo_um] = Efn_bulk # fermi levels are flat in the bulk region
    # Math for depletion region's fermi level:
    Efn[x < -xpo_um] = np.interp( # interpolation, adjusts the fermi line between two points
        x[x < -xpo_um], # the target
        [-xpo_um - decay_length_um, -xpo_um], #x coordinates
        [Efp_bulk, Efn_bulk] # y coordinates
    )

    Efp[x <= xno_um] = Efp_bulk
    Efp[x > xno_um] = np.interp( 
        x[x > xno_um],           
        [xno_um, xno_um + decay_length_um], 
        [Efp_bulk, Efn_bulk] 
    )
    
    return Efn, Efp

# the main physics part
def compute_junction(Va: float):
    Va = min(Va, V0 - 0.001) # Clamp to prevent crashes for forward bias

    # Biased depletion widths
    W   = np.sqrt(2 * epss * (V0 - Va) / q * (1/Na + 1/Nd)) # total depletion width
    xno = W * (Na / (Na + Nd)) # extension into n side 
    xpo = W * (Nd / (Na + Nd)) # extention into p side

    # Spatial grid increased 3 times the depletion width for better visualization
    margin = 3 * W0
    x    = np.linspace(-xpo - margin, xno + margin, 2000) # 2000 points in linear space
    x_um = x * 1e4    # Convert to µm for plotting
    xpo_um = xpo * 1e4
    xno_um = xno * 1e4

    # ----- Analytical potential φ(x) -------------------------------
    phi = np.zeros_like(x) # Initialize the same number of values as x

    phi[x < -xpo] = 0 # Equals to 0 on p side bulk
    # P side depletion region
    mask_p = (x >= -xpo) & (x < 0) # Identifying the points in p side depletion region
    phi[mask_p] = (q * Na / (2 * epss)) * (x[mask_p] + xpo)**2 # Poisson's Equation for a step junction

    phi[x >= xno] = V0 - Va # Equals to V0-Va (since it curves down from the bulk)
    # N side depletion region
    mask_n = (x >= 0) & (x < xno)
    phi[mask_n] = (V0 - Va) - (q * Nd / (2 * epss)) * (x[mask_n] - xno)**2

    # ---------- Band edges ----------------------------------------
    Ec = (V0 - Va) - phi # downward energy curve
    Ev = Ec - Eg # same shape as Ec but down

    Efn_bulk = Vt * np.log(Nd/ni) # Fermi level in n side relative to intrinsic level Ei
    Efp_bulk = Efn_bulk - Va # Difference between Efn and Efp is Va (applied voltage)

    # Calculating quasi fermi levels
    Efn, Efp = _quasi_fermi_levels(x_um, xpo_um, xno_um, Efn_bulk, Efp_bulk)

    # Returning all necessary values
    return x_um, xpo_um, xno_um, Ec, Ev, Efn, Efp


# Models relationship between Va and resulting current I
def shockley(Va_array):
    I0 = 1e-12    # reverse saturation current [A/cm²]
    return I0 * (np.exp(Va_array / Vt) - 1) # return calculated I

# ------- Initial computation ----------------------------------
Va_init = 0.0 # Va starts as 0
x_um, xpo_um, xno_um, Ec, Ev, Efn, Efp = compute_junction(Va_init)

# Compute x and y limits at maximum reverse bias for scaling the plot
_, xno_worst, xpo_worst, Ec_worst, Ev_worst, _, _ = compute_junction(-2.0)
y_bottom = float(Ev_worst.min()) - 0.3
y_top    = float(Ec_worst.max()) + 0.4
x_limit = float(max(xno_worst, xpo_worst)) * 2
x_left = float(-x_limit)
x_right = float(x_limit * 1.5)

Va_sweep = np.linspace(-2, V0 * 0.95, 500) # 500 evenly spaced points from -2 till 95% of V0

iv_x = Va_sweep.tolist() # voltage as x-axis for IV curve graph
iv_y = (shockley(Va_sweep) * 1e3).tolist() # current as y-axis for IV curve graph

# ------- Dash App -------------------------------------------
app = Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

# Layout
app.layout = html.Div([
    # Main title
    html.H1("PN Junction - Interactive Bias Control",
            style={"textAlign": "center", "fontFamily": font_title,
                   "marginBottom": "4px"}),
    # First graph
    dcc.Graph(id = "band-diagram",
              style = {"height": "480px"}),
    # Slider div
    html.Div([
        html.Label("Applied bias Vₐ [V]",
                   style={"fontFamily": font_text, "fontWeight": "bold"}),
        dcc.Slider( # Inititalizing the slider
            id = "va-slider",
            min = -2.0,
            max = round(V0*0.95, 3),
            step = 0.01,
            value = 0.0,
            marks = {v: f"{v:.1f} V"
                     for v in np.arange(-2.0, round(V0 * 0.95, 1) + 0.1, 0.2)},
            tooltip = {"placement": "bottom", "always visible": True},        
        ),
    ], style = {"width": "70%", "margin": "0 auto 20px auto"}),
    # Second graph
    dcc.Graph(id = "iv-curve",
              style = {"height": "400px"}),
], style={"maxWidth": "900px", "margin": "0 auto", "padding":"10px"})

# Callback setup
@app.callback(
    Output("band-diagram", "figure"),
    Output("iv-curve", "figure"),
    Input("va-slider", "value"),
)

# --------- Update function -------------------------------------
def update(Va):
    # In case of None, Va = 0
    if Va is None:
        Va = 0.0
        
    x_um, xpo_um, xno_um, Ec, Ev, Efn, Efp = compute_junction(Va)

    # ----- Band diagram --------------------------
    band_fig = make_subplots(rows=1, cols=1) # Graph skeleton initialized

    # Ec line + shadow
    band_fig.add_trace(go.Scatter(
        x = x_um.tolist(), y=Ec.tolist(),
        name = "E꜀", line = dict(color = colors["blue"], width=2)))
    band_fig.add_trace(go.Scatter(
        x=np.concatenate([x_um, x_um[::-1]]).tolist(), # Loop through x
        y=np.concatenate([Ec + 0.15, Ec[::-1]]).tolist(), # Loop through y
        fill = "toself", fillcolor = fills["blue"] , # fill the loop
        line = dict(width=0), showlegend = False, hoverinfo = "skip"))
    
    # Ev line + shadow
    band_fig.add_trace(go.Scatter(
        x = x_um.tolist(), y=Ev.tolist(),
        name = "Eᵥ", line = dict(color = colors["pink"], width=2)))
    band_fig.add_trace(go.Scatter(
        x=np.concatenate([x_um, x_um[::-1]]).tolist(),
        y=np.concatenate([Ev - 0.15, Ev[::-1]]).tolist(),
        fill = "toself", fillcolor = fills["pink"],
        line = dict(width=0), showlegend = False, hoverinfo = "skip"))
    
    # Quasi Fermi levels
    band_fig.add_trace(go.Scatter(
        x = x_um.tolist(), y=Efn.tolist(),
        name = "Eₙ (quasi)", line = dict(color = colors["orange"], width=1.5, dash="dash")))
    band_fig.add_trace(go.Scatter(
        x = x_um.tolist(), y=Efp.tolist(),
        name = "Eₚ (quasi)", line = dict(color = colors["green"], width=1.5, dash="dash")))
    
    # Texts, legend, depletion region
    band_fig.update_layout(
        title = dict(text = f"Band Diagram - Vₐ = {Va:.2f} V",
                     x = 0.5, font = dict(size = 20, family = font_text)),
        xaxis = dict(title = "Position [µm]", range = [x_left, x_right]),          
        yaxis = dict(title = "Energy [eV]", range = [y_bottom, y_top]),      
        legend = dict(orientation = "h", y = -0.15),
        margin = dict(l = 60, r = 20, t = 50, b = 60),
        annotations = [
            dict(
                x = (xno_um - xpo_um) / 2, 
                y = y_top - 0.1,
                showarrow = False,
                text = "Depletion Region",
                font = dict( family = font_text, size = 12, color = colors["dark_gray"])
            )
        ],
        shapes = [
            # Depletion region
            dict(type = "rect",
                 x0 = -xpo_um, x1 = xno_um,
                 y0 = y_bottom, y1 = y_top,
                 fillcolor = fills["gray"],
                 line = dict(color = "gray", width = 1, dash = "dot")),
            # Line at y = 0
            dict(type = "line",
                 x0 = x_left, x1 = x_right,
                 y0 = 0, y1 = 0,
                 line = dict(color = "black", width = 0.5)),
        ],  
    )

    # ------- I-V curve ------------
    I_op = float(shockley(Va) * 1e3)
    
    iv_fig = go.Figure()
    # Initialize the I-V curve line
    iv_fig.add_trace(go.Scatter(
        x=iv_x, y=iv_y,
        name = "I-V", line = dict(color = colors["blue"], width = 2)))
    # Initialize the op point on the curve
    iv_fig.add_trace(go.Scatter(
        x = [Va], y = [I_op],
        mode = "markers",
        marker = dict(color = colors['orange'], size = 10, symbol = "circle"),
        name = "Operating point"))
    # Title, labels, legend
    iv_fig.update_layout(
        title = dict(text = "I-V Curve (Shockley)", x = 0.5, font = dict(size = 20, family = font_text)),
        xaxis = dict(title = "Applied bias Vₐ [V]", range = [-2, 0.8],
                     zeroline = True, zerolinecolor = "black", zerolinewidth = 0.8),
        yaxis = dict(title = "Current density [mA/cm²]",
                     range = [-0.1, float(shockley(V0 * 0.95) * 1e3) * 1.1],
                     zeroline = True, zerolinecolor = "black", zerolinewidth = 0.8),
        legend = dict(orientation = "h", y = -0.18),
        margin = dict(l=60, r=20, t=50, b=60),          
    )
    
    return band_fig, iv_fig

# Starting the app
if __name__ == "__main__":
    app.run(debug=True)