import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

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
y_bottom = Ev_worst.min() - 0.3
y_top    = Ec_worst.max() + 0.4
x_limit = max(xno_worst, xpo_worst) * 2
x_left = -x_limit
x_right = x_limit * 1.5

# ------- Figure layout -------------------------------------------
fig = plt.figure(figsize=(13, 6))
fig.suptitle("PN junction — interactive bias control", fontsize=13)
ax_band = fig.add_subplot(1, 2, 1)
ax_iv   = fig.add_subplot(1, 2, 2)
plt.subplots_adjust(
    left=0.08,  
    right=0.96,
    top=0.90,   
    bottom=0.25,
    wspace=0.2,  
    hspace=0.2   
)
# -------- Band diagram ------------------------------------------------
band_shade = 0.15 # For visualizing which direction the bands continue towards

# Plotting conduction and valence band lines, as well as the quasi fermi levels
l_ec, = ax_band.plot(x_um, Ec,  color="#4d6bc0", lw=2,   label="$E_c$")
l_ev, = ax_band.plot(x_um, Ev,  color="#ca5dca", lw=2,   label="$E_v$")
l_efn, = ax_band.plot(x_um, Efn, color="#e07b3a", lw=1.5, label="$E_{fn}$", linestyle="--")
l_efp, = ax_band.plot(x_um, Efp, color="#2ca02c", lw=1.5, label="$E_{fp}$", linestyle="--")

# Shade on the conduction and valence band lines
fill_ec = ax_band.fill_between(x_um, Ec, Ec + band_shade, alpha=0.15, color="#4d6bc0")
fill_ev = ax_band.fill_between(x_um, Ev, Ev - band_shade, alpha=0.15, color="#ca5dca")

# Making the depletion region
dep_span = ax_band.axvspan(-xpo_um, xno_um, color="gray", alpha=0.08, label="Depletion region")
dep_left = ax_band.axvline(-xpo_um, color="gray", lw=0.8, ls="--")
dep_right = ax_band.axvline( xno_um, color="gray", lw=0.8, ls="--")
ax_band.axhline(0, color="black", lw=0.5)

# Text at to top left indicating the current Va set
va_text = ax_band.text(0.05, 0.95, f"$V_a$ = {Va_init:.2f} V",
                       transform=ax_band.transAxes,
                       fontsize=10, va="top", color="gray")

# Labels and legend
ax_band.set_xlabel("Position (µm)")
ax_band.set_ylabel("Energy (eV)")
ax_band.set_title("Band diagram")
ax_band.legend(loc="upper right")
ax_band.set_ylim(y_bottom, y_top)
ax_band.set_xlim(x_left, x_right)

# ------- I-V curve ---------------------------------------------------
Va_sweep = np.linspace(-2, V0 * 0.95, 500) # 500 evenly spaced points from -2 till 95% of V0

# Plot the I-V curve line
ax_iv.plot(Va_sweep, shockley(Va_sweep) * 1e3, color="#4d6bc0", lw=2)
ax_iv.axhline(0, color="black", lw=0.5)
ax_iv.axvline(0, color="black", lw=0.5)

# Plot the dot tracking the line as Va is adjusted
op_dot, = ax_iv.plot(Va_init, shockley(Va_init) * 1e3,
                     'o', color="#e07b3a", markersize=8, zorder=5)

# Labels
ax_iv.set_xlabel("Applied bias $V_a$ (V)")
ax_iv.set_ylabel("Current density (mA/cm²)")
ax_iv.set_title("I-V curve (Shockley)")
ax_iv.set_xlim(-2, 0.8) # from the start or Va_Sweep until the maximum for Si diode
ax_iv.set_ylim(-0.1, shockley(V0 * 0.95) * 1e3 * 1.1) # From -1 until shockely at 95% of Va + 10% headroom

# ------- Slider ------------------------------------------
ax_slider = plt.axes([0.2, 0.08, 0.6, 0.03]) # slider dimensions
# Initializing the slider
slider = Slider(ax_slider, "$V_a$ (V)", -2.0, V0 * 0.95,
                   valinit=Va_init, color="#7f77dd")

# --------- Update function --------------------------------
def update(val):
    Va = val
    x_um, xpo_um, xno_um, Ec, Ev, Efn, Efp = compute_junction(Va)

    # Update band lines
    l_ec.set_xdata(x_um);  l_ec.set_ydata(Ec)
    l_ev.set_xdata(x_um);  l_ev.set_ydata(Ev)
    l_efn.set_xdata(x_um); l_efn.set_ydata(Efn)
    l_efp.set_xdata(x_um); l_efp.set_ydata(Efp)

    # Update band shading
    global fill_ec, fill_ev
    fill_ec.remove()
    fill_ev.remove()
    fill_ec = ax_band.fill_between(x_um, Ec, Ec + band_shade, alpha=0.15, color="#4d6bc0")
    fill_ev = ax_band.fill_between(x_um, Ev, Ev - band_shade, alpha=0.15, color="#ca5dca")

    # Update depletion region shading
    global dep_span, dep_left, dep_right
    dep_span.remove()
    dep_left.remove()
    dep_right.remove()
    dep_span  = ax_band.axvspan(-xpo_um, xno_um, color="gray", alpha=0.08)
    dep_left  = ax_band.axvline(-xpo_um, color="gray", lw=0.8, ls="--")
    dep_right = ax_band.axvline( xno_um, color="gray", lw=0.8, ls="--")

    # Update I-V operating point and Va label
    op_dot.set_data([Va], [shockley(Va) * 1e3])
    va_text.set_text(f"$V_a$ = {Va:.2f} V")

    fig.canvas.draw_idle()

slider.on_changed(update)
plt.show()

