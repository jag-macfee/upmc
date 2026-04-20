import time
from solver import solve_flat_plate_unsteady, solve_flat_plate_upm_init
import numpy as np
import matplotlib.pyplot as plt
import sys
from plotter import (
    plot_gamma_distributions_by_timestep,
    plot_gamma_surface_over_time,
    plot_normalised_lift_components_over_time,
    plot_lift_over_time,
    plot_lift_frequency_spectrum,
)
from utils import get_lift_history_components

# Freestream / gust properties
U_inf = 10
v_0 = 1
rho = 1.225

# Airfoil / config properties
c = 1
n_panels = 70
k_max = 10 * n_panels  # Length of the wake in terms of panels
panel_width = c / n_panels
delta_t = panel_width / U_inf


# Have to make sure dirac is shifted for x_n rather than gamma_n
# Note: Due to floating-point precision in calculation of t and delta_t,
# we have to acknowledge their difference to be system eps
def dirac(t):
    if abs(t) <= sys.float_info.epsilon:
        return 1.0
    return 0.0


# Have gusts with magnitude one for now, we will sort out normalising later
# will need to normalise with respect to each xi position
half_chord_per_cycle_omega = 2 * np.pi / (0.5 * c / U_inf)
v_gust = lambda t: dirac(t - 0.5 * delta_t)

start_time = time.time()
gamma, gamma_wake, xi, xi_wake = solve_flat_plate_unsteady(
    c, n_panels, v_gust, k_max, U_inf
)
end_time = time.time()
print(f"Solver execution time: {end_time - start_time:.4f} seconds")

# Preserve raw circulation strengths for force calculations.
gamma_raw = np.copy(gamma)

# normalise
gamma = gamma / panel_width
gamma_wake = gamma_wake / panel_width
xi = xi / c
xi_wake = xi_wake / c

# Lift histories
gamma_init, xi_init, x_init = solve_flat_plate_upm_init(v_0, c, n_panels)
circulatory_history, non_circulatory_history, total_lift_history = (
    get_lift_history_components(k_max, rho, U_inf, delta_t, panel_width, gamma_raw)
)

# Save lift histories to a CSV file
lift_data = np.column_stack(
    (circulatory_history, non_circulatory_history, total_lift_history)
)
np.savetxt(
    "lift_history.csv",
    lift_data,
    delimiter=",",
    header="Circulatory,Non-Circulatory,Total",
    comments="",
)

# Plot
plot_gamma_surface_over_time(gamma, gamma_wake, xi, xi_wake, k_max, show=False)
plot_gamma_distributions_by_timestep(gamma, gamma_wake, xi, xi_wake, k_max, show=False)
plot_normalised_lift_components_over_time(
    circulatory_history,
    non_circulatory_history,
    total_lift_history,
    rho,
    c,
    v_0,
    U_inf,
    delta_t,
    show=False,
)
plot_lift_over_time(
    circulatory_history,
    rho,
    c,
    v_0,
    U_inf,
    delta_t,
    label="Circulatory",
    title="Normalised Circulatory Lift over Time",
    show=False,
)
plot_lift_frequency_spectrum(
    circulatory_history,
    rho,
    c,
    v_0,
    U_inf,
    delta_t,
    label="Circulatory",
    title="Circulatory Lift Frequency Spectrum",
    show=False,
)
plot_lift_frequency_spectrum(
    total_lift_history,
    rho,
    c,
    v_0,
    U_inf,
    delta_t,
    label="Total",
    title="Total Lift Frequency Spectrum",
    show=False,
)


# Show all generated figures together.
plt.show()
