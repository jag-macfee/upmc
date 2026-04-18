import time
from solver import solve_flat_plate_unsteady
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys

U_inf = 4
c = 2
n_panels = 50
k_max = n_panels  # Length of the wake in terms of panels
panel_width = c / n_panels
delta_t = panel_width / U_inf


# Have to make sure dirac is shifted for x_n rather than gamma_n
# Note: Due to floating-point precision in calculation of t and delta_t,
# we have to acknowledge their difference to be system eps
def dirac(t):
    if abs(t) <= sys.float_info.epsilon:
        return 1
    return 0


# Have gusts with magnitude one for now, we will sort out normalising later
# will need to normalise with respect to each xi position
half_chord_per_cycle_omega = 2 * np.pi / (0.5 * c / U_inf)
v_gust = (
    lambda t: (np.e / (np.e**2 - 1)) * np.exp(np.cos(half_chord_per_cycle_omega * t))
    - 1 / np.e
)

start_time = time.time()
gamma, gamma_wake, xi, xi_wake = solve_flat_plate_unsteady(
    c, n_panels, v_gust, k_max, U_inf
)
end_time = time.time()
print(f"Solver execution time: {end_time - start_time:.4f} seconds")

# normalise
gamma = gamma / panel_width
gamma_wake = gamma_wake / panel_width
xi = xi / c
xi_wake = xi_wake / c

# 3D plot of gamma distribution over time
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Create a meshgrid for the plot
# We need to create a grid of xi and k values
# The number of xi points changes with each timestep due to the wake
max_xi_points = len(xi) + k_max - 1
xi_grid = np.zeros((k_max, max_xi_points))
gamma_grid = np.zeros((k_max, max_xi_points))
k_grid = np.zeros((k_max, max_xi_points))

for k in range(k_max):
    # Combine bound and wake gammas
    gamma_k = gamma[k, :]
    gamma_wake_k = gamma_wake[k, :k]
    full_gamma = np.concatenate((gamma_k, gamma_wake_k))

    # Combine bound and wake xis
    xi_wake_k = xi_wake[k, :k]
    full_xi = np.concatenate((xi, xi_wake_k))

    # Sort by xi for a clean plot
    sort_indices = np.argsort(full_xi)
    full_xi_sorted = full_xi[sort_indices]
    full_gamma_sorted = full_gamma[sort_indices]

    num_points = len(full_xi_sorted)
    xi_grid[k, :num_points] = full_xi_sorted
    gamma_grid[k, :num_points] = full_gamma_sorted
    k_grid[k, :num_points] = k + 1  # Use 1-indexed timesteps for plotting

    # Pad the rest of the row with the last value to make a rectangular grid for the surface plot
    if num_points < max_xi_points:
        xi_grid[k, num_points:] = full_xi_sorted[-1]
        gamma_grid[k, num_points:] = 0  # Assume gamma is zero beyond the wake
        k_grid[k, num_points:] = k + 1


# Plot the surface
ax.plot_surface(xi_grid, k_grid, gamma_grid, cmap="viridis")

ax.set_xlabel("xi (Position x/c)")
ax.set_ylabel("k (Timestep)")
ax.set_zlabel("gamma / panel_width")
ax.set_title("Gamma Distribution over Time")

plt.show()

# Plot gamma distributions
plt.figure()

# Timesteps to plot
timesteps_to_plot = [
    0,
    k_max // 4,
    k_max // 2,
    3 * k_max // 4,
    k_max - 1,
]
labels = ["first", "quarter", "half", "three-quarter", "last"]

for i, k in enumerate(timesteps_to_plot):
    # Combine bound and wake gammas
    gamma_k = gamma[k, :]
    gamma_wake_k = gamma_wake[k, :k]
    full_gamma = np.concatenate((gamma_k, gamma_wake_k))

    # Combine bound and wake xis
    xi_wake_k = xi_wake[k, :k]
    full_xi = np.concatenate((xi, xi_wake_k))

    # Sort by xi for a clean plot
    sort_indices = np.argsort(full_xi)
    full_xi_sorted = full_xi[sort_indices]
    full_gamma_sorted = full_gamma[sort_indices]

    plt.plot(
        full_xi_sorted,
        full_gamma_sorted,
        label=f"k = {k+1} ({labels[i]})",
    )


plt.xlabel("xi (x/c)")
plt.ylabel("Circulation density")
plt.title("Gamma Distribution at Different Timesteps")
plt.legend()
plt.grid(True)
plt.show()
