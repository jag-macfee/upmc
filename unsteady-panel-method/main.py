import time
from solver import solve_flat_plate_unsteady
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

U_inf = 1
c = 1
n_panels = 100
k_max = 100
panel_width = c / n_panels

v_gust = lambda t: 1


start_time = time.time()
gamma, gamma_wake, zeta, zeta_wake = solve_flat_plate_unsteady(
    c, n_panels, v_gust, k_max, U_inf
)
end_time = time.time()
print(f"Solver execution time: {end_time - start_time:.4f} seconds")

# 3D plot of gamma distribution over time
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Create a meshgrid for the plot
# We need to create a grid of zeta and k values
# The number of zeta points changes with each timestep due to the wake
max_zeta_points = len(zeta) + k_max - 1
zeta_grid = np.zeros((k_max, max_zeta_points))
gamma_grid = np.zeros((k_max, max_zeta_points))
k_grid = np.zeros((k_max, max_zeta_points))

for k in range(k_max):
    # Combine bound and wake gammas
    gamma_k = gamma[k, :]
    gamma_wake_k = gamma_wake[k, :k]
    full_gamma = np.concatenate((gamma_k, gamma_wake_k))

    # Combine bound and wake zetas
    zeta_wake_k = zeta_wake[k, :k]
    full_zeta = np.concatenate((zeta, zeta_wake_k))

    # Sort by zeta for a clean plot
    sort_indices = np.argsort(full_zeta)
    full_zeta_sorted = full_zeta[sort_indices]
    full_gamma_sorted = full_gamma[sort_indices]

    num_points = len(full_zeta_sorted)
    zeta_grid[k, :num_points] = full_zeta_sorted
    gamma_grid[k, :num_points] = full_gamma_sorted / panel_width
    k_grid[k, :num_points] = k + 1  # Use 1-indexed timesteps for plotting

    # Pad the rest of the row with the last value to make a rectangular grid for the surface plot
    if num_points < max_zeta_points:
        zeta_grid[k, num_points:] = full_zeta_sorted[-1]
        gamma_grid[k, num_points:] = 0  # Assume gamma is zero beyond the wake
        k_grid[k, num_points:] = k + 1


# Plot the surface
ax.plot_surface(zeta_grid, k_grid, gamma_grid, cmap="viridis")

ax.set_xlabel("zeta (Position)")
ax.set_ylabel("k (Timestep)")
ax.set_zlabel("gamma / panel_width")
ax.set_title("Gamma Distribution over Time")

plt.show()

# # Plot gamma distributions
# plt.figure()

# # Timesteps to plot
# timesteps_to_plot = [
#     0,
#     k_max // 4,
#     k_max // 2,
#     3 * k_max // 4,
#     k_max - 1,
# ]
# labels = ["first", "quarter", "half", "three-quarter", "last"]

# for i, k in enumerate(timesteps_to_plot):
#     # Combine bound and wake gammas
#     gamma_k = gamma[k, :]
#     gamma_wake_k = gamma_wake[k, :k]
#     full_gamma = np.concatenate((gamma_k, gamma_wake_k))

#     # Combine bound and wake zetas
#     zeta_wake_k = zeta_wake[k, :k]
#     full_zeta = np.concatenate((zeta, zeta_wake_k))

#     # Sort by zeta for a clean plot
#     sort_indices = np.argsort(full_zeta)
#     full_zeta_sorted = full_zeta[sort_indices]
#     full_gamma_sorted = full_gamma[sort_indices]

#     plt.plot(
#         full_zeta_sorted,
#         full_gamma_sorted / panel_width,
#         label=f"k = {k+1} ({labels[i]})",
#     )


# plt.xlabel("zeta")
# plt.ylabel("Circulation density")
# plt.title("Gamma Distribution at Different Timesteps")
# plt.legend()
# plt.grid(True)
# plt.show()
