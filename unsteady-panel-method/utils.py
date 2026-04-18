import numpy as np


# Theoretical solution for circulation per unit length
# For (u)nsteady (p)anel (m)ethod flat plate
# Normalised for v_o
def upm_flat_plate_gamma(x, c):
    return 2 * (1 - 2 * x / c) / np.sqrt(1 - (1 - 2 * x / c) ** 2)


# Gets the A matrix associated with the removed kutta condition steady state solution
# Lysak p.18
def A_flat_plate_no_kutta_init(xi, x):
    n_vortices = len(xi)
    num_panels = len(x)

    A_mat = np.array([calculate_A_row(xi, x, i) for i in range(1, num_panels + 1)])
    A_mat = np.vstack([A_mat, np.full(n_vortices, 1)])

    return (1 / (2 * np.pi)) * A_mat


# Returns the nth row of the A matrix (helper for fn above)
# Takes in xi, x, and n
# Note: n is 1-indexed in input, so we adjust it in the calc
def calculate_A_row(xi, x, n):
    num_vortices = len(xi)

    return [1 / (xi[i] - x[n - 1]) for i in range(0, num_vortices)]


# Gets xi of the ith vortex in the wake
def wake_vortex_position(xi_last, panel_width, i):
    return xi_last + panel_width * i
