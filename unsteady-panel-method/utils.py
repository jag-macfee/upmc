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


# Calculate the potential jump across panel n at time step k given the xi array
# Both n and k are indexed as 1 in the input
def delta_phi_k_n(k, n, gamma):
    # Get the right time step solution
    gamma_k = gamma[k - 1]

    delta_phi_n = sum([gamma_k[i] for i in range(0, n)])
    return delta_phi_n


# Calculate the pressure difference across panel n at time step k
# n and k are 1-indexed again
def delta_p_k_n(rho, k, n, delta_t, gamma):
    d_phi_prev = 0 if k == 1 else delta_phi_k_n(k - 1, n - 1, gamma)
    return rho * (delta_phi_k_n(k, n, gamma) - d_phi_prev) / delta_t


# Calculate the total lift at time step k, accounting for both steady and unsteady effects
# k and n are 1-indexed
# Note that gamma here is the plate gamma solution (not inclusive of the wake)
def lift_components_k(rho, U_inf, k, delta_t, delta_x, gamma):
    gamma_k = gamma[k - 1]
    num_panels = len(gamma_k) - 1

    # Calculate circulatory component
    circulatory_lift = rho * U_inf * sum(gamma_k)

    # Calculate non-circulatory component
    non_circulatory_lift = (
        sum([delta_p_k_n(rho, k, n, delta_t, gamma) for n in range(1, num_panels + 1)])
    ) * delta_x

    total_lift = circulatory_lift + non_circulatory_lift

    return circulatory_lift, non_circulatory_lift, total_lift


def get_lift_history_components(k_max, rho, U_inf, delta_t, delta_x, gamma):
    circulatory_history = []
    non_circulatory_history = []
    total_history = []

    for k in range(1, k_max + 1):
        circulatory_lift, non_circulatory_lift, total_lift = lift_components_k(
            rho, U_inf, k, delta_t, delta_x, gamma
        )
        circulatory_history.append(circulatory_lift)
        non_circulatory_history.append(non_circulatory_lift)
        total_history.append(total_lift)

    return circulatory_history, non_circulatory_history, total_history
