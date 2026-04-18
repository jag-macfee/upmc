# Helper methods for unsteady panel method analyses

import numpy as np
from utils import A_flat_plate_no_kutta_init, wake_vortex_position


# Solve for gamma vector (vortex strengths) for the flat plate case with no kutta condition
# From Lysak p.18:
#   In a flow where no Horizantal freestream is present, this serves as a valid
#   steady-state solution.
#   However, in situations where U_inf > 0, vortices will be carried into the wake
#   and thus this serves as an initial condition solution
# Takes in:
# v_o - the gust velocity (velocity component of flow parallel to y axis) in m/s
# c - cord length in m
# n_panels - number of panels (and thus control points)
def solve_flat_plate_upm_init(v_o, c, n_panels):
    # Vortex positions xi_n as denoted in Lysak
    xi = np.array([(c * (n - 1) / n_panels) for n in range(1, n_panels + 2)])

    # Control points, denoted by x_n in Lysak
    x = np.array([c * (n - 0.5) / n_panels for n in range(1, n_panels + 1)])

    # A matrix
    A_mat = A_flat_plate_no_kutta_init(xi, x)

    # b vector
    b = np.full(n_panels, -v_o)
    b = np.append(b, 0)

    # solve
    gamma = np.linalg.solve(A_mat, b)

    return gamma, xi, x


# Solve for the gamma vector for k_max number of time steps for flat plate in unsteady flow
# Takes in:
# c - the chord length of the flat plate
# n_panels - number of panels (and thus control points)
# v_gust - a function describing the speed of the gust (purely vertical) over time
#   Expressed as v(t), where t = 0 marks the point where the gust hits the LE
# k_max - the maximum number of time steps to calculate
# U_inf - horizontal freestream velocity
def solve_flat_plate_unsteady(c, n_panels, v_gust, k_max, U_inf):
    # Initial vortex positions xi_n as denoted in Lysak
    xi = np.array([(c * (n - 1) / n_panels) for n in range(1, n_panels + 2)])

    # Control points, denoted by x_n in Lysak
    x = np.array([c * (n - 0.5) / n_panels for n in range(1, n_panels + 1)])

    # A matrix stays the same as before
    A_mat = A_flat_plate_no_kutta_init(xi, x)

    # Set up time and space resolution variables
    panel_width = c / n_panels  # delta x, as denoted in Lysak
    time_step_len = panel_width / U_inf  # time resolution, or delta t

    # Set up v matrix, where v[n - 1][k - 1] denotes the gust strength at control point n at time step k
    # Points are indexed by 1 in Lysak but arrays need to be indexed by 0
    v = np.zeros((n_panels, k_max))
    for n in range(1, n_panels + 1):
        for k in range(1, k_max + 1):
            # Calculate t from number of time steps
            # Note: Time steps are also indexed at 1, thus time step 1 indicates the initial state
            # Hence the k - 1 term
            t = time_step_len * (k - 1)

            # t_n refers to the time it takes for the gust to reach control point n
            # Proof layed out in "Justification for v<k>_n.pdf"
            t_n = x[n - 1] / U_inf

            v[n - 1][k - 1] = v_gust(t - t_n)

    # gamma[k - 1] provides the full gamma vector at time step k
    gamma = np.empty((0, n_panels + 1))

    # gamma_wake[k - 1] will provide a vector of all wake vortex gammas at each time step
    # Note: Shorter gamma_wake lists from earlier time steps will be padded with zeroes
    gamma_wake = np.array([])

    # xi_wake[k - 1] stores wake vortex x-positions corresponding to gamma_wake[k - 1]
    # Same shape as gamma_wake: one row per timestep and padded with zeroes
    xi_wake = np.array([])

    # Initial solution
    b = -np.copy(v[:, 0])
    b = np.append(b, 0)
    gamma = np.vstack((gamma, np.linalg.solve(A_mat, b)))
    gamma_wake = np.zeros((1, k_max - 1))
    xi_wake = np.zeros((1, k_max - 1))

    # March through timesteps k = 2 all the way to k = k_max
    xi_last = xi[-1]  # end of airfoil
    for k in range(2, k_max + 1):
        gamma, gamma_wake, xi_wake = _advance_flat_plate_unsteady_timestep(
            k,
            n_panels,
            x,
            v,
            gamma,
            gamma_wake,
            xi_wake,
            xi_last,
            panel_width,
            A_mat,
        )

    # Return timestepped solutions
    return gamma, gamma_wake, xi, xi_wake


# Helper for function above
# Advances the unsteady solution by one time step
def _advance_flat_plate_unsteady_timestep(
    k, n_panels, x, v, gamma, gamma_wake, xi_wake, xi_last, panel_width, A_mat
):
    # Construct b vector row by row for this timestep
    b = np.array([])
    for n in range(1, n_panels + 1):
        wake_vortex_contribution = 0
        for i in range(1, k):
            xi_i_wake = wake_vortex_position(xi_last, panel_width, i)
            wake_vortex_contribution += gamma[k - i - 1, -1] / (
                2 * np.pi * (xi_i_wake - x[n - 1])
            )

        b = np.append(b, -v[n - 1][k - 1] - wake_vortex_contribution)

    wake_vortex_sum = sum([gamma[k - i - 1, -1] for i in range(1, k)])
    b = np.append(b, -wake_vortex_sum)

    # Solve current timestep and update wake strengths/positions
    last_TE_gamma = gamma[-1, -1]
    gamma = np.vstack((gamma, np.linalg.solve(A_mat, b)))

    prev_wake_row = gamma_wake[-1]
    new_wake_row = np.empty_like(prev_wake_row)
    new_wake_row[0] = last_TE_gamma
    new_wake_row[1:] = prev_wake_row[:-1]
    gamma_wake = np.vstack((gamma_wake, new_wake_row))

    prev_xi_row = xi_wake[-1]
    new_xi_row = np.copy(prev_xi_row)

    # k is 1-indexed for timesteps, and we are in timestep k.
    # The new vortex is the (k-1)th vortex shed into the wake.
    # Its position is fixed relative to the trailing edge.
    new_xi_row[k - 2] = wake_vortex_position(xi_last, panel_width, k - 1)
    xi_wake = np.vstack((xi_wake, new_xi_row))

    return gamma, gamma_wake, xi_wake
