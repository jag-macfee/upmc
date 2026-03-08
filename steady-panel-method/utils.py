# Helper methods for steady panel method analyses

import numpy as np


##############################################################################################
##############################################################################################
############################## FLAT PLATE FUNCTIONS ##########################################
##############################################################################################
##############################################################################################


# Theoretical solution for circulation per unit length
# For (s)teady (p)anel (m)ethod flat plate
# Normalised for v_o
def spm_flat_plate_gamma(x, c):
    return 2 * np.sqrt((1 - x / c) / (x / c))


# Gives the magnitude of the velocity vector induced by a point vertex of strength gamma
def vort_vel(gamma, r):
    return gamma / (2 * np.pi * r)


# Gives the A matrix in the flat plate steady panel case
# Takes in zeta and x vectors
def A_flat_plate_steady(zeta, x):
    num_panels = len(x)
    A_mat = np.array([calculate_A_row(zeta, x, i) for i in range(1, num_panels + 1)])

    return (1 / (2 * np.pi)) * A_mat


# Returns the nth row of the A matrix (helper for fn above)
# Takes in zeta, x, and n
# Note: n is 1-indexed in input, so we adjust it in the calc
def calculate_A_row(zeta, x, n):
    num_vortices = len(zeta)

    # Kutta condition allows us to ignore last vortex in zeta arr
    return [1 / (zeta[i] - x[n - 1]) for i in range(0, num_vortices - 1)]


##############################################################################################
##############################################################################################
############################## MEAN CAMBER FUNCTIONS #########################################
##############################################################################################
##############################################################################################


# Given a vortex strength gamma, the coordinates of the vortex (x_j, z_j),
# calculates the velocity (u, w) of a point P(x, z) induced by the vortex
# Katz p. 223, matrix equation
# As used in Katz, will call the method vor2D
def vor2D(gamma, x, z, x_j, z_j):
    r_squared = (x - x_j) ** 2 + (z - z_j) ** 2

    # Might move this to global var to avoid redeclaring
    rhs_mat = np.array([[0, 1], [-1, 0]])
    rhs_vec = np.array([x - x_j, z - z_j])

    v = gamma / (2 * np.pi * r_squared) * (rhs_mat @ rhs_vec)
    return v


# Given a panel defined by two endpoints (vortex positions), return a unit vector
# normal to it
def normal_vec_to_panel(p1: np.array, p2: np.array):
    (x1, z1) = p1
    (x2, z2) = p2

    panel_vec = [x2 - x1, z2 - z1]
    n = np.array([-panel_vec[1], panel_vec[0]])

    return 1 / (np.sqrt(n.dot(n))) * n


# Given a set of vortex and collocation point coordinate lists,
# create the 'a' matrix as used in Katz
# This is the 'influence coefficient matrix' which will be used to solve for the linear system
# and represents summations of induced vortex velocities at each collocation point
def construct_a_matrix(zeta_points, collocation_points):
    n_vortices = len(zeta_points)

    a_mat = np.array(
        [
            construct_a_matrix_row(i, zeta_points, collocation_points)
            for i in range(1, n_vortices + 1)
        ]
    )

    return a_mat


# Returns the ith row of the 'a' matrix
# Note: i is passed in 1-indexed so we adjust for appropriate position
def construct_a_matrix_row(i, zeta_points, collocation_points):
    # Panels can be defined with 2 points: vortex i and collocation point i
    vortex_i = zeta_points[i - 1]
    collocation_point_i = collocation_points[i - 1]

    n_i = normal_vec_to_panel(vortex_i, collocation_point_i)

    # Construct row left to right, iterating over all vortices (in zeta) j
    x_ci = collocation_point_i[0]
    z_ci = collocation_point_i[1]

    row = []
    for vortex_point in zeta_points:
        x_0j = vortex_point[0]
        z_0j = vortex_point[1]

        induced_vel_per_circulation = vor2D(1, x_ci, z_ci, x_0j, z_0j)
        a_ij = np.dot(induced_vel_per_circulation, n_i)

        row.append(a_ij)

    return np.array(row)


# Constructs the RHS vector in the steady panel method matrix eq.
# Represents the component of the freestream velocity normal to each panel at
# each collocation point i
def construct_rhs_vector_b(Q_inf, zeta_points, collocation_points):
    # Construct entries one by one
    b = []

    n_collocation_points = len(collocation_points)
    for i in range(1, n_collocation_points + 1):
        # Panels can be defined with 2 points: vortex i and collocation point i
        vortex_i = zeta_points[i - 1]
        collocation_point_i = collocation_points[i - 1]

        n_i = normal_vec_to_panel(vortex_i, collocation_point_i)
        dot_prod = np.dot(-Q_inf, n_i)

        b.append(dot_prod)

    return np.array(b)
