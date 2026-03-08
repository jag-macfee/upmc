# Helper methods for steady panel method analyses

import numpy as np


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
def normal_vec_to_panel(p1: tuple, p2: tuple):
    (x1, z1) = p1
    (x2, z2) = p2

    panel_vec = [x2 - x1, z2 - z1]
    n = np.array([-panel_vec[1], panel_vec[0]])

    return 1 / (np.sqrt(n.dot(n))) * n
