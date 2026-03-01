# Helper methods for steady panel method analyses

import numpy as np


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
