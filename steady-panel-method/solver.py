# Helper methods for steady panel method analyses

import numpy as np
from utils import A_flat_plate_steady


# Solve for gamma vector (vortex strengths) for the flat plate case
# Takes in:
# v_o - the gust velocity (velocity component of flow parallel to y axis) in m/s
# c - cord length in m
# n_panels - number of panels (and thus control points)
def solve_flat_plate_spm(v_o, c, n_panels):
    # Vortex positions zeta_n as denoted in Lysak
    zeta = np.array([(c * (n - 1) / n_panels) for n in range(1, n_panels + 2)])

    # Control points, denoted by x_n in Lysak
    x = np.array([c * (n - 0.5) / n_panels for n in range(1, n_panels + 1)])

    # A matrix
    A_mat = A_flat_plate_steady(zeta, x)

    # b vector
    b = np.full(n_panels, -v_o)

    # solve
    gamma = np.linalg.solve(A_mat, b)

    return gamma, zeta, x
