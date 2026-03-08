# Helper methods for steady panel method analyses

import numpy as np
from utils import A_flat_plate_steady, construct_a_matrix, construct_rhs_vector_b
import sys
import os

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)  # to allow naca.py import to work
from naca import naca


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


# Solve for the gamma vector for a NACA airfoil mean camber line (thin airfoil assumption)
# Uses Katz's method
# Takes in:
# airfoil_digits - a string representing the digits of a NACA airfoil (ie. a NACA2412 airfoil will set this to "2412")
# c - the chord length in m
# n_panels - number of panels (and thus collocation points)
# Q_inf - a vector representing the freestream velocity
#
# Note: Katz's solution includes N collocation points and N vortices. This is different to Lysak's flat plate solution,
# which includes a vortex at the trailing edge (TE).
# Since the Kutta condition simply enforces the TE vortex to have gamma = 0, including it as part of the gamma solution will not affect
# any calculations before it.
def solve_naca_airfoil_camberline(airfoil_digits: str, c, n_panels, Q_inf: np.array):
    # Create Airfoil class for use
    airfoil = naca(airfoil_digits, c)

    # Vortex positions
    # Ignore TE vortex for now
    zeta_x = np.array([(c * (n - 1) / n_panels) for n in range(1, n_panels + 1)])
    zeta_z = np.array(airfoil.camber(zeta_x))

    zeta_points = np.column_stack((zeta_x, zeta_z))

    # Collocation point positions
    collocation_x = np.array([c * (n - 0.5) / n_panels for n in range(1, n_panels + 1)])
    collocation_z = np.array(airfoil.camber(collocation_x))

    collocation_points = np.column_stack((collocation_x, collocation_z))

    # Construct 'a' matrix (influence coeff. matrix)
    a_mat = construct_a_matrix(zeta_points, collocation_points)
    rhs_b = construct_rhs_vector_b(Q_inf, zeta_points, collocation_points)

    gamma_sol = np.linalg.solve(a_mat, rhs_b)

    # Add last vortex on the end, with 0 circulation for Kutta condition.
    # Preserve (N, 2) shape when appending the TE point.
    gamma_sol = np.append(gamma_sol, 0)
    zeta_points = np.vstack((zeta_points, np.array([c, 0.0])))

    return gamma_sol, zeta_points, collocation_points
