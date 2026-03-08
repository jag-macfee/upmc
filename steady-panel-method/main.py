from solver import solve_naca_airfoil_camberline
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)  # to allow naca.py import to work
from naca import naca

# --- Running solver and plotting results ---
airfoil_code = "4412"
n_panels = 100
chord_length = 1.0
panel_width = chord_length / n_panels

# Angle of attack of 5 degrees
alpha = np.deg2rad(5)
Q_inf_magnitude = 1.0
Q_inf = np.array([Q_inf_magnitude * np.cos(alpha), Q_inf_magnitude * np.sin(alpha)])

# Get the solution from the solver
gamma, zeta_points, collocation_points = solve_naca_airfoil_camberline(
    airfoil_code, chord_length, n_panels, Q_inf
)

# --- Plot 1: Collocation and Vortex Points ---
airfoil = naca(airfoil_code, chord_length)
x_theory = np.linspace(0, chord_length, 200)
y_theory = airfoil.camber(x_theory)

plt.figure(figsize=(12, 8))
plt.plot(x_theory, y_theory, "k-", label=f"NACA {airfoil_code} Camber Line")
plt.scatter(
    zeta_points[:, 0],
    zeta_points[:, 1],
    marker="o",
    color="blue",
    label="Vortex Points (zeta)",
    s=80,
)
plt.scatter(
    collocation_points[:, 0],
    collocation_points[:, 1],
    marker="x",
    color="red",
    label="Collocation Points",
    s=80,
)
plt.title(f"Vortex and Collocation Points for NACA {airfoil_code}")
plt.xlabel("x/c")
plt.ylabel("y/c")
plt.legend()
plt.grid(True)
plt.axis("equal")
plt.show()

# --- Plot 2: Gamma Distribution ---
plt.figure(figsize=(10, 6))
plt.plot(
    zeta_points[:, 0],
    gamma / panel_width,
    "o-",
    label="Numerical gamma / panel width",
)
plt.title(f"Gamma Distribution for NACA {airfoil_code}")
plt.xlabel("x/c (Vortex Position)")
plt.ylabel("Circulation per unit length")
plt.grid(True)
plt.legend()
plt.show()
