from solver import solve_flat_plate_spm
import matplotlib.pyplot as plt
import numpy as np
from utils import spm_flat_plate_gamma
import sys
import os

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)  # to allow naca.py import to work
from naca import naca

v_o = 1
c = 1
N = 1000
panel_width = c / N

# gamma, zeta, x_abs = solve_flat_plate_spm(v_o, c, N)

# # Convert to circulation per unit length
# gamma_adjusted = gamma / (v_o * panel_width)
# rel_vortex_positions = zeta[:-1] / c  # Remove last zeta as we are ignoring

# eps = 1e-6 * c
# x = np.linspace(eps, c, 500)  # includes endpoint by default
# gamma_theoret = spm_flat_plate_gamma(x, c)
# x_theoret_rel = x / c

# # Plot both on the same axes
# plt.scatter(rel_vortex_positions, gamma_adjusted, label="Numerical solution")
# plt.plot(x_theoret_rel, gamma_theoret, label="Theoretical solution")

# plt.xlabel("Position x/c")
# plt.ylabel("Circulation density")
# plt.xlim(left=0)
# plt.ylim(bottom=0, top=50)
# plt.legend()
# plt.show()

# --- Plotting NACA airfoil camber lines ---
airfoil_codes = ["2412", "2415", "4412"]
c = 1.0  # Chord length
x_points = np.linspace(0, c, 200)

plt.figure(figsize=(10, 6))

for code in airfoil_codes:
    airfoil = naca(code, c)
    y_camber = [airfoil.camber(x) for x in x_points]
    plt.plot(x_points, y_camber, label=f"NACA {code}")

plt.title("NACA 4-Digit Airfoil Mean Camber Lines")
plt.xlabel("x/c (Position along chord)")
plt.ylabel("y/c (Camber)")
plt.legend()
plt.grid(True)
plt.axis("equal")
plt.show()
