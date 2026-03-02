from solver import solve_flat_plate_spm
import matplotlib.pyplot as plt
import numpy as np
from utils import spm_flat_plate_gamma

v_o = 1
c = 1
N = 1000
panel_width = c / N

gamma, zeta, x_abs = solve_flat_plate_spm(v_o, c, N)

# Convert to circulation per unit length
gamma_adjusted = gamma / (v_o * panel_width)
rel_vortex_positions = zeta[:-1] / c  # Remove last zeta as we are ignoring

eps = 1e-6 * c
x = np.linspace(eps, c, 500)  # includes endpoint by default
gamma_theoret = spm_flat_plate_gamma(x, c)
x_theoret_rel = x / c

# Plot both on the same axes
plt.scatter(rel_vortex_positions, gamma_adjusted, label="Numerical solution")
plt.plot(x_theoret_rel, gamma_theoret, label="Theoretical solution")

plt.xlabel("Position x/c")
plt.ylabel("Circulation density")
plt.xlim(left=0)
plt.ylim(bottom=0, top=50)
plt.legend()
plt.show()
