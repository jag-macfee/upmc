from solver import solve_flat_plate_upm
import matplotlib.pyplot as plt
import numpy as np
from utils import upm_flat_plate_gamma

v_o = 1
c = 1
N = 1000
panel_width = c / N

gamma, zeta, x_abs = solve_flat_plate_upm(v_o, c, N)

# Convert to circulation per unit length
gamma_adjusted = gamma / (v_o * panel_width)
relative_vortex_positions = zeta / c


eps = 1e-6 * c
x = np.linspace(eps, c - eps, 500)  # don't include endpoints (div by zero)
gamma_theoret = upm_flat_plate_gamma(x, c)
x_theoret_rel = x / c

# Plot both on the same axes
plt.scatter(relative_vortex_positions, gamma_adjusted, label="Numerical solution")
plt.plot(x_theoret_rel, gamma_theoret, label="Theoretical solution")

plt.xlabel("Position x/c")
plt.ylabel("Circulation density")
plt.xlim(left=0)
plt.ylim(bottom=-50, top=50)
plt.legend()
plt.show()
