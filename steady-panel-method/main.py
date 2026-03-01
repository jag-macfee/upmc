from solver import solve_flat_plate
import matplotlib.pyplot as plt
import numpy as np

v_o = 1
c = 1
N = 50
panel_width = c / N

gamma, zeta, x_abs = solve_flat_plate(v_o, c, N)

# Convert to circulation per unit length
gamma_adjusted = gamma / (v_o * panel_width)
x_rel = x_abs / c

# Theoretical solution
gamma_theoretical = lambda x: 2 * np.sqrt((1 - x / c) / (x / c))

eps = 1e-6 * c
x = np.linspace(eps, c, 500)  # includes endpoint by default
gamma_theoret = gamma_theoretical(x) / v_o
x_theoret_rel = x / c

# Plot both on the same axes
plt.scatter(x_rel, gamma_adjusted, label="Numerical solution")
plt.plot(x_theoret_rel, gamma_theoret, label="Theoretical solution")

plt.xlabel("Position x/c")
plt.ylabel("Circulation density")
plt.xlim(left=0)
plt.ylim(bottom=0, top=50)
plt.legend()
plt.show()
