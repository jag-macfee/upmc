import matplotlib.pyplot as plt
import numpy as np


def _get_sorted_gamma_and_xi(gamma, gamma_wake, xi, xi_wake, k):
    """Return sorted combined bound and wake distributions for timestep k."""
    gamma_k = gamma[k, :]
    gamma_wake_k = gamma_wake[k, :k]
    full_gamma = np.concatenate((gamma_k, gamma_wake_k))

    xi_wake_k = xi_wake[k, :k]
    full_xi = np.concatenate((xi, xi_wake_k))

    sort_indices = np.argsort(full_xi)
    return full_xi[sort_indices], full_gamma[sort_indices]


def plot_gamma_surface_over_time(gamma, gamma_wake, xi, xi_wake, k_max, show=True):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    max_xi_points = len(xi) + k_max - 1
    xi_grid = np.zeros((k_max, max_xi_points))
    gamma_grid = np.zeros((k_max, max_xi_points))
    k_grid = np.zeros((k_max, max_xi_points))

    for k in range(k_max):
        full_xi_sorted, full_gamma_sorted = _get_sorted_gamma_and_xi(
            gamma, gamma_wake, xi, xi_wake, k
        )

        num_points = len(full_xi_sorted)
        xi_grid[k, :num_points] = full_xi_sorted
        gamma_grid[k, :num_points] = full_gamma_sorted
        k_grid[k, :num_points] = k + 1

        if num_points < max_xi_points:
            xi_grid[k, num_points:] = full_xi_sorted[-1]
            gamma_grid[k, num_points:] = 0
            k_grid[k, num_points:] = k + 1

    ax.plot_surface(xi_grid, k_grid, gamma_grid, cmap="viridis")

    ax.set_xlabel("xi (Position x/c)")
    ax.set_ylabel("k (Timestep)")
    ax.set_zlabel("gamma / panel_width")
    ax.set_title("Gamma Distribution over Time")

    if show:
        plt.show()


def plot_gamma_distributions_by_timestep(
    gamma, gamma_wake, xi, xi_wake, k_max, show=True
):
    plt.figure()

    timesteps_to_plot = [
        0,
        k_max // 4,
        k_max // 2,
        3 * k_max // 4,
        k_max - 1,
    ]
    labels = ["first", "quarter", "half", "three-quarter", "last"]

    for i, k in enumerate(timesteps_to_plot):
        full_xi_sorted, full_gamma_sorted = _get_sorted_gamma_and_xi(
            gamma, gamma_wake, xi, xi_wake, k
        )

        plt.plot(
            full_xi_sorted,
            full_gamma_sorted,
            label=f"k = {k+1} ({labels[i]})",
        )

    plt.xlabel("xi (x/c)")
    plt.ylabel("Circulation density")
    plt.title("Gamma Distribution at Different Timesteps")
    plt.legend()
    plt.grid(True)
    if show:
        plt.show()


def plot_normalised_lift_components_over_time(
    circulatory_history,
    non_circulatory_history,
    total_history,
    rho,
    c,
    v_0,
    U_inf,
    delta_t,
    show=True,
):
    ss_lift = np.pi * rho * c * v_0 * U_inf
    circulatory_normalised = np.array(circulatory_history) / ss_lift
    non_circulatory_normalised = np.array(non_circulatory_history) / ss_lift
    total_normalised = np.array(total_history) / ss_lift

    num_time_steps = len(total_normalised)
    normalised_time = []
    for k in range(1, num_time_steps + 1):
        t = delta_t * (k - 1)
        chord_lengths_per_second = U_inf / c
        normalised_time.append(t * chord_lengths_per_second)

    plt.figure()
    plt.plot(normalised_time, circulatory_normalised, label="Circulatory")
    plt.plot(normalised_time, non_circulatory_normalised, label="Non-circulatory")
    plt.plot(normalised_time, total_normalised, label="Total")
    plt.xlabel("Normalised time (t U_inf / c)")
    plt.ylabel("Normalised lift")
    plt.title("Normalised Lift Components over Time")
    plt.legend()
    plt.grid(True)
    if show:
        plt.show()


def plot_lift_over_time(
    lift_history,
    rho,
    c,
    v_0,
    U_inf,
    delta_t,
    label="Lift",
    title="Normalised Lift over Time",
    show=True,
):
    ss_lift = np.pi * rho * c * v_0 * U_inf
    lift_normalised = np.array(lift_history) / ss_lift

    num_time_steps = len(lift_normalised)
    normalised_time = []
    for k in range(1, num_time_steps + 1):
        t = delta_t * (k - 1)
        chord_lengths_per_second = U_inf / c
        normalised_time.append(t * chord_lengths_per_second)

    plt.figure()
    plt.plot(normalised_time, lift_normalised, label=label)
    plt.xlabel("Normalised time (t U_inf / c)")
    plt.ylabel("Normalised lift")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    if show:
        plt.show()


def plot_lift_frequency_spectrum(
    h_k,
    rho,
    c,
    v_0,
    U_inf,
    delta_t,
    label="Lift",
    title="Lift Frequency Spectrum",
    show=True,
    show_sears_solution=True,
    sears_fn=None,
):
    h_k = np.asarray(h_k, dtype=float)

    ss_lift = np.pi * rho * c * v_0 * U_inf

    H = np.fft.rfft(h_k)
    freqs = np.fft.rfftfreq(len(h_k), d=delta_t)

    normalised_frequency = freqs * c / U_inf

    H_squared_normalised = (np.abs(H) ** 2) / (ss_lift**2)

    # Remove DC for log axis
    normalised_frequency = normalised_frequency[1:]
    H_squared_normalised = H_squared_normalised[1:]

    valid = (normalised_frequency > 0) & (H_squared_normalised > 0)
    normalised_frequency = normalised_frequency[valid]
    H_squared_normalised = H_squared_normalised[valid]

    # Save only the numerical solution points
    filename = title.lower().replace(" ", "_") + ".csv"
    data_to_save = np.column_stack((normalised_frequency, H_squared_normalised))
    np.savetxt(
        filename,
        data_to_save,
        delimiter=",",
        header="Normalised Frequency,Normalised H squared",
        comments="",
    )

    plt.figure()
    plt.plot(
        normalised_frequency,
        H_squared_normalised,
        "o",
        markersize=3,
        label=label,
    )

    if show_sears_solution:
        # Default: Lysak approximation to |S(k)|^2
        if sears_fn is None:

            def sears_mag_sq_approx(k):
                return (1 / (1 + (2 * np.pi * k) ** 1.3)) ** (1 / 1.3)

            sears_fn = sears_mag_sq_approx

        # Plot Sears all the way to 100 regardless of trusted numerical limit
        x_min_plot = 1e-2
        x_max_plot = 100.0

        sears_freqs = np.logspace(np.log10(x_min_plot), np.log10(x_max_plot), 500)

        # Sears argument is reduced frequency k = pi * f c / U_inf
        k_vals = np.pi * sears_freqs
        sears_squared = sears_fn(k_vals)

        plt.plot(sears_freqs, sears_squared, label="Sears solution")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Normalised frequency (f c / U_inf)")
    plt.ylabel(r"$|H|^2 / (\pi \rho c v_0 U_\infty)^2$")
    plt.xlim(1e-2, 100.0)
    plt.title(title)
    plt.grid(True, which="both")
    plt.legend()

    if show:
        plt.show()
