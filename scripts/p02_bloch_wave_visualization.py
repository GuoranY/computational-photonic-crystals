import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Physical and numerical parameters
# =============================================================================

lattice_constant = 1.0
number_of_cells = 8
points_per_cell = 200

x_min = 0.0
x_max = number_of_cells * lattice_constant

number_of_points = number_of_cells * points_per_cell + 1
x = np.linspace(x_min, x_max, number_of_points)


# =============================================================================
# Bloch wave parameters
# =============================================================================

reciprocal_lattice_vector = 2.0 * np.pi / lattice_constant

wave_vector = 0.6 * np.pi / lattice_constant
periodic_modulation_amplitude = 0.35


# --------------------------------------------------
# Construct the Bloch wave
# --------------------------------------------------

periodic_part = (
    1.0
    + periodic_modulation_amplitude
    * np.cos(reciprocal_lattice_vector * x)
)

plane_wave = np.exp(1j * wave_vector * x)

bloch_wave = periodic_part * plane_wave


# --------------------------------------------------
# Visualize the Bloch wave
# --------------------------------------------------

figure, axes = plt.subplots(
    nrows=3,
    ncols=1,
    figsize=(10, 9),
    sharex=True
)

axes[0].plot(
    x,
    periodic_part,
    label=r"$u_k(x)$"
)

axes[0].set_ylabel(r"$u_k(x)$")
axes[0].set_title("Periodic Part of the Bloch Wave")
axes[0].grid(True, alpha=0.3)
axes[0].legend()

axes[1].plot(
    x,
    np.real(plane_wave),
    label=r"$\mathrm{Re}[e^{ikx}]$"
)

axes[1].set_ylabel("Amplitude")
axes[1].set_title("Real Part of the Plane Wave")
axes[1].grid(True, alpha=0.3)
axes[1].legend()

axes[2].plot(
    x,
    np.real(bloch_wave),
    label=r"$\mathrm{Re}[E_k(x)]$"
)

axes[2].set_xlabel(r"Position $x$")
axes[2].set_ylabel("Amplitude")
axes[2].set_title("Real Part of the Bloch Wave")
axes[2].grid(True, alpha=0.3)
axes[2].legend()

cell_boundaries = np.arange(
    x_min,
    x_max + lattice_constant,
    lattice_constant
)

for axis in axes:
    for boundary in cell_boundaries:
        axis.axvline(
            boundary,
            linestyle="--",
            linewidth=0.8,
            alpha=0.4
        )

plt.tight_layout()

plt.savefig(
    "../figures/p02_bloch_wave_visualization.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
