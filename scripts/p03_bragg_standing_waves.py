"""
Project:
    Computational Photonic Crystals

Module:
    P03 - Bragg Scattering and Standing Waves

Description:
    Visualize Bragg scattering at the boundary of the first Brillouin zone
    and the formation of two standing-wave modes from coupled forward and
    backward propagating waves.
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Physical and structural parameters
# =============================================================================

n1 = 1.0
n2 = 3.5
lattice_constant = 1.0
number_of_cells = 6
fill_fraction = 0.5


# =============================================================================
# Spatial domain
# =============================================================================

x = np.linspace(0, number_of_cells * lattice_constant, 2000)


# =============================================================================
# Brillouin-zone boundary wave vector
# =============================================================================

k_bragg = np.pi / lattice_constant


# =============================================================================
# Two standing-wave solutions
# =============================================================================

electric_field_cos = np.cos(k_bragg * x)
electric_field_sin = np.sin(k_bragg * x)


# =============================================================================
# Standing-wave intensities
# =============================================================================

intensity_cos = electric_field_cos ** 2
intensity_sin = electric_field_sin ** 2

# =============================================================================
# Plot standing-wave fields and intensities
# =============================================================================

fig, axes = plt.subplots(
    4,
    1,
    figsize=(10, 11),
    sharex=True,
)

axes[0].plot(
    x,
    electric_field_cos,
    linewidth=2,
    label=r"$\cos(k_{\mathrm{B}}x)$",
)

axes[1].plot(
    x,
    intensity_cos,
    linewidth=2,
    label=r"$\cos^2(k_{\mathrm{B}}x)$",
)

axes[2].plot(
    x,
    electric_field_sin,
    linewidth=2,
    label=r"$\sin(k_{\mathrm{B}}x)$",
)

axes[3].plot(
    x,
    intensity_sin,
    linewidth=2,
    label=r"$\sin^2(k_{\mathrm{B}}x)$",
)

# =============================================================================
# Add the periodic dielectric structure as background
# =============================================================================

for axis in axes:
    for cell_index in range(number_of_cells):
        cell_start = cell_index * lattice_constant
        interface_position = (
            cell_start + fill_fraction * lattice_constant
        )
        cell_end = (cell_index + 1) * lattice_constant

        axis.axvspan(
            cell_start,
            interface_position,
            alpha=0.15,
            label=rf"$n_1={n1}$" if cell_index == 0 else None,
        )

        axis.axvspan(
            interface_position,
            cell_end,
            alpha=0.3,
            label=rf"$n_2={n2}$" if cell_index == 0 else None,
        )


# =============================================================================
# Format the subplots
# =============================================================================

axes[0].set_ylabel(r"$E_{\cos}(x)$")
axes[0].set_title(
    r"Cosine standing wave at $k=\pi/a$"
)

axes[1].set_ylabel(r"$|E_{\cos}(x)|^2$")
axes[1].set_title(
    "Cosine-mode intensity"
)

axes[2].set_ylabel(r"$E_{\sin}(x)$")
axes[2].set_title(
    r"Sine standing wave at $k=\pi/a$"
)

axes[3].set_xlabel(r"Position $x/a$")
axes[3].set_ylabel(r"$|E_{\sin}(x)|^2$")
axes[3].set_title(
    "Sine-mode intensity"
)

for axis in axes:
    axis.set_xlim(0, number_of_cells * lattice_constant)
    axis.grid(alpha=0.25)
    axis.legend(loc="upper right")

axes[0].axhline(0, linewidth=0.8)
axes[2].axhline(0, linewidth=0.8)

axes[0].set_ylim(-1.2, 1.2)
axes[2].set_ylim(-1.2, 1.2)

axes[1].set_ylim(-0.05, 1.15)
axes[3].set_ylim(-0.05, 1.15)

fig.suptitle(
    "Bragg Standing Waves and Field Intensities "
    "in a One-Dimensional Photonic Crystal",
    y=0.975,
)

plt.tight_layout(rect=[0, 0, 1, 0.98])

plt.savefig(
    "../figures/p03_bragg_standing_waves.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()