"""
Project:
    Computational Photonic Crystals

Module:
    P07 - Photonic Band Structure

Description:
    Calculate and visualize the photonic band structure of a
    one-dimensional periodic dielectric using the Bloch condition.
"""

import numpy as np
import matplotlib.pyplot as plt

from utils.bloch_utils import bloch_wavevector


# =============================================================================
# Refractive indices
# =============================================================================

n_1 = 1.0
n_2 = 3.5


# =============================================================================
# Lattice parameters
# =============================================================================

lattice_constant = 1.0
fill_fraction = 0.5


# =============================================================================
# Layer thicknesses
# =============================================================================

d_1 = fill_fraction * lattice_constant
d_2 = (1.0 - fill_fraction) * lattice_constant


# =============================================================================
# Frequency range
# =============================================================================

normalized_frequencies = np.linspace(
    0.0,
    1.0,
    4000,
)

bloch_wavevectors = np.array(
    [
        bloch_wavevector(
            normalized_frequency=frequency,
            n_1=n_1,
            n_2=n_2,
            d_1=d_1,
            d_2=d_2,
            lattice_constant=lattice_constant,
        )
        for frequency in normalized_frequencies
    ]
)

normalized_wavevectors = (
    bloch_wavevectors
    * lattice_constant
    / np.pi
)


# =============================================================================
# Visualization
# =============================================================================

fig, ax = plt.subplots(figsize=(7, 6))

ax.plot(
    normalized_wavevectors,
    normalized_frequencies,
    linewidth=1.5,
)

ax.plot(
    -normalized_wavevectors,
    normalized_frequencies,
    linewidth=1.5,
)

ax.set_xlabel(r"Normalized Bloch wavevector $ka / \pi$")
ax.set_ylabel(r"Normalized frequency $\omega a / (2\pi c)$")
ax.set_title("Photonic Band Structure of a 1D Photonic Crystal")

ax.set_xlim(-1.0, 1.0)
ax.set_ylim(
    normalized_frequencies[0],
    normalized_frequencies[-1],
)

ax.set_xticks(
    [-1.0, 0.0, 1.0],
    [r"$X$", r"$\Gamma$", r"$X$"],
)

ax.grid(alpha=0.3)

fig.tight_layout()

fig.savefig(
    "../figures/p07_photonic_band_structure.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
