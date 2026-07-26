"""
Project:
    Computational Photonic Crystals

Module:
    P04 - Transfer Matrix and Transmission Spectrum

Description:
    Calculate and visualize the transmission spectrum of a finite
    one-dimensional photonic crystal using the transfer matrix method.
"""

import numpy as np
import matplotlib.pyplot as plt

from utils.transmission_utils import transmission_spectrum

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
number_of_cells = 10


# =============================================================================
# Layer thicknesses
# =============================================================================

d_1 = fill_fraction * lattice_constant
d_2 = (1.0 - fill_fraction) * lattice_constant


# =============================================================================
# Surrounding media
# =============================================================================

n_incident = 1.0
n_exit = 1.0


if __name__ == "__main__":
    # =========================================================================
    # Frequency range
    # =========================================================================

    normalized_frequencies = np.linspace(0.0, 1.0, 2000)

    transmissions = transmission_spectrum(
        normalized_frequencies=normalized_frequencies,
        n_1=n_1,
        n_2=n_2,
        d_1=d_1,
        d_2=d_2,
        lattice_constant=lattice_constant,
        number_of_cells=number_of_cells,
        n_incident=n_incident,
        n_exit=n_exit,
    )

    # ============================================================================
    # Visualization
    # ============================================================================

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(
        normalized_frequencies,
        transmissions,
        linewidth=1.5,
    )

    ax.set_xlabel(r"Normalized frequency $\omega a / (2\pi c)$")
    ax.set_ylabel("Transmission")
    ax.set_title("Transmission Spectrum of a 1D Photonic Crystal")

    ax.set_xlim(normalized_frequencies[0], normalized_frequencies[-1])
    ax.set_ylim(0.0, 1.05)

    ax.grid(alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        "../figures/p04_transmission_spectrum.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
