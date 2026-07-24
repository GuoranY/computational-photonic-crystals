"""
Project:
    Computational Photonic Crystals

Module:
    P05 - Bloch Condition and Band Classification

Description:
    Apply the Bloch condition to the unit-cell transfer matrix
    and identify allowed bands and forbidden frequency gaps.
"""

import numpy as np
import matplotlib.pyplot as plt

from utils.bloch_utils import bloch_function


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

normalized_frequencies = np.linspace(0.0, 1.0, 2000)

bloch_values = np.array(
    [
        bloch_function(
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


# =============================================================================
# Band classification
# =============================================================================

allowed_mask = np.abs(bloch_values) <= 1.0
forbidden_mask = np.abs(bloch_values) > 1.0


# =============================================================================
# Visualization
# =============================================================================

fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(
    normalized_frequencies,
    bloch_values,
    linewidth=1.5,
    label=r"$F(\omega)$",
)

ax.axhline(
    1.0,
    linestyle="--",
    linewidth=1.0,
    label=r"$F=1$",
)

ax.axhline(
    -1.0,
    linestyle="--",
    linewidth=1.0,
    label=r"$F=-1$",
)

ax.set_xlabel(r"Normalized frequency $\omega a / (2\pi c)$")
ax.set_ylabel(r"$F(\omega)$")
ax.set_title("Bloch Function of a 1D Photonic Crystal")

ax.set_xlim(
    normalized_frequencies[0],
    normalized_frequencies[-1],
)

ax.grid(alpha=0.3)
ax.legend()

fig.tight_layout()

fig.savefig(
    "../figures/p05_bloch_function.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()

fig, ax = plt.subplots(figsize=(9, 2.5))

ax.fill_between(
    normalized_frequencies,
    0.0,
    1.0,
    color="tab:blue",
    where=allowed_mask,
    alpha=0.5,
    label="Allowed bands",
)

ax.fill_between(
    normalized_frequencies,
    0.0,
    1.0,
    color="silver",
    where=forbidden_mask,
    alpha=0.5,
    label="Forbidden gaps",
)

ax.set_xlim(
    normalized_frequencies[0],
    normalized_frequencies[-1],
)

ax.set_ylim(0.0, 1.0)

ax.set_xlabel(r"Normalized frequency $\omega a / (2\pi c)$")
ax.set_yticks([])

ax.set_title("Allowed Bands and Forbidden Gaps")

ax.legend()

fig.tight_layout()

fig.savefig(
    "../figures/p05_allowed_forbidden_bands.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
