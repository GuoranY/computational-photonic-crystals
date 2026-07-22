"""
Project:
    Computational Photonic Crystals

Module:
    P01 - One-Dimensional Periodic Dielectric

Description:
    Construct and visualize a one-dimensional periodic dielectric
    structure as the foundation for photonic crystal simulations.
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Physical parameters
# =============================================================================

n1 = 1.0                  # Refractive index of material A
n2 = 3.5                  # Refractive index of material B

lattice_constant = 1.0    # Lattice constant (normalized units)
num_periods = 10          # Number of unit cells

fill_fraction = 0.5       # Fraction of each unit cell occupied by material A

num_points = 1000         # Spatial sampling points

# =============================================================================
# Spatial grid
# =============================================================================

total_length = lattice_constant * num_periods

x = np.linspace(
    start=0.0,
    stop=total_length,
    num=num_points,
    endpoint=False,
)

# =============================================================================
# Periodic refractive-index profile
# =============================================================================

position_in_cell = np.mod(x, lattice_constant)

refractive_index = np.where(
    position_in_cell < fill_fraction * lattice_constant,
    n1,
    n2,
)

# =============================================================================
# Visualization
# =============================================================================

figure, axis = plt.subplots(
    figsize=(10, 3),
    dpi=300,
)

axis.step(
    x,
    refractive_index,
    where="post",
    linewidth=2.0,
)

axis.set_xlabel(
    r"Position ($x/a$)"
)
axis.set_ylabel(
    r"Refractive Index $n(x)$"
)

axis.set_title(
    "Periodic Refractive-Index Profile"
)

axis.grid(
    visible=True,
)

figure.tight_layout()

figure.savefig(
    "../figures/p01_periodic_dielectric.png",
    dpi=300,
    bbox_inches="tight",
)

plt.show()
