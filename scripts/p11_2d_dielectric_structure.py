"""
Project:
    Computational Photonic Crystals

Module:
    P08 - Two-Dimensional Periodic Dielectric

Description:
    Construct and visualize a two-dimensional square lattice of
    dielectric rods in air as the foundation for two-dimensional
    photonic crystal band-structure calculations.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ============================================================================
# Physical parameters
# ============================================================================

lattice_constant = 1.0
rod_radius = 0.2 * lattice_constant

background_index = 1.0
rod_index = 3.5

# ============================================================================
# Structure parameters
# ============================================================================

number_of_cells_x = 6
number_of_cells_y = 6
points_per_cell = 100

# ============================================================================
# Spatial domain
# ============================================================================
x_min = 0.0
x_max = number_of_cells_x * lattice_constant

y_min = 0.0
y_max = number_of_cells_y * lattice_constant

number_of_points_x = number_of_cells_x * points_per_cell
number_of_points_y = number_of_cells_y * points_per_cell

x = np.linspace(
    x_min,
    x_max,
    number_of_points_x,
    endpoint=False,
)

y = np.linspace(
    y_min,
    y_max,
    number_of_points_y,
    endpoint=False,
)

X, Y = np.meshgrid(x, y)

# ============================================================================
# Initialize the entire structure as background material
# ============================================================================

refractive_index = background_index * np.ones_like(X)


# ============================================================================
# Place one dielectric rod at the center of each unit cell
# ============================================================================

for cell_x in range(number_of_cells_x):
    for cell_y in range(number_of_cells_y):
        rod_center_x = (cell_x + 0.5) * lattice_constant
        rod_center_y = (cell_y + 0.5) * lattice_constant

        distance_squared = (
            (X - rod_center_x) ** 2
            + (Y - rod_center_y) ** 2
        )

        inside_rod = distance_squared <= rod_radius ** 2
        refractive_index[inside_rod] = rod_index

# ============================================================================
# Convert refractive index to relative permittivity
# ============================================================================

relative_permittivity = refractive_index ** 2

# ============================================================================
# Plot the two-dimensional refractive-index distribution
# ============================================================================

fig, ax = plt.subplots(figsize=(7, 6))

image = ax.pcolormesh(
    X,
    Y,
    refractive_index,
    shading="auto",
)

ax.set_xlabel(r"$x/a$")
ax.set_ylabel(r"$y/a$")
ax.set_title("Two-Dimensional Periodic Dielectric")

ax.set_aspect("equal")
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

colorbar = fig.colorbar(image, ax=ax)
colorbar.set_label("Refractive Index")

fig.tight_layout()

project_root = Path(__file__).resolve().parents[1]
figure_directory = project_root / "figures"
figure_directory.mkdir(parents=True, exist_ok=True)

output_path = figure_directory / "p11_2d_dielectric_structure.png"
fig.savefig(output_path, dpi=300, bbox_inches="tight")

plt.show()
