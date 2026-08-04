"""
Project:
    Computational Photonic Crystals

Module:
    P18 - Photonic Crystal Waveguide

Description:
    Construct a two-dimensional line-defect photonic-crystal
    waveguide, calculate its TM projected band structure, and
    identify a guided mode from its transverse field confinement.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utils.field_reconstruction_utils import (
    field_confinement_factor,
    reconstruct_tm_field,
)
from utils.waveguide_utils import (
    calculate_tm_waveguide_band_structure,
    rectangular_reciprocal_vector_set,
    square_lattice_waveguide_positions,
    supercell_convolution_matrix,
    waveguide_reciprocal_vectors,
    waveguide_supercell_rod_positions,
    waveguide_supercell_vectors,
    waveguide_wavevector_path,
)


# ============================================================================
# Structural parameters
# ============================================================================

lattice_constant = 1.0

background_refractive_index = 1.0
rod_refractive_index = 3.5

background_permittivity = (
    background_refractive_index**2
)

rod_permittivity = (
    rod_refractive_index**2
)

rod_radius = (
    0.2 * lattice_constant
)

number_of_rows = 7
number_of_columns = 7


# ============================================================================
# Plane-wave expansion parameters
# ============================================================================

reciprocal_index_limit_x = 3
reciprocal_index_limit_y = 10

number_of_wavevector_points = 41
number_of_bands = 24


# ============================================================================
# Guided-mode search parameters
# ============================================================================

gap_lower = 0.278188
gap_upper = 0.415995

confinement_width = (
    lattice_constant
)

selected_wavevector_index = (
    number_of_wavevector_points // 2
)


# ============================================================================
# Field-reconstruction grid
# ============================================================================

x_grid = np.linspace(
    -3.0 * lattice_constant,
    3.0 * lattice_constant,
    300
)

y_grid = np.linspace(
    -3.5 * lattice_constant,
    3.5 * lattice_constant,
    350
)


# ============================================================================
# Construct waveguide geometry for visualization
# ============================================================================

visualization_rod_positions = (
    square_lattice_waveguide_positions(
        lattice_constant,
        number_of_rows,
        number_of_columns
    )
)


# ============================================================================
# Construct computational waveguide supercell
# ============================================================================

real_vector_1, real_vector_2 = (
    waveguide_supercell_vectors(
        lattice_constant,
        number_of_rows
    )
)

reciprocal_vector_1, reciprocal_vector_2 = (
    waveguide_reciprocal_vectors(
        lattice_constant,
        number_of_rows
    )
)

supercell_rod_positions = (
    waveguide_supercell_rod_positions(
        lattice_constant,
        number_of_rows
    )
)

supercell_area = (
    np.linalg.norm(real_vector_1)
    * np.linalg.norm(real_vector_2)
)


# ============================================================================
# Construct wavevector path along the waveguide
# ============================================================================

wavevectors, normalized_wavevectors = (
    waveguide_wavevector_path(
        lattice_constant,
        number_of_wavevector_points
    )
)


# ============================================================================
# Construct reciprocal-vector basis
# ============================================================================

reciprocal_vectors = (
    rectangular_reciprocal_vector_set(
        reciprocal_vector_1,
        reciprocal_vector_2,
        reciprocal_index_limit_x,
        reciprocal_index_limit_y
    )
)


# ============================================================================
# Construct permittivity convolution matrix
# ============================================================================

permittivity_convolution_matrix = (
    supercell_convolution_matrix(
        reciprocal_vectors,
        supercell_rod_positions,
        rod_radius,
        rod_permittivity,
        background_permittivity,
        supercell_area
    )
)

hermitian_error = np.max(
    np.abs(
        permittivity_convolution_matrix
        - permittivity_convolution_matrix.conj().T
    )
)


# ============================================================================
# Calculate TM projected band structure
# ============================================================================

tm_frequencies, tm_eigenvectors = (
    calculate_tm_waveguide_band_structure(
        wavevectors,
        reciprocal_vectors,
        permittivity_convolution_matrix,
        lattice_constant,
        number_of_bands
    )
)


# ============================================================================
# Identify candidate modes inside the selected frequency interval
# ============================================================================

selected_wavevector = (
    wavevectors[
        selected_wavevector_index
    ]
)

candidate_bands = []

for band_index in range(
    number_of_bands
):
    frequency = (
        tm_frequencies[
            selected_wavevector_index,
            band_index
        ]
    )

    if (
        gap_lower
        <
        frequency
        <
        gap_upper
    ):
        candidate_bands.append(
            band_index
        )

if not candidate_bands:
    raise RuntimeError(
        "No candidate TM waveguide mode was found "
        "inside the selected frequency interval."
    )


# ============================================================================
# Evaluate transverse confinement of candidate modes
# ============================================================================

candidate_confinements = []

for band_index in candidate_bands:

    candidate_coefficients = (
        tm_eigenvectors[
            selected_wavevector_index,
            :,
            band_index
        ]
    )

    candidate_field = reconstruct_tm_field(
        candidate_coefficients,
        reciprocal_vectors,
        selected_wavevector,
        x_grid,
        y_grid
    )

    candidate_intensity = (
        np.abs(candidate_field)**2
    )

    confinement = (
        field_confinement_factor(
            candidate_intensity,
            y_grid,
            confinement_width
        )
    )

    candidate_confinements.append(
        confinement
    )


# ============================================================================
# Select the most strongly confined candidate
# ============================================================================

best_candidate_index = int(
    np.argmax(
        candidate_confinements
    )
)

selected_band = (
    candidate_bands[
        best_candidate_index
    ]
)

selected_frequency = (
    tm_frequencies[
        selected_wavevector_index,
        selected_band
    ]
)

selected_confinement = (
    candidate_confinements[
        best_candidate_index
    ]
)


# ============================================================================
# Reconstruct selected guided mode
# ============================================================================

selected_mode_coefficients = (
    tm_eigenvectors[
        selected_wavevector_index,
        :,
        selected_band
    ]
)

selected_field = reconstruct_tm_field(
    selected_mode_coefficients,
    reciprocal_vectors,
    selected_wavevector,
    x_grid,
    y_grid
)

selected_intensity = (
    np.abs(selected_field)**2
)

maximum_intensity = np.max(
    selected_intensity
)

if maximum_intensity > 0.0:
    selected_intensity /= (
        maximum_intensity
    )


# ============================================================================
# Print computational results
# ============================================================================

print(
    "Waveguide supercell vectors:"
)

print(
    "a1 =",
    real_vector_1
)

print(
    "a2 =",
    real_vector_2
)

print(
    "Waveguide reciprocal vectors:"
)

print(
    "b1 =",
    reciprocal_vector_1
)

print(
    "b2 =",
    reciprocal_vector_2
)

print(
    "Number of supercell rods =",
    supercell_rod_positions.shape[0]
)

print(
    "Number of plane waves =",
    reciprocal_vectors.shape[0]
)

print(
    "Permittivity-matrix Hermitian error =",
    hermitian_error
)

print(
    "Candidate waveguide bands =",
    candidate_bands
)

for band_index, confinement in zip(
    candidate_bands,
    candidate_confinements
):
    print(
        "Band",
        band_index,
        "frequency =",
        tm_frequencies[
            selected_wavevector_index,
            band_index
        ],
        "confinement =",
        confinement
    )

print(
    "Selected band =",
    selected_band
)

print(
    "Selected frequency =",
    selected_frequency
)

print(
    "Selected confinement =",
    selected_confinement
)


# ============================================================================
# Plot waveguide structure
# ============================================================================

structure_figure, structure_axis = (
    plt.subplots(
        figsize=(7.0, 4.5)
    )
)

structure_axis.scatter(
    visualization_rod_positions[:, 0]
    / lattice_constant,
    visualization_rod_positions[:, 1]
    / lattice_constant,
    s=120
)

structure_axis.axhspan(
    -0.5,
    0.5,
    alpha=0.10,
    label="Line defect"
)

structure_axis.set_aspect(
    "equal"
)

structure_axis.set_xlabel(
    r"$x/a$"
)

structure_axis.set_ylabel(
    r"$y/a$"
)

structure_axis.set_title(
    "Two-Dimensional Photonic-Crystal Waveguide"
)

structure_axis.grid(
    alpha=0.25
)

structure_figure.tight_layout()

structure_output_path = Path(
    "../figures/"
    "p18_waveguide_structure.png"
)

structure_figure.savefig(
    structure_output_path,
    dpi=300,
    bbox_inches="tight"
)


# ============================================================================
# Plot TM projected band structure
# ============================================================================

band_figure, band_axis = (
    plt.subplots(
        figsize=(7.0, 5.0)
    )
)

band_axis.axhspan(
    gap_lower,
    gap_upper,
    alpha=0.12,
    label="Selected TM gap interval"
)

for band_index in range(
    number_of_bands
):

    if band_index == selected_band:
        band_axis.plot(
            normalized_wavevectors,
            tm_frequencies[
                :,
                band_index
            ],
            linewidth=2.5,
            label="Guided-mode band"
        )

    else:
        band_axis.plot(
            normalized_wavevectors,
            tm_frequencies[
                :,
                band_index
            ],
            linewidth=0.9,
            alpha=0.50
        )

band_axis.set_xlim(
    0.0,
    1.0
)

band_axis.set_ylim(
    0.0,
    0.65
)

band_axis.set_xticks(
    [
        0.0,
        1.0
    ]
)

band_axis.set_xticklabels(
    [
        r"$\Gamma_{\mathrm{w}}$",
        r"$X_{\mathrm{w}}$"
    ]
)

band_axis.set_xlabel(
    r"Wavevector $k_xa/\pi$"
)

band_axis.set_ylabel(
    r"Normalized frequency "
    r"$\omega a/2\pi c$"
)

band_axis.set_title(
    "TM Projected Band Structure of a "
    "Photonic-Crystal Waveguide"
)

band_axis.grid(
    alpha=0.25
)

band_axis.legend()

band_figure.tight_layout()

band_output_path = Path(
    "../figures/"
    "p18_waveguide_projected_band_structure.png"
)

band_figure.savefig(
    band_output_path,
    dpi=300,
    bbox_inches="tight"
)


# ============================================================================
# Plot selected guided-mode intensity
# ============================================================================

field_figure, field_axis = (
    plt.subplots(
        figsize=(7.0, 5.0)
    )
)

field_image = field_axis.imshow(
    selected_intensity,
    extent=[
        x_grid[0] / lattice_constant,
        x_grid[-1] / lattice_constant,
        y_grid[0] / lattice_constant,
        y_grid[-1] / lattice_constant
    ],
    origin="lower",
    aspect="equal"
)

field_axis.axhline(
    -0.5,
    linewidth=0.8,
    alpha=0.60
)

field_axis.axhline(
    0.5,
    linewidth=0.8,
    alpha=0.60
)

field_axis.set_xlabel(
    r"$x/a$"
)

field_axis.set_ylabel(
    r"$y/a$"
)

field_axis.set_title(
    rf"TM Guided Mode, "
    rf"$\omega a/2\pi c="
    rf"{selected_frequency:.3f}$"
)

field_figure.colorbar(
    field_image,
    ax=field_axis,
    label=r"Normalized $|E_z|^2$"
)

field_figure.tight_layout()

field_output_path = Path(
    "../figures/"
    "p18_waveguide_mode_field.png"
)

field_figure.savefig(
    field_output_path,
    dpi=300,
    bbox_inches="tight"
)


# ============================================================================
# Display figures
# ============================================================================

plt.show()
