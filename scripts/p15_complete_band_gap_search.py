"""
Project:
    Computational Photonic Crystals

Module:
    P15 - Complete Photonic Band-Gap Search

Description:
    Calculate the TE and TM photonic band structures of a
    two-dimensional square lattice and identify frequency
    intervals that are forbidden for both polarizations.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utils.band_gap_utils import (
    find_complete_band_gaps,
    print_band_gaps,
    print_complete_band_gaps,
)
from utils.plane_wave_expansion_utils import (
    calculate_te_tm_band_structure,
    dielectric_convolution_matrix,
    reciprocal_lattice_vector_set,
    square_lattice_wavevector_path,
)
from utils.reciprocal_lattice_utils import (
    reciprocal_lattice_vectors,
)


# ============================================================================
# Structural parameters
# ============================================================================

background_refractive_index = 1.0
rod_refractive_index = 3.5

background_permittivity = (
    background_refractive_index**2
)

rod_permittivity = (
    rod_refractive_index**2
)

lattice_constant = 1.0

rod_radius = (
    0.2 * lattice_constant
)


# ============================================================================
# Plane-wave expansion parameters
# ============================================================================

reciprocal_index_limit = 3

points_per_segment = 30

number_of_bands = 8


# ============================================================================
# Real-space and reciprocal-space basis vectors
# ============================================================================

real_vector_1 = np.array(
    [
        lattice_constant,
        0.0,
    ]
)

real_vector_2 = np.array(
    [
        0.0,
        lattice_constant,
    ]
)

(
    reciprocal_vector_1,
    reciprocal_vector_2,
) = reciprocal_lattice_vectors(
    real_vector_1=real_vector_1,
    real_vector_2=real_vector_2,
)


# ============================================================================
# Output directory
# ============================================================================

project_root = Path(
    __file__
).resolve().parents[1]

figure_directory = (
    project_root / "figures"
)

figure_directory.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================================
# Main calculation
# ============================================================================

if __name__ == "__main__":
    # ------------------------------------------------------------------------
    # Generate the truncated reciprocal-lattice basis
    # ------------------------------------------------------------------------

    (
        reciprocal_vectors,
        reciprocal_indices,
    ) = reciprocal_lattice_vector_set(
        reciprocal_vector_1=reciprocal_vector_1,
        reciprocal_vector_2=reciprocal_vector_2,
        index_limit=reciprocal_index_limit,
    )

    number_of_plane_waves = len(
        reciprocal_vectors
    )

    # ------------------------------------------------------------------------
    # Construct the Gamma-X-M-Gamma wavevector path
    # ------------------------------------------------------------------------

    (
        wavevectors,
        symmetry_positions,
        symmetry_labels,
    ) = square_lattice_wavevector_path(
        lattice_constant=lattice_constant,
        points_per_segment=points_per_segment,
    )

    # ------------------------------------------------------------------------
    # Construct material convolution matrices
    # ------------------------------------------------------------------------

    dielectric_matrix = (
        dielectric_convolution_matrix(
            reciprocal_indices=reciprocal_indices,
            lattice_constant=lattice_constant,
            rod_radius=rod_radius,
            rod_permittivity=rod_permittivity,
            background_permittivity=(
                background_permittivity
            ),
            inverse=False,
        )
    )

    inverse_dielectric_matrix = (
        dielectric_convolution_matrix(
            reciprocal_indices=reciprocal_indices,
            lattice_constant=lattice_constant,
            rod_radius=rod_radius,
            rod_permittivity=rod_permittivity,
            background_permittivity=(
                background_permittivity
            ),
            inverse=True,
        )
    )

    # ------------------------------------------------------------------------
    # Calculate TE and TM band structures
    # ------------------------------------------------------------------------

    te_bands, tm_bands = (
        calculate_te_tm_band_structure(
            wavevectors=wavevectors,
            reciprocal_vectors=reciprocal_vectors,
            dielectric_matrix=dielectric_matrix,
            inverse_dielectric_matrix=(
                inverse_dielectric_matrix
            ),
            lattice_constant=lattice_constant,
            number_of_bands=number_of_bands,
        )
    )

    # ------------------------------------------------------------------------
    # Search for TE, TM, and complete band gaps
    # ------------------------------------------------------------------------

    (
        te_band_gaps,
        tm_band_gaps,
        complete_band_gaps,
    ) = find_complete_band_gaps(
        te_frequencies=te_bands,
        tm_frequencies=tm_bands,
    )

    # ------------------------------------------------------------------------
    # Print calculation information
    # ------------------------------------------------------------------------

    print(
        "Complete photonic band-gap search"
    )

    print(
        f"Reciprocal index limit: "
        f"{reciprocal_index_limit}"
    )

    print(
        f"Number of plane waves: "
        f"{number_of_plane_waves}"
    )

    print(
        f"Number of wavevectors: "
        f"{len(wavevectors)}"
    )

    print(
        f"Number of calculated bands: "
        f"{number_of_bands}"
    )

    print_band_gaps(
        band_gaps=te_band_gaps,
        polarization="TE",
    )

    print_band_gaps(
        band_gaps=tm_band_gaps,
        polarization="TM",
    )

    print_complete_band_gaps(
        complete_band_gaps=complete_band_gaps
    )

    # ------------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------------

    path_coordinates = np.arange(
        len(wavevectors)
    )

    figure, axis = plt.subplots(
        figsize=(
            8.5,
            5.6,
        )
    )

    for band_index in range(
        number_of_bands
    ):
        axis.plot(
            path_coordinates,
            te_bands[:, band_index],
            linewidth=1.4,
            linestyle="--",
            label=(
                "TE bands"
                if band_index == 0
                else None
            ),
        )

        axis.plot(
            path_coordinates,
            tm_bands[:, band_index],
            linewidth=1.4,
            linestyle="-",
            label=(
                "TM bands"
                if band_index == 0
                else None
            ),
        )

    for symmetry_position in (
        symmetry_positions
    ):
        axis.axvline(
            symmetry_position,
            linewidth=0.8,
            alpha=0.45,
        )

    for gap_index, gap in enumerate(
        complete_band_gaps
    ):
        axis.axhspan(
            gap["lower_frequency"],
            gap["upper_frequency"],
            alpha=0.22,
            label=(
                "Complete band gap"
                if gap_index == 0
                else None
            ),
        )

    if not complete_band_gaps:
        axis.text(
            0.5,
            0.96,
            "No complete TE-TM band gap found",
            transform=axis.transAxes,
            horizontalalignment="center",
            verticalalignment="top",
        )

    axis.set_xlim(
        symmetry_positions[0],
        symmetry_positions[-1],
    )

    axis.set_xticks(
        symmetry_positions
    )

    axis.set_xticklabels(
        symmetry_labels
    )

    axis.set_xlabel(
        "Wavevector path"
    )

    axis.set_ylabel(
        r"Normalized frequency "
        r"$\omega a/(2\pi c)$"
    )

    axis.set_title(
        "Complete Photonic Band-Gap Search"
    )

    axis.grid(
        alpha=0.22
    )

    axis.legend()

    figure.tight_layout()

    output_path = (
        figure_directory
        / "p15_complete_band_gap_search.png"
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
