"""
Project:
    Computational Photonic Crystals

Module:
    P16 - Two-Dimensional Point-Defect Cavity

Description:
    Construct a square photonic-crystal supercell containing a
    missing dielectric rod, calculate the TM supercell band
    structure, identify a defect mode inside the original
    photonic band gap, and reconstruct its localized electric-field
    intensity.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from utils.band_gap_utils import (
    find_band_gaps,
)
from utils.plane_wave_expansion_utils import (
    dielectric_convolution_matrix,
    reciprocal_lattice_vector_set,
    square_lattice_wavevector_path,
    tm_eigenfrequencies,
)
from utils.reciprocal_lattice_utils import (
    reciprocal_lattice_vectors,
)
from utils.supercell_utils import (
    calculate_tm_supercell_band_structure,
    field_localization_fraction,
    normalized_field_intensity,
    reconstruct_tm_field,
    square_supercell_rod_positions,
    square_supercell_wavevector_path,
    supercell_convolution_matrix,
    supercell_reciprocal_vector_set,
    tm_supercell_eigensystem,
)


# ============================================================================
# Primitive photonic-crystal parameters
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
# Primitive-cell plane-wave expansion parameters
# ============================================================================

primitive_reciprocal_index_limit = 3

primitive_points_per_segment = 30

primitive_number_of_bands = 8


# ============================================================================
# Defect-supercell parameters
# ============================================================================

number_of_supercell_cells = 5

supercell_length = (
    number_of_supercell_cells
    * lattice_constant
)

supercell_reciprocal_index_limit = 7

supercell_points_per_segment = 10

supercell_number_of_bands = 32


# ============================================================================
# Field-reconstruction parameters
# ============================================================================

number_of_field_points = 181

localization_radius = (
    0.75 * lattice_constant
)


# ============================================================================
# Primitive real-space and reciprocal-space basis vectors
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
    # Calculate the TM band structure of the perfect primitive crystal
    # ------------------------------------------------------------------------

    (
        primitive_reciprocal_vectors,
        primitive_reciprocal_indices,
    ) = reciprocal_lattice_vector_set(
        reciprocal_vector_1=(
            reciprocal_vector_1
        ),
        reciprocal_vector_2=(
            reciprocal_vector_2
        ),
        index_limit=(
            primitive_reciprocal_index_limit
        ),
    )

    (
        primitive_wavevectors,
        _,
        _,
    ) = square_lattice_wavevector_path(
        lattice_constant=lattice_constant,
        points_per_segment=(
            primitive_points_per_segment
        ),
    )

    primitive_dielectric_matrix = (
        dielectric_convolution_matrix(
            reciprocal_indices=(
                primitive_reciprocal_indices
            ),
            lattice_constant=lattice_constant,
            rod_radius=rod_radius,
            rod_permittivity=(
                rod_permittivity
            ),
            background_permittivity=(
                background_permittivity
            ),
            inverse=False,
        )
    )

    primitive_tm_bands = np.zeros(
        (
            len(
                primitive_wavevectors
            ),
            primitive_number_of_bands,
        )
    )

    for (
        wavevector_index,
        wavevector,
    ) in enumerate(
        primitive_wavevectors
    ):
        primitive_tm_bands[
            wavevector_index
        ] = tm_eigenfrequencies(
            wavevector=wavevector,
            reciprocal_vectors=(
                primitive_reciprocal_vectors
            ),
            dielectric_matrix=(
                primitive_dielectric_matrix
            ),
            lattice_constant=(
                lattice_constant
            ),
            number_of_bands=(
                primitive_number_of_bands
            ),
        )

    primitive_tm_band_gaps = find_band_gaps(
        frequencies=primitive_tm_bands
    )

    if not primitive_tm_band_gaps:
        raise RuntimeError(
            "No TM band gap was found for the "
            "primitive photonic crystal."
        )

    target_band_gap = (
        primitive_tm_band_gaps[0]
    )

    lower_gap_frequency = (
        target_band_gap[
            "lower_frequency"
        ]
    )

    upper_gap_frequency = (
        target_band_gap[
            "upper_frequency"
        ]
    )

    midgap_frequency = (
        target_band_gap[
            "midgap_frequency"
        ]
    )

    # ------------------------------------------------------------------------
    # Construct the point-defect supercell
    # ------------------------------------------------------------------------

    perfect_rod_positions = (
        square_supercell_rod_positions(
            number_of_cells=(
                number_of_supercell_cells
            ),
            lattice_constant=(
                lattice_constant
            ),
            remove_center_rod=False,
        )
    )

    defect_rod_positions = (
        square_supercell_rod_positions(
            number_of_cells=(
                number_of_supercell_cells
            ),
            lattice_constant=(
                lattice_constant
            ),
            remove_center_rod=True,
        )
    )

    (
        supercell_reciprocal_vectors,
        supercell_reciprocal_indices,
    ) = supercell_reciprocal_vector_set(
        supercell_length=(
            supercell_length
        ),
        reciprocal_index_limit=(
            supercell_reciprocal_index_limit
        ),
    )

    defect_dielectric_matrix = (
        supercell_convolution_matrix(
            reciprocal_indices=(
                supercell_reciprocal_indices
            ),
            supercell_length=(
                supercell_length
            ),
            rod_radius=rod_radius,
            rod_positions=(
                defect_rod_positions
            ),
            rod_value=rod_permittivity,
            background_value=(
                background_permittivity
            ),
        )
    )

    # ------------------------------------------------------------------------
    # Calculate the defect-supercell band structure
    # ------------------------------------------------------------------------

    (
        supercell_wavevectors,
        supercell_symmetry_positions,
        supercell_symmetry_labels,
    ) = square_supercell_wavevector_path(
        supercell_length=(
            supercell_length
        ),
        points_per_segment=(
            supercell_points_per_segment
        ),
    )

    defect_tm_bands = (
        calculate_tm_supercell_band_structure(
            wavevectors=(
                supercell_wavevectors
            ),
            reciprocal_vectors=(
                supercell_reciprocal_vectors
            ),
            dielectric_matrix=(
                defect_dielectric_matrix
            ),
            normalization_length=(
                lattice_constant
            ),
            number_of_bands=(
                supercell_number_of_bands
            ),
        )
    )

    # ------------------------------------------------------------------------
    # Calculate Gamma-point supercell eigenmodes
    # ------------------------------------------------------------------------

    gamma_wavevector = np.array(
        [
            0.0,
            0.0,
        ]
    )

    (
        gamma_frequencies,
        gamma_eigenvectors,
    ) = tm_supercell_eigensystem(
        wavevector=gamma_wavevector,
        reciprocal_vectors=(
            supercell_reciprocal_vectors
        ),
        dielectric_matrix=(
            defect_dielectric_matrix
        ),
        normalization_length=(
            lattice_constant
        ),
        number_of_bands=(
            supercell_number_of_bands
        ),
    )

    candidate_mode_indices = np.where(
        (
            gamma_frequencies
            > lower_gap_frequency
        )
        & (
            gamma_frequencies
            < upper_gap_frequency
        )
    )[0]

    if len(
        candidate_mode_indices
    ) == 0:
        raise RuntimeError(
            "No Gamma-point defect mode was found "
            "inside the selected TM band gap. "
            "Increase supercell_reciprocal_index_limit "
            "or supercell_number_of_bands."
        )

    # ------------------------------------------------------------------------
    # Reconstruct every candidate and select the most localized mode
    # ------------------------------------------------------------------------

    half_supercell_length = (
        0.5 * supercell_length
    )

    x_coordinates = np.linspace(
        -half_supercell_length,
        half_supercell_length,
        number_of_field_points,
    )

    y_coordinates = np.linspace(
        -half_supercell_length,
        half_supercell_length,
        number_of_field_points,
    )

    candidate_results = []

    for mode_index in candidate_mode_indices:
        current_field = reconstruct_tm_field(
            x_coordinates=x_coordinates,
            y_coordinates=y_coordinates,
            wavevector=gamma_wavevector,
            reciprocal_vectors=(
                supercell_reciprocal_vectors
            ),
            eigenvector=(
                gamma_eigenvectors[
                    :,
                    mode_index
                ]
            ),
        )

        current_intensity = (
            normalized_field_intensity(
                field=current_field
            )
        )

        current_localization_fraction = (
            field_localization_fraction(
                intensity=current_intensity,
                x_coordinates=(
                    x_coordinates
                ),
                y_coordinates=(
                    y_coordinates
                ),
                localization_radius=(
                    localization_radius
                ),
            )
        )

        candidate_results.append(
            {
                "mode_index": int(
                    mode_index
                ),
                "frequency": float(
                    gamma_frequencies[
                        mode_index
                    ]
                ),
                "field": current_field,
                "intensity": (
                    current_intensity
                ),
                "localization_fraction": float(
                    current_localization_fraction
                ),
            }
        )

    selected_mode = max(
        candidate_results,
        key=lambda result: result[
            "localization_fraction"
        ],
    )

    selected_mode_index = (
        selected_mode[
            "mode_index"
        ]
    )

    selected_mode_frequency = (
        selected_mode[
            "frequency"
        ]
    )

    selected_intensity = (
        selected_mode[
            "intensity"
        ]
    )

    selected_localization_fraction = (
        selected_mode[
            "localization_fraction"
        ]
    )

    # ------------------------------------------------------------------------
    # Print calculation information
    # ------------------------------------------------------------------------

    print(
        "Two-dimensional point-defect cavity"
    )

    print(
        f"Supercell size: "
        f"{number_of_supercell_cells} x "
        f"{number_of_supercell_cells}"
    )

    print(
        f"Number of rods in perfect supercell: "
        f"{len(perfect_rod_positions)}"
    )

    print(
        f"Number of rods in defect supercell: "
        f"{len(defect_rod_positions)}"
    )

    print(
        f"Supercell reciprocal index limit: "
        f"{supercell_reciprocal_index_limit}"
    )

    print(
        f"Number of supercell plane waves: "
        f"{len(supercell_reciprocal_vectors)}"
    )

    print()
    print(
        "Selected primitive-crystal TM band gap"
    )

    print(
        f"Lower frequency: "
        f"{lower_gap_frequency:.6f}"
    )

    print(
        f"Upper frequency: "
        f"{upper_gap_frequency:.6f}"
    )

    print(
        f"Midgap frequency: "
        f"{midgap_frequency:.6f}"
    )

    print()
    print(
        "Gamma-point modes inside the TM band gap"
    )

    for result in candidate_results:
        print(
            f"Band {result['mode_index'] + 1}: "
            f"frequency = "
            f"{result['frequency']:.6f}, "
            f"localization fraction = "
            f"{100.0 * result['localization_fraction']:.3f}%"
        )

    print()
    print(
        "Selected defect-cavity mode"
    )

    print(
        f"Band index: "
        f"{selected_mode_index + 1}"
    )

    print(
        f"Normalized frequency: "
        f"{selected_mode_frequency:.6f}"
    )

    print(
        f"Localization fraction inside "
        f"r <= {localization_radius / lattice_constant:.2f}a: "
        f"{100.0 * selected_localization_fraction:.3f}%"
    )

    # ------------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------------

    path_coordinates = np.arange(
        len(
            supercell_wavevectors
        )
    )

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(
            16.0,
            5.2,
        ),
    )

    # ------------------------------------------------------------------------
    # Plot the defect-supercell geometry
    # ------------------------------------------------------------------------

    geometry_axis = axes[0]

    geometry_axis.set_facecolor(
        "white"
    )

    for rod_position in defect_rod_positions:
        rod = Circle(
            xy=(
                rod_position[0],
                rod_position[1],
            ),
            radius=rod_radius,
            fill=True,
            alpha=0.78,
        )

        geometry_axis.add_patch(
            rod
        )

    defect_marker = Circle(
        xy=(
            0.0,
            0.0,
        ),
        radius=rod_radius,
        fill=False,
        linestyle="--",
        linewidth=1.6,
        alpha=0.85,
    )

    geometry_axis.add_patch(
        defect_marker
    )

    geometry_axis.axhline(
        0.0,
        linewidth=0.7,
        alpha=0.25,
    )

    geometry_axis.axvline(
        0.0,
        linewidth=0.7,
        alpha=0.25,
    )

    geometry_axis.set_xlim(
        -half_supercell_length,
        half_supercell_length,
    )

    geometry_axis.set_ylim(
        -half_supercell_length,
        half_supercell_length,
    )

    geometry_axis.set_aspect(
        "equal"
    )

    geometry_axis.set_xlabel(
        r"$x/a$"
    )

    geometry_axis.set_ylabel(
        r"$y/a$"
    )

    geometry_axis.set_title(
        "Point-Defect Supercell"
    )

    geometry_axis.grid(
        alpha=0.18
    )

    # ------------------------------------------------------------------------
    # Plot the defect-supercell bands
    # ------------------------------------------------------------------------

    band_axis = axes[1]

    for band_index in range(
        supercell_number_of_bands
    ):
        is_selected_band = (
            band_index
            == selected_mode_index
        )

        band_axis.plot(
            path_coordinates,
            defect_tm_bands[
                :,
                band_index
            ],
            linewidth=(
                2.4
                if is_selected_band
                else 1.0
            ),
            alpha=(
                1.0
                if is_selected_band
                else 0.68
            ),
            label=(
                "Defect band"
                if is_selected_band
                else None
            ),
        )

    band_axis.axhspan(
        lower_gap_frequency,
        upper_gap_frequency,
        alpha=0.20,
        label=(
            "Perfect-crystal TM gap"
        ),
    )

    band_axis.axhline(
        selected_mode_frequency,
        linestyle="--",
        linewidth=1.1,
        alpha=0.75,
    )

    for symmetry_position in (
        supercell_symmetry_positions
    ):
        band_axis.axvline(
            symmetry_position,
            linewidth=0.8,
            alpha=0.40,
        )

    band_axis.set_xlim(
        supercell_symmetry_positions[0],
        supercell_symmetry_positions[-1],
    )

    vertical_margin = (
        0.18
        * (
            upper_gap_frequency
            - lower_gap_frequency
        )
    )

    band_axis.set_ylim(
        lower_gap_frequency
        - vertical_margin,
        upper_gap_frequency
        + vertical_margin,
    )

    band_axis.set_xticks(
        supercell_symmetry_positions
    )

    band_axis.set_xticklabels(
        supercell_symmetry_labels
    )

    band_axis.set_xlabel(
        "Supercell wavevector path"
    )

    band_axis.set_ylabel(
        r"Normalized frequency "
        r"$\omega a/(2\pi c)$"
    )

    band_axis.set_title(
        "TM Defect-Supercell Bands"
    )

    band_axis.grid(
        alpha=0.20
    )

    band_axis.legend(
        fontsize=8
    )

    # ------------------------------------------------------------------------
    # Plot the localized defect-mode intensity
    # ------------------------------------------------------------------------

    field_axis = axes[2]

    field_image = field_axis.imshow(
        selected_intensity,
        origin="lower",
        extent=[
            -half_supercell_length,
            half_supercell_length,
            -half_supercell_length,
            half_supercell_length,
        ],
        aspect="equal",
        interpolation="bilinear",
    )

    for rod_position in defect_rod_positions:
        rod_outline = Circle(
            xy=(
                rod_position[0],
                rod_position[1],
            ),
            radius=rod_radius,
            fill=False,
            linewidth=0.65,
            alpha=0.55,
        )

        field_axis.add_patch(
            rod_outline
        )

    localization_circle = Circle(
        xy=(
            0.0,
            0.0,
        ),
        radius=localization_radius,
        fill=False,
        linestyle="--",
        linewidth=1.1,
        alpha=0.75,
    )

    field_axis.add_patch(
        localization_circle
    )

    field_axis.set_xlim(
        -half_supercell_length,
        half_supercell_length,
    )

    field_axis.set_ylim(
        -half_supercell_length,
        half_supercell_length,
    )

    field_axis.set_xlabel(
        r"$x/a$"
    )

    field_axis.set_ylabel(
        r"$y/a$"
    )

    field_axis.set_title(
        (
            r"Localized TM Mode "
            f"\n"
            r"$\omega a/(2\pi c)="
            f"{selected_mode_frequency:.4f}"
            r"$"
        )
    )

    colorbar = figure.colorbar(
        field_image,
        ax=field_axis,
        fraction=0.046,
        pad=0.04,
    )

    colorbar.set_label(
        r"Normalized $|E_z|^2$"
    )

    figure.suptitle(
        "Two-Dimensional Point-Defect Cavity",
        fontsize=14,
    )

    figure.tight_layout()

    output_path = (
        figure_directory
        / "p16_2d_tm_defect_cavity.png"
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    print()
    print(
        f"Figure saved to: "
        f"{output_path}"
    )

    plt.show()
