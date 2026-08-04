"""
Project:
    Computational Photonic Crystals

Module:
    P17 - Two-Dimensional TE Point-Defect Cavity

Description:
    Construct a square photonic-crystal supercell containing a
    missing dielectric rod, calculate the TE defect-supercell band
    structure, search for modes inside the perfect-crystal TE band
    gaps, and visualize the most localized TE mode candidate.

    If no genuine in-gap TE defect mode exists, the program selects
    the Gamma-point TE mode with the strongest central localization
    and labels it as the best TE candidate rather than as a genuine
    localized defect mode.
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
    te_eigenfrequencies,
)
from utils.reciprocal_lattice_utils import (
    reciprocal_lattice_vectors,
)
from utils.supercell_utils import (
    calculate_te_supercell_band_structure,
    field_localization_fraction,
    normalized_field_intensity,
    reconstruct_tm_field,
    square_supercell_rod_positions,
    square_supercell_wavevector_path,
    supercell_convolution_matrix,
    supercell_reciprocal_vector_set,
    te_supercell_eigensystem,
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

background_inverse_permittivity = (
    1.0
    / background_permittivity
)

rod_inverse_permittivity = (
    1.0
    / rod_permittivity
)

lattice_constant = 1.0

rod_radius = (
    0.2
    * lattice_constant
)


# ============================================================================
# Primitive-cell plane-wave expansion parameters
# ============================================================================

primitive_reciprocal_index_limit = 5

primitive_points_per_segment = 40

primitive_number_of_bands = 10


# ============================================================================
# Defect-supercell parameters
# ============================================================================

number_of_supercell_cells = 5

supercell_length = (
    number_of_supercell_cells
    * lattice_constant
)

supercell_reciprocal_index_limit = 7

supercell_points_per_segment = 12

supercell_number_of_bands = 40


# ============================================================================
# Field-reconstruction and localization parameters
# ============================================================================

number_of_field_points = 181

localization_radius = (
    0.75
    * lattice_constant
)


# ============================================================================
# Primitive real-space and reciprocal-space lattice vectors
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
    project_root
    / "figures"
)

figure_directory.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================================
# Main calculation
# ============================================================================

if __name__ == "__main__":
    # ========================================================================
    # Perfect primitive-crystal TE band structure
    # ========================================================================

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
        lattice_constant=(
            lattice_constant
        ),
        points_per_segment=(
            primitive_points_per_segment
        ),
    )

    primitive_inverse_dielectric_matrix = (
        dielectric_convolution_matrix(
            reciprocal_indices=(
                primitive_reciprocal_indices
            ),
            lattice_constant=(
                lattice_constant
            ),
            rod_radius=(
                rod_radius
            ),
            rod_permittivity=(
                rod_permittivity
            ),
            background_permittivity=(
                background_permittivity
            ),
            inverse=True,
        )
    )

    primitive_te_bands = np.zeros(
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
        primitive_te_bands[
            wavevector_index
        ] = te_eigenfrequencies(
            wavevector=(
                wavevector
            ),
            reciprocal_vectors=(
                primitive_reciprocal_vectors
            ),
            inverse_dielectric_matrix=(
                primitive_inverse_dielectric_matrix
            ),
            lattice_constant=(
                lattice_constant
            ),
            number_of_bands=(
                primitive_number_of_bands
            ),
        )

    primitive_te_band_gaps = find_band_gaps(
        frequencies=(
            primitive_te_bands
        )
    )

    print(
        "Perfect-crystal TE band gaps"
    )

    if primitive_te_band_gaps:
        for (
            gap_index,
            band_gap,
        ) in enumerate(
            primitive_te_band_gaps,
            start=1,
        ):
            print(
                f"Gap {gap_index}: "
                f"{band_gap['lower_frequency']:.6f} "
                f"to "
                f"{band_gap['upper_frequency']:.6f}"
            )

    else:
        print(
            "No TE band gaps were found."
        )

    # ========================================================================
    # Construct the missing-rod supercell
    # ========================================================================

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

    defect_inverse_dielectric_matrix = (
        supercell_convolution_matrix(
            reciprocal_indices=(
                supercell_reciprocal_indices
            ),
            supercell_length=(
                supercell_length
            ),
            rod_radius=(
                rod_radius
            ),
            rod_positions=(
                defect_rod_positions
            ),
            rod_value=(
                rod_inverse_permittivity
            ),
            background_value=(
                background_inverse_permittivity
            ),
        )
    )

    # ========================================================================
    # Supercell high-symmetry wavevector path
    # ========================================================================

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

    defect_te_bands = (
        calculate_te_supercell_band_structure(
            wavevectors=(
                supercell_wavevectors
            ),
            reciprocal_vectors=(
                supercell_reciprocal_vectors
            ),
            inverse_dielectric_matrix=(
                defect_inverse_dielectric_matrix
            ),
            normalization_length=(
                lattice_constant
            ),
            number_of_bands=(
                supercell_number_of_bands
            ),
        )
    )

    # ========================================================================
    # Gamma-point TE eigensystem
    # ========================================================================

    gamma_wavevector = np.array(
        [
            0.0,
            0.0,
        ]
    )

    (
        gamma_frequencies,
        gamma_eigenvectors,
    ) = te_supercell_eigensystem(
        wavevector=(
            gamma_wavevector
        ),
        reciprocal_vectors=(
            supercell_reciprocal_vectors
        ),
        inverse_dielectric_matrix=(
            defect_inverse_dielectric_matrix
        ),
        normalization_length=(
            lattice_constant
        ),
        number_of_bands=(
            supercell_number_of_bands
        ),
    )

    # ========================================================================
    # Real-space field grid
    # ========================================================================

    half_supercell_length = (
        0.5
        * supercell_length
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

    # ========================================================================
    # Search for Gamma-point modes inside perfect-crystal TE gaps
    # ========================================================================

    candidate_results = []

    for (
        gap_index,
        band_gap,
    ) in enumerate(
        primitive_te_band_gaps
    ):
        lower_gap_frequency = (
            band_gap[
                "lower_frequency"
            ]
        )

        upper_gap_frequency = (
            band_gap[
                "upper_frequency"
            ]
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

        for mode_index in candidate_mode_indices:
            current_field = reconstruct_tm_field(
                x_coordinates=(
                    x_coordinates
                ),
                y_coordinates=(
                    y_coordinates
                ),
                wavevector=(
                    gamma_wavevector
                ),
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
                    field=(
                        current_field
                    )
                )
            )

            current_localization_fraction = (
                field_localization_fraction(
                    intensity=(
                        current_intensity
                    ),
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
                    "gap_index": int(
                        gap_index
                    ),
                    "mode_index": int(
                        mode_index
                    ),
                    "frequency": float(
                        gamma_frequencies[
                            mode_index
                        ]
                    ),
                    "field": (
                        current_field
                    ),
                    "intensity": (
                        current_intensity
                    ),
                    "localization_fraction": float(
                        current_localization_fraction
                    ),
                    "lower_gap_frequency": float(
                        lower_gap_frequency
                    ),
                    "upper_gap_frequency": float(
                        upper_gap_frequency
                    ),
                }
            )

    # ========================================================================
    # Select an in-gap mode or use the best fallback candidate
    # ========================================================================

    genuine_in_gap_mode_found = bool(
        candidate_results
    )

    if genuine_in_gap_mode_found:
        selected_mode = max(
            candidate_results,
            key=lambda result: result[
                "localization_fraction"
            ],
        )

        selected_gap_index = (
            selected_mode[
                "gap_index"
            ]
        )

        selected_lower_gap_frequency = (
            selected_mode[
                "lower_gap_frequency"
            ]
        )

        selected_upper_gap_frequency = (
            selected_mode[
                "upper_gap_frequency"
            ]
        )

    else:
        print()
        print(
            "No Gamma-point TE supercell mode was found "
            "inside any perfect-crystal TE band gap."
        )

        print(
            "Searching all calculated Gamma-point TE modes "
            "for the strongest central localization."
        )

        fallback_results = []

        # Skip the first zero-frequency Gamma-point mode.
        for mode_index in range(
            1,
            supercell_number_of_bands,
        ):
            current_field = reconstruct_tm_field(
                x_coordinates=(
                    x_coordinates
                ),
                y_coordinates=(
                    y_coordinates
                ),
                wavevector=(
                    gamma_wavevector
                ),
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
                    field=(
                        current_field
                    )
                )
            )

            current_localization_fraction = (
                field_localization_fraction(
                    intensity=(
                        current_intensity
                    ),
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

            fallback_results.append(
                {
                    "mode_index": int(
                        mode_index
                    ),
                    "frequency": float(
                        gamma_frequencies[
                            mode_index
                        ]
                    ),
                    "field": (
                        current_field
                    ),
                    "intensity": (
                        current_intensity
                    ),
                    "localization_fraction": float(
                        current_localization_fraction
                    ),
                }
            )

        selected_mode = max(
            fallback_results,
            key=lambda result: result[
                "localization_fraction"
            ],
        )

        selected_gap_index = None

        selected_lower_gap_frequency = None

        selected_upper_gap_frequency = None

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

    # ========================================================================
    # Print selected-mode information
    # ========================================================================

    print()
    print(
        "Two-dimensional TE point-defect cavity"
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
        f"Number of supercell plane waves: "
        f"{len(supercell_reciprocal_vectors)}"
    )

    print()
    print(
        "Selected TE mode"
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
        f"r <= "
        f"{localization_radius / lattice_constant:.2f}a: "
        f"{selected_localization_fraction:.4f}"
    )

    if genuine_in_gap_mode_found:
        print(
            f"Perfect-crystal TE gap: "
            f"{selected_gap_index + 1}"
        )

        print(
            f"Gap interval: "
            f"{selected_lower_gap_frequency:.6f} "
            f"to "
            f"{selected_upper_gap_frequency:.6f}"
        )

        print(
            "Classification: genuine in-gap "
            "TE defect-mode candidate"
        )

    else:
        print(
            "Classification: best localized TE candidate, "
            "but not an in-gap defect mode"
        )

    # ========================================================================
    # Visualization
    # ========================================================================

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

    # ========================================================================
    # Plot the point-defect supercell
    # ========================================================================

    geometry_axis = axes[0]

    for rod_position in defect_rod_positions:
        rod = Circle(
            xy=(
                rod_position[0],
                rod_position[1],
            ),
            radius=(
                rod_radius
            ),
            alpha=0.78,
        )

        geometry_axis.add_patch(
            rod
        )

    missing_rod_outline = Circle(
        xy=(
            0.0,
            0.0,
        ),
        radius=(
            rod_radius
        ),
        fill=False,
        linestyle="--",
        linewidth=1.5,
        alpha=0.85,
    )

    geometry_axis.add_patch(
        missing_rod_outline
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
        alpha=0.22
    )

    # ========================================================================
    # Plot the TE defect-supercell bands
    # ========================================================================

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
            defect_te_bands[
                :,
                band_index
            ],
            linewidth=(
                2.4
                if is_selected_band
                else 0.9
            ),
            alpha=(
                1.0
                if is_selected_band
                else 0.68
            ),
            label=(
                "Selected candidate band"
                if is_selected_band
                else None
            ),
        )

    if genuine_in_gap_mode_found:
        band_axis.axhspan(
            selected_lower_gap_frequency,
            selected_upper_gap_frequency,
            alpha=0.16,
            label="Perfect-crystal TE gap",
        )

        band_axis.axhline(
            selected_mode_frequency,
            linestyle="--",
            linewidth=1.1,
            alpha=0.80,
            label="Selected-mode frequency",
        )

    for symmetry_position in (
            supercell_symmetry_positions
    ):
        band_axis.axvline(
            symmetry_position,
            linewidth=0.8,
            alpha=0.35,
        )

    band_axis.set_xlim(
        supercell_symmetry_positions[0],
        supercell_symmetry_positions[-1],
    )

    if genuine_in_gap_mode_found:
        selected_gap_width = (
                selected_upper_gap_frequency
                - selected_lower_gap_frequency
        )

        vertical_margin = max(
            0.18
            * selected_gap_width,
            0.01,
        )

        band_axis.set_ylim(
            selected_lower_gap_frequency
            - vertical_margin,
            selected_upper_gap_frequency
            + vertical_margin,
        )

    else:
        frequency_margin = 0.08

        band_axis.set_ylim(
            max(
                0.0,
                selected_mode_frequency
                - frequency_margin,
            ),
            selected_mode_frequency
            + frequency_margin,
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
        "TE Defect-Supercell Bands"
    )

    band_axis.grid(
        alpha=0.20
    )

    band_axis.legend(
        fontsize=8
    )

    # ========================================================================
    # Plot the selected TE field intensity
    # ========================================================================

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
        vmin=0.0,
        vmax=1.0,
    )

    for rod_position in defect_rod_positions:
        rod_outline = Circle(
            xy=(
                rod_position[0],
                rod_position[1],
            ),
            radius=(
                rod_radius
            ),
            fill=False,
            linewidth=0.65,
            alpha=0.50,
        )

        field_axis.add_patch(
            rod_outline
        )

    localization_circle = Circle(
        xy=(
            0.0,
            0.0,
        ),
        radius=(
            localization_radius
        ),
        fill=False,
        linestyle="--",
        linewidth=1.1,
        alpha=0.80,
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

    if genuine_in_gap_mode_found:
        field_title = (
            "Localized TE Mode\n"
            r"$\omega a/(2\pi c)="
            f"{selected_mode_frequency:.4f}"
            r"$, "
            r"\eta="
            f"{selected_localization_fraction:.2f}"
            r"$"
        )

    else:
        field_title = (
            "TE Best Candidate\n"
            r"$\omega a/(2\pi c)="
            f"{selected_mode_frequency:.4f}"
            r"$, "
            r"\eta="
            f"{selected_localization_fraction:.2f}"
            r"$"
        )

    field_axis.set_title(
        field_title
    )

    colorbar = figure.colorbar(
        field_image,
        ax=field_axis,
        fraction=0.046,
        pad=0.04,
    )

    colorbar.set_label(
        r"Normalized $|H_z|^2$"
    )

    figure.suptitle(
        "Two-Dimensional TE Point-Defect Cavity",
        fontsize=14,
    )

    figure.tight_layout()

    output_path = (
        figure_directory
        / "p17_2d_te_defect_cavity.png"
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
