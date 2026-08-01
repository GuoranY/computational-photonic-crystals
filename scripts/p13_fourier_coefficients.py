"""
Project:
    Computational Photonic Crystals

Module:
    P13 - Fourier Coefficients of the Dielectric Functions

Description:
    Calculate and visualize the reciprocal-space Fourier coefficients
    of the dielectric function and inverse dielectric function for a
    two-dimensional square lattice of circular dielectric rods.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(
        0,
        str(PROJECT_ROOT)
    )


from utils.fourier_utils import (
    dielectric_fourier_coefficients,
    inverse_dielectric_fourier_coefficients
)

from utils.reciprocal_lattice_utils import (
    reciprocal_lattice_vectors,
    generate_reciprocal_lattice_points
)

from utils.plotting_utils import (
    plot_fourier_coefficient_maps
)


# ============================================================================
# Main calculation
# ============================================================================

def main() -> None:
    """
    Calculate and visualize the Fourier coefficients of the dielectric
    function and inverse dielectric function for a square lattice of
    circular dielectric rods.
    """

    # ------------------------------------------------------------------------
    # Physical parameters
    # ------------------------------------------------------------------------

    lattice_constant = 1.0

    background_refractive_index = 1.0
    rod_refractive_index = 3.5

    rod_radius = 0.2 * lattice_constant

    reciprocal_index_limit = 5

    # ------------------------------------------------------------------------
    # Construct real-space lattice vectors
    # ------------------------------------------------------------------------

    real_vector_1 = np.array(
        [
            lattice_constant,
            0.0
        ]
    )

    real_vector_2 = np.array(
        [
            0.0,
            lattice_constant
        ]
    )

    # ------------------------------------------------------------------------
    # Calculate reciprocal-lattice basis vectors
    # ------------------------------------------------------------------------

    (
        reciprocal_vector_1,
        reciprocal_vector_2
    ) = reciprocal_lattice_vectors(
        real_vector_1=real_vector_1,
        real_vector_2=real_vector_2
    )

    # ------------------------------------------------------------------------
    # Generate reciprocal-lattice vectors
    # ------------------------------------------------------------------------

    reciprocal_points = generate_reciprocal_lattice_points(
        reciprocal_vector_1=reciprocal_vector_1,
        reciprocal_vector_2=reciprocal_vector_2,
        index_limit=reciprocal_index_limit
    )

    number_of_indices = (
        2 * reciprocal_index_limit + 1
    )

    reciprocal_vector_grid = reciprocal_points.reshape(
        number_of_indices,
        number_of_indices,
        2
    )

    # ------------------------------------------------------------------------
    # Calculate Fourier coefficients
    # ------------------------------------------------------------------------

    dielectric_coefficients = dielectric_fourier_coefficients(
        reciprocal_vectors=reciprocal_vector_grid,
        lattice_constant=lattice_constant,
        radius=rod_radius,
        rod_refractive_index=rod_refractive_index,
        background_refractive_index=background_refractive_index
    )

    inverse_dielectric_coefficients = (
        inverse_dielectric_fourier_coefficients(
            reciprocal_vectors=reciprocal_vector_grid,
            lattice_constant=lattice_constant,
            radius=rod_radius,
            rod_refractive_index=rod_refractive_index,
            background_refractive_index=background_refractive_index
        )
    )

    # ------------------------------------------------------------------------
    # Calculate reference quantities
    # ------------------------------------------------------------------------

    filling_fraction = (
        np.pi
        * rod_radius ** 2
        / lattice_constant ** 2
    )

    zero_index = reciprocal_index_limit

    zero_order_dielectric_coefficient = (
        dielectric_coefficients[
            zero_index,
            zero_index
        ]
    )

    zero_order_inverse_dielectric_coefficient = (
        inverse_dielectric_coefficients[
            zero_index,
            zero_index
        ]
    )

    # ------------------------------------------------------------------------
    # Print calculated values
    # ------------------------------------------------------------------------

    print(
        "First reciprocal-lattice vector:"
    )

    print(
        np.array2string(
            reciprocal_vector_1,
            precision=8
        )
    )

    print(
        "\nSecond reciprocal-lattice vector:"
    )

    print(
        np.array2string(
            reciprocal_vector_2,
            precision=8
        )
    )

    print(
        f"\nFilling fraction: "
        f"{filling_fraction:.6f}"
    )

    print(
        "Zero-order dielectric coefficient: "
        f"{zero_order_dielectric_coefficient:.6f}"
    )

    print(
        "Zero-order inverse-dielectric coefficient: "
        f"{zero_order_inverse_dielectric_coefficient:.6f}"
    )

    # ------------------------------------------------------------------------
    # Plot Fourier-coefficient maps
    # ------------------------------------------------------------------------

    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(12.5, 5.2),
        constrained_layout=True
    )

    plot_fourier_coefficient_maps(
        axes=axes,
        dielectric_coefficients=dielectric_coefficients,
        inverse_dielectric_coefficients=(
            inverse_dielectric_coefficients
        ),
        maximum_index=reciprocal_index_limit
    )

    figure.suptitle(
        "Fourier Coefficients of the Dielectric Functions",
        fontsize=14
    )

    # ------------------------------------------------------------------------
    # Save and display figure
    # ------------------------------------------------------------------------

    output_path = (
        PROJECT_ROOT
        / "figures"
        / "p13_fourier_coefficients.png"
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


# ============================================================================
# Script entry point
# ============================================================================

if __name__ == "__main__":
    main()
