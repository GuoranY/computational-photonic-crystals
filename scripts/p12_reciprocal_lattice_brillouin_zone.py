"""
Project:
    Computational Photonic Crystals

Module:
    P12 - Reciprocal Lattice and Brillouin Zone

Description:
    Construct the reciprocal lattice of a two-dimensional square lattice,
    visualize the first Brillouin zone, identify the high-symmetry points,
    and display the standard Gamma-X-M-Gamma path used for photonic
    band-structure calculations.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from utils.plotting_utils import (
    plot_first_brillouin_zone,
    plot_reciprocal_lattice
)
from utils.reciprocal_lattice_utils import (
    generate_reciprocal_lattice_points,
    reciprocal_lattice_vectors,
    square_brillouin_zone_vertices,
    square_lattice_high_symmetry_points
)


# ============================================================================
# Project paths
# ============================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

FIGURE_PATH = (
    PROJECT_ROOT
    / "figures"
    / "p12_reciprocal_lattice_brillouin_zone.png"
)


# ============================================================================
# Main program
# ============================================================================

def main() -> None:
    """
    Construct and visualize the reciprocal lattice and first
    Brillouin zone of a two-dimensional square lattice.
    """

    lattice_constant = 1.0
    reciprocal_index_limit = 2

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

    reciprocal_vector_1, reciprocal_vector_2 = (
        reciprocal_lattice_vectors(
            real_vector_1=real_vector_1,
            real_vector_2=real_vector_2
        )
    )

    reciprocal_points = generate_reciprocal_lattice_points(
        reciprocal_vector_1=reciprocal_vector_1,
        reciprocal_vector_2=reciprocal_vector_2,
        index_limit=reciprocal_index_limit
    )

    brillouin_zone_vertices = square_brillouin_zone_vertices(
        lattice_constant=lattice_constant
    )

    high_symmetry_points = square_lattice_high_symmetry_points(
        lattice_constant=lattice_constant
    )

    figure: Figure
    axes: np.ndarray

    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(12.0, 5.4)
    )

    plot_reciprocal_lattice(
        axis=axes[0],
        reciprocal_points=reciprocal_points,
        reciprocal_vector_1=reciprocal_vector_1,
        reciprocal_vector_2=reciprocal_vector_2,
        lattice_constant=lattice_constant
    )

    plot_first_brillouin_zone(
        axis=axes[1],
        brillouin_zone_vertices=brillouin_zone_vertices,
        high_symmetry_points=high_symmetry_points,
        lattice_constant=lattice_constant
    )

    figure.suptitle(
        "Reciprocal Lattice and Brillouin Zone of a Square Lattice",
        fontsize=14
    )

    figure.tight_layout()

    FIGURE_PATH.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    figure.savefig(
        FIGURE_PATH,
        dpi=300,
        bbox_inches="tight"
    )

    print(
        "Real-space lattice vectors:"
    )

    print(
        f"a1 = {real_vector_1}"
    )

    print(
        f"a2 = {real_vector_2}"
    )

    print(
        "\nReciprocal-lattice vectors:"
    )

    print(
        f"b1 = {reciprocal_vector_1}"
    )

    print(
        f"b2 = {reciprocal_vector_2}"
    )

    print(
        "\nHigh-symmetry points:"
    )

    for point_name, point_coordinates in high_symmetry_points.items():
        print(
            f"{point_name} = {point_coordinates}"
        )

    print(
        f"\nFigure saved to:\n{FIGURE_PATH}"
    )

    plt.show()


if __name__ == "__main__":
    main()
