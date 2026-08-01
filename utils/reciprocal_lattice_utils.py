"""
Utility functions for reciprocal-lattice calculations,
Brillouin-zone construction, and high-symmetry points.
"""

import numpy as np


def reciprocal_lattice_vectors(
    real_vector_1: np.ndarray,
    real_vector_2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the reciprocal-lattice basis vectors corresponding
    to two two-dimensional real-space lattice vectors.

    The reciprocal vectors satisfy

        a_i · b_j = 2π δ_ij.

    Parameters
    ----------
    real_vector_1:
        First real-space lattice vector.

    real_vector_2:
        Second real-space lattice vector.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The reciprocal-lattice vectors b_1 and b_2.
    """

    real_lattice_matrix = np.column_stack(
        (
            real_vector_1,
            real_vector_2
        )
    )

    reciprocal_lattice_matrix = (
        2.0
        * np.pi
        * np.linalg.inv(real_lattice_matrix).T
    )

    reciprocal_vector_1 = reciprocal_lattice_matrix[:, 0]
    reciprocal_vector_2 = reciprocal_lattice_matrix[:, 1]

    return reciprocal_vector_1, reciprocal_vector_2


def generate_reciprocal_lattice_points(
    reciprocal_vector_1: np.ndarray,
    reciprocal_vector_2: np.ndarray,
    index_limit: int
) -> np.ndarray:
    """
    Generate reciprocal-lattice points of the form

        G = m b_1 + n b_2,

    where m and n are integers between -index_limit and index_limit.

    Parameters
    ----------
    reciprocal_vector_1:
        First reciprocal-lattice basis vector.

    reciprocal_vector_2:
        Second reciprocal-lattice basis vector.

    index_limit:
        Maximum absolute value of the reciprocal-lattice indices.

    Returns
    -------
    np.ndarray
        Array containing the generated reciprocal-lattice points.
    """

    if index_limit < 0:
        raise ValueError(
            "index_limit must be non-negative."
        )

    reciprocal_points = []

    for first_index in range(-index_limit, index_limit + 1):
        for second_index in range(-index_limit, index_limit + 1):
            reciprocal_point = (
                first_index * reciprocal_vector_1
                + second_index * reciprocal_vector_2
            )

            reciprocal_points.append(
                reciprocal_point
            )

    return np.asarray(
        reciprocal_points
    )


def square_brillouin_zone_vertices(
    lattice_constant: float
) -> np.ndarray:
    """
    Return the four vertices of the first Brillouin zone
    of a two-dimensional square lattice.

    The first Brillouin zone is bounded by

        -π/a <= k_x <= π/a,
        -π/a <= k_y <= π/a.

    Parameters
    ----------
    lattice_constant:
        Real-space lattice constant.

    Returns
    -------
    np.ndarray
        Vertices of the square first Brillouin zone.
    """

    if lattice_constant <= 0.0:
        raise ValueError(
            "lattice_constant must be positive."
        )

    zone_boundary = np.pi / lattice_constant

    return np.array(
        [
            [-zone_boundary, -zone_boundary],
            [zone_boundary, -zone_boundary],
            [zone_boundary, zone_boundary],
            [-zone_boundary, zone_boundary]
        ]
    )


def square_lattice_high_symmetry_points(
    lattice_constant: float
) -> dict[str, np.ndarray]:
    """
    Return the standard high-symmetry points of the first
    Brillouin zone of a square lattice.

    The points are

        Gamma = (0, 0),
        X     = (π/a, 0),
        M     = (π/a, π/a).

    Parameters
    ----------
    lattice_constant:
        Real-space lattice constant.

    Returns
    -------
    dict[str, np.ndarray]
        Dictionary containing the Gamma, X, and M points.
    """

    if lattice_constant <= 0.0:
        raise ValueError(
            "lattice_constant must be positive."
        )

    zone_boundary = np.pi / lattice_constant

    return {
        "Gamma": np.array(
            [
                0.0,
                0.0
            ]
        ),
        "X": np.array(
            [
                zone_boundary,
                0.0
            ]
        ),
        "M": np.array(
            [
                zone_boundary,
                zone_boundary
            ]
        )
    }