"""
Utility functions for constructing
two-dimensional photonic-crystal waveguides.
"""

import numpy as np
from scipy.special import j1
from scipy.linalg import eigh


def square_lattice_waveguide_positions(
    lattice_constant: float,
    number_of_rows: int,
    number_of_columns: int
) -> np.ndarray:
    """
    Generate dielectric rod positions for a line-defect
    photonic crystal waveguide.

    The central row of rods is removed to create
    a waveguide defect.

    Parameters
    ----------
    lattice_constant:
        Lattice constant a.

    number_of_rows:
        Number of rows perpendicular to waveguide direction.

    number_of_columns:
        Number of rods along waveguide direction.

    Returns
    -------
    positions:
        Array containing rod coordinates.
    """

    positions = []

    center_row = number_of_rows // 2

    for iy in range(number_of_rows):

        # remove central row
        if iy == center_row:
            continue

        y = (
            iy - center_row
        ) * lattice_constant

        for ix in range(number_of_columns):
            x = (
                        ix - number_of_columns // 2
                ) * lattice_constant

            positions.append(
                [
                    x,
                    y
                ]
            )


    return np.array(positions)


def waveguide_supercell_vectors(
    lattice_constant: float,
    number_of_rows: int
):
    """
    Return real-space lattice vectors
    of the rectangular waveguide supercell.
    """

    a1 = np.array(
        [
            lattice_constant,
            0.0
        ]
    )

    a2 = np.array(
        [
            0.0,
            number_of_rows * lattice_constant
        ]
    )

    return a1, a2


def waveguide_reciprocal_vectors(
    lattice_constant: float,
    number_of_rows: int
):
    """
    Calculate reciprocal lattice vectors
    for the rectangular waveguide supercell.
    """

    b1 = np.array(
        [
            2*np.pi/lattice_constant,
            0.0
        ]
    )


    b2 = np.array(
        [
            0.0,
            2*np.pi /
            (
                number_of_rows*lattice_constant
            )
        ]
    )


    return b1, b2


def waveguide_supercell_rod_positions(
    lattice_constant: float,
    number_of_rows: int
) -> np.ndarray:
    """
    Generate rod positions inside one rectangular waveguide
    supercell.

    The supercell has width a along the waveguide direction and
    height number_of_rows * a perpendicular to the waveguide.
    The central rod is removed to form the line defect.

    Parameters
    ----------
    lattice_constant:
        Lattice constant a.

    number_of_rows:
        Number of primitive-cell rows contained in the
        waveguide supercell. This value must be odd.

    Returns
    -------
    np.ndarray
        Rod-center coordinates with shape
        (number_of_rows - 1, 2).
    """

    if number_of_rows < 3:
        raise ValueError(
            "number_of_rows must be at least 3."
        )

    if number_of_rows % 2 == 0:
        raise ValueError(
            "number_of_rows must be odd so that a unique "
            "central row can be removed."
        )

    center_row = number_of_rows // 2

    positions = []

    for row_index in range(number_of_rows):

        if row_index == center_row:
            continue

        y = (
            row_index - center_row
        ) * lattice_constant

        positions.append(
            [
                0.0,
                y
            ]
        )

    return np.asarray(
        positions,
        dtype=float
    )


def waveguide_wavevector_path(
    lattice_constant: float,
    number_of_points: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Construct the Gamma-to-X wavevector path along the
    waveguide direction.

    The line defect extends along x, so ky remains zero and
    kx varies from 0 to pi / a.

    Parameters
    ----------
    lattice_constant:
        Lattice constant a.

    number_of_points:
        Number of wavevectors along the path.

    Returns
    -------
    wavevectors:
        Array of shape (number_of_points, 2).

    normalized_wavevectors:
        Dimensionless coordinate k_x a / pi, ranging from
        zero at Gamma to one at X.
    """

    if number_of_points < 2:
        raise ValueError(
            "number_of_points must be at least 2."
        )

    kx_values = np.linspace(
        0.0,
        np.pi / lattice_constant,
        number_of_points
    )

    wavevectors = np.column_stack(
        (
            kx_values,
            np.zeros_like(kx_values)
        )
    )

    normalized_wavevectors = (
        kx_values * lattice_constant / np.pi
    )

    return wavevectors, normalized_wavevectors


def circular_rod_supercell_fourier_coefficient(
    reciprocal_vector: np.ndarray,
    rod_positions: np.ndarray,
    rod_radius: float,
    rod_value: float,
    background_value: float,
    supercell_area: float
) -> complex:
    """
    Calculate one Fourier coefficient of a piecewise-constant
    material function for identical circular rods at arbitrary
    positions inside a supercell.

    The material function equals rod_value inside each rod and
    background_value outside the rods.

    Parameters
    ----------
    reciprocal_vector:
        Reciprocal vector G.

    rod_positions:
        Array containing the center of every rod inside the
        supercell.

    rod_radius:
        Radius of each circular rod.

    rod_value:
        Material value inside each rod.

    background_value:
        Material value in the background.

    supercell_area:
        Area of the rectangular supercell.

    Returns
    -------
    complex
        Fourier coefficient at reciprocal vector G.
    """

    reciprocal_vector = np.asarray(
        reciprocal_vector,
        dtype=float
    )

    rod_positions = np.asarray(
        rod_positions,
        dtype=float
    )

    reciprocal_magnitude = np.linalg.norm(
        reciprocal_vector
    )

    single_rod_fill_fraction = (
        np.pi * rod_radius**2 / supercell_area
    )

    number_of_rods = rod_positions.shape[0]

    if np.isclose(
        reciprocal_magnitude,
        0.0
    ):
        total_fill_fraction = (
            number_of_rods
            * single_rod_fill_fraction
        )

        return (
            background_value
            + (
                rod_value - background_value
            )
            * total_fill_fraction
        )

    shape_factor = (
        2.0
        * j1(
            reciprocal_magnitude * rod_radius
        )
        / (
            reciprocal_magnitude * rod_radius
        )
    )

    phase_sum = np.sum(
        np.exp(
            -1j
            * (
                rod_positions
                @ reciprocal_vector
            )
        )
    )

    return (
        (
            rod_value - background_value
        )
        * single_rod_fill_fraction
        * shape_factor
        * phase_sum
    )


def supercell_convolution_matrix(
    reciprocal_vectors: np.ndarray,
    rod_positions: np.ndarray,
    rod_radius: float,
    rod_value: float,
    background_value: float,
    supercell_area: float
) -> np.ndarray:
    """
    Construct the Fourier-space convolution matrix for a
    circular-rod supercell.

    Matrix element (i, j) equals the Fourier coefficient at
    G_i - G_j.
    """

    reciprocal_vectors = np.asarray(
        reciprocal_vectors,
        dtype=float
    )

    number_of_vectors = reciprocal_vectors.shape[0]

    convolution_matrix = np.empty(
        (
            number_of_vectors,
            number_of_vectors
        ),
        dtype=complex
    )

    for row_index in range(number_of_vectors):
        for column_index in range(number_of_vectors):

            reciprocal_difference = (
                reciprocal_vectors[row_index]
                - reciprocal_vectors[column_index]
            )

            convolution_matrix[
                row_index,
                column_index
            ] = (
                circular_rod_supercell_fourier_coefficient(
                    reciprocal_difference,
                    rod_positions,
                    rod_radius,
                    rod_value,
                    background_value,
                    supercell_area
                )
            )

    return convolution_matrix


def rectangular_reciprocal_vector_set(
    reciprocal_vector_1: np.ndarray,
    reciprocal_vector_2: np.ndarray,
    index_limit_1: int,
    index_limit_2: int
) -> np.ndarray:
    """
    Generate reciprocal-lattice vectors for a rectangular
    supercell.

    The reciprocal vectors are

        G = m b_1 + n b_2,

    where m and n are independently truncated.

    Parameters
    ----------
    reciprocal_vector_1:
        First reciprocal-lattice basis vector.

    reciprocal_vector_2:
        Second reciprocal-lattice basis vector.

    index_limit_1:
        Maximum absolute reciprocal index along b_1.

    index_limit_2:
        Maximum absolute reciprocal index along b_2.

    Returns
    -------
    np.ndarray
        Reciprocal vectors with shape
        ((2 * index_limit_1 + 1)
        * (2 * index_limit_2 + 1), 2).
    """

    if index_limit_1 < 0 or index_limit_2 < 0:
        raise ValueError(
            "Reciprocal index limits must be non-negative."
        )

    reciprocal_vectors = []

    for index_1 in range(
        -index_limit_1,
        index_limit_1 + 1
    ):
        for index_2 in range(
            -index_limit_2,
            index_limit_2 + 1
        ):
            reciprocal_vector = (
                index_1 * reciprocal_vector_1
                + index_2 * reciprocal_vector_2
            )

            reciprocal_vectors.append(
                reciprocal_vector
            )

    return np.asarray(
        reciprocal_vectors,
        dtype=float
    )


def calculate_tm_waveguide_band_structure(
    wavevectors: np.ndarray,
    reciprocal_vectors: np.ndarray,
    permittivity_convolution_matrix: np.ndarray,
    lattice_constant: float,
    number_of_bands: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the TM projected band structure of a
    two-dimensional photonic-crystal waveguide.

    The generalized eigenvalue problem is

        K E = lambda epsilon E,

    where

        K_GG' = |k + G|^2 delta_GG'

    and

        lambda = (omega / c)^2.

    Parameters
    ----------
    wavevectors:
        Bloch wavevectors along the waveguide direction.

    reciprocal_vectors:
        Truncated reciprocal-lattice vector set.

    permittivity_convolution_matrix:
        Fourier-space permittivity convolution matrix.

    lattice_constant:
        Primitive square-lattice constant a.

    number_of_bands:
        Number of lowest-frequency bands to retain.

    Returns
    -------
    frequencies:
        Normalized frequencies omega * a / (2 pi c), with shape
        (number_of_wavevectors, number_of_bands).

    eigenvectors:
        Plane-wave coefficients with shape
        (number_of_wavevectors, number_of_plane_waves,
        number_of_bands).
    """

    wavevectors = np.asarray(
        wavevectors,
        dtype=float
    )

    reciprocal_vectors = np.asarray(
        reciprocal_vectors,
        dtype=float
    )

    permittivity_convolution_matrix = np.asarray(
        permittivity_convolution_matrix,
        dtype=complex
    )

    number_of_plane_waves = reciprocal_vectors.shape[0]

    if number_of_bands < 1:
        raise ValueError(
            "number_of_bands must be positive."
        )

    number_of_bands = min(
        number_of_bands,
        number_of_plane_waves
    )

    frequencies = np.empty(
        (
            wavevectors.shape[0],
            number_of_bands
        ),
        dtype=float
    )

    eigenvectors = np.empty(
        (
            wavevectors.shape[0],
            number_of_plane_waves,
            number_of_bands
        ),
        dtype=complex
    )

    hermitian_permittivity_matrix = (
        0.5
        * (
            permittivity_convolution_matrix
            + permittivity_convolution_matrix.conj().T
        )
    )

    for wavevector_index, wavevector in enumerate(
        wavevectors
    ):
        shifted_vectors = (
            reciprocal_vectors + wavevector
        )

        squared_magnitudes = np.sum(
            shifted_vectors**2,
            axis=1
        )

        kinetic_matrix = np.diag(
            squared_magnitudes
        )

        eigenvalues, mode_coefficients = eigh(
            kinetic_matrix,
            hermitian_permittivity_matrix,
            subset_by_index=(
                0,
                number_of_bands - 1
            )
        )

        eigenvalues = np.maximum(
            np.real(eigenvalues),
            0.0
        )

        frequencies[
            wavevector_index,
            :
        ] = (
            lattice_constant
            * np.sqrt(eigenvalues)
            / (
                2.0 * np.pi
            )
        )

        eigenvectors[
            wavevector_index,
            :,
            :
        ] = mode_coefficients

    return frequencies, eigenvectors
