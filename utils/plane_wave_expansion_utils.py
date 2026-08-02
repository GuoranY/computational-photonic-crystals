"""
Utility functions for two-dimensional photonic-band calculations
using the plane-wave expansion method.
"""

import numpy as np
from scipy.linalg import eigh

from utils.fourier_utils import circular_inclusion_fourier_coefficients


# ============================================================================
# Reciprocal-lattice vector generation
# ============================================================================

def reciprocal_lattice_vector_set(
    reciprocal_vector_1: np.ndarray,
    reciprocal_vector_2: np.ndarray,
    index_limit: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate a square truncated set of reciprocal-lattice vectors.

    The reciprocal vectors are

        G = m b_1 + n b_2,

    where

        -index_limit <= m, n <= index_limit.

    Parameters
    ----------
    reciprocal_vector_1:
        First reciprocal-lattice basis vector.

    reciprocal_vector_2:
        Second reciprocal-lattice basis vector.

    index_limit:
        Maximum absolute reciprocal-lattice index.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        reciprocal_vectors:
            Array with shape (number_of_vectors, 2).

        reciprocal_indices:
            Integer index pairs (m, n), with shape
            (number_of_vectors, 2).
    """

    reciprocal_vectors = []
    reciprocal_indices = []

    for m_index in range(
        -index_limit,
        index_limit + 1,
    ):
        for n_index in range(
            -index_limit,
            index_limit + 1,
        ):
            reciprocal_vector = (
                m_index * reciprocal_vector_1
                + n_index * reciprocal_vector_2
            )

            reciprocal_vectors.append(
                reciprocal_vector
            )

            reciprocal_indices.append(
                [
                    m_index,
                    n_index,
                ]
            )

    return (
        np.asarray(
            reciprocal_vectors,
            dtype=float,
        ),
        np.asarray(
            reciprocal_indices,
            dtype=int,
        ),
    )


# ============================================================================
# High-symmetry wavevector path
# ============================================================================

def interpolate_wavevector_segment(
    start_point: np.ndarray,
    end_point: np.ndarray,
    number_of_points: int,
    include_endpoint: bool = False,
) -> np.ndarray:
    """
    Generate equally spaced wavevectors along one path segment.
    """

    return np.linspace(
        start_point,
        end_point,
        number_of_points,
        endpoint=include_endpoint,
    )


def square_lattice_wavevector_path(
    lattice_constant: float,
    points_per_segment: int,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Construct the square-lattice high-symmetry path

        Gamma -> X -> M -> Gamma.

    Parameters
    ----------
    lattice_constant:
        Real-space lattice constant.

    points_per_segment:
        Number of sampled wavevectors on each path segment.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, list[str]]
        wavevectors:
            Sampled wavevectors along the complete path.

        symmetry_positions:
            Horizontal-axis positions of Gamma, X, M, and Gamma.

        symmetry_labels:
            Corresponding high-symmetry-point labels.
    """

    gamma_point = np.array(
        [
            0.0,
            0.0,
        ]
    )

    x_point = np.array(
        [
            np.pi / lattice_constant,
            0.0,
        ]
    )

    m_point = np.array(
        [
            np.pi / lattice_constant,
            np.pi / lattice_constant,
        ]
    )

    gamma_to_x = interpolate_wavevector_segment(
        start_point=gamma_point,
        end_point=x_point,
        number_of_points=points_per_segment,
        include_endpoint=False,
    )

    x_to_m = interpolate_wavevector_segment(
        start_point=x_point,
        end_point=m_point,
        number_of_points=points_per_segment,
        include_endpoint=False,
    )

    m_to_gamma = interpolate_wavevector_segment(
        start_point=m_point,
        end_point=gamma_point,
        number_of_points=points_per_segment + 1,
        include_endpoint=True,
    )

    wavevectors = np.vstack(
        [
            gamma_to_x,
            x_to_m,
            m_to_gamma,
        ]
    )

    symmetry_positions = np.array(
        [
            0,
            points_per_segment,
            2 * points_per_segment,
            3 * points_per_segment,
        ],
        dtype=int,
    )

    symmetry_labels = [
        r"$\Gamma$",
        r"$X$",
        r"$M$",
        r"$\Gamma$",
    ]

    return (
        wavevectors,
        symmetry_positions,
        symmetry_labels,
    )


# ============================================================================
# Dielectric convolution matrices
# ============================================================================

def dielectric_convolution_matrix(
    reciprocal_indices: np.ndarray,
    lattice_constant: float,
    rod_radius: float,
    rod_permittivity: float,
    background_permittivity: float,
    inverse: bool = False,
) -> np.ndarray:
    """
    Construct the reciprocal-space convolution matrix containing

        epsilon_(G-G')

    or

        (1 / epsilon)_(G-G').
    """

    number_of_vectors = len(
        reciprocal_indices
    )

    index_differences = (
        reciprocal_indices[:, np.newaxis, :]
        - reciprocal_indices[np.newaxis, :, :]
    )

    reciprocal_differences = (
        2.0
        * np.pi
        / lattice_constant
        * index_differences
    )

    flattened_differences = (
        reciprocal_differences.reshape(
            -1,
            2,
        )
    )

    if inverse:
        inclusion_value = (
            1.0 / rod_permittivity
        )

        background_value = (
            1.0 / background_permittivity
        )

    else:
        inclusion_value = rod_permittivity
        background_value = background_permittivity

    coefficients = (
        circular_inclusion_fourier_coefficients(
            reciprocal_vectors=(
                flattened_differences
            ),
            lattice_constant=lattice_constant,
            radius=rod_radius,
            inclusion_value=inclusion_value,
            background_value=background_value,
        )
    )

    convolution_matrix = coefficients.reshape(
        number_of_vectors,
        number_of_vectors,
    )

    return convolution_matrix


# ============================================================================
# TE and TM plane-wave matrices
# ============================================================================

def te_plane_wave_matrix(
    wavevector: np.ndarray,
    reciprocal_vectors: np.ndarray,
    inverse_dielectric_matrix: np.ndarray,
) -> np.ndarray:
    """
    Construct the TE-polarization plane-wave matrix.

    For TE polarization, H_z is expanded in plane waves and

        A_GG'
        =
        (k + G) . (k + G')
        (1 / epsilon)_(G-G').
    """

    shifted_vectors = (
        reciprocal_vectors
        + wavevector[np.newaxis, :]
    )

    dot_product_matrix = (
        shifted_vectors
        @ shifted_vectors.T
    )

    matrix = (
        dot_product_matrix
        * inverse_dielectric_matrix
    )

    return 0.5 * (
        matrix + matrix.T
    )


def tm_plane_wave_matrices(
    wavevector: np.ndarray,
    reciprocal_vectors: np.ndarray,
    dielectric_matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Construct the TM generalized-eigenvalue matrices.

    For TM polarization, E_z satisfies

        A E = (omega / c)^2 B E,

    where

        A_GG' = |k + G|^2 delta_GG'

    and

        B_GG' = epsilon_(G-G').
    """

    shifted_vectors = (
        reciprocal_vectors
        + wavevector[np.newaxis, :]
    )

    squared_magnitudes = np.sum(
        shifted_vectors**2,
        axis=1,
    )

    operator_matrix = np.diag(
        squared_magnitudes
    )

    return (
        operator_matrix,
        dielectric_matrix,
    )


# ============================================================================
# Eigenfrequency calculation
# ============================================================================

def normalized_frequencies_from_eigenvalues(
    eigenvalues: np.ndarray,
    lattice_constant: float,
) -> np.ndarray:
    """
    Convert eigenvalues (omega / c)^2 into normalized frequencies

        omega a / (2 pi c).
    """

    cleaned_eigenvalues = np.maximum(
        np.real(eigenvalues),
        0.0,
    )

    return (
        lattice_constant
        * np.sqrt(cleaned_eigenvalues)
        / (2.0 * np.pi)
    )


def te_eigenfrequencies(
    wavevector: np.ndarray,
    reciprocal_vectors: np.ndarray,
    inverse_dielectric_matrix: np.ndarray,
    lattice_constant: float,
    number_of_bands: int,
) -> np.ndarray:
    """
    Calculate the lowest TE normalized eigenfrequencies.
    """

    matrix = te_plane_wave_matrix(
        wavevector=wavevector,
        reciprocal_vectors=reciprocal_vectors,
        inverse_dielectric_matrix=(
            inverse_dielectric_matrix
        ),
    )

    eigenvalues = eigh(
        matrix,
        eigvals_only=True,
        check_finite=False,
    )

    frequencies = (
        normalized_frequencies_from_eigenvalues(
            eigenvalues=eigenvalues,
            lattice_constant=lattice_constant,
        )
    )

    return frequencies[:number_of_bands]


def tm_eigenfrequencies(
    wavevector: np.ndarray,
    reciprocal_vectors: np.ndarray,
    dielectric_matrix: np.ndarray,
    lattice_constant: float,
    number_of_bands: int,
) -> np.ndarray:
    """
    Calculate the lowest TM normalized eigenfrequencies.
    """

    operator_matrix, material_matrix = (
        tm_plane_wave_matrices(
            wavevector=wavevector,
            reciprocal_vectors=reciprocal_vectors,
            dielectric_matrix=dielectric_matrix,
        )
    )

    eigenvalues = eigh(
        operator_matrix,
        material_matrix,
        eigvals_only=True,
        check_finite=False,
    )

    frequencies = (
        normalized_frequencies_from_eigenvalues(
            eigenvalues=eigenvalues,
            lattice_constant=lattice_constant,
        )
    )

    return frequencies[:number_of_bands]


def calculate_te_tm_band_structure(
    wavevectors: np.ndarray,
    reciprocal_vectors: np.ndarray,
    dielectric_matrix: np.ndarray,
    inverse_dielectric_matrix: np.ndarray,
    lattice_constant: float,
    number_of_bands: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate TE and TM photonic bands along a wavevector path.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        TE and TM band arrays, each with shape

            (number_of_wavevectors, number_of_bands).
    """

    number_of_wavevectors = len(
        wavevectors
    )

    te_bands = np.zeros(
        (
            number_of_wavevectors,
            number_of_bands,
        )
    )

    tm_bands = np.zeros_like(
        te_bands
    )

    for wavevector_index, wavevector in enumerate(
        wavevectors
    ):
        te_bands[
            wavevector_index
        ] = te_eigenfrequencies(
            wavevector=wavevector,
            reciprocal_vectors=reciprocal_vectors,
            inverse_dielectric_matrix=(
                inverse_dielectric_matrix
            ),
            lattice_constant=lattice_constant,
            number_of_bands=number_of_bands,
        )

        tm_bands[
            wavevector_index
        ] = tm_eigenfrequencies(
            wavevector=wavevector,
            reciprocal_vectors=reciprocal_vectors,
            dielectric_matrix=dielectric_matrix,
            lattice_constant=lattice_constant,
            number_of_bands=number_of_bands,
        )

    return (
        te_bands,
        tm_bands,
    )
