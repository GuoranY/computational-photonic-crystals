"""
Utility functions for two-dimensional photonic-crystal supercells,
point-defect cavities, and localized TM-mode reconstruction.
"""

import numpy as np
from scipy.linalg import eigh
from scipy.special import j1


# ============================================================================
# Supercell geometry
# ============================================================================

def square_supercell_rod_positions(
    number_of_cells: int,
    lattice_constant: float,
    remove_center_rod: bool = False,
) -> np.ndarray:
    """
    Generate the positions of dielectric rods in a square supercell.

    The supercell contains

        number_of_cells x number_of_cells

    primitive unit cells. The rod positions are centered around the
    origin.

    Parameters
    ----------
    number_of_cells:
        Number of primitive cells along each spatial direction.
        An odd value is recommended so that one rod lies exactly
        at the origin.

    lattice_constant:
        Primitive-lattice constant.

    remove_center_rod:
        If True, remove the rod located at the origin to create
        a point-defect cavity.

    Returns
    -------
    np.ndarray
        Rod positions with shape (number_of_rods, 2).
    """

    if number_of_cells <= 0:
        raise ValueError(
            "number_of_cells must be positive."
        )

    if lattice_constant <= 0.0:
        raise ValueError(
            "lattice_constant must be positive."
        )

    if (
        remove_center_rod
        and number_of_cells % 2 == 0
    ):
        raise ValueError(
            "An odd number_of_cells is required "
            "to remove a unique center rod."
        )

    half_width = (
        number_of_cells - 1
    ) / 2.0

    coordinate_indices = np.arange(
        number_of_cells,
        dtype=float,
    ) - half_width

    coordinates = (
        coordinate_indices
        * lattice_constant
    )

    rod_positions = []

    for x_position in coordinates:
        for y_position in coordinates:
            position = np.array(
                [
                    x_position,
                    y_position,
                ]
            )

            is_center_rod = np.allclose(
                position,
                0.0,
            )

            if (
                remove_center_rod
                and is_center_rod
            ):
                continue

            rod_positions.append(
                position
            )

    return np.asarray(
        rod_positions,
        dtype=float,
    )


# ============================================================================
# Reciprocal-lattice vector generation
# ============================================================================

def supercell_reciprocal_vector_set(
    supercell_length: float,
    reciprocal_index_limit: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate a square truncated reciprocal-lattice basis for
    a square supercell.

    The reciprocal vectors are

        G = m B_1 + n B_2,

    where

        B_1 = (2 pi / L, 0),
        B_2 = (0, 2 pi / L),

    and L is the supercell side length.

    Parameters
    ----------
    supercell_length:
        Side length of the square supercell.

    reciprocal_index_limit:
        Maximum absolute reciprocal-lattice index.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        reciprocal_vectors:
            Reciprocal vectors with shape (number_of_vectors, 2).

        reciprocal_indices:
            Integer reciprocal indices with shape
            (number_of_vectors, 2).
    """

    if supercell_length <= 0.0:
        raise ValueError(
            "supercell_length must be positive."
        )

    if reciprocal_index_limit < 0:
        raise ValueError(
            "reciprocal_index_limit must be nonnegative."
        )

    reciprocal_vector_1 = np.array(
        [
            2.0 * np.pi / supercell_length,
            0.0,
        ]
    )

    reciprocal_vector_2 = np.array(
        [
            0.0,
            2.0 * np.pi / supercell_length,
        ]
    )

    reciprocal_vectors = []
    reciprocal_indices = []

    for m_index in range(
        -reciprocal_index_limit,
        reciprocal_index_limit + 1,
    ):
        for n_index in range(
            -reciprocal_index_limit,
            reciprocal_index_limit + 1,
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
# Supercell Fourier coefficients
# ============================================================================

def circular_rod_form_factor(
    reciprocal_magnitudes: np.ndarray,
    rod_radius: float,
) -> np.ndarray:
    """
    Calculate the circular-inclusion form factor

        2 J_1(|G| r) / (|G| r).

    The zero-wavevector value is defined by its limiting value 1.
    """

    reciprocal_magnitudes = np.asarray(
        reciprocal_magnitudes,
        dtype=float,
    )

    argument = (
        reciprocal_magnitudes
        * rod_radius
    )

    form_factor = np.ones_like(
        reciprocal_magnitudes,
        dtype=float,
    )

    nonzero_mask = ~np.isclose(
        reciprocal_magnitudes,
        0.0,
    )

    form_factor[nonzero_mask] = (
        2.0
        * j1(
            argument[nonzero_mask]
        )
        / argument[nonzero_mask]
    )

    return form_factor


def supercell_material_fourier_coefficients(
    reciprocal_vectors: np.ndarray,
    supercell_length: float,
    rod_radius: float,
    rod_positions: np.ndarray,
    rod_value: float,
    background_value: float,
) -> np.ndarray:
    """
    Calculate Fourier coefficients of a piecewise-constant material
    function for multiple circular rods inside a square supercell.

    The coefficient is

        f_G =
            f_background delta_G0
            + (f_rod - f_background)
              (A_rod / A_supercell)
              F(G)
              sum_j exp(-i G . r_j),

    where F(G) is the circular-rod form factor.

    Parameters
    ----------
    reciprocal_vectors:
        Reciprocal vectors with shape (..., 2).

    supercell_length:
        Side length of the square supercell.

    rod_radius:
        Radius of every dielectric rod.

    rod_positions:
        Rod-center positions with shape (number_of_rods, 2).

    rod_value:
        Material value inside each rod.

    background_value:
        Material value in the surrounding background.

    Returns
    -------
    np.ndarray
        Complex Fourier coefficients with shape
        reciprocal_vectors.shape[:-1].
    """

    reciprocal_vectors = np.asarray(
        reciprocal_vectors,
        dtype=float,
    )

    rod_positions = np.asarray(
        rod_positions,
        dtype=float,
    )

    if reciprocal_vectors.shape[-1] != 2:
        raise ValueError(
            "reciprocal_vectors must have shape (..., 2)."
        )

    if (
        rod_positions.ndim != 2
        or rod_positions.shape[1] != 2
    ):
        raise ValueError(
            "rod_positions must have shape "
            "(number_of_rods, 2)."
        )

    if supercell_length <= 0.0:
        raise ValueError(
            "supercell_length must be positive."
        )

    if rod_radius <= 0.0:
        raise ValueError(
            "rod_radius must be positive."
        )

    supercell_area = (
        supercell_length**2
    )

    rod_area = (
        np.pi * rod_radius**2
    )

    reciprocal_magnitudes = np.linalg.norm(
        reciprocal_vectors,
        axis=-1,
    )

    form_factor = circular_rod_form_factor(
        reciprocal_magnitudes=(
            reciprocal_magnitudes
        ),
        rod_radius=rod_radius,
    )

    phase_arguments = np.tensordot(
        reciprocal_vectors,
        rod_positions,
        axes=(
            [-1],
            [1],
        ),
    )

    structure_factor = np.sum(
        np.exp(
            -1.0j * phase_arguments
        ),
        axis=-1,
    )

    coefficients = (
        (
            rod_value
            - background_value
        )
        * rod_area
        / supercell_area
        * form_factor
        * structure_factor
    )

    zero_wavevector_mask = np.isclose(
        reciprocal_magnitudes,
        0.0,
    )

    coefficients = np.asarray(
        coefficients,
        dtype=complex,
    )

    coefficients[
        zero_wavevector_mask
    ] += background_value

    return coefficients


def supercell_convolution_matrix(
    reciprocal_indices: np.ndarray,
    supercell_length: float,
    rod_radius: float,
    rod_positions: np.ndarray,
    rod_value: float,
    background_value: float,
) -> np.ndarray:
    """
    Construct the reciprocal-space convolution matrix

        f_(G-G')

    for a square photonic-crystal supercell.
    """

    reciprocal_indices = np.asarray(
        reciprocal_indices,
        dtype=int,
    )

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
        / supercell_length
        * index_differences
    )

    flattened_differences = (
        reciprocal_differences.reshape(
            -1,
            2,
        )
    )

    coefficients = (
        supercell_material_fourier_coefficients(
            reciprocal_vectors=(
                flattened_differences
            ),
            supercell_length=(
                supercell_length
            ),
            rod_radius=rod_radius,
            rod_positions=rod_positions,
            rod_value=rod_value,
            background_value=(
                background_value
            ),
        )
    )

    convolution_matrix = coefficients.reshape(
        number_of_vectors,
        number_of_vectors,
    )

    return 0.5 * (
        convolution_matrix
        + convolution_matrix.conj().T
    )


# ============================================================================
# Supercell Brillouin-zone path
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


def square_supercell_wavevector_path(
    supercell_length: float,
    points_per_segment: int,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """
    Construct the square-supercell high-symmetry path

        Gamma -> X_s -> M_s -> Gamma.
    """

    gamma_point = np.array(
        [
            0.0,
            0.0,
        ]
    )

    x_point = np.array(
        [
            np.pi / supercell_length,
            0.0,
        ]
    )

    m_point = np.array(
        [
            np.pi / supercell_length,
            np.pi / supercell_length,
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
        number_of_points=(
            points_per_segment + 1
        ),
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
        r"$X_{\mathrm{s}}$",
        r"$M_{\mathrm{s}}$",
        r"$\Gamma$",
    ]

    return (
        wavevectors,
        symmetry_positions,
        symmetry_labels,
    )


# ============================================================================
# TM supercell eigensystem
# ============================================================================

def normalized_frequencies_from_eigenvalues(
    eigenvalues: np.ndarray,
    normalization_length: float,
) -> np.ndarray:
    """
    Convert eigenvalues (omega / c)^2 into normalized frequencies

        omega a / (2 pi c),

    where normalization_length is normally the primitive-lattice
    constant a.
    """

    cleaned_eigenvalues = np.maximum(
        np.real(eigenvalues),
        0.0,
    )

    return (
        normalization_length
        * np.sqrt(
            cleaned_eigenvalues
        )
        / (2.0 * np.pi)
    )


def tm_supercell_eigensystem(
    wavevector: np.ndarray,
    reciprocal_vectors: np.ndarray,
    dielectric_matrix: np.ndarray,
    normalization_length: float,
    number_of_bands: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate TM eigenfrequencies and plane-wave eigenvectors
    for a photonic-crystal supercell.

    For TM polarization, the out-of-plane electric field E_z
    satisfies

        A E = (omega / c)^2 B E,

    where

        A_GG' = |k + G|^2 delta_GG'

    and

        B_GG' = epsilon_(G-G').
    """

    number_of_plane_waves = len(
        reciprocal_vectors
    )

    if number_of_bands > number_of_plane_waves:
        raise ValueError(
            "number_of_bands cannot exceed the "
            "number of plane waves."
        )

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

    eigenvalues, eigenvectors = eigh(
        operator_matrix,
        dielectric_matrix,
        subset_by_index=[
            0,
            number_of_bands - 1,
        ],
        check_finite=False,
    )

    frequencies = (
        normalized_frequencies_from_eigenvalues(
            eigenvalues=eigenvalues,
            normalization_length=(
                normalization_length
            ),
        )
    )

    return (
        frequencies,
        eigenvectors,
    )


def calculate_tm_supercell_band_structure(
    wavevectors: np.ndarray,
    reciprocal_vectors: np.ndarray,
    dielectric_matrix: np.ndarray,
    normalization_length: float,
    number_of_bands: int,
) -> np.ndarray:
    """
    Calculate TM bands along a supercell wavevector path.
    """

    number_of_wavevectors = len(
        wavevectors
    )

    frequencies = np.zeros(
        (
            number_of_wavevectors,
            number_of_bands,
        )
    )

    for wavevector_index, wavevector in enumerate(
        wavevectors
    ):
        (
            current_frequencies,
            _,
        ) = tm_supercell_eigensystem(
            wavevector=wavevector,
            reciprocal_vectors=(
                reciprocal_vectors
            ),
            dielectric_matrix=(
                dielectric_matrix
            ),
            normalization_length=(
                normalization_length
            ),
            number_of_bands=(
                number_of_bands
            ),
        )

        frequencies[
            wavevector_index
        ] = current_frequencies

    return frequencies


# ============================================================================
# Field reconstruction and localization
# ============================================================================

def reconstruct_tm_field(
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    wavevector: np.ndarray,
    reciprocal_vectors: np.ndarray,
    eigenvector: np.ndarray,
    chunk_size: int = 4000,
) -> np.ndarray:
    """
    Reconstruct the real-space TM field

        E_z(r) =
            sum_G E_G exp[i(k + G).r].

    The calculation is performed in chunks to reduce memory usage.
    """

    x_grid, y_grid = np.meshgrid(
        x_coordinates,
        y_coordinates,
    )

    positions = np.column_stack(
        [
            x_grid.ravel(),
            y_grid.ravel(),
        ]
    )

    shifted_vectors = (
        reciprocal_vectors
        + wavevector[np.newaxis, :]
    )

    flattened_field = np.zeros(
        len(positions),
        dtype=complex,
    )

    for start_index in range(
        0,
        len(positions),
        chunk_size,
    ):
        end_index = min(
            start_index + chunk_size,
            len(positions),
        )

        current_positions = positions[
            start_index:end_index
        ]

        phase_matrix = np.exp(
            1.0j
            * (
                current_positions
                @ shifted_vectors.T
            )
        )

        flattened_field[
            start_index:end_index
        ] = (
            phase_matrix
            @ eigenvector
        )

    field = flattened_field.reshape(
        x_grid.shape
    )

    return field


def normalized_field_intensity(
    field: np.ndarray,
) -> np.ndarray:
    """
    Calculate and normalize the field intensity |E|^2.
    """

    intensity = np.abs(
        field
    )**2

    maximum_intensity = np.max(
        intensity
    )

    if maximum_intensity <= 0.0:
        return intensity

    return (
        intensity
        / maximum_intensity
    )


def field_localization_fraction(
    intensity: np.ndarray,
    x_coordinates: np.ndarray,
    y_coordinates: np.ndarray,
    localization_radius: float,
) -> float:
    """
    Calculate the fraction of the total sampled field intensity
    contained inside a circle centered on the point defect.
    """

    x_grid, y_grid = np.meshgrid(
        x_coordinates,
        y_coordinates,
    )

    radial_distance = np.sqrt(
        x_grid**2
        + y_grid**2
    )

    localization_mask = (
        radial_distance
        <= localization_radius
    )

    total_intensity = np.sum(
        intensity
    )

    if total_intensity <= 0.0:
        return 0.0

    localized_intensity = np.sum(
        intensity[
            localization_mask
        ]
    )

    return float(
        localized_intensity
        / total_intensity
    )


# ============================================================================
# TE supercell eigensystem
# ============================================================================

def te_supercell_eigensystem(
    wavevector: np.ndarray,
    reciprocal_vectors: np.ndarray,
    inverse_dielectric_matrix: np.ndarray,
    normalization_length: float,
    number_of_bands: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate TE eigenfrequencies and plane-wave eigenvectors
    for a two-dimensional photonic-crystal supercell.

    For TE polarization, the out-of-plane magnetic field H_z
    satisfies the plane-wave eigenvalue equation

        sum_G' [
            (k + G) . (k + G')
            epsilon_inverse_(G-G')
        ] H_G'
        =
        (omega / c)^2 H_G.

    Parameters
    ----------
    wavevector:
        Bloch wavevector with shape (2,).

    reciprocal_vectors:
        Truncated reciprocal-vector set with shape
        (number_of_vectors, 2).

    inverse_dielectric_matrix:
        Convolution matrix of the inverse permittivity.

    normalization_length:
        Length used to define the normalized frequency.
        This should normally be the primitive lattice constant.

    number_of_bands:
        Number of eigenmodes to return.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Normalized eigenfrequencies and the corresponding
        plane-wave eigenvectors.
    """

    wavevector = np.asarray(
        wavevector,
        dtype=float,
    )

    reciprocal_vectors = np.asarray(
        reciprocal_vectors,
        dtype=float,
    )

    inverse_dielectric_matrix = np.asarray(
        inverse_dielectric_matrix,
        dtype=complex,
    )

    number_of_plane_waves = len(
        reciprocal_vectors
    )

    if number_of_bands <= 0:
        raise ValueError(
            "number_of_bands must be positive."
        )

    if number_of_bands > number_of_plane_waves:
        raise ValueError(
            "number_of_bands cannot exceed the "
            "number of plane waves."
        )

    shifted_vectors = (
        reciprocal_vectors
        + wavevector[np.newaxis, :]
    )

    dot_product_matrix = (
        shifted_vectors
        @ shifted_vectors.T
    )

    operator_matrix = (
        dot_product_matrix
        * inverse_dielectric_matrix
    )

    operator_matrix = 0.5 * (
        operator_matrix
        + operator_matrix.conj().T
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        operator_matrix
    )

    eigenvalues = np.real(
        eigenvalues
    )

    eigenvalues = np.maximum(
        eigenvalues,
        0.0,
    )

    frequencies = (
        normalization_length
        * np.sqrt(
            eigenvalues
        )
        / (
            2.0 * np.pi
        )
    )

    return (
        frequencies[
            :number_of_bands
        ],
        eigenvectors[
            :,
            :number_of_bands
        ],
    )


def calculate_te_supercell_band_structure(
    wavevectors: np.ndarray,
    reciprocal_vectors: np.ndarray,
    inverse_dielectric_matrix: np.ndarray,
    normalization_length: float,
    number_of_bands: int,
) -> np.ndarray:
    """
    Calculate the TE band structure of a photonic-crystal
    supercell along a supplied wavevector path.
    """

    number_of_wavevectors = len(
        wavevectors
    )

    frequencies = np.zeros(
        (
            number_of_wavevectors,
            number_of_bands,
        )
    )

    for (
        wavevector_index,
        wavevector,
    ) in enumerate(
        wavevectors
    ):
        (
            current_frequencies,
            _,
        ) = te_supercell_eigensystem(
            wavevector=wavevector,
            reciprocal_vectors=(
                reciprocal_vectors
            ),
            inverse_dielectric_matrix=(
                inverse_dielectric_matrix
            ),
            normalization_length=(
                normalization_length
            ),
            number_of_bands=(
                number_of_bands
            ),
        )

        frequencies[
            wavevector_index
        ] = current_frequencies

    return frequencies
