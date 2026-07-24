"""
Transfer-matrix utilities for one-dimensional photonic crystals.
"""

import numpy as np


# ============================================================================
# Propagation matrix
# ============================================================================

def propagation_matrix(
    refractive_index: float,
    thickness: float,
    normalized_frequency: float,
    lattice_constant: float,
) -> np.ndarray:
    """
    Construct a propagation matrix that maps the wave amplitudes
    on the right side of a dielectric layer to those on the left side.
    """
    phase = (
        2.0
        * np.pi
        * normalized_frequency
        * refractive_index
        * thickness
        / lattice_constant
    )

    return np.array(
        [
            [np.exp(-1j * phase), 0.0],
            [0.0, np.exp(1j * phase)],
        ],
        dtype=complex,
    )


# ============================================================================
# Interface matrix
# ============================================================================

def interface_matrix(
    refractive_index_left: float,
    refractive_index_right: float,
) -> np.ndarray:
    """
    Construct an interface matrix that maps the wave amplitudes
    on the right side of an interface to those on the left side.
    """
    ratio = refractive_index_right / refractive_index_left

    return 0.5 * np.array(
        [
            [1.0 + ratio, 1.0 - ratio],
            [1.0 - ratio, 1.0 + ratio],
        ],
        dtype=complex,
    )


# ============================================================================
# Unit-cell matrix
# ============================================================================

def unit_cell_matrix(
    normalized_frequency: float,
    n_1: float,
    n_2: float,
    d_1: float,
    d_2: float,
    lattice_constant: float,
) -> np.ndarray:
    """
    Construct the transfer matrix of one periodic unit cell.
    """
    propagation_1 = propagation_matrix(
        refractive_index=n_1,
        thickness=d_1,
        normalized_frequency=normalized_frequency,
        lattice_constant=lattice_constant,
    )

    propagation_2 = propagation_matrix(
        refractive_index=n_2,
        thickness=d_2,
        normalized_frequency=normalized_frequency,
        lattice_constant=lattice_constant,
    )

    interface_12 = interface_matrix(
        refractive_index_left=n_1,
        refractive_index_right=n_2,
    )

    interface_21 = interface_matrix(
        refractive_index_left=n_2,
        refractive_index_right=n_1,
    )

    return (
        propagation_1
        @ interface_12
        @ propagation_2
        @ interface_21
    )
