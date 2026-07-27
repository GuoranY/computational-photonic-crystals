"""
Utilities for transfer-matrix calculations using the electromagnetic
field-state representation [E, H]^T.
"""

import numpy as np


def field_layer_matrix(
    refractive_index: float,
    thickness: float,
    normalized_frequency: float,
    lattice_constant: float
) -> np.ndarray:
    """
    Construct the transfer matrix of one homogeneous dielectric layer.

    The electromagnetic state is represented by

        [E, H]^T,

    and the matrix propagates the state from the left boundary of the
    layer to the right boundary.

    Parameters
    ----------
    refractive_index:
        Refractive index of the dielectric layer.

    thickness:
        Physical thickness of the dielectric layer.

    normalized_frequency:
        Normalized frequency

            nu = omega a / (2 pi c).

    lattice_constant:
        Lattice constant a used in the normalized-frequency definition.

    Returns
    -------
    np.ndarray
        The 2 x 2 complex transfer matrix of the layer.
    """

    if refractive_index <= 0.0:
        raise ValueError(
            "refractive_index must be positive."
        )

    if thickness <= 0.0:
        raise ValueError(
            "thickness must be positive."
        )

    if normalized_frequency < 0.0:
        raise ValueError(
            "normalized_frequency cannot be negative."
        )

    if lattice_constant <= 0.0:
        raise ValueError(
            "lattice_constant must be positive."
        )

    phase = (
        2.0
        * np.pi
        * refractive_index
        * normalized_frequency
        * thickness
        / lattice_constant
    )

    return np.array(
        [
            [
                np.cos(phase),
                1j * np.sin(phase) / refractive_index
            ],
            [
                1j * refractive_index * np.sin(phase),
                np.cos(phase)
            ]
        ],
        dtype=complex
    )


def field_unit_cell_matrix(
    normalized_frequency: float,
    n_1: float,
    d_1: float,
    n_2: float,
    d_2: float,
    lattice_constant: float
) -> np.ndarray:
    """
    Construct the field-state transfer matrix of one AB unit cell.

    The physical order of the unit cell is

        material 1 followed by material 2.

    Therefore, the state propagates according to

        state_after_cell
        = M_2 M_1 state_before_cell.

    Parameters
    ----------
    normalized_frequency:
        Normalized frequency

            nu = omega a / (2 pi c).

    n_1:
        Refractive index of material 1.

    d_1:
        Thickness of material 1.

    n_2:
        Refractive index of material 2.

    d_2:
        Thickness of material 2.

    lattice_constant:
        Lattice constant a.

    Returns
    -------
    np.ndarray
        The 2 x 2 field-state transfer matrix of one unit cell.
    """

    layer_1_matrix = field_layer_matrix(
        refractive_index=n_1,
        thickness=d_1,
        normalized_frequency=normalized_frequency,
        lattice_constant=lattice_constant
    )

    layer_2_matrix = field_layer_matrix(
        refractive_index=n_2,
        thickness=d_2,
        normalized_frequency=normalized_frequency,
        lattice_constant=lattice_constant
    )

    return layer_2_matrix @ layer_1_matrix


# ============================================================================
# Total transfer matrix
# ============================================================================

def field_total_transfer_matrix(
    layers: list[tuple[float, float]],
    normalized_frequency: float,
    lattice_constant: float
) -> np.ndarray:
    """
    Multiply the field-state layer matrices in their physical
    propagation order.

    Parameters
    ----------
    layers:
        Sequence of dielectric layers represented as

            (refractive_index, thickness).

    normalized_frequency:
        Normalized frequency

            nu = omega a / (2 pi c).

    lattice_constant:
        Lattice constant a.

    Returns
    -------
    np.ndarray
        The 2 x 2 total field-state transfer matrix.
    """

    total_matrix = np.eye(
        2,
        dtype=complex
    )

    for refractive_index, thickness in layers:
        total_matrix = (
            field_layer_matrix(
                refractive_index=refractive_index,
                thickness=thickness,
                normalized_frequency=normalized_frequency,
                lattice_constant=lattice_constant
            )
            @ total_matrix
        )

    return total_matrix
