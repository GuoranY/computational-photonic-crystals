"""
Transmission utilities for finite one-dimensional photonic crystals.
"""

import numpy as np

from utils.transfer_matrix_utils import (
    interface_matrix,
    propagation_matrix,
)


# ============================================================================
# Finite-crystal matrix
# ============================================================================

def finite_crystal_matrix(
    normalized_frequency: float,
    n_1: float,
    n_2: float,
    d_1: float,
    d_2: float,
    lattice_constant: float,
    number_of_cells: int,
    n_incident: float,
    n_exit: float,
) -> np.ndarray:
    """
    Construct the transfer matrix of the complete finite crystal,
    including the incident and exit interfaces.
    """
    if number_of_cells < 1:
        raise ValueError("number_of_cells must be at least 1.")

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

    interface_incident_1 = interface_matrix(
        refractive_index_left=n_incident,
        refractive_index_right=n_1,
    )

    interface_12 = interface_matrix(
        refractive_index_left=n_1,
        refractive_index_right=n_2,
    )

    interface_21 = interface_matrix(
        refractive_index_left=n_2,
        refractive_index_right=n_1,
    )

    interface_2_exit = interface_matrix(
        refractive_index_left=n_2,
        refractive_index_right=n_exit,
    )

    total_matrix = interface_incident_1

    for cell_index in range(number_of_cells):
        total_matrix = (
            total_matrix
            @ propagation_1
            @ interface_12
            @ propagation_2
        )

        if cell_index < number_of_cells - 1:
            total_matrix = total_matrix @ interface_21

    return total_matrix @ interface_2_exit


# ============================================================================
# Transmission spectrum
# ============================================================================

def transmission_spectrum(
    normalized_frequencies: np.ndarray,
    n_1: float,
    n_2: float,
    d_1: float,
    d_2: float,
    lattice_constant: float,
    number_of_cells: int,
    n_incident: float,
    n_exit: float,
) -> np.ndarray:
    """
    Calculate the power transmission spectrum of a finite
    one-dimensional photonic crystal.

    Parameters
    ----------
    normalized_frequencies:
        One-dimensional array of normalized frequencies.
    n_1:
        Refractive index of the first dielectric layer.
    n_2:
        Refractive index of the second dielectric layer.
    d_1:
        Thickness of the first dielectric layer.
    d_2:
        Thickness of the second dielectric layer.
    lattice_constant:
        Lattice constant of the periodic structure.
    number_of_cells:
        Number of unit cells in the finite crystal.
    n_incident:
        Refractive index of the incident medium.
    n_exit:
        Refractive index of the exit medium.

    Returns
    -------
    np.ndarray
        Power transmission coefficients corresponding to the input
        normalized frequencies.
    """
    frequencies = np.asarray(
        normalized_frequencies,
        dtype=float,
    )

    if frequencies.ndim != 1:
        raise ValueError(
            "normalized_frequencies must be a one-dimensional array."
        )

    if n_incident <= 0.0 or n_exit <= 0.0:
        raise ValueError(
            "n_incident and n_exit must be positive."
        )

    transmissions = np.empty_like(frequencies)

    for index, frequency in enumerate(frequencies):
        total_matrix = finite_crystal_matrix(
            normalized_frequency=frequency,
            n_1=n_1,
            n_2=n_2,
            d_1=d_1,
            d_2=d_2,
            lattice_constant=lattice_constant,
            number_of_cells=number_of_cells,
            n_incident=n_incident,
            n_exit=n_exit,
        )

        transmission_amplitude = 1.0 / total_matrix[0, 0]

        transmissions[index] = (
            n_exit
            / n_incident
            * np.abs(transmission_amplitude) ** 2
        )

    return transmissions
