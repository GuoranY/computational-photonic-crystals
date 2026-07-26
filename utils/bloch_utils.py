"""
Bloch-analysis utilities for one-dimensional photonic crystals.
"""

import numpy as np

from utils.transfer_matrix_utils import unit_cell_matrix


# =============================================================================
# Bloch function
# =============================================================================

def bloch_function(
    normalized_frequency: float,
    n_1: float,
    n_2: float,
    d_1: float,
    d_2: float,
    lattice_constant: float,
) -> float:
    """
    Calculate the Bloch function from the unit-cell transfer matrix.
    """
    cell_matrix = unit_cell_matrix(
        normalized_frequency=normalized_frequency,
        n_1=n_1,
        n_2=n_2,
        d_1=d_1,
        d_2=d_2,
        lattice_constant=lattice_constant,
    )

    return float(np.real(np.trace(cell_matrix) / 2.0))


# =============================================================================
# Bloch wavevector
# =============================================================================

def bloch_wavevector(
    normalized_frequency: float,
    n_1: float,
    n_2: float,
    d_1: float,
    d_2: float,
    lattice_constant: float,
    tolerance: float = 1e-10,
) -> float:
    """
    Calculate the real Bloch wavevector in an allowed frequency band.

    Returns np.nan inside a forbidden frequency gap.
    """
    bloch_value = bloch_function(
        normalized_frequency=normalized_frequency,
        n_1=n_1,
        n_2=n_2,
        d_1=d_1,
        d_2=d_2,
        lattice_constant=lattice_constant,
    )

    if np.abs(bloch_value) > 1.0 + tolerance:
        return np.nan

    bloch_value = np.clip(bloch_value, -1.0, 1.0)

    return float(np.arccos(bloch_value) / lattice_constant)


def find_band_gaps(
    normalized_frequencies: np.ndarray,
    n_1: float,
    n_2: float,
    d_1: float,
    d_2: float,
    lattice_constant: float,
) -> list[tuple[float, float]]:
    """
    Find photonic band gaps from the unit-cell transfer matrix.

    A frequency belongs to a forbidden gap when

        |Tr(M_cell) / 2| > 1.

    Returns
    -------
    list[tuple[float, float]]
        Lower and upper frequency boundaries of each detected band gap.
    """
    bloch_functions = np.empty_like(normalized_frequencies)

    for index, frequency in enumerate(normalized_frequencies):
        cell_matrix = unit_cell_matrix(
            normalized_frequency=frequency,
            n_1=n_1,
            n_2=n_2,
            d_1=d_1,
            d_2=d_2,
            lattice_constant=lattice_constant,
        )

        bloch_functions[index] = np.real(
            0.5 * np.trace(cell_matrix)
        )

    forbidden = np.abs(bloch_functions) > 1.0

    transitions = np.diff(forbidden.astype(int))

    gap_starts = np.where(transitions == 1)[0] + 1
    gap_ends = np.where(transitions == -1)[0]

    if forbidden[0]:
        gap_starts = np.insert(gap_starts, 0, 0)

    if forbidden[-1]:
        gap_ends = np.append(
            gap_ends,
            len(normalized_frequencies) - 1,
        )

    return [
        (
            float(normalized_frequencies[start]),
            float(normalized_frequencies[end]),
        )
        for start, end in zip(gap_starts, gap_ends)
    ]
