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
