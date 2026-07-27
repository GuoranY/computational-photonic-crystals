"""
Bloch-analysis utilities for one-dimensional photonic crystals.
"""

from collections.abc import Callable

import numpy as np

from utils.field_transfer_utils import field_unit_cell_matrix
from utils.transfer_matrix_utils import unit_cell_matrix


# ============================================================================
# Bloch function
# ============================================================================

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

    return float(
        np.real(
            0.5 * np.trace(cell_matrix)
        )
    )


# ============================================================================
# Bloch wavevector
# ============================================================================

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

    bloch_value = np.clip(
        bloch_value,
        -1.0,
        1.0,
    )

    return float(
        np.arccos(bloch_value)
        / lattice_constant
    )


def find_band_gaps(
    normalized_frequencies: np.ndarray,
    n_1: float,
    n_2: float,
    d_1: float,
    d_2: float,
    lattice_constant: float,
) -> list[tuple[float, float]]:
    """
    Find approximate photonic band gaps from sampled frequencies.

    A sampled frequency belongs to a forbidden gap when

        |Tr(M_cell) / 2| > 1.

    The returned boundaries correspond to points on the supplied
    frequency grid and are therefore approximate.

    Returns
    -------
    list[tuple[float, float]]
        Approximate lower and upper frequency boundaries of each
        detected band gap.
    """
    if normalized_frequencies.ndim != 1:
        raise ValueError(
            "normalized_frequencies must be a one-dimensional array."
        )

    if normalized_frequencies.size == 0:
        return []

    bloch_functions = np.empty_like(
        normalized_frequencies,
        dtype=float,
    )

    for index, frequency in enumerate(normalized_frequencies):
        bloch_functions[index] = bloch_function(
            normalized_frequency=float(frequency),
            n_1=n_1,
            n_2=n_2,
            d_1=d_1,
            d_2=d_2,
            lattice_constant=lattice_constant,
        )

    forbidden = np.abs(bloch_functions) > 1.0

    transitions = np.diff(
        forbidden.astype(int)
    )

    gap_starts = np.where(
        transitions == 1
    )[0] + 1

    gap_ends = np.where(
        transitions == -1
    )[0]

    if forbidden[0]:
        gap_starts = np.insert(
            gap_starts,
            0,
            0,
        )

    if forbidden[-1]:
        gap_ends = np.append(
            gap_ends,
            normalized_frequencies.size - 1,
        )

    return [
        (
            float(normalized_frequencies[start]),
            float(normalized_frequencies[end]),
        )
        for start, end in zip(
            gap_starts,
            gap_ends,
        )
    ]


# ============================================================================
# Bloch trace function
# ============================================================================

def field_bloch_trace(
    normalized_frequency: float,
    n_1: float,
    d_1: float,
    n_2: float,
    d_2: float,
    lattice_constant: float,
) -> float:
    """
    Evaluate the Bloch trace function

        F(nu) = 1/2 Tr(M_cell)

    using the field-state transfer matrix.
    """
    matrix = field_unit_cell_matrix(
        normalized_frequency=normalized_frequency,
        n_1=n_1,
        d_1=d_1,
        n_2=n_2,
        d_2=d_2,
        lattice_constant=lattice_constant,
    )

    return float(
        np.real(
            0.5 * np.trace(matrix)
        )
    )


# ============================================================================
# Find a root using bisection
# ============================================================================

def bisection_root(
    function: Callable[[float], float],
    left_boundary: float,
    right_boundary: float,
    target_value: float,
    tolerance: float = 1e-12,
    maximum_iterations: int = 200,
) -> float:
    """
    Find a solution of

        function(x) = target_value

    within a specified interval using the bisection method.
    """
    if left_boundary >= right_boundary:
        raise ValueError(
            "left_boundary must be smaller than right_boundary."
        )

    if tolerance <= 0.0:
        raise ValueError(
            "tolerance must be positive."
        )

    if maximum_iterations <= 0:
        raise ValueError(
            "maximum_iterations must be positive."
        )

    left_value = (
        function(left_boundary)
        - target_value
    )

    right_value = (
        function(right_boundary)
        - target_value
    )

    if abs(left_value) < tolerance:
        return left_boundary

    if abs(right_value) < tolerance:
        return right_boundary

    if left_value * right_value > 0.0:
        raise ValueError(
            "The supplied interval does not bracket a root."
        )

    for _ in range(maximum_iterations):
        midpoint = 0.5 * (
            left_boundary
            + right_boundary
        )

        midpoint_value = (
            function(midpoint)
            - target_value
        )

        if abs(midpoint_value) < tolerance:
            return midpoint

        if left_value * midpoint_value <= 0.0:
            right_boundary = midpoint
        else:
            left_boundary = midpoint
            left_value = midpoint_value

    return 0.5 * (
        left_boundary
        + right_boundary
    )


# ============================================================================
# Obtain the band-edge Bloch state
# ============================================================================

def field_band_edge_state(
    normalized_frequency: float,
    n_1: float,
    d_1: float,
    n_2: float,
    d_2: float,
    lattice_constant: float,
    band_edge_tolerance: float = 1e-8,
) -> np.ndarray:
    """
    Calculate the electromagnetic field state at a photonic band edge.

    At a band edge, the Bloch eigenvalue is either +1 or -1:

        M_cell v = lambda v.

    The corresponding state vector is obtained from the approximate
    null space of

        M_cell - lambda I.

    Returns
    -------
    np.ndarray
        Complex electromagnetic state vector [E, H]^T.
    """
    matrix = field_unit_cell_matrix(
        normalized_frequency=normalized_frequency,
        n_1=n_1,
        d_1=d_1,
        n_2=n_2,
        d_2=d_2,
        lattice_constant=lattice_constant,
    )

    trace_value = float(
        np.real(
            0.5 * np.trace(matrix)
        )
    )

    if not np.isclose(
        abs(trace_value),
        1.0,
        atol=band_edge_tolerance,
        rtol=0.0,
    ):
        raise ValueError(
            "The supplied frequency is not sufficiently close "
            "to a photonic band edge."
        )

    bloch_eigenvalue = (
        1.0
        if trace_value >= 0.0
        else -1.0
    )

    shifted_matrix = (
        matrix
        - bloch_eigenvalue
        * np.eye(
            2,
            dtype=complex,
        )
    )

    _, _, conjugate_transpose = np.linalg.svd(
        shifted_matrix
    )

    state = (
        conjugate_transpose
        .conj()
        .T[:, -1]
    )

    # Remove the arbitrary global complex phase using the component
    # with the largest magnitude as the phase reference.
    reference_index = int(
        np.argmax(
            np.abs(state)
        )
    )

    state *= np.exp(
        -1j
        * np.angle(
            state[reference_index]
        )
    )

    # Fix the otherwise arbitrary amplitude of the state vector.
    state /= np.linalg.norm(state)

    return state
