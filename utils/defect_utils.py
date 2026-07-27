"""
Utilities for constructing one-dimensional photonic-crystal layer sequences.
"""

# ============================================================================
# Construct layer sequences
# ============================================================================

def perfect_crystal_layers(
    n_1: float,
    d_1: float,
    n_2: float,
    d_2: float,
    number_of_mirror_cells: int
) -> list[tuple[float, float]]:
    """
    Construct a finite perfect crystal containing twice the specified
    number of unit cells.

    The structure is

        (AB)^(2N).
    """

    if number_of_mirror_cells <= 0:
        raise ValueError(
            "number_of_mirror_cells must be positive."
        )

    return [
        (n_1, d_1),
        (n_2, d_2),
    ] * (2 * number_of_mirror_cells)


def defect_crystal_layers(
    n_1: float,
    d_1: float,
    n_2: float,
    d_2: float,
    number_of_mirror_cells: int,
    defect_index: float,
    defect_thickness: float
) -> list[tuple[float, float]]:
    """
    Construct the symmetric defect structure

        (AB)^N D (BA)^N.
    """

    if number_of_mirror_cells <= 0:
        raise ValueError(
            "number_of_mirror_cells must be positive."
        )

    if defect_index <= 0.0:
        raise ValueError(
            "defect_index must be positive."
        )

    if defect_thickness <= 0.0:
        raise ValueError(
            "defect_thickness must be positive."
        )

    left_mirror = [
                      (n_1, d_1),
                      (n_2, d_2),
                  ] * number_of_mirror_cells

    right_mirror = [
                       (n_2, d_2),
                       (n_1, d_1),
                   ] * number_of_mirror_cells

    return (
            left_mirror
            + [(defect_index, defect_thickness)]
            + right_mirror
    )
