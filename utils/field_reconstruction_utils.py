"""
Utility functions for reconstructing
electromagnetic fields from plane-wave
expansion coefficients.
"""


import numpy as np



def reconstruct_tm_field(
    coefficients: np.ndarray,
    reciprocal_vectors: np.ndarray,
    wavevector: np.ndarray,
    x_grid: np.ndarray,
    y_grid: np.ndarray
) -> np.ndarray:
    """
    Reconstruct the TM electric field Ez(x,y)
    from plane-wave expansion coefficients.

    Ez(r)=sum_G E_G exp(i(k+G)·r)

    Parameters
    ----------
    coefficients:
        Plane-wave coefficients of one eigenmode.

    reciprocal_vectors:
        Reciprocal lattice vectors.

    wavevector:
        Bloch wavevector.

    x_grid:
        x coordinates.

    y_grid:
        y coordinates.

    Returns
    -------
    np.ndarray
        Complex electric field Ez(x,y).
    """


    field = np.zeros(
        (
            len(y_grid),
            len(x_grid)
        ),
        dtype=complex
    )


    shifted_vectors = (
        reciprocal_vectors
        + wavevector
    )


    for coefficient, vector in zip(
        coefficients,
        shifted_vectors
    ):

        phase = (
            vector[0]
            * x_grid[np.newaxis, :]
            +
            vector[1]
            * y_grid[:, np.newaxis]
        )


        field += (
            coefficient
            *
            np.exp(
                1j * phase
            )
        )


    return field


def field_confinement_factor(
    intensity: np.ndarray,
    y_grid: np.ndarray,
    confinement_width: float
) -> float:
    """
    Calculate fraction of field energy
    inside the waveguide channel.
    """

    total_energy = np.sum(
        intensity
    )

    channel_mask = (
        np.abs(y_grid)
        <= confinement_width / 2
    )

    channel_energy = np.sum(
        intensity[channel_mask, :]
    )

    return (
        channel_energy
        /
        total_energy
    )
