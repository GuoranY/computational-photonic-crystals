"""
Project:
    Computational Photonic Crystals

Module:
    P04 - Transfer Matrix and Transmission Spectrum

Description:
    Calculate and visualize the transmission spectrum of a finite
    one-dimensional photonic crystal using the transfer matrix method.
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Refractive indices
# =============================================================================

n_1 = 1.0
n_2 = 3.5


# =============================================================================
# Lattice parameters
# =============================================================================

lattice_constant = 1.0
fill_fraction = 0.5
number_of_cells = 10


# =============================================================================
# Layer thicknesses
# =============================================================================

d_1 = fill_fraction * lattice_constant
d_2 = (1.0 - fill_fraction) * lattice_constant


# =============================================================================
# Surrounding media
# =============================================================================

n_incident = 1.0
n_exit = 1.0


# ============================================================================
# Propagation matrix
# ============================================================================

def propagation_matrix(
    refractive_index: float,
    thickness: float,
    normalized_frequency: float,
) -> np.ndarray:
    """
    Construct a propagation matrix that maps the wave amplitudes
    on the right side of a dielectric layer to those on the left side.

    Parameters
    ----------
    refractive_index:
        Refractive index of the dielectric layer.
    thickness:
        Thickness of the dielectric layer.
    normalized_frequency:
        Normalized frequency omega * a / (2 * pi * c).

    Returns
    -------
    np.ndarray
        The 2 x 2 propagation matrix.
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

    Parameters
    ----------
    refractive_index_left:
        Refractive index on the left side of the interface.
    refractive_index_right:
        Refractive index on the right side of the interface.

    Returns
    -------
    np.ndarray
        The 2 x 2 interface matrix.
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

def unit_cell_matrix(normalized_frequency: float) -> np.ndarray:
    """
    Construct the transfer matrix of one unit cell.

    Parameters
    ----------
    normalized_frequency:
        Normalized frequency omega * a / (2 * pi * c).

    Returns
    -------
    np.ndarray
        The 2 x 2 transfer matrix of one unit cell.
    """
    propagation_1 = propagation_matrix(
        refractive_index=n_1,
        thickness=d_1,
        normalized_frequency=normalized_frequency,
    )

    propagation_2 = propagation_matrix(
        refractive_index=n_2,
        thickness=d_2,
        normalized_frequency=normalized_frequency,
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


# ============================================================================
# Crystal matrix
# ============================================================================

def crystal_matrix(normalized_frequency: float) -> np.ndarray:
    """
    Construct the transfer matrix of the complete finite crystal,
    including the incident and exit interfaces.

    Parameters
    ----------
    normalized_frequency:
        Normalized frequency omega * a / (2 * pi * c).

    Returns
    -------
    np.ndarray
        The 2 x 2 transfer matrix of the complete crystal.
    """
    propagation_1 = propagation_matrix(
        refractive_index=n_1,
        thickness=d_1,
        normalized_frequency=normalized_frequency,
    )

    propagation_2 = propagation_matrix(
        refractive_index=n_2,
        thickness=d_2,
        normalized_frequency=normalized_frequency,
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

    total_matrix = total_matrix @ interface_2_exit

    return total_matrix


# ============================================================================
# Transmission
# ============================================================================

def transmission(normalized_frequency: float) -> float:
    """
    Calculate the power transmission coefficient of the finite crystal.

    Parameters
    ----------
    normalized_frequency:
        Normalized frequency omega * a / (2 * pi * c).

    Returns
    -------
    float
        Power transmission coefficient.
    """
    total_matrix = crystal_matrix(normalized_frequency)

    transmission_amplitude = 1.0 / total_matrix[0, 0]

    return (
        n_exit
        / n_incident
        * np.abs(transmission_amplitude) ** 2
    )


if __name__ == "__main__":
# ============================================================================
# Frequency range
# ============================================================================

    normalized_frequencies = np.linspace(0.0, 1.0, 2000)

    transmissions = np.array(
        [
            transmission(frequency)
            for frequency in normalized_frequencies
        ]
    )


# ============================================================================
# Visualization
# ============================================================================

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(
        normalized_frequencies,
        transmissions,
        linewidth=1.5,
    )

    ax.set_xlabel(r"Normalized frequency $\omega a / (2\pi c)$")
    ax.set_ylabel("Transmission")
    ax.set_title("Transmission Spectrum of a 1D Photonic Crystal")

    ax.set_xlim(normalized_frequencies[0], normalized_frequencies[-1])
    ax.set_ylim(0.0, 1.05)

    ax.grid(alpha=0.3)

    fig.tight_layout()

    fig.savefig(
        "../figures/p04_transmission_spectrum.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
