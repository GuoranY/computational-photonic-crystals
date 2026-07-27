"""
Utilities for reconstructing electric-field profiles in one-dimensional
multilayer dielectric structures.

This module provides functions for calculating the electric field inside
individual homogeneous layers and across complete finite multilayer systems.

The electromagnetic state is represented by

    [E, H]^T,

and the field inside each dielectric layer is reconstructed as a
superposition of forward- and backward-propagating waves.

The utilities are used to visualize band-edge modes, defect modes, and
other spatial field distributions in one-dimensional photonic crystals.
"""

import numpy as np

from utils.field_transfer_utils import field_layer_matrix
from utils.transmission_utils import scattering_amplitudes

# ============================================================================
# Field inside one layer
# ============================================================================

def field_inside_layer(
    initial_state: np.ndarray,
    refractive_index: float,
    thickness: float,
    normalized_frequency: float,
    lattice_constant: float,
    points_per_layer: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct the electric field inside one homogeneous layer.
    """

    if refractive_index <= 0.0:
        raise ValueError(
            "refractive_index must be positive."
        )

    if thickness <= 0.0:
        raise ValueError(
            "thickness must be positive."
        )

    if lattice_constant <= 0.0:
        raise ValueError(
            "lattice_constant must be positive."
        )

    if points_per_layer < 2:
        raise ValueError(
            "points_per_layer must be at least 2."
        )

    initial_electric_field = initial_state[0]
    initial_magnetic_field = initial_state[1]

    forward_amplitude = 0.5 * (
        initial_electric_field
        + initial_magnetic_field / refractive_index
    )

    backward_amplitude = 0.5 * (
        initial_electric_field
        - initial_magnetic_field / refractive_index
    )

    local_positions = np.linspace(
        0.0,
        thickness,
        points_per_layer,
        endpoint=False
    )

    wave_number = (
        2.0
        * np.pi
        * refractive_index
        * normalized_frequency
        / lattice_constant
    )

    electric_field = (
        forward_amplitude
        * np.exp(
            1j * wave_number * local_positions
        )
        + backward_amplitude
        * np.exp(
            -1j * wave_number * local_positions
        )
    )

    return (
        local_positions,
        electric_field
    )


# ===========================================================================
# Field profile through the defect crystal
# ===========================================================================

def calculate_field_profile(
    layers: list[tuple[float, float]],
    normalized_frequency: float,
    lattice_constant: float,
    points_per_layer: int,
    incident_index: float,
    exit_index: float,
    defect_layer_index: int
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    tuple[float, float]
]:
    """
    Reconstruct the electric field throughout a multilayer structure.

    Parameters
    ----------
    layers:
        Dielectric layer sequence represented by

            (refractive_index, thickness).

    normalized_frequency:
        Normalized frequency

            nu = omega a / (2 pi c).

    lattice_constant:
        Lattice constant a.

    points_per_layer:
        Number of spatial sample points inside each layer.

    incident_index:
        Refractive index of the incident medium.

    exit_index:
        Refractive index of the exit medium.

    defect_layer_index:
        Position of the defect layer in the layer sequence.

    Returns
    -------
    positions:
        Spatial coordinates across the multilayer structure.

    electric_field:
        Complex electric-field values.

    refractive_index_profile:
        Refractive index at each spatial point.

    defect_region:
        Start and end positions of the defect layer.
    """

    if not layers:
        raise ValueError(
            "layers cannot be empty."
        )

    if not 0 <= defect_layer_index < len(layers):
        raise ValueError(
            "defect_layer_index is outside the layer sequence."
        )

    reflection_amplitude, _ = scattering_amplitudes(
        layers=layers,
        normalized_frequency=normalized_frequency,
        lattice_constant=lattice_constant,
        incident_index=incident_index,
        exit_index=exit_index
    )

    current_state = np.array(
        [
            1.0 + reflection_amplitude,
            incident_index * (
                1.0 - reflection_amplitude
            )
        ],
        dtype=complex
    )

    position_segments = []
    field_segments = []
    index_segments = []

    current_position = 0.0

    defect_start = None
    defect_end = None

    for layer_index, (
        refractive_index,
        thickness
    ) in enumerate(layers):

        if layer_index == defect_layer_index:
            defect_start = current_position
            defect_end = (
                current_position
                + thickness
            )

        local_positions, layer_field = field_inside_layer(
            initial_state=current_state,
            refractive_index=refractive_index,
            thickness=thickness,
            normalized_frequency=normalized_frequency,
            lattice_constant=lattice_constant,
            points_per_layer=points_per_layer
        )

        position_segments.append(
            current_position
            + local_positions
        )

        field_segments.append(
            layer_field
        )

        index_segments.append(
            np.full(
                local_positions.shape,
                refractive_index
            )
        )

        current_state = (
            field_layer_matrix(
                refractive_index=refractive_index,
                thickness=thickness,
                normalized_frequency=normalized_frequency,
                lattice_constant=lattice_constant
            )
            @ current_state
        )

        current_position += thickness

    if defect_start is None or defect_end is None:
        raise RuntimeError(
            "The defect layer could not be identified."
        )

    return (
        np.concatenate(position_segments),
        np.concatenate(field_segments),
        np.concatenate(index_segments),
        (
            defect_start,
            defect_end
        )
    )


# ============================================================================
# Reconstruct the field over several unit cells
# ============================================================================

def calculate_band_edge_field_profile(
    normalized_frequency: float,
    initial_state: np.ndarray,
    n_1: float,
    d_1: float,
    n_2: float,
    d_2: float,
    number_of_cells: int,
    lattice_constant: float,
    points_per_layer: int
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray
]:
    """
    Reconstruct a band-edge field profile over several periodic
    unit cells.
    """

    if number_of_cells <= 0:
        raise ValueError(
            "number_of_cells must be positive."
        )

    current_state = initial_state.copy()

    position_segments = []
    field_segments = []
    index_segments = []
    label_segments = []

    current_position = 0.0

    materials = [
        (n_1, d_1, 0),
        (n_2, d_2, 1),
    ]

    for _ in range(number_of_cells):
        for refractive_index, thickness, material_label in materials:
            local_positions, layer_field = field_inside_layer(
                initial_state=current_state,
                refractive_index=refractive_index,
                thickness=thickness,
                normalized_frequency=normalized_frequency,
                lattice_constant=lattice_constant,
                points_per_layer=points_per_layer,
            )

            position_segments.append(
                current_position + local_positions
            )

            field_segments.append(
                layer_field
            )

            index_segments.append(
                np.full(
                    local_positions.shape,
                    refractive_index,
                )
            )

            label_segments.append(
                np.full(
                    local_positions.shape,
                    material_label,
                    dtype=int,
                )
            )

            current_state = (
                    field_layer_matrix(
                        refractive_index=refractive_index,
                        thickness=thickness,
                        normalized_frequency=normalized_frequency,
                        lattice_constant=lattice_constant,
                    )
                    @ current_state
            )

            current_position += thickness

    positions = np.concatenate(
        position_segments
    )

    electric_field = np.concatenate(
        field_segments
    )

    refractive_index_profile = np.concatenate(
        index_segments
    )

    material_labels = np.concatenate(
        label_segments
    )

    return (
        positions,
        electric_field,
        refractive_index_profile,
        material_labels
    )
