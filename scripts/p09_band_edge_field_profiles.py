"""
Project:
    Computational Photonic Crystals

Module:
    P09 - One-Dimensional Band-Edge Field Profiles

Description:
    Locate the edges of the first photonic band gap and visualize
    the corresponding electric-field intensity profiles.

    The two band-edge modes form different standing-wave patterns:
    the lower-frequency mode is concentrated mainly in the
    high-dielectric region, whereas the higher-frequency mode is
    concentrated mainly in the low-dielectric region.

    Their different spatial overlap with the dielectric structure
    explains the frequency splitting and the opening of the
    photonic band gap.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

n1 = 1.0
n2 = 3.5

lattice_constant = 1.0
fill_fraction = 0.5

d1 = fill_fraction * lattice_constant
d2 = (1.0 - fill_fraction) * lattice_constant

number_of_cells = 6
points_per_layer = 200

minimum_frequency = 0.0
maximum_frequency = 0.5
number_of_frequencies = 5000


# ------------------------------------------------------------
# Layer transfer matrix
# ------------------------------------------------------------

def layer_matrix(
    refractive_index: float,
    thickness: float,
    normalized_frequency: float
) -> np.ndarray:
    """
    Construct the transfer matrix of a homogeneous dielectric layer.

    The electromagnetic state is represented by

        [E, H]^T

    where E is the electric-field amplitude and H is the normalized
    magnetic-field amplitude.

    Parameters
    ----------
    refractive_index:
        Refractive index of the layer.

    thickness:
        Thickness of the layer.

    normalized_frequency:
        Normalized frequency nu = omega a / (2 pi c).

    Returns
    -------
    np.ndarray
        The 2 x 2 layer transfer matrix.
    """

    phase = (
        2.0
        * np.pi
        * refractive_index
        * normalized_frequency
        * thickness
        / lattice_constant
    )

    matrix = np.array(
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

    return matrix


# ------------------------------------------------------------
# Unit-cell transfer matrix
# ------------------------------------------------------------

def unit_cell_matrix(
    normalized_frequency: float
) -> np.ndarray:
    """
    Construct the transfer matrix of one unit cell.

    Each unit cell contains:

        material 1 followed by material 2.
    """

    material_1_matrix = layer_matrix(
        refractive_index=n1,
        thickness=d1,
        normalized_frequency=normalized_frequency
    )

    material_2_matrix = layer_matrix(
        refractive_index=n2,
        thickness=d2,
        normalized_frequency=normalized_frequency
    )

    return material_2_matrix @ material_1_matrix


# ------------------------------------------------------------
# Bloch trace function
# ------------------------------------------------------------

def bloch_trace(
    normalized_frequency: float
) -> float:
    """
    Evaluate

        F(nu) = 1/2 Tr(M_cell).

    Allowed bands satisfy |F| <= 1, whereas photonic band gaps
    satisfy |F| > 1.
    """

    matrix = unit_cell_matrix(normalized_frequency)

    return float(
        np.real(
            0.5 * np.trace(matrix)
        )
    )


# ------------------------------------------------------------
# Find a root using bisection
# ------------------------------------------------------------

def bisection_root(
    function,
    left_boundary: float,
    right_boundary: float,
    target_value: float,
    tolerance: float = 1e-12,
    maximum_iterations: int = 200
) -> float:
    """
    Find a solution of

        function(x) = target_value

    within a specified interval using the bisection method.
    """

    left_value = (
        function(left_boundary)
        - target_value
    )

    right_value = (
        function(right_boundary)
        - target_value
    )

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
            right_value = midpoint_value
        else:
            left_boundary = midpoint
            left_value = midpoint_value

    return 0.5 * (
        left_boundary
        + right_boundary
    )


# ------------------------------------------------------------
# Locate the first photonic band gap
# ------------------------------------------------------------

def find_first_band_gap(
    frequencies: np.ndarray,
    trace_values: np.ndarray
) -> tuple[float, float]:
    """
    Locate the lower and upper edges of the first photonic band gap.

    The first gap is identified as the first continuous frequency
    interval satisfying

        |F(nu)| > 1.
    """

    gap_mask = np.abs(trace_values) > 1.0

    transitions = np.where(
        np.diff(gap_mask.astype(int)) != 0
    )[0]

    if len(transitions) < 2:
        raise RuntimeError(
            "No complete photonic band gap was found in the "
            "selected frequency range."
        )

    lower_transition = transitions[0]
    upper_transition = transitions[1]

    lower_left = frequencies[lower_transition]
    lower_right = frequencies[lower_transition + 1]

    upper_left = frequencies[upper_transition]
    upper_right = frequencies[upper_transition + 1]

    lower_target = np.sign(
        trace_values[lower_transition + 1]
    )

    upper_target = np.sign(
        trace_values[upper_transition]
    )

    lower_edge = bisection_root(
        function=bloch_trace,
        left_boundary=lower_left,
        right_boundary=lower_right,
        target_value=lower_target
    )

    upper_edge = bisection_root(
        function=bloch_trace,
        left_boundary=upper_left,
        right_boundary=upper_right,
        target_value=upper_target
    )

    return lower_edge, upper_edge


# ------------------------------------------------------------
# Obtain the band-edge Bloch state
# ------------------------------------------------------------

def band_edge_state(
    normalized_frequency: float
) -> np.ndarray:
    """
    Calculate the electromagnetic state at a band edge.

    At a band edge, the Bloch eigenvalue is either +1 or -1:

        M_cell v = lambda v.

    The state vector v is found from the approximate null space of

        M_cell - lambda I.
    """

    matrix = unit_cell_matrix(normalized_frequency)

    trace_value = bloch_trace(
        normalized_frequency
    )

    bloch_eigenvalue = (
        1.0
        if trace_value >= 0.0
        else -1.0
    )

    shifted_matrix = (
        matrix
        - bloch_eigenvalue * np.eye(2)
    )

    _, _, conjugate_transpose = np.linalg.svd(
        shifted_matrix
    )

    state = conjugate_transpose.conj().T[:, -1]

    # Remove an arbitrary global complex phase.
    if abs(state[0]) > 1e-14:
        state *= np.exp(
            -1j * np.angle(state[0])
        )

    state /= np.linalg.norm(state)

    return state


# ------------------------------------------------------------
# Field inside one dielectric layer
# ------------------------------------------------------------

def field_inside_layer(
    initial_state: np.ndarray,
    refractive_index: float,
    thickness: float,
    normalized_frequency: float,
    number_of_points: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstruct the electric field inside one homogeneous layer.

    The field is written as

        E(x) = E_plus exp(i q x)
             + E_minus exp(-i q x).
    """

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
        number_of_points,
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

    return local_positions, electric_field


# ------------------------------------------------------------
# Reconstruct the field over several unit cells
# ------------------------------------------------------------

def calculate_field_profile(
    normalized_frequency: float
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray
]:
    """
    Calculate the band-edge field profile across several unit cells.

    Returns
    -------
    positions:
        Spatial coordinates.

    electric_field:
        Complex electric-field values.

    refractive_index_profile:
        Refractive index at each position.

    material_labels:
        Integer labels identifying the two materials.
    """

    current_state = band_edge_state(
        normalized_frequency
    )

    position_segments = []
    field_segments = []
    index_segments = []
    label_segments = []

    current_position = 0.0

    for _ in range(number_of_cells):

        # Material 1
        local_positions, electric_field = (
            field_inside_layer(
                initial_state=current_state,
                refractive_index=n1,
                thickness=d1,
                normalized_frequency=normalized_frequency,
                number_of_points=points_per_layer
            )
        )

        position_segments.append(
            current_position + local_positions
        )

        field_segments.append(electric_field)

        index_segments.append(
            np.full(
                local_positions.shape,
                n1
            )
        )

        label_segments.append(
            np.zeros(
                local_positions.shape,
                dtype=int
            )
        )

        current_state = (
            layer_matrix(
                refractive_index=n1,
                thickness=d1,
                normalized_frequency=normalized_frequency
            )
            @ current_state
        )

        current_position += d1

        # Material 2
        local_positions, electric_field = (
            field_inside_layer(
                initial_state=current_state,
                refractive_index=n2,
                thickness=d2,
                normalized_frequency=normalized_frequency,
                number_of_points=points_per_layer
            )
        )

        position_segments.append(
            current_position + local_positions
        )

        field_segments.append(electric_field)

        index_segments.append(
            np.full(
                local_positions.shape,
                n2
            )
        )

        label_segments.append(
            np.ones(
                local_positions.shape,
                dtype=int
            )
        )

        current_state = (
            layer_matrix(
                refractive_index=n2,
                thickness=d2,
                normalized_frequency=normalized_frequency
            )
            @ current_state
        )

        current_position += d2

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


# ------------------------------------------------------------
# Frequency scan
# ------------------------------------------------------------

frequencies = np.linspace(
    minimum_frequency,
    maximum_frequency,
    number_of_frequencies
)

trace_values = np.array(
    [
        bloch_trace(frequency)
        for frequency in frequencies
    ]
)


# ------------------------------------------------------------
# Find the first gap edges
# ------------------------------------------------------------

lower_edge_frequency, upper_edge_frequency = (
    find_first_band_gap(
        frequencies=frequencies,
        trace_values=trace_values
    )
)


# ------------------------------------------------------------
# Calculate the two band-edge field profiles
# ------------------------------------------------------------

(
    lower_positions,
    lower_field,
    lower_index_profile,
    lower_material_labels
) = calculate_field_profile(
    normalized_frequency=lower_edge_frequency
)

(
    upper_positions,
    upper_field,
    upper_index_profile,
    upper_material_labels
) = calculate_field_profile(
    normalized_frequency=upper_edge_frequency
)


# ------------------------------------------------------------
# Normalize the field intensities
# ------------------------------------------------------------

lower_intensity = np.abs(lower_field) ** 2
upper_intensity = np.abs(upper_field) ** 2

lower_intensity /= np.max(lower_intensity)
upper_intensity /= np.max(upper_intensity)


# ------------------------------------------------------------
# Compare field localization
# ------------------------------------------------------------

lower_intensity_in_material_1 = np.mean(
    lower_intensity[
        lower_material_labels == 0
    ]
)

lower_intensity_in_material_2 = np.mean(
    lower_intensity[
        lower_material_labels == 1
    ]
)

upper_intensity_in_material_1 = np.mean(
    upper_intensity[
        upper_material_labels == 0
    ]
)

upper_intensity_in_material_2 = np.mean(
    upper_intensity[
        upper_material_labels == 1
    ]
)


print(
    "First photonic band gap:"
)

print(
    f"Lower edge frequency = "
    f"{lower_edge_frequency:.6f}"
)

print(
    f"Upper edge frequency = "
    f"{upper_edge_frequency:.6f}"
)

print()

print(
    "Lower-edge mode:"
)

print(
    f"Mean intensity in n = {n1:.1f} region: "
    f"{lower_intensity_in_material_1:.6f}"
)

print(
    f"Mean intensity in n = {n2:.1f} region: "
    f"{lower_intensity_in_material_2:.6f}"
)

print()

print(
    "Upper-edge mode:"
)

print(
    f"Mean intensity in n = {n1:.1f} region: "
    f"{upper_intensity_in_material_1:.6f}"
)

print(
    f"Mean intensity in n = {n2:.1f} region: "
    f"{upper_intensity_in_material_2:.6f}"
)


# ------------------------------------------------------------
# Plotting function for the spatial field profiles
# ------------------------------------------------------------

def plot_field_profile(
    axis,
    positions: np.ndarray,
    intensity: np.ndarray,
    refractive_index_profile: np.ndarray,
    frequency: float,
    title: str
) -> None:
    """
    Plot a normalized electric-field intensity profile together
    with the refractive-index distribution.
    """

    axis.plot(
        positions,
        intensity,
        linewidth=2.0,
        label=r"Normalized $|E(x)|^2$"
    )

    axis.set_ylabel(
        r"Normalized $|E(x)|^2$"
    )

    axis.set_ylim(
        0.0,
        1.12
    )

    axis.grid(alpha=0.25)

    index_axis = axis.twinx()

    index_axis.step(
        positions,
        refractive_index_profile,
        where="post",
        linewidth=1.2,
        alpha=0.55,
        label=r"$n(x)$"
    )

    index_axis.set_ylabel(
        r"Refractive index $n(x)$"
    )

    index_axis.set_ylim(
        0.0,
        1.25 * n2
    )

    axis.set_title(
        title
        + "\n"
        + rf"$\nu = {frequency:.6f}$"
    )

    field_lines, field_labels = (
        axis.get_legend_handles_labels()
    )

    index_lines, index_labels = (
        index_axis.get_legend_handles_labels()
    )

    axis.legend(
        field_lines + index_lines,
        field_labels + index_labels,
        loc="upper right"
    )


# ------------------------------------------------------------
# Create figure
# ------------------------------------------------------------

fig, axes = plt.subplots(
    nrows=3,
    ncols=1,
    figsize=(12, 11)
)


# ------------------------------------------------------------
# Top panel: gap opening in frequency space
# ------------------------------------------------------------

axes[0].plot(
    frequencies,
    trace_values,
    linewidth=1.8,
    label=r"$F(\nu)=\frac{1}{2}\mathrm{Tr}(M_{\mathrm{cell}})$"
)

axes[0].axhline(
    1.0,
    linestyle="--",
    linewidth=1.0
)

axes[0].axhline(
    -1.0,
    linestyle="--",
    linewidth=1.0
)

axes[0].axvspan(
    lower_edge_frequency,
    upper_edge_frequency,
    alpha=0.2,
    label="First photonic band gap"
)

axes[0].axvline(
    lower_edge_frequency,
    linestyle=":",
    linewidth=1.5,
    label="Lower band edge"
)

axes[0].axvline(
    upper_edge_frequency,
    linestyle="-.",
    linewidth=1.5,
    label="Upper band edge"
)

axes[0].set_title(
    "First Photonic Band Gap"
)

axes[0].set_xlabel(
    r"Normalized frequency $\nu=\omega a/(2\pi c)$"
)

axes[0].set_ylabel(
    r"$F(\nu)$"
)

axes[0].grid(alpha=0.25)
axes[0].legend()


# ------------------------------------------------------------
# Middle panel: lower-frequency band-edge mode
# ------------------------------------------------------------

plot_field_profile(
    axis=axes[1],
    positions=lower_positions,
    intensity=lower_intensity,
    refractive_index_profile=lower_index_profile,
    frequency=lower_edge_frequency,
    title=(
        "Lower Band-Edge Mode: "
        "Field Concentrated in the High-Index Regions"
    )
)


# ------------------------------------------------------------
# Bottom panel: upper-frequency band-edge mode
# ------------------------------------------------------------

plot_field_profile(
    axis=axes[2],
    positions=upper_positions,
    intensity=upper_intensity,
    refractive_index_profile=upper_index_profile,
    frequency=upper_edge_frequency,
    title=(
        "Upper Band-Edge Mode: "
        "Field Concentrated in the Low-Index Regions"
    )
)

axes[2].set_xlabel(
    r"Position $x/a$"
)


# ------------------------------------------------------------
# Figure formatting and output
# ------------------------------------------------------------

fig.suptitle(
    "Band-Edge Field Profiles and Photonic Band-Gap Opening",
    fontsize=15
)

fig.tight_layout(
    rect=[0.0, 0.0, 1.0, 0.97]
)


project_root = Path(__file__).resolve().parents[1]

output_directory = project_root / "figures"

output_directory.mkdir(
    parents=True,
    exist_ok=True
)

output_path = (
    output_directory
    / "p09_band_edge_field_profiles.png"
)

fig.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
