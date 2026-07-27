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


from utils.bloch_utils import (
    bisection_root,
    field_band_edge_state,
    field_bloch_trace,
    find_band_gaps,
)

from utils.field_profile_utils import calculate_band_edge_field_profile

from utils.plotting_utils import (
    plot_field_intensity_with_index_profile,
)

# ============================================================================
# Parameters
# ============================================================================

n_1 = 1.0
n_2 = 3.5

lattice_constant = 1.0
fill_fraction = 0.5

d_1 = fill_fraction * lattice_constant
d_2 = (1.0 - fill_fraction) * lattice_constant

number_of_cells = 6
points_per_layer = 200

minimum_frequency = 0.0
maximum_frequency = 0.5
number_of_frequencies = 5000

# ============================================================================
# Frequency scan
# ============================================================================

frequencies = np.linspace(
    minimum_frequency,
    maximum_frequency,
    number_of_frequencies
)

trace_values = np.array(
    [
        field_bloch_trace(
            normalized_frequency=frequency,
            n_1=n_1,
            d_1=d_1,
            n_2=n_2,
            d_2=d_2,
            lattice_constant=lattice_constant
        )
        for frequency in frequencies
    ]
)

# ============================================================================
# Find the first gap edges
# ============================================================================

band_gaps = find_band_gaps(
    normalized_frequencies=frequencies,
    n_1=n_1,
    n_2=n_2,
    d_1=d_1,
    d_2=d_2,
    lattice_constant=lattice_constant
)

if len(band_gaps) == 0:
    raise RuntimeError(
        "No photonic band gap was found."
    )

approximate_lower_edge, approximate_upper_edge = band_gaps[0]

frequency_step = frequencies[1] - frequencies[0]

trace_function = lambda frequency: field_bloch_trace(
    normalized_frequency=frequency,
    n_1=n_1,
    d_1=d_1,
    n_2=n_2,
    d_2=d_2,
    lattice_constant=lattice_constant,
)

lower_edge_frequency = bisection_root(
    function=trace_function,
    left_boundary=(
        approximate_lower_edge
        - frequency_step
    ),
    right_boundary=approximate_lower_edge,
    target_value=-1.0,
)

upper_edge_frequency = bisection_root(
    function=trace_function,
    left_boundary=approximate_upper_edge,
    right_boundary=(
        approximate_upper_edge
        + frequency_step
    ),
    target_value=-1.0,
)

# ============================================================================
# Calculate the two band-edge field profiles
# ============================================================================

lower_initial_state = field_band_edge_state(
    normalized_frequency=lower_edge_frequency,
    n_1=n_1,
    d_1=d_1,
    n_2=n_2,
    d_2=d_2,
    lattice_constant=lattice_constant
)

upper_initial_state = field_band_edge_state(
    normalized_frequency=upper_edge_frequency,
    n_1=n_1,
    d_1=d_1,
    n_2=n_2,
    d_2=d_2,
    lattice_constant=lattice_constant
)

(
    lower_positions,
    lower_field,
    lower_index_profile,
    lower_material_labels
) = calculate_band_edge_field_profile(
    normalized_frequency=lower_edge_frequency,
    initial_state=lower_initial_state,
    n_1=n_1,
    d_1=d_1,
    n_2=n_2,
    d_2=d_2,
    number_of_cells=number_of_cells,
    lattice_constant=lattice_constant,
    points_per_layer=points_per_layer
)

(
    upper_positions,
    upper_field,
    upper_index_profile,
    upper_material_labels
) = calculate_band_edge_field_profile(
    normalized_frequency=upper_edge_frequency,
    initial_state=upper_initial_state,
    n_1=n_1,
    d_1=d_1,
    n_2=n_2,
    d_2=d_2,
    number_of_cells=number_of_cells,
    lattice_constant=lattice_constant,
    points_per_layer=points_per_layer
)

# ============================================================================
# Normalize the field intensities
# ============================================================================

lower_intensity = np.abs(lower_field) ** 2
upper_intensity = np.abs(upper_field) ** 2

lower_intensity /= np.max(lower_intensity)
upper_intensity /= np.max(upper_intensity)

# ============================================================================
# Compare field localization
# ============================================================================

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
    f"Mean intensity in n = {n_1:.1f} region: "
    f"{lower_intensity_in_material_1:.6f}"
)

print(
    f"Mean intensity in n = {n_2:.1f} region: "
    f"{lower_intensity_in_material_2:.6f}"
)

print()

print(
    "Upper-edge mode:"
)

print(
    f"Mean intensity in n = {n_1:.1f} region: "
    f"{upper_intensity_in_material_1:.6f}"
)

print(
    f"Mean intensity in n = {n_2:.1f} region: "
    f"{upper_intensity_in_material_2:.6f}"
)

# ============================================================================
# Create figure
# ============================================================================

fig, axes = plt.subplots(
    nrows=3,
    ncols=1,
    figsize=(12, 11)
)

# ============================================================================
# Top panel: gap opening in frequency space
# ============================================================================

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

# ============================================================================
# Middle panel: lower-frequency band-edge mode
# ============================================================================

plot_field_intensity_with_index_profile(
    axis=axes[1],
    positions=lower_positions,
    intensity=lower_intensity,
    refractive_index_profile=lower_index_profile,
    frequency=lower_edge_frequency,
    title=(
        "Lower Band-Edge Mode: "
        "Field Concentrated in the High-Index Regions"
    ),
    maximum_refractive_index=max(n_1, n_2)
)

# ============================================================================
# Bottom panel: upper-frequency band-edge mode
# ============================================================================

plot_field_intensity_with_index_profile(
    axis=axes[2],
    positions=upper_positions,
    intensity=upper_intensity,
    refractive_index_profile=upper_index_profile,
    frequency=upper_edge_frequency,
    title=(
        "Upper Band-Edge Mode: "
        "Field Concentrated in the Low-Index Regions"
    ),
    maximum_refractive_index=max(n_1, n_2)
)

# ============================================================================
# Figure formatting and output
# ============================================================================

fig.suptitle(
    "Band-Edge Field Profiles and Photonic Band-Gap Opening",
    fontsize=15
)

fig.tight_layout(
    rect=(0.0, 0.0, 1.0, 0.97)
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
