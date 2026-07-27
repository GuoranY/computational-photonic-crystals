"""
Project:
    Computational Photonic Crystals

Module:
    P10 - One-Dimensional Defect Mode

Description:
    Introduce a central defect layer into a finite one-dimensional
    photonic crystal and calculate the resulting transmission spectrum.

    A localized defect mode appears as a narrow transmission resonance
    inside the first photonic band gap. The corresponding electric-field
    intensity is concentrated near the defect and decays into the
    surrounding periodic mirrors.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utils.bloch_utils import (
    find_band_gaps,
)

from utils.defect_utils import (
    defect_crystal_layers,
    perfect_crystal_layers
)

from utils.transmission_utils import (
    transmission,
)

from utils.field_profile_utils import calculate_field_profile

# ============================================================================
# Parameters
# ============================================================================

n1 = 1.0
n2 = 3.5

incident_index = 1.0
exit_index = 1.0

lattice_constant = 1.0
fill_fraction = 0.5

d1 = fill_fraction * lattice_constant
d2 = (1.0 - fill_fraction) * lattice_constant

number_of_mirror_cells = 6

defect_index = n1
defect_thickness = lattice_constant

minimum_frequency = 0.10
maximum_frequency = 0.32
number_of_frequencies = 8000

points_per_layer = 120

# ============================================================================
# Frequency scan
# ============================================================================

frequencies = np.linspace(
    minimum_frequency,
    maximum_frequency,
    number_of_frequencies
)

band_gaps = find_band_gaps(
    normalized_frequencies=frequencies,
    n_1=n1,
    n_2=n2,
    d_1=d1,
    d_2=d2,
    lattice_constant=lattice_constant
)

if len(band_gaps) == 0:
    raise RuntimeError(
        "No photonic band gap was found."
    )

lower_gap_edge, upper_gap_edge = band_gaps[0]

# ============================================================================
# Build the two structures
# ============================================================================

perfect_layers = perfect_crystal_layers(
    n_1=n1,
    d_1=d1,
    n_2=n2,
    d_2=d2,
    number_of_mirror_cells=number_of_mirror_cells
)

defect_layers = defect_crystal_layers(
    n_1=n1,
    d_1=d1,
    n_2=n2,
    d_2=d2,
    number_of_mirror_cells=number_of_mirror_cells,
    defect_index=defect_index,
    defect_thickness=defect_thickness
)

# ============================================================================
# Calculate transmission spectra
# ============================================================================

perfect_transmission = np.array(
    [
        transmission(
            layers=perfect_layers,
            normalized_frequency=frequency,
            lattice_constant=lattice_constant,
            incident_index=incident_index,
            exit_index=exit_index
        )
        for frequency in frequencies
    ]
)

defect_transmission = np.array(
    [
        transmission(
            layers=defect_layers,
            normalized_frequency=frequency,
            lattice_constant=lattice_constant,
            incident_index=incident_index,
            exit_index=exit_index
        )
        for frequency in frequencies
    ]
)

# ============================================================================
# Locate the defect resonance inside the gap
# ============================================================================

gap_mask = (
    (frequencies > lower_gap_edge)
    & (frequencies < upper_gap_edge)
)

gap_frequencies = frequencies[
    gap_mask
]

gap_transmission = defect_transmission[
    gap_mask
]

resonance_index = np.argmax(
    gap_transmission
)

defect_frequency = gap_frequencies[
    resonance_index
]

defect_peak_transmission = gap_transmission[
    resonance_index
]

# ============================================================================
# Calculate the resonant field profile
# ============================================================================

defect_layer_index = (
    2 * number_of_mirror_cells
)

(
    positions,
    electric_field,
    refractive_index_profile,
    defect_region
) = calculate_field_profile(
    layers=defect_layers,
    normalized_frequency=defect_frequency,
    lattice_constant=lattice_constant,
    points_per_layer=points_per_layer,
    incident_index=incident_index,
    exit_index=exit_index,
    defect_layer_index=defect_layer_index
)

field_intensity = (
    np.abs(electric_field) ** 2
)

normalized_intensity = (
    field_intensity
    / np.max(field_intensity)
)

# ============================================================================
# Print numerical results
# ============================================================================

print("First photonic band gap:")

print(
    f"Lower edge frequency = "
    f"{lower_gap_edge:.6f}"
)

print(
    f"Upper edge frequency = "
    f"{upper_gap_edge:.6f}"
)

print()

print("Defect resonance:")

print(
    f"Defect-mode frequency = "
    f"{defect_frequency:.6f}"
)

print(
    f"Peak transmission = "
    f"{defect_peak_transmission:.6f}"
)

# ============================================================================
# Create figure
# ============================================================================

fig, axes = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(12, 9)
)

# ============================================================================
# Top panel: transmission spectra
# ============================================================================

axes[0].plot(
    frequencies,
    perfect_transmission,
    linewidth=1.8,
    label="Perfect crystal"
)

axes[0].plot(
    frequencies,
    defect_transmission,
    linewidth=1.8,
    label="Crystal with central defect"
)

axes[0].axvspan(
    lower_gap_edge,
    upper_gap_edge,
    alpha=0.18,
    label="First photonic band gap"
)

axes[0].axvline(
    defect_frequency,
    linestyle="--",
    linewidth=1.5,
    label=(
        rf"Defect mode: "
        rf"$\nu={defect_frequency:.6f}$"
    )
)

axes[0].set_title(
    "Defect Resonance inside the First Photonic Band Gap"
)

axes[0].set_xlabel(
    r"Normalized frequency "
    r"$\nu=\omega a/(2\pi c)$"
)

axes[0].set_ylabel(
    "Power transmission"
)

axes[0].set_xlim(
    minimum_frequency,
    maximum_frequency
)

axes[0].set_ylim(
    -0.03,
    1.08
)

axes[0].grid(alpha=0.25)
axes[0].legend()

# ============================================================================
# Bottom panel: localized field profile
# ============================================================================

axes[1].plot(
    positions,
    normalized_intensity,
    linewidth=2.0,
    label=r"Normalized $|E(x)|^2$"
)

axes[1].axvspan(
    defect_region[0],
    defect_region[1],
    alpha=0.25,
    label="Defect layer"
)

axes[1].set_title(
    "Localized Electric-Field Intensity at the Defect Resonance"
)

axes[1].set_xlabel(
    r"Position $x/a$"
)

axes[1].set_ylabel(
    r"Normalized $|E(x)|^2$"
)

axes[1].set_ylim(
    0.0,
    1.10
)

axes[1].grid(alpha=0.25)

index_axis = axes[1].twinx()

index_axis.step(
    positions,
    refractive_index_profile,
    where="post",
    linewidth=1.1,
    alpha=0.50,
    label=r"$n(x)$"
)

index_axis.set_ylabel(
    r"Refractive index $n(x)$"
)

index_axis.set_ylim(
    0.0,
    1.25 * n2
)

field_lines, field_labels = (
    axes[1].get_legend_handles_labels()
)

index_lines, index_labels = (
    index_axis.get_legend_handles_labels()
)

axes[1].legend(
    field_lines + index_lines,
    field_labels + index_labels,
    loc="upper right"
)

# ============================================================================
# Figure formatting and output
# ============================================================================

fig.suptitle(
    "One-Dimensional Photonic-Crystal Defect Mode",
    fontsize=15
)

fig.tight_layout(
    rect=(
        0.0,
        0.0,
        1.0,
        0.97
    )
)

project_root = Path(__file__).resolve().parents[1]

output_directory = (
    project_root
    / "figures"
)

output_directory.mkdir(
    parents=True,
    exist_ok=True
)

output_path = (
    output_directory
    / "p10_1d_defect_mode.png"
)

fig.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight"
)

plt.show()
