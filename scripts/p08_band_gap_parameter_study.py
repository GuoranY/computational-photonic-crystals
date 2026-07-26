"""
Project:
    Computational Photonic Crystals

Module:
    P08 - Band-Gap Parameter Study

Description:
    Investigate how the photonic band gap of a one-dimensional periodic
    dielectric depends on refractive-index contrast, fill fraction, and
    the number of unit cells.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utils.bloch_utils import find_band_gaps
from utils.transmission_utils import transmission_spectrum


# ============================================================================
# Base structural parameters
# ============================================================================

n_1 = 1.0
n_2 = 3.5
lattice_constant = 1.0
fill_fraction = 0.5
number_of_cells = 10

n_incident = 1.0
n_exit = 1.0

d_1 = fill_fraction * lattice_constant
d_2 = (1.0 - fill_fraction) * lattice_constant


# ============================================================================
# Frequency sampling
# ============================================================================

normalized_frequencies = np.linspace(1e-4, 1.0, 4000)


# ============================================================================
# Output directory
# ============================================================================

project_root = Path(__file__).resolve().parents[1]
figure_directory = project_root / "figures"
figure_directory.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Refractive-index contrast study
# ============================================================================

n_2_values = [
    1.5,
    2.0,
    2.5,
    3.5,
    4.5,
]

# ============================================================================
# Fill-fraction study
# ============================================================================

fill_fraction_values = [
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.85,
    0.9,
    0.95,
]

# ============================================================================
# Number-of-cells study
# ============================================================================

number_of_cells_values = [
    2,
    4,
    6,
    10,
    20,
]

if __name__ == "__main__":
    transmission_spectra = []

    first_gap_lower_edges = []
    first_gap_upper_edges = []
    first_gap_widths = []
    relative_gap_widths = []

    for n_2_value in n_2_values:
        transmissions = transmission_spectrum(
            normalized_frequencies=normalized_frequencies,
            n_1=n_1,
            n_2=n_2_value,
            d_1=d_1,
            d_2=d_2,
            lattice_constant=lattice_constant,
            number_of_cells=number_of_cells,
            n_incident=n_incident,
            n_exit=n_exit,
        )

        transmission_spectra.append(transmissions)

        band_gaps = find_band_gaps(
            normalized_frequencies=normalized_frequencies,
            n_1=n_1,
            n_2=n_2_value,
            d_1=d_1,
            d_2=d_2,
            lattice_constant=lattice_constant,
        )

        if not band_gaps:
            raise RuntimeError(
                f"No band gap was found for n_2 = {n_2_value}."
            )

        first_gap_lower, first_gap_upper = band_gaps[0]

        gap_width = first_gap_upper - first_gap_lower
        gap_center = 0.5 * (
            first_gap_lower + first_gap_upper
        )

        first_gap_lower_edges.append(first_gap_lower)
        first_gap_upper_edges.append(first_gap_upper)
        first_gap_widths.append(gap_width)
        relative_gap_widths.append(
            gap_width / gap_center
        )

        print(
            f"n_2 = {n_2_value:.1f}: "
            f"gap = [{first_gap_lower:.4f}, "
            f"{first_gap_upper:.4f}], "
            f"width = {gap_width:.4f}, "
            f"relative width = "
            f"{gap_width / gap_center:.4f}"
        )

    index_contrasts = (
        np.array(n_2_values, dtype=float) / n_1
    )

    # ============================================================================
    # Refractive-index-contrast transmission spectra
    # ============================================================================

    fig, axes = plt.subplots(
        nrows=len(n_2_values),
        ncols=1,
        figsize=(9, 10),
        sharex=True,
        sharey=True,
    )

    for ax, current_n_2, current_transmissions in zip(
            axes,
            n_2_values,
            transmission_spectra,
    ):
        ax.plot(
            normalized_frequencies,
            current_transmissions,
            linewidth=1.3,
        )

        ax.set_ylabel("Transmission")

        ax.text(
            0.95,
            0.78,
            rf"$n_2={current_n_2}$",
            transform=ax.transAxes,
            horizontalalignment="right",
            verticalalignment="top",
        )

        ax.grid(alpha=0.3)

    axes[0].set_title(
        "Effect of Refractive-Index Contrast on the First Band Gap"
    )

    axes[-1].set_xlabel(
        r"Normalized frequency $\omega a/(2\pi c)$"
    )

    axes[-1].set_xlim(0.0, 0.45)
    axes[-1].set_ylim(0.0, 1.05)

    fig.tight_layout()

    output_path = (
        figure_directory
        / "p08_index_contrast_transmission.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ============================================================================
    # First band-gap edges
    # ============================================================================

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        index_contrasts,
        first_gap_lower_edges,
        marker="o",
        label="Lower edge",
    )

    ax.plot(
        index_contrasts,
        first_gap_upper_edges,
        marker="o",
        label="Upper edge",
    )

    ax.fill_between(
        index_contrasts,
        first_gap_lower_edges,
        first_gap_upper_edges,
        alpha=0.2,
    )

    ax.set(
        xlabel=r"Refractive-index contrast $n_2/n_1$",
        ylabel=r"Normalized frequency $\nu$",
        title="First Band-Gap Edges versus Refractive-Index Contrast",
    )

    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()

    output_path = (
        figure_directory
        / "p08_index_contrast_gap_edges.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()


    # ============================================================================
    # First band-gap width
    # ============================================================================

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        index_contrasts,
        first_gap_widths,
        marker="o",
        linewidth=1.5,
    )

    ax.set(
        xlabel=r"Refractive-index contrast $n_2/n_1$",
        ylabel=r"First band-gap width $\Delta\nu$",
        title=(
            "First Band-Gap Width versus "
            "Refractive-Index Contrast"
        ),
    )

    ax.grid(alpha=0.3)

    fig.tight_layout()

    output_path = (
        figure_directory
        / "p08_index_contrast_gap_width.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()


    # ============================================================================
    # Relative band-gap width
    # ============================================================================

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        index_contrasts,
        relative_gap_widths,
        marker="o",
        linewidth=1.5,
    )

    ax.set(
        xlabel=r"Refractive-index contrast $n_2/n_1$",
        ylabel=(
            r"Relative gap width "
            r"$\Delta\nu/\nu_{\mathrm{mid}}$"
        ),
        title=(
            "Relative Band-Gap Width versus "
            "Refractive-Index Contrast"
        ),
    )

    ax.grid(alpha=0.3)

    fig.tight_layout()

    output_path = (
        figure_directory
        / "p08_index_contrast_relative_gap_width.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ============================================================================
    # Fill-fraction calculations
    # ============================================================================

    fill_fraction_transmission_spectra = []

    fill_fraction_lower_edges = []
    fill_fraction_upper_edges = []
    fill_fraction_gap_widths = []
    fill_fraction_relative_gap_widths = []

    for current_fill_fraction in fill_fraction_values:
        current_d_1 = (
            current_fill_fraction * lattice_constant
        )

        current_d_2 = (
            (1.0 - current_fill_fraction)
            * lattice_constant
        )

        transmissions = transmission_spectrum(
            normalized_frequencies=normalized_frequencies,
            n_1=n_1,
            n_2=n_2,
            d_1=current_d_1,
            d_2=current_d_2,
            lattice_constant=lattice_constant,
            number_of_cells=number_of_cells,
            n_incident=n_incident,
            n_exit=n_exit,
        )

        fill_fraction_transmission_spectra.append(
            transmissions
        )

        band_gaps = find_band_gaps(
            normalized_frequencies=normalized_frequencies,
            n_1=n_1,
            n_2=n_2,
            d_1=current_d_1,
            d_2=current_d_2,
            lattice_constant=lattice_constant,
        )

        if not band_gaps:
            raise RuntimeError(
                "No band gap was found for "
                f"fill fraction = {current_fill_fraction}."
            )

        first_gap_lower, first_gap_upper = band_gaps[0]

        gap_width = (
                first_gap_upper - first_gap_lower
        )

        gap_center = 0.5 * (
                first_gap_lower + first_gap_upper
        )

        fill_fraction_lower_edges.append(
            first_gap_lower
        )

        fill_fraction_upper_edges.append(
            first_gap_upper
        )

        fill_fraction_gap_widths.append(
            gap_width
        )

        fill_fraction_relative_gap_widths.append(
            gap_width / gap_center
        )

        print(
            f"fill fraction = "
            f"{current_fill_fraction:.2f}: "
            f"gap = [{first_gap_lower:.4f}, "
            f"{first_gap_upper:.4f}], "
            f"width = {gap_width:.4f}, "
            f"relative width = "
            f"{gap_width / gap_center:.4f}"
        )

    # ============================================================================
    # Fill-fraction transmission spectra
    # ============================================================================

    fig, axes = plt.subplots(
        nrows=len(fill_fraction_values),
        ncols=1,
        figsize=(9, 13),
        sharex=True,
        sharey=True,
    )

    for (
            ax,
            current_fill_fraction,
            current_transmissions,
    ) in zip(
        axes,
        fill_fraction_values,
        fill_fraction_transmission_spectra,
    ):
        ax.plot(
            normalized_frequencies,
            current_transmissions,
            linewidth=1.3,
        )

        ax.set_ylabel("Transmission")

        ax.text(
            0.95,
            0.78,
            rf"$f={current_fill_fraction}$",
            transform=ax.transAxes,
            horizontalalignment="right",
            verticalalignment="top",
        )

        ax.grid(alpha=0.3)

    axes[0].set_title(
        "Effect of Fill Fraction on the First Band Gap"
    )

    axes[-1].set_xlabel(
        r"Normalized frequency $\omega a/(2\pi c)$"
    )

    axes[-1].set_xlim(0.0, 0.55)
    axes[-1].set_ylim(0.0, 1.05)

    fig.tight_layout()

    output_path = (
            figure_directory
            / "p08_fill_fraction_transmission.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ============================================================================
    # Fill-fraction band-gap edges
    # ============================================================================

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        fill_fraction_values,
        fill_fraction_lower_edges,
        marker="o",
        label="Lower edge",
    )

    ax.plot(
        fill_fraction_values,
        fill_fraction_upper_edges,
        marker="o",
        label="Upper edge",
    )

    ax.fill_between(
        fill_fraction_values,
        fill_fraction_lower_edges,
        fill_fraction_upper_edges,
        alpha=0.2,
    )

    ax.set(
        xlabel=r"Fill fraction $f=d_1/a$",
        ylabel=r"Normalized frequency $\nu$",
        title=(
            "First Band-Gap Edges versus "
            "Fill Fraction"
        ),
    )

    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()

    output_path = (
            figure_directory
            / "p08_fill_fraction_gap_edges.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ============================================================================
    # Fill-fraction band-gap width
    # ============================================================================

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        fill_fraction_values,
        fill_fraction_gap_widths,
        marker="o",
        linewidth=1.5,
    )

    ax.set(
        xlabel=r"Fill fraction $f=d_1/a$",
        ylabel=r"First band-gap width $\Delta\nu$",
        title=(
            "First Band-Gap Width versus "
            "Fill Fraction"
        ),
    )

    ax.grid(alpha=0.3)

    fig.tight_layout()

    output_path = (
            figure_directory
            / "p08_fill_fraction_gap_width.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ============================================================================
    # Fill-fraction relative band-gap width
    # ============================================================================

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(
        fill_fraction_values,
        fill_fraction_relative_gap_widths,
        marker="o",
        linewidth=1.5,
    )

    ax.set(
        xlabel=r"Fill fraction $f=d_1/a$",
        ylabel=(
            r"Relative gap width "
            r"$\Delta\nu/\nu_{\mathrm{mid}}$"
        ),
        title=(
            "Relative Band-Gap Width versus "
            "Fill Fraction"
        ),
    )

    ax.grid(alpha=0.3)

    fig.tight_layout()

    output_path = (
            figure_directory
            / "p08_fill_fraction_relative_gap_width.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ============================================================================
    # Number-of-cells calculations
    # ============================================================================

    reference_band_gaps = find_band_gaps(
        normalized_frequencies=normalized_frequencies,
        n_1=n_1,
        n_2=n_2,
        d_1=d_1,
        d_2=d_2,
        lattice_constant=lattice_constant,
    )

    if not reference_band_gaps:
        raise RuntimeError(
            "No band gap was found for the reference unit cell."
        )

    first_gap_lower, first_gap_upper = reference_band_gaps[0]
    first_gap_center = 0.5 * (
            first_gap_lower + first_gap_upper
    )

    cell_number_transmission_spectra = []
    center_transmissions = []
    minimum_gap_transmissions = []

    for current_number_of_cells in number_of_cells_values:
        transmissions = transmission_spectrum(
            normalized_frequencies=normalized_frequencies,
            n_1=n_1,
            n_2=n_2,
            d_1=d_1,
            d_2=d_2,
            lattice_constant=lattice_constant,
            number_of_cells=current_number_of_cells,
            n_incident=n_incident,
            n_exit=n_exit,
        )

        cell_number_transmission_spectra.append(
            transmissions
        )

        center_index = np.argmin(
            np.abs(
                normalized_frequencies
                - first_gap_center
            )
        )

        gap_mask = (
                (normalized_frequencies >= first_gap_lower)
                & (normalized_frequencies <= first_gap_upper)
        )

        center_transmission = transmissions[center_index]
        minimum_gap_transmission = np.min(
            transmissions[gap_mask]
        )

        center_transmissions.append(
            center_transmission
        )

        minimum_gap_transmissions.append(
            minimum_gap_transmission
        )

        print(
            f"number of cells = "
            f"{current_number_of_cells:2d}: "
            f"center transmission = "
            f"{center_transmission:.6e}, "
            f"minimum gap transmission = "
            f"{minimum_gap_transmission:.6e}"
        )

    # ============================================================================
    # Number-of-cells transmission spectra
    # ============================================================================

    fig, axes = plt.subplots(
        nrows=len(number_of_cells_values),
        ncols=1,
        figsize=(9, 10),
        sharex=True,
        sharey=True,
    )

    for (
            ax,
            current_number_of_cells,
            current_transmissions,
    ) in zip(
        axes,
        number_of_cells_values,
        cell_number_transmission_spectra,
    ):
        ax.plot(
            normalized_frequencies,
            current_transmissions,
            linewidth=1.3,
        )

        ax.axvspan(
            first_gap_lower,
            first_gap_upper,
            alpha=0.15,
        )

        ax.set_ylabel("Transmission")

        ax.text(
            0.95,
            0.78,
            rf"$N={current_number_of_cells}$",
            transform=ax.transAxes,
            horizontalalignment="right",
            verticalalignment="top",
        )

        ax.grid(alpha=0.3)

    axes[0].set_title(
        "Effect of the Number of Unit Cells on Transmission"
    )

    axes[-1].set_xlabel(
        r"Normalized frequency $\omega a/(2\pi c)$"
    )

    axes[-1].set_xlim(0.05, 0.40)
    axes[-1].set_ylim(0.0, 1.05)

    fig.tight_layout()

    output_path = (
        figure_directory
        / "p08_cell_number_transmission.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()

    # ============================================================================
    # Band-gap-center transmission
    # ============================================================================

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.semilogy(
        number_of_cells_values,
        center_transmissions,
        marker="o",
        linewidth=1.5,
    )

    ax.set(
        xlabel="Number of unit cells",
        ylabel="Transmission at band-gap center",
        title=(
            "Band-Gap Suppression versus "
            "Number of Unit Cells"
        ),
    )

    ax.grid(
        alpha=0.3,
        which="both",
    )

    fig.tight_layout()

    output_path = (
            figure_directory
            / "p08_cell_number_gap_suppression.png"
    )

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.show()
