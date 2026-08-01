"""
Utilities for plotting field distributions and refractive-index profiles
in photonic-crystal simulations.
"""

from matplotlib.axes import Axes
from matplotlib.patches import Polygon
import numpy as np


def plot_field_intensity_with_index_profile(
    axis: Axes,
    positions: np.ndarray,
    intensity: np.ndarray,
    refractive_index_profile: np.ndarray,
    frequency: float,
    title: str,
    maximum_refractive_index: float
) -> None:
    """
    Plot a normalized electric-field intensity profile together with
    the corresponding refractive-index distribution.

    Parameters
    ----------
    axis:
        Matplotlib axis used for the electric-field intensity.

    positions:
        Spatial coordinates.

    intensity:
        Normalized electric-field intensity values.

    refractive_index_profile:
        Refractive index at each spatial coordinate.

    frequency:
        Normalized frequency displayed in the plot title.

    title:
        Main title of the subplot.

    maximum_refractive_index:
        Maximum refractive index used to set the secondary-axis limit.
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
        1.25 * maximum_refractive_index
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


def plot_reciprocal_lattice(
    axis: Axes,
    reciprocal_points: np.ndarray,
    reciprocal_vector_1: np.ndarray,
    reciprocal_vector_2: np.ndarray,
    lattice_constant: float
) -> None:
    """
    Plot reciprocal-lattice points together with the two
    reciprocal-lattice basis vectors.
    """

    reciprocal_scale = (
        2.0
        * np.pi
        / lattice_constant
    )

    normalized_points = (
        reciprocal_points
        / reciprocal_scale
    )

    normalized_vector_1 = (
        reciprocal_vector_1
        / reciprocal_scale
    )

    normalized_vector_2 = (
        reciprocal_vector_2
        / reciprocal_scale
    )

    axis.scatter(
        normalized_points[:, 0],
        normalized_points[:, 1],
        s=42,
        zorder=3
    )

    axis.quiver(
        0.0,
        0.0,
        normalized_vector_1[0],
        normalized_vector_1[1],
        angles="xy",
        scale_units="xy",
        scale=1.0,
        width=0.008,
        zorder=4
    )

    axis.quiver(
        0.0,
        0.0,
        normalized_vector_2[0],
        normalized_vector_2[1],
        angles="xy",
        scale_units="xy",
        scale=1.0,
        width=0.008,
        zorder=4
    )

    axis.text(
        normalized_vector_1[0] + 0.08,
        normalized_vector_1[1] - 0.12,
        r"$\mathbf{b}_1$",
        fontsize=12
    )

    axis.text(
        normalized_vector_2[0] + 0.08,
        normalized_vector_2[1] + 0.02,
        r"$\mathbf{b}_2$",
        fontsize=12
    )

    axis.axhline(
        0.0,
        linewidth=0.8,
        alpha=0.5
    )

    axis.axvline(
        0.0,
        linewidth=0.8,
        alpha=0.5
    )

    axis.set_xlabel(
        r"$G_x/(2\pi/a)$"
    )

    axis.set_ylabel(
        r"$G_y/(2\pi/a)$"
    )

    axis.set_title(
        "Reciprocal Lattice"
    )

    axis.set_aspect(
        "equal"
    )

    axis.set_xlim(
        -2.5,
        2.5
    )

    axis.set_ylim(
        -2.5,
        2.5
    )

    axis.grid(
        alpha=0.25
    )


def plot_first_brillouin_zone(
    axis: Axes,
    brillouin_zone_vertices: np.ndarray,
    high_symmetry_points: dict[str, np.ndarray],
    lattice_constant: float
) -> None:
    """
    Plot the first Brillouin zone of a square lattice together
    with the Gamma-X-M-Gamma high-symmetry path.
    """

    normalization_scale = (
        np.pi
        / lattice_constant
    )

    normalized_vertices = (
        brillouin_zone_vertices
        / normalization_scale
    )

    normalized_gamma = (
        high_symmetry_points["Gamma"]
        / normalization_scale
    )

    normalized_x = (
        high_symmetry_points["X"]
        / normalization_scale
    )

    normalized_m = (
        high_symmetry_points["M"]
        / normalization_scale
    )

    brillouin_zone_patch = Polygon(
        normalized_vertices,
        closed=True,
        fill=False,
        linewidth=2.0
    )

    axis.add_patch(
        brillouin_zone_patch
    )

    symmetry_path = np.array(
        [
            normalized_gamma,
            normalized_x,
            normalized_m,
            normalized_gamma
        ]
    )

    axis.plot(
        symmetry_path[:, 0],
        symmetry_path[:, 1],
        marker="o",
        linewidth=2.0,
        markersize=6,
        label=(
            r"$\Gamma\rightarrow X"
            r"\rightarrow M\rightarrow\Gamma$"
        )
    )

    axis.text(
        normalized_gamma[0] - 0.13,
        normalized_gamma[1] - 0.16,
        r"$\Gamma$",
        fontsize=13
    )

    axis.text(
        normalized_x[0] + 0.06,
        normalized_x[1] - 0.08,
        r"$X$",
        fontsize=13
    )

    axis.text(
        normalized_m[0] + 0.06,
        normalized_m[1] + 0.04,
        r"$M$",
        fontsize=13
    )

    axis.axhline(
        0.0,
        linewidth=0.8,
        alpha=0.5
    )

    axis.axvline(
        0.0,
        linewidth=0.8,
        alpha=0.5
    )

    axis.set_xlabel(
        r"$k_x/(\pi/a)$"
    )

    axis.set_ylabel(
        r"$k_y/(\pi/a)$"
    )

    axis.set_title(
        "First Brillouin Zone"
    )

    axis.set_aspect(
        "equal"
    )

    axis.set_xlim(
        -1.25,
        1.35
    )

    axis.set_ylim(
        -1.25,
        1.25
    )

    axis.grid(
        alpha=0.25
    )

    axis.legend(
        loc="lower left"
    )


def plot_fourier_coefficient_maps(
    axes,
    dielectric_coefficients: np.ndarray,
    inverse_dielectric_coefficients: np.ndarray,
    maximum_index: int
) -> None:
    """
    Plot Fourier-coefficient maps for the dielectric function and
    inverse dielectric function.

    Parameters
    ----------
    axes:
        Sequence containing two Matplotlib axes.

    dielectric_coefficients:
        Two-dimensional array containing epsilon_G.

    inverse_dielectric_coefficients:
        Two-dimensional array containing (1 / epsilon)_G.

    maximum_index:
        Maximum reciprocal-lattice index displayed in the maps.
    """

    extent = (
        -maximum_index - 0.5,
        maximum_index + 0.5,
        -maximum_index - 0.5,
        maximum_index + 0.5
    )

    dielectric_image = axes[0].imshow(
        np.real(dielectric_coefficients),
        origin="lower",
        extent=extent,
        interpolation="nearest"
    )

    axes[0].set_title(
        r"Dielectric coefficients $\varepsilon_{\mathbf{G}}$"
    )

    axes[0].set_xlabel(
        r"Reciprocal index $m$"
    )

    axes[0].set_ylabel(
        r"Reciprocal index $n$"
    )

    axes[0].set_xticks(
        np.arange(
            -maximum_index,
            maximum_index + 1
        )
    )

    axes[0].set_yticks(
        np.arange(
            -maximum_index,
            maximum_index + 1
        )
    )

    axes[0].figure.colorbar(
        dielectric_image,
        ax=axes[0],
        label=r"$\varepsilon_{\mathbf{G}}$"
    )

    inverse_dielectric_image = axes[1].imshow(
        np.real(inverse_dielectric_coefficients),
        origin="lower",
        extent=extent,
        interpolation="nearest"
    )

    axes[1].set_title(
        r"Inverse-dielectric coefficients "
        r"$(1/\varepsilon)_{\mathbf{G}}$"
    )

    axes[1].set_xlabel(
        r"Reciprocal index $m$"
    )

    axes[1].set_ylabel(
        r"Reciprocal index $n$"
    )

    axes[1].set_xticks(
        np.arange(
            -maximum_index,
            maximum_index + 1
        )
    )

    axes[1].set_yticks(
        np.arange(
            -maximum_index,
            maximum_index + 1
        )
    )

    axes[1].figure.colorbar(
        inverse_dielectric_image,
        ax=axes[1],
        label=r"$(1/\varepsilon)_{\mathbf{G}}$"
    )

    for axis in axes:
        axis.grid(
            alpha=0.2
        )
