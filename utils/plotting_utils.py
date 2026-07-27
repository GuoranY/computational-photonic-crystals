"""
Utilities for plotting field distributions and refractive-index profiles
in photonic-crystal simulations.
"""

from matplotlib.axes import Axes
import numpy as np

# ============================================================================
# Plotting function for the spatial field profiles
# ============================================================================

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
