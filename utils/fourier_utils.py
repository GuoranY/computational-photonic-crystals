"""
Fourier utilities for two-dimensional photonic crystals.

This module calculates the reciprocal-space Fourier coefficients
of piecewise-constant material functions for a square lattice of
circular dielectric rods.
"""

import numpy as np
from scipy.special import j1


def circular_inclusion_fourier_coefficients(
    reciprocal_vectors: np.ndarray,
    lattice_constant: float,
    radius: float,
    inclusion_value: float,
    background_value: float
) -> np.ndarray:
    """
    Calculate Fourier coefficients for a circular inclusion embedded
    in a uniform background.

    The circular inclusion is assumed to be centered at the origin
    of the unit cell.

    Parameters
    ----------
    reciprocal_vectors:
        Array of reciprocal-lattice vectors with shape (..., 2).

    lattice_constant:
        Real-space lattice constant.

    radius:
        Radius of the circular inclusion.

    inclusion_value:
        Value of the material function inside the circular inclusion.

    background_value:
        Value of the material function outside the inclusion.

    Returns
    -------
    np.ndarray
        Fourier coefficients with shape reciprocal_vectors.shape[:-1].
    """

    reciprocal_vectors = np.asarray(
        reciprocal_vectors,
        dtype=float
    )

    if reciprocal_vectors.shape[-1] != 2:
        raise ValueError(
            "reciprocal_vectors must have shape (..., 2)."
        )

    if lattice_constant <= 0.0:
        raise ValueError(
            "lattice_constant must be positive."
        )

    if radius <= 0.0:
        raise ValueError(
            "radius must be positive."
        )

    if radius >= lattice_constant / 2.0:
        raise ValueError(
            "radius must be smaller than half the lattice constant."
        )

    reciprocal_magnitudes = np.linalg.norm(
        reciprocal_vectors,
        axis=-1
    )

    filling_fraction = (
        np.pi * radius ** 2
        / lattice_constant ** 2
    )

    argument = reciprocal_magnitudes * radius

    form_factor = np.ones_like(
        reciprocal_magnitudes,
        dtype=float
    )

    nonzero_mask = ~np.isclose(
        reciprocal_magnitudes,
        0.0
    )

    form_factor[nonzero_mask] = (
        2.0
        * j1(argument[nonzero_mask])
        / argument[nonzero_mask]
    )

    coefficients = (
        inclusion_value - background_value
    ) * filling_fraction * form_factor

    zero_mask = np.isclose(
        reciprocal_magnitudes,
        0.0
    )

    coefficients[zero_mask] += background_value

    return coefficients


def dielectric_fourier_coefficients(
    reciprocal_vectors: np.ndarray,
    lattice_constant: float,
    radius: float,
    rod_refractive_index: float,
    background_refractive_index: float
) -> np.ndarray:
    """
    Calculate Fourier coefficients of the dielectric function epsilon(r).

    Parameters
    ----------
    reciprocal_vectors:
        Array of reciprocal-lattice vectors with shape (..., 2).

    lattice_constant:
        Real-space lattice constant.

    radius:
        Radius of the dielectric rods.

    rod_refractive_index:
        Refractive index inside the rods.

    background_refractive_index:
        Refractive index of the background.

    Returns
    -------
    np.ndarray
        Fourier coefficients epsilon_G.
    """

    rod_permittivity = rod_refractive_index ** 2
    background_permittivity = background_refractive_index ** 2

    return circular_inclusion_fourier_coefficients(
        reciprocal_vectors=reciprocal_vectors,
        lattice_constant=lattice_constant,
        radius=radius,
        inclusion_value=rod_permittivity,
        background_value=background_permittivity
    )


def inverse_dielectric_fourier_coefficients(
    reciprocal_vectors: np.ndarray,
    lattice_constant: float,
    radius: float,
    rod_refractive_index: float,
    background_refractive_index: float
) -> np.ndarray:
    """
    Calculate Fourier coefficients of the inverse dielectric function
    1 / epsilon(r).

    Parameters
    ----------
    reciprocal_vectors:
        Array of reciprocal-lattice vectors with shape (..., 2).

    lattice_constant:
        Real-space lattice constant.

    radius:
        Radius of the dielectric rods.

    rod_refractive_index:
        Refractive index inside the rods.

    background_refractive_index:
        Refractive index of the background.

    Returns
    -------
    np.ndarray
        Fourier coefficients eta_G = (1 / epsilon)_G.
    """

    rod_inverse_permittivity = (
        1.0 / rod_refractive_index ** 2
    )

    background_inverse_permittivity = (
        1.0 / background_refractive_index ** 2
    )

    return circular_inclusion_fourier_coefficients(
        reciprocal_vectors=reciprocal_vectors,
        lattice_constant=lattice_constant,
        radius=radius,
        inclusion_value=rod_inverse_permittivity,
        background_value=background_inverse_permittivity
    )
