"""
Utility functions for identifying band gaps in computed
photonic band structures.
"""

import numpy as np


def find_band_gaps(
    frequencies: np.ndarray,
    tolerance: float = 1.0e-10
) -> list[dict]:
    """
    Find frequency gaps between adjacent photonic bands.

    A band gap exists between band n and band n + 1 when the
    minimum frequency of the upper band is greater than the
    maximum frequency of the lower band over all sampled
    wavevectors.

    Parameters
    ----------
    frequencies:
        Two-dimensional frequency array with shape

            (number_of_k_points, number_of_bands).

        Each column represents one photonic band.

    tolerance:
        Minimum positive gap width required for a gap to be
        retained.

    Returns
    -------
    list[dict]
        Information about each detected band gap. Each dictionary
        contains

        - lower_band
        - upper_band
        - lower_frequency
        - upper_frequency
        - gap_width
        - midgap_frequency
        - gap_to_midgap_ratio
    """

    frequencies = np.asarray(
        frequencies,
        dtype=float
    )

    if frequencies.ndim != 2:
        raise ValueError(
            "frequencies must be a two-dimensional array."
        )

    if frequencies.shape[1] < 2:
        return []

    band_gaps = []

    number_of_bands = frequencies.shape[1]

    for lower_band_index in range(number_of_bands - 1):
        upper_band_index = lower_band_index + 1

        lower_frequency = float(
            np.max(
                frequencies[:, lower_band_index]
            )
        )

        upper_frequency = float(
            np.min(
                frequencies[:, upper_band_index]
            )
        )

        gap_width = (
            upper_frequency
            - lower_frequency
        )

        if gap_width > tolerance:
            midgap_frequency = 0.5 * (
                lower_frequency
                + upper_frequency
            )

            gap_to_midgap_ratio = (
                gap_width
                / midgap_frequency
                if midgap_frequency > 0.0
                else np.nan
            )

            band_gaps.append(
                {
                    "lower_band": lower_band_index + 1,
                    "upper_band": upper_band_index + 1,
                    "lower_frequency": lower_frequency,
                    "upper_frequency": upper_frequency,
                    "gap_width": gap_width,
                    "midgap_frequency": midgap_frequency,
                    "gap_to_midgap_ratio": gap_to_midgap_ratio
                }
            )

    return band_gaps


def find_overlapping_band_gaps(
    first_band_gaps: list[dict],
    second_band_gaps: list[dict],
    tolerance: float = 1.0e-10
) -> list[dict]:
    """
    Find overlapping frequency intervals between two sets
    of photonic band gaps.

    This function can be used to identify complete photonic
    band gaps shared by TE and TM polarizations.

    Parameters
    ----------
    first_band_gaps:
        Band-gap information for the first polarization.

    second_band_gaps:
        Band-gap information for the second polarization.

    tolerance:
        Minimum overlap width required for an interval to be
        retained.

    Returns
    -------
    list[dict]
        Information about every overlapping band gap.
    """

    overlapping_gaps = []

    for first_gap in first_band_gaps:
        for second_gap in second_band_gaps:
            lower_frequency = max(
                first_gap["lower_frequency"],
                second_gap["lower_frequency"]
            )

            upper_frequency = min(
                first_gap["upper_frequency"],
                second_gap["upper_frequency"]
            )

            gap_width = (
                upper_frequency
                - lower_frequency
            )

            if gap_width > tolerance:
                midgap_frequency = 0.5 * (
                    lower_frequency
                    + upper_frequency
                )

                gap_to_midgap_ratio = (
                    gap_width
                    / midgap_frequency
                    if midgap_frequency > 0.0
                    else np.nan
                )

                overlapping_gaps.append(
                    {
                        "first_lower_band":
                            first_gap["lower_band"],

                        "first_upper_band":
                            first_gap["upper_band"],

                        "second_lower_band":
                            second_gap["lower_band"],

                        "second_upper_band":
                            second_gap["upper_band"],

                        "lower_frequency":
                            lower_frequency,

                        "upper_frequency":
                            upper_frequency,

                        "gap_width":
                            gap_width,

                        "midgap_frequency":
                            midgap_frequency,

                        "gap_to_midgap_ratio":
                            gap_to_midgap_ratio
                    }
                )

    return overlapping_gaps


def find_complete_band_gaps(
    te_frequencies: np.ndarray,
    tm_frequencies: np.ndarray,
    tolerance: float = 1.0e-10
) -> tuple[list[dict], list[dict], list[dict]]:
    """
    Find TE, TM, and complete photonic band gaps.

    Parameters
    ----------
    te_frequencies:
        TE frequency array with shape

            (number_of_k_points, number_of_bands).

    tm_frequencies:
        TM frequency array with shape

            (number_of_k_points, number_of_bands).

    tolerance:
        Minimum positive gap or overlap width.

    Returns
    -------
    tuple[list[dict], list[dict], list[dict]]
        TE band gaps, TM band gaps, and complete band gaps.
    """

    te_band_gaps = find_band_gaps(
        frequencies=te_frequencies,
        tolerance=tolerance
    )

    tm_band_gaps = find_band_gaps(
        frequencies=tm_frequencies,
        tolerance=tolerance
    )

    complete_band_gaps = find_overlapping_band_gaps(
        first_band_gaps=te_band_gaps,
        second_band_gaps=tm_band_gaps,
        tolerance=tolerance
    )

    return (
        te_band_gaps,
        tm_band_gaps,
        complete_band_gaps
    )


def print_band_gaps(
    band_gaps: list[dict],
    polarization: str
) -> None:
    """
    Print detected photonic band gaps.
    """

    print()
    print(
        f"{polarization} band gaps"
    )
    print(
        "-" * 72
    )

    if not band_gaps:
        print(
            "No band gaps were found."
        )
        return

    for gap_number, gap in enumerate(
        band_gaps,
        start=1
    ):
        print(
            f"Gap {gap_number}: "
            f"bands {gap['lower_band']}-"
            f"{gap['upper_band']}"
        )

        print(
            "  Frequency interval: "
            f"{gap['lower_frequency']:.6f} "
            f"to {gap['upper_frequency']:.6f}"
        )

        print(
            f"  Gap width: "
            f"{gap['gap_width']:.6f}"
        )

        print(
            f"  Midgap frequency: "
            f"{gap['midgap_frequency']:.6f}"
        )

        print(
            "  Gap-to-midgap ratio: "
            f"{100.0 * gap['gap_to_midgap_ratio']:.3f}%"
        )


def print_complete_band_gaps(
    complete_band_gaps: list[dict]
) -> None:
    """
    Print complete photonic band gaps shared by TE and TM modes.
    """

    print()
    print(
        "Complete TE-TM photonic band gaps"
    )
    print(
        "-" * 72
    )

    if not complete_band_gaps:
        print(
            "No complete photonic band gap was found."
        )
        return

    for gap_number, gap in enumerate(
        complete_band_gaps,
        start=1
    ):
        print(
            f"Complete gap {gap_number}"
        )

        print(
            "  TE bands: "
            f"{gap['first_lower_band']}-"
            f"{gap['first_upper_band']}"
        )

        print(
            "  TM bands: "
            f"{gap['second_lower_band']}-"
            f"{gap['second_upper_band']}"
        )

        print(
            "  Frequency interval: "
            f"{gap['lower_frequency']:.6f} "
            f"to {gap['upper_frequency']:.6f}"
        )

        print(
            f"  Gap width: "
            f"{gap['gap_width']:.6f}"
        )

        print(
            f"  Midgap frequency: "
            f"{gap['midgap_frequency']:.6f}"
        )

        print(
            "  Gap-to-midgap ratio: "
            f"{100.0 * gap['gap_to_midgap_ratio']:.3f}%"
        )
