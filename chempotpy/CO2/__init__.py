"""Geometry ordering helpers for the CO2 five-state RKHS family."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("O", "C", "O")
SURFACES = (
    "CO2_1Ap_RKHS",
    "CO2_1App_RKHS",
    "CO2_21Ap_RKHS",
    "CO2_3Ap_RKHS",
    "CO2_3App_RKHS",
)


def intro():
    """Print the surfaces currently registered for this system."""

    print("CO2 potential energy surfaces: {}".format(", ".join(SURFACES)))


def intro_detail():
    """Print the detailed surface list for this system."""

    print(
        "Five single-state reactive RKHS surfaces for the 1A', 1A\"), "
        "(2)1A', 3A', and 3A\" states of CO2. Each provides energy and "
        "analytic Cartesian gradients. Reference: J. C. San Vicente Veliz "
        "et al., Phys. Chem. Chem. Phys. 23, 11251-11263 (2021)."
    )


def check(system, surface, geometry):
    """Validate and reorder one geometry into O-C-O order."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown CO2 surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("CO2 requires {} atoms".format(len(ATOM_SYMBOLS)))
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("CO2 requires atoms O, C, O")
    carbon = supplied_symbols.index("C")
    oxygens = [index for index, symbol in enumerate(supplied_symbols) if symbol == "O"]
    ordered_indices = [oxygens[0], carbon, oxygens[1]]
    coordinates = np.asarray(
        [[float(value) for value in geometry[index][1:4]] for index in ordered_indices],
        dtype=float,
    )
    rank = np.asarray(
        [(ATOM_SYMBOLS[index], original) for index, original in enumerate(ordered_indices)],
        dtype=object,
    )
    return coordinates, rank


def reverse_order(energy, gradient, coupling, geometry, rank_ordered):
    """Restore Cartesian derivatives to the caller's atom order."""

    del geometry
    reordered_gradient = np.zeros_like(gradient)
    reordered_coupling = np.zeros_like(coupling)
    for ordered, (_, original) in enumerate(rank_ordered):
        reordered_gradient[:, int(original), :] = gradient[:, ordered, :]
        reordered_coupling[:, :, int(original), :] = coupling[:, :, ordered, :]
    return energy, reordered_gradient, reordered_coupling


def reverse_order_dpem(matrix, matrix_gradient, geometry, rank_ordered):
    """Restore diabatic derivatives to the caller's atom order."""

    del geometry
    reordered_gradient = np.zeros_like(matrix_gradient)
    for ordered, (_, original) in enumerate(rank_ordered):
        reordered_gradient[:, :, int(original), :] = matrix_gradient[:, :, ordered, :]
    return matrix, reordered_gradient
