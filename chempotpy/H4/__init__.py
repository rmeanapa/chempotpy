"""Geometry helpers for the global BMKP H4 potential energy surface."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("H", "H", "H", "H")
SURFACES = ("H4_BMKP_2002",)


def intro():
    """Print the available H4 potential energy surfaces."""

    print("H4 potential energy surfaces: H4_BMKP_2002 (energy and gradient)")


def intro_detail():
    """Print the scientific scope of the BMKP H4 surface."""

    print(
        "H4_BMKP_2002: global ground-state H4 surface for H2 + H2 dynamics, "
        "with analytic Cartesian gradients. A. I. Boothroyd, P. G. Martin, "
        "W. J. Keogh, and M. R. Peterson, J. Chem. Phys. 116, 666 (2002)."
    )


def check(system, surface, geometry):
    """Validate an H4 geometry while preserving its equivalent-atom order."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown H4 surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("H4 requires four atoms")
    symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("H4 requires four H atoms")
    coordinates = np.asarray(
        [[float(value) for value in atom[1:4]] for atom in geometry],
        dtype=float,
    )
    rank = np.asarray([(symbol, index) for index, symbol in enumerate(symbols)], dtype=object)
    return coordinates, rank


def reverse_order(energy, gradient, coupling, geometry, rank_ordered):
    """Restore adiabatic derivatives to the caller's atom order."""

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
