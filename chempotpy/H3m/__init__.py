"""Geometry ordering helpers for the H3- system."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("H", "H", "H")
SURFACES = ("H3m_AS_2004",)


def intro():
    """Print the surfaces currently registered for this system."""

    print("H3m potential energy surfaces: H3m_AS_2004 (energy)")


def intro_detail():
    """Print the detailed surface list for this system."""

    print(
        "H3m_AS_2004: single-state many-body surface with an asymptotic "
        "long-range correction; energy only. A. N. Panda and N. "
        "Sathyamurthy, J. Chem. Phys. 121, 9343 (2004)."
    )


def check(system, surface, geometry):
    """Validate one H3- geometry; the evaluator is permutation-invariant."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown H3m surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("H3m requires three atoms")
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("H3m requires three H atoms")
    coordinates = np.asarray(
        [[float(value) for value in atom[1:4]] for atom in geometry],
        dtype=float,
    )
    rank = np.asarray(
        [(symbol, index) for index, symbol in enumerate(supplied_symbols)],
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
