"""Geometry ordering helpers for the formic-acid dimer surface."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("H", "H", "H", "H", "O", "O", "O", "O", "C", "C")
SURFACES = ("C2H4O4_FAD_PIP",)


def intro():
    """Print the surfaces currently registered for this system."""

    print("C2H4O4 potential energy surfaces: C2H4O4_FAD_PIP (energy)")


def intro_detail():
    """Print the detailed surface list for this system."""

    print(
        "C2H4O4_FAD_PIP: single-state permutationally invariant polynomial "
        "surface for the formic acid dimer; energy only. C. Qu and J. M. "
        "Bowman, Phys. Chem. Chem. Phys. 18, 24835 (2016)."
    )


def check(system, surface, geometry):
    """Validate and reorder one geometry into H4-O4-C2 order."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown C2H4O4 surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("C2H4O4 requires {} atoms".format(len(ATOM_SYMBOLS)))
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("C2H4O4 requires four H, four O, and two C atoms")
    available = {}
    for index, symbol in enumerate(supplied_symbols):
        available.setdefault(symbol, []).append(index)
    ordered_indices = [available[symbol].pop(0) for symbol in ATOM_SYMBOLS]
    coordinates = np.asarray(
        [
            [float(value) for value in geometry[index][1:4]]
            for index in ordered_indices
        ],
        dtype=float,
    )
    rank = np.asarray(
        [
            (ATOM_SYMBOLS[index], original)
            for index, original in enumerate(ordered_indices)
        ],
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
