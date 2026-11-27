"""Geometry ordering and metadata for the Klos Kr--NO interaction potential."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("Kr", "N", "O")
SURFACES = ("KrNO_Klos_2009",)


def intro():
    """Print the surfaces currently registered for KrNO."""

    print("KrNO potential energy surfaces: {}".format(", ".join(SURFACES)))


def intro_detail():
    """Print scientific scope and output ordering for the Kr--NO model."""

    print(
        "KrNO_Klos_2009: rigid-NO interaction potential; states are "
        "A-double-prime and A-prime"
    )


def check(system, surface, geometry):
    """Validate and reorder a Kr--NO geometry to Kr, N, O atom order."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown KrNO surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("KrNO requires 3 atoms")
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("KrNO requires atoms Kr, N, O")
    available = {}
    for index, symbol in enumerate(supplied_symbols):
        available.setdefault(symbol, []).append(index)
    ordered_indices = [available[symbol].pop(0) for symbol in ATOM_SYMBOLS]
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
    """Restore derivatives to the caller's atom ordering."""

    del geometry
    reordered_gradient = np.zeros_like(gradient)
    reordered_coupling = np.zeros_like(coupling)
    for ordered, (_, original) in enumerate(rank_ordered):
        reordered_gradient[:, int(original), :] = gradient[:, ordered, :]
        reordered_coupling[:, :, int(original), :] = coupling[:, :, ordered, :]
    return energy, reordered_gradient, reordered_coupling


def reverse_order_dpem(matrix, matrix_gradient, geometry, rank_ordered):
    """Restore diabatic derivatives to the caller's atom ordering."""

    del geometry
    reordered_gradient = np.zeros_like(matrix_gradient)
    for ordered, (_, original) in enumerate(rank_ordered):
        reordered_gradient[:, :, int(original), :] = matrix_gradient[:, :, ordered, :]
    return matrix, reordered_gradient
