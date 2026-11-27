"""Geometry ordering and metadata for the Klos He--BeO interaction potential."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("He", "Be", "O")
SURFACES = ("HeBeO_Klos_2013",)


def intro():
    """Print the surfaces currently registered for HeBeO."""

    print("HeBeO potential energy surfaces: {}".format(", ".join(SURFACES)))


def intro_detail():
    """Print scientific scope and output ordering for the He--BeO model."""

    print(
        "HeBeO_Klos_2013: rigid-BeO CCSD(T)/AV5Z RKHS interaction potential"
    )


def check(system, surface, geometry):
    """Validate and reorder a He--BeO geometry to He, Be, O atom order."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown HeBeO surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("HeBeO requires 3 atoms")
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("HeBeO requires atoms He, Be, O")
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
