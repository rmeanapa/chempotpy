"""Geometry ordering and metadata for the rigid-rotor He--NO surfaces."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("He", "N", "O")
SURFACES = ("HeNO_Klos_2000",)


def intro():
    """Print the surfaces currently registered for HeNO."""

    print("HeNO potential energy surfaces: {}".format(", ".join(SURFACES)))


def intro_detail():
    """Print scientific scope and output ordering for the He--NO model."""

    print(
        "HeNO_Klos_2000: rigid-NO RCCSD(T) interaction potential; "
        "states are A-double-prime and A-prime"
    )


def check(system, surface, geometry):
    """Validate and reorder a He--NO geometry to He, N, O atom order."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown HeNO surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("HeNO requires 3 atoms")
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("HeNO requires atoms He, N, O")
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
