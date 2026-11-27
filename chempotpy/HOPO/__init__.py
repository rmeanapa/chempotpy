"""Geometry ordering helpers for metaphosphorous acid, HOPO."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("P", "O", "O", "H")
SURFACES = ("HOPO_Tew",)


def intro():
    """Print the available HOPO potential energy surfaces."""

    print("HOPO potential energy surfaces: HOPO_Tew (energy and gradient)")


def intro_detail():
    """Print the scientific scope of the available HOPO surface."""

    print(
        "HOPO_Tew: single-state analytic Morse/Gaussian surface for the cis "
        "and trans rotamers of metaphosphorous acid, with analytic gradients."
    )


def check(system, surface, geometry):
    """Validate HOPO and order P, hydroxyl O, terminal O, H."""

    del system, surface
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("HOPO requires four atoms")
    symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("HOPO requires atoms P, O, O, H")
    coordinates = np.asarray(
        [[float(value) for value in atom[1:4]] for atom in geometry],
        dtype=float,
    )
    phosphorus = symbols.index("P")
    hydrogen = symbols.index("H")
    oxygens = [index for index, symbol in enumerate(symbols) if symbol == "O"]
    hydroxyl_oxygen = min(
        oxygens,
        key=lambda index: np.linalg.norm(coordinates[index] - coordinates[hydrogen]),
    )
    terminal_oxygen = next(index for index in oxygens if index != hydroxyl_oxygen)
    order = [phosphorus, hydroxyl_oxygen, terminal_oxygen, hydrogen]
    rank = np.asarray([(symbols[index], index) for index in order], dtype=object)
    return coordinates[order, :], rank


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
