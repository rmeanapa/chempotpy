"""Geometry helpers for the CH2+ three-state coupled diabatic potential."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("H", "H", "C")
SURFACES = ("CH2p_3x3_DPEM_2025",)


def intro():
    """Print the available CH2+ potential energy surfaces."""

    print("CH2p potential energy surfaces: CH2p_3x3_DPEM_2025 (DPEM and gradient)")


def intro_detail():
    """Print the scientific scope of the CH2+ coupled diabatic surface."""

    print(
        "CH2p_3x3_DPEM_2025: three-state coupled diabatic neural-network "
        "surface for H + CH+ and C+ + H2 dynamics. O. Roncero and "
        "P. del Mazo-Sevillano, Phys. Chem. Chem. Phys. 27, 14292 (2025), "
        "doi:10.1039/D5CP01718A."
    )


def check(system, surface, geometry):
    """Validate and reorder a CH2+ geometry into the H, H, C fit order."""

    del system
    if surface not in SURFACES:
        raise ValueError("unknown CH2p surface: {}".format(surface))
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("CH2p requires three atoms")
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("CH2p requires two H atoms and one C atom")
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
        [(ATOM_SYMBOLS[index], original) for index, original in enumerate(ordered_indices)],
        dtype=object,
    )
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
