"""Geometry ordering helpers for the H3O2m system."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("O", "O", "H", "H", "H")
SURFACES = ("H3O2m_sGDML_2024",)


def intro():
    """Print the surfaces currently registered for this system."""

    print("H3O2m potential energy surfaces: {}".format(", ".join(SURFACES)))


def intro_detail():
    """Print the detailed surface list for this system."""

    print(
        "H3O2m_sGDML_2024: single-state symmetric gradient-domain "
        "machine-learning surface with analytic gradients. M. Arandhara "
        "and S. G. Ramesh, J. Phys. Chem. A 128, 1600-1610 (2024), "
        "doi:10.1021/acs.jpca.3c08027."
    )


def check(system, surface, geometry):
    """Validate and reorder one geometry into the system's canonical order."""

    del system, surface
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("H3O2m requires {} atoms".format(len(ATOM_SYMBOLS)))
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError(
            "H3O2m requires atoms {}".format(", ".join(ATOM_SYMBOLS))
        )
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
