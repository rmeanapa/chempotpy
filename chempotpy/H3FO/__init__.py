"""Geometry ordering helpers for the water--hydrogen-fluoride complex."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("O", "H", "H", "H", "F")
SURFACES = ("H3FO_Tew_2022",)


def intro():
    """Print the available H3FO potential energy surfaces."""

    print("H3FO potential energy surfaces: H3FO_Tew_2022 (energy and gradient)")


def intro_detail():
    """Print the scientific scope of the H2O--HF surface."""

    print(
        "H3FO_Tew_2022: full-dimensional H2O--HF surface fitted to "
        "counterpoise-corrected CCSD(T)(F12*)/aug-cc-pVTZ energies, with "
        "analytic gradients."
    )


def check(system, surface, geometry):
    """Validate H3FO and order O, the two water H atoms, the HF H, and F."""

    del system, surface
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("H3FO requires five atoms")
    symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("H3FO requires atoms O, H, H, H, F")
    coordinates = np.asarray(
        [[float(value) for value in atom[1:4]] for atom in geometry],
        dtype=float,
    )
    oxygen = symbols.index("O")
    fluorine = symbols.index("F")
    hydrogens = [index for index, symbol in enumerate(symbols) if symbol == "H"]
    hf_hydrogen = min(
        hydrogens,
        key=lambda index: np.linalg.norm(coordinates[index] - coordinates[fluorine]),
    )
    water_hydrogens = [index for index in hydrogens if index != hf_hydrogen]
    order = [oxygen, *water_hydrogens, hf_hydrogen, fluorine]
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
