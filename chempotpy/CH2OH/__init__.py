"""Geometry ordering helpers for the CH2OH multistate surfaces."""

from collections import Counter

import numpy as np


ATOM_SYMBOLS = ("O", "C", "H", "H", "H")
SURFACES = (
    "CH2OH_Surfgen_2state_2018",
    "CH2OH_Surfgen_2state_2018_DPEM",
    "CH2OH_Surfgen_3state_2017",
    "CH2OH_Surfgen_3state_2017_DPEM",
)


def intro():
    """Print the surfaces currently registered for this system."""

    print("CH2OH potential energy surfaces: {}".format(", ".join(SURFACES)))


def intro_detail():
    """Print the detailed surface list for this system."""

    print(
        "CH2OH_Surfgen_2state_2018: two-state quasi-diabatic Surfgen model "
        "with energies, gradients, and derivative couplings; a DPEM interface "
        "is also available. C. Xie, C. L. Malbon, D. R. Yarkony, D. Xie, "
        "and H. Guo (2018).\n"
        "CH2OH_Surfgen_3state_2017: three-state quasi-diabatic Surfgen model "
        "with the same adiabatic and diabatic interfaces. C. L. Malbon and "
        "D. R. Yarkony, J. Chem. Phys. 146, 134302 (2017)."
    )


def check(system, surface, geometry):
    """Validate CH2OH and place the hydroxyl hydrogen last in the canonical order."""

    del system, surface
    if len(geometry) != len(ATOM_SYMBOLS):
        raise ValueError("CH2OH requires {} atoms".format(len(ATOM_SYMBOLS)))
    supplied_symbols = tuple(str(atom[0]) for atom in geometry)
    if Counter(supplied_symbols) != Counter(ATOM_SYMBOLS):
        raise ValueError("CH2OH requires one O, one C, and three H atoms")
    oxygen = supplied_symbols.index("O")
    carbon = supplied_symbols.index("C")
    hydrogens = [index for index, symbol in enumerate(supplied_symbols) if symbol == "H"]
    raw_coordinates = np.asarray(
        [[float(value) for value in atom[1:4]] for atom in geometry], dtype=float
    )
    hydroxyl_hydrogen = min(
        hydrogens,
        key=lambda index: np.linalg.norm(raw_coordinates[index]-raw_coordinates[oxygen]),
    )
    carbon_hydrogens = [index for index in hydrogens if index != hydroxyl_hydrogen]
    ordered_indices = [oxygen, carbon, *carbon_hydrogens, hydroxyl_hydrogen]
    coordinates = raw_coordinates[ordered_indices, :]
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
