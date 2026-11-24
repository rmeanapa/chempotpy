"""Geometry ordering helpers for CH4 above a fixed Ni(111) environment."""

import numpy as np


ATOM_SYMBOLS = ("H", "H", "H", "H", "C")
SURFACES = ("CH4Ni_EAM6_J1_1996",)


def intro():
    """Print the surfaces currently registered for this system."""

    print("CH4/Ni potential energy surfaces: " + ", ".join(SURFACES))


def intro_detail():
    """Print the detailed surface list for this system."""

    print(
        "CH4Ni_EAM6_J1_1996: four active H atoms and one active C atom "
        "above the canonical fixed 397-atom Ni(111) POTLIB environment; "
        "EAM/DIM plus modified-J1 energy with a centered numerical gradient."
    )


def check(system, surface, geometry):
    """Validate and order an active CH4 geometry as H,H,H,H,C."""

    del system, surface
    symbols = [str(atom[0]) for atom in geometry]
    if len(symbols) != 5 or symbols.count("H") != 4 or symbols.count("C") != 1:
        raise ValueError("CH4Ni requires four H atoms and one C atom")
    ordered_indices = [index for index, symbol in enumerate(symbols) if symbol == "H"]
    ordered_indices.extend(index for index, symbol in enumerate(symbols) if symbol == "C")
    coordinates = np.asarray(
        [[float(value) for value in geometry[index][1:4]] for index in ordered_indices],
        dtype=float,
    )
    rank = np.asarray([(symbols[index], index) for index in ordered_indices], dtype=object)
    return coordinates, rank


def reverse_order(energy, gradient, coupling, geometry, rank_ordered):
    """Restore gradient and coupling atom axes to the caller's order."""

    del geometry
    original_indices = np.asarray(rank_ordered[:, 1], dtype=int)
    restored_gradient = np.zeros_like(gradient)
    restored_coupling = np.zeros_like(coupling)
    restored_gradient[:, original_indices, :] = gradient
    restored_coupling[:, :, original_indices, :] = coupling
    return energy, restored_gradient, restored_coupling


def reverse_order_dpem(matrix, matrix_gradient, geometry, rank_ordered):
    """Restore diabatic-gradient atom axes to the caller's order."""

    del geometry
    original_indices = np.asarray(rank_ordered[:, 1], dtype=int)
    restored_gradient = np.zeros_like(matrix_gradient)
    restored_gradient[:, :, original_indices, :] = matrix_gradient
    return matrix, restored_gradient
