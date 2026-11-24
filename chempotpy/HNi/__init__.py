"""Geometry ordering helpers for active H atoms on fixed Ni environments."""

import numpy as np


ATOM_SYMBOLS = ("H",)
SURFACES = ("HNi_EAM5_1995", "HNi_EAM6_1996", "HNi_EAM7_1998")


def intro():
    """Print the surfaces currently registered for this system."""

    print("H/Ni potential energy surfaces: " + ", ".join(SURFACES))


def intro_detail():
    """Print the detailed surface list for this system."""

    print(
        "HNi_EAM5_1995, HNi_EAM6_1996, and HNi_EAM7_1998: one active H "
        "atom above the canonical fixed 698-atom Ni(111) POTLIB environment; "
        "analytic energy and gradient."
    )


def check(system, surface, geometry):
    """Validate the one-active-H geometry expected by the fixed-slab adapter."""

    del system, surface
    if len(geometry) != 1 or str(geometry[0][0]) != "H":
        raise ValueError("HNi requires one active H atom")
    coordinates = np.asarray([[float(value) for value in geometry[0][1:4]]], dtype=float)
    rank = np.asarray([("H", 0)], dtype=object)
    return coordinates, rank


def reverse_order(energy, gradient, coupling, geometry, rank_ordered):
    """Return results unchanged because the active system contains one atom."""

    del geometry, rank_ordered
    return energy, gradient, coupling


def reverse_order_dpem(matrix, matrix_gradient, geometry, rank_ordered):
    """Return diabatic results unchanged for the one-atom ordering."""

    del geometry, rank_ordered
    return matrix, matrix_gradient
