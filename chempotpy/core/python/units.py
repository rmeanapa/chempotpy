"""Unit definitions used at CHEMPOTPY's public and kernel boundaries."""

from __future__ import annotations

import numpy as np

from .constants import (
    AVOGADRO_CONSTANT_PER_MOLE,
    BOHR_RADIUS_ANGSTROM,
    ELEMENTARY_CHARGE_COULOMB,
    HARTREE_ENERGY_ELECTRONVOLT,
    HARTREE_ENERGY_JOULE,
    PLANCK_CONSTANT_JOULE_SECOND,
    RYDBERG_CONSTANT_PER_METER,
    SPEED_OF_LIGHT_METER_PER_SECOND,
)


JOULE_PER_THERMOCHEMICAL_KILOCALORIE = 4184.0
BOHR_TO_ANGSTROM = BOHR_RADIUS_ANGSTROM
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
HARTREE_TO_EV = HARTREE_ENERGY_ELECTRONVOLT
EV_TO_HARTREE = 1.0 / HARTREE_TO_EV
HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM = HARTREE_TO_EV / BOHR_TO_ANGSTROM
EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR = (
    1.0 / HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM
)
HARTREE_TO_JOULE = HARTREE_ENERGY_JOULE
JOULE_TO_HARTREE = 1.0 / HARTREE_TO_JOULE
HARTREE_TO_KILOJOULE_PER_MOLE = (
    HARTREE_ENERGY_JOULE * AVOGADRO_CONSTANT_PER_MOLE / 1000.0
)
HARTREE_TO_KILOCALORIE_PER_MOLE = (
    HARTREE_ENERGY_JOULE
    * AVOGADRO_CONSTANT_PER_MOLE
    / JOULE_PER_THERMOCHEMICAL_KILOCALORIE
)
ELECTRONVOLT_TO_KILOJOULE_PER_MOLE = (
    ELEMENTARY_CHARGE_COULOMB * AVOGADRO_CONSTANT_PER_MOLE / 1000.0
)
ELECTRONVOLT_TO_KILOCALORIE_PER_MOLE = (
    ELEMENTARY_CHARGE_COULOMB
    * AVOGADRO_CONSTANT_PER_MOLE
    / JOULE_PER_THERMOCHEMICAL_KILOCALORIE
)
ELECTRONVOLT_TO_INVERSE_CENTIMETER = ELEMENTARY_CHARGE_COULOMB / (
    PLANCK_CONSTANT_JOULE_SECOND
    * SPEED_OF_LIGHT_METER_PER_SECOND
    * 100.0
)
HARTREE_TO_INVERSE_CENTIMETER = 2.0 * RYDBERG_CONSTANT_PER_METER / 100.0

# The current Python API uses these rounded values. They are named explicitly
# so compatibility behavior is never mistaken for the canonical contract.
CHEMPOTPY_LEGACY_BOHR_TO_ANGSTROM = 0.529177211
CHEMPOTPY_LEGACY_HARTREE_TO_EV = 27.211386
CHEMPOTPY_LEGACY_HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM = 51.422067


def angstrom_to_bohr( coordinates: np.ndarray ) -> np.ndarray:
    """Convert Cartesian coordinates from angstrom to bohr."""

    return np.asarray(coordinates, dtype=np.float64) * ANGSTROM_TO_BOHR


def bohr_to_angstrom( coordinates: np.ndarray ) -> np.ndarray:
    """Convert Cartesian coordinates from bohr to angstrom."""

    return np.asarray(coordinates, dtype=np.float64) * BOHR_TO_ANGSTROM


def hartree_to_ev( energy: np.ndarray ) -> np.ndarray:
    """Convert energy values from hartree to electronvolts."""

    return np.asarray(energy, dtype=np.float64) * HARTREE_TO_EV


def ev_to_hartree( energy: np.ndarray ) -> np.ndarray:
    """Convert energy values from electronvolts to hartree."""

    return np.asarray(energy, dtype=np.float64) * EV_TO_HARTREE


def hartree_per_bohr_to_ev_per_angstrom( gradient: np.ndarray ) -> np.ndarray:
    """Convert gradients from hartree/bohr to eV/angstrom."""

    return (
        np.asarray(gradient, dtype=np.float64)
        * HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM
    )


def ev_per_angstrom_to_hartree_per_bohr( gradient: np.ndarray ) -> np.ndarray:
    """Convert gradients from eV/angstrom to hartree/bohr."""

    return (
        np.asarray(gradient, dtype=np.float64)
        * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR
    )
