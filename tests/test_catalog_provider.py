"""Tests for catalog-wide conversion to the modern PES contract."""

import unittest

import numpy as np

from chempotpy.core.catalog import CatalogSurfaceRecord
from chempotpy.core.contract import PESRequest, Representation
from chempotpy.core.providers import (
    CatalogPotentialEnergySurface,
    evaluate_state,
)
from chempotpy.core.units import (
    EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR,
    EV_TO_HARTREE,
)


class FakeRuntime:
    def prepare(self, system_name, surface_name, geometry):
        self.geometry = geometry
        return geometry

    def adiabatic(self, prepared, derivative):
        energy = np.array([2.0, 3.0])
        gradient = np.full((2, 2, 3), 4.0) if derivative else None
        coupling = np.zeros((2, 2, 2, 3)) if derivative == 2 else None
        return energy, gradient, coupling

    def numerical_adiabatic(self, prepared):
        return np.array([2.0, 3.0]), np.full((2, 2, 3), 5.0)

    def diabatic(self, prepared, derivative):
        matrix = np.array([[2.0, 0.0], [0.0, 3.0]])
        gradient = np.zeros((2, 2, 2, 3))
        gradient[0, 0] = 4.0
        gradient[1, 1] = 6.0
        return matrix, gradient if derivative else None


def record(entry_point="pes", analytic_gradient=True):
    return CatalogSurfaceRecord(
        system_name="H2",
        surface_name="test",
        entry_point=entry_point,
        atom_count=2,
        state_count=2,
        analytic_gradient=analytic_gradient,
        nonadiabatic_coupling=False,
        variable_atom_count=False,
    )


class CatalogProviderTests(unittest.TestCase):
    def setUp(self) -> None:
        self.coordinates = np.array([[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]])

    def test_converts_legacy_energy_and_gradient_to_atomic_units(self) -> None:
        surface = CatalogPotentialEnergySurface(
            record(), ("H", "H"), runtime=FakeRuntime()
        )
        result = surface.evaluate(self.coordinates, PESRequest(gradient=True))
        np.testing.assert_allclose(result.energy, [2.0, 3.0] * np.array(EV_TO_HARTREE))
        np.testing.assert_allclose(
            result.gradient,
            4.0 * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR,
        )

    def test_uses_numerical_fallback_for_energy_only_surface(self) -> None:
        surface = CatalogPotentialEnergySurface(
            record(analytic_gradient=False),
            ("H", "H"),
            runtime=FakeRuntime(),
        )
        result = surface.evaluate(self.coordinates, PESRequest(gradient=True))
        np.testing.assert_allclose(
            result.gradient,
            5.0 * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR,
        )

    def test_diagonalizes_dpem_for_state_bridge(self) -> None:
        surface = CatalogPotentialEnergySurface(
            record(entry_point="dpem"),
            ("H", "H"),
            runtime=FakeRuntime(),
        )
        state = evaluate_state(surface, self.coordinates, state_index=1)
        self.assertAlmostEqual(state.energy, 3.0 * EV_TO_HARTREE)
        np.testing.assert_allclose(
            state.gradient,
            6.0 * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR,
        )


if __name__ == "__main__":
    unittest.main()
