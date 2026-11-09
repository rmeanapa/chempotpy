"""Tests for the unified runtime used by all CHEMPOTPY surfaces."""

import types
import unittest
from pathlib import Path

import numpy as np

from chempotpy.core.runtime import LegacyRuntime, SurfaceRuntimeError


class RuntimeTests( unittest.TestCase ):
    def setUp( self ) -> None:
        self.geometry = [
            ["H", 1.0, 2.0, 3.0],
            ["Al", 4.0, 5.0, 6.0],
        ]
        self.calls = []

    def _system_module( self ):
        coordinates = np.array([row[1:] for row in self.geometry])

        def check( system_name, surface_name, geometry ):
            return coordinates.copy(), [0, 1]

        def reverse_order( energy, gradient, coupling, geometry, rank ):
            return energy, gradient, coupling

        def reverse_order_dpem( matrix, gradient, geometry, rank ):
            return matrix, gradient

        return types.SimpleNamespace(
            check=check,
            reverse_order=reverse_order,
            reverse_order_dpem=reverse_order_dpem,
            nAlnH=lambda *args: (1, 1),
            nSinOnH=lambda *args: (1, 0, 1),
        )

    def _pes( self, derivative_index ):
        def pes( *args ):
            self.calls.append(args)
            coordinates = np.asarray(args[0])[:2]
            energy = np.array([np.sum(coordinates * coordinates)])
            gradient = np.array([2.0 * coordinates])
            coupling = np.zeros((1, 1, 2, 3))
            return energy, gradient, coupling

        return pes

    def _runtime( self, surface_name, extension, **categories ):
        system_module = self._system_module()
        modules = {
            "chempotpy.Test": system_module,
            "chempotpy.Test.{}".format(surface_name): extension,
        }
        runtime = LegacyRuntime(
            importer=modules.__getitem__,
            data_root=Path("/test/data"),
            **categories
        )
        return runtime, runtime.prepare("Test", surface_name, self.geometry)

    def test_dispatches_every_historical_pes_signature( self ) -> None:
        cases = (
            ("default", {}, 2),
            ("data", {"requires_data_path": ("data",)}, 3),
            ("variable", {"variable_atom_surfaces": ("variable",)}, 3),
            (
                "almhn",
                {"aluminium_hydride_surfaces": ("almhn",)},
                5,
            ),
            (
                "sioh",
                {"silicon_oxygen_hydride_surfaces": ("sioh",)},
                6,
            ),
        )
        for surface_name, categories, argument_count in cases:
            with self.subTest(surface=surface_name):
                self.calls.clear()
                extension = types.SimpleNamespace(pes=self._pes(-1))
                runtime, prepared = self._runtime(
                    surface_name, extension, **categories
                )
                runtime.adiabatic(prepared, 1)
                self.assertEqual(len(self.calls[-1]), argument_count)
                if surface_name in {"variable", "almhn", "sioh"}:
                    self.assertEqual(self.calls[-1][0].shape, (10000, 3))

    def test_numerical_adiabatic_gradient( self ) -> None:
        extension = types.SimpleNamespace(pes=self._pes(-1))
        runtime, prepared = self._runtime("quadratic", extension)
        _, gradient = runtime.numerical_adiabatic(prepared)
        np.testing.assert_allclose(
            gradient,
            2.0 * prepared.coordinates[np.newaxis, :, :],
            rtol=0.0,
            atol=1.0e-10,
        )

    def test_numerical_diabatic_gradient( self ) -> None:
        weights = np.array([[1.0, 2.0], [2.0, 3.0]])

        def dpem( coordinates, derivative_order ):
            value = np.sum(coordinates * coordinates)
            return value * weights, np.zeros((2, 2, 2, 3))

        extension = types.SimpleNamespace(dpem=dpem)
        runtime, prepared = self._runtime("diabatic", extension)
        _, gradient = runtime.numerical_diabatic(prepared)
        expected = (
            weights[:, :, np.newaxis, np.newaxis]
            * 2.0
            * prepared.coordinates[np.newaxis, np.newaxis, :, :]
        )
        np.testing.assert_allclose(
            gradient, expected, rtol=0.0, atol=1.0e-10
        )

    def test_dispatches_legacy_diabatic_workspace( self ) -> None:
        def dpem(
            coordinates,
            derivative_order,
            energy,
            gradient,
            coupling,
            transformation,
        ):
            self.calls.append((
                coordinates,
                derivative_order,
                energy,
                gradient,
                coupling,
                transformation,
            ))
            return (
                np.zeros((3, 3)),
                np.zeros((3, 3, 2, 3)),
            )

        extension = types.SimpleNamespace(dpem=dpem)
        runtime, prepared = self._runtime(
            "workspace",
            extension,
            diabatic_workspace_state_counts={"workspace": 3},
        )
        runtime.diabatic(prepared, 0)
        self.assertEqual(len(self.calls[-1]), 6)
        self.assertEqual(self.calls[-1][2].shape, (3,))
        self.assertEqual(self.calls[-1][3].shape, (3, 2, 3))
        self.assertEqual(self.calls[-1][4].shape, (3, 3, 2, 3))
        self.assertEqual(self.calls[-1][5].shape, (3, 3))

    def test_rejects_missing_entry_point( self ) -> None:
        runtime, prepared = self._runtime("missing", types.SimpleNamespace())
        with self.assertRaises(SurfaceRuntimeError):
            runtime.adiabatic(prepared, 0)
        with self.assertRaises(SurfaceRuntimeError):
            runtime.diabatic(prepared, 0)


if __name__ == "__main__":
    unittest.main()
