"""Unified runtime dispatch for every historical CHEMPOTPY extension module."""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Optional

import numpy as np


NUMERICAL_STEP_ANGSTROM = 1.0e-3
MAX_VARIABLE_ATOMS = 10000


class SurfaceRuntimeError( RuntimeError ):
    """Report an unavailable entry point or invalid runtime configuration."""


@dataclass
class PreparedSurface:
    """Store one validated geometry and its imported CHEMPOTPY modules."""

    system_name: str
    surface_name: str
    geometry: Any
    coordinates: np.ndarray
    rank_ordered: Any
    system_module: Any
    extension_module: Any


class LegacyRuntime:
    """Route all public APIs through one compatibility implementation."""

    def __init__(
        self,
        requires_data_path: Iterable[str] = (),
        variable_atom_surfaces: Iterable[str] = (),
        aluminium_hydride_surfaces: Iterable[str] = (),
        silicon_oxygen_hydride_surfaces: Iterable[str] = (),
        diabatic_workspace_state_counts: Optional[
            Mapping[str, int]
        ] = None,
        importer: Callable[[str], Any] = importlib.import_module,
        data_root: Optional[Path] = None,
    ) -> None:
        self._requires_data_path = frozenset(requires_data_path)
        self._variable_atom_surfaces = frozenset(variable_atom_surfaces)
        self._aluminium_hydride_surfaces = frozenset(
            aluminium_hydride_surfaces
        )
        self._silicon_oxygen_hydride_surfaces = frozenset(
            silicon_oxygen_hydride_surfaces
        )
        self._diabatic_workspace_state_counts = dict(
            diabatic_workspace_state_counts or {}
        )
        self._importer = importer
        self._data_root = data_root or Path(__file__).resolve().parents[2]

    def prepare(
        self, system_name: str, surface_name: str, geometry: Any
    ) -> PreparedSurface:
        """Validate atom ordering and import one compiled surface extension."""

        system_module = self._importer("chempotpy.{}".format(system_name))
        coordinates, rank_ordered = system_module.check(
            system_name, surface_name, geometry
        )
        extension_module = self._importer(
            "chempotpy.{}.{}".format(system_name, surface_name)
        )
        return PreparedSurface(
            system_name=system_name,
            surface_name=surface_name,
            geometry=geometry,
            coordinates=np.asarray(coordinates, dtype=np.float64),
            rank_ordered=rank_ordered,
            system_module=system_module,
            extension_module=extension_module,
        )

    def adiabatic(
        self, prepared: PreparedSurface, derivative_order: int
    ):
        """Evaluate and restore the caller's atom order for adiabatic output."""

        energy, gradient, coupling = self._call_pes(
            prepared, prepared.coordinates, derivative_order
        )
        return prepared.system_module.reverse_order(
            energy,
            gradient,
            coupling,
            prepared.geometry,
            prepared.rank_ordered,
        )

    def numerical_adiabatic(
        self,
        prepared: PreparedSurface,
        step: float = NUMERICAL_STEP_ANGSTROM,
    ):
        """Calculate an eV/angstrom gradient by central finite differences."""

        if step <= 0.0:
            raise ValueError("finite-difference step must be positive")
        energy, _, coupling = self._call_pes(
            prepared, prepared.coordinates, 0
        )
        gradient = np.zeros(
            (len(energy), len(prepared.coordinates), 3), dtype=np.float64
        )
        for atom_index in range(len(prepared.coordinates)):
            for axis in range(3):
                forward = prepared.coordinates.copy()
                backward = prepared.coordinates.copy()
                forward[atom_index, axis] += step
                backward[atom_index, axis] -= step
                energy_forward, _, _ = self._call_pes(prepared, forward, 0)
                energy_backward, _, _ = self._call_pes(prepared, backward, 0)
                gradient[:, atom_index, axis] = (
                    np.asarray(energy_forward) - np.asarray(energy_backward)
                ) / (2.0 * step)
        energy, gradient, _ = prepared.system_module.reverse_order(
            energy,
            gradient,
            coupling,
            prepared.geometry,
            prepared.rank_ordered,
        )
        return energy, gradient

    def diabatic(
        self, prepared: PreparedSurface, derivative_order: int
    ):
        """Evaluate and restore atom order for a diabatic matrix and gradient."""

        matrix, gradient = self._call_dpem(
            prepared, prepared.coordinates, derivative_order
        )
        return prepared.system_module.reverse_order_dpem(
            matrix,
            gradient,
            prepared.geometry,
            prepared.rank_ordered,
        )

    def numerical_diabatic(
        self,
        prepared: PreparedSurface,
        step: float = NUMERICAL_STEP_ANGSTROM,
    ):
        """Calculate a diabatic eV/angstrom gradient by central differences."""

        if step <= 0.0:
            raise ValueError("finite-difference step must be positive")
        matrix, _ = self._call_dpem(prepared, prepared.coordinates, 0)
        state_count = len(matrix)
        gradient = np.zeros(
            (
                state_count,
                state_count,
                len(prepared.coordinates),
                3,
            ),
            dtype=np.float64,
        )
        for atom_index in range(len(prepared.coordinates)):
            for axis in range(3):
                forward = prepared.coordinates.copy()
                backward = prepared.coordinates.copy()
                forward[atom_index, axis] += step
                backward[atom_index, axis] -= step
                matrix_forward, _ = self._call_dpem(prepared, forward, 0)
                matrix_backward, _ = self._call_dpem(prepared, backward, 0)
                gradient[:, :, atom_index, axis] = (
                    np.asarray(matrix_forward) - np.asarray(matrix_backward)
                ) / (2.0 * step)
        return prepared.system_module.reverse_order_dpem(
            matrix,
            gradient,
            prepared.geometry,
            prepared.rank_ordered,
        )

    def _call_pes(
        self,
        prepared: PreparedSurface,
        coordinates: np.ndarray,
        derivative_order: int,
    ):
        """Invoke one of the historical adiabatic wrapper signatures."""

        if not hasattr(prepared.extension_module, "pes"):
            raise SurfaceRuntimeError(
                "surface {} has no adiabatic pes entry point".format(
                    prepared.surface_name
                )
            )
        kernel = prepared.extension_module.pes
        surface_name = prepared.surface_name
        if surface_name in self._requires_data_path:
            return kernel(
                coordinates, derivative_order, str(self._data_root)
            )
        if surface_name in self._variable_atom_surfaces:
            atom_count = len(prepared.geometry)
            return kernel(
                self._expanded_coordinates(coordinates),
                atom_count,
                derivative_order,
            )
        if surface_name in self._aluminium_hydride_surfaces:
            atom_count = len(prepared.geometry)
            hydrogen_count, aluminium_count = (
                prepared.system_module.nAlnH(
                    prepared.system_name,
                    surface_name,
                    prepared.geometry,
                )
            )
            return kernel(
                self._expanded_coordinates(coordinates),
                hydrogen_count,
                aluminium_count,
                atom_count,
                derivative_order,
            )
        if surface_name in self._silicon_oxygen_hydride_surfaces:
            atom_count = len(prepared.geometry)
            silicon_count, oxygen_count, hydrogen_count = (
                prepared.system_module.nSinOnH(
                    prepared.system_name,
                    surface_name,
                    prepared.geometry,
                )
            )
            return kernel(
                self._expanded_coordinates(coordinates),
                silicon_count,
                oxygen_count,
                hydrogen_count,
                atom_count,
                derivative_order,
            )
        return kernel(coordinates, derivative_order)

    def _call_dpem(
        self,
        prepared: PreparedSurface,
        coordinates: np.ndarray,
        derivative_order: int,
    ):
        """Invoke a historical diabatic wrapper signature."""

        if not hasattr(prepared.extension_module, "dpem"):
            raise SurfaceRuntimeError(
                "surface {} has no diabatic dpem entry point".format(
                    prepared.surface_name
                )
            )
        kernel = prepared.extension_module.dpem
        if prepared.surface_name in self._requires_data_path:
            return kernel(
                coordinates, derivative_order, str(self._data_root)
            )
        state_count = self._diabatic_workspace_state_counts.get(
            prepared.surface_name
        )
        if state_count is not None:
            atom_count = len(coordinates)
            energy = np.zeros(state_count, dtype=np.float64)
            gradient = np.zeros(
                (state_count, atom_count, 3), dtype=np.float64
            )
            coupling = np.zeros(
                (state_count, state_count, atom_count, 3),
                dtype=np.float64,
            )
            transformation = np.zeros(
                (state_count, state_count), dtype=np.float64
            )
            return kernel(
                coordinates,
                derivative_order,
                energy,
                gradient,
                coupling,
                transformation,
            )
        return kernel(coordinates, derivative_order)

    @staticmethod
    def _expanded_coordinates( coordinates: np.ndarray ) -> np.ndarray:
        """Copy a variable-size geometry into the historical work array."""

        if len(coordinates) > MAX_VARIABLE_ATOMS:
            raise SurfaceRuntimeError(
                "geometry exceeds the historical 10000-atom work array"
            )
        expanded = np.zeros((MAX_VARIABLE_ATOMS, 3), dtype=np.float64)
        expanded[: len(coordinates), :] = coordinates
        return expanded
