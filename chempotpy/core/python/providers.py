"""Catalog-wide adapters from historical extensions to the PES contract."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np

from .catalog import (
    CatalogSurfaceRecord,
    SurfaceCatalogError,
    get_surface_record,
    resolve_atom_symbols,
)
from .contract import (
    PESCapabilities,
    PESContractError,
    PESMetadata,
    PESRequest,
    PESResult,
    Representation,
    validate_coordinates,
)
from .units import (
    BOHR_TO_ANGSTROM,
    EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR,
    EV_TO_HARTREE,
)


GRADIENT_AUTO = "auto"
GRADIENT_ANALYTIC = "analytic"
GRADIENT_NUMERICAL = "numerical"
_GRADIENT_STRATEGIES = {
    GRADIENT_AUTO,
    GRADIENT_ANALYTIC,
    GRADIENT_NUMERICAL,
}


@dataclass(frozen=True)
class StateEvaluation:
    """Store one adiabatic state's energy and gradient in atomic units."""

    energy: float
    gradient: np.ndarray


class CatalogPotentialEnergySurface:
    """Expose any installed historical extension through the modern contract."""

    def __init__(
        self,
        record: CatalogSurfaceRecord,
        atom_symbols: Iterable[str],
        gradient_strategy: str = GRADIENT_AUTO,
        runtime=None,
    ) -> None:
        if gradient_strategy not in _GRADIENT_STRATEGIES:
            raise SurfaceCatalogError(
                "unknown gradient strategy {!r}".format(gradient_strategy)
            )
        if gradient_strategy == GRADIENT_ANALYTIC and not record.analytic_gradient:
            raise SurfaceCatalogError(
                "{} does not publish an analytic gradient".format(record.key)
            )
        self.record = record
        self.gradient_strategy = gradient_strategy
        self._runtime_override = runtime
        self._metadata = PESMetadata(
            system_name=record.system_name,
            surface_name=record.surface_name,
            atom_symbols=tuple(atom_symbols),
            state_count=record.state_count,
            capabilities=PESCapabilities(
                gradient=record.entry_point == "pes",
                nonadiabatic_coupling=(
                    record.entry_point == "pes"
                    and record.nonadiabatic_coupling
                ),
                diabatic_matrix=record.entry_point == "dpem",
                diabatic_gradient=record.entry_point == "dpem",
            ),
        )

    @property
    def metadata( self ) -> PESMetadata:
        """Return immutable catalog metadata for this configured surface."""

        return self._metadata

    def evaluate(
        self, coordinates_bohr: np.ndarray, request: PESRequest
    ) -> PESResult:
        """Evaluate a catalog surface and convert its historical public units."""

        coordinates = validate_coordinates(coordinates_bohr, self.metadata)
        geometry = [
            [symbol, *point]
            for symbol, point in zip(
                self.metadata.atom_symbols,
                coordinates * BOHR_TO_ANGSTROM,
            )
        ]
        runtime = self._runtime()
        prepared = runtime.prepare(
            self.record.system_name,
            self.record.surface_name,
            geometry,
        )
        if request.representation is Representation.ADIABATIC:
            result = self._evaluate_adiabatic(runtime, prepared, request)
        else:
            result = self._evaluate_diabatic(runtime, prepared, request)
        result.validate(self.metadata, request)
        return result

    def _runtime( self ):
        if self._runtime_override is not None:
            return self._runtime_override
        import chempotpy

        return chempotpy._runtime

    def _use_analytic_gradient( self ) -> bool:
        if self.gradient_strategy == GRADIENT_NUMERICAL:
            return False
        if self.gradient_strategy == GRADIENT_ANALYTIC:
            return True
        return self.record.analytic_gradient

    def _evaluate_adiabatic( self, runtime, prepared, request ) -> PESResult:
        if self.record.entry_point != "pes":
            raise PESContractError(
                "{} exposes a diabatic matrix, not adiabatic energies".format(
                    self.record.key
                )
            )
        if (
            request.nonadiabatic_coupling
            and not self.record.nonadiabatic_coupling
        ):
            raise PESContractError(
                "{} does not provide nonadiabatic coupling".format(
                    self.record.key
                )
            )
        if request.nonadiabatic_coupling:
            energy, gradient, coupling = runtime.adiabatic(prepared, 2)
        elif request.gradient and self._use_analytic_gradient():
            energy, gradient, coupling = runtime.adiabatic(prepared, 1)
        elif request.gradient:
            energy, gradient = runtime.numerical_adiabatic(prepared)
            coupling = None
        else:
            energy, _, _ = runtime.adiabatic(prepared, 0)
            gradient = None
            coupling = None
        result = PESResult(
            energy=np.asarray(energy, dtype=np.float64) * EV_TO_HARTREE
        )
        if request.gradient:
            result.gradient = (
                np.asarray(gradient, dtype=np.float64)
                * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR
            )
        if request.nonadiabatic_coupling:
            result.nonadiabatic_coupling = (
                np.asarray(coupling, dtype=np.float64) * BOHR_TO_ANGSTROM
            )
        return result

    def _evaluate_diabatic( self, runtime, prepared, request ) -> PESResult:
        if self.record.entry_point != "dpem":
            raise PESContractError(
                "{} exposes adiabatic energies, not a diabatic matrix".format(
                    self.record.key
                )
            )
        if request.nonadiabatic_coupling:
            raise PESContractError(
                "nonadiabatic coupling requires an adiabatic PES entry"
            )
        if request.gradient and self._use_analytic_gradient():
            matrix, gradient = runtime.diabatic(prepared, 1)
        elif request.gradient:
            matrix, gradient = runtime.numerical_diabatic(prepared)
        else:
            matrix, _ = runtime.diabatic(prepared, 0)
            gradient = None
        result = PESResult(
            diabatic_matrix=(
                np.asarray(matrix, dtype=np.float64) * EV_TO_HARTREE
            )
        )
        if request.gradient:
            result.diabatic_gradient = (
                np.asarray(gradient, dtype=np.float64)
                * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR
            )
        return result


def create_surface(
    system_name: str,
    surface_name: str,
    atom_symbols: Optional[Iterable[str]] = None,
    gradient_strategy: str = GRADIENT_AUTO,
    runtime=None,
) -> CatalogPotentialEnergySurface:
    """Create a modern provider for any callable catalog entry."""

    record = get_surface_record(system_name, surface_name)
    symbols = resolve_atom_symbols(record, atom_symbols)
    return CatalogPotentialEnergySurface(
        record,
        symbols,
        gradient_strategy=gradient_strategy,
        runtime=runtime,
    )


def evaluate_state(
    surface: CatalogPotentialEnergySurface,
    coordinates_bohr: np.ndarray,
    state_index: int = 0,
) -> StateEvaluation:
    """Return one adiabatic state from either a PES or a DPEM provider."""

    if not 0 <= state_index < surface.metadata.state_count:
        raise PESContractError(
            "state_index {} is outside 0..{}".format(
                state_index, surface.metadata.state_count - 1
            )
        )
    if surface.metadata.capabilities.diabatic_matrix:
        request = PESRequest(
            representation=Representation.DIABATIC,
            gradient=True,
        )
        result = surface.evaluate(coordinates_bohr, request)
        energies, vectors = np.linalg.eigh(result.diabatic_matrix)
        state_vector = vectors[:, state_index]
        gradient = np.einsum(
            "i,ijab,j->ab",
            state_vector,
            result.diabatic_gradient,
            state_vector,
            optimize=True,
        )
        return StateEvaluation(float(energies[state_index]), gradient)
    result = surface.evaluate(coordinates_bohr, PESRequest(gradient=True))
    return StateEvaluation(
        float(result.energy[state_index]),
        np.asarray(result.gradient[state_index], dtype=np.float64),
    )
