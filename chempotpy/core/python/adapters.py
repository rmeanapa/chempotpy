"""Compatibility adapters from historical f2py wrappers to the PES contract."""

from __future__ import annotations

from typing import Callable, Tuple

import numpy as np

from .contract import (
    PESContractError,
    PESMetadata,
    PESRequest,
    PESResult,
    Representation,
    validate_coordinates,
)
from .units import (
    CHEMPOTPY_LEGACY_BOHR_TO_ANGSTROM,
    CHEMPOTPY_LEGACY_HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM,
    CHEMPOTPY_LEGACY_HARTREE_TO_EV,
)


LegacyAdiabaticKernel = Callable[
    [np.ndarray, int], Tuple[np.ndarray, np.ndarray, np.ndarray]
]


class LegacyAdiabaticF2PyAdapter:
    """Expose a historical adiabatic f2py ``pes`` wrapper in atomic units.

    Historical CHEMPOTPY wrappers accept angstrom and return eV, eV/angstrom,
    and 1/angstrom. This first adapter intentionally uses their rounded
    conversion factors so it reproduces released behavior exactly. New direct
    kernel adapters should use the canonical constants in ``core.units``.
    """

    def __init__(
        self, kernel: LegacyAdiabaticKernel, metadata: PESMetadata
    ) -> None:
        self._kernel = kernel
        self._metadata = metadata

    @property
    def metadata( self ) -> PESMetadata:
        """Return immutable metadata for the wrapped surface."""

        return self._metadata

    def evaluate(
        self, coordinates_bohr: np.ndarray, request: PESRequest
    ) -> PESResult:
        """Evaluate an ordered Cartesian geometry using the atomic-unit contract."""

        coordinates = validate_coordinates(coordinates_bohr, self.metadata)
        if request.representation is not Representation.ADIABATIC:
            raise PESContractError(
                "legacy adiabatic adapter does not provide a diabatic matrix"
            )
        if request.nonadiabatic_coupling:
            raise PESContractError(
                "this compatibility adapter does not provide nonadiabatic coupling"
            )
        if request.gradient and not self.metadata.capabilities.gradient:
            raise PESContractError("surface does not support gradients")

        coordinates_angstrom = (
            coordinates * CHEMPOTPY_LEGACY_BOHR_TO_ANGSTROM
        )
        derivative_order = 1 if request.gradient else 0
        energy_ev, gradient_ev_per_angstrom, _ = self._kernel(
            coordinates_angstrom, derivative_order
        )

        result = PESResult(
            energy=(
                np.asarray(energy_ev, dtype=np.float64)
                / CHEMPOTPY_LEGACY_HARTREE_TO_EV
            )
        )
        if request.gradient:
            result.gradient = (
                np.asarray(gradient_ev_per_angstrom, dtype=np.float64)
                / CHEMPOTPY_LEGACY_HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM
            )
        result.validate(self.metadata, request)
        return result
