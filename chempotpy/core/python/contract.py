"""Language-neutral data model for modern PES adapters.

The low-level contract uses atomic units and Cartesian arrays. The historical
Python API remains free to present eV and angstrom after calling an adapter.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional, Protocol, Tuple

import numpy as np


class PESContractError( ValueError ):
    """Report a violation of the shared PES array or capability contract."""


class Representation( str, Enum ):
    """Electronic representation requested from a surface."""

    ADIABATIC = "adiabatic"
    DIABATIC = "diabatic"


@dataclass(frozen=True)
class PESCapabilities:
    """Describe the derivative and representation outputs a PES supports."""

    gradient: bool = False
    nonadiabatic_coupling: bool = False
    diabatic_matrix: bool = False
    diabatic_gradient: bool = False


@dataclass(frozen=True)
class PESMetadata:
    """Describe one callable potential-energy surface."""

    system_name: str
    surface_name: str
    atom_symbols: Tuple[str, ...]
    state_count: int = 1
    capabilities: PESCapabilities = PESCapabilities()

    def __post_init__( self ) -> None:
        if not self.system_name:
            raise PESContractError("system_name must not be empty")
        if not self.surface_name:
            raise PESContractError("surface_name must not be empty")
        if not self.atom_symbols:
            raise PESContractError("atom_symbols must contain at least one atom")
        if self.state_count < 1:
            raise PESContractError("state_count must be positive")

    @property
    def atom_count( self ) -> int:
        """Return the number of atoms required by the surface."""

        return len(self.atom_symbols)


@dataclass(frozen=True)
class PESRequest:
    """Select the electronic representation and derivative outputs."""

    representation: Representation = Representation.ADIABATIC
    gradient: bool = False
    nonadiabatic_coupling: bool = False


@dataclass
class PESResult:
    """Store adiabatic or diabatic PES outputs in atomic units."""

    energy: Optional[np.ndarray] = None
    gradient: Optional[np.ndarray] = None
    nonadiabatic_coupling: Optional[np.ndarray] = None
    diabatic_matrix: Optional[np.ndarray] = None
    diabatic_gradient: Optional[np.ndarray] = None

    def validate( self, metadata: PESMetadata, request: PESRequest ) -> None:
        """Validate output presence and shapes against metadata and a request."""

        states = metadata.state_count
        atoms = metadata.atom_count
        capabilities = metadata.capabilities

        if request.representation is Representation.ADIABATIC:
            self.energy = _array_with_shape(self.energy, (states,), "energy")
            if request.gradient:
                if not capabilities.gradient:
                    raise PESContractError("surface does not support gradients")
                self.gradient = _array_with_shape(
                    self.gradient, (states, atoms, 3), "gradient"
                )
            if request.nonadiabatic_coupling:
                if not capabilities.nonadiabatic_coupling:
                    raise PESContractError(
                        "surface does not support nonadiabatic coupling"
                    )
                self.nonadiabatic_coupling = _array_with_shape(
                    self.nonadiabatic_coupling,
                    (states, states, atoms, 3),
                    "nonadiabatic_coupling",
                )
        else:
            if not capabilities.diabatic_matrix:
                raise PESContractError("surface does not support a diabatic matrix")
            self.diabatic_matrix = _array_with_shape(
                self.diabatic_matrix, (states, states), "diabatic_matrix"
            )
            if request.gradient:
                if not capabilities.diabatic_gradient:
                    raise PESContractError(
                        "surface does not support diabatic gradients"
                    )
                self.diabatic_gradient = _array_with_shape(
                    self.diabatic_gradient,
                    (states, states, atoms, 3),
                    "diabatic_gradient",
                )


class PotentialEnergySurface( Protocol ):
    """Protocol implemented by Python, Fortran, or C-backed PES adapters."""

    @property
    def metadata( self ) -> PESMetadata:
        """Return immutable surface metadata."""

    def evaluate(
        self, coordinates_bohr: np.ndarray, request: PESRequest
    ) -> PESResult:
        """Evaluate Cartesian coordinates shaped ``(natoms, 3)`` in bohr."""


def validate_coordinates(
    coordinates_bohr: np.ndarray, metadata: PESMetadata
) -> np.ndarray:
    """Return a floating Cartesian array after validating its contract shape."""

    coordinates = np.asarray(coordinates_bohr, dtype=np.float64)
    expected = (metadata.atom_count, 3)
    if coordinates.shape != expected:
        raise PESContractError(
            f"coordinates have shape {coordinates.shape}; expected {expected}"
        )
    if not np.all(np.isfinite(coordinates)):
        raise PESContractError("coordinates contain non-finite values")
    return coordinates


def _array_with_shape(
    value: Optional[np.ndarray], expected: Tuple[int, ...], name: str
) -> np.ndarray:
    if value is None:
        raise PESContractError(f"{name} was requested but not returned")
    array = np.asarray(value, dtype=np.float64)
    if array.shape != expected:
        raise PESContractError(
            f"{name} has shape {array.shape}; expected {expected}"
        )
    if not np.all(np.isfinite(array)):
        raise PESContractError(f"{name} contains non-finite values")
    return array
