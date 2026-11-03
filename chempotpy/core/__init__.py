"""Shared interfaces for CHEMPOTPY potential-energy surfaces."""

from .adapters import LegacyAdiabaticF2PyAdapter
from .contract import (
    PESCapabilities,
    PESContractError,
    PESMetadata,
    PESRequest,
    PESResult,
    PotentialEnergySurface,
    Representation,
)
from .registry import SurfaceSource, discover_surface_sources
from .runtime import LegacyRuntime, SurfaceRuntimeError

__all__ = [
    "LegacyAdiabaticF2PyAdapter",
    "LegacyRuntime",
    "PESCapabilities",
    "PESContractError",
    "PESMetadata",
    "PESRequest",
    "PESResult",
    "PotentialEnergySurface",
    "Representation",
    "SurfaceSource",
    "SurfaceRuntimeError",
    "discover_surface_sources",
]
