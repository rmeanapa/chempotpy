"""Shared CHEMPOTPY core with compatibility imports for its Python API."""

import sys

from .python import adapters, contract, registry, runtime, units
from .python.adapters import LegacyAdiabaticF2PyAdapter
from .python.contract import (
    PESCapabilities,
    PESContractError,
    PESMetadata,
    PESRequest,
    PESResult,
    PotentialEnergySurface,
    Representation,
)
from .python.registry import SurfaceSource, discover_surface_sources
from .python.runtime import LegacyRuntime, SurfaceRuntimeError

# Preserve the established chempotpy.core.<module> import paths while keeping
# their implementations together under chempotpy/core/python.
sys.modules[__name__ + ".adapters"] = adapters
sys.modules[__name__ + ".contract"] = contract
sys.modules[__name__ + ".registry"] = registry
sys.modules[__name__ + ".runtime"] = runtime
sys.modules[__name__ + ".units"] = units

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
