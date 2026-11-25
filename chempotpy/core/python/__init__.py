"""Python implementation of the shared CHEMPOTPY core."""

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
from .catalog import (
    CatalogSurfaceRecord,
    SurfaceCatalogError,
    get_surface_record,
    resolve_atom_symbols,
    surface_catalog,
)
from .ipi import IPIDriver, IPIProtocolError
from .mdi import MDIBridgeError, MDIEngine, initialize_mdi
from .providers import (
    CatalogPotentialEnergySurface,
    StateEvaluation,
    create_surface,
    evaluate_state,
)
from .registry import SurfaceSource, discover_surface_sources
from .runtime import LegacyRuntime, SurfaceRuntimeError

__all__ = [
    "LegacyAdiabaticF2PyAdapter",
    "LegacyRuntime",
    "CatalogPotentialEnergySurface",
    "CatalogSurfaceRecord",
    "IPIDriver",
    "IPIProtocolError",
    "MDIBridgeError",
    "MDIEngine",
    "PESCapabilities",
    "PESContractError",
    "PESMetadata",
    "PESRequest",
    "PESResult",
    "PotentialEnergySurface",
    "Representation",
    "SurfaceSource",
    "SurfaceRuntimeError",
    "StateEvaluation",
    "SurfaceCatalogError",
    "create_surface",
    "discover_surface_sources",
    "evaluate_state",
    "get_surface_record",
    "initialize_mdi",
    "resolve_atom_symbols",
    "surface_catalog",
]
