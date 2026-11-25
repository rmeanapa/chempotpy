"""Shared CHEMPOTPY core with compatibility imports for its Python API."""

import sys

from .python import (
    adapters,
    catalog,
    constants,
    contract,
    ipi,
    mdi,
    providers,
    registry,
    runtime,
    socket_transport,
    units,
)
from .python.catalog import (
    CatalogSurfaceRecord,
    SurfaceCatalogError,
    get_surface_record,
    resolve_atom_symbols,
    surface_catalog,
)
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
from .python.ipi import IPIDriver, IPIProtocolError
from .python.mdi import MDIBridgeError, MDIEngine, initialize_mdi
from .python.providers import (
    CatalogPotentialEnergySurface,
    StateEvaluation,
    create_surface,
    evaluate_state,
)
from .python.socket_transport import (
    PESWireProtocolError,
    RemotePESError,
    SocketPESClient,
    SocketPESServer,
    serve_connection,
)

# Preserve the established chempotpy.core.<module> import paths while keeping
# their implementations together under chempotpy/core/python.
sys.modules[__name__ + ".adapters"] = adapters
sys.modules[__name__ + ".catalog"] = catalog
sys.modules[__name__ + ".constants"] = constants
sys.modules[__name__ + ".contract"] = contract
sys.modules[__name__ + ".ipi"] = ipi
sys.modules[__name__ + ".mdi"] = mdi
sys.modules[__name__ + ".providers"] = providers
sys.modules[__name__ + ".registry"] = registry
sys.modules[__name__ + ".runtime"] = runtime
sys.modules[__name__ + ".socket_transport"] = socket_transport
sys.modules[__name__ + ".units"] = units

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
    "PESWireProtocolError",
    "PotentialEnergySurface",
    "RemotePESError",
    "Representation",
    "SocketPESClient",
    "SocketPESServer",
    "SurfaceSource",
    "StateEvaluation",
    "SurfaceCatalogError",
    "SurfaceRuntimeError",
    "discover_surface_sources",
    "create_surface",
    "evaluate_state",
    "get_surface_record",
    "initialize_mdi",
    "resolve_atom_symbols",
    "serve_connection",
    "surface_catalog",
]
