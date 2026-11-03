"""Atomic-unit contract adapters for CH5 potential-energy surfaces."""

from __future__ import annotations

from types import ModuleType
from typing import Optional

from chempotpy.core import (
    LegacyAdiabaticF2PyAdapter,
    PESCapabilities,
    PESMetadata,
)


CH5_J1_METADATA = PESMetadata(
    system_name="CH5",
    surface_name="CH5_GEN_J1_1987",
    atom_symbols=("H", "H", "H", "H", "C", "H"),
    state_count=1,
    capabilities=PESCapabilities(gradient=True),
)


def create_ch5_j1_surface(
    extension: Optional[ModuleType] = None,
) -> LegacyAdiabaticF2PyAdapter:
    """Create the CH5 J1 adapter, importing its f2py extension on demand."""

    if extension is None:
        from . import CH5_GEN_J1_1987 as extension

    return LegacyAdiabaticF2PyAdapter(extension.pes, CH5_J1_METADATA)
