"""Installed metadata catalog for constructing transport-neutral PES providers."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import importlib
import json
from pathlib import Path
import re
from typing import Iterable, Optional, Tuple


CATALOG_SCHEMA_VERSION = 1
_CATALOG_PATH = Path(__file__).resolve().parents[1] / "provider_catalog.json"
_SPECIAL_ATOM_SYMBOLS = {
    "ALC": ("A", "L", "C"),
    "CH4_H2O_H2O": ("C", "H", "H", "H", "H", "O", "H", "H", "O", "H", "H"),
    "COCO": ("C", "O", "C", "O"),
    "HOCO_Ar": ("H", "O", "C", "O", "Ar"),
    "MCH": ("M", "C", "H"),
    "MXH": ("M", "X", "H"),
    "YRH": ("Y", "R", "H"),
    "malonaldehyde": ("C", "C", "O", "C", "O", "H", "H", "H", "H"),
}
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)([0-9]*)")


class SurfaceCatalogError( ValueError ):
    """Report a missing surface or incomplete catalog configuration."""


@dataclass(frozen=True)
class CatalogSurfaceRecord:
    """Describe one callable PES or DPEM entry in the installed catalog."""

    system_name: str
    surface_name: str
    entry_point: str
    atom_count: int
    state_count: int
    analytic_gradient: bool
    nonadiabatic_coupling: bool
    variable_atom_count: bool

    @property
    def key( self ) -> str:
        """Return the stable ``system/surface`` catalog key."""

        return "{}/{}".format(self.system_name, self.surface_name)


@lru_cache(maxsize=1)
def surface_catalog() -> Tuple[CatalogSurfaceRecord, ...]:
    """Load and validate the installed provider catalog once."""

    payload = json.loads(_CATALOG_PATH.read_text())
    if payload.get("schema_version") != CATALOG_SCHEMA_VERSION:
        raise SurfaceCatalogError("unsupported provider catalog schema")
    records = tuple(
        CatalogSurfaceRecord(**record) for record in payload.get("surfaces", ())
    )
    if not records:
        raise SurfaceCatalogError("provider catalog is empty")
    keys = [record.key for record in records]
    if len(keys) != len(set(keys)):
        raise SurfaceCatalogError("provider catalog contains duplicate surfaces")
    return records


@lru_cache(maxsize=None)
def get_surface_record(
    system_name: str, surface_name: str
) -> CatalogSurfaceRecord:
    """Return one catalog record or raise a descriptive error."""

    for record in surface_catalog():
        if (
            record.system_name == system_name
            and record.surface_name == surface_name
        ):
            return record
    raise SurfaceCatalogError(
        "unknown CHEMPOTPY surface {}/{}".format(system_name, surface_name)
    )


def resolve_atom_symbols(
    record: CatalogSurfaceRecord,
    atom_symbols: Optional[Iterable[str]] = None,
) -> Tuple[str, ...]:
    """Return caller ordering for a fixed surface or explicit variable system."""

    if atom_symbols is not None:
        resolved = tuple(str(symbol) for symbol in atom_symbols)
        if not resolved or any(not symbol for symbol in resolved):
            raise SurfaceCatalogError("atom_symbols must not contain empty values")
        if not record.variable_atom_count and len(resolved) != record.atom_count:
            raise SurfaceCatalogError(
                "{} requires {} atoms, received {}".format(
                    record.key, record.atom_count, len(resolved)
                )
            )
        return resolved
    if record.variable_atom_count:
        raise SurfaceCatalogError(
            "{} has a variable composition; atom_symbols are required".format(
                record.key
            )
        )
    system_module = importlib.import_module(
        "chempotpy.{}".format(record.system_name)
    )
    module_symbols = getattr(system_module, "ATOM_SYMBOLS", None)
    if module_symbols is not None:
        resolved = tuple(str(symbol) for symbol in module_symbols)
    elif record.system_name in _SPECIAL_ATOM_SYMBOLS:
        resolved = _SPECIAL_ATOM_SYMBOLS[record.system_name]
    else:
        resolved = _symbols_from_formula(record.system_name)
    if len(resolved) != record.atom_count:
        raise SurfaceCatalogError(
            "catalog composition for {} has {} atoms; expected {}".format(
                record.key, len(resolved), record.atom_count
            )
        )
    return resolved


def _symbols_from_formula( system_name: str ) -> Tuple[str, ...]:
    formula = system_name[:-1] if system_name[-1:] in {"p", "m"} else system_name
    matches = tuple(_FORMULA_TOKEN.finditer(formula))
    if "".join(match.group(0) for match in matches) != formula:
        raise SurfaceCatalogError(
            "cannot infer atom symbols for system {!r}".format(system_name)
        )
    return tuple(
        symbol
        for match in matches
        for symbol in (match.group(1),) * int(match.group(2) or 1)
    )
