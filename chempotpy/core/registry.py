"""Discover PES sources without relying on a hard-coded system list."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Tuple

SOURCE_SUFFIXES = {".f", ".f90", ".F", ".F90"}
ENTRY_POINT_PATTERN = re.compile(
    r"^\s*(?:recursive\s+)?subroutine\s+(pes|dpem)\s*\(",
    re.IGNORECASE | re.MULTILINE,
)
SOURCE_FORM_OVERRIDES = {
    ("CH5p", "CH5p_PIP"): "fixed",
}


@dataclass(frozen=True)
class SurfaceSource:
    """Identify one Fortran source belonging to a CHEMPOTPY system."""

    system_name: str
    surface_name: str
    path: Path
    source_form: str
    entry_point: Optional[str]

    @property
    def is_surface( self ) -> bool:
        """Return whether the source exposes a CHEMPOTPY entry point."""

        return self.entry_point is not None


def discover_surface_sources(
    package_root: Optional[Path] = None,
) -> Tuple[SurfaceSource, ...]:
    """Return every visible Fortran PES source in deterministic order."""

    root = package_root or Path(__file__).resolve().parents[1]
    sources = []
    for system_directory in sorted(root.iterdir()):
        if not system_directory.is_dir():
            continue
        if system_directory.name.startswith((".", "_")):
            continue
        if not (system_directory / "__init__.py").is_file():
            continue
        for path in sorted(system_directory.iterdir()):
            if path.is_file() and path.suffix in SOURCE_SUFFIXES:
                match = ENTRY_POINT_PATTERN.search(
                    path.read_text(errors="replace")
                )
                sources.append(
                    SurfaceSource(
                        system_name=system_directory.name,
                        surface_name=path.stem,
                        path=path,
                        source_form=SOURCE_FORM_OVERRIDES.get(
                            (system_directory.name, path.stem),
                            "free" if path.suffix.lower() == ".f90" else "fixed",
                        ),
                        entry_point=(match.group(1).lower() if match else None),
                    )
                )
    return tuple(sources)


def duplicate_surface_sources(
    sources: Iterable[SurfaceSource],
) -> Dict[Tuple[str, str], Tuple[SurfaceSource, ...]]:
    """Return surfaces that have multiple competing source files."""

    grouped: Dict[Tuple[str, str], List[SurfaceSource]] = {}
    for source in sources:
        grouped.setdefault(
            (source.system_name, source.surface_name), []
        ).append(source)
    return {
        key: tuple(value)
        for key, value in grouped.items()
        if len(value) > 1
    }
