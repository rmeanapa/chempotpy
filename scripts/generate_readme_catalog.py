#!/usr/bin/env python3
"""Generate the complete PES availability table in the project README."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys
from typing import Dict, Iterable


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import SurfaceSource, discover_surface_sources


README = REPOSITORY_ROOT / "README.md"
FIXTURE = REPOSITORY_ROOT / "tests" / "data" / "all_surface_references.json"
START_MARKER = "<!-- BEGIN GENERATED PES CATALOG -->"
END_MARKER = "<!-- END GENERATED PES CATALOG -->"
AVAILABILITY_TOKENS = {"P", "G", "D", "U", "UG"}


def curated_availability(source: SurfaceSource) -> str:
    """Return the availability code maintained in a system's short catalog."""

    initialization = source.path.parent / "__init__.py"
    if not initialization.is_file():
        return ""
    surface_pattern = re.compile(
        r"\b{}\s*:".format(re.escape(source.surface_name))
    )
    availability_pattern = re.compile(r"([A-Z]+(?:/[A-Z]+)*)\s*$")
    for line in initialization.read_text(errors="replace").splitlines():
        if surface_pattern.search(line) is None:
            continue
        match = availability_pattern.search(line.strip())
        if match is None:
            continue
        tokens = match.group(1).split("/")
        if tokens and set(tokens) <= AVAILABILITY_TOKENS:
            return "/".join(tokens)
    return ""


def supports_first_derivatives(source: SurfaceSource) -> bool:
    """Return whether source metadata declares first or higher derivatives."""

    match = re.search(
        r"^\s*!\s*Number\s+of\s+derivatives\s*:\s*(.*?)\s*$",
        source.path.read_text(errors="replace"),
        flags=re.IGNORECASE | re.MULTILINE,
    )
    if match is None:
        raise ValueError(
            "{} has neither curated availability nor derivative metadata".format(
                source.path.relative_to(REPOSITORY_ROOT)
            )
        )
    value = match.group(1).strip().lower()
    number = re.match(r"(\d+)", value)
    if number is not None:
        return int(number.group(1)) >= 1
    if value.startswith("first"):
        return True
    raise ValueError(
        "unsupported derivative metadata {!r} in {}".format(
            match.group(1), source.path.relative_to(REPOSITORY_ROOT)
        )
    )


def pes_availability(source: SurfaceSource, state_count: int) -> str:
    """Return a PES availability code from curated or source metadata."""

    curated = curated_availability(source)
    if curated:
        return curated
    if not supports_first_derivatives(source):
        return "P"
    return "P/G/D" if state_count > 1 else "P/G"


def catalog_availability(
    sources: Iterable[SurfaceSource], records: Dict[str, dict]
) -> Dict[str, str]:
    """Return availability codes for every callable catalog source."""

    source_map = {
        "{}/{}".format(source.system_name, source.surface_name): source
        for source in sources
    }
    availability: Dict[str, str] = {}
    for key, source in source_map.items():
        if source.entry_point == "pes":
            availability[key] = pes_availability(
                source, int(records[key]["nstates"])
            )
    for key, source in source_map.items():
        if source.entry_point != "dpem":
            continue
        curated = curated_availability(source)
        if curated:
            availability[key] = curated
            continue
        pes_key = key[:-len("_DPEM")]
        if pes_key not in availability:
            raise ValueError("{} has no paired PES availability".format(key))
        availability[key] = "U/UG" if "G" in availability[pes_key] else "U"
    return availability


def generated_catalog() -> str:
    """Render the complete Markdown availability table."""

    sources = tuple(sorted(
        (
            source for source in discover_surface_sources()
            if source.is_surface
        ),
        key=lambda source: (
            source.system_name.casefold(), source.surface_name.casefold()
        ),
    ))
    records = json.loads(FIXTURE.read_text())["surfaces"]
    keys = {
        "{}/{}".format(source.system_name, source.surface_name)
        for source in sources
    }
    if keys != set(records):
        raise ValueError("README catalog and numerical fixture differ")
    availability = catalog_availability(sources, records)
    systems = {source.system_name for source in sources}
    lines = [
        "The catalog currently exposes **{} callable surfaces** across "
        "**{} chemical systems**. Availability codes are: **P** potential "
        "energy, **G** Cartesian gradient, **D** nonadiabatic coupling vector, "
        "**U** diabatic potential matrix, and **UG** Cartesian gradient of "
        "the diabatic matrix.".format(len(sources), len(systems)),
        "",
        "| System | Surface | Interface | States | Available |",
        "| :--- | :--- | :---: | ---: | :--- |",
    ]
    for source in sources:
        key = "{}/{}".format(source.system_name, source.surface_name)
        lines.append(
            "| `{}` | `{}` | {} | {} | **{}** |".format(
                source.system_name,
                source.surface_name,
                source.entry_point.upper(),
                records[key]["nstates"],
                availability[key],
            )
        )
    return "\n".join(lines)


def updated_readme(original: str) -> str:
    """Return README text with the generated catalog replaced."""

    if original.count(START_MARKER) != 1 or original.count(END_MARKER) != 1:
        raise ValueError("README must contain one PES catalog marker pair")
    start = original.index(START_MARKER)+len(START_MARKER)
    end = original.index(END_MARKER)
    if start > end:
        raise ValueError("README PES catalog markers are reversed")
    return (
        original[:start]
        + "\n\n"
        + generated_catalog()
        + "\n\n"
        + original[end:]
    )


def main() -> int:
    """Write or verify the README PES catalog."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    original = README.read_text()
    updated = updated_readme(original)
    if arguments.check:
        if original != updated:
            print("outdated: README.md PES catalog")
            return 1
    else:
        README.write_text(updated)
        print("README catalog surfaces: {}".format(
            len(json.loads(FIXTURE.read_text())["surfaces"])
        ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
