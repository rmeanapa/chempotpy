#!/usr/bin/env python3
"""Generate the installed catalog used by external PES providers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import discover_surface_sources


REFERENCE_CATALOG = (
    REPOSITORY_ROOT / "tests" / "data" / "all_surface_references.json"
)
PACKAGE_CATALOG = (
    REPOSITORY_ROOT / "chempotpy" / "core" / "provider_catalog.json"
)
LEGACY_API = REPOSITORY_ROOT / "chempotpy" / "__init__.py"
AVAILABILITY_PATTERN = re.compile(
    r"(?m)^\s*\d+\.\s+(?:[SMZ]\d+\.\s+)?"
    r"([A-Za-z0-9_]+):?\s+.*?\b(P(?:/G)?(?:/D)?|U(?:/UG)?)\s*$"
)
DERIVATIVE_PATTERN = re.compile(
    r"(?im)^!\s*Number of derivatives\s*:\s*(\d+)"
)
VARIABLE_ATOM_SYSTEMS = {"AlmHn", "Aln", "SiOH"}


def availability_map() -> dict[str, str]:
    """Return historical capability annotations keyed by surface name."""

    return dict(AVAILABILITY_PATTERN.findall(LEGACY_API.read_text()))


def generated_catalog() -> dict:
    """Return the compact installed provider catalog."""

    references = json.loads(REFERENCE_CATALOG.read_text())["surfaces"]
    availability = availability_map()
    sources = {
        "{}/{}".format(source.system_name, source.surface_name): source
        for source in discover_surface_sources()
        if source.is_surface
    }
    if set(sources) != set(references):
        raise ValueError("provider catalog and numerical fixture differ")
    records = []
    for key in sorted(sources, key=str.casefold):
        source = sources[key]
        reference = references[key]
        source_text = source.path.read_text(errors="replace")[:5000]
        derivative_match = DERIVATIVE_PATTERN.search(source_text)
        derivative_count = (
            int(derivative_match.group(1)) if derivative_match else 0
        )
        annotation = availability.get(source.surface_name, "")
        analytic_gradient = (
            "/G" in annotation
            or "/UG" in annotation
            or derivative_count >= 1
        )
        records.append({
            "system_name": source.system_name,
            "surface_name": source.surface_name,
            "entry_point": reference["entry_point"],
            "atom_count": reference["natoms"],
            "state_count": reference["nstates"],
            "analytic_gradient": analytic_gradient,
            "nonadiabatic_coupling": "/D" in annotation,
            "variable_atom_count": (
                source.system_name in VARIABLE_ATOM_SYSTEMS
            ),
        })
    return {"schema_version": 1, "surfaces": records}


def main() -> int:
    """Write the catalog or verify that the committed file is current."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    rendered = json.dumps(generated_catalog(), indent=2) + "\n"
    if arguments.check:
        if not PACKAGE_CATALOG.is_file() or PACKAGE_CATALOG.read_text() != rendered:
            print("outdated: chempotpy/core/provider_catalog.json")
            return 1
    else:
        PACKAGE_CATALOG.write_text(rendered)
        print("provider catalog written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
