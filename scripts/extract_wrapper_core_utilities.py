#!/usr/bin/env python3
"""Move repeated public-wrapper geometry and unit literals to the core API."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


WRAPPER_START = re.compile(
    r"^\s*subroutine\s+(?P<name>pes|dpem)\s*\(",
    re.IGNORECASE | re.MULTILINE,
)
DISTANCE = re.compile(
    r"sqrt\s*\(\s*"
    r"\(\s*x\(\s*(?P<first>\d+)\s*,\s*1\s*\)\s*-\s*"
    r"x\(\s*(?P<second>\d+)\s*,\s*1\s*\)\s*\)\s*\*\*\s*2\s*\+\s*"
    r"\(\s*x\(\s*(?P=first)\s*,\s*2\s*\)\s*-\s*"
    r"x\(\s*(?P=second)\s*,\s*2\s*\)\s*\)\s*\*\*\s*2\s*\+\s*"
    r"\(\s*x\(\s*(?P=first)\s*,\s*3\s*\)\s*-\s*"
    r"x\(\s*(?P=second)\s*,\s*3\s*\)\s*\)\s*\*\*\s*2\s*\)",
    re.IGNORECASE,
)
LEGACY_BOHR = re.compile(
    r"(?<![a-z0-9_.])0\.529177211(?![a-z0-9_.])",
    re.IGNORECASE,
)
LEGACY_HARTREE = re.compile(
    r"(?<![a-z0-9_.])27\.211386(?![a-z0-9_.])",
    re.IGNORECASE,
)
TRIATOMIC_MATRIX = re.compile(
    r"(?P<indent>^[ \t]*)(?P<name>[a-z][a-z0-9_]*)\(1,\s*1\)\s*=\s*"
    r"cartesian_distance\(x\(1,\s*:\),\s*x\(2,\s*:\)\)\s*/\s*"
    r"legacy_bohr_to_angstrom\s*\n"
    r"(?P=indent)(?P=name)\(1,\s*2\)\s*=\s*"
    r"cartesian_distance\(x\(2,\s*:\),\s*x\(3,\s*:\)\)\s*/\s*"
    r"legacy_bohr_to_angstrom\s*\n"
    r"(?P=indent)(?P=name)\(1,\s*3\)\s*=\s*"
    r"cartesian_distance\(x\(1,\s*:\),\s*x\(3,\s*:\)\)\s*/\s*"
    r"legacy_bohr_to_angstrom",
    re.IGNORECASE | re.MULTILINE,
)
TRIATOMIC_VECTOR = re.compile(
    r"(?P<indent>^[ \t]*)(?P<name>[a-z][a-z0-9_]*)\(1\)\s*=\s*"
    r"cartesian_distance\(x\(1,\s*:\),\s*x\(2,\s*:\)\)\s*/\s*"
    r"legacy_bohr_to_angstrom\s*\n"
    r"(?P=indent)(?P=name)\(2\)\s*=\s*"
    r"cartesian_distance\(x\(2,\s*:\),\s*x\(3,\s*:\)\)\s*/\s*"
    r"legacy_bohr_to_angstrom\s*\n"
    r"(?P=indent)(?P=name)\(3\)\s*=\s*"
    r"cartesian_distance\(x\(1,\s*:\),\s*x\(3,\s*:\)\)\s*/\s*"
    r"legacy_bohr_to_angstrom",
    re.IGNORECASE | re.MULTILINE,
)


def wrapper_span(source: str):
    """Return the source offsets for the public PES or DPEM wrapper."""

    start = WRAPPER_START.search(source)
    if start is None:
        return None
    name = start.group("name")
    subroutine_end = re.compile(
        r"^\s*end\s*subroutine(?:\s+{})?\s*$".format(re.escape(name)),
        re.IGNORECASE | re.MULTILINE,
    ).search(source, start.end())
    if subroutine_end is not None:
        return start.start(), subroutine_end.end()
    bare_end = re.compile(
        r"^\s*end\s*$", re.IGNORECASE | re.MULTILINE
    ).search(source, start.end())
    if bare_end is None:
        raise ValueError("end of {} wrapper not found".format(name))
    return start.start(), bare_end.end()


def add_import(wrapper: str, module_name: str, names: list[str]) -> str:
    """Insert one dependency import immediately before IMPLICIT NONE."""

    if re.search(
        r"^\s*use\s+{}\b".format(re.escape(module_name)),
        wrapper,
        re.IGNORECASE | re.MULTILINE,
    ):
        raise ValueError(
            "{} already imported; merge its only-list manually".format(
                module_name
            )
        )
    implicit = re.search(
        r"^(?P<indent>\s*)implicit\s+none\b",
        wrapper,
        re.IGNORECASE | re.MULTILINE,
    )
    if implicit is None:
        raise ValueError("wrapper has no implicit none")
    statement = "{}use {}, only: {}\n".format(
        implicit.group("indent"), module_name, ", ".join(names)
    )
    return wrapper[:implicit.start()] + statement + wrapper[implicit.start():]


def update_import(
    wrapper: str,
    module_name: str,
    add_names: tuple[str, ...] = (),
    remove_names: tuple[str, ...] = (),
) -> str:
    """Update one single-line ONLY import, removing it when empty."""

    pattern = re.compile(
        r"^(?P<indent>\s*)use\s+{}\s*,\s*only\s*:\s*"
        r"(?P<names>[^\n]*)\n?".format(re.escape(module_name)),
        re.IGNORECASE | re.MULTILINE,
    )
    match = pattern.search(wrapper)
    if match is None:
        if not add_names:
            return wrapper
        return add_import(wrapper, module_name, list(add_names))
    names = [name.strip() for name in match.group("names").split(",")]
    removed = {name.lower() for name in remove_names}
    names = [name for name in names if name.lower() not in removed]
    existing = {name.lower() for name in names}
    names.extend(name for name in add_names if name.lower() not in existing)
    if names:
        replacement = "{}use {}, only: {}\n".format(
            match.group("indent"), module_name, ", ".join(names)
        )
    else:
        replacement = ""
    return wrapper[:match.start()] + replacement + wrapper[match.end():]


def used_outside_import(wrapper: str, name: str) -> bool:
    """Return whether NAME occurs outside a USE statement."""

    without_imports = re.sub(
        r"^\s*use\b[^\n]*$", "", wrapper, flags=re.IGNORECASE | re.MULTILINE
    )
    return re.search(r"\b{}\b".format(re.escape(name)), without_imports, re.I) is not None


def migrate_source(source: str) -> tuple[str, int, int, int, int]:
    """Return one source with wrapper-only core substitutions."""

    span = wrapper_span(source)
    if span is None:
        return source, 0, 0, 0, 0
    start, end = span
    wrapper = source[start:end]

    def distance_call(match) -> str:
        return "cartesian_distance(x({}, :), x({}, :))".format(
            match.group("first"), match.group("second")
        )

    wrapper, distance_count = DISTANCE.subn(distance_call, wrapper)
    wrapper, coordinate_count = LEGACY_BOHR.subn(
        "legacy_bohr_to_angstrom", wrapper
    )
    wrapper, energy_count = LEGACY_HARTREE.subn(
        "legacy_hartree_to_electronvolt", wrapper
    )
    if distance_count:
        wrapper = add_import(
            wrapper, "chempotpy_geometry", ["cartesian_distance"]
        )
    unit_names = []
    if coordinate_count:
        unit_names.append("legacy_bohr_to_angstrom")
    if energy_count:
        unit_names.append("legacy_hartree_to_electronvolt")
    if unit_names:
        wrapper = add_import(
            wrapper, "chempotpy_legacy_units", unit_names
        )

    def matrix_call(match) -> str:
        return "{}call legacy_triatomic_distances_to_bohr( x, {}(1, :) )".format(
            match.group("indent"), match.group("name")
        )

    def vector_call(match) -> str:
        return "{}call legacy_triatomic_distances_to_bohr( x, {} )".format(
            match.group("indent"), match.group("name")
        )

    wrapper, matrix_count = TRIATOMIC_MATRIX.subn(matrix_call, wrapper)
    wrapper, vector_count = TRIATOMIC_VECTOR.subn(vector_call, wrapper)
    triatomic_count = matrix_count+vector_count
    if triatomic_count:
        wrapper = update_import(
            wrapper,
            "chempotpy_legacy_coordinates",
            add_names=("legacy_triatomic_distances_to_bohr",),
        )
        if not used_outside_import(wrapper, "cartesian_distance"):
            wrapper = update_import(
                wrapper,
                "chempotpy_geometry",
                remove_names=("cartesian_distance",),
            )
        if not used_outside_import(wrapper, "legacy_bohr_to_angstrom"):
            wrapper = update_import(
                wrapper,
                "chempotpy_legacy_units",
                remove_names=("legacy_bohr_to_angstrom",),
            )
    return (
        source[:start] + wrapper + source[end:],
        distance_count,
        coordinate_count,
        energy_count,
        triatomic_count,
    )


def main() -> int:
    """Migrate requested sources and print deterministic totals."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    distances = 0
    coordinates = 0
    energies = 0
    triatomic_blocks = 0
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        try:
            migrated, distance_count, coordinate_count, energy_count, triatomic_count = (
                migrate_source(original)
            )
        except ValueError as error:
            raise ValueError("{}: {}".format(path, error)) from error
        changed += migrated != original
        distances += distance_count
        coordinates += coordinate_count
        energies += energy_count
        triatomic_blocks += triatomic_count
        if arguments.write and migrated != original:
            path.write_text(migrated)
    print("fortran_sources_changed: {}".format(changed))
    print("distance_expressions_replaced: {}".format(distances))
    print("coordinate_literals_replaced: {}".format(coordinates))
    print("energy_literals_replaced: {}".format(energies))
    print("triatomic_blocks_replaced: {}".format(triatomic_blocks))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
