#!/usr/bin/env python3
"""Replace duplicated POTLIB CARTOU bodies with a bundled core call."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


CARTOU_START = re.compile(r"^\s*subroutine\s+cartou\s*$", re.IGNORECASE)
CARTOU_END = re.compile(
    r"^\s*end\s+subroutine\s+cartou\s*$", re.IGNORECASE
)
CARTTOR_START = re.compile(r"^\s*subroutine\s+carttor\s*$", re.IGNORECASE)
CARTTOR_END = re.compile(
    r"^\s*end\s+subroutine\s+carttor\s*$", re.IGNORECASE
)
BODY_START = re.compile(
    r"^\s*if\s*\(\s*ireorder\s*\.eq\.\s*1\s*\)\s*then\s*$",
    re.IGNORECASE,
)
BODY_END = re.compile(r"^\s*end\s*if\s*$", re.IGNORECASE)
CORE_USE = (
    "    use chempotpy_legacy_coordinates, only: convert_cartesian_input"
)
CORE_CALL = (
    "    call convert_cartesian_input( cart, nulbl, natoms, cnvrtd, "
    "ireorder == 1, cartnu )"
)
CARTTOR_USE = (
    "    use chempotpy_legacy_coordinates, only: "
    "cartesian_to_potlib_coordinates"
)
CARTTOR_STATUS = "    integer :: coordinate_status"
CARTTOR_CALL = (
    "    call cartesian_to_potlib_coordinates( cartnu, natoms, icartr, r, "
    "coordinate_status )"
)
CARTTOR_ERROR = [
    "    if (coordinate_status /= 0) then",
    "        write(nflag(18),'(2x,''WRONG ICARTR FOR CARTNU; ICARTR ='',i5//)') icartr",
    "        stop",
    "    end if",
]


def migrate_source(source: str) -> tuple[str, int]:
    """Return source with each canonical CARTOU implementation extracted."""

    lines = source.splitlines()
    migrated = 0
    index = 0
    while index < len(lines):
        if CARTOU_START.match(lines[index]) is None:
            index += 1
            continue
        end = next(
            candidate
            for candidate in range(index + 1, len(lines))
            if CARTOU_END.match(lines[candidate])
        )
        if any("chempotpy_legacy_coordinates" in line.lower()
               for line in lines[index:end]):
            index = end + 1
            continue
        body_start = next(
            candidate
            for candidate in range(index + 1, end)
            if BODY_START.match(lines[candidate])
        )
        body_end = next(
            candidate
            for candidate in range(end - 1, body_start, -1)
            if BODY_END.match(lines[candidate])
        )
        lines.insert(index + 1, CORE_USE)
        body_start += 1
        body_end += 1
        lines[body_start:body_end + 1] = [CORE_CALL]
        migrated += 1
        index = body_start + 2
    return "\n".join(lines) + "\n", migrated


def migrate_carttor(source: str) -> tuple[str, int]:
    """Return source with each canonical CARTTOR implementation extracted."""

    lines = source.splitlines()
    migrated = 0
    index = 0
    while index < len(lines):
        if CARTTOR_START.match(lines[index]) is None:
            index += 1
            continue
        end = next(
            candidate
            for candidate in range(index + 1, len(lines))
            if CARTTOR_END.match(lines[candidate])
        )
        if any("cartesian_to_potlib_coordinates" in line.lower()
               for line in lines[index:end]):
            index = end + 1
            continue
        body_start = next(
            candidate for candidate in range(index + 1, end)
            if re.match(r"^\s*if\s*\(\s*icartr", lines[candidate], re.IGNORECASE)
        )
        body_end = next(
            candidate for candidate in range(end - 1, body_start, -1)
            if BODY_END.match(lines[candidate])
        )
        lines.insert(index + 1, CARTTOR_USE)
        body_start += 1
        body_end += 1
        implicit_none = next(
            candidate for candidate in range(index + 1, body_start)
            if re.match(r"^\s*implicit\s+none\s*$", lines[candidate], re.IGNORECASE)
        )
        lines.insert(implicit_none + 1, CARTTOR_STATUS)
        body_start += 1
        body_end += 1
        lines[body_start:body_end + 1] = [CARTTOR_CALL, *CARTTOR_ERROR]
        migrated += 1
        index = body_start + len(CARTTOR_ERROR) + 2
    return "\n".join(lines) + "\n", migrated


def main() -> int:
    """Migrate requested sources and report the number of extracted copies."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    source_count = 0
    routine_count = 0
    carttor_count = 0
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        migrated, count = migrate_source(original)
        migrated, carttor_migrated = migrate_carttor(migrated)
        source_count += migrated != original
        routine_count += count
        carttor_count += carttor_migrated
        if arguments.write and migrated != original:
            path.write_text(migrated)
    print("fortran_sources_changed: {}".format(source_count))
    print("cartou_routines_extracted: {}".format(routine_count))
    print("carttor_routines_extracted: {}".format(carttor_count))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
