#!/usr/bin/env python3
"""Enforce CHEMPOTPY style on Fortran sources marked as modernized."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPOSITORY_ROOT / "chempotpy"
MODERN_MARKER = "!> @modern"
PROCEDURE_PATTERN = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:[a-z][a-z0-9_]*(?:\([^)]*\))?\s+)?"
    r"(subroutine|function)\s+([a-z][a-z0-9_]*)\s*\(",
    flags=re.IGNORECASE,
)
MODULE_PATTERN = re.compile(
    r"^\s*module\s+(?!procedure\b)([a-z][a-z0-9_]*)",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True)
class StyleViolation:
    """Describe one source-style violation."""

    path: str
    line: int
    rule: str
    message: str


def modern_sources() -> Sequence[Path]:
    """Return free-form sources explicitly opted into the modern style gate."""

    candidates = sorted(PACKAGE_ROOT.rglob("*.f90"))
    return tuple(
        path for path in candidates
        if MODERN_MARKER in path.read_text(errors="replace")
    )


def style_inventory() -> dict:
    """Return repository-wide procedure documentation and spacing counts."""

    source_suffixes = {".f", ".f03", ".f08", ".f90", ".f95", ".pyf"}
    procedure_count = 0
    spaced_count = 0
    brief_count = 0
    files_with_procedures = set()
    for path in sorted(PACKAGE_ROOT.rglob("*")):
        if path.suffix.lower() not in source_suffixes:
            continue
        lines = path.read_text(errors="replace").splitlines()
        for index, line in enumerate(lines):
            match = PROCEDURE_PATTERN.match(line)
            if match is None:
                continue
            procedure_count += 1
            files_with_procedures.add(path)
            opening = line.find("(", match.start(2))
            closing = line.rfind(")")
            opening_spaced = opening >= 0 and line[opening: opening+2] == "( "
            closing_spaced = (
                closing <= opening or line[closing-1: closing+1] == " )"
            )
            if opening_spaced and closing_spaced:
                spaced_count += 1
            previous = index-1
            while previous >= 0 and not lines[previous].strip():
                previous -= 1
            if previous >= 0 and lines[previous].lstrip().startswith(
                "!> @brief"
            ):
                brief_count += 1
    return {
        "fortran_files_with_procedures": len(files_with_procedures),
        "procedure_declarations": procedure_count,
        "spaced_declarations": spaced_count,
        "briefed_procedures": brief_count,
    }


def check_source(path: Path) -> Sequence[StyleViolation]:
    """Return all style violations in one marked modern Fortran source."""

    lines = path.read_text(errors="replace").splitlines()
    relative_path = str(path.relative_to(REPOSITORY_ROOT))
    violations: List[StyleViolation] = []
    modules = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.lstrip()
        if "\t" in line:
            violations.append(StyleViolation(
                relative_path, line_number, "indent-tabs",
                "tabs are not permitted",
            ))
        if stripped and not stripped.startswith("!"):
            indentation = len(line)-len(stripped)
            if indentation % 4 != 0:
                violations.append(StyleViolation(
                    relative_path, line_number, "indent-four",
                    "indentation must be a multiple of four spaces",
                ))
        procedure_match = PROCEDURE_PATTERN.match(line)
        if procedure_match is not None:
            opening = line.find("(", procedure_match.start(2))
            if opening >= 0 and line[opening: opening+2] != "( ":
                violations.append(StyleViolation(
                    relative_path, line_number, "argument-spacing",
                    "procedure arguments must start with '( '",
                ))
            closing = line.rfind(")")
            if closing > opening and line[closing-1: closing+1] != " )":
                violations.append(StyleViolation(
                    relative_path, line_number, "argument-spacing",
                    "procedure arguments must end with ' )'",
                ))
            previous = line_number-2
            while previous >= 0 and not lines[previous].strip():
                previous -= 1
            if previous < 0 or not lines[previous].lstrip().startswith(
                "!> @brief"
            ):
                violations.append(StyleViolation(
                    relative_path, line_number, "procedure-brief",
                    "procedure requires an immediately preceding !> @brief",
                ))
        module_match = MODULE_PATTERN.match(line)
        if module_match is not None:
            modules.append((line_number, module_match.group(1)))
        if re.match(r"^\s*use\s+", line, flags=re.IGNORECASE):
            if not re.search(r",\s*only\s*:", line, flags=re.IGNORECASE):
                violations.append(StyleViolation(
                    relative_path, line_number, "use-only",
                    "module imports must use only:",
                ))

    lower_lines = [line.strip().lower() for line in lines]
    for line_number, module_name in modules:
        end_index = next(
            (
                index for index in range(line_number, len(lines))
                if lower_lines[index].startswith("contains")
                or lower_lines[index].startswith("end module")
            ),
            len(lines),
        )
        specification = lower_lines[line_number:end_index]
        if "implicit none" not in specification:
            violations.append(StyleViolation(
                relative_path, line_number, "implicit-none",
                "module {} requires implicit none".format(module_name),
            ))
        if "private" not in specification:
            violations.append(StyleViolation(
                relative_path, line_number, "private-default",
                "module {} must be private by default".format(module_name),
            ))
    return tuple(violations)


def main() -> int:
    """Check all opted-in sources and fail when any violation is present."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    inventory = style_inventory()
    sources = modern_sources()
    violations = tuple(
        violation
        for source in sources
        for violation in check_source(source)
    )
    if args.json:
        print(json.dumps({
            "inventory": inventory,
            "sources": [str(path.relative_to(REPOSITORY_ROOT)) for path in sources],
            "violations": [asdict(violation) for violation in violations],
        }, indent=2, sort_keys=True))
    else:
        for key, value in inventory.items():
            print("{}: {}".format(key, value))
        print("modern_fortran_sources_checked: {}".format(len(sources)))
        print("modern_fortran_style_violations: {}".format(len(violations)))
        for violation in violations:
            print("{}:{}: {}: {}".format(
                violation.path,
                violation.line,
                violation.rule,
                violation.message,
            ))
    return 1 if violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
