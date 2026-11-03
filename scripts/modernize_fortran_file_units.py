#!/usr/bin/env python3
"""Replace hard-coded OPEN units with modern NEWUNIT variables."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PROCEDURE_START = re.compile(
    r"^\s*(?!end\b)(?:(?:recursive|pure|elemental|impure|module)\s+)*"
    r"(?:(?:double\s+precision|double\s+complex|real(?:\s*\([^)]*\))?|"
    r"integer(?:\s*\([^)]*\))?|logical(?:\s*\([^)]*\))?|"
    r"complex(?:\s*\([^)]*\))?|character(?:\s*\([^)]*\))?)\s+)?"
    r"(?:subroutine|function)\b",
    re.IGNORECASE,
)
PROCEDURE_END = re.compile(
    r"^\s*(?:end\s*(?:subroutine|function)\b.*|end)\s*$", re.IGNORECASE
)
OPEN_UNIT = re.compile(
    r"\bopen\s*\(\s*(?:unit\s*=\s*)?(\d+)\b", re.IGNORECASE
)


def code_without_comment(line: str) -> str:
    """Return code before a comment marker outside a character literal."""

    quote = ""
    index = 0
    while index < len(line):
        character = line[index]
        if quote:
            if character == quote:
                if index + 1 < len(line) and line[index + 1] == quote:
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            quote = character
        elif character == "!":
            return line[:index]
        index += 1
    return line


def replace_unit(code: str, unit: str, variable: str) -> str:
    """Replace active I/O references to one unit in a physical code line."""

    comment = code[len(code_without_comment(code)):]
    active = code_without_comment(code)
    active = re.sub(
        rf"\bopen\s*\(\s*unit\s*=\s*{unit}\b",
        "open(newunit={}".format(variable),
        active,
        flags=re.IGNORECASE,
    )
    active = re.sub(
        rf"\bopen\s*\(\s*{unit}\s*,",
        "open(newunit={},".format(variable),
        active,
        flags=re.IGNORECASE,
    )
    for keyword in ("read", "write", "close", "rewind", "backspace", "endfile"):
        active = re.sub(
            rf"\b{keyword}\s*\(\s*unit\s*=\s*{unit}\b",
            "{}(unit={}".format(keyword, variable),
            active,
            flags=re.IGNORECASE,
        )
        active = re.sub(
            rf"\b{keyword}\s*\(\s*{unit}\b",
            "{}({}".format(keyword, variable),
            active,
            flags=re.IGNORECASE,
        )
    return active + comment


def migrate_source(source: str) -> tuple[str, int]:
    """Migrate every procedure containing an active hard-coded OPEN."""

    lines = source.splitlines()
    changed_procedures = 0
    index = 0
    while index < len(lines):
        if PROCEDURE_START.match(lines[index]) is None:
            index += 1
            continue
        end = next(
            candidate for candidate in range(index + 1, len(lines))
            if PROCEDURE_END.match(lines[candidate])
        )
        units = sorted({
            match.group(1)
            for line in lines[index:end + 1]
            for match in OPEN_UNIT.finditer(code_without_comment(line))
        }, key=int)
        if not units:
            index = end + 1
            continue
        for unit in units:
            variable = "file_unit_{}".format(unit)
            for line_index in range(index, end + 1):
                lines[line_index] = replace_unit(
                    lines[line_index], unit, variable
                )
        implicit_none = next((
            candidate for candidate in range(index + 1, end)
            if re.match(
                r"^\s*implicit\s+none\s*$", lines[candidate], re.IGNORECASE
            )
        ), None)
        if implicit_none is not None:
            insertion = implicit_none + 1
        else:
            insertion = index + 1
            while insertion < end and code_without_comment(
                lines[insertion-1]
            ).rstrip().endswith("&"):
                insertion += 1
            for candidate in range(insertion, end):
                if re.match(r"^\s*use(?:\s|,)", lines[candidate], re.IGNORECASE):
                    insertion = candidate + 1
        declarations = [
            "    integer :: {}".format("file_unit_{}".format(unit))
            for unit in units
        ]
        lines[insertion:insertion] = declarations
        end += len(declarations)
        changed_procedures += 1
        index = end + 1
    return "\n".join(lines) + "\n", changed_procedures


def main() -> int:
    """Migrate requested sources and report updated procedure counts."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    procedures = 0
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        modernized, count = migrate_source(original)
        changed += modernized != original
        procedures += count
        if arguments.write and modernized != original:
            path.write_text(modernized)
    print("fortran_sources_changed: {}".format(changed))
    print("procedures_using_newunit: {}".format(procedures))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
