#!/usr/bin/env python3
"""Give Fortran modules a private default and explicit public API."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, List, Sequence

try:
    from scripts.explicit_fortran_types import (
        SYMBOL_PATTERN,
        compiler_dump,
        declaration_lines,
    )
except ModuleNotFoundError:  # Direct execution from the scripts directory.
    from explicit_fortran_types import (
        SYMBOL_PATTERN,
        compiler_dump,
        declaration_lines,
    )


MODULE_START = re.compile(
    r"^(?P<indent>\s*)module\s+(?!procedure\b)(?P<name>[a-z][a-z0-9_]*)",
    flags=re.IGNORECASE,
)
PROCEDURE_START = re.compile(
    r"^\s*(?!end\b)(?:(?:recursive|pure|elemental|impure|module)\s+)*"
    r"(?:(?:double\s+precision|real(?:\s*\([^)]*\)|\*\d+)?|integer|logical|"
    r"complex|character(?:\s*\([^)]*\)|\*\d+)?)\s+)?"
    r"(?:subroutine|function)\s+([a-z][a-z0-9_]*)",
    flags=re.IGNORECASE,
)


def use_associations(dump: str) -> DefaultDict[str, set[str]]:
    """Return compiler-resolved public names imported from each module."""

    associations: DefaultDict[str, set[str]] = defaultdict(set)
    for match in SYMBOL_PATTERN.finditer(dump):
        symbol = match.group(1).lower()
        module_match = re.search(
            r"USE-ASSOC\(([^)]+)\)", match.group("body"), re.IGNORECASE
        )
        if module_match is None:
            continue
        module_name = module_match.group(1).lower()
        if symbol != module_name and re.fullmatch(r"[a-z][a-z0-9_]*", symbol):
            associations[module_name].add(symbol)
    return associations


def module_ranges(lines: Sequence[str]) -> List[tuple[str, int, int, str]]:
    """Locate module specification and body ranges."""

    modules = []
    for start, line in enumerate(lines):
        match = MODULE_START.match(line)
        if match is None:
            continue
        name = match.group("name").lower()
        end = next(
            (
                index for index in range(start + 1, len(lines))
                if re.match(
                    r"^\s*end\s+module(?:\s+{})?\b".format(re.escape(name)),
                    lines[index],
                    re.IGNORECASE,
                )
            ),
            len(lines),
        )
        modules.append((name, start, end, match.group("indent")))
    return modules


def explicitly_public_names(specification: Sequence[str]) -> set[str]:
    """Return names already carrying an explicit public attribute."""

    logical_lines: List[str] = []
    pending = ""
    for source_line in specification:
        code = source_line.split("!", 1)[0].strip()
        if not code:
            continue
        if pending:
            code = pending + " " + code.lstrip("& ")
        if code.endswith("&"):
            pending = code[:-1].rstrip()
            continue
        logical_lines.append(code)
        pending = ""
    if pending:
        logical_lines.append(pending)

    names: set[str] = set()
    for line in logical_lines:
        if not re.search(r"\bpublic\b", line, re.IGNORECASE):
            continue
        items_match = re.search(r"::\s*(.*)$", line)
        if items_match is None:
            items_match = re.match(r"\s*public\s+(.*)$", line, re.IGNORECASE)
        if items_match is None:
            continue
        for item in items_match.group(1).split(","):
            name_match = re.match(r"\s*([a-z][a-z0-9_]*)", item, re.IGNORECASE)
            if name_match is not None:
                names.add(name_match.group(1).lower())
    return names


def rewrite_source(source_text: str, associations: dict[str, set[str]]) -> str:
    """Return source with private module defaults and explicit public names."""

    lines = source_text.splitlines()
    edits: List[tuple[int, List[str]]] = []
    for name, start, end, indent in module_ranges(lines):
        contains = next(
            (
                index for index in range(start + 1, end)
                if re.match(r"^\s*contains\b", lines[index], re.IGNORECASE)
            ),
            end,
        )
        specification = lines[start + 1:contains]
        if any(re.match(r"^\s*private\s*$", line, re.IGNORECASE) for line in specification):
            continue
        public_names = set(associations.get(name, set()))
        for line in lines[contains + 1:end]:
            procedure = PROCEDURE_START.match(line)
            if procedure is not None and procedure.group(1).lower() in {"pes", "dpem"}:
                public_names.add(procedure.group(1).lower())
        public_names.difference_update(explicitly_public_names(specification))
        implicit_index = next(
            (
                index for index in range(start + 1, contains)
                if re.match(r"^\s*implicit\s+none\b", lines[index], re.IGNORECASE)
            ),
            None,
        )
        if implicit_index is None:
            raise RuntimeError("module {} lacks implicit none".format(name))
        body_indent = indent + "    "
        inserted = [body_indent + "private"]
        if public_names:
            inserted.extend(declaration_lines(
                "public", public_names, body_indent
            ))
        edits.append((implicit_index + 1, inserted + [""]))
    for index, inserted in sorted(edits, reverse=True):
        lines[index:index] = inserted
    return "\n".join(lines) + "\n"


def main() -> int:
    """Privatize modules in requested source files."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    for source in arguments.sources:
        original = source.read_text(errors="strict")
        if not any(MODULE_START.match(line) for line in original.splitlines()):
            continue
        associations = use_associations(compiler_dump(source, arguments.compiler))
        rewritten = rewrite_source(original, associations)
        if rewritten != original:
            changed += 1
            if arguments.write:
                source.write_text(rewritten)
    print("fortran_module_sources_changed: {}".format(changed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
