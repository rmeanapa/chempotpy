#!/usr/bin/env python3
"""Fail when a catalog PES reintroduces a forbidden legacy construct."""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = REPOSITORY_ROOT / "chempotpy"


@dataclass(frozen=True)
class Rule:
    name: str
    pattern: re.Pattern[str]
    message: str


RULES = (
    Rule("data", re.compile(r"^\s*data\b", re.I), "replace DATA with parameters or explicit assignments"),
    Rule("save", re.compile(r"\bsave\b", re.I), "remove hidden saved procedure state"),
    Rule("common", re.compile(r"^\s*common\b", re.I), "encapsulate shared state in a private module"),
    Rule("entry", re.compile(r"^\s*entry\b", re.I), "use one explicit procedure per interface"),
    Rule("block-data", re.compile(r"^\s*block\s+data\b", re.I), "use module parameters or module-owned data"),
    Rule("goto", re.compile(r"\b(?:go\s*to|goto)\b", re.I), "use structured control flow"),
    Rule("arithmetic-if", re.compile(r"^\s*if\s*\(.*\)\s*\d+\s*,\s*\d+\s*,\s*\d+", re.I), "use IF/ELSEIF/ENDIF"),
    Rule("numeric-label", re.compile(r"^\s*\d+\s+\S", re.I), "remove obsolete statement labels"),
    Rule("continue", re.compile(r"^\s*continue\s*$", re.I), "use ENDDO or remove the no-op"),
    Rule("standalone-format", re.compile(r"^\s*format\s*\(", re.I), "place the format in the WRITE statement"),
    Rule("implicit-typing", re.compile(r"^\s*implicit\s+(?!none\b)", re.I), "declare every variable explicitly"),
    Rule("prepot", re.compile(r"^\s*(?:call\s+|subroutine\s+)\w*prepot\w*\b", re.I), "evaluate through the owning PES module without a setup entry point"),
    Rule("hollerith-format", re.compile(r"\b\d+H[A-Za-z]", re.I), "use a quoted character edit descriptor"),
    Rule("adifor", re.compile(r"ADIFOR", re.I), "remove generated ADIFOR boilerplate"),
)


def active_code(line: str) -> str:
    """Return the statement portion of a free-form Fortran source line."""
    quote = None
    index = 0
    while index < len(line):
        character = line[index]
        if quote:
            if character == quote:
                if index + 1 < len(line) and line[index + 1] == quote:
                    index += 2
                    continue
                quote = None
        elif character in {"'", '"'}:
            quote = character
        elif character == "!":
            return line[:index]
        index += 1
    return line


def violations(source_root: Path = SOURCE_ROOT) -> list[str]:
    """Return formatted violations for all modern free-form PES sources."""
    found: list[str] = []
    legacy_files = sorted(
        path for path in source_root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".f", ".for", ".f77"}
    )
    for path in legacy_files:
        found.append(f"{path.relative_to(REPOSITORY_ROOT)}: legacy fixed-form source remains")
    for path in sorted(source_root.rglob("*.f90")):
        for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            code = active_code(line)
            for rule in RULES:
                if rule.pattern.search(code):
                    found.append(
                        f"{path.relative_to(REPOSITORY_ROOT)}:{line_number}: "
                        f"{rule.name}: {rule.message}: {code.strip()}"
                    )
    return found


def main() -> int:
    found = violations()
    if found:
        print("modernization violations:")
        print("\n".join(f"  {item}" for item in found))
        print(f"modernization_violations: {len(found)}")
        return 1
    source_count = sum(1 for _ in SOURCE_ROOT.rglob("*.f90"))
    print(f"modernization_sources_checked: {source_count}")
    print("modernization_violations: 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
