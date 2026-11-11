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
    Rule("split-elseif", re.compile(r"^\s*else\s+if\s*\(", re.I), "use the compact ELSEIF keyword"),
    Rule("numeric-label", re.compile(r"^\s*\d+\s+\S", re.I), "remove obsolete statement labels"),
    Rule("continue", re.compile(r"^\s*continue\s*$", re.I), "use ENDDO or remove the no-op"),
    Rule("standalone-format", re.compile(r"^\s*format\s*\(", re.I), "place the format in the WRITE statement"),
    Rule("standalone-dimension", re.compile(r"^\s*dimension\b", re.I), "put array shapes in typed declarations"),
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
        lines = path.read_text(errors="replace").splitlines()
        first_system_line = next(
            (
                line_number
                for line_number, line in enumerate(lines, 1)
                if re.match(r"^\s*!\s*System\s*:", line, re.I)
            ),
            None,
        )
        if first_system_line is not None:
            first_content_line = next(
                (
                    line_number
                    for line_number, line in enumerate(lines, 1)
                    if line.strip()
                ),
                None,
            )
            if first_system_line != first_content_line:
                found.append(
                    f"{path.relative_to(REPOSITORY_ROOT)}:{first_system_line}: "
                    "metadata-position: put the PES metadata block at the "
                    "beginning of the file"
                )
        for line_number, line in enumerate(lines, 1):
            code = active_code(line)
            for rule in RULES:
                if rule.pattern.search(code):
                    found.append(
                        f"{path.relative_to(REPOSITORY_ROOT)}:{line_number}: "
                        f"{rule.name}: {rule.message}: {code.strip()}"
                    )
            if re.fullmatch(r"\s*contains\s*", code, re.I) is not None:
                before = line_number >= 2 and not lines[line_number-2].strip()
                after = line_number < len(lines) and not lines[line_number].strip()
                extra_before = line_number >= 3 and not lines[line_number-3].strip()
                extra_after = line_number+1 < len(lines) and not lines[line_number+1].strip()
                if not before or not after or extra_before or extra_after:
                    found.append(
                        f"{path.relative_to(REPOSITORY_ROOT)}:{line_number}: "
                        "contains-spacing: put exactly one blank line before "
                        "and after CONTAINS"
                    )
            if re.fullmatch(
                r"\s*end\s+(?:subroutine|function)(?:\s+[a-z][a-z0-9_]*)?\s*",
                code,
                re.I,
            ) is not None:
                next_code = line_number
                while next_code < len(lines) and (
                    not lines[next_code].strip()
                    or lines[next_code].lstrip().startswith("!")
                ):
                    next_code += 1
                if (
                    next_code < len(lines)
                    and re.match(
                        r"^\s*(?:(?:recursive|pure|elemental|impure|module)\s+)*"
                        r"(?:(?:double\s+precision|real|integer|logical|complex|character)"
                        r"(?:\s*\([^)]*\)|\*\d+)?\s+)?(?:subroutine|function)\b",
                        active_code(lines[next_code]),
                        re.I,
                    )
                    is not None
                ):
                    one_blank = (
                        line_number < len(lines) and not lines[line_number].strip()
                    )
                    extra_blank = (
                        line_number+1 < len(lines)
                        and not lines[line_number+1].strip()
                    )
                    if not one_blank or extra_blank:
                        found.append(
                            f"{path.relative_to(REPOSITORY_ROOT)}:{line_number}: "
                            "procedure-spacing: put exactly one blank line "
                            "between adjacent procedures"
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
