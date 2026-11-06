#!/usr/bin/env python3
"""Enforce CHEMPOTPY style on every free-form Fortran source."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Sequence

try:
    from scripts.format_fortran import format_source_text
except ModuleNotFoundError:  # Direct execution from the scripts directory.
    from format_fortran import format_source_text


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPOSITORY_ROOT / "chempotpy"
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
LEGACY_DECLARATION_PATTERN = re.compile(
    r"^\s*(?:double\s+precision|double\s+complex|"
    r"real(?:\s*\([^)]*\))?|integer(?:\s*\([^)]*\))?|"
    r"logical(?:\s*\([^)]*\))?|complex(?:\s*\([^)]*\))?|"
    r"character(?:\s*\([^)]*\))?)(?:\s*,[^:]*)?\s+"
    r"(?!function\b|subroutine\b)(?!.*::)",
    flags=re.IGNORECASE,
)
FLATTEN_LOOP = re.compile(
    r"(?P<indent>^[ \t]*)do\s+(?P<atom>\w+)\s*=\s*1\s*,\s*(?P<count>\w+)\s*\n"
    r"[ \t]*do\s+(?P<direction>\w+)\s*=\s*1\s*,\s*3\s*\n"
    r"[ \t]*(?P<index>\w+)\s*=\s*\(\s*(?P=atom)\s*-\s*1\s*\)\s*\*\s*3\s*\+\s*(?P=direction)\s*\n"
    r"[ \t]*(?P<flat>\w+)\s*\(\s*(?P=index)\s*\)\s*=\s*"
    r"(?P<matrix>\w+)\s*\(\s*(?P=atom)\s*,\s*(?P=direction)\s*\)\s*\n"
    r"[ \t]*enddo\s*\n[ \t]*enddo",
    flags=re.IGNORECASE | re.MULTILINE,
)
CHAIN_RULE_LOOP_HEADER = (
    r"[ \t]*do\s+(?P<cart_index>\w+)\s*=\s*1\s*,[^\n]+\n"
    r"[ \t]*do\s+(?P<internal_index>\w+)\s*=\s*1\s*,[^\n]+\n"
)
CHAIN_RULE_LOOP_END = r"[ \t]*enddo\s*\n[ \t]*enddo"
CHAIN_RULE_PATTERNS = (
    re.compile(
        r"(?P<indent>^[ \t]*)(?P<tensor_output>\w+)\s*=\s*0\.d0\s*\n"
        r"[ \t]*(?P<matrix_output>\w+)\s*=\s*0\.d0\s*\n"
        + CHAIN_RULE_LOOP_HEADER
        + r"[ \t]*(?P=tensor_output)\(:,:\s*,\s*(?P=cart_index)\)\s*=\s*"
        r"(?P=tensor_output)\(:,:\s*,\s*(?P=cart_index)\)\s*\+\s*"
        r"(?P<tensor_input>\w+)\(:,:\s*,\s*(?P=internal_index)\)\s*\*\s*"
        r"(?P<jacobian>\w+)\((?P=internal_index)\s*,\s*(?P=cart_index)\)\s*\n"
        r"[ \t]*(?P=matrix_output)\(:,\s*(?P=cart_index)\)\s*=\s*"
        r"(?P=matrix_output)\(:,\s*(?P=cart_index)\)\s*\+\s*"
        r"(?P<matrix_input>\w+)\(:,\s*(?P=internal_index)\)\s*\*\s*"
        r"(?P=jacobian)\((?P=internal_index)\s*,\s*(?P=cart_index)\)\s*\n"
        + CHAIN_RULE_LOOP_END,
        flags=re.IGNORECASE | re.MULTILINE,
    ),
    re.compile(
        r"(?P<indent>^[ \t]*)(?P<output>\w+)\s*=\s*0\.d0\s*\n"
        + CHAIN_RULE_LOOP_HEADER
        + r"[ \t]*(?P=output)\((?P=cart_index)\)\s*=\s*(?P=output)\((?P=cart_index)\)\s*\+\s*"
        r"(?P<input>\w+)\((?P=internal_index)\)\s*\*\s*"
        r"(?P<jacobian>\w+)\((?P=internal_index)\s*,\s*(?P=cart_index)\)\s*\n"
        + CHAIN_RULE_LOOP_END,
        flags=re.IGNORECASE | re.MULTILINE,
    ),
    re.compile(
        r"(?P<indent>^[ \t]*)(?P<output>\w+)\s*=\s*0\.d0\s*\n"
        + CHAIN_RULE_LOOP_HEADER
        + r"[ \t]*(?P=output)\(:,\s*(?P=cart_index)\)\s*=\s*"
        r"(?P=output)\(:,\s*(?P=cart_index)\)\s*\+\s*"
        r"(?P<input>\w+)(?:\(:,\s*|\()(?P=internal_index)\)\s*\*\s*"
        r"(?P<jacobian>\w+)\((?P=internal_index)\s*,\s*(?P=cart_index)\)\s*\n"
        + CHAIN_RULE_LOOP_END,
        flags=re.IGNORECASE | re.MULTILINE,
    ),
    re.compile(
        r"(?P<indent>^[ \t]*)(?P<output>\w+)\s*=\s*0\.d0\s*\n"
        + CHAIN_RULE_LOOP_HEADER
        + r"[ \t]*(?P=output)\(:,:\s*,\s*(?P=cart_index)\)\s*=\s*"
        r"(?P=output)\(:,:\s*,\s*(?P=cart_index)\)\s*\+\s*"
        r"(?P<input>\w+)\(:,:\s*,\s*(?P=internal_index)\)\s*\*\s*"
        r"(?P<jacobian>\w+)\((?P=internal_index)\s*,\s*(?P=cart_index)\)\s*\n"
        + CHAIN_RULE_LOOP_END,
        flags=re.IGNORECASE | re.MULTILINE,
    ),
)


def has_duplicated_gradient_chain_rule(source: str) -> bool:
    """Return whether a wrapper still contains a migrated chain-rule loop."""

    if "call flatten_cartesian_coordinates" not in source.lower():
        return False
    wrapper = re.search(
        r"(?ims)^subroutine\s+(?:pes|dpem)\b.*?^end\s+subroutine(?:\s+(?:pes|dpem))?\s*$",
        source,
    )
    return wrapper is not None and any(pattern.search(wrapper.group(0)) for pattern in CHAIN_RULE_PATTERNS)


@dataclass(frozen=True)
class StyleViolation:
    """Describe one source-style violation."""

    path: str
    line: int
    rule: str
    message: str


def modern_sources() -> Sequence[Path]:
    """Return every free-form source governed by the modern style gate."""

    return tuple(sorted(PACKAGE_ROOT.rglob("*.f90")))


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
            closing = line.find(")", opening)
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
    """Return all style violations in one modern Fortran source."""

    lines = path.read_text(errors="replace").splitlines()
    relative_path = str(path.relative_to(REPOSITORY_ROOT))
    violations: List[StyleViolation] = []
    source_text = path.read_text(errors="replace")
    if format_source_text(source_text) != source_text:
        violations.append(StyleViolation(
            relative_path, 1, "structural-indent",
            "source does not match four-space structural indentation",
        ))
    for match in FLATTEN_LOOP.finditer(source_text):
        violations.append(StyleViolation(
            relative_path,
            source_text.count("\n", 0, match.start())+1,
            "duplicated-cartesian-flattening",
            "use chempotpy_geometry flatten_cartesian_coordinates",
        ))
    if has_duplicated_gradient_chain_rule(source_text):
        violations.append(StyleViolation(
            relative_path,
            1,
            "duplicated-gradient-chain-rule",
            "use chempotpy_geometry transform_internal_gradient_to_cartesian",
        ))
    modules = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.lstrip()
        code = line.split("!", 1)[0]
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
            closing = line.find(")", opening)
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
        if re.match(r"^\s*use(?:\s|,)", line, flags=re.IGNORECASE):
            if not re.search(r",\s*only\s*:", line, flags=re.IGNORECASE):
                violations.append(StyleViolation(
                    relative_path, line_number, "use-only",
                    "module imports must use only:",
                ))
        if re.match(r"^\s*implicit\b", line, flags=re.IGNORECASE) and not re.match(
            r"^\s*implicit\s+none\b", line, flags=re.IGNORECASE
        ):
            violations.append(StyleViolation(
                relative_path, line_number, "implicit-typing",
                "implicit typing maps are forbidden",
            ))
        if re.search(
            r"\b(?:real|integer|complex|logical|character)\s*\*",
            code,
            flags=re.IGNORECASE,
        ):
            violations.append(StyleViolation(
                relative_path, line_number, "star-kind",
                "legacy star-sized intrinsic declarations are forbidden",
            ))
        if re.search(r"\b(?:ADIFOR|DISCLAIMER)\b", line, flags=re.IGNORECASE):
            violations.append(StyleViolation(
                relative_path, line_number, "generated-disclaimer",
                "remove obsolete generator disclaimers and provenance boilerplate",
            ))
        if re.match(r"^\s*!+\s*\*{5,}\s*!?\s*$", line):
            violations.append(StyleViolation(
                relative_path, line_number, "decorative-separator",
                "decorative star separator comments are forbidden",
            ))
        if re.match(r"^\s*intrinsic\b", code, flags=re.IGNORECASE):
            violations.append(StyleViolation(
                relative_path, line_number, "intrinsic-declaration",
                "do not redeclare standard intrinsic procedures",
            ))
        if re.match(r"^\s*entry\b", code, flags=re.IGNORECASE):
            violations.append(StyleViolation(
                relative_path, line_number, "entry-statement",
                "replace ENTRY with an explicit module procedure",
            ))
        if re.match(r"^\s*(?:[a-z][a-z0-9_]*\s*:\s*)?do\s+\d+\b", code, flags=re.IGNORECASE):
            violations.append(StyleViolation(
                relative_path, line_number, "labeled-do",
                "use structured do ... enddo instead of a labeled loop terminator",
            ))
        if re.search(r"\b(?:go\s*to|goto)\b", code, flags=re.IGNORECASE):
            violations.append(StyleViolation(
                relative_path, line_number, "goto-statement",
                "replace GOTO with structured control flow",
            ))
        if LEGACY_DECLARATION_PATTERN.match(code):
            violations.append(StyleViolation(
                relative_path, line_number, "declaration-separator",
                "variable declarations must use the :: separator",
            ))
        if re.match(r"^\s*parameter\s*\(", code, flags=re.IGNORECASE):
            violations.append(StyleViolation(
                relative_path, line_number, "legacy-parameter",
                "declare constants with a type, parameter, and ::",
            ))
        if re.match(
            r"^\s*character(?:\s*\([^)]*\))?.*::.*"
            r"\b[a-z][a-z0-9_]*\s*\*\s*\d+",
            code,
            flags=re.IGNORECASE,
        ):
            violations.append(StyleViolation(
                relative_path, line_number, "legacy-character-length",
                "put character length in character(len=...), not the entity",
            ))
        if re.search(
            r"\bopen\s*\(\s*(?:unit\s*=\s*)?\d+",
            code,
            flags=re.IGNORECASE,
        ):
            violations.append(StyleViolation(
                relative_path, line_number, "hardcoded-open-unit",
                "OPEN statements must acquire a unit with newunit=",
            ))
        if re.match(r"^\s*\d+\s+format\s*\(", line, flags=re.IGNORECASE):
            violations.append(StyleViolation(
                relative_path, line_number, "labeled-format",
                "put the format specification directly in READ or WRITE",
            ))
        if re.match(r"^\s*return\s*(?:!.*)?$", line, flags=re.IGNORECASE):
            next_index = line_number
            while next_index < len(lines) and (
                not lines[next_index].strip()
                or lines[next_index].lstrip().startswith("!")
            ):
                next_index += 1
            if next_index < len(lines) and re.match(
                r"^\s*(?:end\s*(?:subroutine|function)\b.*|end)\s*$",
                lines[next_index],
                flags=re.IGNORECASE,
            ):
                violations.append(StyleViolation(
                    relative_path, line_number, "terminal-return",
                    "RETURN immediately before a procedure end is redundant",
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
    """Check all free-form sources and fail when any violation is present."""

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
