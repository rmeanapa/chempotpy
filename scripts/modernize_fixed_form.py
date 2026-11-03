#!/usr/bin/env python3
"""Convert CHEMPOTPY fixed-form sources to behavior-preserving free form."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


PROCEDURE_PATTERN = re.compile(
    r"^(?P<indent>\s*)(?:(?:recursive|pure|elemental)\s+)*"
    r"(?:(?:double\s+precision|real(?:\*\d+)?|integer|logical|complex|"
    r"character(?:\([^)]*\)|\*\d+)?)\s+)?"
    r"(?P<kind>subroutine|function)\s+(?P<name>[a-z][a-z0-9_]*)"
    r"(?P<arguments>\s*\()?",
    flags=re.IGNORECASE,
)
MODULE_PATTERN = re.compile(
    r"^\s*module\s+(?!procedure\b)", flags=re.IGNORECASE
)


@dataclass
class Procedure:
    """Track the external procedure currently being converted."""

    kind: str
    name: str


def comment_index(line: str) -> int:
    """Return the first free-form comment marker outside a string literal."""

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
            return index
        index += 1
    return -1


def has_unclosed_literal(line: str) -> bool:
    """Return whether the code portion ends inside a character literal."""

    stop = comment_index(line)
    code = line if stop < 0 else line[:stop]
    quote = ""
    index = 0
    while index < len(code):
        character = code[index]
        if quote:
            if character == quote:
                if index + 1 < len(code) and code[index + 1] == quote:
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            quote = character
        index += 1
    return bool(quote)


def continue_line(
    line: str, inside_literal: bool, join_boundary_token: bool
) -> str:
    """Append a free-form continuation marker without swallowing comments."""

    if inside_literal or join_boundary_token:
        return line.rstrip() + "&"
    stop = comment_index(line)
    if stop < 0:
        return line.rstrip() + " &"
    return line[:stop].rstrip() + " &  " + line[stop:]


def procedure_argument_spacing(line: str, declaration: bool) -> str:
    """Apply spaces inside procedure declarations and after their commas."""

    if not declaration:
        return line
    line = re.sub(r",\s*", ", ", line)
    match = PROCEDURE_PATTERN.match(line)
    if match is not None and match.group("arguments"):
        opening = line.find("(", match.start("name"))
        if line[opening:opening + 2] != "( ":
            line = line[:opening + 1] + " " + line[opening + 1:]
    closing = line.rfind(")")
    if closing >= 0 and line[closing - 1:closing] != " ":
        line = line[:closing] + " " + line[closing:]
    return line


def normalize_fixed_form_blanks(line: str) -> str:
    """Remove blanks that fixed form permits inside numbers and operators."""

    stop = comment_index(line)
    code = line if stop < 0 else line[:stop]
    comment = "" if stop < 0 else line[stop:]
    pieces = re.split(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")", code)
    operators = (
        "eqv", "neqv", "and", "not", "or", "eq", "ne", "lt", "le",
        "gt", "ge", "true", "false",
    )
    operator_pattern = re.compile(
        r"\.\s*(" + "|".join(
            r"\s*".join(name) for name in operators
        ) + r")\s*\.",
        flags=re.IGNORECASE,
    )
    for index in range(0, len(pieces), 2):
        piece = pieces[index]
        previous = None
        while previous != piece:
            previous = piece
            piece = re.sub(r"(?<=\d)\s+(?=\d)", "", piece)
            piece = re.sub(
                r"(?<=[dDeE])\s+([+-]?)\s*(?=\d)",
                lambda match: match.group(1),
                piece,
            )
        piece = operator_pattern.sub(
            lambda match: "." + re.sub(r"\s+", "", match.group(1)) + ".",
            piece,
        )
        pieces[index] = piece
    return "".join(pieces) + comment


def converted_statement(
    expanded: str,
    continuation: bool,
    literal_continuation: bool,
    join_boundary_token: bool,
) -> str:
    """Convert one fixed-form statement field to a free-form physical line."""

    padded = expanded.ljust(6)
    label = padded[:5].strip()
    statement = padded[6:72].rstrip()
    if not statement:
        return ""
    if continuation:
        if literal_continuation:
            return "        &" + statement
        if join_boundary_token:
            return "        &" + statement.lstrip()
        return "        & " + statement.lstrip()
    if label:
        return label + " " + statement.lstrip()
    match = PROCEDURE_PATTERN.match(statement)
    indentation = 0 if match is not None else 4
    leading = len(statement) - len(statement.lstrip())
    if leading:
        indentation += 4*(leading//4)
    return " "*indentation + statement.lstrip()


def convert_source_text(source_text: str, surface_name: str) -> str:
    """Convert fixed-form text while retaining its statement order and tokens."""

    original_lines = source_text.splitlines()
    contains_module = any(
        MODULE_PATTERN.match(line.expandtabs(8)[6:72])
        for line in original_lines
        if len(line.expandtabs(8)) >= 6
    )
    output: List[str] = []
    if not contains_module:
        output.append(
            "!> @modern CHEMPOTPY modern Fortran style is enforced for this file."
        )
    last_code_index: Optional[int] = None
    procedure: Optional[Procedure] = None
    declaration_continues = False
    previous_ends_token = False
    for original in original_lines:
        expanded = original.expandtabs(8)
        if not expanded.strip():
            output.append("")
            continue
        if expanded[0] in "cC*!":
            output.append("!" + expanded[1:72].rstrip())
            continue
        padded = expanded.ljust(6)
        continuation = padded[5] not in " 0"
        literal_continuation = (
            continuation
            and last_code_index is not None
            and has_unclosed_literal(output[last_code_index])
        )
        statement = padded[6:72].rstrip()
        current_begins_token = bool(
            statement.lstrip()
            and re.match(r"[a-z0-9_.+-]", statement.lstrip(), re.IGNORECASE)
        )
        join_boundary_token = (
            continuation
            and previous_ends_token
            and current_begins_token
            and not literal_continuation
        )
        if continuation and last_code_index is not None:
            output[last_code_index] = continue_line(
                output[last_code_index],
                literal_continuation,
                join_boundary_token,
            )
        converted = converted_statement(
            expanded,
            continuation,
            literal_continuation,
            join_boundary_token,
        )
        if not converted:
            output.append("")
            continue
        converted = normalize_fixed_form_blanks(converted)
        match = PROCEDURE_PATTERN.match(converted)
        if match is not None:
            procedure = Procedure(
                kind=match.group("kind").lower(),
                name=match.group("name"),
            )
            output.append(
                "!> @brief Evaluate {} for the {} potential.".format(
                    procedure.name, surface_name
                )
            )
            declaration_continues = (
                match.group("arguments") is not None and ")" not in converted
            )
            converted = procedure_argument_spacing(converted, True)
        elif declaration_continues:
            converted = procedure_argument_spacing(converted, True)
            declaration_continues = ")" not in converted
        if procedure is not None and converted.strip().lower() in {
            "end", "endsubroutine", "end subroutine", "endfunction",
            "end function",
        }:
            converted = "end {} {}".format(procedure.kind, procedure.name)
            procedure = None
        elif re.match(
            r"^\s*end\s+(?:subroutine|function)\b",
            converted,
            flags=re.IGNORECASE,
        ):
            procedure = None
        output.append(converted)
        last_code_index = len(output) - 1
        previous_ends_token = bool(
            statement
            and re.search(r"[a-z0-9_.+-]$", statement, re.IGNORECASE)
        )
    while output and not output[-1]:
        output.pop()
    return "\n".join(output) + "\n"


def convert_source(source: Path, overwrite: bool = False) -> Path:
    """Write the free-form sibling of one fixed-form source."""

    destination = source.with_suffix(".f90")
    if destination.exists() and not overwrite:
        raise FileExistsError("destination already exists: {}".format(destination))
    destination.write_text(convert_source_text(
        source.read_text(errors="strict"), source.stem
    ))
    return destination


def main() -> int:
    """Convert all requested fixed-form sources without removing originals."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    for source in args.sources:
        destination = convert_source(source, overwrite=args.overwrite)
        print("{} -> {}".format(source, destination))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
