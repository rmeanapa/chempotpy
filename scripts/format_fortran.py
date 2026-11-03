#!/usr/bin/env python3
"""Apply CHEMPOTPY's four-space structural indentation to free-form Fortran."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence


UNIT_START = re.compile(
    r"^(?:module\s+(?!procedure\b)|program\b|block\s+data\b|"
    r"(?:(?:recursive|pure|elemental|impure|module)\s+)*"
    r"(?:(?:double\s+precision|double\s+complex|"
    r"real(?:\s*\([^)]*\)|\*\d+)?|integer(?:\s*\([^)]*\)|\*\d+)?|"
    r"logical(?:\s*\([^)]*\)|\*\d+)?|complex(?:\s*\([^)]*\)|\*\d+)?|"
    r"character(?:\s*\([^)]*\)|\*\d+)?)\s+)?"
    r"(?:subroutine|function)\b)",
    flags=re.IGNORECASE,
)
UNIT_END = re.compile(
    r"^end\s*(?:module|program|block\s+data|subroutine|function)?\b",
    flags=re.IGNORECASE,
)
PROCEDURE_DECLARATION = re.compile(
    r"^(?P<indent>\s*)(?!end\b)"
    r"(?:(?:recursive|pure|elemental|impure|module)\s+)*"
    r"(?:(?:double\s+precision|double\s+complex|"
    r"real(?:\s*\([^)]*\)|\*\d+)?|integer(?:\s*\([^)]*\)|\*\d+)?|"
    r"logical(?:\s*\([^)]*\)|\*\d+)?|complex(?:\s*\([^)]*\)|\*\d+)?|"
    r"character(?:\s*\([^)]*\)|\*\d+)?)\s+)?"
    r"(?P<kind>subroutine|function)\s+(?P<name>[a-z][a-z0-9_]*)"
    r"(?P<arguments>\s*\()?",
    flags=re.IGNORECASE,
)
NUMBERED_LINE = re.compile(r"^(\d+)\s+(.*)$")
DECLARATION = re.compile(
    r"^(?P<indent>\s*)(?P<type>double\s+precision|double\s+complex|"
    r"real(?:\s*\([^)]*\)|\*\d+)?|integer(?:\s*\([^)]*\)|\*\d+)?|"
    r"logical(?:\s*\([^)]*\)|\*\d+)?|complex(?:\s*\([^)]*\)|\*\d+)?|"
    r"character(?:\s*\([^)]*\)|\*\d+)?)"
    r"(?P<attributes>\s*,[^:]*)?\s*::\s*(?P<entities>.*)$",
    flags=re.IGNORECASE,
)
ASSIGNMENT = re.compile(
    r"^(?P<indent>\s*)(?P<target>[a-z][a-z0-9_]*(?:%[a-z][a-z0-9_]*)*"
    r"(?:\([^=]*\))?)\s*=\s*(?P<expression>[^=>].*)$",
    flags=re.IGNORECASE,
)
PARAMETER_ASSIGNMENT = re.compile(
    r"^(?P<indent>\s*)parameter\s*\(\s*"
    r"(?P<name>[a-z][a-z0-9_]*)\s*=\s*(?P<value>[^,]+?)\s*\)\s*$",
    flags=re.IGNORECASE,
)
LABELED_DO = re.compile(
    r"^(?:[a-z][a-z0-9_]*\s*:\s*)?do\s+(\d+)\s+",
    flags=re.IGNORECASE,
)


@dataclass
class Block:
    """Represent one indentation-producing Fortran construct."""

    kind: str
    end_label: str | None = None


def code_without_comment(line: str) -> str:
    """Return code before the first comment marker outside a literal."""

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


def lowercase_fortran_code(source_text: str) -> str:
    """Lowercase Fortran code while preserving literals and comments exactly."""

    output = []
    quote = ""
    comment = False
    preprocessor = False
    line_start = True
    index = 0
    while index < len(source_text):
        character = source_text[index]
        if character == "\n":
            output.append(character)
            comment = False
            preprocessor = False
            line_start = True
            index += 1
            continue
        if line_start and character in " \t":
            output.append(character)
            index += 1
            continue
        if line_start:
            preprocessor = character == "#"
            line_start = False
        if comment or preprocessor:
            output.append(character)
        elif quote:
            output.append(character)
            if character == quote:
                if index + 1 < len(source_text) and source_text[index + 1] == quote:
                    output.append(source_text[index + 1])
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            quote = character
            output.append(character)
        elif character == "!":
            comment = True
            output.append(character)
        else:
            output.append(character.lower())
        index += 1
    return "".join(output)


def normalized_code(line: str) -> tuple[str | None, str]:
    """Return an optional statement label and normalized lowercase code."""

    stripped = code_without_comment(line).strip()
    label_match = NUMBERED_LINE.match(stripped)
    if label_match is None:
        return None, stripped.lower()
    return label_match.group(1), label_match.group(2).strip().lower()


def close_kind(code: str) -> str | None:
    """Return the block kind closed by *code*."""

    compact = re.sub(r"\s+", " ", code)
    endings = (
        (r"^end\s*do\b", "do"),
        (r"^end\s*if\b", "if"),
        (r"^end\s*select\b", "select"),
        (r"^end\s*where\b", "where"),
        (r"^end\s*forall\b", "forall"),
        (r"^end\s*associate\b", "associate"),
        (r"^end\s*block\b", "block"),
        (r"^end\s*critical\b", "critical"),
        (r"^end\s*interface\b", "interface"),
        (r"^end\s*type\b", "type"),
        (r"^end\s*enum\b", "enum"),
        (r"^end\s*team\b", "team"),
    )
    for pattern, kind in endings:
        if re.match(pattern, compact, re.IGNORECASE):
            return kind
    if UNIT_END.match(compact):
        return "unit"
    return None


def branch_statement(code: str) -> bool:
    """Return whether *code* aligns with its owning block header."""

    return bool(re.match(
        r"^(?:else(?:\s*if)?\b|elseif\b|elsewhere\b|"
        r"case\b|class\s+(?:is|default)\b|type\s+is\b|"
        r"rank\s+(?:default|\()|end\s*select\b)",
        code,
        re.IGNORECASE,
    ))


def opening_block(code: str) -> Block | None:
    """Return the block opened by one complete statement, if any."""

    stripped = re.sub(r"^[a-z][a-z0-9_]*\s*:\s*", "", code, flags=re.IGNORECASE)
    label_match = LABELED_DO.match(code)
    if label_match is not None:
        return Block("do", label_match.group(1))
    if re.match(r"^do(?:\s|$)", stripped, re.IGNORECASE):
        return Block("do")
    if re.match(r"^if\s*\(.*\)\s*then\b", stripped, re.IGNORECASE):
        return Block("if")
    if re.match(r"^select\s+(?:case|type|rank)\b", stripped, re.IGNORECASE):
        return Block("select")
    if re.match(r"^associate\s*\(", stripped, re.IGNORECASE):
        return Block("associate")
    if re.match(r"^critical(?:\s|$)", stripped, re.IGNORECASE):
        return Block("critical")
    if re.match(r"^(?:abstract\s+)?interface(?:\s|$)", stripped, re.IGNORECASE):
        return Block("interface")
    if re.match(r"^enum\s*,", stripped, re.IGNORECASE):
        return Block("enum")
    if re.match(r"^(?:change\s+)?team\s*\(", stripped, re.IGNORECASE):
        return Block("team")
    if re.match(r"^block\s*$", stripped, re.IGNORECASE):
        return Block("block")
    if re.match(
        r"^type\s*(?:,.*)?::\s*[a-z][a-z0-9_]*\s*$",
        stripped,
        re.IGNORECASE,
    ):
        return Block("type")
    if re.match(r"^where\s*\(.*\)\s*$", stripped, re.IGNORECASE):
        return Block("where")
    if re.match(r"^forall\s*\(.*\)\s*$", stripped, re.IGNORECASE):
        return Block("forall")
    if UNIT_START.match(stripped):
        return Block("unit")
    return None


def pop_block(stack: List[Block], kind: str) -> None:
    """Close the nearest matching block while tolerating legacy constructs."""

    for index in range(len(stack) - 1, -1, -1):
        if stack[index].kind == kind:
            del stack[index:]
            return


def matching_argument_close(
    lines: Sequence[str], start: int, name_end: int
) -> tuple[int, int] | None:
    """Locate the closing parenthesis of a procedure dummy argument list."""

    opening = lines[start].find("(", name_end)
    if opening < 0:
        return None
    depth = 0
    quote = ""
    for line_index in range(start, len(lines)):
        code = code_without_comment(lines[line_index])
        column = opening if line_index == start else 0
        while column < len(code):
            character = code[column]
            if quote:
                if character == quote:
                    quote = ""
            elif character in "'\"":
                quote = character
            elif character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0:
                    return line_index, column
            column += 1
    return None


def normalize_procedures(source_text: str) -> str:
    """Add procedure briefs and normalize dummy-argument spacing."""

    lines = source_text.splitlines()
    index = 0
    while index < len(lines):
        match = PROCEDURE_DECLARATION.match(lines[index])
        if match is None:
            index += 1
            continue
        previous = index - 1
        while previous >= 0 and not lines[previous].strip():
            previous -= 1
        if previous < 0 or not lines[previous].lstrip().startswith("!> @brief"):
            lines.insert(
                index,
                match.group("indent")
                + "!> @brief Evaluate {}.".format(match.group("name")),
            )
            index += 1
            match = PROCEDURE_DECLARATION.match(lines[index])

        if match.group("arguments") is None:
            index += 1
            continue
        opening = lines[index].find("(", match.start("name"))
        if lines[index][opening + 1:opening + 2] != " ":
            lines[index] = lines[index][:opening + 1] + " " + lines[index][opening + 1:]
        close = matching_argument_close(lines, index, match.end("name"))
        if close is None:
            index += 1
            continue
        close_line, _ = close
        for header_index in range(index, close_line + 1):
            lines[header_index] = re.sub(r",\s*", ", ", lines[header_index])
        close = matching_argument_close(lines, index, match.end("name"))
        if close is not None:
            close_line, close_column = close
            if lines[close_line][close_column - 1:close_column] != " ":
                lines[close_line] = (
                    lines[close_line][:close_column]
                    + " "
                    + lines[close_line][close_column:]
                )
        index = close_line + 1
    return "\n".join(lines) + "\n"


def align_declarations(lines: Sequence[str]) -> List[str]:
    """Align declaration separators in consecutive declaration groups."""

    aligned = list(lines)
    index = 0
    while index < len(aligned):
        first = DECLARATION.match(aligned[index])
        if first is None:
            index += 1
            continue
        group = [(index, first)]
        next_index = index + 1
        while next_index < len(aligned):
            match = DECLARATION.match(aligned[next_index])
            if match is None or match.group("indent") != first.group("indent"):
                break
            group.append((next_index, match))
            next_index += 1
        prefixes = {
            line_index: (
                match.group("type") + (match.group("attributes") or "")
            ).rstrip()
            for line_index, match in group
        }
        width = max(len(prefix) for prefix in prefixes.values())
        for line_index, match in group:
            prefix = prefixes[line_index]
            aligned[line_index] = "{}{}{} :: {}".format(
                match.group("indent"),
                prefix,
                " " * (width-len(prefix)),
                match.group("entities"),
            )
        index = next_index
    return aligned


def declaration_end(lines: Sequence[str], index: int) -> int:
    """Return the first line after one possibly continued declaration."""

    end = index + 1
    while (
        end < len(lines)
        and code_without_comment(lines[end - 1]).rstrip().endswith("&")
    ):
        end += 1
    return end


def declaration_rank(declaration) -> int:
    """Order constants, dummy arguments, and local variables."""

    attributes = declaration.group("attributes") or ""
    if re.search(r"(?:^|,)\s*parameter\s*(?:,|$)", attributes, re.I):
        return 0
    if re.search(r"(?:^|,)\s*intent\s*\(", attributes, re.I):
        return 1
    return 2


def order_initial_declarations(lines: Sequence[str]) -> List[str]:
    """Place intended dummy arguments before local variable declarations."""

    ordered = list(lines)
    index = 0
    while index < len(ordered):
        if not re.match(r"^\s*implicit\s+none\b", ordered[index], re.I):
            index += 1
            continue

        first = index + 1
        while first < len(ordered) and (
            not ordered[first].strip()
            or ordered[first].lstrip().startswith("!")
        ):
            first += 1
        if first >= len(ordered) or DECLARATION.match(ordered[first]) is None:
            index += 1
            continue

        entries = []
        pending = []
        cursor = first
        while cursor < len(ordered):
            if (
                not ordered[cursor].strip()
                or ordered[cursor].lstrip().startswith("!")
            ):
                pending.append(ordered[cursor])
                cursor += 1
                continue
            declaration = DECLARATION.match(ordered[cursor])
            if declaration is None:
                break
            end = declaration_end(ordered, cursor)
            comments = [
                line for line in pending if line.lstrip().startswith("!")
            ]
            entries.append((
                declaration_rank(declaration),
                comments + ordered[cursor:end],
            ))
            pending = []
            cursor = end

        ranks = [rank for rank, _ in entries]
        if 1 not in ranks or ranks == sorted(ranks):
            index = cursor
            continue

        replacement = []
        previous_rank = None
        for rank, chunk in sorted(entries, key=lambda item: item[0]):
            if previous_rank is not None and rank != previous_rank:
                replacement.append("")
            replacement.extend(chunk)
            previous_rank = rank
        if pending:
            if replacement and replacement[-1].strip():
                replacement.append("")
            replacement.extend(
                line for line in pending if line.strip()
            )
        ordered[first:cursor] = replacement
        index = first + len(replacement)
    return ordered


def collapse_blank_lines(lines: Sequence[str]) -> List[str]:
    """Collapse redundant runs while retaining one visual separator."""

    collapsed = []
    for line in lines:
        if not line.strip() and (not collapsed or not collapsed[-1].strip()):
            continue
        collapsed.append(line)
    while collapsed and not collapsed[-1].strip():
        collapsed.pop()
    return collapsed


def remove_terminal_returns(lines: Sequence[str]) -> List[str]:
    """Remove unlabeled RETURN statements immediately before procedure ends."""

    redundant = set()
    for index, line in enumerate(lines):
        if not re.match(r"^\s*return\s*(?:!.*)?$", line, re.IGNORECASE):
            continue
        next_index = index + 1
        while next_index < len(lines) and (
            not lines[next_index].strip()
            or lines[next_index].lstrip().startswith("!")
        ):
            next_index += 1
        if next_index < len(lines) and re.match(
            r"^\s*(?:end\s*(?:subroutine|function)\b.*|end)\s*$",
            lines[next_index],
            re.IGNORECASE,
        ):
            redundant.add(index)
    return [line for index, line in enumerate(lines) if index not in redundant]


def repair_assignment_text(line: str) -> str:
    """Repair legacy token splits in one complete assignment statement."""

    output = []
    code = []
    quote = ""
    index = 0

    def flush_code() -> None:
        text = "".join(code)
        text = re.sub(
            r"(?<=[a-z0-9_]) +(?=[a-z0-9_])",
            "",
            text,
            flags=re.IGNORECASE,
        )
        text = re.sub(r"(?<=\d) +(?=\.\d)", "", text)
        text = re.sub(r"\* +\*", "**", text)
        text = re.sub(r"/ +/", "//", text)
        text = re.sub(r"([<>=/]) +=", r"\1=", text)
        text = re.sub(r": +:", "::", text)
        output.append(text)
        code.clear()

    while index < len(line):
        character = line[index]
        if quote:
            output.append(character)
            if character == quote:
                if index + 1 < len(line) and line[index + 1] == quote:
                    output.append(line[index + 1])
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            flush_code()
            quote = character
            output.append(character)
        elif character == "!":
            flush_code()
            output.append(line[index:])
            index = len(line)
            continue
        else:
            code.append(character)
        index += 1
    flush_code()
    return "".join(output)


def assignment_breaks(expression: str) -> List[int]:
    """Return safe source-continuation offsets in an assignment expression."""

    breaks = []
    index = 0
    quote = ""
    while index < len(expression):
        character = expression[index]
        if quote:
            if character == quote:
                if index + 1 < len(expression) and expression[index + 1] == quote:
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            quote = character
        elif character in ",":
            breaks.append(index + 1)
        elif character in "+-":
            exponent_sign = (
                index > 0
                and expression[index - 1] in "dDeE"
                and index + 1 < len(expression)
                and expression[index + 1].isdigit()
            )
            if not exponent_sign:
                breaks.append(index + 1)
        elif character in "*/":
            if index + 1 < len(expression) and expression[index + 1] == character:
                breaks.append(index + 2)
                index += 1
            else:
                breaks.append(index + 1)
        index += 1
    return breaks


def wrap_assignment(
    statement: str, indentation: str, maximum_width: int
) -> List[str]:
    """Wrap one complete assignment deterministically at safe operators."""

    content = statement.strip()
    pieces = []
    first = True
    while True:
        prefix_width = len(indentation) + (0 if first else 6)
        final_width = maximum_width-prefix_width
        if len(content) <= final_width:
            pieces.append(content)
            break
        available = final_width-2  # Reserve the trailing " &".
        candidates = [
            position for position in assignment_breaks(content)
            if position <= available
        ]
        if candidates:
            split = candidates[-1]
        else:
            later = assignment_breaks(content)
            if not later:
                pieces.append(content)
                break
            split = later[0]
        if split >= len(content):
            pieces.append(content)
            break
        pieces.append(content[:split].rstrip())
        content = content[split:].lstrip()
        first = False

    wrapped = []
    for piece_index, piece in enumerate(pieces):
        prefix = indentation if piece_index == 0 else indentation + "    & "
        suffix = " &" if piece_index < len(pieces)-1 else ""
        wrapped.append(prefix + piece + suffix)
    return wrapped


def combine_short_assignment_continuations(
    lines: Sequence[str], maximum_width: int = 120
) -> List[str]:
    """Combine unnecessarily short physical lines in continued assignments."""

    combined: List[str] = []
    index = 0
    while index < len(lines):
        first = lines[index]
        if (
            ASSIGNMENT.match(first) is None
            or not code_without_comment(first).rstrip().endswith("&")
            or code_without_comment(first) != first
        ):
            combined.append(first)
            index += 1
            continue

        statement = [first]
        next_index = index + 1
        valid = True
        while code_without_comment(statement[-1]).rstrip().endswith("&"):
            if next_index >= len(lines):
                valid = False
                break
            continuation = lines[next_index]
            if (
                not continuation.strip()
                or continuation.lstrip().startswith("!")
                or code_without_comment(continuation) != continuation
            ):
                valid = False
                break
            statement.append(continuation)
            next_index += 1
        if not valid or len(statement) == 1:
            combined.append(first)
            index += 1
            continue

        segments = []
        for segment_index, line in enumerate(statement):
            content = line.strip()
            if segment_index:
                content = re.sub(r"^&", "", content)
            if segment_index < len(statement)-1:
                content = re.sub(r"&\s*$", "", content)
            segments.append(content)

        logical = segments[0].rstrip()
        for segment in segments[1:]:
            left = logical.rstrip()
            right = segment.lstrip()
            last_character = left[-1:] or " "
            first_character = right[:1] or " "
            token_continuation = (
                (
                    re.match(r"[a-z0-9_]", last_character, re.IGNORECASE)
                    and re.match(r"[a-z0-9_]", first_character, re.IGNORECASE)
                )
                or last_character + first_character
                in {"**", "//", "<=", ">=", "==", "/=", "::"}
            )
            separator = "" if token_continuation else " "
            logical = left + separator + right

        indentation = re.match(r"^\s*", first).group(0)
        logical = repair_assignment_text(logical)
        combined.extend(wrap_assignment(logical, indentation, maximum_width))
        index = next_index
    return combined


def repair_split_assignment_tokens(
    lines: Sequence[str], maximum_width: int = 120
) -> List[str]:
    """Repair split tokens and wrap long assignments at safe operators."""

    repaired = []
    for line in lines:
        if ASSIGNMENT.match(line) is None:
            repaired.append(line)
            continue
        statement = repair_assignment_text(line)
        if (
            len(statement) <= maximum_width
            or code_without_comment(statement) != statement
        ):
            repaired.append(statement)
            continue
        indentation = re.match(r"^\s*", statement).group(0)
        repaired.extend(wrap_assignment(
            statement, indentation, maximum_width
        ))
    return repaired


def align_parameters(lines: Sequence[str]) -> List[str]:
    """Align equals signs in consecutive single-name PARAMETER statements."""

    aligned = list(lines)
    index = 0
    while index < len(aligned):
        first = PARAMETER_ASSIGNMENT.match(aligned[index])
        if first is None:
            index += 1
            continue
        group = [(index, first)]
        next_index = index + 1
        while next_index < len(aligned):
            match = PARAMETER_ASSIGNMENT.match(aligned[next_index])
            if match is None or match.group("indent") != first.group("indent"):
                break
            group.append((next_index, match))
            next_index += 1
        width = max(len(match.group("name")) for _, match in group)
        for line_index, match in group:
            name = match.group("name")
            aligned[line_index] = "{}parameter ({}{} = {})".format(
                match.group("indent"),
                name,
                " " * (width-len(name)),
                match.group("value").strip(),
            )
        index = next_index
    return aligned


def typed_parameter_assignment(line: str):
    """Return the top-level assignment in one typed parameter declaration."""

    declaration = DECLARATION.match(line)
    if declaration is None or not re.search(
        r"(?:^|,)\s*parameter\s*(?:,|$)",
        declaration.group("attributes") or "",
        re.IGNORECASE,
    ):
        return None
    entities = declaration.group("entities")
    depth = 0
    quote = ""
    index = 0
    while index < len(entities):
        character = entities[index]
        if quote:
            if character == quote:
                if index + 1 < len(entities) and entities[index + 1] == quote:
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            quote = character
        elif character in "([":
            depth += 1
        elif character in ")]":
            depth -= 1
        elif character == "," and depth == 0:
            return None
        elif character == "=" and depth == 0:
            prefix = line[:declaration.start("entities")]
            return prefix + entities[:index].rstrip(), entities[index + 1:].strip()
        index += 1
    return None


def align_typed_parameters(lines: Sequence[str]) -> List[str]:
    """Align equals signs in consecutive typed parameter declarations."""

    aligned = list(lines)
    index = 0
    while index < len(aligned):
        first = typed_parameter_assignment(aligned[index])
        if first is None:
            index += 1
            continue
        indentation = re.match(r"^\s*", aligned[index]).group(0)
        group = [(index, first)]
        next_index = index + 1
        while next_index < len(aligned):
            assignment = typed_parameter_assignment(aligned[next_index])
            if (
                assignment is None
                or re.match(r"^\s*", aligned[next_index]).group(0) != indentation
            ):
                break
            group.append((next_index, assignment))
            next_index += 1
        width = max(len(lhs) for _, (lhs, _) in group)
        for line_index, (lhs, value) in group:
            aligned[line_index] = "{}{} = {}".format(
                lhs,
                " " * (width-len(lhs)),
                value,
            )
        index = next_index
    return aligned


def align_assignments(lines: Sequence[str]) -> List[str]:
    """Align equals signs in consecutive simple assignment groups."""

    aligned = list(lines)
    index = 0
    while index < len(aligned):
        first = ASSIGNMENT.match(aligned[index])
        if first is None or aligned[index].rstrip().endswith("&"):
            index += 1
            continue
        group = [(index, first)]
        next_index = index + 1
        while next_index < len(aligned):
            match = ASSIGNMENT.match(aligned[next_index])
            if (
                match is None
                or match.group("indent") != first.group("indent")
                or aligned[next_index].rstrip().endswith("&")
            ):
                break
            group.append((next_index, match))
            next_index += 1
        width = max(len(match.group("target")) for _, match in group)
        for line_index, match in group:
            aligned[line_index] = "{}{}{} = {}".format(
                match.group("indent"),
                match.group("target"),
                " " * (width-len(match.group("target"))),
                match.group("expression").strip(),
            )
        index = next_index
    return aligned


def format_source_text(source_text: str) -> str:
    """Return structurally indented free-form Fortran source."""

    source_text = lowercase_fortran_code(source_text)
    source_text = source_text.expandtabs(4)
    source_text = normalize_procedures(source_text)
    source_text = "\n".join(
        remove_terminal_returns(source_text.splitlines())
    ) + "\n"
    source_text = "\n".join(
        combine_short_assignment_continuations(source_text.splitlines())
    ) + "\n"
    source_text = "\n".join(
        repair_split_assignment_tokens(source_text.splitlines())
    ) + "\n"
    output: List[str] = []
    stack: List[Block] = []
    continuing = False
    for original in source_text.splitlines():
        if not original.strip():
            output.append("")
            continue
        stripped = original.strip()
        if stripped.startswith("#"):
            output.append(stripped)
            continue
        if stripped.startswith("!"):
            if continuing:
                output.append(" " * (4 * (len(stack) + 1)) + stripped)
            elif stripped.startswith("!> @brief"):
                output.append(" " * (4 * len(stack)) + stripped)
            else:
                output.append(" " * (4 * len(stack)) + stripped)
            continue

        label, code = normalized_code(original)
        if continuing:
            indent_level = len(stack) + 1
            content = stripped
            output.append(" " * (4 * indent_level) + content)
            continuing = code_without_comment(original).rstrip().endswith("&")
            continue

        if label is not None:
            while stack and stack[-1].end_label == label:
                stack.pop()

        ending = close_kind(code)
        if ending is not None:
            pop_block(stack, ending)

        indent_level = len(stack)
        if branch_statement(code) and ending is None:
            indent_level = max(0, indent_level - 1)

        if label is not None:
            content = "{} {}".format(label, NUMBERED_LINE.match(stripped).group(2).strip())
            output.append(" " * (4 * indent_level) + content)
        else:
            output.append(" " * (4 * indent_level) + stripped)

        continuing = code_without_comment(original).rstrip().endswith("&")
        if ending is None:
            opened = opening_block(code)
            if opened is not None:
                stack.append(opened)

    output = order_initial_declarations(output)
    output = align_declarations(output)
    output = align_assignments(output)
    output = align_parameters(output)
    output = align_typed_parameters(output)
    output = collapse_blank_lines(output)
    return "\n".join(output) + "\n"


def format_source(path: Path, write: bool) -> bool:
    """Format one source and optionally write it."""

    original = path.read_text(errors="strict")
    formatted = format_source_text(original)
    changed = formatted != original
    if write and changed:
        path.write_text(formatted)
    return changed


def main() -> int:
    """Format requested free-form sources."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = sum(format_source(path, arguments.write) for path in arguments.sources)
    print("fortran_sources_changed: {}".format(changed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
