#!/usr/bin/env python3
"""Inline labeled Fortran FORMAT statements at their READ/WRITE call sites."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


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
FORMAT_START = re.compile(
    r"^\s*(?P<label>\d+)\s+format\s*\(", re.IGNORECASE
)


@dataclass(frozen=True)
class SourceSpan:
    """Represent a half-open source-text span."""

    start: int
    end: int


@dataclass(frozen=True)
class FormatStatement:
    """Represent one labeled FORMAT statement and its descriptor."""

    label: str
    descriptor: str
    span: SourceSpan


@dataclass(frozen=True)
class ConversionResult:
    """Summarize one source conversion."""

    text: str
    formats: int
    references: int
    skipped_labels: tuple[str, ...]


def code_mask(source: str) -> str:
    """Mask literals and comments while preserving source offsets."""

    masked = list(source)
    quote = ""
    index = 0
    while index < len(source):
        character = source[index]
        if quote:
            masked[index] = " "
            if character == quote:
                if index + 1 < len(source) and source[index + 1] == quote:
                    masked[index + 1] = " "
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            quote = character
            masked[index] = " "
        elif character == "!":
            while index < len(source) and source[index] != "\n":
                masked[index] = " "
                index += 1
            continue
        index += 1
    return "".join(masked)


def code_without_comment(line: str) -> str:
    """Return source before a comment marker outside a character literal."""

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


def line_offsets(source: str) -> tuple[list[str], list[int]]:
    """Return physical lines and their absolute starting offsets."""

    lines = source.splitlines(keepends=True)
    offsets = []
    offset = 0
    for line in lines:
        offsets.append(offset)
        offset += len(line)
    return lines, offsets


def has_trailing_continuation(line: str) -> bool:
    """Return whether a physical free-form line ends in an ampersand."""

    return code_without_comment(line).rstrip().endswith("&")


def logical_format(lines: Sequence[str]) -> str:
    """Join the physical lines of one FORMAT statement."""

    parts = []
    for index, line in enumerate(lines):
        code = code_without_comment(line.rstrip("\r\n"))
        code = code.rstrip()
        if code.endswith("&"):
            code = code[:-1]
        code = code.strip()
        if index and code.startswith("&"):
            code = code[1:].lstrip()
        parts.append(code)
    return " ".join(parts)


def find_procedure_spans(source: str) -> tuple[SourceSpan, ...]:
    """Return spans for external and module procedures in a source file."""

    lines, offsets = line_offsets(source)
    spans = []
    start = None
    depth = 0
    for index, line in enumerate(lines):
        if PROCEDURE_START.match(line):
            if start is None:
                start = offsets[index]
            depth += 1
        if start is not None and PROCEDURE_END.match(line):
            depth -= 1
            if depth == 0:
                spans.append(SourceSpan(start, offsets[index] + len(line)))
                start = None
    return tuple(spans)


def find_formats(procedure: str) -> tuple[FormatStatement, ...]:
    """Return every labeled FORMAT statement in a procedure."""

    lines, offsets = line_offsets(procedure)
    formats = []
    index = 0
    while index < len(lines):
        match = FORMAT_START.match(lines[index])
        if match is None:
            index += 1
            continue
        end_index = index
        while end_index + 1 < len(lines) and has_trailing_continuation(
            lines[end_index]
        ):
            end_index += 1
            while (
                end_index + 1 < len(lines)
                and not code_without_comment(lines[end_index]).strip()
            ):
                end_index += 1
        logical = logical_format(lines[index:end_index + 1])
        descriptor_match = re.match(
            r"^\s*\d+\s+format\s*(\(.*\))\s*$",
            logical,
            re.IGNORECASE | re.DOTALL,
        )
        if descriptor_match is None:
            index = end_index + 1
            continue
        formats.append(FormatStatement(
            match.group("label"),
            descriptor_match.group(1),
            SourceSpan(offsets[index], offsets[end_index] + len(lines[end_index])),
        ))
        index = end_index + 1
    return tuple(formats)


def matching_parenthesis(mask: str, opening: int) -> int | None:
    """Find the closing parenthesis matching *opening* in masked code."""

    depth = 0
    for index in range(opening, len(mask)):
        if mask[index] == "(":
            depth += 1
        elif mask[index] == ")":
            depth -= 1
            if depth == 0:
                return index
    return None


def control_items(mask: str, start: int, end: int) -> tuple[SourceSpan, ...]:
    """Split an I/O control list into top-level comma-separated items."""

    items = []
    depth = 0
    item_start = start
    for index in range(start, end):
        if mask[index] == "(":
            depth += 1
        elif mask[index] == ")":
            depth -= 1
        elif mask[index] == "," and depth == 0:
            items.append(SourceSpan(item_start, index))
            item_start = index + 1
    items.append(SourceSpan(item_start, end))
    return tuple(items)


def label_token(item: str, label: str, named: bool) -> SourceSpan | None:
    """Locate *label* when an I/O control item is a format reference."""

    if named:
        match = re.match(
            rf"^\s*fmt\s*=\s*(?P<label>{re.escape(label)})\s*$",
            item,
            re.IGNORECASE,
        )
    else:
        match = re.match(rf"^\s*(?P<label>{re.escape(label)})\s*$", item)
    if match is None:
        return None
    return SourceSpan(match.start("label"), match.end("label"))


def format_expression(descriptor: str, column: int) -> str:
    """Return a safely wrapped character expression for a format descriptor."""

    pieces = []
    remaining = descriptor
    limit = max(16, 116 - column)
    while remaining:
        used = 0
        count = 0
        for character in remaining:
            encoded = 2 if character == "'" else 1
            if count and used + encoded > limit:
                break
            used += encoded
            count += 1
        pieces.append(remaining[:count])
        remaining = remaining[count:]
        limit = 84

    encoded = ["'{}'".format(piece.replace("'", "''")) for piece in pieces]
    if len(encoded) == 1:
        return encoded[0]
    continuation_indent = " " * (column + 4)
    return (" // &\n" + continuation_indent).join(encoded)


def reference_spans(procedure: str, label: str) -> tuple[SourceSpan, ...]:
    """Find direct READ, WRITE, and PRINT references to a format label."""

    mask = code_mask(procedure)
    references = []
    io_pattern = re.compile(r"\b(?:read|write)\s*\(", re.IGNORECASE)
    for match in io_pattern.finditer(mask):
        opening = mask.find("(", match.start())
        closing = matching_parenthesis(mask, opening)
        if closing is None:
            continue
        items = control_items(mask, opening + 1, closing)
        for item_index, item in enumerate(items):
            text = mask[item.start:item.end]
            token = label_token(text, label, named=True)
            if token is None and item_index == 1:
                token = label_token(text, label, named=False)
            if token is not None:
                references.append(SourceSpan(
                    item.start + token.start, item.start + token.end
                ))

    for pattern in (
        re.compile(rf"\bprint\s+(?P<label>{re.escape(label)})\s*,", re.IGNORECASE),
        re.compile(rf"\bread\s+(?P<label>{re.escape(label)})\s*,", re.IGNORECASE),
    ):
        for match in pattern.finditer(mask):
            references.append(SourceSpan(match.start("label"), match.end("label")))
    return tuple(sorted(set(references), key=lambda span: span.start))


def has_assigned_format(procedure: str, label: str) -> bool:
    """Return whether *label* participates in legacy assigned-format I/O."""

    mask = code_mask(procedure)
    return bool(re.search(
        rf"\bassign\s+{re.escape(label)}\s+to\b", mask, re.IGNORECASE
    ))


def apply_replacements(source: str, replacements: Iterable[tuple[SourceSpan, str]]) -> str:
    """Apply nonoverlapping replacements from the end of the source backward."""

    result = source
    ordered = sorted(replacements, key=lambda item: item[0].start, reverse=True)
    for span, replacement in ordered:
        result = result[:span.start] + replacement + result[span.end:]
    return result


def convert_procedure(procedure: str) -> ConversionResult:
    """Inline all directly referenced formats in one procedure."""

    formats = find_formats(procedure)
    replacements = []
    reference_count = 0
    converted_count = 0
    skipped = []
    for statement in formats:
        if has_assigned_format(procedure, statement.label):
            skipped.append(statement.label)
            continue
        references = reference_spans(procedure, statement.label)
        for reference in references:
            line_start = procedure.rfind("\n", 0, reference.start) + 1
            column = reference.start - line_start
            replacements.append((
                reference,
                format_expression(statement.descriptor, column),
            ))
        replacements.append((statement.span, ""))
        reference_count += len(references)
        converted_count += 1
    return ConversionResult(
        apply_replacements(procedure, replacements),
        converted_count,
        reference_count,
        tuple(skipped),
    )


def convert_source(source: str) -> ConversionResult:
    """Inline labeled formats in every procedure in a source file."""

    replacements = []
    format_count = 0
    reference_count = 0
    skipped = []
    for span in find_procedure_spans(source):
        result = convert_procedure(source[span.start:span.end])
        replacements.append((span, result.text))
        format_count += result.formats
        reference_count += result.references
        skipped.extend(result.skipped_labels)
    return ConversionResult(
        apply_replacements(source, replacements),
        format_count,
        reference_count,
        tuple(skipped),
    )


def main() -> int:
    """Convert the requested sources and print a migration summary."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    formats = 0
    references = 0
    skipped = []
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        result = convert_source(original)
        changed += result.text != original
        formats += result.formats
        references += result.references
        skipped.extend("{}:{}".format(path, label) for label in result.skipped_labels)
        if arguments.write and result.text != original:
            path.write_text(result.text)
    print("fortran_sources_changed: {}".format(changed))
    print("format_statements_inlined: {}".format(formats))
    print("format_references_rewritten: {}".format(references))
    print("assigned_formats_skipped: {}".format(len(skipped)))
    for item in skipped:
        print("skipped: {}".format(item))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
