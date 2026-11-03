#!/usr/bin/env python3
"""Convert simple legacy PARAMETER statements to typed declarations."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

try:
    from scripts.format_fortran import DECLARATION
except ModuleNotFoundError:  # Direct execution from the scripts directory.
    from format_fortran import DECLARATION


SIMPLE_PARAMETER = re.compile(
    r"^(?P<indent>\s*)parameter\s*\(\s*"
    r"(?P<name>[a-z][a-z0-9_]*)\s*=\s*(?P<value>[^,]+?)\s*\)\s*$",
    re.IGNORECASE,
)
PARAMETER_START = re.compile(r"^(?P<indent>\s*)parameter\s*\(", re.IGNORECASE)
UNIT_START = re.compile(
    r"^\s*(?!end\b)(?:module\s+(?!procedure\b)|program\b|block\s+data\b|"
    r"(?:(?:recursive|pure|elemental|impure|module)\s+)*"
    r"(?:(?:double\s+precision|double\s+complex|real(?:\s*\([^)]*\))?|"
    r"integer(?:\s*\([^)]*\))?|logical(?:\s*\([^)]*\))?|"
    r"complex(?:\s*\([^)]*\))?|character(?:\s*\([^)]*\))?)\s+)?"
    r"(?:subroutine|function)\b)",
    re.IGNORECASE,
)


def split_entities(entities: str) -> list[str]:
    """Split a declaration entity list at top-level commas."""

    result = []
    start = 0
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
            result.append(entities[start:index].strip())
            start = index + 1
        index += 1
    result.append(entities[start:].strip())
    return result


def entity_name(entity: str) -> str:
    """Return the identifier declared by one entity."""

    match = re.match(r"([a-z][a-z0-9_]*)", entity, re.IGNORECASE)
    return match.group(1).lower() if match else ""


def declaration_at(lines: list[str], index: int):
    """Return one complete continued declaration beginning at index."""

    declaration = DECLARATION.match(lines[index])
    if declaration is None:
        return None
    entity_parts = []
    content = declaration.group("entities").strip()
    continued = content.endswith("&")
    if continued:
        content = content[:-1].rstrip()
    entity_parts.append(content)
    end = index + 1
    while continued and end < len(lines):
        content = lines[end].strip()
        if content.startswith("&"):
            content = content[1:].lstrip()
        continued = content.endswith("&")
        if continued:
            content = content[:-1].rstrip()
        entity_parts.append(content)
        end += 1
    if continued:
        return None
    return declaration, split_entities(" ".join(entity_parts)), end


def render_declaration(declaration, entities: list[str]) -> list[str]:
    """Render a declaration entity list with deterministic comma wrapping."""

    if not entities:
        return []
    prefix = (
        declaration.group("type") + (declaration.group("attributes") or "")
    ).rstrip()
    head = "{}{} :: ".format(declaration.group("indent"), prefix)
    continuation = declaration.group("indent") + "    & "
    rendered = []
    current = head + entities[0]
    for entity in entities[1:]:
        addition = ", " + entity
        if len(current) + len(addition) <= 120:
            current += addition
        else:
            rendered.append(current + ", &")
            current = continuation + entity
    rendered.append(current)
    return rendered


def split_comment(line: str) -> tuple[str, str]:
    """Split a source line at a comment marker outside character literals."""

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
            return line[:index], line[index:]
        index += 1
    return line, ""


def parameter_at(lines: list[str], index: int):
    """Return every assignment in one complete legacy PARAMETER statement."""

    start = PARAMETER_START.match(lines[index])
    if start is None:
        return None
    pieces = []
    comments = []
    depth = 0
    opened = False
    quote = ""
    end = index
    while end < len(lines):
        code, comment = split_comment(lines[end])
        if comment:
            comments.append(comment.strip())
        content = code.strip()
        if end > index and content.startswith("&"):
            content = content[1:].lstrip()
        if content.endswith("&"):
            content = content[:-1].rstrip()
        pieces.append(content)
        for character in content:
            if quote:
                if character == quote:
                    quote = ""
            elif character in "'\"":
                quote = character
            elif character == "(":
                depth += 1
                opened = True
            elif character == ")":
                depth -= 1
        end += 1
        if opened and depth == 0:
            break
    if not opened or depth != 0 or len(comments) > 1:
        return None
    logical = " ".join(pieces)
    opening = logical.find("(")
    closing = logical.rfind(")")
    if opening < 0 or closing <= opening:
        return None
    assignments = []
    for entity in split_entities(logical[opening + 1:closing]):
        match = re.match(
            r"^\s*([a-z][a-z0-9_]*)\s*=\s*(.*?)\s*$",
            entity,
            re.IGNORECASE,
        )
        if match is None:
            return None
        assignments.append((match.group(1).lower(), match.group(2)))
    return start.group("indent"), assignments, end, comments[0] if comments else ""


def migrate_source(source: str) -> tuple[str, int]:
    """Convert legacy parameters whose explicit types are locally available."""

    lines = source.splitlines()
    converted = 0
    unit_start = 0
    index = 0
    while index < len(lines):
        if UNIT_START.match(lines[index]):
            unit_start = index
        parameter = parameter_at(lines, index)
        if parameter is None:
            index += 1
            continue
        indentation, assignments, parameter_end, comment = parameter
        declarations = {}
        candidate = unit_start + 1
        while candidate < index:
            complete = declaration_at(lines, candidate)
            if complete is None:
                candidate += 1
                continue
            match, entities, end = complete
            for entity in entities:
                declarations[entity_name(entity)] = (
                    candidate,
                    end,
                    match,
                    entity,
                    entities,
                )
            candidate = end
        if any(name not in declarations for name, _ in assignments):
            index += 1
            continue

        parameter_lines = []
        names = {name for name, _ in assignments}
        for name, value in assignments:
            _, _, declaration, entity, _ = declarations[name]
            prefix = (
                declaration.group("type")
                + (declaration.group("attributes") or "")
            ).rstrip()
            parameter_lines.append("{}{}, parameter :: {} = {}".format(
                indentation,
                prefix,
                entity.split("=", 1)[0].strip(),
                value.strip(),
            ))
        if comment:
            parameter_lines[-1] += "  " + comment

        declaration_groups = {}
        for name in names:
            start, end, declaration, _, entities = declarations[name]
            declaration_groups[start] = (end, declaration, entities)
        total_delta = 0
        for start in sorted(declaration_groups, reverse=True):
            end, declaration, entities = declaration_groups[start]
            remaining = [
                entity for entity in entities if entity_name(entity) not in names
            ]
            replacement = render_declaration(declaration, remaining)
            lines[start:end] = replacement
            total_delta += len(replacement) - (end-start)
        index += total_delta
        parameter_end += total_delta
        lines[index:parameter_end] = parameter_lines
        converted += len(assignments)
        index += len(parameter_lines)
    return "\n".join(lines) + "\n", converted


def main() -> int:
    """Convert requested sources and print conversion counts."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    converted = 0
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        modernized, count = migrate_source(original)
        changed += modernized != original
        converted += count
        if arguments.write and modernized != original:
            path.write_text(modernized)
    print("fortran_sources_changed: {}".format(changed))
    print("parameters_modernized: {}".format(converted))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
