#!/usr/bin/env python3
"""Replace legacy star-sized intrinsic declarations with standard syntax."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


LEGACY_DECLARATION = re.compile(
    r"^(?P<indent>\s*)(?P<type>double\s+precision|double\s+complex|"
    r"real(?:\s*\([^)]*\))?|integer(?:\s*\([^)]*\))?|"
    r"logical(?:\s*\([^)]*\))?|complex(?:\s*\([^)]*\))?|"
    r"character(?:\s*\([^)]*\))?|type\s*\([^)]*\)|class\s*\([^)]*\))"
    r"(?P<attributes>\s*,[^:]*)?\s+(?P<entities>.+)$",
    flags=re.IGNORECASE,
)


def modernize_line(code: str) -> str:
    """Modernize declaration syntax outside literals and comments on one line."""

    output = []
    segment = []
    quote = ""
    index = 0

    def flush_segment() -> None:
        text = "".join(segment)
        text = re.sub(
            r"\breal\s*\*\s*8\b", "double precision", text,
            flags=re.IGNORECASE,
        )
        text = re.sub(
            r"\breal\s*\*\s*4\b", "real", text,
            flags=re.IGNORECASE,
        )
        text = re.sub(
            r"\bcharacter\s*\*\s*\(\s*\*\s*\)",
            "character(len=*)", text, flags=re.IGNORECASE,
        )
        text = re.sub(
            r"\bcharacter\s*\*\s*(\d+)",
            r"character(len=\1)", text, flags=re.IGNORECASE,
        )
        output.append(text)
        segment.clear()

    while index < len(code):
        character = code[index]
        if quote:
            output.append(character)
            if character == quote:
                if index + 1 < len(code) and code[index + 1] == quote:
                    output.append(code[index + 1])
                    index += 2
                    continue
                quote = ""
        elif character in "'\"":
            flush_segment()
            quote = character
            output.append(character)
        elif character == "!":
            flush_segment()
            output.append(code[index:])
            index = len(code)
            continue
        else:
            segment.append(character)
        index += 1
    flush_segment()
    modernized = "".join(output)
    declaration_code = modernized.split("!", 1)[0].rstrip("\r\n")
    entity_character = re.match(
        r"^(?P<indent>\s*)character\s+(?P<entities>[^:]+?)\s*$",
        declaration_code,
        re.IGNORECASE,
    )
    if entity_character is not None and "::" not in declaration_code:
        grouped = []
        valid = True
        for entity in entity_character.group("entities").split(","):
            match = re.match(
                r"^\s*([a-z][a-z0-9_]*)(?:\s*\*\s*(\d+))?\s*$",
                entity,
                re.IGNORECASE,
            )
            if match is None:
                valid = False
                break
            length = match.group(2) or "1"
            if grouped and grouped[-1][0] == length:
                grouped[-1][1].append(match.group(1))
            else:
                grouped.append((length, [match.group(1)]))
        if valid and grouped:
            indent = entity_character.group("indent")
            replacements = []
            for length, entities in grouped:
                type_spec = (
                    "character" if length == "1"
                    else "character(len={})".format(length)
                )
                replacements.append("{}{} :: {}".format(
                    indent, type_spec, ", ".join(entities)
                ))
            suffix = modernized[len(declaration_code):]
            if suffix.startswith("!"):
                replacements[-1] += " "
            return "\n".join(replacements) + suffix
    declaration = LEGACY_DECLARATION.match(declaration_code)
    if (
        declaration is not None
        and "::" not in declaration_code
        and not re.match(
            r"^(?:function|subroutine)\b",
            declaration.group("entities"),
            re.IGNORECASE,
        )
    ):
        replacement = "{}{}{} :: {}".format(
            declaration.group("indent"),
            declaration.group("type"),
            declaration.group("attributes") or "",
            declaration.group("entities").strip(),
        )
        suffix = modernized[len(declaration_code):]
        if suffix.startswith("!"):
            replacement += " "
        modernized = replacement + suffix
    return modernized


def modernize_code(code: str) -> str:
    """Modernize declaration syntax outside literals and comments."""

    return "".join(modernize_line(line) for line in code.splitlines(keepends=True))


def main() -> int:
    """Modernize requested sources and print a deterministic summary."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        modernized = modernize_code(original)
        changed += modernized != original
        if arguments.write and modernized != original:
            path.write_text(modernized)
    print("fortran_sources_changed: {}".format(changed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
