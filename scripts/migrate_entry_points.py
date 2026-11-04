#!/usr/bin/env python3
"""Replace legacy Fortran ENTRY pairs with module-owned dispatch procedures."""

from __future__ import annotations

import argparse
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

SUBROUTINE_RE = re.compile(
    r"^(?P<indent>\s*)subroutine\s+(?P<name>[a-z]\w*)\s*(?:\((?P<args>[^\n]*)\))?\s*$",
    re.IGNORECASE,
)
END_SUBROUTINE_RE = re.compile(
    r"^\s*(?:end\s+subroutine(?:\s+[a-z]\w*)?|end)\s*$",
    re.IGNORECASE,
)
ENTRY_RE = re.compile(
    r"^\s*entry\s+(?P<name>[a-z]\w*)\s*(?:\((?P<args>[^\n]*)\))?\s*$",
    re.IGNORECASE,
)
RETURN_RE = re.compile(r"^\s*return\s*$", re.IGNORECASE)
DECLARATION_RE = re.compile(
    r"^(?:use\b|implicit\b|integer\b|double\s+precision\b|real\b|logical\b|"
    r"character\b|complex\b|dimension\b|common\b|save\b|data\b|parameter\b|"
    r"equivalence\b|external\b|intrinsic\b|allocatable\b|pointer\b|target\b|"
    r"optional\b|intent\b|namelist\b)",
    re.IGNORECASE,
)
TYPE_RE = re.compile(
    r"^(?P<type>double\s+precision|integer|real(?:\s*\([^)]*\))?|logical|"
    r"character(?:\s*\([^)]*\))?|complex(?:\s*\([^)]*\))?)(?![a-z0-9_])(?P<rest>.*)$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class EntryUnit:
    start: int
    end: int
    entry_index: int
    return_index: int
    executable_index: int
    preserve_return: bool
    initializer_name: str
    initializer_args: tuple[str, ...]
    entry_name: str
    entry_args: tuple[str, ...]


def split_arguments(text: str | None) -> tuple[str, ...]:
    if not text or not text.strip():
        return ()
    return tuple(part.strip().lower() for part in text.split(",") if part.strip())


def normalize_identifier_continuations(lines: list[str]) -> list[str]:
    normalized: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        while line.rstrip().endswith("&") and index + 1 < len(lines):
            continuation = lines[index + 1].lstrip()
            if not re.match(r"&\S", continuation):
                break
            line = line.rstrip()[:-1] + continuation[1:]
            index += 1
        normalized.append(line)
        index += 1
    return normalized


def logical_statements(lines: list[str], start: int, end: int) -> list[tuple[int, int, str]]:
    statements: list[tuple[int, int, str]] = []
    index = start
    while index < end:
        line = lines[index]
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            index += 1
            continue
        first = index
        parts = [stripped]
        separators: list[str] = []
        while parts[-1].endswith("&") and index + 1 < end:
            index += 1
            continuation = lines[index].strip()
            separator = " "
            if continuation.startswith("&"):
                continuation = continuation[1:]
                if continuation and not continuation[0].isspace():
                    separator = ""
                continuation = continuation.lstrip()
            separators.append(separator)
            parts.append(continuation)
        statement = parts[0]
        for separator, part in zip(separators, parts[1:]):
            statement = statement[:-1].rstrip() + separator + part
        statements.append((first, index, statement))
        index += 1
    return statements


def find_units(lines: list[str]) -> list[EntryUnit]:
    units: list[EntryUnit] = []
    index = 0
    while index < len(lines):
        match = SUBROUTINE_RE.match(lines[index])
        if not match or lines[index].lstrip().lower().startswith("end "):
            index += 1
            continue
        start = index
        end = start + 1
        while end < len(lines) and not END_SUBROUTINE_RE.match(lines[end]):
            end += 1
        if end == len(lines):
            raise ValueError(f"unterminated subroutine {match.group('name')}")
        entries = [(i, ENTRY_RE.match(lines[i])) for i in range(start + 1, end) if ENTRY_RE.match(lines[i])]
        if not entries:
            index = end + 1
            continue
        if len(entries) != 1:
            raise ValueError(f"{match.group('name')} has {len(entries)} ENTRY statements; expected one")
        entry_index, entry_match = entries[0]
        assert entry_match is not None
        return_index = next(
            (candidate for candidate in range(entry_index - 1, start, -1) if RETURN_RE.match(lines[candidate])),
            -1,
        )
        if return_index < 0:
            raise ValueError(f"ENTRY {entry_match.group('name')} has no preceding RETURN")
        statements = logical_statements(lines, start + 1, return_index)
        executable_index = -1
        for first, _, statement in statements:
            if not DECLARATION_RE.match(statement):
                executable_index = first
                break
        if executable_index < 0:
            executable_index = return_index
        preserve_return = any(
            stripped
            and not stripped.startswith("!")
            and not re.match(r"^(?:end\s*if|endif)\b", stripped, re.IGNORECASE)
            for stripped in (line.strip() for line in lines[return_index + 1 : entry_index])
        )
        units.append(
            EntryUnit(
                start=start,
                end=end,
                entry_index=entry_index,
                return_index=return_index,
                executable_index=executable_index,
                preserve_return=preserve_return,
                initializer_name=match.group("name").lower(),
                initializer_args=split_arguments(match.group("args")),
                entry_name=entry_match.group("name").lower(),
                entry_args=split_arguments(entry_match.group("args")),
            )
        )
        index = end + 1
    return units


def split_entities(text: str) -> list[str]:
    entities: list[str] = []
    level = 0
    start = 0
    for index, character in enumerate(text):
        if character == "(":
            level += 1
        elif character == ")":
            level -= 1
        elif character == "," and level == 0:
            entities.append(text[start:index].strip())
            start = index + 1
    entities.append(text[start:].strip())
    return [entity for entity in entities if entity]


def argument_declarations(lines: list[str], unit: EntryUnit) -> dict[str, str]:
    arguments = set(unit.initializer_args + unit.entry_args)
    types: dict[str, str] = {}
    shapes: dict[str, str] = {}
    for _, _, statement in logical_statements(lines, unit.start + 1, unit.executable_index):
        clean = statement.strip()
        dimension_match = re.match(r"dimension\s+(.*)$", clean, re.IGNORECASE)
        if dimension_match:
            for entity in split_entities(dimension_match.group(1)):
                match = re.match(r"([a-z]\w*)\s*(\([^=]*\))", entity, re.IGNORECASE)
                if match and match.group(1).lower() in arguments:
                    shapes[match.group(1).lower()] = match.group(2).replace(" ", "")
            continue
        type_match = TYPE_RE.match(clean)
        if not type_match:
            continue
        type_spec = re.sub(r"\s+", " ", type_match.group("type").strip().lower())
        rest = type_match.group("rest")
        entity_text = rest.split("::", 1)[1] if "::" in rest else rest
        for entity in split_entities(entity_text):
            match = re.match(r"([a-z]\w*)\s*(\([^=]*\))?", entity.strip(), re.IGNORECASE)
            if not match:
                continue
            name = match.group(1).lower()
            if name not in arguments:
                continue
            types[name] = type_spec
            if match.group(2):
                shapes[name] = match.group(2).replace(" ", "")
    missing = sorted(arguments - types.keys())
    if missing:
        raise ValueError(f"missing explicit types for arguments: {', '.join(missing)}")
    return {name: f"{types[name]} :: {name}{shapes.get(name, '')}" for name in arguments}


def module_name(path: Path, unit: EntryUnit) -> str:
    base = re.sub(r"[^a-z0-9]+", "_", f"chempotpy_{path.stem}_{unit.entry_name}_legacy_state".lower()).strip("_")
    relative_path = path.resolve().relative_to(REPOSITORY_ROOT)
    digest = hashlib.sha1(relative_path.as_posix().encode()).hexdigest()[:8]
    return f"{base[:53]}_{digest}"


def signature(name: str, arguments: tuple[str, ...]) -> str:
    if arguments:
        return f"subroutine {name}( {', '.join(arguments)} )"
    return f"subroutine {name}"


def wrapper(
    name: str,
    arguments: tuple[str, ...],
    declarations: dict[str, str],
    state_module: str,
    initialize_only: bool,
) -> list[str]:
    result = [f"!> @brief Compatibility interface for {name}.", signature(name, arguments)]
    result.append(f"    use {state_module}, only: dispatch")
    result.append("    implicit none")
    if any(re.search(r"\(\s*dp\s*\)", declarations[argument], re.IGNORECASE) for argument in arguments):
        result.append("    integer, parameter :: dp = kind(0.d0)")
    for argument in arguments:
        result.append(f"    {declarations[argument]}")
    size_argument = next((argument for argument in reversed(arguments) if argument in {"n", "nt"}), None)
    keywords: list[str] = []
    for argument in arguments:
        actual = argument
        shape_match = re.search(r"\(([^)]*)\)\s*$", declarations[argument])
        if shape_match and "*" in shape_match.group(1):
            if size_argument is None:
                raise ValueError(f"cannot section assumed-size argument {argument} without n or nt")
            dimensions = [dimension.strip() for dimension in shape_match.group(1).split(",")]
            sections = [f"1:{size_argument}" if dimension == "*" else f"1:{dimension}" for dimension in dimensions]
            actual = f"{argument}({','.join(sections)})"
        keywords.append(f"{argument}={actual}")
    keywords.append(f"initialize_only={'.true.' if initialize_only else '.false.'}")
    result.append(f"    call dispatch({', '.join(keywords)})")
    result.append(f"end subroutine {name}")
    return result


def transform_unit(lines: list[str], path: Path, unit: EntryUnit) -> list[str]:
    declarations = argument_declarations(lines, unit)
    state_module = module_name(path, unit)
    all_arguments = tuple(dict.fromkeys(unit.initializer_args + unit.entry_args))
    original = lines[unit.start : unit.end + 1]
    start_relative = 0
    executable_relative = unit.executable_index - unit.start
    return_relative = unit.return_index - unit.start
    entry_relative = unit.entry_index - unit.start
    end_relative = unit.end - unit.start
    original[start_relative] = signature("dispatch", all_arguments + ("initialize_only",))
    implicit_relative = next(
        index for index, line in enumerate(original) if re.match(r"\s*implicit\s+none\s*$", line, re.IGNORECASE)
    )
    insertions = ["    logical, intent(in) :: initialize_only"]
    if all_arguments:
        insertions.append(f"    optional :: {', '.join(all_arguments)}")
    original[implicit_relative + 1 : implicit_relative + 1] = insertions
    offset = len(insertions)
    executable_relative += offset
    return_relative += offset
    entry_relative += offset
    end_relative += offset
    for argument in all_arguments:
        declaration = declarations[argument]
        shape_match = re.search(r"(\([^)]*\))\s*$", declaration)
        if not shape_match:
            continue
        rank = shape_match.group(1).count(",") + 1
        assumed_shape = "(" + ",".join(":" for _ in range(rank)) + ")"
        pattern = re.compile(rf"\b{re.escape(argument)}\s*\([^)]*\)", re.IGNORECASE)
        for line_index in range(implicit_relative + 1, executable_relative):
            original[line_index] = pattern.sub(f"{argument}{assumed_shape}", original[line_index])
    original.insert(executable_relative, "    if( initialize_only )then")
    return_relative += 1
    entry_relative += 1
    end_relative += 1
    if unit.preserve_return:
        pass
    elif unit.executable_index == unit.return_index:
        original[return_relative] = "        continue"
    else:
        del original[return_relative]
        entry_relative -= 1
        end_relative -= 1
    original[entry_relative] = "    else"
    original.insert(end_relative, "    endif")
    end_relative += 1
    original[end_relative] = "end subroutine dispatch"
    result = [f"module {state_module}", "    implicit none", "    private", "    public :: dispatch", "contains"]
    result.extend(original)
    result.append(f"end module {state_module}")
    result.append("")
    result.extend(wrapper(unit.initializer_name, unit.initializer_args, declarations, state_module, True))
    result.append("")
    result.extend(wrapper(unit.entry_name, unit.entry_args, declarations, state_module, False))
    return result


def migrate(path: Path, write: bool) -> bool:
    original_text = path.read_text()
    lines = normalize_identifier_continuations(original_text.splitlines())
    units = find_units(lines)
    if not units:
        return False
    for unit in reversed(units):
        lines[unit.start : unit.end + 1] = transform_unit(lines, path, unit)
    new_text = "\n".join(lines) + "\n"
    if new_text == original_text:
        return False
    if write:
        path.write_text(new_text)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    changed = 0
    for path in args.paths:
        if migrate(path, args.write):
            changed += 1
            print(path)
    print(f"fortran_sources_changed: {changed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
