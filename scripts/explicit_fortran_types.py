#!/usr/bin/env python3
"""Replace implicit Fortran typing with compiler-resolved declarations.

The legacy sources rely on implicit maps such as ``A-H,O-Z`` being double
precision.  GFortran's parsed namespace dump is used as the source of truth so
that this migration preserves the compiler-resolved type of every symbol.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import DefaultDict, Iterable, List, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.build_pes import core_fortran_dependencies


UNIT_PATTERN = re.compile(
    r"^(?P<indent>\s*)(?:"
    r"module\s+(?!procedure\b)(?P<module>[a-z][a-z0-9_]*)|"
    r"block\s+data(?:\s+(?P<block>[a-z][a-z0-9_]*))?|"
    r"program\s+(?P<program>[a-z][a-z0-9_]*)|"
    r"(?:(?:recursive|pure|elemental|impure|module)\s+)*"
    r"(?:(?:double\s+precision|real(?:\s*\([^)]*\)|\*\d+)?|"
    r"integer(?:\s*\([^)]*\)|\*\d+)?|logical(?:\s*\([^)]*\)|\*\d+)?|"
    r"complex(?:\s*\([^)]*\)|\*\d+)?|"
    r"character(?:\s*\([^)]*\)|\*\d+)?)\s+)?"
    r"(?P<kind>subroutine|function)\s+(?P<procedure>[a-z][a-z0-9_]*)"
    r")",
    flags=re.IGNORECASE,
)
IMPLICIT_PATTERN = re.compile(r"^\s*implicit\b", flags=re.IGNORECASE)
IMPLICIT_NONE_PATTERN = re.compile(
    r"^\s*implicit\s+none\b", flags=re.IGNORECASE
)
SYMBOL_PATTERN = re.compile(
    r"^\s*symtree:\s*'([^']+)'[^\n]*\n"
    r"(?P<body>.*?)(?=^\s*symtree:|^\s*code:|^Namespace:|\Z)",
    flags=re.MULTILINE | re.DOTALL,
)
TYPE_PATTERN = re.compile(r"^\s*type spec\s*:\s*\(([^)]*)\)", re.MULTILINE)
ATTRIBUTES_PATTERN = re.compile(
    r"^\s*attributes:\s*\(([^)]*)\)", re.MULTILINE
)
PROCEDURE_NAME_PATTERN = re.compile(
    r"^procedure name\s*=\s*(\S+)", re.MULTILINE
)


@dataclass
class Namespace:
    """Compiler-resolved implicit symbols for one source unit."""

    name: str
    declarations: DefaultDict[str, List[str]] = field(
        default_factory=lambda: defaultdict(list)
    )


@dataclass
class SourceUnit:
    """Location of one source unit declaration."""

    name: str
    line_index: int
    next_line_index: int
    indent: str


def compiler_dump(
    source: Path, compiler: str, source_text: str | None = None
) -> str:
    """Return GFortran's resolved namespace dump for *source*."""

    with tempfile.TemporaryDirectory(prefix="chempotpy-types-") as directory:
        source = source.resolve()
        include_directory = source.parent
        if source_text is None:
            source_text = source.read_text(errors="strict")
            compile_source = source
        else:
            compile_source = Path(directory) / source.name
            compile_source.write_text(source_text)
        dependencies = core_fortran_dependencies(
            source_text
        )
        for dependency in dependencies:
            dependency_command = [
                compiler,
                "-c",
                "-std=f2008",
                "-J",
                directory,
                "-I",
                directory,
                str(dependency.resolve()),
            ]
            dependency_result = subprocess.run(
                dependency_command,
                cwd=directory,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            if dependency_result.returncode != 0:
                raise RuntimeError(
                    "dependency compile failed for {}:\n{}".format(
                        dependency, dependency_result.stdout
                    )
                )
        command = [
            compiler,
            "-fsyntax-only",
            "-fallow-argument-mismatch",
            "-std=legacy",
            "-ffree-form",
            "-ffree-line-length-none",
            "-fdump-fortran-original",
            "-J",
            directory,
            "-I",
            directory,
            "-I",
            str(include_directory),
            str(compile_source),
        ]
        result = subprocess.run(
            command,
            cwd=directory,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(
            "compiler namespace dump failed for {}:\n{}\n{}".format(
                source, result.stdout, result.stderr
            )
        )
    return result.stdout


def declaration_type(type_spec: str) -> str | None:
    """Translate a GFortran type specification to a declaration type."""

    fields = type_spec.split()
    if not fields:
        return None
    category = fields[0].upper()
    if category == "REAL":
        return "double precision" if fields[1:2] == ["8"] else "real"
    if category == "INTEGER":
        return "integer" if fields[1:2] in ([], ["4"]) else "integer*{}".format(fields[1])
    if category == "LOGICAL":
        return "logical" if fields[1:2] in ([], ["4"]) else "logical*{}".format(fields[1])
    if category == "COMPLEX":
        return "double complex" if fields[1:2] == ["16"] else "complex"
    if category == "CHARACTER":
        length = fields[1].split("_")[0] if len(fields) > 1 else "1"
        return "character(len={})".format(length)
    return None


def normalized_namespace_name(name: str) -> str:
    """Map compiler-generated entry-master names to source procedure names."""

    lowered = name.lower()
    if lowered.startswith("master."):
        return lowered.rsplit(".", 1)[-1]
    return lowered


def parse_namespaces(dump: str) -> List[Namespace]:
    """Extract implicitly typed, non-intrinsic symbols from a compiler dump."""

    starts = [match.start() for match in re.finditer(r"^Namespace:", dump, re.MULTILINE)]
    namespaces: List[Namespace] = []
    for position, start in enumerate(starts):
        end = starts[position + 1] if position + 1 < len(starts) else len(dump)
        block = dump[start:end]
        name_match = PROCEDURE_NAME_PATTERN.search(block)
        if name_match is None:
            continue
        namespace = Namespace(normalized_namespace_name(name_match.group(1)))
        for symbol_match in SYMBOL_PATTERN.finditer(block):
            symbol = symbol_match.group(1).lower()
            if re.fullmatch(r"[a-z][a-z0-9_]*", symbol) is None:
                continue
            body = symbol_match.group("body")
            attributes_match = ATTRIBUTES_PATTERN.search(body)
            type_match = TYPE_PATTERN.search(body)
            if attributes_match is None or type_match is None:
                continue
            attributes = attributes_match.group(1).upper()
            if "IMPLICIT-TYPE" not in attributes:
                continue
            if "INTRINSIC-PROC" in attributes:
                continue
            declared_type = declaration_type(type_match.group(1))
            if declared_type is None:
                continue
            namespace.declarations[declared_type].append(symbol)
        namespaces.append(namespace)
    return namespaces


def parse_source_units(lines: Sequence[str]) -> List[SourceUnit]:
    """Locate source-unit declarations and their following source unit."""

    starts: List[tuple[str, int, str]] = []
    for line_index, line in enumerate(lines):
        if re.match(r"^\s*end\b", line, re.IGNORECASE):
            continue
        match = UNIT_PATTERN.match(line)
        if match is None:
            continue
        name = next(
            (
                match.group(group)
                for group in ("module", "block", "program", "procedure")
                if match.group(group)
            ),
            "__unnamed_block_data__",
        )
        starts.append((name.lower(), line_index, match.group("indent")))
    units: List[SourceUnit] = []
    for position, (name, line_index, indent) in enumerate(starts):
        next_index = starts[position + 1][1] if position + 1 < len(starts) else len(lines)
        units.append(SourceUnit(name, line_index, next_index, indent))
    return units


def declaration_lines(
    declared_type: str, names: Iterable[str], indent: str, width: int = 100
) -> List[str]:
    """Format one explicit declaration using four-space continuations."""

    unique_names = sorted(set(names))
    prefix = indent + declared_type + " :: "
    continuation_prefix = indent + "    & "
    result: List[str] = []
    current = prefix
    for name in unique_names:
        addition = name if current.endswith(":: ") or current.endswith("& ") else ", " + name
        if len(current) + len(addition) <= width:
            current += addition
            continue
        result.append(current + ", &")
        current = continuation_prefix + name
    result.append(current)
    return result


def header_end(lines: Sequence[str], start: int) -> int:
    """Return the first line after a continued source-unit declaration."""

    index = start
    while index < len(lines) and lines[index].rstrip().endswith("&"):
        index += 1
    return min(index + 1, len(lines))


def implicit_insertion_index(lines: Sequence[str], unit: SourceUnit) -> int:
    """Find the valid location after a unit header and any USE/IMPORT lines."""

    index = header_end(lines, unit.line_index)
    while index < unit.next_line_index:
        stripped = lines[index].strip().lower()
        if not stripped or stripped.startswith("!"):
            index += 1
            continue
        if stripped.startswith("use ") or stripped == "use" or stripped.startswith("import"):
            while index < unit.next_line_index and lines[index].rstrip().endswith("&"):
                index += 1
            index += 1
            continue
        break
    return index


def specification_statements(
    lines: Sequence[str], start: int, end: int
) -> Iterable[str]:
    """Yield logical statements from a free-form specification region."""

    current = ""
    for line in lines[start:end]:
        code = line.split("!", 1)[0].strip()
        if not code:
            continue
        if code.startswith("&"):
            code = code[1:].lstrip()
        continues = code.endswith("&")
        if continues:
            code = code[:-1].rstrip()
        current = (current + " " + code).strip()
        if not continues:
            yield current
            current = ""
    if current:
        yield current


def explicitly_declared_names(
    lines: Sequence[str], unit: SourceUnit
) -> set[str]:
    """Return names already carrying an explicit type in a source unit."""

    type_pattern = re.compile(
        r"^(?:double\s+precision|double\s+complex|"
        r"real(?:\s*\([^)]*\)|\*\d+)?|"
        r"integer(?:\s*\([^)]*\)|\*\d+)?|"
        r"logical(?:\s*\([^)]*\)|\*\d+)?|"
        r"complex(?:\s*\([^)]*\)|\*\d+)?|"
        r"character(?:\s*\([^)]*\)|\*\d+)?)\b(?P<rest>.*)$",
        flags=re.IGNORECASE,
    )
    names: set[str] = set()
    for statement in specification_statements(
        lines, unit.line_index + 1, unit.next_line_index
    ):
        match = type_pattern.match(statement)
        if match is None or re.search(r"\bfunction\b", statement, re.IGNORECASE):
            continue
        rest = match.group("rest")
        if "::" in rest:
            rest = rest.split("::", 1)[1]
        elif rest.lstrip().startswith(","):
            continue
        depth = 0
        entity = ""
        entities: List[str] = []
        for character in rest:
            if character == "(":
                depth += 1
            elif character == ")":
                depth = max(0, depth - 1)
            if character == "," and depth == 0:
                entities.append(entity)
                entity = ""
            else:
                entity += character
        entities.append(entity)
        for item in entities:
            name_match = re.match(r"\s*([a-z][a-z0-9_]*)", item, re.IGNORECASE)
            if name_match is not None:
                names.add(name_match.group(1).lower())
    return names


def rewrite_source(source_text: str, namespaces: Sequence[Namespace]) -> tuple[str, int]:
    """Return source text with implicit typing replaced by explicit declarations."""

    lines = source_text.splitlines()
    units = parse_source_units(lines)
    by_name: DefaultDict[str, List[SourceUnit]] = defaultdict(list)
    for unit in units:
        by_name[unit.name].append(unit)

    edits: List[tuple[int, int, List[str]]] = []
    changed_units = 0
    used_units: set[int] = set()
    namespace_names = {namespace.name for namespace in namespaces}
    for namespace in namespaces:
        candidates = by_name.get(namespace.name, [])
        unit = next((candidate for candidate in candidates if id(candidate) not in used_units), None)
        if unit is None:
            if not namespace.declarations:
                continue
            raise RuntimeError(
                "cannot match compiler namespace {!r} to a source unit".format(
                    namespace.name
                )
            )
        used_units.add(id(unit))
        implicit_indexes = [
            index
            for index in range(unit.line_index + 1, unit.next_line_index)
            if IMPLICIT_PATTERN.match(lines[index])
        ]
        declaration_indent = unit.indent + "    "
        explicit_names = explicitly_declared_names(lines, unit)
        generated: List[str] = []
        for declared_type in sorted(namespace.declarations):
            undeclared_names = [
                name for name in namespace.declarations[declared_type]
                if name not in explicit_names
            ]
            if not undeclared_names:
                continue
            generated.extend(
                declaration_lines(
                    declared_type,
                    undeclared_names,
                    declaration_indent,
                )
            )
        if implicit_indexes:
            first = implicit_indexes[0]
            replacement = [declaration_indent + "implicit none", ""] + generated
            edits.append((first, first + 1, replacement))
            for duplicate in implicit_indexes[1:]:
                edits.append((duplicate, duplicate + 1, []))
        else:
            insertion = implicit_insertion_index(lines, unit)
            edits.append(
                (
                    insertion,
                    insertion,
                    [declaration_indent + "implicit none", ""] + generated,
                )
            )
        changed_units += 1

    unmatched_units = [
        unit.name for unit in units
        if id(unit) not in used_units and unit.name in namespace_names
    ]
    if unmatched_units:
        raise RuntimeError(
            "compiler namespaces were not uniquely matched for: {}".format(
                ", ".join(sorted(set(unmatched_units)))
            )
        )

    for start, end, replacement in sorted(edits, reverse=True):
        lines[start:end] = replacement
    return "\n".join(lines) + "\n", changed_units


def git_baseline(source: Path, revision: str) -> str:
    """Return a tracked source's content at *revision*."""

    relative = source.resolve().relative_to(REPOSITORY_ROOT)
    result = subprocess.run(
        ["git", "show", "{}:{}".format(revision, relative.as_posix())],
        cwd=REPOSITORY_ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "cannot read {} at {}: {}".format(
                relative, revision, result.stderr.strip()
            )
        )
    return result.stdout


def looks_fixed_form(source_text: str) -> bool:
    """Return whether a nominal .f90 source still uses fixed continuations."""

    return any(
        len(line) >= 6
        and line[5] == "$"
        and not line.lstrip().startswith("!")
        for line in source_text.splitlines()
    )


def modernize_source(
    source: Path,
    compiler: str,
    write: bool,
    baseline_revision: str | None = None,
) -> int:
    """Modernize implicit declarations in one free-form source."""

    original = (
        git_baseline(source, baseline_revision)
        if baseline_revision is not None
        else source.read_text(errors="strict")
    )
    if looks_fixed_form(original):
        from scripts.modernize_fixed_form import convert_source_text

        original = convert_source_text(original, source.stem)
    dump = compiler_dump(source, compiler, original)
    namespaces = parse_namespaces(dump)
    rewritten, changed_units = rewrite_source(original, namespaces)
    if write and rewritten != original:
        source.write_text(rewritten)
    return changed_units


def main() -> int:
    """Modernize the requested free-form sources."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--baseline-revision")
    arguments = parser.parse_args()

    total = 0
    for source in arguments.sources:
        changed_units = modernize_source(
            source,
            arguments.compiler,
            arguments.write,
            arguments.baseline_revision,
        )
        total += changed_units
        print("{}: {} source units".format(source, changed_units))
    print("resolved explicit declarations for {} source units".format(total))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
