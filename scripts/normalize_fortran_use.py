#!/usr/bin/env python3
"""Replace broad Fortran USE statements with compiler-resolved ONLY lists."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List

try:
    from scripts.explicit_fortran_types import compiler_dump
    from scripts.privatize_fortran_modules import use_associations
except ModuleNotFoundError:  # Direct execution from the scripts directory.
    from explicit_fortran_types import compiler_dump
    from privatize_fortran_modules import use_associations


USE_PATTERN = re.compile(
    r"^(?P<indent>\s*)use(?P<nature>\s*,\s*(?:non_)?intrinsic)?\s*"
    r"(?:::\s*)?(?P<module>[a-z][a-z0-9_]*)\s*(?:!.*)?$",
    flags=re.IGNORECASE,
)


def use_lines(module: str, names: set[str], indent: str) -> List[str]:
    """Format a USE statement with a wrapped ONLY list."""

    ordered = sorted(names)
    prefix = indent + "use {}, only: ".format(module)
    continuation = indent + "    & "
    result: List[str] = []
    current = prefix
    for name in ordered:
        addition = name if current.endswith(": ") or current.endswith("& ") else ", " + name
        if len(current) + len(addition) <= 100:
            current += addition
        else:
            result.append(current + ", &")
            current = continuation + name
    result.append(current)
    return result


def rewrite_source(source_text: str, associations: dict[str, set[str]]) -> tuple[str, int]:
    """Return source with broad USE imports made explicit."""

    output: List[str] = []
    changed = 0
    for line in source_text.splitlines():
        match = USE_PATTERN.match(line)
        if match is None:
            output.append(line)
            continue
        module = match.group("module").lower()
        names = associations.get(module, set())
        if not names:
            output.append(line)
            continue
        nature = match.group("nature") or ""
        if nature:
            prefix_module = nature.strip(" ,") + " :: " + module
        else:
            prefix_module = module
        output.extend(use_lines(prefix_module, names, match.group("indent")))
        changed += 1
    return "\n".join(output) + "\n", changed


def main() -> int:
    """Normalize broad USE imports in requested sources."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    unresolved = 0
    for source in arguments.sources:
        original = source.read_text(errors="strict")
        broad = [match for line in original.splitlines() if (match := USE_PATTERN.match(line))]
        if not broad:
            continue
        associations = use_associations(compiler_dump(source, arguments.compiler))
        rewritten, source_changed = rewrite_source(original, associations)
        changed += source_changed
        unresolved += len(broad) - source_changed
        if arguments.write and rewritten != original:
            source.write_text(rewritten)
    print("broad_use_statements_changed: {}".format(changed))
    print("broad_use_statements_unresolved: {}".format(unresolved))
    return 1 if unresolved else 0


if __name__ == "__main__":
    raise SystemExit(main())
