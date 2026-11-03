#!/usr/bin/env python3
"""Remove unused generated rank metadata and its one-off COMMON block."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


COMMON = re.compile(
    r"^\s*common\s*/\s*com_para\s*/\s*g_myrank\s*,\s*g_nprocs\s*$",
    re.IGNORECASE,
)
DECLARATION = re.compile(
    r"^\s*integer(?:\s*::)?\s*g_myrank\s*,\s*g_nprocs\s*$",
    re.IGNORECASE,
)


def remove_dead_state(source: str) -> tuple[str, int]:
    """Remove an unused generated rank block when it is the only reference."""

    if len(re.findall(r"\bg_myrank\b|\bg_nprocs\b", source, re.IGNORECASE)) != 4:
        return source, 0
    lines = source.splitlines()
    common_lines = [index for index, line in enumerate(lines) if COMMON.match(line)]
    declaration_lines = [
        index for index, line in enumerate(lines) if DECLARATION.match(line)
    ]
    if len(common_lines) != 1 or len(declaration_lines) != 1:
        return source, 0
    removed = set(common_lines + declaration_lines)
    return "\n".join(
        line for index, line in enumerate(lines) if index not in removed
    ) + "\n", 1


def main() -> int:
    """Remove matching blocks from requested sources."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    removed = 0
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        modernized, count = remove_dead_state(original)
        changed += modernized != original
        removed += count
        if arguments.write and modernized != original:
            path.write_text(modernized)
    print("fortran_sources_changed: {}".format(changed))
    print("dead_common_blocks_removed: {}".format(removed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
