#!/usr/bin/env python3
"""Remove only generated CHEMPOTPY build products from known locations."""

from __future__ import annotations

import shutil
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
GENERATED_SUFFIXES = (".log", ".mod", ".o", ".pyd", ".so", ".swp")


def main() -> int:
    """Remove generated files below system directories and the build tree."""

    package_root = REPOSITORY_ROOT / "chempotpy"
    removed = 0
    for system_directory in package_root.iterdir():
        if not system_directory.is_dir() or system_directory.name == "core":
            continue
        for path in system_directory.iterdir():
            if path.is_file() and path.name.endswith(GENERATED_SUFFIXES):
                path.unlink()
                removed += 1
    build_root = REPOSITORY_ROOT / "build"
    if build_root.exists():
        shutil.rmtree(build_root)
    print("generated_files_removed: {}".format(removed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
