#!/usr/bin/env python3
"""Verify built PES coverage and absence of external numerical libraries."""

from __future__ import annotations

import argparse
import platform
import re
import subprocess
import sys
from pathlib import Path
from typing import Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import discover_surface_sources  # noqa: E402


def extension_artifacts(output_root: Path, system: str, surface: str) -> Sequence[Path]:
    """Return extension files whose module name exactly matches one surface."""

    directory = output_root / system
    artifacts = sorted(directory.glob(surface + ".*.so"))
    artifacts.extend(sorted(directory.glob(surface + ".*.pyd")))
    return tuple(artifacts)


def linked_libraries(path: Path) -> str:
    """Return the platform linker dependency report for one extension."""

    if platform.system() == "Darwin":
        command = ["otool", "-L", str(path)]
    elif platform.system() == "Linux":
        command = ["ldd", str(path)]
    else:
        raise RuntimeError("binary linkage checking is unsupported on this platform")
    return subprocess.run(
        command, check=True, capture_output=True, text=True
    ).stdout.lower()


def defined_symbols(path: Path) -> str:
    """Return externally visible symbols defined by one extension."""

    command = ["nm", "-g"]
    if platform.system() == "Darwin":
        command.append("-U")
    else:
        command.append("--defined-only")
    command.append(str(path))
    return subprocess.run(
        command, check=True, capture_output=True, text=True
    ).stdout


def main() -> int:
    """Fail for missing extensions or external MKL, LAPACK, or BLAS links."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    output_root = args.output_root.resolve()
    missing = []
    external_numerical_links = []
    missing_internal_dsyev = []
    checked = 0
    dsyev_dependent = 0
    for source in discover_surface_sources():
        if not source.is_surface:
            continue
        artifacts = extension_artifacts(
            output_root, source.system_name, source.surface_name
        )
        if len(artifacts) != 1:
            missing.append(
                "{}/{} ({})".format(
                    source.system_name, source.surface_name, len(artifacts)
                )
            )
            continue
        checked += 1
        artifact = artifacts[0]
        linkage = linked_libraries(artifact)
        if any(name in linkage for name in ("mkl", "lapack", "blas")):
            external_numerical_links.append(str(artifact))
        source_text = source.path.read_text(errors="replace")
        if re.search(r"\bcall\s+dsyev\b", source_text, flags=re.IGNORECASE):
            dsyev_dependent += 1
            symbols = defined_symbols(artifact)
            if re.search(r"\bT\s+_?dsyev_$", symbols, flags=re.MULTILINE) is None:
                missing_internal_dsyev.append(str(artifact))

    print("pes_extensions_checked: {}".format(checked))
    print("dsyev_dependent_extensions: {}".format(dsyev_dependent))
    print("missing_or_ambiguous_extensions: {}".format(len(missing)))
    print("external_mkl_lapack_blas_links: {}".format(
        len(external_numerical_links)
    ))
    print("extensions_without_internal_dsyev: {}".format(
        len(missing_internal_dsyev)
    ))
    failures = missing + external_numerical_links + missing_internal_dsyev
    for failure in failures:
        print("FAIL {}".format(failure))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
