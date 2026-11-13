#!/usr/bin/env python3
"""Compile-check every CHEMPOTPY Fortran source in isolated directories."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import SurfaceSource, discover_surface_sources
from scripts.build_pes import SUPPORT_SOURCES, core_fortran_dependencies


@dataclass(frozen=True)
class CheckResult:
    """Store the compiler result for one isolated Fortran source."""

    system_name: str
    surface_name: str
    path: str
    entry_point: Optional[str]
    passed: bool
    command: List[str]
    output: str


def compiler_command(
    source: SurfaceSource, module_directory: Path, compiler: str
) -> List[str]:
    """Return the strict compatibility command for one source form."""

    command = [
        compiler,
        "-fsyntax-only",
        "-fallow-argument-mismatch",
        "-std=legacy",
        "-J",
        str(module_directory),
        "-I",
        str(source.path.parent),
    ]
    if source.source_form == "fixed":
        command.extend(("-ffixed-form", "-ffixed-line-length-72"))
    else:
        command.extend((
            "-ffree-form",
            "-ffree-line-length-none",
            "-fimplicit-none",
        ))
    source_text = (
        source.path.read_text(errors="replace")
        if source.path.exists()
        else ""
    )
    support_sources = tuple(
        source.path.parent / support_name
        for support_name in SUPPORT_SOURCES.get(
            (source.system_name, source.surface_name), ()
        )
    )
    command.extend(
        str(path) for path in core_fortran_dependencies(
            source_text, support_sources
        )
    )
    command.extend(str(path) for path in support_sources)
    command.append(str(source.path))
    return command


def check_source(
    source: SurfaceSource, work_root: Path, compiler: str
) -> CheckResult:
    """Run the compiler syntax gate for one source."""

    module_directory = work_root / source.system_name / source.surface_name
    module_directory.mkdir(parents=True, exist_ok=True)
    command = compiler_command(source, module_directory, compiler)
    process = subprocess.run(
        command,
        cwd=module_directory,
        capture_output=True,
        text=True,
    )
    output = "\n".join(
        part.rstrip() for part in (process.stdout, process.stderr) if part
    )
    return CheckResult(
        system_name=source.system_name,
        surface_name=source.surface_name,
        path=str(source.path.relative_to(REPOSITORY_ROOT)),
        entry_point=source.entry_point,
        passed=process.returncode == 0,
        command=command,
        output=output,
    )


def select_sources(
    system_name: Optional[str],
    surface_name: Optional[str],
    include_support: bool,
) -> Sequence[SurfaceSource]:
    """Select deterministic sources for a requested validation scope."""

    sources = discover_surface_sources()
    return tuple(
        source
        for source in sources
        if (include_support or source.is_surface)
        and (system_name is None or source.system_name == system_name)
        and (surface_name is None or source.surface_name == surface_name)
    )


def main() -> int:
    """Check selected sources and return failure when any source fails."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--system")
    parser.add_argument("--surface")
    parser.add_argument("--include-support", action="store_true")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")

    sources = select_sources(
        args.system, args.surface, args.include_support
    )
    if not sources:
        parser.error("no Fortran sources matched the requested scope")

    with tempfile.TemporaryDirectory(prefix="chempotpy-fortran-check-") as tmp:
        work_root = Path(tmp)
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.jobs
        ) as executor:
            futures = [
                executor.submit(check_source, source, work_root, args.compiler)
                for source in sources
            ]
            results = [future.result() for future in futures]

    results.sort(key=lambda item: (item.system_name, item.surface_name))
    failures = [result for result in results if not result.passed]
    print("fortran_sources_checked: {}".format(len(results)))
    print("fortran_sources_passed: {}".format(len(results) - len(failures)))
    print("fortran_sources_failed: {}".format(len(failures)))
    for failure in failures:
        print("FAIL {}/{}".format(
            failure.system_name, failure.surface_name
        ))
        if failure.output:
            print(failure.output)

    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(
                [asdict(result) for result in results],
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
