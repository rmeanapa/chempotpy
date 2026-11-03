#!/usr/bin/env python3
"""Build discovered CHEMPOTPY f2py extensions with structured results."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import SurfaceSource, discover_surface_sources


SUPPORT_SOURCES = {
    ("CH4", "CH4_GEN_SP_2001"): ("ch4T8.F",),
}
CORE_FORTRAN_ROOT = REPOSITORY_ROOT / "chempotpy" / "core" / "fortran"
CORE_MODULE_SOURCES = {
    "chempotpy_geometry": CORE_FORTRAN_ROOT / "chempotpy_geometry.f90",
    "chempotpy_kinds": CORE_FORTRAN_ROOT / "chempotpy_kinds.f90",
    "chempotpy_legacy_coordinates": (
        CORE_FORTRAN_ROOT / "chempotpy_legacy_coordinates.f90"
    ),
    "chempotpy_legacy_units": (
        CORE_FORTRAN_ROOT / "chempotpy_legacy_units.f90"
    ),
    "chempotpy_neural_network_io": (
        CORE_FORTRAN_ROOT / "chempotpy_neural_network_io.f90"
    ),
    "chempotpy_status": CORE_FORTRAN_ROOT / "chempotpy_status.f90",
    "chempotpy_units": CORE_FORTRAN_ROOT / "chempotpy_units.f90",
}
INTERNAL_DSYEV = (
    CORE_FORTRAN_ROOT / "chempotpy_dsyev.f90"
)
INTERNAL_BLAS = CORE_FORTRAN_ROOT / "chempotpy_blas.f90"
USE_PATTERN = re.compile(
    r"^\s*use(?:\s*,[^:]*)?\s*(?:::\s*)?([a-z][a-z0-9_]*)",
    flags=re.IGNORECASE | re.MULTILINE,
)
F2PY_KIND_MAP = {
    "real": {
        "real32": "float",
        "real64": "double",
    }
}


def copy_ascii_source(source: Path, destination: Path) -> None:
    """Copy Fortran source while neutralizing non-ASCII comment bytes.

    NumPy's legacy f2py parser assumes ASCII unless an optional character-set
    detector is installed. Replacing each non-ASCII byte with one space keeps
    fixed-form columns stable and removes that optional build dependency.
    """

    content = source.read_bytes()
    destination.write_bytes(bytes(byte if byte < 128 else 32 for byte in content))


def write_f2py_kind_map(build_directory: Path) -> Path:
    """Write the explicit modern-Fortran-to-C kind mapping used by f2py."""

    mapping_path = build_directory / ".f2py_f2cmap"
    mapping_path.write_text(repr(F2PY_KIND_MAP) + "\n")
    return mapping_path


@dataclass(frozen=True)
class BuildResult:
    """Store one extension build outcome and its diagnostic log path."""

    system_name: str
    surface_name: str
    source_path: str
    entry_point: str
    passed: bool
    skipped: bool
    artifact_path: Optional[str]
    log_path: str
    command: List[str]
    duration_seconds: float


def select_surfaces(
    system_name: Optional[str], surface_name: Optional[str]
) -> Sequence[SurfaceSource]:
    """Return callable surfaces in deterministic build order."""

    return tuple(
        source
        for source in discover_surface_sources()
        if source.is_surface
        and (system_name is None or source.system_name == system_name)
        and (surface_name is None or source.surface_name == surface_name)
    )


def core_fortran_dependencies(source_text: str) -> Sequence[Path]:
    """Return referenced bundled modules in dependency-first build order."""

    ordered = []
    resolved = set()
    resolving = set()

    def resolve(module_name: str) -> None:
        normalized_name = module_name.lower()
        source_path = CORE_MODULE_SOURCES.get(normalized_name)
        if source_path is None or normalized_name in resolved:
            return
        if normalized_name in resolving:
            raise ValueError(
                "circular bundled Fortran module dependency: {}".format(
                    normalized_name
                )
            )
        resolving.add(normalized_name)
        dependency_text = source_path.read_text(errors="replace")
        for dependency in USE_PATTERN.findall(dependency_text):
            resolve(dependency)
        resolving.remove(normalized_name)
        resolved.add(normalized_name)
        ordered.append(source_path)

    for referenced_module in USE_PATTERN.findall(source_text):
        resolve(referenced_module)
    return tuple(ordered)


def f2py_command(
    source: SurfaceSource,
    local_sources: Sequence[Path],
    python: str,
    backend: str,
    debug: bool = False,
) -> List[str]:
    """Return the f2py command for one normalized source set."""

    command = [
        python,
        "-m",
        "numpy.f2py",
        "-c",
        "--backend",
        backend,
        "-m",
        source.surface_name,
        "--f77flags=-fallow-argument-mismatch -std=legacy",
        "--f90flags=-fallow-argument-mismatch -std=legacy",
    ]
    if debug:
        command.append("--debug")
    command.extend(str(path) for path in local_sources)
    if not any(path.suffix == ".pyf" for path in local_sources):
        command.extend(["only:", source.entry_point or "", ":"])
    return command


def normalized_sources(
    source: SurfaceSource, build_directory: Path
) -> Sequence[Path]:
    """Copy sources with suffixes matching their actual Fortran form."""

    suffix = ".f" if source.source_form == "fixed" else ".f90"
    local_source = build_directory / (source.surface_name + suffix)
    copy_ascii_source(source.path, local_source)
    signature_source = source.path.with_suffix(".pyf")
    source_text = source.path.read_text(errors="replace")
    local_sources = []
    if signature_source.exists():
        local_signature = build_directory / signature_source.name
        copy_ascii_source(signature_source, local_signature)
        local_sources.append(local_signature)
    for core_source in core_fortran_dependencies(source_text):
        local_core = build_directory / core_source.name
        copy_ascii_source(core_source, local_core)
        local_sources.append(local_core)
    local_sources.append(local_source)
    for support_name in SUPPORT_SOURCES.get(
        (source.system_name, source.surface_name), ()
    ):
        support_source = source.path.parent / support_name
        local_support = build_directory / support_name
        copy_ascii_source(support_source, local_support)
        local_sources.append(local_support)
    if re.search(r"\bcall\s+dsyev\b", source_text, flags=re.IGNORECASE):
        local_dsyev = build_directory / INTERNAL_DSYEV.name
        copy_ascii_source(INTERNAL_DSYEV, local_dsyev)
        local_sources.append(local_dsyev)
    if re.search(
        r"\b(?:ddot|dgemv|dgemm)\b", source_text, flags=re.IGNORECASE
    ):
        local_blas = build_directory / INTERNAL_BLAS.name
        copy_ascii_source(INTERNAL_BLAS, local_blas)
        local_sources.append(local_blas)
    return tuple(local_sources)


def build_surface(
    source: SurfaceSource,
    work_root: Path,
    output_root: Path,
    log_root: Path,
    python: str,
    backend: str,
    debug: bool,
) -> BuildResult:
    """Build one surface in isolation and copy out only its extension."""

    started = time.monotonic()
    build_directory = work_root / source.system_name / source.surface_name
    build_directory.mkdir(parents=True, exist_ok=True)
    write_f2py_kind_map(build_directory)
    local_sources = normalized_sources(source, build_directory)
    command = f2py_command(
        source, local_sources, python, backend, debug
    )
    environment = os.environ.copy()
    environment["PATH"] = "{}{}{}".format(
        Path(python).absolute().parent,
        os.pathsep,
        environment.get("PATH", ""),
    )
    process = subprocess.run(
        command,
        cwd=build_directory,
        capture_output=True,
        text=True,
        env=environment,
    )
    log_path = log_root / source.system_name / (source.surface_name + ".log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "$ {}\n\n{}{}".format(
            shlex.join(command), process.stdout, process.stderr
        )
    )

    artifact_path = None
    candidates = sorted(build_directory.glob(source.surface_name + ".*.so"))
    candidates.extend(
        sorted(build_directory.glob(source.surface_name + ".*.pyd"))
    )
    passed = process.returncode == 0 and len(candidates) == 1
    if passed:
        destination_directory = output_root / source.system_name
        destination_directory.mkdir(parents=True, exist_ok=True)
        destination = destination_directory / candidates[0].name
        shutil.copy2(candidates[0], destination)
        artifact_path = str(destination)

    return BuildResult(
        system_name=source.system_name,
        surface_name=source.surface_name,
        source_path=str(source.path.relative_to(REPOSITORY_ROOT)),
        entry_point=source.entry_point or "",
        passed=passed,
        skipped=False,
        artifact_path=artifact_path,
        log_path=str(log_path),
        command=command,
        duration_seconds=time.monotonic()-started,
    )


def existing_result(
    source: SurfaceSource, output_root: Path, log_root: Path
) -> Optional[BuildResult]:
    """Return a successful result when exactly one prior artifact exists."""

    output_directory = output_root / source.system_name
    candidates = sorted(output_directory.glob(source.surface_name + ".*.so"))
    candidates.extend(
        sorted(output_directory.glob(source.surface_name + ".*.pyd"))
    )
    if len(candidates) != 1:
        return None
    return BuildResult(
        system_name=source.system_name,
        surface_name=source.surface_name,
        source_path=str(source.path.relative_to(REPOSITORY_ROOT)),
        entry_point=source.entry_point or "",
        passed=True,
        skipped=True,
        artifact_path=str(candidates[0]),
        log_path=str(
            log_root / source.system_name / (source.surface_name + ".log")
        ),
        command=[],
        duration_seconds=0.0,
    )


def main() -> int:
    """Build selected surfaces and fail when any extension does not build."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--backend", choices=("meson", "distutils"), default="meson")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--system")
    parser.add_argument("--surface")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--log-root", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")

    surfaces = select_surfaces(args.system, args.surface)
    if not surfaces:
        parser.error("no callable surfaces matched the requested scope")
    output_root = (args.output_root or (REPOSITORY_ROOT / "chempotpy")).resolve()
    log_root = (args.log_root or (REPOSITORY_ROOT / "build" / "pes-logs")).resolve()
    report_path = args.report or (REPOSITORY_ROOT / "build" / "pes-build.json")

    existing_results = []
    if args.resume:
        for source in surfaces:
            result = existing_result(source, output_root, log_root)
            if result is not None:
                existing_results.append(result)
        existing_keys = {
            (result.system_name, result.surface_name)
            for result in existing_results
        }
        surfaces = tuple(
            source for source in surfaces
            if (source.system_name, source.surface_name) not in existing_keys
        )

    with tempfile.TemporaryDirectory(prefix="chempotpy-f2py-build-") as tmp:
        work_root = Path(tmp)
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.jobs
        ) as executor:
            futures = [
                executor.submit(
                    build_surface,
                    source,
                    work_root,
                    output_root,
                    log_root,
                    args.python,
                    args.backend,
                    args.debug,
                )
                for source in surfaces
            ]
            results = existing_results + [
                future.result() for future in futures
            ]

    results.sort(key=lambda item: (item.system_name, item.surface_name))
    failures = [result for result in results if not result.passed]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(
            [asdict(result) for result in results], indent=2, sort_keys=True
        )
        + "\n"
    )
    print("pes_extensions_built: {}".format(len(results) - len(failures)))
    print("pes_extensions_skipped: {}".format(
        sum(result.skipped for result in results)
    ))
    print("pes_extensions_failed: {}".format(len(failures)))
    print("build_report: {}".format(report_path))
    for failure in failures:
        print("FAIL {}/{} ({})".format(
            failure.system_name,
            failure.surface_name,
            failure.log_path,
        ))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
