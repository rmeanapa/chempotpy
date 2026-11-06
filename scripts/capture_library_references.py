#!/usr/bin/env python3
"""Capture deterministic outputs from the remaining CHEMPOTPY interfaces."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.build_pes import (
    INTERNAL_BLAS,
    INTERNAL_DSYEV,
    SUPPORT_SOURCES,
    core_fortran_dependencies,
)
ENTRY_PATTERN = re.compile(
    r"^\s*(?:recursive\s+)?subroutine\s+(pes|dpem)\s*\((.*?)\)",
    flags=re.IGNORECASE | re.MULTILINE | re.DOTALL,
)
INTEGER_PATTERN = r"\b{}\s*=\s*(\d+)"


@dataclass(frozen=True)
class CaptureResult:
    """Hold one source's reference record or diagnostic."""

    key: str
    record: Optional[dict]
    error: str


def interface_signature(source_text: str) -> tuple[str, str]:
    """Return the entry point and normalized historical argument signature."""

    match = ENTRY_PATTERN.search(source_text)
    if match is None:
        raise ValueError("no pes or dpem entry point")
    arguments = re.sub(r"[&\s]+", "", match.group(2)).lower()
    return match.group(1).lower(), arguments


def integer_parameter(source_text: str, name: str) -> int:
    """Read one integer parameter from a historical wrapper."""

    match = re.search(INTEGER_PATTERN.format(name), source_text, re.IGNORECASE)
    if match is None:
        raise ValueError("missing {} parameter".format(name))
    return int(match.group(1))


def first_integer_parameter(source_text: str, *names: str) -> int:
    """Return the first named integer parameter present in a wrapper."""

    for name in names:
        try:
            return integer_parameter(source_text, name)
        except ValueError:
            continue
    raise ValueError("missing {} parameter".format(" or ".join(names)))


def coordinates(
    natoms: int,
    geometry_index: int,
    geometry_variant: str = "",
) -> tuple[tuple[float, ...], ...]:
    """Return a deterministic, non-symmetric test geometry in angstroms."""

    points = []
    for index in range(natoms):
        if geometry_index == 1:
            point = (
                0.93*(index % 4) + 0.07*(index//4),
                1.07*((index//4) % 4) + 0.11*(index % 3),
                0.89*(index//16) + 0.13*((2*index + 1) % 5),
            )
        else:
            point = (
                1.01*(index % 3) + 0.17*((index//3) % 2),
                0.97*((index//3) % 3) + 0.09*((3*index + 2) % 4),
                1.13*(index//9) + 0.15*((index + 3) % 6),
            )
            if geometry_variant == "torsion":
                point = (point[0], point[1], point[2] + 0.017*index*index)
        points.append(point)
    return tuple(points)


def assignments(
    natoms: int,
    geometry_index: int,
    array_name: str,
    geometry_variant: str = "",
) -> str:
    """Return coordinate assignments for a Fortran driver."""

    def literal(value: float) -> str:
        return "{:.17e}".format(value).replace("e", "d")

    return "\n".join(
        "    {}({},:) = [{}, {}, {}]".format(
            array_name,
            index,
            literal(point[0]),
            literal(point[1]),
            literal(point[2]),
        ).replace("e+", "e+").replace("e-", "e-").replace("ed", "d")
        for index, point in enumerate(
            coordinates(natoms, geometry_index, geometry_variant), 1
        )
    )


def standard_driver(
    entry_point: str,
    signature: str,
    natoms: int,
    nstates: int,
    data_path: Path,
    system_name: str,
) -> str:
    """Create a driver for one of the historical public signatures."""

    arguments = signature.split(",")
    if entry_point == "pes" and len(arguments) == 5:
        declarations = """    real(8) :: x({natoms},3), p({nstates})
    real(8) :: g({nstates},{natoms},3)
    real(8) :: d({nstates},{nstates},{natoms},3)""".format(
            natoms=natoms, nstates=nstates
        )
        call = "call pes( x, 0, p, g, d )"
        writes = "\n".join(
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'CHEMPOTPY_REF', {0}, p({0})".format(index)
            for index in range(1, nstates + 1)
        )
        coordinate_name = "x"
    elif (
        entry_point == "pes"
        and len(arguments) == 6
        and arguments[2] == "path"
    ):
        declarations = """    real(8) :: x({natoms},3), p({nstates})
    real(8) :: g({nstates},{natoms},3)
    real(8) :: d({nstates},{nstates},{natoms},3)
    character(len=1024) :: data_path""".format(
            natoms=natoms, nstates=nstates
        )
        call = "call pes( x, 0, data_path, p, g, d )"
        writes = "\n".join(
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'CHEMPOTPY_REF', {0}, p({0})".format(index)
            for index in range(1, nstates + 1)
        )
        coordinate_name = "x"
    elif entry_point == "pes" and len(arguments) == 6:
        natoms = 5
        declarations = """    real(8) :: x(10000,3), p({nstates})
    real(8) :: g({nstates},{natoms},3)
    real(8) :: d({nstates},{nstates},{natoms},3)""".format(
            natoms=natoms, nstates=nstates
        )
        call = "call pes( x, 5, 0, p, g, d )"
        writes = "\n".join(
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'CHEMPOTPY_REF', {0}, p({0})".format(index)
            for index in range(1, nstates + 1)
        )
        coordinate_name = "x"
    elif entry_point == "dpem" and len(arguments) == 4:
        declarations = """    real(8) :: x({natoms},3)
    real(8) :: u({nstates},{nstates})
    real(8) :: ug({nstates},{nstates},{natoms},3)""".format(
            natoms=natoms, nstates=nstates
        )
        call = "call dpem( x, 0, u, ug )"
        writes = "\n".join(
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'CHEMPOTPY_REF', {0}, u({1},{2})".format(
                (row - 1)*nstates + column, row, column
            )
            for row in range(1, nstates + 1)
            for column in range(1, nstates + 1)
        )
        coordinate_name = "x"
    elif entry_point == "dpem" and len(arguments) == 8:
        declarations = """    real(8) :: x({natoms},3), v({nstates})
    real(8) :: tg({nstates},{natoms},3)
    real(8) :: td({nstates},{nstates},{natoms},3)
    real(8) :: t({nstates},{nstates}), u({nstates},{nstates})
    real(8) :: ug({nstates},{nstates},{natoms},3)""".format(
            natoms=natoms, nstates=nstates
        )
        call = "call dpem( x, 0, v, tg, td, t, u, ug )"
        writes = "\n".join(
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'CHEMPOTPY_REF', {0}, u({1},{2})".format(
                (row - 1)*nstates + column, row, column
            )
            for row in range(1, nstates + 1)
            for column in range(1, nstates + 1)
        )
        coordinate_name = "x"
    elif entry_point == "pes" and len(arguments) == 8:
        natoms = 3
        declarations = """    real(8) :: xyz(10000,3), p({nstates})
    real(8) :: g({nstates},{natoms},3)
    real(8) :: d({nstates},{nstates},{natoms},3)""".format(
            natoms=natoms, nstates=nstates
        )
        call = "call pes( xyz, 1, 2, 3, 0, p, g, d )"
        writes = "\n".join(
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'CHEMPOTPY_REF', {0}, p({0})".format(index)
            for index in range(1, nstates + 1)
        )
        coordinate_name = "xyz"
    elif entry_point == "pes" and len(arguments) == 9:
        natoms = 3
        declarations = """    real(8) :: xyz(10000,3), p({nstates})
    real(8) :: g({nstates},{natoms},3)
    real(8) :: d({nstates},{nstates},{natoms},3)""".format(
            natoms=natoms, nstates=nstates
        )
        call = "call pes( xyz, 1, 1, 1, 3, 0, p, g, d )"
        writes = "\n".join(
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'CHEMPOTPY_REF', {0}, p({0})".format(index)
            for index in range(1, nstates + 1)
        )
        coordinate_name = "xyz"
    else:
        raise ValueError("unsupported signature {}({}".format(
            entry_point, signature
        ))

    blocks = []
    for geometry_index in (1, 2):
        geometry_variant = (
            "torsion"
            if system_name == "C6H6O" and geometry_index == 2
            else ""
        )
        blocks.append(
            "    {name} = 0.0d0\n{coordinates}\n    {call}\n{writes}".format(
                name=coordinate_name,
                coordinates=assignments(
                    natoms,
                    geometry_index,
                    coordinate_name,
                    geometry_variant,
                ),
                call=call,
                writes=writes,
            )
        )
    path_assignment = ""
    if "path" in arguments:
        path_assignment = "    data_path = {!r}\n".format(
            str(data_path.resolve()) + os.sep
        )
    return """
program capture_library_reference
    implicit none
{declarations}
{path_assignment}{blocks}
end program capture_library_reference
""".format(
        declarations=declarations,
        path_assignment=path_assignment,
        blocks="\n".join(blocks),
    )


def compile_dependencies(source: Path, source_text: str) -> tuple[Path, ...]:
    """Return extra Fortran sources needed to link one historical surface."""

    dependencies = list(core_fortran_dependencies(source_text))
    if re.search(r"\bcall\s+dsyev\b", source_text, re.IGNORECASE):
        dependencies.append(INTERNAL_DSYEV)
    if re.search(r"\b(?:ddot|dgemv|dgemm)\b", source_text, re.IGNORECASE):
        dependencies.append(INTERNAL_BLAS)
    for support_name in SUPPORT_SOURCES.get(
        (source.parent.name, source.stem), ()
    ):
        dependencies.append(source.parent / support_name)
    return tuple(dependencies)


def capture_source(compiler: str, source: Path, work_root: Path) -> CaptureResult:
    """Compile and run one source, returning a serializable reference record."""

    relative = source.relative_to(REPOSITORY_ROOT)
    key = "{}/{}".format(source.parent.name, source.stem)
    try:
        source_text = source.read_text(errors="replace")
        entry_point, signature = interface_signature(source_text)
        arguments = signature.split(",")
        if len(arguments) in {8, 9}:
            natoms = 3
            nstates = first_integer_parameter(
                source_text, "nstates", "state_count"
            )
        elif entry_point == "pes" and len(arguments) == 6 and arguments[2] != "path":
            natoms = 5
            nstates = first_integer_parameter(
                source_text, "nstates", "state_count"
            )
        else:
            natoms = first_integer_parameter(
                source_text, "natoms", "atom_count"
            )
            nstates = first_integer_parameter(
                source_text, "nstates", "state_count"
            )
        work_directory = work_root / source.parent.name / source.stem
        work_directory.mkdir(parents=True)
        driver = work_directory / "driver.f90"
        executable = work_directory / "capture"
        driver.write_text(standard_driver(
            entry_point,
            signature,
            natoms,
            nstates,
            REPOSITORY_ROOT / "chempotpy",
            source.parent.name,
        ))
        form_flags = (
            ["-ffixed-form", "-ffixed-line-length-72"]
            if source.suffix.lower() == ".f"
            else ["-ffree-form", "-ffree-line-length-none"]
        )
        common_flags = [
            compiler,
            "-O0",
            "-fallow-argument-mismatch",
            "-std=legacy",
            "-J",
            str(work_directory),
            "-I",
            str(source.parent),
        ]
        objects = []
        compile_inputs = [
            (
                dependency,
                ["-ffixed-form", "-ffixed-line-length-72"]
                if dependency.suffix.lower() == ".f"
                else ["-ffree-form", "-ffree-line-length-none"],
            )
            for dependency in compile_dependencies(source, source_text)
        ]
        compile_inputs.extend(((source, form_flags), (driver, [
            "-ffree-form", "-ffree-line-length-none"
        ])))
        for index, (compile_source, compile_flags) in enumerate(compile_inputs):
            object_path = work_directory / "source_{}.o".format(index)
            compiled = subprocess.run(
                [
                    *common_flags,
                    *compile_flags,
                    "-c",
                    str(compile_source),
                    "-o",
                    str(object_path),
                ],
                cwd=work_directory,
                capture_output=True,
                text=True,
            )
            if compiled.returncode:
                raise RuntimeError(
                    "compile failed for {}:\n{}".format(
                        compile_source, compiled.stderr
                    )
                )
            objects.append(object_path)
        linked = subprocess.run(
            [compiler, *(str(path) for path in objects), "-o", str(executable)],
            cwd=work_directory,
            capture_output=True,
            text=True,
        )
        if linked.returncode:
            raise RuntimeError("link failed:\n{}".format(linked.stderr))
        completed = subprocess.run(
            [str(executable)],
            cwd=source.parent,
            capture_output=True,
            text=True,
            timeout=30,
        )
        if completed.returncode:
            raise RuntimeError(
                "execution failed ({}):\n{}{}".format(
                    completed.returncode, completed.stdout, completed.stderr
                )
            )
        values = [
            line.split()[-1]
            for line in completed.stdout.splitlines()
            if line.startswith("CHEMPOTPY_REF ")
        ]
        expected_count = 2*(nstates*nstates if entry_point == "dpem" else nstates)
        if len(values) != expected_count:
            raise RuntimeError(
                "captured {} of {} values; output:\n{}".format(
                    len(values), expected_count, completed.stdout
                )
            )
        return CaptureResult(
            key=key,
            record={
                "entry_point": entry_point,
                "source_sha256": hashlib.sha256(
                    source.read_bytes()
                ).hexdigest(),
                "natoms": natoms,
                "nstates": nstates,
                "signature": signature,
                "values": values,
            },
            error="",
        )
    except Exception as error:  # Report every surface in one batch.
        return CaptureResult(
            key=key,
            record=None,
            error="{}: {}".format(relative, error),
        )


def main() -> int:
    """Capture all requested sources and write one deterministic JSON fixture."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="*", type=Path)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.all and args.sources:
        parser.error("sources and --all are mutually exclusive")
    if not args.all and not args.sources:
        parser.error("provide sources or --all")
    if args.all:
        from chempotpy.core.registry import discover_surface_sources

        args.sources = [
            item.path for item in discover_surface_sources() if item.is_surface
        ]
    compiler = shutil.which("gfortran")
    if compiler is None:
        parser.error("gfortran is required")
    sources = tuple(sorted(path.resolve() for path in args.sources))
    with tempfile.TemporaryDirectory(prefix="chempotpy-library-reference-") as tmp:
        work_root = Path(tmp)
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.jobs
        ) as executor:
            results = list(executor.map(
                lambda source: capture_source(compiler, source, work_root),
                sources,
            ))
    failures = [result for result in results if result.error]
    for failure in failures:
        print("FAIL {}: {}".format(failure.key, failure.error))
    document = {
        "surfaces": {
            result.key: result.record
            for result in sorted(results, key=lambda item: item.key)
            if result.record is not None
        }
    }
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(serialized, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print("captured: {}".format(len(results) - len(failures)))
    print("failed: {}".format(len(failures)))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
