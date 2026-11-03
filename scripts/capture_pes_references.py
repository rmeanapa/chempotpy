#!/usr/bin/env python3
"""Compile PES sources and capture deterministic energy reference values."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

try:
    from scripts.build_pes import core_fortran_dependencies
except ModuleNotFoundError:  # Direct execution adds scripts/, not the repository.
    from build_pes import core_fortran_dependencies


ATOM_COUNTS = (2, 3, 5, 8)
GEOMETRIES = {
    "cluster_a": (
        (0.0, 0.0, 0.0),
        (2.8, 0.2, 0.1),
        (0.4, 3.1, 0.3),
        (0.2, 0.6, 3.3),
        (3.0, 3.2, 0.4),
        (3.2, 0.5, 3.4),
        (0.7, 3.5, 3.6),
        (3.7, 3.8, 3.9),
    ),
    "cluster_b": (
        (0.15, -0.20, 0.35),
        (2.45, 0.65, -0.10),
        (-0.55, 2.70, 0.80),
        (0.90, -0.75, 3.05),
        (2.85, 2.25, 0.55),
        (-1.15, 1.40, 3.45),
        (3.25, -0.95, 2.80),
        (1.35, 3.60, 3.15),
    ),
}


def fortran_number(value: float) -> str:
    """Return one double-precision Fortran literal."""

    return "{:.17e}".format(value).replace("e", "d")


def driver_source() -> str:
    """Return the common Fortran driver used for reference capture."""

    blocks = []
    for geometry_name, coordinates in GEOMETRIES.items():
        assignments = "\n".join(
            "    coordinates({0},:) = [{1}]".format(
                index,
                ", ".join(fortran_number(value) for value in point),
            )
            for index, point in enumerate(coordinates, start=1)
        )
        calls = "\n".join(
            "    call pes( coordinates, {0}, 0, energy, gradient, coupling )\n"
            "    write(*,'(A,1X,I0,1X,ES26.17E3)') "
            "'{1}', {0}, energy(1)".format(atom_count, geometry_name)
            for atom_count in ATOM_COUNTS
        )
        blocks.append("{}\n{}".format(assignments, calls))
    return """
program capture_pes_reference
    implicit none
    real(8) :: coordinates(10000,3)
    real(8) :: energy(1), gradient(1,8,3), coupling(1,1,8,3)

    coordinates = 0.0d0
{blocks}
end program capture_pes_reference
""".format(blocks="\n".join(blocks))


def source_dependencies(source: Path) -> Iterable[Path]:
    """Return bundled core dependencies needed by a free-form source."""

    if source.suffix.lower() == ".f":
        return ()
    return core_fortran_dependencies(source.read_text())


def capture_source(compiler: str, source: Path, directory: Path) -> dict:
    """Compile and evaluate one PES source."""

    surface_name = source.stem
    work_directory = directory / surface_name
    work_directory.mkdir()
    driver = work_directory / "driver.f90"
    executable = work_directory / "capture"
    driver.write_text(driver_source())
    subprocess.run(
        [
            compiler,
            *(str(path) for path in source_dependencies(source)),
            str(source),
            str(driver),
            "-o",
            str(executable),
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=work_directory,
    )
    completed = subprocess.run(
        [str(executable)],
        check=True,
        capture_output=True,
        text=True,
    )
    energies = {name: {} for name in GEOMETRIES}
    for line in completed.stdout.splitlines():
        geometry_name, atom_count, energy = line.split()
        energies[geometry_name][atom_count] = float(
            energy.replace("D", "E").replace("d", "e")
        )
    expected_count = len(GEOMETRIES)*len(ATOM_COUNTS)
    actual_count = sum(len(values) for values in energies.values())
    if actual_count != expected_count:
        raise RuntimeError(
            "{} returned {} of {} reference values".format(
                surface_name, actual_count, expected_count
            )
        )
    return {
        "legacy_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "energies_ev": energies,
    }


def main() -> int:
    """Capture reference energies for all requested PES sources."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    compiler = shutil.which("gfortran")
    if compiler is None:
        parser.error("gfortran is required")
    sources = tuple(sorted(path.resolve() for path in args.sources))
    with tempfile.TemporaryDirectory(prefix="chempotpy-reference-") as tmp:
        directory = Path(tmp)
        references = {
            source.stem: capture_source(compiler, source, directory)
            for source in sources
        }
    document = {
        "atom_counts": ATOM_COUNTS,
        "geometries_angstrom": GEOMETRIES,
        "surfaces": references,
    }
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(serialized, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
