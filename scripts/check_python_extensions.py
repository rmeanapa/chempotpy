#!/usr/bin/env python3
"""Import and numerically validate every built CHEMPOTPY extension."""

from __future__ import annotations

import argparse
import concurrent.futures
import importlib.util
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import discover_surface_sources  # noqa: E402
from scripts.capture_library_references import coordinates  # noqa: E402


RESULT_PREFIX = "CHEMPOTPY_PYTHON_EXTENSION_RESULT "


def extension_artifact(
    output_root: Path, system_name: str, surface_name: str
) -> Path:
    """Return the unique compiled extension for one surface."""

    directory = output_root / system_name
    candidates = sorted(directory.glob(surface_name + ".*.so"))
    candidates.extend(sorted(directory.glob(surface_name + ".*.pyd")))
    if len(candidates) != 1:
        raise RuntimeError(
            "expected one extension for {}/{}, found {}".format(
                system_name, surface_name, len(candidates)
            )
        )
    return candidates[0]


def load_extension(path: Path, surface_name: str):
    """Load one extension directly from the requested output tree."""

    specification = importlib.util.spec_from_file_location(
        surface_name, path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("could not create an import specification")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def call_surface(
    module,
    record: dict,
    system_name: str,
    geometry_index: int,
    data_root: Path,
) -> Sequence[float]:
    """Evaluate one fixture geometry through its historical Python wrapper."""

    entry_point = record["entry_point"]
    arguments = record["signature"].split(",")
    atom_count = int(record["natoms"])
    state_count = int(record["nstates"])
    geometry_variant = (
        "torsion"
        if system_name == "C6H6O" and geometry_index == 2
        else ""
    )
    if entry_point == "pes" and len(arguments) == 6 and arguments[2] != "path":
        atom_count = 5
    elif entry_point == "pes" and len(arguments) in (8, 9):
        atom_count = 3
    coordinate_array = np.asarray(
        coordinates(
            atom_count,
            geometry_index,
            geometry_variant,
            system_name,
        ),
        dtype=np.float64,
    )
    if entry_point == "pes" and len(arguments) == 5:
        result = module.pes(coordinate_array, 0)
        selected = result[0]
    elif (
        entry_point == "pes"
        and len(arguments) == 6
        and arguments[2] == "path"
    ):
        result = module.pes(
            coordinate_array, 0, str(data_root.resolve()) + os.sep
        )
        selected = result[0]
    elif entry_point == "pes" and len(arguments) == 6:
        expanded = np.zeros((10000, 3), dtype=np.float64)
        expanded[:atom_count, :] = coordinate_array
        result = module.pes(expanded, atom_count, 0)
        selected = result[0]
    elif entry_point == "dpem" and len(arguments) == 4:
        result = module.dpem(coordinate_array, 0)
        selected = result[0]
    elif entry_point == "dpem" and len(arguments) == 8:
        energy = np.zeros(state_count, dtype=np.float64)
        gradient = np.zeros(
            (state_count, atom_count, 3), dtype=np.float64
        )
        coupling = np.zeros(
            (state_count, state_count, atom_count, 3),
            dtype=np.float64,
        )
        transformation = np.zeros(
            (state_count, state_count), dtype=np.float64
        )
        result = module.dpem(
            coordinate_array,
            0,
            energy,
            gradient,
            coupling,
            transformation,
        )
        selected = result[0]
    elif entry_point == "pes" and len(arguments) == 8:
        expanded = np.zeros((10000, 3), dtype=np.float64)
        expanded[:atom_count, :] = coordinate_array
        result = module.pes(expanded, 1, 2, 3, 0)
        selected = result[0]
    elif entry_point == "pes" and len(arguments) == 9:
        expanded = np.zeros((10000, 3), dtype=np.float64)
        expanded[:atom_count, :] = coordinate_array
        result = module.pes(expanded, 1, 1, 1, 3, 0)
        selected = result[0]
    else:
        raise RuntimeError(
            "unsupported signature {}({})".format(
                entry_point, record["signature"]
            )
        )
    expected_size = (
        state_count*state_count
        if entry_point == "dpem"
        else state_count
    )
    values = np.asarray(selected, dtype=np.float64).reshape(-1)
    if values.size != expected_size:
        raise RuntimeError(
            "expected {} values, received {}".format(
                expected_size, values.size
            )
        )
    return tuple(float(value) for value in values)


def validate_one(output_root: Path, key: str, record: dict) -> dict:
    """Import and evaluate one extension in the current process."""

    system_name, surface_name = key.split("/", 1)
    artifact = extension_artifact(
        output_root, system_name, surface_name
    )
    module = load_extension(artifact, surface_name)
    unexpected_modules = sorted(
        name for name in dir(module) if name.startswith("chempotpy_")
    )
    if unexpected_modules:
        raise RuntimeError(
            "internal Fortran modules are public: {}".format(
                ", ".join(unexpected_modules)
            )
        )
    actual = []
    for geometry_index in (1, 2):
        actual.extend(call_surface(
            module,
            record,
            system_name,
            geometry_index,
            output_root,
        ))
    return {"key": key, "values": actual}


def values_match(expected: str, actual: float) -> bool:
    """Return whether one Python result matches its Fortran reference."""

    expected_number = float(expected)
    if math.isnan(expected_number):
        return math.isnan(actual)
    if expected_number == actual:
        return True
    if abs(expected_number) < 1.0e-300 and abs(actual) < 1.0e-300:
        return True
    return math.isclose(
        actual,
        expected_number,
        rel_tol=1.0e-13,
        abs_tol=1.0e-13,
    )


def child_command(
    output_root: Path, fixture: Path, key: str
) -> List[str]:
    """Return the isolated validation command for one surface."""

    return [
        sys.executable,
        str(Path(__file__).resolve()),
        "--output-root",
        str(output_root),
        "--fixture",
        str(fixture),
        "--surface",
        key,
    ]


def validate_child(
    output_root: Path, fixture: Path, key: str
) -> tuple[str, Optional[dict], str]:
    """Run one surface in a subprocess so native crashes are isolated."""

    try:
        process = subprocess.run(
            child_command(output_root, fixture, key),
            cwd=REPOSITORY_ROOT,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        return key, None, "validation timed out after 300 seconds"
    if process.returncode != 0:
        return key, None, "exit {}: {}{}".format(
            process.returncode, process.stdout, process.stderr
        ).strip()
    result_lines = [
        line[len(RESULT_PREFIX):]
        for line in process.stdout.splitlines()
        if line.startswith(RESULT_PREFIX)
    ]
    if len(result_lines) != 1:
        return key, None, "missing isolated result: {}{}".format(
            process.stdout, process.stderr
        ).strip()
    return key, json.loads(result_lines[0]), ""


def main() -> int:
    """Validate one requested surface or orchestrate the full catalog."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--fixture",
        type=Path,
        default=REPOSITORY_ROOT / "tests" / "data"
        / "all_surface_references.json",
    )
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--surface", help=argparse.SUPPRESS)
    args = parser.parse_args()
    references = json.loads(args.fixture.read_text())["surfaces"]
    output_root = args.output_root.resolve()
    if args.surface:
        try:
            result = validate_one(
                output_root, args.surface, references[args.surface]
            )
        except Exception as error:
            print("{}{}".format(
                RESULT_PREFIX,
                json.dumps({"key": args.surface, "error": str(error)}),
            ))
            return 1
        print(RESULT_PREFIX + json.dumps(result))
        return 0
    catalog = {
        "{}/{}".format(source.system_name, source.surface_name)
        for source in discover_surface_sources()
        if source.is_surface
    }
    failures = []
    if catalog != set(references):
        failures.append("fixture and callable catalog differ")
    keys = sorted(catalog & set(references))
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=args.jobs
    ) as executor:
        results = executor.map(
            lambda key: validate_child(
                output_root, args.fixture.resolve(), key
            ),
            keys,
        )
        for key, result, error in results:
            if error:
                failures.append("{}: {}".format(key, error))
                continue
            expected = references[key]["values"]
            actual = result["values"]
            if len(actual) != len(expected):
                failures.append(
                    "{}: expected {} values, received {}".format(
                        key, len(expected), len(actual)
                    )
                )
                continue
            for index, (expected_value, actual_value) in enumerate(
                zip(expected, actual), 1
            ):
                if not values_match(expected_value, actual_value):
                    failures.append(
                        "{} value {}: expected {}, received {:.17e}".format(
                            key, index, expected_value, actual_value
                        )
                    )
    print("python_extensions_checked: {}".format(len(keys)))
    print("python_extension_evaluations: {}".format(2*len(keys)))
    print("python_extension_failures: {}".format(len(failures)))
    for failure in failures:
        print("FAIL {}".format(failure))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
