#!/usr/bin/env python3
"""Generate immutable NH3 Surfgen parameters from the released contract."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests" / "data" / "external_parameters" / "NH3"
OUTPUT = ROOT / "chempotpy" / "NH3" / "chempotpy_nh3_surfgen_parameters.f90"


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one provenance file."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def integer_array(name: str, values: list[int], chunk: int = 8) -> str:
    """Render one integer parameter array."""

    lines = [f"    integer, parameter :: {name}({len(values)}) = [ &"]
    for start in range(0, len(values), chunk):
        group = ", ".join(str(value) for value in values[start:start+chunk])
        suffix = ", &" if start+chunk < len(values) else " ]"
        lines.append("        & "+group+suffix)
    return "\n".join(lines)


def real_array(name: str, values: list[str], chunk: int = 4) -> str:
    """Render one real parameter array with binary64-kind literals."""

    lines = [f"    real(wp), parameter :: {name}({len(values)}) = [ &"]
    for start in range(0, len(values), chunk):
        group = ", ".join(value.replace("D", "e").replace("d", "e")+"_wp" for value in values[start:start+chunk])
        suffix = ", &" if start+chunk < len(values) else " ]"
        lines.append("        & "+group+suffix)
    return "\n".join(lines)


def generated_document() -> str:
    """Parse the flattened released Surfgen contract and render its module."""

    contract = SOURCE / "nh3_surfgen_contract.txt"
    lines = contract.read_text().splitlines()
    header = lines[0].split()
    if header != ["HEADER", "25", "4", "2"]:
        raise ValueError("unexpected NH3 Surfgen contract header")
    coordinate_types: list[int] = []
    scaling_types: list[int] = []
    atom_indices: list[int] = []
    coordinate_parameters: list[str] = []
    recurrence: dict[int, tuple[list[int], list[int], list[int]]] = {}
    coefficient_order: list[int] = []
    coefficient_row: list[int] = []
    coefficient_column: list[int] = []
    coefficient_term: list[int] = []
    coefficient_value: list[str] = []
    index = 1
    while index < len(lines):
        fields = lines[index].split()
        label = fields[0]
        if label == "COORD":
            atom_count = int(fields[4])
            coordinate_types.append(int(fields[2]))
            scaling_types.append(int(fields[3]))
            atoms = [int(value) for value in fields[5:5+atom_count]]
            atom_indices.extend(atoms+[0]*(4-atom_count))
            coordinate_parameters.extend(fields[5+atom_count:7+atom_count])
            index += 1
        elif label == "RECURRENCE":
            order = int(fields[1])
            count = int(fields[2])
            parents: list[int] = []
            directions: list[int] = []
            divisors: list[int] = []
            for expected in range(1, count+1):
                index += 1
                term, parent, direction, divisor = map(int, lines[index].split())
                if term != expected:
                    raise ValueError("nonsequential Surfgen recurrence")
                parents.append(parent)
                directions.append(direction)
                divisors.append(divisor)
            recurrence[order] = (parents, directions, divisors)
            index += 1
        elif label == "COEF":
            coefficient_order.append(int(fields[1]))
            coefficient_row.append(int(fields[2]))
            coefficient_column.append(int(fields[3]))
            coefficient_term.append(int(fields[4]))
            coefficient_value.append(fields[5])
            index += 1
        else:
            raise ValueError(f"unknown contract record: {label}")
    output = [
        "!> @brief Store the immutable contracted Surfgen model for the two-state NH3 surface.",
        "!",
        "! Generated from the released Zhu-Yarkony files and a flattened Surfgen 2.8.5 contract.",
    ]
    for name in ("coord.in", "hd.data", "irrep.in", "surfgen.in", "nh3_surfgen_contract.txt"):
        output.append(f"! {name} SHA-256: {sha256(SOURCE/name)}")
    output.extend([
        "module chempotpy_nh3_surfgen_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: atom_indices, coordinate_parameters, coordinate_types, scaling_types",
        "    public :: second_order_direction, second_order_divisor, second_order_parent",
        "    public :: third_order_direction, third_order_divisor, third_order_parent",
        "    public :: fourth_order_direction, fourth_order_divisor, fourth_order_parent",
        "    public :: coefficient_column, coefficient_order, coefficient_row, coefficient_term, coefficient_value",
        integer_array("coordinate_types", coordinate_types),
        integer_array("scaling_types", scaling_types),
        "    integer, parameter :: atom_indices(4,25) = reshape([ &",
    ])
    for start in range(0, len(atom_indices), 8):
        group = ", ".join(str(value) for value in atom_indices[start:start+8])
        suffix = ", &" if start+8 < len(atom_indices) else " ], [4,25])"
        output.append("        & "+group+suffix)
    output.append("    real(wp), parameter :: coordinate_parameters(2,25) = reshape([ &")
    for start in range(0, len(coordinate_parameters), 4):
        group = ", ".join(value+"_wp" for value in coordinate_parameters[start:start+4])
        suffix = ", &" if start+4 < len(coordinate_parameters) else " ], [2,25])"
        output.append("        & "+group+suffix)
    for order, prefix in ((2, "second"), (3, "third"), (4, "fourth")):
        parents, directions, divisors = recurrence[order]
        output.extend([
            integer_array(prefix+"_order_parent", parents),
            integer_array(prefix+"_order_direction", directions),
            integer_array(prefix+"_order_divisor", divisors),
        ])
    output.extend([
        integer_array("coefficient_order", coefficient_order),
        integer_array("coefficient_row", coefficient_row),
        integer_array("coefficient_column", coefficient_column),
        integer_array("coefficient_term", coefficient_term),
        real_array("coefficient_value", coefficient_value),
        "end module chempotpy_nh3_surfgen_parameters",
        "",
    ])
    return "\n".join(output)


def main() -> int:
    """Write or verify the immutable Surfgen parameter module."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    content = generated_document()
    if arguments.check:
        if not OUTPUT.is_file() or OUTPUT.read_text() != content:
            print(f"outdated: {OUTPUT.relative_to(ROOT)}")
            return 1
    else:
        OUTPUT.write_text(content)
        print(OUTPUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
