#!/usr/bin/env python3
"""Generate miscellaneous immutable PES parameter modules."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.format_fortran import format_source_text


PROVENANCE_ROOT = REPOSITORY_ROOT / "tests/data/external_parameters"


def real_literal(value: str) -> str:
    """Return one numeric token as an explicit real64 literal."""

    normalized = value.lower().replace("d", "e")
    if "." not in normalized and "e" not in normalized:
        normalized += ".0"
    return normalized+"_real64"


def stable_format(source: str) -> str:
    """Format generated Fortran until another pass makes no changes."""

    for _ in range(4):
        formatted = format_source_text(source)
        if formatted == source:
            return source
        source = formatted
    raise RuntimeError("generated Fortran formatting did not converge")


def f2h2_module() -> str:
    """Render the released F2H2 JBKKL coefficient vector."""

    source = PROVENANCE_ROOT / "F2H2/bjkkl.par"
    released_values = source.read_text().split()
    if len(released_values) != 360:
        raise ValueError("F2H2/bjkkl.par must contain 72 five-value records")
    # The released reader requested 359 values and intentionally left the
    # final field of the last record unused.
    values = released_values[:359]
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store immutable coefficients for F2H2_GEN_JBKKL_1990.",
        "!",
        "! Generated from the released bjkkl.par model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_f2h2_jbkkl_parameters",
        "    use iso_fortran_env, only: real64",
        "    implicit none",
        "    private",
        "    public :: jbkkl_coefficients",
        "    real(real64), parameter :: jbkkl_coefficients(359) = [ &",
    ]
    literals = [real_literal(value) for value in values]
    for index in range(0, len(literals), 3):
        chunk = ", ".join(literals[index:index+3])
        suffix = ", &" if index+3 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    lines.extend([
        "        ]",
        "end module chempotpy_f2h2_jbkkl_parameters",
    ])
    return stable_format("\n".join(lines)+"\n")


def real_array(name: str, values: list[str]) -> list[str]:
    """Render one real64 parameter vector."""

    lines = [
        "    real(real64), parameter :: {}({}) = [ &".format(
            name, len(values)
        )
    ]
    literals = [real_literal(value) for value in values]
    for index in range(0, len(literals), 3):
        chunk = ", ".join(literals[index:index+3])
        suffix = ", &" if index+3 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    lines.append("        ]")
    return lines


def real_shaped_array(
    name: str, dimensions: str, values: list[str]
) -> list[str]:
    """Render one multidimensional real64 parameter array."""

    lines = [
        "    real(real64), parameter :: {}({}) = reshape([ &".format(
            name, dimensions
        )
    ]
    literals = [real_literal(value) for value in values]
    for index in range(0, len(literals), 3):
        chunk = ", ".join(literals[index:index+3])
        suffix = ", &" if index+3 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    lines.extend(["        ], [{}])".format(dimensions)])
    return lines


def hobr_module() -> str:
    """Render the released HOBr RKHS parameters."""

    source = PROVENANCE_ROOT / "HOBr/o_hbr_a_pp.par"
    values = source.read_text().split()
    if len(values) != 14+14+4*1110:
        raise ValueError("HOBr parameter record count is inconsistent")
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    records = values[28:]
    arrays = {
        "phi_oh": values[:14],
        "phi_hbr": values[14:28],
        "rkhs_x_1": records[0::4],
        "rkhs_x_2": records[1::4],
        "rkhs_y_1": records[2::4],
        "rkhs_coefficients": records[3::4],
    }
    lines = [
        "!> @brief Store immutable RKHS parameters for HOBr_3App_PIP.",
        "!",
        "! Generated from the released o_hbr_a_pp.par model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_hobr_3app_parameters",
        "    use iso_fortran_env, only: real64",
        "    implicit none",
        "    private",
        "    public :: phi_oh, phi_hbr, rkhs_x_1, rkhs_x_2, rkhs_y_1, &",
        "        & rkhs_coefficients",
    ]
    for name, array in arrays.items():
        lines.extend(real_array(name, array))
    lines.append("end module chempotpy_hobr_3app_parameters")
    return stable_format("\n".join(lines)+"\n")


def coco_module() -> str:
    """Render the released COCO neural-network parameters."""

    source = PROVENANCE_ROOT / "COCO/nn_ococ_w20.txt"
    tokens = source.read_text().split()
    values = tokens[:2492]
    if len(values) != 2492:
        raise ValueError("COCO neural-network data is incomplete")
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    offset = 0

    def take(count: int) -> list[str]:
        nonlocal offset
        result = values[offset:offset+count]
        offset += count
        return result

    arrays = (
        ("input_weights", "7,45", take(7*45)),
        ("hidden_1_biases", "45", take(45)),
        ("hidden_weights", "45,45", take(45*45)),
        ("hidden_2_biases", "45", take(45)),
        ("output_weights", "45", take(45)),
    )
    output_bias = take(1)[0]
    input_ranges = take(2*7)
    output_range = take(2)
    if offset != len(values):
        raise ValueError("COCO neural-network parameter count is inconsistent")
    lines = [
        "!> @brief Store immutable neural-network parameters for COCO_PIPNN.",
        "!",
        "! Generated from the released nn_ococ_w20.txt model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_coco_pipnn_parameters",
        "    use iso_fortran_env, only: real64",
        "    implicit none",
        "    private",
        "    public :: input_weights, hidden_1_biases, hidden_weights, &",
        "        & hidden_2_biases, output_weights, output_bias, input_ranges, &",
        "        & output_range",
    ]
    for name, dimensions, array in arrays:
        if "," in dimensions:
            lines.extend(real_shaped_array(name, dimensions, array))
        else:
            lines.extend(real_array(name, array))
    lines.append(
        "    real(real64), parameter :: output_bias = {}".format(
            real_literal(output_bias)
        )
    )
    lines.extend(real_shaped_array("input_ranges", "2,7", input_ranges))
    lines.extend(real_array("output_range", output_range))
    lines.append("end module chempotpy_coco_pipnn_parameters")
    return stable_format("\n".join(lines)+"\n")


def c3h7no_module() -> str:
    """Render the released C3H7NO PIP coefficient vector."""

    source = PROVENANCE_ROOT / "C3H7NO/coeff.dat"
    values = source.read_text().split()
    if len(values) != 8040:
        raise ValueError("C3H7NO coefficient count is inconsistent")
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store immutable coefficients for C3H7NO_PIP.",
        "!",
        "! Generated from the released coeff.dat model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_c3h7no_pip_parameters",
        "    use iso_fortran_env, only: real64",
        "    implicit none",
        "    private",
        "    public :: pip_coefficients",
    ]
    lines.extend(real_array("pip_coefficients", values))
    lines.append("end module chempotpy_c3h7no_pip_parameters")
    return stable_format("\n".join(lines)+"\n")


def generated_documents() -> dict[Path, str]:
    """Return every generated miscellaneous PES parameter module."""

    return {
        REPOSITORY_ROOT / "chempotpy/F2H2/chempotpy_f2h2_jbkkl_parameters.f90": (
            f2h2_module()
        ),
        REPOSITORY_ROOT / "chempotpy/HOBr/chempotpy_hobr_3app_parameters.f90": (
            hobr_module()
        ),
        REPOSITORY_ROOT / "chempotpy/COCO/chempotpy_coco_pipnn_parameters.f90": (
            coco_module()
        ),
        REPOSITORY_ROOT / "chempotpy/C3H7NO/chempotpy_c3h7no_pip_parameters.f90": (
            c3h7no_module()
        ),
    }


def main() -> int:
    """Write or verify the generated modules."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    failures = []
    for path, content in generated_documents().items():
        if arguments.check:
            if not path.is_file() or path.read_text() != content:
                failures.append(path)
        else:
            path.write_text(content)
            print(path.relative_to(REPOSITORY_ROOT))
    for path in failures:
        print("outdated: {}".format(path.relative_to(REPOSITORY_ROOT)))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
