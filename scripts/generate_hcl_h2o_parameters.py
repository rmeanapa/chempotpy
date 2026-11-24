#!/usr/bin/env python3
"""Generate immutable HCl-H2O PIP parameters from the POTLIB release."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests" / "data" / "external_parameters" / "HClH2O"
OUTPUT = ROOT / "chempotpy" / "H3ClO" / "chempotpy_hcl_h2o_pip_parameters.f90"


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one provenance file."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(name: str) -> list[str]:
    """Read the coefficient payload from one released PCF file."""

    lines = (SOURCE/name).read_text().splitlines()
    count = int(lines[1].split()[0])
    values = " ".join(lines[2:]).replace("D", "e").replace("d", "e").split()
    if len(values) != count:
        raise ValueError(f"{name} declares {count} coefficients but contains {len(values)}")
    return values


def real_array(name: str, values: list[str], chunk: int = 4) -> str:
    """Render one real parameter array."""

    lines = [f"    real(wp), parameter :: {name}({len(values)}) = [ &"]
    for start in range(0, len(values), chunk):
        group = ", ".join(value+"_wp" for value in values[start:start+chunk])
        suffix = ", &" if start+chunk < len(values) else " ]"
        lines.append("        & "+group+suffix)
    return "\n".join(lines)


def generated_document() -> str:
    """Return the surface-owned immutable coefficient module."""

    files = (
        "pcf-x1.dat", "pcf-y1.dat", "pcf-z1.dat", "pcf-x2.dat",
        "pcf-x1y1.dat", "pcf-x1z1.dat", "pcf-y1z1.dat", "pcf-x3y1z1.dat",
    )
    arrays = {
        "hh_coefficients": coefficients("pcf-x2.dat"),
        "ho_coefficients": coefficients("pcf-x1y1.dat"),
        "hcl_coefficients": coefficients("pcf-x1z1.dat"),
        "ocl_coefficients": coefficients("pcf-y1z1.dat"),
        "five_body_coefficients": coefficients("pcf-x3y1z1.dat"),
    }
    output = [
        "!> @brief Store immutable Bowman-Mancini HCl-H2O PIP coefficients.",
        "!",
        "! Generated from the released PCF files; do not edit numeric values by hand.",
    ]
    output.extend(f"! {name} SHA-256: {sha256(SOURCE/name)}" for name in files)
    output.extend([
        "module chempotpy_hcl_h2o_pip_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: five_body_coefficients, hcl_coefficients, hh_coefficients, ho_coefficients, ocl_coefficients",
        "    real(wp), parameter :: pair_cutoff_radius = 15.0_wp",
        "    real(wp), parameter :: pair_morse_length = 2.0_wp",
        "    real(wp), parameter :: five_body_cutoff_radius = 1.0e8_wp",
        "    public :: five_body_cutoff_radius, pair_cutoff_radius, pair_morse_length",
    ])
    output.extend(real_array(name, values) for name, values in arrays.items())
    output.extend(["end module chempotpy_hcl_h2o_pip_parameters", ""])
    return "\n".join(output)


def main() -> int:
    """Write or verify the surface-owned immutable coefficient module."""

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
