#!/usr/bin/env python3
"""Generate immutable H3O2- sGDML parameters from the released model."""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
import hashlib
from pathlib import Path
import struct
import zipfile


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = (
    REPOSITORY_ROOT
    / "tests/data/external_parameters/H3O2m/sgdml_model.npz"
)
OUTPUT_PATH = (
    REPOSITORY_ROOT
    / "chempotpy/H3O2m/chempotpy_h3o2m_sgdml_parameters.f90"
)
CHUNK_SIZE = 300


@dataclass(frozen=True)
class NpyArray:
    """Store one numeric NPY member without requiring NumPy."""

    shape: tuple[int, ...]
    fortran_order: bool
    values: tuple[float | int, ...]

    def at(self, *indices: int) -> float | int:
        """Return one element using the member's declared storage order."""

        if len(indices) != len(self.shape):
            raise IndexError("incorrect NPY index rank")
        offset = 0
        stride = 1
        dimensions = range(len(self.shape))
        if not self.fortran_order:
            dimensions = reversed(tuple(dimensions))
        for dimension in dimensions:
            index = indices[dimension]
            if index < 0 or index >= self.shape[dimension]:
                raise IndexError("NPY index is outside the array")
            offset += index*stride
            stride *= self.shape[dimension]
        return self.values[offset]


def read_npy_member(archive: zipfile.ZipFile, name: str) -> NpyArray:
    """Read one little-endian float64 or int64 array from an NPZ archive."""

    payload = archive.read(name+".npy")
    if payload[:6] != b"\x93NUMPY":
        raise ValueError("{} is not an NPY member".format(name))
    major, minor = payload[6], payload[7]
    del minor
    if major == 1:
        header_length = struct.unpack_from("<H", payload, 8)[0]
        header_start = 10
    elif major in {2, 3}:
        header_length = struct.unpack_from("<I", payload, 8)[0]
        header_start = 12
    else:
        raise ValueError("unsupported NPY version {}".format(major))
    header = ast.literal_eval(
        payload[header_start:header_start+header_length].decode("latin1")
    )
    descriptor = header["descr"]
    formats = {"<f8": "d", "<i8": "q", "|i8": "q"}
    if descriptor not in formats:
        raise ValueError("unsupported NPY descriptor {}".format(descriptor))
    shape = tuple(int(value) for value in header["shape"])
    count = 1
    for dimension in shape:
        count *= dimension
    data_start = header_start+header_length
    values = struct.unpack_from(
        "<{}{}".format(count, formats[descriptor]), payload, data_start
    )
    return NpyArray(shape, bool(header["fortran_order"]), values)


def wp_literal(value: float) -> str:
    """Return a round-trip-safe working-precision literal."""

    return "{:.17e}_wp".format(value)


def render_chunks(prefix: str, values: list[float]) -> tuple[list[str], list[str]]:
    """Render standard-conforming parameter chunks and return their names."""

    lines: list[str] = []
    names: list[str] = []
    for chunk_index, start in enumerate(range(0, len(values), CHUNK_SIZE), 1):
        chunk = values[start:start+CHUNK_SIZE]
        name = "{}_{:03d}".format(prefix, chunk_index)
        names.append(name)
        lines.append(
            "    real(wp), parameter :: {}({}) = [ &".format(
                name, len(chunk)
            )
        )
        literals = [wp_literal(value) for value in chunk]
        for index in range(0, len(literals), 3):
            record = ", ".join(literals[index:index+3])
            suffix = ", &" if index+3 < len(literals) else " &"
            lines.append("        & "+record+suffix)
        lines.append("        ]")
    return lines, names


def render_joined_array(
    name: str,
    dimensions: str,
    chunk_names: list[str],
) -> list[str]:
    """Join parameter chunks into one shaped public array."""

    lines = [
        "    real(wp), parameter :: {}({}) = reshape([ &".format(
            name, dimensions
        )
    ]
    for index, chunk_name in enumerate(chunk_names):
        suffix = ", &" if index+1 < len(chunk_names) else " &"
        lines.append("        & "+chunk_name+suffix)
    lines.append("        ], [{}])".format(dimensions))
    return lines


def generated_source() -> str:
    """Return the immutable parameter module generated from the release."""

    with zipfile.ZipFile(MODEL_PATH) as archive:
        descriptors = read_npy_member(archive, "R_desc")
        coefficients = read_npy_member(archive, "R_d_desc_alpha")
        permutations = read_npy_member(archive, "tril_perms_lin")
        sigma = read_npy_member(archive, "sig")
        energy_offset = read_npy_member(archive, "c")
        energy_scale = read_npy_member(archive, "std")
    if descriptors.shape != (10, 3000):
        raise ValueError("unexpected H3O2- descriptor dimensions")
    if coefficients.shape != (3000, 10):
        raise ValueError("unexpected H3O2- coefficient dimensions")
    if permutations.shape != (120,):
        raise ValueError("unexpected H3O2- permutation dimensions")
    descriptor_values = [
        float(descriptors.at(descriptor, training))
        for training in range(3000)
        for descriptor in range(10)
    ]
    coefficient_values = [
        float(coefficients.at(training, descriptor))
        for training in range(3000)
        for descriptor in range(10)
    ]
    permutation_values = [
        int(permutations.at(descriptor*12+permutation)) % 10+1
        for permutation in range(12)
        for descriptor in range(10)
    ]
    descriptor_lines, descriptor_chunks = render_chunks(
        "descriptor_chunk", descriptor_values
    )
    coefficient_lines, coefficient_chunks = render_chunks(
        "coefficient_chunk", coefficient_values
    )
    digest = hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store the immutable released H3O2- sGDML model parameters.",
        "!",
        "! Generated from the released authors' h3o2m_traindata-unknown-train3000-sym12.npz.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_h3o2m_sgdml_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: descriptor_coefficients, descriptor_permutations, energy_offset, energy_scale",
        "    public :: kernel_sigma, training_descriptors",
        "    integer, parameter :: descriptor_count = 10",
        "    integer, parameter :: permutation_count = 12",
        "    integer, parameter :: training_count = 3000",
        "    real(wp), parameter :: kernel_sigma = {}".format(
            wp_literal(float(sigma.values[0]))
        ),
        "    real(wp), parameter :: energy_offset = {}".format(
            wp_literal(float(energy_offset.values[0]))
        ),
        "    real(wp), parameter :: energy_scale = {}".format(
            wp_literal(float(energy_scale.values[0]))
        ),
        "    integer, parameter :: descriptor_permutations(descriptor_count,permutation_count) = reshape([ &",
    ]
    for index in range(0, len(permutation_values), 12):
        values = ", ".join(
            str(value) for value in permutation_values[index:index+12]
        )
        suffix = ", &" if index+12 < len(permutation_values) else " &"
        lines.append("        & "+values+suffix)
    lines.append("        ], [descriptor_count, permutation_count])")
    lines.extend(descriptor_lines)
    lines.extend(coefficient_lines)
    lines.extend(render_joined_array(
        "training_descriptors",
        "descriptor_count,training_count",
        descriptor_chunks,
    ))
    lines.extend(render_joined_array(
        "descriptor_coefficients",
        "descriptor_count,training_count",
        coefficient_chunks,
    ))
    lines.append("end module chempotpy_h3o2m_sgdml_parameters")
    return "\n".join(lines)+"\n"


def main() -> int:
    """Write the generated module or verify its committed contents."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    rendered = generated_source()
    if arguments.check:
        if not OUTPUT_PATH.is_file() or OUTPUT_PATH.read_text() != rendered:
            print("outdated: {}".format(OUTPUT_PATH.relative_to(REPOSITORY_ROOT)))
            return 1
    else:
        OUTPUT_PATH.write_text(rendered)
        print(OUTPUT_PATH.relative_to(REPOSITORY_ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
