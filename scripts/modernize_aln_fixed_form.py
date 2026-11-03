#!/usr/bin/env python3
"""Convert the remaining Aln fixed-form kernels to free-form sources."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


PROCEDURE_PATTERN = re.compile(
    r"^(?P<indent>\s*)(?:(?:double\s+precision|real(?:\*\d+)?|integer|"
    r"logical|complex)(?:\s+kind\s*=\s*\d+)?\s+)?"
    r"(?P<kind>subroutine|function)\s+(?P<name>[a-z][a-z0-9_]*)\s*\(",
    flags=re.IGNORECASE,
)


@dataclass
class Procedure:
    """Track the external procedure currently being converted."""

    kind: str
    name: str


def modern_wrapper() -> str:
    """Return the common core-aware Aln PES wrapper."""

    return """!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Evaluate an aluminum-cluster potential through the Python contract.
subroutine pes( coordinates_angstrom, atom_count, derivative_order, &
    energy_ev, gradient_ev_per_angstrom, coupling_per_angstrom )
    use chempotpy_legacy_units, only: legacy_coordinates_to_bohr, &
        legacy_energy_to_electronvolt
    use chempotpy_kinds, only: real64
    implicit none

    integer, parameter :: state_count = 1
    real(real64), intent(in) :: coordinates_angstrom(10000,3)
    integer, intent(in) :: atom_count, derivative_order
    real(real64), intent(out) :: energy_ev(state_count)
    real(real64), intent(out) :: gradient_ev_per_angstrom( &
        state_count,atom_count,3 )
    real(real64), intent(out) :: coupling_per_angstrom( &
        state_count,state_count,atom_count,3 )
    real(real64) :: coordinates_bohr(atom_count,3), energy_hartree

    energy_ev = 0.0d0
    gradient_ev_per_angstrom = 0.0d0
    coupling_per_angstrom = 0.0d0
    call legacy_coordinates_to_bohr( &
        coordinates_angstrom(1:atom_count,:), coordinates_bohr )
    call pot( coordinates_bohr(:,1), coordinates_bohr(:,2), &
        coordinates_bohr(:,3), energy_hartree, atom_count, atom_count )
    if( derivative_order == 0 ) then
        energy_ev(1) = legacy_energy_to_electronvolt(energy_hartree)
    else
        write(*,*) 'Only energy is available'
    endif
end subroutine pes
"""


def first_procedure_end(lines: List[str]) -> int:
    """Return the index following the released PES wrapper."""

    for index, line in enumerate(lines):
        if re.match(r"^\s{6}end(?:subroutine)?\s*$", line, re.IGNORECASE):
            return index + 1
    raise ValueError("released PES wrapper has no end statement")


def fixed_line_to_free(line: str) -> tuple[str, bool, bool]:
    """Convert one fixed-form physical line without changing its tokens."""

    expanded = line.expandtabs(8).rstrip("\n")
    if not expanded.strip():
        return "", False, False
    if expanded[0] in "cC*!":
        return "!" + expanded[1:].rstrip(), False, True
    padded = expanded.ljust(6)
    label = padded[:5].strip()
    continuation = padded[5] not in " 0"
    statement = padded[6:72].rstrip()
    if not statement:
        return "", continuation, False
    leading = len(statement) - len(statement.lstrip())
    if label:
        converted = label + " " + statement.lstrip()
    elif continuation:
        converted = "        & " + statement.lstrip()
    else:
        indentation = 0 if PROCEDURE_PATTERN.match(statement) else 4
        if leading:
            indentation += 4*(leading//4)
        converted = " "*indentation + statement.lstrip()
    return converted, continuation, False


def space_procedure_arguments(line: str, declaration: bool) -> str:
    """Apply the project argument spacing to a procedure declaration line."""

    if not declaration:
        return line
    line = re.sub(r",\s*", ", ", line)
    match = PROCEDURE_PATTERN.match(line)
    if match is not None:
        opening = line.find("(", match.start("name"))
        if line[opening:opening + 2] != "( ":
            line = line[:opening + 1] + " " + line[opening + 1:]
    closing = line.rfind(")")
    if closing >= 0 and line[closing - 1:closing] != " ":
        line = line[:closing] + " " + line[closing:]
    return line


def continue_free_form_line(line: str) -> str:
    """Add a continuation marker before any trailing inline comment."""

    comment = line.find("!")
    if comment < 0:
        return line.rstrip() + " &"
    code = line[:comment].rstrip()
    trailing_comment = line[comment:]
    return code + " &  " + trailing_comment


def convert_kernel(lines: List[str], surface_name: str) -> str:
    """Convert all helper procedures following the released wrapper."""

    output: List[str] = []
    last_code_index: Optional[int] = None
    procedure: Optional[Procedure] = None
    declaration_continues = False
    for original in lines:
        converted, continuation, comment = fixed_line_to_free(original)
        if continuation and last_code_index is not None:
            output[last_code_index] = continue_free_form_line(
                output[last_code_index]
            )
        if not converted:
            output.append("")
            continue
        match = PROCEDURE_PATTERN.match(converted)
        if match is not None:
            procedure = Procedure(
                kind=match.group("kind").lower(),
                name=match.group("name"),
            )
            output.append(
                "!> @brief Evaluate {} for the {} potential.".format(
                    procedure.name, surface_name
                )
            )
            declaration_continues = ")" not in converted
            converted = space_procedure_arguments(converted, True)
        elif declaration_continues:
            converted = space_procedure_arguments(converted, True)
            declaration_continues = ")" not in converted
        if procedure is not None and converted.strip().lower() in {
            "end", "endsubroutine", "end subroutine", "endfunction",
            "end function",
        }:
            converted = "end {} {}".format(procedure.kind, procedure.name)
            procedure = None
        output.append(converted)
        if not comment:
            last_code_index = len(output) - 1
    while output and not output[-1]:
        output.pop()
    return "\n".join(output) + "\n"


def convert_source(source: Path) -> Path:
    """Create a free-form sibling of one released Aln source."""

    lines = source.read_text(errors="strict").splitlines()
    wrapper_end = first_procedure_end(lines)
    destination = source.with_suffix(".f90")
    destination.write_text(
        modern_wrapper()
        + "\n"
        + convert_kernel(lines[wrapper_end:], source.stem)
    )
    return destination


def main() -> int:
    """Convert all requested fixed-form Aln sources."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    args = parser.parse_args()
    for source in args.sources:
        destination = convert_source(source)
        print("{} -> {}".format(source, destination))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
