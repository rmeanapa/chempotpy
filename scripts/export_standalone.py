#!/usr/bin/env python3
"""Export one PES and its bundled dependencies as standalone Fortran."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import SurfaceSource, discover_surface_sources
from scripts.capture_library_references import (
    assignments,
    compile_dependencies,
    first_integer_parameter,
    interface_signature,
)


def select_surface(selector: str) -> SurfaceSource:
    """Resolve ``SYSTEM/SURFACE`` or one unambiguous surface name."""

    normalized = selector.strip()
    if normalized.endswith(".f90"):
        normalized = normalized[:-4]
    matches = []
    for source in discover_surface_sources():
        if not source.is_surface:
            continue
        qualified = "{}/{}".format(source.system_name, source.surface_name)
        if normalized in {qualified, source.surface_name}:
            matches.append(source)
    if not matches:
        raise ValueError("unknown callable surface: {}".format(selector))
    if len(matches) > 1:
        choices = ", ".join(
            "{}/{}".format(item.system_name, item.surface_name)
            for item in matches
        )
        raise ValueError(
            "surface name is ambiguous; use SYSTEM/SURFACE: {}".format(choices)
        )
    return matches[0]


def surface_interface(source: SurfaceSource) -> tuple[str, str, int, int]:
    """Return entry point, signature, atom count, and state count."""

    source_text = source.path.read_text(errors="replace")
    entry_point, signature = interface_signature(source_text)
    arguments = signature.split(",")
    if len(arguments) in {8, 9}:
        atom_count = 3
    elif entry_point == "pes" and len(arguments) == 6 and arguments[2] != "path":
        atom_count = 5
    else:
        atom_count = first_integer_parameter(source_text, "natoms", "atom_count")
    state_count = first_integer_parameter(source_text, "nstates", "state_count")
    return entry_point, signature, atom_count, state_count


def standalone_sources(source: SurfaceSource) -> tuple[Path, ...]:
    """Return dependency-first sources required by one PES."""

    source_text = source.path.read_text(errors="replace")
    ordered = [*compile_dependencies(source.path, source_text), source.path]
    unique = []
    seen = set()
    for path in ordered:
        resolved = path.resolve()
        if resolved not in seen:
            unique.append(path)
            seen.add(resolved)
    return tuple(unique)


def combined_source(source: SurfaceSource) -> str:
    """Return one dependency-first free-form source file for a PES."""

    sections = []
    for path in standalone_sources(source):
        relative = path.relative_to(REPOSITORY_ROOT)
        sections.append("! Standalone source: {}\n{}".format(
            relative, path.read_text(errors="strict").rstrip()
        ))
    return "\n\n".join(sections) + "\n"


def evaluator_body(entry_point: str, signature: str) -> str:
    """Return declarations and a call for the historical public interface."""

    arguments = signature.split(",")
    if entry_point == "pes" and len(arguments) == 5:
        return """\
        real(8) :: energy(nstates), gradient(nstates,natoms,3)
        real(8) :: coupling(nstates,nstates,natoms,3)
        call pes( coordinates, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 6 and arguments[2] == "path":
        return """\
        real(8) :: energy(nstates), gradient(nstates,natoms,3)
        real(8) :: coupling(nstates,nstates,natoms,3)
        call pes( coordinates, 0, data_path, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 6:
        return """\
        real(8) :: workspace(10000,3), energy(nstates)
        real(8) :: gradient(nstates,natoms,3)
        real(8) :: coupling(nstates,nstates,natoms,3)
        workspace = 0.0d0; workspace(1:natoms,:) = coordinates
        call pes( workspace, natoms, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 8:
        return """\
        real(8) :: workspace(10000,3), energy(nstates)
        real(8) :: gradient(nstates,natoms,3)
        real(8) :: coupling(nstates,nstates,natoms,3)
        workspace = 0.0d0; workspace(1:natoms,:) = coordinates
        call pes( workspace, 1, 2, natoms, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 9:
        return """\
        real(8) :: workspace(10000,3), energy(nstates)
        real(8) :: gradient(nstates,natoms,3)
        real(8) :: coupling(nstates,nstates,natoms,3)
        workspace = 0.0d0; workspace(1:natoms,:) = coordinates
        call pes( workspace, 1, 1, 1, natoms, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "dpem" and len(arguments) == 4:
        return """\
        real(8) :: matrix(nstates,nstates)
        real(8) :: matrix_gradient(nstates,nstates,natoms,3)
        call dpem( coordinates, 0, matrix, matrix_gradient )
        values = reshape(matrix, [nvalues])"""
    if entry_point == "dpem" and len(arguments) == 8:
        return """\
        real(8) :: energy(nstates), gradient(nstates,natoms,3)
        real(8) :: coupling(nstates,nstates,natoms,3)
        real(8) :: transform(nstates,nstates), matrix(nstates,nstates)
        real(8) :: matrix_gradient(nstates,nstates,natoms,3)
        call dpem( coordinates, 0, energy, gradient, coupling, transform, matrix, matrix_gradient )
        values = reshape(matrix, [nvalues])"""
    raise ValueError("unsupported standalone interface: {}({})".format(
        entry_point, signature
    ))


def example_driver(
    source: SurfaceSource,
    geometry_count: int,
    difference_step: float,
) -> str:
    """Return a program reporting values and finite-difference derivatives."""

    entry_point, signature, atom_count, state_count = surface_interface(source)
    value_count = state_count if entry_point == "pes" else state_count*state_count
    path_argument = "path" in signature.split(",")
    geometry_blocks = []
    for geometry_index in range(1, geometry_count + 1):
        reference_index = 1 if geometry_index % 2 else 2
        geometry_blocks.append(
            "    coordinates = 0.0d0\n{}\n"
            "    call report_geometry({}, coordinates)".format(
                assignments(atom_count, reference_index, "coordinates"),
                geometry_index,
            )
        )
    data_declaration = "    character(len=1024) :: data_path\n" if path_argument else ""
    data_setup = ""
    if path_argument:
        data_setup = """\
    data_path = '.'
    if( command_argument_count() >= 1 ) call get_command_argument(1, data_path)
    if( len_trim(data_path) > 0 )then
        if( data_path(len_trim(data_path):len_trim(data_path)) /= '/' ) data_path = trim(data_path)//'/'
    endif
"""
    value_label = "energy" if entry_point == "pes" else "diabatic_matrix_element"
    return """\
program chempotpy_standalone_example
    implicit none
    integer, parameter :: natoms  = {atom_count}
    integer, parameter :: nstates = {state_count}
    integer, parameter :: nvalues = {value_count}
    real(8), parameter :: difference_step = {difference_step:.17e}
    real(8)            :: coordinates(natoms,3)
{data_declaration}{data_setup}{geometry_blocks}
    contains
    subroutine evaluate_values( coordinates, values )
        implicit none
        real(8), intent(in)  :: coordinates(natoms,3)
        real(8), intent(out) :: values(nvalues)
{evaluator_body}
    end subroutine evaluate_values

    subroutine report_geometry( geometry_index, coordinates )
        implicit none
        integer, intent(in) :: geometry_index
        real(8), intent(in) :: coordinates(natoms,3)
        real(8)             :: base(nvalues), forward(nvalues), backward(nvalues)
        real(8)             :: mixed_pp(nvalues), mixed_pm(nvalues)
        real(8)             :: mixed_mp(nvalues), mixed_mm(nvalues)
        real(8)             :: numerical_gradient(nvalues,3*natoms)
        real(8)             :: numerical_hessian(nvalues,3*natoms,3*natoms)
        real(8)             :: displaced(natoms,3)
        integer             :: atom_i, atom_j, coordinate_i, coordinate_j
        integer             :: direction_i, direction_j, value_index
        call evaluate_values(coordinates, base)
        numerical_gradient = 0.0d0; numerical_hessian = 0.0d0
        do coordinate_i=1, 3*natoms
            atom_i      = (coordinate_i-1)/3+1
            direction_i = mod(coordinate_i-1, 3)+1
            displaced = coordinates
            displaced(atom_i,direction_i) = displaced(atom_i,direction_i)+difference_step
            call evaluate_values(displaced, forward)
            displaced(atom_i,direction_i) = displaced(atom_i,direction_i)-2.0d0*difference_step
            call evaluate_values(displaced, backward)
            numerical_gradient(:,coordinate_i) = (forward-backward)/(2.0d0*difference_step)
            numerical_hessian(:,coordinate_i,coordinate_i) = &
                (forward-2.0d0*base+backward)/(difference_step*difference_step)
            do coordinate_j=coordinate_i+1, 3*natoms
                atom_j      = (coordinate_j-1)/3+1
                direction_j = mod(coordinate_j-1, 3)+1
                displaced = coordinates
                displaced(atom_i,direction_i) = displaced(atom_i,direction_i)+difference_step
                displaced(atom_j,direction_j) = displaced(atom_j,direction_j)+difference_step
                call evaluate_values(displaced, mixed_pp)
                displaced(atom_j,direction_j) = displaced(atom_j,direction_j)-2.0d0*difference_step
                call evaluate_values(displaced, mixed_pm)
                displaced(atom_i,direction_i) = displaced(atom_i,direction_i)-2.0d0*difference_step
                call evaluate_values(displaced, mixed_mm)
                displaced(atom_j,direction_j) = displaced(atom_j,direction_j)+2.0d0*difference_step
                call evaluate_values(displaced, mixed_mp)
                numerical_hessian(:,coordinate_i,coordinate_j) = &
                    (mixed_pp-mixed_pm-mixed_mp+mixed_mm)/(4.0d0*difference_step*difference_step)
                numerical_hessian(:,coordinate_j,coordinate_i) = numerical_hessian(:,coordinate_i,coordinate_j)
            enddo
        enddo
        write(*,'(A,1X,I0)') 'geometry', geometry_index
        do value_index=1, nvalues
            write(*,'(A,1X,I0,1X,ES26.17E3)') '{value_label}', value_index, base(value_index)
            do coordinate_i=1, 3*natoms
                write(*,'(A,2(1X,I0),1X,ES26.17E3)') &
                    'numerical_gradient', value_index, coordinate_i, numerical_gradient(value_index,coordinate_i)
                do coordinate_j=1, 3*natoms
                    write(*,'(A,3(1X,I0),1X,ES26.17E3)') 'numerical_hessian', value_index, coordinate_i, &
                        coordinate_j, numerical_hessian(value_index,coordinate_i,coordinate_j)
                enddo
            enddo
        enddo
    end subroutine report_geometry
end program chempotpy_standalone_example
""".format(
        atom_count=atom_count,
        state_count=state_count,
        value_count=value_count,
        difference_step=difference_step,
        data_declaration=data_declaration,
        data_setup=data_setup,
        geometry_blocks="\n".join(geometry_blocks),
        evaluator_body=evaluator_body(entry_point, signature),
        value_label=value_label,
    )


def export_surface(
    source: SurfaceSource,
    output_directory: Path,
    with_driver: bool,
    geometry_count: int,
    difference_step: float,
) -> tuple[Path, ...]:
    """Write the requested standalone artifacts and return their paths."""

    output_directory.mkdir(parents=True, exist_ok=True)
    source_output = output_directory / (source.surface_name + "_standalone.f90")
    source_output.write_text(combined_source(source))
    outputs = [source_output]
    if with_driver:
        driver_output = output_directory / (source.surface_name + "_example.f90")
        driver_output.write_text(example_driver(
            source, geometry_count, difference_step
        ))
        outputs.append(driver_output)
    return tuple(outputs)


def main() -> int:
    """Run the standalone exporter."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--surface", required=True)
    parser.add_argument("--output", type=Path, default=Path("standalone"))
    parser.add_argument("--with-driver", action="store_true")
    parser.add_argument("--geometries", type=int, default=2)
    parser.add_argument("--difference-step", type=float, default=1.0e-3)
    args = parser.parse_args()
    if args.geometries < 1:
        parser.error("--geometries must be positive")
    if args.difference_step <= 0.0:
        parser.error("--difference-step must be positive")
    try:
        source = select_surface(args.surface)
        outputs = export_surface(
            source,
            args.output,
            args.with_driver,
            args.geometries,
            args.difference_step,
        )
    except ValueError as error:
        parser.error(str(error))
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
