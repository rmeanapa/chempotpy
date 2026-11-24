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
    REGRESSION_GEOMETRY_COUNT,
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


def standalone_sources(
    source: SurfaceSource,
    include_validation: bool = False,
) -> tuple[Path, ...]:
    """Return dependency-first sources required by one PES."""

    source_text = source.path.read_text(errors="replace")
    if include_validation:
        source_text = (
            "use chempotpy_finite_difference, "
            "only: central_finite_difference\n"+source_text
        )
    ordered = [*compile_dependencies(source.path, source_text), source.path]
    unique = []
    seen = set()
    for path in ordered:
        resolved = path.resolve()
        if resolved not in seen:
            unique.append(path)
            seen.add(resolved)
    return tuple(unique)


def combined_source(
    source: SurfaceSource,
    include_validation: bool = False,
) -> str:
    """Return one dependency-first free-form source file for a PES."""

    sections = []
    for path in standalone_sources(source, include_validation):
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
        real(wp) :: energy(nstates), gradient(nstates,natoms,3)
        real(wp) :: coupling(nstates,nstates,natoms,3)
        call pes( coordinates, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 6 and arguments[2] == "path":
        return """\
        real(wp) :: energy(nstates), gradient(nstates,natoms,3)
        real(wp) :: coupling(nstates,nstates,natoms,3)
        call pes( coordinates, 0, data_path, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 6:
        return """\
        real(wp) :: workspace(10000,3), energy(nstates)
        real(wp) :: gradient(nstates,natoms,3)
        real(wp) :: coupling(nstates,nstates,natoms,3)
        workspace = 0.0_wp; workspace(1:natoms,:) = coordinates
        call pes( workspace, natoms, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 8:
        return """\
        real(wp) :: workspace(10000,3), energy(nstates)
        real(wp) :: gradient(nstates,natoms,3)
        real(wp) :: coupling(nstates,nstates,natoms,3)
        workspace = 0.0_wp; workspace(1:natoms,:) = coordinates
        call pes( workspace, 1, 2, natoms, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "pes" and len(arguments) == 9:
        return """\
        real(wp) :: workspace(10000,3), energy(nstates)
        real(wp) :: gradient(nstates,natoms,3)
        real(wp) :: coupling(nstates,nstates,natoms,3)
        workspace = 0.0_wp; workspace(1:natoms,:) = coordinates
        call pes( workspace, 1, 1, 1, natoms, 0, energy, gradient, coupling )
        values = energy"""
    if entry_point == "dpem" and len(arguments) == 4:
        return """\
        real(wp) :: matrix(nstates,nstates)
        real(wp) :: matrix_gradient(nstates,nstates,natoms,3)
        call dpem( coordinates, 0, matrix, matrix_gradient )
        values = reshape(matrix, [nvalues])"""
    if entry_point == "dpem" and len(arguments) == 8:
        return """\
        real(wp) :: energy(nstates), gradient(nstates,natoms,3)
        real(wp) :: coupling(nstates,nstates,natoms,3)
        real(wp) :: transform(nstates,nstates), matrix(nstates,nstates)
        real(wp) :: matrix_gradient(nstates,nstates,natoms,3)
        call dpem( coordinates, 0, energy, gradient, coupling, transform, matrix, matrix_gradient )
        values = reshape(matrix, [nvalues])"""
    raise ValueError("unsupported standalone interface: {}({})".format(
        entry_point, signature
    ))


def evaluator_interface(entry_point: str, signature: str) -> str:
    """Return an explicit interface for the exported external evaluator."""

    arguments = signature.split(",")
    if entry_point == "pes" and len(arguments) == 5:
        return """\
        interface
            subroutine pes( x, igrad, energy, gradient, coupling )
                import :: wp, natoms, nstates
                integer, intent(in)   :: igrad
                real(wp), intent(in)  :: x(natoms,3)
                real(wp), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
                real(wp), intent(out) :: coupling(nstates,nstates,natoms,3)
            end subroutine pes
        end interface"""
    if entry_point == "pes" and len(arguments) == 6 and arguments[2] == "path":
        return """\
        interface
            subroutine pes( x, igrad, path, energy, gradient, coupling )
                import :: wp, natoms, nstates
                integer, intent(in)          :: igrad
                real(wp), intent(in)         :: x(natoms,3)
                character(len=*), intent(in) :: path
                real(wp), intent(out)        :: energy(nstates), gradient(nstates,natoms,3)
                real(wp), intent(out)        :: coupling(nstates,nstates,natoms,3)
            end subroutine pes
        end interface"""
    if entry_point == "pes" and len(arguments) == 6:
        return """\
        interface
            subroutine pes( workspace, atom_count, igrad, energy, gradient, coupling )
                import :: wp, natoms, nstates
                integer, intent(in)   :: atom_count, igrad
                real(wp), intent(in)  :: workspace(10000,3)
                real(wp), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
                real(wp), intent(out) :: coupling(nstates,nstates,natoms,3)
            end subroutine pes
        end interface"""
    if entry_point == "pes" and len(arguments) in {8, 9}:
        species_arguments = len(arguments)-6
        species_names = ["species_{}".format(index) for index in range(1, species_arguments+1)]
        names = ", ".join(species_names)
        return """\
        interface
            subroutine pes( workspace, {names}, atom_count, igrad, energy, gradient, coupling )
                import :: wp, natoms, nstates
                integer, intent(in)   :: {names}, atom_count, igrad
                real(wp), intent(in)  :: workspace(10000,3)
                real(wp), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
                real(wp), intent(out) :: coupling(nstates,nstates,natoms,3)
            end subroutine pes
        end interface""".format(names=names)
    if entry_point == "dpem" and len(arguments) == 4:
        return """\
        interface
            subroutine dpem( x, igrad, matrix, matrix_gradient )
                import :: wp, natoms, nstates
                integer, intent(in)   :: igrad
                real(wp), intent(in)  :: x(natoms,3)
                real(wp), intent(out) :: matrix(nstates,nstates)
                real(wp), intent(out) :: matrix_gradient(nstates,nstates,natoms,3)
            end subroutine dpem
        end interface"""
    if entry_point == "dpem" and len(arguments) == 8:
        return """\
        interface
            subroutine dpem( x, igrad, energy, gradient, coupling, transform, matrix, matrix_gradient )
                import :: wp, natoms, nstates
                integer, intent(in)   :: igrad
                real(wp), intent(in)  :: x(natoms,3)
                real(wp), intent(out) :: energy(nstates), gradient(nstates,natoms,3)
                real(wp), intent(out) :: coupling(nstates,nstates,natoms,3)
                real(wp), intent(out) :: transform(nstates,nstates), matrix(nstates,nstates)
                real(wp), intent(out) :: matrix_gradient(nstates,nstates,natoms,3)
            end subroutine dpem
        end interface"""
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
        reference_index = (
            (geometry_index-1) % REGRESSION_GEOMETRY_COUNT
        )+1
        geometry_blocks.append(
            "    coordinates = 0.0_wp\n{}\n"
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
    use chempotpy_finite_difference, only: central_finite_difference
    use chempotpy_kinds,             only: wp
    use chempotpy_status,            only: chempotpy_status_t
    implicit none
    integer, parameter :: natoms  = {atom_count}
    integer, parameter :: nstates = {state_count}
    integer, parameter :: nvalues = {value_count}
    real(wp), parameter :: difference_step = {difference_step:.17e}_wp
    real(wp)            :: coordinates(natoms,3)
{data_declaration}{data_setup}{geometry_blocks}
    contains
    subroutine evaluate_values( coordinates, values, status )
        implicit none
        real(wp), intent(in)                    :: coordinates(:,:)
        real(wp), intent(out)                   :: values(:)
        type(chempotpy_status_t), intent(inout) :: status
{evaluator_interface}
{evaluator_body}
        if( status%failed() ) values = 0.0_wp
    end subroutine evaluate_values

    subroutine report_geometry( geometry_index, coordinates )
        implicit none
        integer, intent(in)                    :: geometry_index
        real(wp), intent(in)                   :: coordinates(natoms,3)
        real(wp)                               :: base(nvalues)
        real(wp)                               :: numerical_gradient(nvalues,3*natoms)
        real(wp)                               :: numerical_hessian(nvalues,3*natoms,3*natoms)
        type(chempotpy_status_t)               :: status
        integer                                :: coordinate_i, coordinate_j, value_index
        call central_finite_difference(evaluate_values, coordinates, difference_step, base, numerical_gradient, &
            numerical_hessian, status)
        if( status%failed() )then
            write(*,'(A,1X,I0,2A)') 'geometry', geometry_index, ' failed: ', trim(status%message)
            return
        endif
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
        evaluator_interface=evaluator_interface(entry_point, signature),
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
    source_output.write_text(combined_source(
        source, include_validation=with_driver
    ))
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
