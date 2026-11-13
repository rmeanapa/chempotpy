"""Tests for the maintained modern Fortran formatter."""

import unittest
from pathlib import Path

from scripts.format_fortran import (
    fold_standalone_dimensions,
    format_source_text,
    hoist_pes_metadata,
    normalize_contains_blank_lines,
    normalize_procedure_separator_blank_lines,
    remove_dead_commented_code,
    specialize_pes_documentation,
)


class FortranFormatMigrationTests(unittest.TestCase):
    def test_standalone_dimensions_are_folded_into_declarations(self) -> None:
        source = """\
subroutine switching( switch, dswitch )
    implicit none
    double precision :: dswitch, switch
    integer :: indices, count
    dimension switch(4), dswitch(4, 4)
    dimension indices(2), count(3, 2)
end subroutine switching
"""
        converted = format_source_text(source)
        self.assertIn(
            "double precision :: dswitch(4,4), switch(4)", converted
        )
        self.assertIn("integer          :: indices(2), count(3,2)", converted)
        self.assertNotRegex(converted, r"(?im)^\s*dimension\b")
        focused = fold_standalone_dimensions(source)
        self.assertIn(
            "    double precision :: dswitch(4,4), switch(4)", focused
        )
        self.assertIn("    integer :: indices(2), count(3,2)", focused)

        parameter_shape = """\
subroutine evaluate
    implicit none
    double precision :: values
    integer, parameter :: value_count = 3
    dimension values(value_count)
end subroutine evaluate
"""
        focused = fold_standalone_dimensions(parameter_shape)
        self.assertLess(
            focused.index("integer, parameter :: value_count = 3"),
            focused.index("double precision :: values(value_count)"),
        )

    def test_dimension_attributes_are_folded_into_entities(self) -> None:
        source = """\
subroutine evaluate( x, work )
    implicit none
    integer, parameter, dimension(0:2) :: bounds = [0, 2, 4]
    real(wp), dimension(natoms, 3), intent(in) :: x
    real(wp), dimension(:, :), allocatable :: work, scratch
    real(wp), dimension(size(x,2)*(size(x,2)-1)/2) :: distances
end subroutine evaluate
"""
        converted = fold_standalone_dimensions(source)
        self.assertIn(
            "integer, parameter :: bounds(0:2) = [0, 2, 4]", converted
        )
        self.assertIn("real(wp), intent(in) :: x(natoms,3)", converted)
        self.assertIn(
            "real(wp), allocatable :: work(:,:), scratch(:,:)", converted
        )
        self.assertIn(
            "real(wp) :: distances(size(x,2)*(size(x,2)-1)/2)", converted
        )
        self.assertNotRegex(converted, r"(?i)\bdimension\s*\(")

    def test_dimension_initializer_keeps_physical_continuations(self) -> None:
        source = """\
subroutine evaluate
    implicit none
    real(wp), dimension(0:3) :: coefficients = (/ &
        & 1.0, &
        & 2.0, &
        & 3.0, &
        & 4.0 /)
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("real(wp) :: coefficients(0:3) = (/ &", converted)
        self.assertIn("\n        & 1.0, &\n", converted)
        self.assertIn("\n        & 2.0, &\n", converted)
        self.assertIn("\n        & 3.0, &\n", converted)
        self.assertIn("\n        & 4.0 /)\n", converted)

    def test_redundant_terminal_return_is_removed(self) -> None:
        source = """\
subroutine report( value )
    implicit none
    double precision :: value
    value = value+1.d0
    return
end subroutine report
"""
        converted = format_source_text(source)
        self.assertNotRegex(
            converted.lower(), r"\breturn\s*\n\s*end subroutine"
        )

    def test_early_return_is_preserved(self) -> None:
        source = """\
subroutine evaluate( valid )
    implicit none
    logical :: valid
    if (.not. valid) return
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("if( .not. valid ) return", converted)

    def test_intent_declarations_precede_local_variables(self) -> None:
        source = """\
subroutine evaluate( values, count )
    implicit none

    integer :: index
    integer, parameter :: maximum_count = 10
    integer, intent(in) :: count
    real, intent(out) :: values(maximum_count)

    values = 0.0
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertLess(
            converted.index("integer, intent(in)"),
            converted.index(":: index"),
        )
        self.assertLess(
            converted.index("integer, parameter"),
            converted.index("integer, intent(in)"),
        )
        self.assertNotIn("maximum_count = 10\n\n", converted)
        self.assertNotIn("values(maximum_count)\n\n", converted)

    def test_blank_lines_are_removed_inside_procedures(self) -> None:
        source = """\
subroutine evaluate


    implicit none



end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertNotIn("\n\n", converted)

    def test_uninformative_comment_lines_are_removed(self) -> None:
        source = """\
subroutine evaluate
    implicit none
    !
    !--------------------
    ! ================== !
    ! Keep this scientific explanation.
    ! \\       /
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertNotIn("!--------------------", converted)
        self.assertNotIn("! ================== !", converted)
        self.assertIn("! Keep this scientific explanation.", converted)
        self.assertIn("! \\       /", converted)

    def test_commented_legacy_code_is_removed(self) -> None:
        source = """\
module example
contains
    subroutine evaluate
        implicit none
        real :: value
        value = 1.0
        !      if( value > 0.0 )then
        !          value = old_value
        !          call old_evaluator(value, old_value)
        !     +        old_argument
        !     *        another_old_argument
        !      endif
        !      write(*,*) value
        !      call old_no_argument_path
        !      if( value > 0.0 ) write(*,*) value
        !      parameter(old_threshold=2.0)
        !d     write(7,*) value
        !c     call old_fixed_form_path(value)
        !C old_value = 2.0
        !$     w() = w()+old_disabled_term
        ! C
        ! Original parameters
    end subroutine evaluate
end module example
"""
        converted = remove_dead_commented_code(source)
        self.assertNotIn("old_value", converted)
        self.assertNotIn("old_evaluator", converted)
        self.assertNotIn("another_old_argument", converted)
        self.assertNotIn("old_fixed_form_path", converted)
        self.assertNotIn("old_no_argument_path", converted)
        self.assertNotIn("old_threshold", converted)
        self.assertNotIn("old_value = 2.0", converted)
        self.assertNotIn("old_disabled_term", converted)
        self.assertNotIn("! C", converted)
        self.assertNotIn("write(*,*)", converted)
        self.assertNotIn("Original parameters", converted)

    def test_openmp_sentinel_is_preserved(self) -> None:
        source = """\
subroutine evaluate
    implicit none
    !$omp parallel do
    do index=1, 4
        values(index) = index
    enddo
    !$omp end parallel do
end subroutine evaluate
"""
        converted = remove_dead_commented_code(source)
        self.assertIn("!$omp parallel do", converted)
        self.assertIn("!$omp end parallel do", converted)

    def test_commented_common_block_is_removed_from_header(self) -> None:
        source = """\
! System: H3
! COMMON /PT1CM/ R(3), ENERGY
module example
end module example
"""
        converted = remove_dead_commented_code(source)
        self.assertIn("! System: H3", converted)
        self.assertNotIn("COMMON /PT1CM/", converted)

    def test_scientific_comment_equations_are_preserved(self) -> None:
        source = """\
module example
contains
    subroutine evaluate
        ! R(1) = R(O-H) ! first bond distance
        implicit none
        real :: value
        value = 1.0
        ! d(tanh) = (sech)**2 = 1-(tanh)**2
        ! CPEM(1) = average scaling * CPEM(0)
    end subroutine evaluate
end module example
"""
        converted = remove_dead_commented_code(source)
        self.assertIn("R(1) = R(O-H)", converted)
        self.assertIn("d(tanh) = (sech)**2", converted)
        self.assertIn("CPEM(1) = average scaling", converted)

    def test_pes_documentation_names_surface_and_system(self) -> None:
        source = """\
!> @modern CHEMPOTPY modern Fortran style is enforced for this file.
!> @brief Own the public evaluator for this potential energy surface.
module example
contains
    !> @brief Evaluate this surface through the modern module contract.
    subroutine evaluate
    end subroutine evaluate
end module example
"""
        converted = specialize_pes_documentation(
            source, Path("chempotpy/H3/H3_GEN_A2_2002.f90")
        )
        self.assertNotIn("@modern", converted)
        self.assertIn(
            "!> @brief Implement H3_GEN_A2_2002, a potential energy surface for H3.",
            converted,
        )
        self.assertIn(
            "!> @brief Evaluate the H3_GEN_A2_2002 potential energy surface for H3.",
            converted,
        )

    def test_dpem_documentation_identifies_diabatic_matrix(self) -> None:
        source = """\
!> @brief Own the public evaluator for this potential energy surface.
module example
end module example
"""
        converted = specialize_pes_documentation(
            source, Path("chempotpy/BrH2/BrH2_LEPS_LTBZ2_1995_DPEM.f90")
        )
        self.assertIn(
            "!> @brief Implement BrH2_LEPS_LTBZ2_1995, "
            "a diabatic potential matrix for BrH2.",
            converted,
        )

    def test_short_initializations_use_semicolons(self) -> None:
        source = """\
subroutine initialize
    implicit none
    real :: energy, gradient, coupling

    ! initialize outputs
    energy = 0.0
    gradient = 0.0
    coupling = 0.0
end subroutine initialize
"""
        converted = format_source_text(source)
        self.assertIn(
            "energy = 0.0; gradient = 0.0; coupling = 0.0",
            converted,
        )

    def test_variables_are_separated_by_spaces(self) -> None:
        source = """\
subroutine evaluate( values, count )
    implicit none
    integer, intent(in) :: count
    real, intent(out) :: values(count,3)
    real :: first,second,third
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("values(count,3)", converted)
        self.assertIn("first, second, third", converted)

    def test_array_slices_are_compact(self) -> None:
        source = """\
subroutine evaluate( r, values )
    implicit none
    real, intent(in) :: values(3, 3)
    real, intent(out) :: r(1, 3)
    r(1, :) = values(1, :)
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("values(3,3)", converted)
        self.assertIn("r(1,3)", converted)
        self.assertIn("r(1,:) = values(1,:)", converted)

    def test_pes_metadata_values_are_aligned(self) -> None:
        source = """\
!   System: H3
!   Functional form:
!   Common name: a2
!   Number of electronic surfaces: 1
!   Interface: 3-1V
subroutine evaluate
    implicit none
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("!   System:                        H3", converted)
        self.assertIn("!   Common name:                   a2", converted)
        self.assertIn("!   Number of electronic surfaces: 1", converted)
        self.assertIn("!   Interface:                     3-1V", converted)

    def test_source_units_have_one_blank_separator(self) -> None:
        source = """\
subroutine first
    implicit none
end subroutine first
!> @brief Evaluate second.
subroutine second
    implicit none
end subroutine second
"""
        converted = format_source_text(source)
        self.assertIn("end subroutine first\n\n!> @brief", converted)
        self.assertNotIn("end subroutine first\n\n\n", converted)

    def test_adjacent_contained_procedures_have_one_blank_separator(self) -> None:
        source = """\
module example
    implicit none
    contains
    subroutine first
        implicit none
    end subroutine first
    !> @brief Evaluate second.

    subroutine second
        implicit none
    end subroutine second
end module example
"""
        converted = format_source_text(source)
        self.assertIn(
            "end subroutine first\n\n    !> @brief Evaluate second.\n"
            "    subroutine second",
            converted,
        )
        self.assertNotIn("end subroutine first\n\n\n", converted)
        focused = normalize_procedure_separator_blank_lines(source)
        self.assertIn(
            "end subroutine first\n\n    !> @brief Evaluate second.\n"
            "    subroutine second",
            focused,
        )

    def test_pes_metadata_is_at_the_beginning_of_the_file(self) -> None:
        source = """\
!> @brief Implement example.
module example
    implicit none
    contains
    subroutine evaluate
        implicit none
    end subroutine evaluate
    ! System: CH4F
    ! Functional form:
    ! Common name: ch4f-noso
    ! Notes: Published PES metadata.
end module example
"""
        converted = hoist_pes_metadata(source)
        self.assertTrue(converted.startswith("! System: CH4F\n"))
        self.assertIn(
            "! Notes: Published PES metadata.\n\n!> @brief Implement example.",
            converted,
        )
        self.assertEqual(converted.count("! System: CH4F"), 1)

    def test_contains_has_one_blank_line_before_and_after(self) -> None:
        source = """\
module example
    implicit none

    contains

    subroutine evaluate
        implicit none

    end subroutine evaluate

end module example
"""
        converted = format_source_text(source)
        self.assertIn("    implicit none\n\n    contains\n\n", converted)
        self.assertIn(
            "    contains\n\n    !> @brief Evaluate evaluate.\n",
            converted,
        )
        self.assertNotIn("\n\n\n", converted)
        focused = normalize_contains_blank_lines(source)
        self.assertIn("    implicit none\n\n    contains\n\n", focused)
        self.assertNotIn("\n\n\n", focused)

    def test_call_parentheses_have_no_inner_edge_spaces(self) -> None:
        source = """\
subroutine evaluate( first, second )
    implicit none
    real :: first, second
    call calculate( first, second )
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("subroutine evaluate( first, second )", converted)
        self.assertIn("call calculate(first, second)", converted)

    def test_if_conditions_are_spaced(self) -> None:
        source = """\
subroutine evaluate( r, re )
    implicit none
    real :: r, re
    if((r-re).ne.0.d0)then
        r = re
    endif
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("if( (r-re) .ne. 0.d0 )then", converted)
        self.assertIn("endif", converted)
        self.assertNotIn("endif\n\n", converted)

    def test_if_spacing_preserves_character_literals(self) -> None:
        source = """\
subroutine evaluate( label )
    implicit none
    character(len=*), intent(in) :: label
    if(label=='a>b')then
        return
    endif
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("if( label == 'a>b' )then", converted)

    def test_else_if_uses_compact_keyword(self) -> None:
        source = """\
subroutine evaluate( value )
    implicit none
    integer, intent(in) :: value
    if( value == 1 )then
        return
    else if( value == 2 )then
        return
    endif
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("elseif( value == 2 )then", converted)
        self.assertNotIn("else if", converted)

    def test_multiline_else_if_uses_compact_keyword(self) -> None:
        source = """\
subroutine evaluate( first, second )
    implicit none
    logical, intent(in) :: first, second
    if( first )then
        return
    else if( second .and. &
        first )then
        return
    endif
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("elseif( second .and. &", converted)
        self.assertNotIn("else if", converted)

    def test_compact_end_statements_are_separated(self) -> None:
        source = """\
subroutine evaluate
    implicit none
    if (.true.) then
    endif
endsubroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("endif", converted)
        self.assertIn("end subroutine evaluate", converted)

    def test_labeled_do_continue_is_structured(self) -> None:
        source = """\
subroutine evaluate( values )
    implicit none
    real, intent(out) :: values(2,2)
    integer :: i, j
    do 10 i=1, 2
        do 10 j=1, 2
            values(i,j) = i+j
    10 continue
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("do i=1, 2", converted)
        self.assertIn("do j=1, 2", converted)
        self.assertIn("        10 enddo\n    enddo", converted)
        self.assertNotRegex(converted, r"(?im)^\s*do\s+10\b")
        self.assertNotIn("continue", converted.lower())

    def test_labeled_do_executable_terminal_is_preserved(self) -> None:
        source = """\
subroutine evaluate( values )
    implicit none
    real, intent(out) :: values(2,2)
    integer :: i, j
    do 20 i=1, 2
        do 20 j=1, 2
    20 values(i,j) = i+j
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("20 values(i,j) = i+j\n        enddo\n    enddo", converted)
        self.assertNotRegex(converted, r"(?im)^\s*do\s+20\b")

    def test_only_clauses_are_aligned(self) -> None:
        source = """\
subroutine evaluate
    use geometry, only: distance
    use legacy_conversion_utilities, only: bohr_to_angstrom
    implicit none
end subroutine evaluate
"""
        converted = format_source_text(source)
        self.assertIn("use geometry,                    only: distance", converted)
        self.assertIn(
            "use legacy_conversion_utilities, only: bohr_to_angstrom",
            converted,
        )


if __name__ == "__main__":
    unittest.main()
