"""Tests for labeled-format and redundant-return modernization."""

import unittest

from scripts.format_fortran import format_source_text
from scripts.inline_fortran_formats import convert_source


class FortranFormatMigrationTests(unittest.TestCase):
    def test_format_is_moved_into_write_statement(self) -> None:
        source = """\
subroutine report( value )
    implicit none
    double precision :: value
    write(6, 100) value
100 format('value = ', f12.6)
    return
end subroutine report
"""
        converted = format_source_text(convert_source(source).text)
        self.assertIn("write(6, '(''value = '', f12.6)') value", converted)
        self.assertNotRegex(converted.lower(), r"\bformat\s*\(")
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

    def test_modules_have_no_internal_blank_lines(self) -> None:
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
        self.assertNotIn("\n\n", converted)

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
