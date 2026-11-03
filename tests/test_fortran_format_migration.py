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
        self.assertIn("if (.not. valid) return", converted)


if __name__ == "__main__":
    unittest.main()
