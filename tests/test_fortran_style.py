"""Tests for the repository-wide modern Fortran style gate."""

import re
import unittest

from scripts.check_fortran_style import check_source, modern_sources


class FortranStyleTests( unittest.TestCase ):
    def test_entry_statements_are_forbidden( self ) -> None:
        offenders = [
            "{}:{}".format(source, line_number)
            for source in modern_sources()
            for line_number, line in enumerate(
                source.read_text(errors="replace").splitlines(),
                start=1,
            )
            if re.match(
                r"^\s*entry\b",
                line.split("!", 1)[0],
                flags=re.IGNORECASE,
            )
        ]
        self.assertEqual(offenders, [])

    def test_all_sources_satisfy_modern_style( self ) -> None:
        sources = modern_sources()
        self.assertGreaterEqual(len(sources), 2)
        violations = [
            violation
            for source in sources
            for violation in check_source(source)
        ]
        self.assertEqual(violations, [])


if __name__ == "__main__":
    unittest.main()
