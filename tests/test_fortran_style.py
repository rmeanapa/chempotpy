"""Tests for the opt-in modern Fortran style gate."""

import unittest

from scripts.check_fortran_style import check_source, modern_sources


class FortranStyleTests( unittest.TestCase ):
    def test_all_marked_sources_satisfy_modern_style( self ) -> None:
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
