"""Tests for Python-level compiled extension validation."""

import unittest

from scripts.check_python_extensions import values_match


class CheckPythonExtensionsTests( unittest.TestCase ):
    def test_expected_nan_matches_nan_only( self ) -> None:
        self.assertTrue(values_match("NaN", float("nan")))
        self.assertFalse(values_match("NaN", 0.0))


if __name__ == "__main__":
    unittest.main()
