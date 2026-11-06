"""Tests for tolerance-guarded reference hash refreshes."""

import unittest

from scripts.refresh_library_reference_hashes import equivalent_values


class RefreshLibraryReferenceHashesTests(unittest.TestCase):
    """Verify that hash refreshes accept roundoff, not numerical changes."""

    def test_identical_values_are_equivalent(self) -> None:
        self.assertTrue(equivalent_values(["1.0"], ["1.0"]))

    def test_roundoff_level_values_are_equivalent(self) -> None:
        self.assertTrue(equivalent_values(["1.0"], ["1.00000000000001"]))

    def test_materially_changed_values_are_rejected(self) -> None:
        self.assertFalse(equivalent_values(["1.0"], ["1.000001"]))

    def test_different_value_counts_are_rejected(self) -> None:
        self.assertFalse(equivalent_values(["1.0"], ["1.0", "2.0"]))


if __name__ == "__main__":
    unittest.main()
