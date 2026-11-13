"""Unit tests for the immutable-baseline numerical regression checker."""

import unittest

from scripts.capture_library_references import coordinates
from scripts.check_numerical_regression import (
    compare_records,
    escape_workflow_command,
    select_baseline,
)


def record(values):
    """Return one minimal PES reference record."""

    return {
        "entry_point": "pes",
        "natoms": 3,
        "nstates": 1,
        "signature": "x,igrad,p,g,d",
        "values": values,
    }


class NumericalRegressionTests(unittest.TestCase):
    """Verify strict numerical and interface comparisons."""

    def test_equal_records_pass(self) -> None:
        diagnostics, count = compare_records(
            "H3/example", record(["1.0", "-2.0"]),
            record(["1.0", "-2.0"])
        )
        self.assertEqual(diagnostics, [])
        self.assertEqual(count, 2)

    def test_roundoff_within_cross_platform_tolerance_passes(self) -> None:
        diagnostics, _ = compare_records(
            "H3/example", record(["1.00000000000001"]), record(["1.0"])
        )
        self.assertEqual(diagnostics, [])

    def test_numerical_change_fails_with_surface_and_index(self) -> None:
        diagnostics, _ = compare_records(
            "H3/example", record(["1.01"]), record(["1.0"])
        )
        self.assertEqual(len(diagnostics), 1)
        self.assertIn("H3/example: value 0 changed", diagnostics[0])

    def test_matching_nan_still_fails(self) -> None:
        diagnostics, _ = compare_records(
            "H3/example", record(["NaN"]), record(["NaN"])
        )
        self.assertEqual(len(diagnostics), 1)
        self.assertIn("is non-finite", diagnostics[0])

    def test_new_non_finite_value_fails(self) -> None:
        diagnostics, _ = compare_records(
            "H3/example", record(["NaN"]), record(["1.0"])
        )
        self.assertEqual(len(diagnostics), 1)
        self.assertIn("is non-finite", diagnostics[0])

    def test_interface_change_fails(self) -> None:
        actual = record(["1.0"])
        expected = record(["1.0"])
        actual["natoms"] = 4
        diagnostics, _ = compare_records("H3/example", actual, expected)
        self.assertEqual(len(diagnostics), 1)
        self.assertIn("natoms changed", diagnostics[0])

    def test_workflow_command_text_is_escaped(self) -> None:
        self.assertEqual(
            escape_workflow_command("100%\r\nfailed"),
            "100%25%0D%0Afailed",
        )

    def test_legacy_fixture_bootstraps_version_one_once(self) -> None:
        current = {"H3/example": record(["1.0"])}
        selected, bootstrapped = select_baseline(
            None, {"H3/example": record(["0.0"])}, 1, current
        )
        self.assertIs(selected, current)
        self.assertTrue(bootstrapped)

    def test_versioned_base_cannot_be_replaced_by_current_fixture(self) -> None:
        baseline = {"H3/example": record(["1.0"])}
        current = {"H3/example": record(["2.0"])}
        selected, bootstrapped = select_baseline(
            1, baseline, 1, current
        )
        self.assertIs(selected, baseline)
        self.assertFalse(bootstrapped)

    def test_problem_systems_use_physical_regression_geometries(self) -> None:
        for system_name, atom_count in (("Aln", 5), ("NH3", 4)):
            for geometry_index in (1, 2):
                points = coordinates(
                    atom_count,
                    geometry_index,
                    system_name=system_name,
                )
                self.assertEqual(len(points), atom_count)
                self.assertEqual(len(set(points)), atom_count)


if __name__ == "__main__":
    unittest.main()
