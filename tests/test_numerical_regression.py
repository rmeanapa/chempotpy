"""Unit tests for the immutable-baseline numerical regression checker."""

import unittest

from scripts.capture_library_references import (
    BASELINE_VERSION,
    REGRESSION_GEOMETRY_COUNT,
    coordinates,
    standard_driver,
)
from scripts.check_numerical_regression import (
    compare_records,
    comparison_record,
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
            "H3/example", record(["1.0000000004"]), record(["1.0"])
        )
        self.assertEqual(diagnostics, [])

    def test_change_outside_cross_platform_tolerance_fails(self) -> None:
        diagnostics, _ = compare_records(
            "H3/example", record(["1.00000005"]), record(["1.0"])
        )
        self.assertEqual(len(diagnostics), 1)

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

    def test_new_surface_bootstraps_values_but_keeps_fixture_metadata(self) -> None:
        actual = record(["2.0"])
        fixture = record(["1.0"])
        fixture["natoms"] = 4
        selected = comparison_record(
            "H4/new", actual, {}, {"H4/new": fixture}
        )
        self.assertEqual(selected["values"], ["2.0"])
        self.assertEqual(selected["natoms"], 4)

    def test_existing_surface_uses_base_revision_values(self) -> None:
        baseline = {"H3/example": record(["1.0"])}
        selected = comparison_record(
            "H3/example",
            record(["2.0"]),
            baseline,
            {"H3/example": record(["3.0"])},
        )
        self.assertIs(selected, baseline["H3/example"])

    def test_legacy_fixture_bootstraps_current_version_once(self) -> None:
        current = {"H3/example": record(["1.0"])}
        selected, bootstrapped = select_baseline(
            None, {"H3/example": record(["0.0"])}, BASELINE_VERSION, current
        )
        self.assertIs(selected, current)
        self.assertTrue(bootstrapped)

    def test_previous_version_uses_approved_changed_record(self) -> None:
        previous = {"H3/example": record(["1.0"])}
        current = {"H3/example": record(["2.0"])}
        selected, bootstrapped = select_baseline(
            BASELINE_VERSION-1, previous, BASELINE_VERSION, current
        )
        self.assertIs(selected["H3/example"], current["H3/example"])
        self.assertTrue(bootstrapped)

    def test_previous_version_keeps_same_platform_unchanged_record(self) -> None:
        previous = {"H3/example": record(["1.0"])}
        current = {"H3/example": record(["1.0"])}
        selected, bootstrapped = select_baseline(
            BASELINE_VERSION-1, previous, BASELINE_VERSION, current
        )
        self.assertIs(selected["H3/example"], previous["H3/example"])
        self.assertTrue(bootstrapped)

    def test_versioned_base_cannot_be_replaced_by_current_fixture(self) -> None:
        baseline = {"H3/example": record(["1.0"])}
        current = {"H3/example": record(["2.0"])}
        selected, bootstrapped = select_baseline(
            BASELINE_VERSION, baseline, BASELINE_VERSION, current
        )
        self.assertIs(selected, baseline)
        self.assertFalse(bootstrapped)

    def test_problem_systems_use_physical_regression_geometries(self) -> None:
        for system_name, atom_count in (("Aln", 5), ("NH3", 4)):
            geometries = []
            for geometry_index in range(1, REGRESSION_GEOMETRY_COUNT+1):
                points = coordinates(
                    atom_count,
                    geometry_index,
                    system_name=system_name,
                )
                self.assertEqual(len(points), atom_count)
                self.assertEqual(len(set(points)), atom_count)
                geometries.append(points)
            self.assertEqual(len(set(geometries)), REGRESSION_GEOMETRY_COUNT)

    def test_each_endpoint_pair_has_its_interpolated_midpoint(self) -> None:
        geometries = [
            coordinates(3, geometry_index)
            for geometry_index in range(1, REGRESSION_GEOMETRY_COUNT+1)
        ]
        self.assertEqual(len(set(geometries)), REGRESSION_GEOMETRY_COUNT)
        for first_index, second_index, midpoint_index in (
            (0, 1, 2),
            (3, 4, 5),
            (6, 7, 8),
        ):
            first = geometries[first_index]
            second = geometries[second_index]
            midpoint = geometries[midpoint_index]
            for atom in range(3):
                for axis in range(3):
                    self.assertAlmostEqual(
                        midpoint[atom][axis],
                        0.5*first[atom][axis]+0.5*second[atom][axis],
                    )

    def test_additional_pairs_change_internal_distances(self) -> None:
        geometries = [
            coordinates(4, geometry_index)
            for geometry_index in range(1, REGRESSION_GEOMETRY_COUNT+1)
        ]
        pair_distances = []
        for geometry in geometries:
            distances = []
            for first_atom in range(3):
                for second_atom in range(first_atom+1, 4):
                    distances.append(sum(
                        (geometry[first_atom][axis]
                         - geometry[second_atom][axis])**2
                        for axis in range(3)
                    )**0.5)
            pair_distances.append(tuple(distances))
        self.assertEqual(
            len(set(pair_distances)),
            REGRESSION_GEOMETRY_COUNT,
        )

    def test_variable_atom_driver_uses_the_public_array_shape(self) -> None:
        driver = standard_driver(
            "pes",
            "xyz,nsi,no,nh,natoms,igrad,p,g,d",
            3,
            1,
            None,
            "SiOH",
        )

        self.assertIn("real(8) :: xyz(3,3)", driver)
        self.assertNotIn("real(8) :: xyz(10000,3)", driver)


if __name__ == "__main__":
    unittest.main()
