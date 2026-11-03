"""Tests for extracting repeated compatibility-wrapper operations."""

import unittest

from scripts.extract_wrapper_core_utilities import migrate_source


class WrapperCoreMigrationTests(unittest.TestCase):
    def test_triatomic_distance_block_uses_one_core_adapter(self) -> None:
        source = """\
subroutine pes( x, energy )
    use chempotpy_geometry, only: cartesian_distance
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, legacy_hartree_to_electronvolt
    implicit none
    double precision, intent(in) :: x(3, 3)
    double precision, intent(out) :: energy
    double precision :: r(1, 3)
    r(1, 1) = cartesian_distance(x(1, :), x(2, :))/legacy_bohr_to_angstrom
    r(1, 2) = cartesian_distance(x(2, :), x(3, :))/legacy_bohr_to_angstrom
    r(1, 3) = cartesian_distance(x(1, :), x(3, :))/legacy_bohr_to_angstrom
    energy = legacy_hartree_to_electronvolt
end subroutine pes
"""
        migrated, distances, coordinates, energies, triatomic = (
            migrate_source(source)
        )
        self.assertEqual((distances, coordinates, energies, triatomic), (0, 0, 0, 1))
        self.assertIn(
            "call legacy_triatomic_distances_to_bohr( x, r(1, :) )",
            migrated,
        )
        self.assertIn(
            "use chempotpy_legacy_coordinates, only: "
            "legacy_triatomic_distances_to_bohr",
            migrated,
        )
        self.assertNotIn("use chempotpy_geometry", migrated)
        self.assertNotIn("legacy_bohr_to_angstrom,", migrated)
        self.assertIn("legacy_hartree_to_electronvolt", migrated)


if __name__ == "__main__":
    unittest.main()
