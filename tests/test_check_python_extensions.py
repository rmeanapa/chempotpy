"""Tests for Python-level compiled extension validation."""

import unittest
from pathlib import Path

import numpy as np

from scripts.check_python_extensions import call_surface, values_match


class CheckPythonExtensionsTests( unittest.TestCase ):
    def test_expected_nan_matches_nan_only( self ) -> None:
        self.assertTrue(values_match("NaN", float("nan")))
        self.assertFalse(values_match("NaN", 0.0))

    def test_cross_toolchain_roundoff_uses_regression_tolerance( self ) -> None:
        self.assertTrue(values_match("6.60457730896920303", 6.60457730896168549))
        self.assertFalse(values_match("6.60457730896920303", 6.6045772))

    def test_sgdml_uses_its_reviewed_reduction_tolerance( self ) -> None:
        self.assertTrue(values_match(
            "5.77890957375201797",
            5.77890956989007432,
            "H3O2m/H3O2m_sGDML_2024",
        ))
        self.assertFalse(values_match(
            "5.77890957375201797",
            5.77890956989007432,
            "H3/example",
        ))

    def test_variable_atom_wrapper_places_optional_count_after_igrad(
        self,
    ) -> None:
        class FakeModule:
            arguments = None

            @classmethod
            def pes(cls, *arguments):
                cls.arguments = arguments
                return (np.asarray([0.0]),)

        record = {
            "entry_point": "pes",
            "signature": "xyz,nsi,no,nh,natoms,igrad,p,g,d",
            "natoms": 3,
            "nstates": 1,
        }

        values = call_surface(FakeModule, record, "SiOH", 1, Path("."))

        self.assertEqual(values, (0.0,))
        self.assertEqual(FakeModule.arguments[0].shape, (3, 3))
        self.assertEqual(FakeModule.arguments[4:], (0, 3))


if __name__ == "__main__":
    unittest.main()
