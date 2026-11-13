"""Tests for the shared PES shapes, capabilities, and units."""

import unittest

import numpy as np

from chempotpy.core.contract import (
    PESCapabilities,
    PESContractError,
    PESMetadata,
    PESRequest,
    PESResult,
    Representation,
    validate_coordinates,
)
from chempotpy.core.units import (
    BOHR_TO_ANGSTROM,
    ELECTRONVOLT_TO_INVERSE_CENTIMETER,
    ELECTRONVOLT_TO_KILOCALORIE_PER_MOLE,
    HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM,
    HARTREE_TO_EV,
    HARTREE_TO_INVERSE_CENTIMETER,
    HARTREE_TO_KILOCALORIE_PER_MOLE,
    angstrom_to_bohr,
    bohr_to_angstrom,
    ev_to_hartree,
    hartree_to_ev,
)


class ContractTests( unittest.TestCase ):
    def setUp( self ) -> None:
        self.metadata = PESMetadata(
            system_name="H3",
            surface_name="example",
            atom_symbols=("H", "H", "H"),
            state_count=2,
            capabilities=PESCapabilities(
                gradient=True,
                nonadiabatic_coupling=True,
                diabatic_matrix=True,
                diabatic_gradient=True,
            ),
        )

    def test_adiabatic_result_shapes( self ) -> None:
        request = PESRequest(gradient=True, nonadiabatic_coupling=True)
        result = PESResult(
            energy=np.zeros(2),
            gradient=np.zeros((2, 3, 3)),
            nonadiabatic_coupling=np.zeros((2, 2, 3, 3)),
        )
        result.validate(self.metadata, request)

    def test_diabatic_result_shapes( self ) -> None:
        request = PESRequest(
            representation=Representation.DIABATIC, gradient=True
        )
        result = PESResult(
            diabatic_matrix=np.zeros((2, 2)),
            diabatic_gradient=np.zeros((2, 2, 3, 3)),
        )
        result.validate(self.metadata, request)

    def test_rejects_wrong_coordinate_shape( self ) -> None:
        with self.assertRaises(PESContractError):
            validate_coordinates(np.zeros((3, 2)), self.metadata)

    def test_unit_round_trips( self ) -> None:
        coordinates = np.array([[0.1, 0.2, 0.3]])
        energy = np.array([1.25])
        np.testing.assert_allclose(
            bohr_to_angstrom(angstrom_to_bohr(coordinates)), coordinates
        )
        np.testing.assert_allclose(
            ev_to_hartree(hartree_to_ev(energy)), energy
        )

    def test_codata_2022_conversion_values( self ) -> None:
        self.assertEqual(BOHR_TO_ANGSTROM, 0.529177210544)
        self.assertEqual(HARTREE_TO_EV, 27.211386245981)
        self.assertAlmostEqual(
            HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM,
            51.42206751119799,
            places=14,
        )
        self.assertAlmostEqual(
            HARTREE_TO_KILOCALORIE_PER_MOLE,
            627.5094740628975,
            places=12,
        )
        self.assertAlmostEqual(
            ELECTRONVOLT_TO_KILOCALORIE_PER_MOLE,
            23.06054783061903,
            places=14,
        )
        self.assertAlmostEqual(
            HARTREE_TO_INVERSE_CENTIMETER,
            219474.63136314,
            places=8,
        )
        self.assertAlmostEqual(
            ELECTRONVOLT_TO_INVERSE_CENTIMETER,
            8065.543937349212,
            places=11,
        )


if __name__ == "__main__":
    unittest.main()
