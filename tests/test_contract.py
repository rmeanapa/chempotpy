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


if __name__ == "__main__":
    unittest.main()
