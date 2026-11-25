"""Tests for the MDI command adapter without requiring pymdi."""

import unittest

import numpy as np

from chempotpy.core.contract import (
    PESCapabilities,
    PESMetadata,
    PESRequest,
    PESResult,
)
from chempotpy.core.mdi import MDIEngine


class QuadraticSurface:
    metadata = PESMetadata(
        system_name="H2",
        surface_name="quadratic",
        atom_symbols=("H", "H"),
        capabilities=PESCapabilities(gradient=True),
    )

    def evaluate(self, coordinates_bohr, request: PESRequest) -> PESResult:
        value = np.sum(coordinates_bohr * coordinates_bohr)
        return PESResult(
            energy=np.array([value]),
            gradient=np.array([2.0 * coordinates_bohr]),
        )


class FakeMDI:
    MDI_ENGINE = 1
    MDI_INT = 2
    MDI_DOUBLE = 3
    MDI_CHAR = 4
    MDI_NAME_LENGTH = 12
    MDI_LABEL_LENGTH = 4

    def __init__(self):
        self.commands = []
        self.sent = []
        self.received = np.array([0.0, 0.0, -0.5, 0.0, 0.0, 0.5])

    def MDI_Get_Role(self):
        return self.MDI_ENGINE

    def MDI_Register_Node(self, node):
        self.node = node

    def MDI_Register_Command(self, node, command):
        self.commands.append(command)

    def MDI_Accept_Communicator(self):
        return 7

    def MDI_Recv(self, count, datatype, communicator):
        return self.received.copy()

    def MDI_Send(self, value, count, datatype, communicator):
        self.sent.append((value, count, datatype, communicator))


class MDIBridgeTests(unittest.TestCase):
    def test_registers_and_executes_standard_energy_force_commands(self) -> None:
        mdi = FakeMDI()
        engine = MDIEngine(QuadraticSurface(), mdi_module=mdi)
        self.assertIn("<ELEMENTS", mdi.commands)
        engine.execute(">COORDS")
        engine.execute("<ENERGY")
        engine.execute("<FORCES")
        self.assertAlmostEqual(mdi.sent[0][0], 0.5)
        np.testing.assert_array_equal(
            mdi.sent[1][0],
            np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]),
        )


if __name__ == "__main__":
    unittest.main()
