"""MolSSI Driver Interface engine bridge for CHEMPOTPY surfaces."""

from __future__ import annotations

import importlib
from typing import Optional

import numpy as np

from .atomic_elements import atomic_numbers
from .contract import PotentialEnergySurface
from .providers import evaluate_state


class MDIBridgeError(RuntimeError):
    """Report missing MDI support or an unsupported MDI command."""


class MDIEngine:
    """Expose one configured PES through standard MDI engine commands."""

    def __init__(
        self,
        surface: PotentialEnergySurface,
        state_index: int = 0,
        mdi_module=None,
    ) -> None:
        self.surface = surface
        self.state_index = state_index
        self.mdi = mdi_module or _import_mdi()
        self.coordinates = np.zeros(
            (surface.metadata.atom_count, 3), dtype=np.float64
        )
        self._evaluation = None
        self._exit_requested = False
        self._atomic_numbers = atomic_numbers(surface.metadata.atom_symbols)
        role = self.mdi.MDI_Get_Role()
        if role != self.mdi.MDI_ENGINE:
            raise MDIBridgeError("MDI role must be ENGINE")
        self._register_commands()
        self.communicator = self.mdi.MDI_Accept_Communicator()

    def run(self) -> None:
        """Receive and execute MDI commands until ``EXIT``."""

        while not self._exit_requested:
            self.execute(self.mdi.MDI_Recv_Command(self.communicator))

    def execute(self, command: str) -> None:
        """Execute one registered command on the active communicator."""

        atom_count = self.surface.metadata.atom_count
        if command == "EXIT":
            self._exit_requested = True
        elif command == "<@":
            self._send("@DEFAULT", self.mdi.MDI_NAME_LENGTH, self.mdi.MDI_CHAR)
        elif command == "<NAME":
            self._send("CHEMPOTPY", self.mdi.MDI_NAME_LENGTH, self.mdi.MDI_CHAR)
        elif command == "<NATOMS":
            self._send(atom_count, 1, self.mdi.MDI_INT)
        elif command == "<LABELS":
            labels = "".join(
                symbol.ljust(self.mdi.MDI_LABEL_LENGTH)
                for symbol in self.surface.metadata.atom_symbols
            )
            self._send(
                labels,
                atom_count * self.mdi.MDI_LABEL_LENGTH,
                self.mdi.MDI_CHAR,
            )
        elif command == "<ELEMENTS" and all(self._atomic_numbers):
            self._send(
                np.asarray(self._atomic_numbers, dtype=np.int32),
                atom_count,
                self.mdi.MDI_INT,
            )
        elif command == "<COORDS":
            self._send(
                np.ascontiguousarray(self.coordinates),
                3 * atom_count,
                self.mdi.MDI_DOUBLE,
            )
        elif command == ">COORDS":
            received = self.mdi.MDI_Recv(
                3 * atom_count, self.mdi.MDI_DOUBLE, self.communicator
            )
            self.coordinates = np.asarray(
                received, dtype=np.float64
            ).reshape(atom_count, 3)
            self._evaluation = None
        elif command == "<ENERGY":
            self._evaluate_if_needed()
            self._send(self._evaluation.energy, 1, self.mdi.MDI_DOUBLE)
        elif command == "<FORCES":
            self._evaluate_if_needed()
            self._send(
                np.ascontiguousarray(-self._evaluation.gradient),
                3 * atom_count,
                self.mdi.MDI_DOUBLE,
            )
        else:
            raise MDIBridgeError(
                "unsupported MDI command {!r}".format(command)
            )

    def _register_commands(self) -> None:
        self.mdi.MDI_Register_Node("@DEFAULT")
        commands = [
            "EXIT", "<@", "<NAME", "<NATOMS", "<LABELS", "<COORDS",
            ">COORDS", "<ENERGY", "<FORCES",
        ]
        if all(self._atomic_numbers):
            commands.append("<ELEMENTS")
        for command in commands:
            self.mdi.MDI_Register_Command("@DEFAULT", command)

    def _evaluate_if_needed(self) -> None:
        if self._evaluation is None:
            self._evaluation = evaluate_state(
                self.surface, self.coordinates, self.state_index
            )

    def _send(self, value, count: int, datatype) -> None:
        self.mdi.MDI_Send(value, count, datatype, self.communicator)


def initialize_mdi(options: str, mdi_module=None):
    """Initialize pymdi with caller-provided engine options."""

    mdi = mdi_module or _import_mdi()
    mdi.MDI_Init(options)
    return mdi


def _import_mdi():
    try:
        return importlib.import_module("mdi")
    except ImportError as error:
        raise MDIBridgeError(
            "MDI support requires the optional pymdi package"
        ) from error
