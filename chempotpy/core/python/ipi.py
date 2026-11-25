"""i-PI driver protocol bridge for any configured CHEMPOTPY surface."""

from __future__ import annotations

import json
import socket
import struct
from typing import Optional

import numpy as np

from .contract import PotentialEnergySurface
from .providers import evaluate_state


_HEADER_BYTES = 12
_INT32 = struct.Struct("=i")
_FLOAT64 = struct.Struct("=d")


class IPIProtocolError(RuntimeError):
    """Report malformed or incompatible i-PI driver traffic."""


class IPIDriver:
    """Serve one PES as a state-selectable i-PI force provider."""

    def __init__(
        self,
        connection: socket.socket,
        surface: PotentialEnergySurface,
        state_index: int = 0,
    ) -> None:
        self._connection = connection
        self._surface = surface
        self._state_index = state_index
        self._evaluation = None
        self._closed = False

    @classmethod
    def connect_tcp(
        cls,
        surface: PotentialEnergySurface,
        host: str,
        port: int,
        state_index: int = 0,
        timeout: Optional[float] = None,
    ) -> "IPIDriver":
        """Connect to an i-PI TCP server and return a configured driver."""

        connection = socket.create_connection((host, port), timeout=timeout)
        return cls(connection, surface, state_index=state_index)

    @classmethod
    def connect_unix(
        cls,
        surface: PotentialEnergySurface,
        path: str,
        state_index: int = 0,
        timeout: Optional[float] = None,
    ) -> "IPIDriver":
        """Connect to an i-PI UNIX-domain server."""

        if not hasattr(socket, "AF_UNIX"):
            raise IPIProtocolError("UNIX-domain sockets are unavailable")
        connection = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        if timeout is not None:
            connection.settimeout(timeout)
        connection.connect(path)
        return cls(connection, surface, state_index=state_index)

    def run(self) -> None:
        """Process i-PI commands until the server sends ``EXIT``."""

        while not self._closed:
            command = self._receive_header()
            if command == "STATUS":
                self._send_header("HAVEDATA" if self._evaluation else "READY")
            elif command == "INIT":
                self._receive_initialization()
            elif command == "POSDATA":
                self._receive_positions_and_evaluate()
            elif command == "GETFORCE":
                self._send_force()
            elif command == "EXIT":
                self.close()
            else:
                raise IPIProtocolError(
                    "unsupported i-PI command {!r}".format(command)
                )

    def close(self) -> None:
        """Close the driver socket idempotently."""

        if not self._closed:
            self._closed = True
            self._connection.close()

    def _receive_initialization(self) -> None:
        self._receive_exact(_INT32.size)
        message_length = _INT32.unpack(self._receive_exact(_INT32.size))[0]
        if message_length < 0:
            raise IPIProtocolError("negative i-PI initialization length")
        self._receive_exact(message_length)
        self._evaluation = None

    def _receive_positions_and_evaluate(self) -> None:
        self._receive_exact(18 * _FLOAT64.size)
        atom_count = _INT32.unpack(self._receive_exact(_INT32.size))[0]
        expected = self._surface.metadata.atom_count
        if atom_count != expected:
            raise IPIProtocolError(
                "i-PI sent {} atoms; surface requires {}".format(
                    atom_count, expected
                )
            )
        coordinate_bytes = self._receive_exact(
            3 * atom_count * _FLOAT64.size
        )
        coordinates = np.frombuffer(
            coordinate_bytes, dtype=np.dtype("=f8")
        ).copy().reshape(atom_count, 3)
        self._evaluation = evaluate_state(
            self._surface, coordinates, self._state_index
        )

    def _send_force(self) -> None:
        if self._evaluation is None:
            raise IPIProtocolError("GETFORCE received before POSDATA")
        atom_count = self._surface.metadata.atom_count
        force = np.ascontiguousarray(-self._evaluation.gradient, dtype="=f8")
        extras = json.dumps(
            {
                "system": self._surface.metadata.system_name,
                "surface": self._surface.metadata.surface_name,
                "state_index": self._state_index,
            },
            separators=(",", ":"),
        ).encode("utf-8")
        self._send_header("FORCEREADY")
        self._send_all(_FLOAT64.pack(self._evaluation.energy))
        self._send_all(_INT32.pack(atom_count))
        self._send_all(force.tobytes())
        self._send_all(np.zeros(9, dtype="=f8").tobytes())
        self._send_all(_INT32.pack(len(extras)))
        self._send_all(extras)
        self._evaluation = None

    def _receive_header(self) -> str:
        payload = self._receive_exact(_HEADER_BYTES)
        try:
            return payload.decode("ascii").strip().upper()
        except UnicodeDecodeError as error:
            raise IPIProtocolError("non-ASCII i-PI command") from error

    def _send_header(self, command: str) -> None:
        encoded = command.encode("ascii")
        if len(encoded) > _HEADER_BYTES:
            raise IPIProtocolError("i-PI command exceeds 12 bytes")
        self._send_all(encoded.ljust(_HEADER_BYTES, b" "))

    def _receive_exact(self, size: int) -> bytes:
        chunks = bytearray()
        while len(chunks) < size:
            chunk = self._connection.recv(size - len(chunks))
            if not chunk:
                raise IPIProtocolError("i-PI connection ended unexpectedly")
            chunks.extend(chunk)
        return bytes(chunks)

    def _send_all(self, payload: bytes) -> None:
        self._connection.sendall(payload)
