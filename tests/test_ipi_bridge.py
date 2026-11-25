"""End-to-end protocol tests for the i-PI catalog bridge."""

import json
import socket
import struct
import threading
import unittest

import numpy as np

from chempotpy.core.contract import (
    PESCapabilities,
    PESMetadata,
    PESRequest,
    PESResult,
    validate_coordinates,
)
from chempotpy.core.ipi import IPIDriver


INT32 = struct.Struct("=i")
FLOAT64 = struct.Struct("=d")


class QuadraticSurface:
    metadata = PESMetadata(
        system_name="H2",
        surface_name="quadratic",
        atom_symbols=("H", "H"),
        capabilities=PESCapabilities(gradient=True),
    )

    def evaluate(self, coordinates_bohr, request: PESRequest) -> PESResult:
        coordinates = validate_coordinates(coordinates_bohr, self.metadata)
        return PESResult(
            energy=np.array([np.sum(coordinates * coordinates)]),
            gradient=np.array([2.0 * coordinates]),
        )


def send_header(connection, command):
    connection.sendall(command.encode("ascii").ljust(12, b" "))


def receive_exact(connection, size):
    payload = bytearray()
    while len(payload) < size:
        payload.extend(connection.recv(size - len(payload)))
    return bytes(payload)


class IPIBridgeTests(unittest.TestCase):
    def test_protocol_round_trip_returns_atomic_energy_and_force(self) -> None:
        server_socket, driver_socket = socket.socketpair()
        driver = IPIDriver(driver_socket, QuadraticSurface())
        thread = threading.Thread(target=driver.run)
        thread.start()
        coordinates = np.array(
            [[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]], dtype="=f8"
        )
        try:
            send_header(server_socket, "STATUS")
            self.assertEqual(receive_exact(server_socket, 12).strip(), b"READY")
            send_header(server_socket, "POSDATA")
            server_socket.sendall(np.eye(3, dtype="=f8").tobytes())
            server_socket.sendall(np.eye(3, dtype="=f8").tobytes())
            server_socket.sendall(INT32.pack(2))
            server_socket.sendall(coordinates.tobytes())
            send_header(server_socket, "STATUS")
            self.assertEqual(receive_exact(server_socket, 12).strip(), b"HAVEDATA")
            send_header(server_socket, "GETFORCE")
            self.assertEqual(
                receive_exact(server_socket, 12).strip(), b"FORCEREADY"
            )
            energy = FLOAT64.unpack(receive_exact(server_socket, 8))[0]
            atom_count = INT32.unpack(receive_exact(server_socket, 4))[0]
            force = np.frombuffer(
                receive_exact(server_socket, 3 * atom_count * 8), dtype="=f8"
            ).reshape(atom_count, 3)
            receive_exact(server_socket, 9 * 8)
            extras_size = INT32.unpack(receive_exact(server_socket, 4))[0]
            extras = json.loads(receive_exact(server_socket, extras_size))
            self.assertAlmostEqual(energy, 0.5)
            np.testing.assert_array_equal(force, -2.0 * coordinates)
            self.assertEqual(extras["surface"], "quadratic")
            send_header(server_socket, "EXIT")
        finally:
            server_socket.close()
            thread.join(timeout=2.0)
        self.assertFalse(thread.is_alive())


if __name__ == "__main__":
    unittest.main()
