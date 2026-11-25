"""Tests for persistent socket access to the atomic-unit PES contract."""

import socket
import tempfile
import threading
import unittest
from pathlib import Path

import numpy as np

from chempotpy.core.contract import (
    PESCapabilities,
    PESMetadata,
    PESRequest,
    PESResult,
    Representation,
    validate_coordinates,
)
from chempotpy.core.socket_transport import (
    PESWireProtocolError,
    RemotePESError,
    SocketPESClient,
    SocketPESServer,
    serve_connection,
)


class QuadraticSurface:
    """Provide a deterministic two-state surface for transport tests."""

    metadata = PESMetadata(
        system_name="H2",
        surface_name="quadratic_test",
        atom_symbols=("H", "H"),
        state_count=2,
        capabilities=PESCapabilities(
            gradient=True,
            nonadiabatic_coupling=True,
            diabatic_matrix=True,
            diabatic_gradient=True,
        ),
    )

    def evaluate(
        self, coordinates_bohr: np.ndarray, request: PESRequest
    ) -> PESResult:
        coordinates = validate_coordinates(coordinates_bohr, self.metadata)
        value = np.sum(coordinates * coordinates)
        if request.representation is Representation.DIABATIC:
            matrix = np.array(
                [[value, 0.5 * value], [0.5 * value, 2.0 * value]]
            )
            result = PESResult(diabatic_matrix=matrix)
            if request.gradient:
                factors = np.array([[1.0, 0.5], [0.5, 2.0]])
                result.diabatic_gradient = (
                    factors[:, :, np.newaxis, np.newaxis]
                    * 2.0
                    * coordinates[np.newaxis, np.newaxis, :, :]
                )
        else:
            result = PESResult(energy=np.array([value, 2.0 * value]))
            if request.gradient:
                factors = np.array([1.0, 2.0])
                result.gradient = (
                    factors[:, np.newaxis, np.newaxis]
                    * 2.0
                    * coordinates[np.newaxis, :, :]
                )
            if request.nonadiabatic_coupling:
                coupling = np.zeros((2, 2, 2, 3), dtype=np.float64)
                coupling[0, 1] = coordinates
                coupling[1, 0] = -coordinates
                result.nonadiabatic_coupling = coupling
        result.validate(self.metadata, request)
        return result


class SocketTransportTests( unittest.TestCase ):
    def setUp( self ) -> None:
        client_socket, server_socket = socket.socketpair()
        self.server_thread = threading.Thread(
            target=serve_connection,
            args=(server_socket, QuadraticSurface()),
        )
        self.server_thread.start()
        self.client = SocketPESClient(client_socket)
        self.coordinates = np.array(
            [[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]], dtype=np.float64
        )

    def tearDown( self ) -> None:
        self.client.close()
        self.server_thread.join(timeout=2.0)
        self.assertFalse(self.server_thread.is_alive())

    def test_handshake_returns_surface_metadata( self ) -> None:
        self.assertEqual(self.client.metadata.system_name, "H2")
        self.assertEqual(self.client.metadata.surface_name, "quadratic_test")
        self.assertEqual(self.client.metadata.atom_symbols, ("H", "H"))
        self.assertTrue(self.client.metadata.capabilities.gradient)

    def test_evaluates_adiabatic_energy_and_gradient( self ) -> None:
        request = PESRequest(gradient=True)
        result = self.client.evaluate(self.coordinates, request)
        value = np.sum(self.coordinates * self.coordinates)
        np.testing.assert_array_equal(result.energy, [value, 2.0 * value])
        np.testing.assert_array_equal(
            result.gradient[0], 2.0 * self.coordinates
        )
        np.testing.assert_array_equal(
            result.gradient[1], 4.0 * self.coordinates
        )

    def test_evaluates_diabatic_matrix_and_gradient( self ) -> None:
        request = PESRequest(
            representation=Representation.DIABATIC,
            gradient=True,
        )
        result = self.client.evaluate(self.coordinates, request)
        value = np.sum(self.coordinates * self.coordinates)
        np.testing.assert_array_equal(
            result.diabatic_matrix,
            [[value, 0.5 * value], [0.5 * value, 2.0 * value]],
        )
        np.testing.assert_array_equal(
            result.diabatic_gradient[0, 1], self.coordinates
        )

    def test_evaluates_nonadiabatic_coupling( self ) -> None:
        request = PESRequest(nonadiabatic_coupling=True)
        result = self.client.evaluate(self.coordinates, request)
        np.testing.assert_array_equal(
            result.nonadiabatic_coupling[0, 1], self.coordinates
        )
        np.testing.assert_array_equal(
            result.nonadiabatic_coupling[1, 0], -self.coordinates
        )

    def test_remote_contract_error_does_not_end_session( self ) -> None:
        with self.assertRaises(RemotePESError) as context:
            self.client.evaluate(np.zeros((1, 3)), PESRequest())
        self.assertEqual(context.exception.error_type, "PESContractError")
        result = self.client.evaluate(self.coordinates, PESRequest())
        self.assertEqual(result.energy.shape, (2,))

    def test_closed_client_rejects_evaluation( self ) -> None:
        self.client.close()
        with self.assertRaises(PESWireProtocolError):
            self.client.evaluate(self.coordinates, PESRequest())


class TCPTransportTests( unittest.TestCase ):
    def test_tcp_listener_serves_the_same_contract( self ) -> None:
        server = SocketPESServer.listen_tcp(QuadraticSurface())
        server_thread = threading.Thread(target=server.serve_once)
        server_thread.start()
        host, port = server.address
        try:
            with SocketPESClient.connect_tcp(host, port, timeout=2.0) as client:
                coordinates = np.array(
                    [[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]],
                    dtype=np.float64,
                )
                result = client.evaluate(coordinates, PESRequest())
                np.testing.assert_array_equal(result.energy, [0.5, 1.0])
        finally:
            server.close()
            server_thread.join(timeout=2.0)
        self.assertFalse(server_thread.is_alive())


@unittest.skipUnless(hasattr(socket, "AF_UNIX"), "UNIX sockets unavailable")
class UnixTransportTests( unittest.TestCase ):
    def test_unix_listener_serves_and_removes_its_socket( self ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            socket_path = Path(temporary_directory) / "chempotpy.sock"
            server = SocketPESServer.listen_unix(
                QuadraticSurface(), str(socket_path)
            )
            server_thread = threading.Thread(target=server.serve_once)
            server_thread.start()
            try:
                with SocketPESClient.connect_unix(
                    str(socket_path), timeout=2.0
                ) as client:
                    coordinates = np.array(
                        [[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]],
                        dtype=np.float64,
                    )
                    result = client.evaluate(coordinates, PESRequest())
                    np.testing.assert_array_equal(result.energy, [0.5, 1.0])
            finally:
                server.close()
                server_thread.join(timeout=2.0)
            self.assertFalse(server_thread.is_alive())
            self.assertFalse(socket_path.exists())


if __name__ == "__main__":
    unittest.main()
