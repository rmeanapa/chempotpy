"""Persistent socket transport for the atomic-unit CHEMPOTPY PES contract.

The transport is deliberately independent of individual surface modules.  A
server owns one :class:`PotentialEnergySurface` implementation and exposes its
metadata and evaluator through a versioned, length-prefixed protocol.  TCP and
UNIX-domain connection setup are kept outside the scientific contract so the
same client can use either kind of connected socket.
"""

from __future__ import annotations

import base64
import json
import math
import os
import socket
import stat
import struct
from typing import Any, Dict, Mapping, Optional

import numpy as np

from .contract import (
    PESCapabilities,
    PESMetadata,
    PESRequest,
    PESResult,
    PotentialEnergySurface,
    Representation,
)


PROTOCOL_NAME = "chempotpy-pes"
PROTOCOL_VERSION = 1
MAX_MESSAGE_BYTES = 64 * 1024 * 1024
COORDINATE_UNIT = "bohr"
ENERGY_UNIT = "hartree"
GRADIENT_UNIT = "hartree/bohr"
COUPLING_UNIT = "1/bohr"
PROTOCOL_UNITS = {
    "coordinates": COORDINATE_UNIT,
    "energy": ENERGY_UNIT,
    "gradient": GRADIENT_UNIT,
    "nonadiabatic_coupling": COUPLING_UNIT,
}

_MESSAGE_SIZE = struct.Struct("!Q")
_WIRE_FLOAT = np.dtype(">f8")


class PESWireProtocolError( RuntimeError ):
    """Report malformed, incompatible, or unexpectedly terminated traffic."""


class RemotePESError( RuntimeError ):
    """Report an evaluation error returned by a remote PES provider."""

    def __init__( self, error_type: str, message: str ) -> None:
        self.error_type = error_type
        super().__init__("{}: {}".format(error_type, message))


class SocketPESClient:
    """Expose a connected socket as a ``PotentialEnergySurface`` provider."""

    def __init__( self, connection: socket.socket ) -> None:
        self._connection = connection
        self._closed = False
        self._next_request_id = 1
        try:
            _send_message(
                self._connection,
                {
                    "operation": "hello",
                    "protocol": PROTOCOL_NAME,
                    "version": PROTOCOL_VERSION,
                },
            )
            response = _receive_message(self._connection)
            _raise_remote_error(response)
            if response.get("protocol") != PROTOCOL_NAME:
                raise PESWireProtocolError(
                    "remote provider returned the wrong protocol"
                )
            if response.get("version") != PROTOCOL_VERSION:
                raise PESWireProtocolError(
                    "remote provider returned an incompatible version"
                )
            units = _required_mapping(response, "units")
            if dict(units) != PROTOCOL_UNITS:
                raise PESWireProtocolError(
                    "remote provider returned an incompatible unit contract"
                )
            self._metadata = _metadata_from_payload(
                _required_mapping(response, "metadata")
            )
        except Exception:
            self._closed = True
            self._connection.close()
            raise

    @classmethod
    def connect_tcp(
        cls,
        host: str,
        port: int,
        timeout: Optional[float] = None,
    ) -> "SocketPESClient":
        """Connect to a TCP PES provider and complete the protocol handshake."""

        connection = socket.create_connection((host, port), timeout=timeout)
        return cls(connection)

    @classmethod
    def connect_unix(
        cls,
        path: str,
        timeout: Optional[float] = None,
    ) -> "SocketPESClient":
        """Connect to a local UNIX-domain PES provider."""

        if not hasattr(socket, "AF_UNIX"):
            raise PESWireProtocolError("UNIX-domain sockets are unavailable")
        connection = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        if timeout is not None:
            connection.settimeout(timeout)
        connection.connect(path)
        return cls(connection)

    @property
    def metadata( self ) -> PESMetadata:
        """Return metadata received during the provider handshake."""

        return self._metadata

    def evaluate(
        self, coordinates_bohr: np.ndarray, request: PESRequest
    ) -> PESResult:
        """Evaluate one Cartesian geometry through the persistent connection."""

        if self._closed:
            raise PESWireProtocolError("socket PES client is closed")
        request_id = self._next_request_id
        self._next_request_id += 1
        _send_message(
            self._connection,
            {
                "operation": "evaluate",
                "request_id": request_id,
                "coordinates": _array_to_payload(coordinates_bohr),
                "request": _request_to_payload(request),
            },
        )
        response = _receive_message(self._connection)
        if response.get("request_id") != request_id:
            raise PESWireProtocolError("remote provider returned the wrong request ID")
        _raise_remote_error(response)
        result = _result_from_payload(_required_mapping(response, "result"))
        result.validate(self.metadata, request)
        return result

    def close( self ) -> None:
        """Request an orderly shutdown and close the local socket."""

        if self._closed:
            return
        try:
            _send_message(self._connection, {"operation": "close"})
            response = _receive_message(self._connection)
            _raise_remote_error(response)
        finally:
            self._closed = True
            self._connection.close()

    def __enter__( self ) -> "SocketPESClient":
        """Return this connected client for context-manager use."""

        return self

    def __exit__( self, exc_type, exc_value, traceback ) -> None:
        """Close the connection when leaving a context manager."""

        self.close()


class SocketPESServer:
    """Own a listening socket for one configured PES provider."""

    def __init__(
        self,
        listener: socket.socket,
        surface: PotentialEnergySurface,
        unix_path: Optional[str] = None,
    ) -> None:
        self._listener = listener
        self._surface = surface
        self._unix_path = unix_path
        self._closed = False

    @classmethod
    def listen_tcp(
        cls,
        surface: PotentialEnergySurface,
        host: str = "127.0.0.1",
        port: int = 0,
        backlog: int = 16,
    ) -> "SocketPESServer":
        """Listen for PES clients on a TCP address, locally by default."""

        listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            listener.bind((host, port))
            listener.listen(backlog)
        except Exception:
            listener.close()
            raise
        return cls(listener, surface)

    @classmethod
    def listen_unix(
        cls,
        surface: PotentialEnergySurface,
        path: str,
        backlog: int = 16,
    ) -> "SocketPESServer":
        """Listen for local PES clients on a new UNIX-domain socket path."""

        if not hasattr(socket, "AF_UNIX"):
            raise PESWireProtocolError("UNIX-domain sockets are unavailable")
        if os.path.lexists(path):
            raise FileExistsError(
                "refusing to replace existing UNIX socket path {!r}".format(path)
            )
        listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        try:
            listener.bind(path)
            listener.listen(backlog)
        except Exception:
            listener.close()
            if os.path.lexists(path):
                os.unlink(path)
            raise
        return cls(listener, surface, unix_path=path)

    @property
    def address( self ):
        """Return the bound TCP address or UNIX-domain socket path."""

        return self._listener.getsockname()

    def serve_once( self ) -> None:
        """Accept and serve one client synchronously."""

        if self._closed:
            raise PESWireProtocolError("socket PES server is closed")
        connection, _ = self._listener.accept()
        serve_connection(connection, self._surface)

    def serve_forever( self ) -> None:
        """Serve clients serially until another thread closes the listener."""

        while not self._closed:
            try:
                self.serve_once()
            except OSError:
                if self._closed:
                    return
                raise

    def close( self ) -> None:
        """Close the listener and remove its owned UNIX-domain socket path."""

        if self._closed:
            return
        self._closed = True
        self._listener.close()
        if self._unix_path is not None and os.path.lexists(self._unix_path):
            path_status = os.stat(self._unix_path, follow_symlinks=False)
            if stat.S_ISSOCK(path_status.st_mode):
                os.unlink(self._unix_path)

    def __enter__( self ) -> "SocketPESServer":
        """Return this listener for context-manager use."""

        return self

    def __exit__( self, exc_type, exc_value, traceback ) -> None:
        """Close the listener when leaving a context manager."""

        self.close()


def serve_connection(
    connection: socket.socket, surface: PotentialEnergySurface
) -> None:
    """Serve one connected client until it requests closure or disconnects."""

    try:
        hello = _receive_message(connection)
        if hello.get("operation") != "hello":
            raise PESWireProtocolError("the first operation must be hello")
        if hello.get("protocol") != PROTOCOL_NAME:
            raise PESWireProtocolError("the client requested an unknown protocol")
        if hello.get("version") != PROTOCOL_VERSION:
            raise PESWireProtocolError("the client requested an incompatible version")
        _send_message(
            connection,
            {
                "status": "ok",
                "protocol": PROTOCOL_NAME,
                "version": PROTOCOL_VERSION,
                "units": PROTOCOL_UNITS,
                "metadata": _metadata_to_payload(surface.metadata),
            },
        )
        while True:
            try:
                message = _receive_message(connection)
            except EOFError:
                return
            operation = message.get("operation")
            if operation == "close":
                _send_message(connection, {"status": "ok"})
                return
            if operation != "evaluate":
                _send_error(connection, message, PESWireProtocolError(
                    "unknown operation {!r}".format(operation)
                ))
                continue
            try:
                request_id = message.get("request_id")
                coordinates = _array_from_payload(
                    _required_mapping(message, "coordinates")
                )
                request = _request_from_payload(
                    _required_mapping(message, "request")
                )
                result = surface.evaluate(coordinates, request)
                result.validate(surface.metadata, request)
                _send_message(
                    connection,
                    {
                        "status": "ok",
                        "request_id": request_id,
                        "result": _result_to_payload(result),
                    },
                )
            except Exception as error:
                _send_error(connection, message, error)
    finally:
        connection.close()


def _send_error(
    connection: socket.socket,
    message: Mapping[str, Any],
    error: Exception,
) -> None:
    _send_message(
        connection,
        {
            "status": "error",
            "request_id": message.get("request_id"),
            "error_type": type(error).__name__,
            "message": str(error),
        },
    )


def _raise_remote_error( response: Mapping[str, Any] ) -> None:
    status = response.get("status")
    if status == "ok":
        return
    if status == "error":
        raise RemotePESError(
            str(response.get("error_type", "RemoteError")),
            str(response.get("message", "remote provider failed")),
        )
    raise PESWireProtocolError("remote provider returned an invalid status")


def _send_message(
    connection: socket.socket, message: Mapping[str, Any]
) -> None:
    try:
        payload = json.dumps(
            message,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise PESWireProtocolError("message is not JSON serializable") from error
    if len(payload) > MAX_MESSAGE_BYTES:
        raise PESWireProtocolError("message exceeds the protocol size limit")
    connection.sendall(_MESSAGE_SIZE.pack(len(payload)))
    connection.sendall(payload)


def _receive_message( connection: socket.socket ) -> Dict[str, Any]:
    header = _receive_exact(connection, _MESSAGE_SIZE.size)
    payload_size = _MESSAGE_SIZE.unpack(header)[0]
    if payload_size > MAX_MESSAGE_BYTES:
        raise PESWireProtocolError("message exceeds the protocol size limit")
    payload = _receive_exact(connection, payload_size)
    try:
        message = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise PESWireProtocolError("message is not valid UTF-8 JSON") from error
    if not isinstance(message, dict):
        raise PESWireProtocolError("message root must be an object")
    return message


def _receive_exact( connection: socket.socket, size: int ) -> bytes:
    chunks = bytearray()
    while len(chunks) < size:
        chunk = connection.recv(size - len(chunks))
        if not chunk:
            if not chunks:
                raise EOFError("socket connection closed")
            raise PESWireProtocolError("socket closed in the middle of a message")
        chunks.extend(chunk)
    return bytes(chunks)


def _array_to_payload( value: np.ndarray ) -> Dict[str, Any]:
    array = np.asarray(value, dtype=_WIRE_FLOAT, order="C")
    return {
        "dtype": "float64",
        "byte_order": "big",
        "shape": list(array.shape),
        "encoding": "base64",
        "data": base64.b64encode(array.tobytes(order="C")).decode("ascii"),
    }


def _array_from_payload( payload: Mapping[str, Any] ) -> np.ndarray:
    if payload.get("dtype") != "float64":
        raise PESWireProtocolError("array dtype must be float64")
    if payload.get("byte_order") != "big":
        raise PESWireProtocolError("array byte order must be big-endian")
    if payload.get("encoding") != "base64":
        raise PESWireProtocolError("array encoding must be base64")
    shape_value = payload.get("shape")
    if not isinstance(shape_value, list) or any(
        not isinstance(length, int) or length < 0 for length in shape_value
    ):
        raise PESWireProtocolError("array shape is invalid")
    data_value = payload.get("data")
    if not isinstance(data_value, str):
        raise PESWireProtocolError("array data must be a base64 string")
    try:
        raw = base64.b64decode(data_value.encode("ascii"), validate=True)
    except (UnicodeEncodeError, ValueError) as error:
        raise PESWireProtocolError("array data is not valid base64") from error
    expected_size = math.prod(shape_value) * _WIRE_FLOAT.itemsize
    if len(raw) != expected_size:
        raise PESWireProtocolError("array byte count does not match its shape")
    return np.frombuffer(raw, dtype=_WIRE_FLOAT).astype(np.float64).reshape(shape_value)


def _metadata_to_payload( metadata: PESMetadata ) -> Dict[str, Any]:
    return {
        "system_name": metadata.system_name,
        "surface_name": metadata.surface_name,
        "atom_symbols": list(metadata.atom_symbols),
        "state_count": metadata.state_count,
        "capabilities": {
            "gradient": metadata.capabilities.gradient,
            "nonadiabatic_coupling": (
                metadata.capabilities.nonadiabatic_coupling
            ),
            "diabatic_matrix": metadata.capabilities.diabatic_matrix,
            "diabatic_gradient": metadata.capabilities.diabatic_gradient,
        },
    }


def _metadata_from_payload( payload: Mapping[str, Any] ) -> PESMetadata:
    capabilities_payload = _required_mapping(payload, "capabilities")
    atom_symbols = payload.get("atom_symbols")
    state_count = payload.get("state_count")
    if not isinstance(atom_symbols, list) or any(
        not isinstance(symbol, str) for symbol in atom_symbols
    ):
        raise PESWireProtocolError("metadata atom_symbols must be strings")
    if isinstance(state_count, bool) or not isinstance(state_count, int):
        raise PESWireProtocolError("metadata state_count must be an integer")
    return PESMetadata(
        system_name=_required_string(payload, "system_name"),
        surface_name=_required_string(payload, "surface_name"),
        atom_symbols=tuple(atom_symbols),
        state_count=state_count,
        capabilities=PESCapabilities(
            gradient=_optional_bool(capabilities_payload, "gradient"),
            nonadiabatic_coupling=_optional_bool(
                capabilities_payload, "nonadiabatic_coupling"
            ),
            diabatic_matrix=_optional_bool(
                capabilities_payload, "diabatic_matrix"
            ),
            diabatic_gradient=_optional_bool(
                capabilities_payload, "diabatic_gradient"
            ),
        ),
    )


def _request_to_payload( request: PESRequest ) -> Dict[str, Any]:
    return {
        "representation": request.representation.value,
        "gradient": request.gradient,
        "nonadiabatic_coupling": request.nonadiabatic_coupling,
    }


def _request_from_payload( payload: Mapping[str, Any] ) -> PESRequest:
    try:
        representation = Representation(payload.get("representation"))
    except ValueError as error:
        raise PESWireProtocolError("unknown electronic representation") from error
    return PESRequest(
        representation=representation,
        gradient=_optional_bool(payload, "gradient"),
        nonadiabatic_coupling=_optional_bool(
            payload, "nonadiabatic_coupling"
        ),
    )


def _result_to_payload( result: PESResult ) -> Dict[str, Any]:
    return {
        "energy": _optional_array_to_payload(result.energy),
        "gradient": _optional_array_to_payload(result.gradient),
        "nonadiabatic_coupling": _optional_array_to_payload(
            result.nonadiabatic_coupling
        ),
        "diabatic_matrix": _optional_array_to_payload(result.diabatic_matrix),
        "diabatic_gradient": _optional_array_to_payload(
            result.diabatic_gradient
        ),
    }


def _result_from_payload( payload: Mapping[str, Any] ) -> PESResult:
    return PESResult(
        energy=_optional_array_from_payload(payload.get("energy")),
        gradient=_optional_array_from_payload(payload.get("gradient")),
        nonadiabatic_coupling=_optional_array_from_payload(
            payload.get("nonadiabatic_coupling")
        ),
        diabatic_matrix=_optional_array_from_payload(
            payload.get("diabatic_matrix")
        ),
        diabatic_gradient=_optional_array_from_payload(
            payload.get("diabatic_gradient")
        ),
    )


def _optional_array_to_payload(
    value: Optional[np.ndarray],
) -> Optional[Dict[str, Any]]:
    if value is None:
        return None
    return _array_to_payload(value)


def _optional_array_from_payload( value: Any ) -> Optional[np.ndarray]:
    if value is None:
        return None
    if not isinstance(value, dict):
        raise PESWireProtocolError("array payload must be an object or null")
    return _array_from_payload(value)


def _required_mapping(
    payload: Mapping[str, Any], name: str
) -> Mapping[str, Any]:
    value = payload.get(name)
    if not isinstance(value, dict):
        raise PESWireProtocolError("{} must be an object".format(name))
    return value


def _required_string( payload: Mapping[str, Any], name: str ) -> str:
    value = payload.get(name)
    if not isinstance(value, str):
        raise PESWireProtocolError("{} must be a string".format(name))
    return value


def _optional_bool(
    payload: Mapping[str, Any], name: str, default: bool = False
) -> bool:
    value = payload.get(name, default)
    if not isinstance(value, bool):
        raise PESWireProtocolError("{} must be a boolean".format(name))
    return value
