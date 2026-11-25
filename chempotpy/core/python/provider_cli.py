"""Command-line entry points for catalog-wide external PES bridges."""

from __future__ import annotations

import argparse

from .catalog import surface_catalog
from .ipi import IPIDriver
from .mdi import MDIEngine, initialize_mdi
from .providers import create_surface
from .socket_transport import SocketPESServer


def _surface(arguments):
    symbols = arguments.atoms.split(",") if arguments.atoms else None
    return create_surface(
        arguments.system,
        arguments.surface,
        atom_symbols=symbols,
        gradient_strategy=arguments.gradient,
    )


def _add_surface_arguments(parser) -> None:
    parser.add_argument("system")
    parser.add_argument("surface")
    parser.add_argument(
        "--atoms",
        help="comma-separated atom symbols (required for variable systems)",
    )
    parser.add_argument(
        "--gradient",
        choices=("auto", "analytic", "numerical"),
        default="auto",
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the bridge command-line parser."""

    parser = argparse.ArgumentParser(prog="chempotpy-provider")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="list all system/surface keys")

    socket_parser = commands.add_parser("socket", help="serve native protocol")
    _add_surface_arguments(socket_parser)
    socket_parser.add_argument("--host", default="127.0.0.1")
    socket_parser.add_argument("--port", type=int, default=8021)

    ipi_parser = commands.add_parser("ipi", help="connect as an i-PI driver")
    _add_surface_arguments(ipi_parser)
    ipi_parser.add_argument("--host", default="127.0.0.1")
    ipi_parser.add_argument("--port", type=int, default=31415)
    ipi_parser.add_argument("--unix")
    ipi_parser.add_argument("--state", type=int, default=0)

    mdi_parser = commands.add_parser("mdi", help="run as an MDI engine")
    _add_surface_arguments(mdi_parser)
    mdi_parser.add_argument("--state", type=int, default=0)
    mdi_parser.add_argument(
        "--mdi-options",
        required=True,
        help="pymdi initialization options, including -role ENGINE",
    )
    return parser


def main(argv=None) -> int:
    """Run the requested catalog or bridge command."""

    arguments = build_parser().parse_args(argv)
    if arguments.command == "list":
        for record in surface_catalog():
            print(record.key)
        return 0
    surface = _surface(arguments)
    if arguments.command == "socket":
        server = SocketPESServer.listen_tcp(
            surface, host=arguments.host, port=arguments.port
        )
        print("CHEMPOTPY socket provider listening on {}:{}".format(*server.address))
        try:
            server.serve_forever()
        finally:
            server.close()
        return 0
    if arguments.command == "ipi":
        if arguments.unix:
            driver = IPIDriver.connect_unix(
                surface, arguments.unix, state_index=arguments.state
            )
        else:
            driver = IPIDriver.connect_tcp(
                surface,
                arguments.host,
                arguments.port,
                state_index=arguments.state,
            )
        try:
            driver.run()
        finally:
            driver.close()
        return 0
    mdi = initialize_mdi(arguments.mdi_options)
    MDIEngine(surface, state_index=arguments.state, mdi_module=mdi).run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
