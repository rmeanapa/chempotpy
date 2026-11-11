#!/usr/bin/env python3
"""Generate immutable Fortran parameters from released neural-network files."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
from pathlib import Path
import sys
from typing import List, Optional, Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.format_fortran import format_source_text


SOURCE_ROOT = REPOSITORY_ROOT / "tests" / "data" / "neural_network_parameters"


@dataclass(frozen=True)
class Network:
    """Store one released feed-forward network in compact layer order."""

    name: str
    nodes: tuple[int, ...]
    activation_function: int
    scale_delta: tuple[str, ...]
    scale_average: tuple[str, ...]
    layer_weights: tuple[tuple[str, ...], ...]
    layer_biases: tuple[tuple[str, ...], ...]


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one provenance file."""

    return hashlib.sha256(path.read_bytes()).hexdigest()


def numeric_tokens(line: str) -> tuple[str, ...]:
    """Return whitespace-separated numeric tokens from one released line."""

    tokens = tuple(line.replace("D", "e").replace("d", "e").split())
    for token in tokens:
        float(token)
    return tokens


def consume_values(
    lines: Sequence[str], index: int, count: int
) -> tuple[tuple[str, ...], int]:
    """Consume exactly ``count`` numeric values from complete input lines."""

    values: List[str] = []
    while len(values) < count:
        tokens = numeric_tokens(lines[index])
        if len(values)+len(tokens) > count:
            raise ValueError("numeric record crosses a model boundary")
        values.extend(tokens)
        index += 1
    return tuple(values), index


def parse_weight_file(path: Path, input_count: int) -> tuple[Network, ...]:
    """Parse three released networks and their scaling/weight records."""

    lines = path.read_text().splitlines()
    index = 0
    partial = []
    for _ in range(3):
        name = lines[index].strip()
        index += 1
        header = tuple(int(value) for value in lines[index].split())
        index += 1
        if len(header) == 2:
            hidden_layer_count, output_count = header
        elif len(header) == 3:
            file_input_count, hidden_layer_count, output_count = header
            if file_input_count != input_count:
                raise ValueError("network input count does not match the PES")
        else:
            raise ValueError("unexpected network header")
        hidden_nodes = tuple(int(value) for value in lines[index].split())
        index += 1
        if len(hidden_nodes) != hidden_layer_count:
            raise ValueError("hidden-layer count does not match its node list")
        nodes = (input_count, *hidden_nodes, output_count)
        activation_function, parameter_count = (
            int(value) for value in lines[index].split()
        )
        index += 1
        scale_count = input_count+output_count
        scale_delta, index = consume_values(lines, index, scale_count)
        scale_average, index = consume_values(lines, index, scale_count)
        layer_weights = []
        weight_count = 0
        for previous_count, current_count in zip(nodes[:-1], nodes[1:]):
            weights, index = consume_values(
                lines, index, previous_count*current_count
            )
            layer_weights.append(weights)
            weight_count += len(weights)
        bias_count = sum(nodes[1:])
        if weight_count+bias_count != parameter_count:
            raise ValueError("released network parameter count is inconsistent")
        partial.append((
            name,
            nodes,
            activation_function,
            scale_delta,
            scale_average,
            tuple(layer_weights),
        ))
    if index != len(lines):
        raise ValueError("unparsed records remain in the weights file")
    return tuple(
        Network(*items, layer_biases=()) for items in partial
    )


def attach_biases(
    path: Path, networks: Sequence[Network]
) -> tuple[Network, ...]:
    """Attach released bias vectors to parsed weight networks."""

    lines = path.read_text().splitlines()
    output = []
    for network in networks:
        bias_label = network.name.replace("weights", "biases").replace(
            "weight", "biases"
        )
        try:
            index = lines.index(bias_label)+1
        except ValueError as error:
            raise ValueError("missing {}".format(bias_label)) from error
        layers = []
        for node_count in network.nodes[1:]:
            values, index = consume_values(lines, index, node_count)
            layers.append(values)
        output.append(Network(
            network.name,
            network.nodes,
            network.activation_function,
            network.scale_delta,
            network.scale_average,
            network.layer_weights,
            tuple(layers),
        ))
    return tuple(output)


def released_vector(path: Path, label: str, count: int) -> tuple[str, ...]:
    """Read a named vector from a released biases/provenance file."""

    lines = path.read_text().splitlines()
    try:
        index = lines.index(label)+1
    except ValueError as error:
        raise ValueError("missing {}".format(label)) from error
    values, _ = consume_values(lines, index, count)
    return values


def real_literal(value: str) -> str:
    """Return one released decimal as an explicit real64 literal."""

    return value.lower().replace("e", "e")+"_real64"


def array_declaration(
    name: str,
    dimensions: str,
    values: Sequence[str],
    shape: Optional[str] = None,
) -> str:
    """Render one compact real64 parameter array."""

    lines = [
        "    real(real64), parameter :: {}({}) = {}[ &".format(
            name, dimensions, "reshape(" if shape else ""
        )
    ]
    literals = [real_literal(value) for value in values]
    for index in range(0, len(literals), 2):
        chunk = ", ".join(literals[index:index+2])
        suffix = ", &" if index+2 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    ending = "        ], [{}])".format(shape) if shape else "        ]"
    lines.append(ending)
    return "\n".join(lines)


def render_module(
    module_name: str,
    networks: Sequence[Network],
    weights_path: Path,
    biases_path: Path,
    system_declarations: Sequence[str],
    extra_arrays: Sequence[tuple[str, str, Sequence[str]]],
) -> str:
    """Render one private-by-default immutable parameter module."""

    nodes = networks[0].nodes
    if any(network.nodes != nodes for network in networks):
        raise ValueError("the three ensemble networks require one architecture")
    activation = networks[0].activation_function
    if any(network.activation_function != activation for network in networks):
        raise ValueError("the three ensemble networks require one activation")
    network_count = len(networks)
    input_count, hidden_1_count, hidden_2_count, output_count = nodes
    scale_count = input_count+output_count
    public_names = [
        "activation_function", "hidden_1_biases", "hidden_2_biases",
        "hidden_1_count", "hidden_2_count", "hidden_hidden_weights",
        "hidden_output_weights", "input_hidden_weights", "network_count",
        "ninput", "nlayer", "nodemax", "nodes", "noutput", "nscale",
        "output_biases", "scale_averages", "scale_deltas",
    ]
    for declaration in system_declarations:
        public_names.append(declaration.split("::", 1)[1].split("=", 1)[0].strip())
    public_names.extend(name for name, _, _ in extra_arrays)
    model_label = module_name
    if model_label.startswith("chempotpy_"):
        model_label = model_label[len("chempotpy_"):]
    if model_label.endswith("_parameters"):
        model_label = model_label[:-len("_parameters")]
    output = [
        "!> @brief Store immutable neural-network parameters for the {} model.".format(
            model_label
        ),
        "!",
        "! Generated from released model data; do not edit numeric values by hand.",
        "! Weights SHA-256: {}".format(sha256(weights_path)),
        "! Biases SHA-256: {}".format(sha256(biases_path)),
        "module {}".format(module_name),
        "    use iso_fortran_env, only: real64",
        "    implicit none",
        "    private",
    ]
    for index in range(0, len(public_names), 5):
        prefix = "    public :: " if index == 0 else "        & "
        suffix = ", &" if index+5 < len(public_names) else ""
        output.append(prefix+", ".join(public_names[index:index+5])+suffix)
    output.extend([
        "    integer, parameter :: network_count       = {}".format(network_count),
        "    integer, parameter :: ninput              = {}".format(input_count),
        "    integer, parameter :: hidden_1_count      = {}".format(hidden_1_count),
        "    integer, parameter :: hidden_2_count      = {}".format(hidden_2_count),
        "    integer, parameter :: noutput             = {}".format(output_count),
        "    integer, parameter :: nlayer              = {}".format(len(nodes)),
        "    integer, parameter :: nscale              = ninput+noutput",
        "    integer, parameter :: nodemax             = {}".format(max(nodes)),
        "    integer, parameter :: activation_function = {}".format(activation),
        "    integer, parameter :: nodes(nlayer) = [{}]".format(
            ", ".join(str(value) for value in nodes)
        ),
        *["    "+declaration for declaration in system_declarations],
    ])
    scale_delta = tuple(
        value for network in networks for value in network.scale_delta
    )
    scale_average = tuple(
        value for network in networks for value in network.scale_average
    )
    output.append(array_declaration(
        "scale_deltas", "nscale,network_count", scale_delta,
        "nscale,network_count"
    ))
    output.append(array_declaration(
        "scale_averages", "nscale,network_count", scale_average,
        "nscale,network_count"
    ))
    layer_names = (
        ("input_hidden_weights", "ninput,hidden_1_count,network_count"),
        ("hidden_hidden_weights", "hidden_1_count,hidden_2_count,network_count"),
        ("hidden_output_weights", "hidden_2_count,noutput,network_count"),
    )
    bias_names = (
        ("hidden_1_biases", "hidden_1_count,network_count"),
        ("hidden_2_biases", "hidden_2_count,network_count"),
        ("output_biases", "noutput,network_count"),
    )
    for layer_index, (name, dimensions) in enumerate(layer_names):
        values = tuple(
            value for network in networks
            for value in network.layer_weights[layer_index]
        )
        output.append(array_declaration(
            name, dimensions, values, dimensions
        ))
    for layer_index, (name, dimensions) in enumerate(bias_names):
        values = tuple(
            value for network in networks
            for value in network.layer_biases[layer_index]
        )
        output.append(array_declaration(
            name, dimensions, values, dimensions
        ))
    for name, dimensions, values in extra_arrays:
        output.append(array_declaration(name, dimensions, values))
    output.append("end module {}".format(module_name))
    source = "\n".join(output)+"\n"
    for _ in range(4):
        formatted = format_source_text(source)
        if formatted == source:
            return formatted
        source = formatted
    raise RuntimeError("generated Fortran formatting did not converge")


def generated_documents() -> dict[Path, str]:
    """Return every generated module and its deterministic contents."""

    n4_root = SOURCE_ROOT / "N4"
    n4_weights = n4_root / "weights.txt-Noconnected"
    n4_biases = n4_root / "biases.txt-Noconnected"
    n4_networks = attach_biases(
        n4_biases, parse_weight_file(n4_weights, input_count=276)
    )
    n4 = render_module(
        "chempotpy_n4_pipnn_parameters",
        n4_networks,
        n4_weights,
        n4_biases,
        (
            "real(real64), parameter :: alpha = 1.0_real64",
            "real(real64), parameter :: pi = 3.1415926_real64",
            "real(real64), parameter :: radian = pi/180.0_real64",
            "integer, parameter :: natom = 4",
            "integer, parameter :: nbasis = 306",
            "integer, parameter :: ndim = 6",
            "integer, parameter :: nnew = 276",
        ),
        (),
    )
    fh2o_root = SOURCE_ROOT / "FH2O"
    fh2o_weights = fh2o_root / "weights.txt"
    fh2o_biases = fh2o_root / "biases.txt"
    fh2o_networks = attach_biases(
        fh2o_biases, parse_weight_file(fh2o_weights, input_count=17)
    )
    fh2o = render_module(
        "chempotpy_fh2o_pipnn_parameters",
        fh2o_networks,
        fh2o_weights,
        fh2o_biases,
        (
            "real(real64), parameter :: alpha = 1.0_real64",
            "real(real64), parameter :: vpesmin = -176.0260035639990_real64",
            "real(real64), parameter :: pi = 3.141592653589793238_real64",
            "real(real64), parameter :: radian = pi/180.0_real64",
            "real(real64), parameter :: bohr = 0.5291772_real64",
            "integer, parameter :: nbasis = 18",
            "integer, parameter :: nrbasissoc = 147",
            "integer, parameter :: npbasissoc = 288",
        ),
        (
            ("r1soccoef", "nrbasissoc", released_vector(
                fh2o_biases, "r1soccoef", 147
            )),
            ("p2soccoef", "npbasissoc", released_vector(
                fh2o_biases, "p2soccoef", 288
            )),
        ),
    )
    return {
        REPOSITORY_ROOT / "chempotpy" / "N4" /
        "chempotpy_n4_pipnn_parameters.f90": n4,
        REPOSITORY_ROOT / "chempotpy" / "FH2O" /
        "chempotpy_fh2o_pipnn_parameters.f90": fh2o,
    }


def main() -> int:
    """Write or verify the embedded parameter modules."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    failures = []
    for path, content in generated_documents().items():
        if arguments.check:
            if not path.is_file() or path.read_text() != content:
                failures.append(path)
        else:
            path.write_text(content)
            print(path.relative_to(REPOSITORY_ROOT))
    if failures:
        for path in failures:
            print("outdated: {}".format(path.relative_to(REPOSITORY_ROOT)))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
