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
    if any(line.strip(" *\t") for line in lines[index:]):
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


def parse_unlabeled_network(path: Path, input_count: int) -> Network:
    """Parse one released network whose file omits an ensemble label."""

    lines = path.read_text().splitlines()
    index = 0
    while index < len(lines):
        try:
            header = tuple(int(value) for value in lines[index].split())
        except ValueError:
            index += 1
            continue
        if len(header) == 2:
            break
        index += 1
    if index == len(lines):
        raise ValueError("network header was not found")
    hidden_layer_count, output_count = (
        value for value in header
    )
    index += 1
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
    for previous_count, current_count in zip(nodes[:-1], nodes[1:]):
        weights, index = consume_values(
            lines, index, previous_count*current_count
        )
        layer_weights.append(weights)
    if sum(len(values) for values in layer_weights)+sum(nodes[1:]) != parameter_count:
        raise ValueError("released network parameter count is inconsistent")
    if any(line.strip(" *\t") for line in lines[index:]):
        raise ValueError("unparsed records remain in the weights file")
    return Network(
        path.stem,
        nodes,
        activation_function,
        scale_delta,
        scale_average,
        tuple(layer_weights),
        (),
    )


def attach_unlabeled_biases(path: Path, network: Network) -> Network:
    """Attach a numeric-only bias file to one released network."""

    lines = path.read_text().splitlines()
    index = 0
    while index < len(lines):
        try:
            numeric_tokens(lines[index])
        except ValueError:
            index += 1
            continue
        break
    layers = []
    for node_count in network.nodes[1:]:
        values, index = consume_values(lines, index, node_count)
        layers.append(values)
    if any(line.strip(" *\t") for line in lines[index:]):
        raise ValueError("unparsed records remain in the biases file")
    return Network(
        network.name,
        network.nodes,
        network.activation_function,
        network.scale_delta,
        network.scale_average,
        network.layer_weights,
        tuple(layers),
    )


def attach_unlabeled_bias_ensemble(
    path: Path, networks: Sequence[Network]
) -> tuple[Network, ...]:
    """Attach sequential bias records to an unlabeled network ensemble."""

    lines = path.read_text().splitlines()
    index = 0
    output = []
    for network in networks:
        while index < len(lines):
            try:
                numeric_tokens(lines[index])
            except ValueError:
                index += 1
                continue
            break
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
    if any(line.strip(" *\t") for line in lines[index:]):
        raise ValueError("unparsed records remain in the biases file")
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


def leading_numeric_values(path: Path, count: int) -> tuple[str, ...]:
    """Return the leading numeric records from a legacy list-directed file."""

    values = []
    for token in path.read_text().split():
        normalized = token.replace("D", "e").replace("d", "e")
        try:
            float(normalized)
        except ValueError:
            break
        values.append(normalized)
        if len(values) == count:
            return tuple(values)
    raise ValueError(
        "{} contains {} leading values; expected {}".format(
            path, len(values), count
        )
    )


def all_numeric_values(path: Path, count: int) -> tuple[str, ...]:
    """Return every numeric token from a fixed-size provenance file."""

    values = tuple(
        token.replace("D", "e").replace("d", "e")
        for token in path.read_text().split()
    )
    for value in values:
        float(value)
    if len(values) != count:
        raise ValueError(
            "{} contains {} values; expected {}".format(
                path, len(values), count
            )
        )
    return values


def natural_spline_second_derivatives(
    coordinates: Sequence[str], values: Sequence[str]
) -> tuple[str, ...]:
    """Reproduce the released Numerical Recipes natural-spline setup."""

    x = tuple(float(value) for value in coordinates)
    y = tuple(float(value) for value in values)
    count = len(x)
    if len(y) != count or count < 2:
        raise ValueError("spline coordinate and value dimensions are inconsistent")
    second = [0.0]*count
    work = [0.0]*count
    for index in range(1, count-1):
        sigma = (x[index]-x[index-1])/(x[index+1]-x[index-1])
        factor = sigma*second[index-1]+2.0
        second[index] = (sigma-1.0)/factor
        slope_change = (
            (y[index+1]-y[index])/(x[index+1]-x[index])
            -(y[index]-y[index-1])/(x[index]-x[index-1])
        )
        work[index] = (
            6.0*slope_change/(x[index+1]-x[index-1])
            -sigma*work[index-1]
        )/factor
    for index in range(count-2, -1, -1):
        second[index] = second[index]*second[index+1]+work[index]
    return tuple("{:.17e}".format(value) for value in second)


def real_literal(value: str) -> str:
    """Return one released decimal as an explicit wp literal."""

    return value.lower().replace("e", "e")+"_wp"


def array_declaration(
    name: str,
    dimensions: str,
    values: Sequence[str],
    shape: Optional[str] = None,
) -> str:
    """Render one compact wp parameter array."""

    lines = [
        "    real(wp), parameter :: {}({}) = {}[ &".format(
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
    weights_path: Path | Sequence[Path],
    biases_path: Path | Sequence[Path],
    system_declarations: Sequence[str],
    extra_arrays: Sequence[tuple[str, str, Sequence[str]]],
    additional_provenance: Sequence[tuple[str, Path]] = (),
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
    weights_paths = (weights_path,) if isinstance(weights_path, Path) else tuple(weights_path)
    biases_paths = (biases_path,) if isinstance(biases_path, Path) else tuple(biases_path)
    output = [
        "!> @brief Store immutable neural-network parameters for the {} model.".format(
            model_label
        ),
        "!",
        "! Generated from released model data; do not edit numeric values by hand.",
    ]
    output.extend(
        "! Weights {} SHA-256: {}".format(index, sha256(path))
        for index, path in enumerate(weights_paths, 1)
    )
    output.extend(
        "! Biases {} SHA-256: {}".format(index, sha256(path))
        for index, path in enumerate(biases_paths, 1)
    )
    output.extend(
        "! {} SHA-256: {}".format(label, sha256(path))
        for label, path in additional_provenance
    )
    output.extend([
        "module {}".format(module_name),
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
    ])
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
        shape = dimensions if "," in dimensions else None
        output.append(array_declaration(name, dimensions, values, shape))
    output.append("end module {}".format(module_name))
    return "\n".join(output)+"\n"


def render_c2n_module(
    module_name: str, parameter_paths: Sequence[Path]
) -> str:
    """Render one immutable three-network C2N parameter module."""

    model_label = module_name
    if model_label.startswith("chempotpy_"):
        model_label = model_label[len("chempotpy_"):]
    if model_label.endswith("_parameters"):
        model_label = model_label[:-len("_parameters")]
    input_count = 12
    hidden_1_count = 30
    hidden_2_count = 30
    network_count = len(parameter_paths)
    parameter_count = (
        input_count*hidden_1_count+hidden_1_count
        + hidden_1_count*hidden_2_count+hidden_2_count
        + hidden_2_count+1+2*input_count+2
    )
    networks = [
        leading_numeric_values(path, parameter_count)
        for path in parameter_paths
    ]
    offsets = [
        0,
        input_count*hidden_1_count,
        input_count*hidden_1_count+hidden_1_count,
        input_count*hidden_1_count+hidden_1_count
        + hidden_1_count*hidden_2_count,
        input_count*hidden_1_count+hidden_1_count
        + hidden_1_count*hidden_2_count+hidden_2_count,
        input_count*hidden_1_count+hidden_1_count
        + hidden_1_count*hidden_2_count+hidden_2_count+hidden_2_count,
        input_count*hidden_1_count+hidden_1_count
        + hidden_1_count*hidden_2_count+hidden_2_count+hidden_2_count+1,
        input_count*hidden_1_count+hidden_1_count
        + hidden_1_count*hidden_2_count+hidden_2_count+hidden_2_count+1
        + 2*input_count,
        parameter_count,
    ]
    names_and_dimensions = (
        ("input_hidden_weights", "ninput,hidden_1_count,network_count"),
        ("hidden_1_biases", "hidden_1_count,network_count"),
        ("hidden_hidden_weights", "hidden_1_count,hidden_2_count,network_count"),
        ("hidden_2_biases", "hidden_2_count,network_count"),
        ("hidden_output_weights", "hidden_2_count,network_count"),
        ("output_biases", "network_count"),
        ("input_ranges", "2,ninput,network_count"),
        ("output_ranges", "2,network_count"),
    )
    output = [
        "!> @brief Store immutable neural-network parameters for the {} model.".format(
            model_label
        ),
        "!",
        "! Generated from released model data; do not edit numeric values by hand.",
    ]
    output.extend(
        "! Network {} SHA-256: {}".format(index, sha256(path))
        for index, path in enumerate(parameter_paths, 1)
    )
    output.extend([
        "module {}".format(module_name),
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: hidden_1_biases, hidden_1_count, hidden_2_biases, hidden_2_count, &",
        "        & hidden_hidden_weights, hidden_output_weights, input_hidden_weights, input_ranges, &",
        "        & network_count, ninput, output_biases, output_ranges",
        "    integer, parameter :: network_count  = {}".format(network_count),
        "    integer, parameter :: ninput         = {}".format(input_count),
        "    integer, parameter :: hidden_1_count = {}".format(hidden_1_count),
        "    integer, parameter :: hidden_2_count = {}".format(hidden_2_count),
    ])
    for array_index, (name, dimensions) in enumerate(names_and_dimensions):
        values = tuple(
            value
            for network in networks
            for value in network[offsets[array_index]:offsets[array_index+1]]
        )
        output.append(array_declaration(name, dimensions, values, dimensions))
    output.append("end module {}".format(module_name))
    return "\n".join(output)+"\n"


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
            "real(wp), parameter :: alpha = 1.0_wp",
            "real(wp), parameter :: pi = 3.1415926_wp",
            "real(wp), parameter :: radian = pi/180.0_wp",
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
            "real(wp), parameter :: alpha = 1.0_wp",
            "real(wp), parameter :: vpesmin = -176.0260035639990_wp",
            "real(wp), parameter :: pi = 3.141592653589793238_wp",
            "real(wp), parameter :: radian = pi/180.0_wp",
            "real(wp), parameter :: bohr = 0.5291772_wp",
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
    h3clo_root = SOURCE_ROOT / "H3ClO"
    h3clo_weights = h3clo_root / "weights.txt"
    h3clo_biases = h3clo_root / "biases.txt"
    h3clo_networks = attach_biases(
        h3clo_biases, parse_weight_file(h3clo_weights, input_count=230)
    )
    h3clo = render_module(
        "chempotpy_h3clo_pipnn_parameters",
        h3clo_networks,
        h3clo_weights,
        h3clo_biases,
        (
            "integer, parameter :: natom = 5",
            "integer, parameter :: nbond = natom*(natom-1)/2",
            "integer, parameter :: nmorse = 110",
        ),
        (),
    )
    c2n_root = SOURCE_ROOT / "C2N"
    c2n_prime_paths = tuple(
        c2n_root / "c2n_2a_prime_para_{}.txt".format(label)
        for label in "abc"
    )
    c2n_dprime_paths = tuple(
        c2n_root / "c2n_2a_dprime_para_{}.txt".format(label)
        for label in "abc"
    )
    c2n_prime = render_c2n_module(
        "chempotpy_c2n_pipnn_ap_parameters", c2n_prime_paths
    )
    c2n_dprime = render_c2n_module(
        "chempotpy_c2n_pipnn_app_parameters", c2n_dprime_paths
    )
    ch4of_root = SOURCE_ROOT / "CH4OF"
    ch4of_weights = ch4of_root / "weights.txt"
    ch4of_biases = ch4of_root / "biases.txt"
    ch4of_network = attach_unlabeled_biases(
        ch4of_biases,
        parse_unlabeled_network(ch4of_weights, input_count=1104),
    )
    ch4of = render_module(
        "chempotpy_ch4of_pipnn_parameters",
        (ch4of_network,),
        ch4of_weights,
        ch4of_biases,
        (
            "real(wp), parameter :: alpha = 1.0_wp",
            "real(wp), parameter :: vpesmin = -215.234262652419_wp",
            "integer, parameter :: natom = 7",
            "integer, parameter :: nbasis = 1105",
            "integer, parameter :: ndim = 21",
        ),
        (),
    )
    single_network_configurations = (
        (
            "CH2O", "chempotpy_ch2o_pipnn_parameters", 17,
            (
                "real(wp), parameter :: alpha = 1.0_wp",
                "real(wp), parameter :: vpesmin = -114.158373618668_wp",
                "integer, parameter :: natom = 4",
                "integer, parameter :: nbasis = 18",
                "integer, parameter :: ndim = 6",
            ),
        ),
        (
            "CH2O2", "chempotpy_ch2o2_pipnn_parameters", 101,
            (
                "real(wp), parameter :: alpha = 1.8_wp",
                "real(wp), parameter :: para2 = 27.2114_wp",
                "real(wp), parameter :: vpesmin = -560.0228149965_wp",
                "integer, parameter :: natom = 5",
                "integer, parameter :: nbasis = 102",
                "integer, parameter :: nbond = 10",
            ),
        ),
        (
            "CH4OCl", "chempotpy_ch4ocl_pipnn_parameters", 1104,
            (
                "real(wp), parameter :: alpha = 1.0_wp",
                "real(wp), parameter :: vpesmin = -575.2960146494_wp",
                "integer, parameter :: natom = 7",
                "integer, parameter :: nbasis = 1105",
                "integer, parameter :: ndim = 21",
            ),
        ),
        (
            "H3S", "chempotpy_h3s_pipnn_parameters", 22,
            (
                "real(wp), parameter :: alpha = 1.0_wp",
                "real(wp), parameter :: para2 = 27.2116_wp",
                "real(wp), parameter :: vpesmin = -399.45863953_wp",
                "integer, parameter :: natom = 4",
                "integer, parameter :: nbasis = 23",
                "integer, parameter :: ndata = 33435",
                "integer, parameter :: ndim = 6",
            ),
        ),
        (
            "O2H3", "chempotpy_o2h3_pipnn_parameters", 138,
            (
                "real(wp), parameter :: alpha = 1.0_wp",
                "real(wp), parameter :: vpesmin = -152.040123400851_wp",
                "integer, parameter :: natom = 5",
                "integer, parameter :: nbasis = 139",
                "integer, parameter :: ndim = 10",
            ),
        ),
        (
            "C2H2N2O", "chempotpy_c2h2n2o_pipnn_parameters", 456,
            (
                "real(wp), parameter :: alpha = 1.8_wp",
                "real(wp), parameter :: para2 = 27.2114_wp",
                "integer, parameter :: natom = 7",
                "integer, parameter :: nbasis = 2390",
                "integer, parameter :: nbond = 21",
            ),
        ),
    )
    additional_single_networks = {}
    for system, module_name, input_count, declarations in single_network_configurations:
        root = SOURCE_ROOT / system
        weights = root / "weights.txt"
        biases_name = "biase.txt" if system == "H3S" else "biases.txt"
        biases = root / biases_name
        if system == "H3S":
            networks = attach_unlabeled_bias_ensemble(
                biases, parse_weight_file(weights, input_count)
            )
        else:
            networks = (
                attach_unlabeled_biases(
                    biases, parse_unlabeled_network(weights, input_count)
                ),
            )
        additional_single_networks[
            REPOSITORY_ROOT / "chempotpy" / system /
            (module_name+".f90")
        ] = render_module(
            module_name, networks, weights, biases, declarations, ()
        )
    ensemble_configurations = (
        (
            "C2H2O", "chempotpy_c2h2o_pipnn_parameters", 101,
            ("real(wp), parameter :: vpescut = 5.0_wp",), "letters",
        ),
        (
            "HClOH", "chempotpy_hcloh_pipnn_parameters", 17,
            ("real(wp), parameter :: vpescut = 5.0_wp",), "letters",
        ),
        (
            "HO3", "chempotpy_ho3_pipnn_parameters", 22,
            ("real(wp), parameter :: vpescut = 6.0_wp",), "letters",
        ),
        (
            "K2Rb2", "chempotpy_k2rb2_pipnn_parameters", 32,
            ("real(wp), parameter :: vpescut = 5000.0_wp",), "letters",
        ),
        (
            "N2HOCp", "chempotpy_n2hocp_pipnn_parameters", 167,
            (
                "real(wp), parameter :: alpha = 1.5_wp",
                "integer, parameter :: natom = 5",
                "integer, parameter :: nbasis = 167",
                "integer, parameter :: ndim = 10",
            ),
            "numbers",
        ),
    )
    additional_ensembles = {}
    for system, module_name, input_count, declarations, naming in ensemble_configurations:
        root = SOURCE_ROOT / system
        if naming == "letters":
            weights = tuple(root / "weights{}.txt".format(label) for label in "abc")
            biases = tuple(root / "biases{}.txt".format(label) for label in "abc")
        else:
            weights = tuple(root / "weights.txt-{}".format(index) for index in range(1, 4))
            biases = tuple(root / "biases.txt-{}".format(index) for index in range(1, 4))
        networks = tuple(
            attach_unlabeled_biases(
                bias_path, parse_unlabeled_network(weight_path, input_count)
            )
            for weight_path, bias_path in zip(weights, biases)
        )
        additional_ensembles[
            REPOSITORY_ROOT / "chempotpy" / system / (module_name+".f90")
        ] = render_module(
            module_name, networks, weights, biases, declarations, ()
        )
    clh2o_root = SOURCE_ROOT / "ClH2O"
    clh2o_weights = tuple(
        clh2o_root / "weights{}.txt".format(label) for label in "abc"
    )
    clh2o_biases = tuple(
        clh2o_root / "biases{}.txt".format(label) for label in "abc"
    )
    clh2o_networks = tuple(
        attach_unlabeled_biases(
            bias_path, parse_unlabeled_network(weight_path, 17)
        )
        for weight_path, bias_path in zip(clh2o_weights, clh2o_biases)
    )
    clh2o_table_path = (
        REPOSITORY_ROOT / "tests" / "data" / "external_parameters" /
        "ClH2O" / "spline_tables.txt"
    )
    clh2o_values = all_numeric_values(clh2o_table_path, 436)
    rxhcl = clh2o_values[0:18]
    rxho = clh2o_values[18:36]
    exhohcl = clh2o_values[36:360]
    r0hcl = clh2o_values[360:378]
    dp0hcl = clh2o_values[378:396]
    r0ho = clh2o_values[396:416]
    dp0ho = clh2o_values[416:436]
    exhohcl2_rows = tuple(
        natural_spline_second_derivatives(
            rxho, tuple(exhohcl[column*18+row] for column in range(18))
        )
        for row in range(18)
    )
    exhohcl2 = tuple(
        exhohcl2_rows[row][column]
        for column in range(18)
        for row in range(18)
    )
    additional_ensembles[
        REPOSITORY_ROOT / "chempotpy" / "ClH2O" /
        "chempotpy_clh2o_pipnn_parameters.f90"
    ] = render_module(
        "chempotpy_clh2o_pipnn_parameters",
        clh2o_networks,
        clh2o_weights,
        clh2o_biases,
        (
            "real(wp), parameter :: alpha = 1.0_wp",
            "real(wp), parameter :: vpesmin = -536.01052952_wp",
            "real(wp), parameter :: vpescut = 5.0_wp",
            "real(wp), parameter :: pi = 3.141592653589793238_wp",
            "real(wp), parameter :: radian = pi/180.0_wp",
            "real(wp), parameter :: bohr = 0.5291772_wp",
            "integer, parameter :: nbasis = 18",
            "integer, parameter :: n0hcl = 18",
            "integer, parameter :: n0ho = 20",
            "integer, parameter :: nxhcl = 18",
            "integer, parameter :: nxho = 18",
        ),
        (
            ("rxhcl", "nxhcl", rxhcl),
            ("rxho", "nxho", rxho),
            ("exhohcl", "nxhcl,nxho", exhohcl),
            ("exhohcl2", "nxhcl,nxho", exhohcl2),
            ("r0hcl", "n0hcl", r0hcl),
            ("dp0hcl", "n0hcl", dp0hcl),
            (
                "dp0hcl2", "n0hcl",
                natural_spline_second_derivatives(r0hcl, dp0hcl),
            ),
            ("r0ho", "n0ho", r0ho),
            ("dp0ho", "n0ho", dp0ho),
            (
                "dp0ho2", "n0ho",
                natural_spline_second_derivatives(r0ho, dp0ho),
            ),
        ),
        (("Spline tables", clh2o_table_path),),
    )
    h3cl_root = SOURCE_ROOT / "H3Cl"
    h3cl_weights = tuple(
        h3cl_root / "weights.txt-{}".format(index) for index in range(1, 4)
    )
    h3cl_biases = tuple(
        h3cl_root / "biases.txt-{}".format(index) for index in range(1, 4)
    )
    h3cl_networks = tuple(
        attach_unlabeled_biases(
            bias_path, parse_unlabeled_network(weight_path, 50)
        )
        for weight_path, bias_path in zip(h3cl_weights, h3cl_biases)
    )
    h2_values = tuple(
        token
        for line in (h3cl_root / "H2.dat").read_text().splitlines()
        for token in numeric_tokens(line)
    )
    hcl_values = tuple(
        token
        for line in (h3cl_root / "HCl.dat").read_text().splitlines()
        for token in numeric_tokens(line)
    )
    long_range_coefficients = tuple(
        token
        for line in (h3cl_root / "c.dat").read_text().splitlines()
        for token in numeric_tokens(line)
    )
    if len(h2_values) != 2*1425 or len(hcl_values) != 2*3200:
        raise ValueError("H3Cl diatomic tables have unexpected dimensions")
    if len(long_range_coefficients) != 644:
        raise ValueError("H3Cl long-range table has an unexpected dimension")
    additional_ensembles[
        REPOSITORY_ROOT / "chempotpy" / "H3Cl" /
        "chempotpy_h3cl_pipnn_parameters.f90"
    ] = render_module(
        "chempotpy_h3cl_pipnn_parameters",
        h3cl_networks,
        h3cl_weights,
        h3cl_biases,
        (
            "real(wp), parameter :: alpha = 1.5_wp",
            "integer, parameter :: natom = 4",
            "integer, parameter :: nbasis = 50",
            "integer, parameter :: ndim = 6",
        ),
        (
            ("h2_distance_grid", "1425", h2_values[0::2]),
            ("h2_potential_grid", "1425", h2_values[1::2]),
            ("hcl_distance_grid", "3200", hcl_values[0::2]),
            ("hcl_potential_grid", "3200", hcl_values[1::2]),
            ("long_range_coefficients", "644", long_range_coefficients),
        ),
    )
    no2_root = SOURCE_ROOT / "NO2"
    for state in ("2Ap", "4Ap", "6Ap"):
        weights = no2_root / "weights_NO2_{}_PIPNN.txt".format(state)
        biases = no2_root / "biases_NO2_{}_PIPNN.txt".format(state)
        networks = attach_unlabeled_bias_ensemble(
            biases, parse_weight_file(weights, 12)
        )
        module_name = "chempotpy_no2_{}_pipnn_parameters".format(state.lower())
        additional_ensembles[
            REPOSITORY_ROOT / "chempotpy" / "NO2" / (module_name+".f90")
        ] = render_module(
            module_name,
            networks,
            weights,
            biases,
            (
                "real(wp), parameter :: alpha = 1.0_wp",
                "real(wp), parameter :: bohr = 0.5291772_wp",
                "integer, parameter :: natom = 3",
                "integer, parameter :: nbond = 3",
                "integer, parameter :: npes = 3",
            ),
            (),
        )
    documents = {
        REPOSITORY_ROOT / "chempotpy" / "N4" /
        "chempotpy_n4_pipnn_parameters.f90": n4,
        REPOSITORY_ROOT / "chempotpy" / "FH2O" /
        "chempotpy_fh2o_pipnn_parameters.f90": fh2o,
        REPOSITORY_ROOT / "chempotpy" / "H3ClO" /
        "chempotpy_h3clo_pipnn_parameters.f90": h3clo,
        REPOSITORY_ROOT / "chempotpy" / "C2N" /
        "chempotpy_c2n_pipnn_ap_parameters.f90": c2n_prime,
        REPOSITORY_ROOT / "chempotpy" / "C2N" /
        "chempotpy_c2n_pipnn_app_parameters.f90": c2n_dprime,
        REPOSITORY_ROOT / "chempotpy" / "CH4OF" /
        "chempotpy_ch4of_pipnn_parameters.f90": ch4of,
    }
    documents.update(additional_single_networks)
    documents.update(additional_ensembles)
    return documents


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
