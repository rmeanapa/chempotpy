#!/usr/bin/env python3
"""Generate miscellaneous immutable PES parameter modules."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

PROVENANCE_ROOT = REPOSITORY_ROOT / "tests/data/external_parameters"


def real_literal(value: str) -> str:
    """Return one numeric token as an explicit wp literal."""

    normalized = value.strip().lower().replace("d", "e")
    if "." not in normalized and "e" not in normalized:
        normalized += ".0"
    return normalized+"_wp"


def stable_format(source: str) -> str:
    """Return generator-controlled source without reparsing large arrays."""

    return source


def f2h2_module() -> str:
    """Render the released F2H2 JBKKL coefficient vector."""

    source = PROVENANCE_ROOT / "F2H2/bjkkl.par"
    released_values = source.read_text().split()
    if len(released_values) != 360:
        raise ValueError("F2H2/bjkkl.par must contain 72 five-value records")
    # The released reader requested 359 values and intentionally left the
    # final field of the last record unused.
    values = released_values[:359]
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store immutable coefficients for F2H2_GEN_JBKKL_1990.",
        "!",
        "! Generated from the released bjkkl.par model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_f2h2_jbkkl_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: jbkkl_coefficients",
        "    real(wp), parameter :: jbkkl_coefficients(359) = [ &",
    ]
    literals = [real_literal(value) for value in values]
    for index in range(0, len(literals), 3):
        chunk = ", ".join(literals[index:index+3])
        suffix = ", &" if index+3 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    lines.extend([
        "        ]",
        "end module chempotpy_f2h2_jbkkl_parameters",
    ])
    return stable_format("\n".join(lines)+"\n")


def real_array(name: str, values: list[str]) -> list[str]:
    """Render one wp parameter vector."""

    lines = [
        "    real(wp), parameter :: {}({}) = [ &".format(
            name, len(values)
        )
    ]
    literals = [real_literal(value) for value in values]
    for index in range(0, len(literals), 3):
        chunk = ", ".join(literals[index:index+3])
        suffix = ", &" if index+3 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    lines.append("        ]")
    return lines


def real_shaped_array(
    name: str, dimensions: str, values: list[str]
) -> list[str]:
    """Render one multidimensional wp parameter array."""

    lines = [
        "    real(wp), parameter :: {}({}) = reshape([ &".format(
            name, dimensions
        )
    ]
    literals = [real_literal(value) for value in values]
    for index in range(0, len(literals), 3):
        chunk = ", ".join(literals[index:index+3])
        suffix = ", &" if index+3 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    lines.extend(["        ], [{}])".format(dimensions)])
    return lines


def integer_shaped_array(
    name: str, dimensions: str, values: list[int]
) -> list[str]:
    """Render one multidimensional integer parameter array."""

    lines = [
        "    integer, parameter :: {}({}) = reshape([ &".format(
            name, dimensions
        )
    ]
    literals = [str(value) for value in values]
    for index in range(0, len(literals), 12):
        chunk = ", ".join(literals[index:index+12])
        suffix = ", &" if index+12 < len(literals) else " &"
        lines.append("        "+chunk+suffix)
    lines.append("        ], [{}])".format(dimensions))
    return lines


def arno_module() -> str:
    """Render the released ArNO spline grid as immutable parameters."""

    source = PROVENANCE_ROOT / "ArNO/spline.dat"
    rows = [line.split() for line in source.read_text().splitlines() if line.split()]
    if any(len(row) != 5 for row in rows):
        raise ValueError("ArNO spline.dat must contain five-value records")
    nroh1, nth, nroh2, state_count = 17, 13, 70, 2
    vmin = "-656.852768990000"
    row_index = 0
    counts = [[0 for _ in range(nth)] for _ in range(nroh1)]
    roh1 = ["0.0" for _ in range(nroh1)]
    thth = ["0.0" for _ in range(nth)]
    roh2 = [[[
        "0.0" for _ in range(nroh2)
    ] for _ in range(nth)] for _ in range(nroh1)]
    energies = [[[ [
        "0.0" for _ in range(state_count)
    ] for _ in range(nroh2)] for _ in range(nth)] for _ in range(nroh1)]
    for i in range(nroh1):
        for j in range(nth):
            reference_r1 = None
            reference_theta = None
            for k in range(nroh2):
                if row_index >= len(rows):
                    raise ValueError("ArNO spline grid ended before its sentinel")
                r1, theta, r2, energy_1, energy_2 = rows[row_index]
                row_index += 1
                for value in (r1, theta, r2, energy_1, energy_2):
                    float(value)
                if float(r1)+float(r2) < 0.1:
                    break
                if reference_r1 is None:
                    reference_r1 = r1
                    reference_theta = theta
                elif abs(float(r1)-float(reference_r1)) + abs(
                    float(theta)-float(reference_theta)
                ) > 0.01:
                    raise ValueError("ArNO spline group has inconsistent coordinates")
                counts[i][j] += 1
                roh1[i] = r1
                thth[j] = theta
                roh2[i][j][k] = r2
                energies[i][j][k] = [
                    format(
                        min((float(energy)-float(vmin))*27.2114, 12.0),
                        ".17e",
                    )
                    for energy in (energy_1, energy_2)
                ]
    if row_index != len(rows):
        raise ValueError("ArNO spline grid contains trailing records")
    count_values = [counts[i][j] for j in range(nth) for i in range(nroh1)]
    roh2_values = [
        roh2[i][j][k]
        for k in range(nroh2)
        for j in range(nth)
        for i in range(nroh1)
    ]
    energy_values = [
        energies[i][j][k][state]
        for k in range(nroh2)
        for j in range(nth)
        for i in range(nroh1)
        for state in range(state_count)
    ]
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store the immutable released spline grid for both ArNO surfaces.",
        "!",
        "! Generated from the released spline.dat model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_arno_spline_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: bvcut_vcut, pesa_ind, pesa_roh1, pesa_roh2, pesa_thth, pesa_ve",
        "    integer, parameter :: nroh1 = 17",
        "    integer, parameter :: nth = 13",
        "    integer, parameter :: nroh2 = 70",
        "    integer, parameter :: state_count = 2",
        "    real(wp), parameter :: bvcut_vcut = 12.0_wp",
    ]
    lines.extend(integer_shaped_array("pesa_ind", "nroh1,nth", count_values))
    lines.extend(real_array("pesa_roh1", roh1))
    lines.extend(real_shaped_array("pesa_roh2", "nroh1,nth,nroh2", roh2_values))
    lines.extend(real_array("pesa_thth", thth))
    lines.extend(real_shaped_array(
        "pesa_ve", "state_count,nroh1,nth,nroh2", energy_values
    ))
    lines.append("end module chempotpy_arno_spline_parameters")
    return stable_format("\n".join(lines)+"\n")


def hobr_module() -> str:
    """Render the released HOBr RKHS parameters."""

    source = PROVENANCE_ROOT / "HOBr/o_hbr_a_pp.par"
    values = source.read_text().split()
    if len(values) != 14+14+4*1110:
        raise ValueError("HOBr parameter record count is inconsistent")
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    records = values[28:]
    arrays = {
        "phi_oh": values[:14],
        "phi_hbr": values[14:28],
        "rkhs_x_1": records[0::4],
        "rkhs_x_2": records[1::4],
        "rkhs_y_1": records[2::4],
        "rkhs_coefficients": records[3::4],
    }
    lines = [
        "!> @brief Store immutable RKHS parameters for HOBr_3App_PIP.",
        "!",
        "! Generated from the released o_hbr_a_pp.par model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_hobr_3app_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: phi_oh, phi_hbr, rkhs_x_1, rkhs_x_2, rkhs_y_1, &",
        "        & rkhs_coefficients",
    ]
    for name, array in arrays.items():
        lines.extend(real_array(name, array))
    lines.append("end module chempotpy_hobr_3app_parameters")
    return stable_format("\n".join(lines)+"\n")


def coco_module() -> str:
    """Render the released COCO neural-network parameters."""

    source = PROVENANCE_ROOT / "COCO/nn_ococ_w20.txt"
    tokens = source.read_text().split()
    values = tokens[:2492]
    if len(values) != 2492:
        raise ValueError("COCO neural-network data is incomplete")
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    offset = 0

    def take(count: int) -> list[str]:
        nonlocal offset
        result = values[offset:offset+count]
        offset += count
        return result

    arrays = (
        ("input_weights", "7,45", take(7*45)),
        ("hidden_1_biases", "45", take(45)),
        ("hidden_weights", "45,45", take(45*45)),
        ("hidden_2_biases", "45", take(45)),
        ("output_weights", "45", take(45)),
    )
    output_bias = take(1)[0]
    input_ranges = take(2*7)
    output_range = take(2)
    if offset != len(values):
        raise ValueError("COCO neural-network parameter count is inconsistent")
    lines = [
        "!> @brief Store immutable neural-network parameters for COCO_PIPNN.",
        "!",
        "! Generated from the released nn_ococ_w20.txt model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_coco_pipnn_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: input_weights, hidden_1_biases, hidden_weights, &",
        "        & hidden_2_biases, output_weights, output_bias, input_ranges, &",
        "        & output_range",
    ]
    for name, dimensions, array in arrays:
        if "," in dimensions:
            lines.extend(real_shaped_array(name, dimensions, array))
        else:
            lines.extend(real_array(name, array))
    lines.append(
        "    real(wp), parameter :: output_bias = {}".format(
            real_literal(output_bias)
        )
    )
    lines.extend(real_shaped_array("input_ranges", "2,7", input_ranges))
    lines.extend(real_array("output_range", output_range))
    lines.append("end module chempotpy_coco_pipnn_parameters")
    return stable_format("\n".join(lines)+"\n")


def c3h7no_module() -> str:
    """Render the released C3H7NO PIP coefficient vector."""

    source = PROVENANCE_ROOT / "C3H7NO/coeff.dat"
    values = source.read_text().split()
    if len(values) != 8040:
        raise ValueError("C3H7NO coefficient count is inconsistent")
    for value in values:
        float(value)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store immutable coefficients for C3H7NO_PIP.",
        "!",
        "! Generated from the released coeff.dat model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_c3h7no_pip_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: pip_coefficients",
    ]
    lines.extend(real_array("pip_coefficients", values))
    lines.append("end module chempotpy_c3h7no_pip_parameters")
    return stable_format("\n".join(lines)+"\n")


def h4o2_module() -> str:
    """Render the released H4O2 four-body PIP coefficients."""

    source = PROVENANCE_ROOT / "H4O2/h4o2.pes4.coeff.dat"
    lines_in = source.read_text().splitlines()
    values = " ".join(lines_in[2:]).split()[:5227]
    if len(values) != 5227:
        raise ValueError("H4O2 four-body coefficient count is inconsistent")
    for value in values:
        float(value.lower().replace("d", "e"))
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store immutable coefficients for H4O2_PIP.",
        "!",
        "! Generated from the released h4o2.pes4.coeff.dat model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_h4o2_pip_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: pip_coefficients",
    ]
    lines.extend(real_array("pip_coefficients", values))
    lines.append("end module chempotpy_h4o2_pip_parameters")
    return stable_format("\n".join(lines)+"\n")


def h2co_module() -> str:
    """Render the released H2CO PIP and HMK parameters."""

    root = PROVENANCE_ROOT / "H2CO"
    coefficient_source = root / "coeff8"
    index_source = root / "index8"
    hmk_source = root / "pothco_hmk.para"
    coefficients = coefficient_source.read_text().splitlines()[1:1562]
    if len(coefficients) != 1561:
        raise ValueError("H2CO coeff8 must contain 1561 fitted coefficients")
    index_rows = [line.split() for line in index_source.read_text().splitlines()]
    if len(index_rows) != 1561 or any(len(row) != 6 for row in index_rows):
        raise ValueError("H2CO index8 dimensions are inconsistent")
    powers = [row[column] for column in range(6) for row in index_rows]
    numeric_rows = [
        line.split() for line in hmk_source.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(numeric_rows) != 712:
        raise ValueError("H2CO HMK parameter record count is inconsistent")
    integer_parameters = [int(value) for row in numeric_rows[:7] for value in row]
    real_parameters = numeric_rows[7]
    offset = 8
    kk_columns = []
    aa_columns = []
    for linear_count, tangent_count in ((230, 132), (214, 128)):
        linear_rows = numeric_rows[offset:offset+linear_count]
        offset += linear_count
        tangent_rows = numeric_rows[offset:offset+tangent_count]
        offset += tangent_count
        if [int(row[0]) for row in linear_rows] != list(range(1, linear_count+1)):
            raise ValueError("H2CO HMK linear indices are inconsistent")
        if [int(row[0]) for row in tangent_rows] != list(range(1, tangent_count+1)):
            raise ValueError("H2CO HMK tangent indices are inconsistent")
        kk_columns.append([row[1] for row in linear_rows]+["0.0"]*(230-linear_count))
        aa_columns.append([row[1] for row in tangent_rows]+["0.0"]*(132-tangent_count))
    if offset != len(numeric_rows):
        raise ValueError("unparsed H2CO HMK records remain")
    names = (
        "gzhc", "gzco", "gzwin", "gnhc", "gnco", "gnwin", "ixhc", "ixco", "ixoh",
        "gwhc", "gwco", "gwwin", "gz2hc", "gn2hc", "ix2hc", "gz2co", "gn2co", "ix2co",
        "gz2oh", "gn2oh", "ix2oh",
    )
    digest_lines = [
        "! coeff8 SHA-256: {}".format(hashlib.sha256(coefficient_source.read_bytes()).hexdigest()),
        "! index8 SHA-256: {}".format(hashlib.sha256(index_source.read_bytes()).hexdigest()),
        "! pothco_hmk.para SHA-256: {}".format(hashlib.sha256(hmk_source.read_bytes()).hexdigest()),
    ]
    lines = [
        "!> @brief Store immutable coefficients and HMK parameters for H2CO_PIP.",
        "!",
        "! Generated from the released H2CO model data.",
        *digest_lines,
        "module chempotpy_h2co_pip_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: aa, easym, fitted_coefficients, gnhc, gnco, gnwin, gwhc, gwco, gwwin, &",
        "        & gn2co, gn2hc, gn2oh, gzhc, gzco, gzwin, gz2co, gz2hc, gz2oh, ixco, ixhc, &",
        "        & ixoh, ix2co, ix2hc, ix2oh, kk, nenadd, powers, wuradd",
        "    integer, parameter :: " + ", ".join(
            "{} = {}".format(name, value) for name, value in zip(names, integer_parameters)
        ),
        "    real(wp), parameter :: wuradd = {}, nenadd = {}, easym = {}".format(
            *(real_literal(value) for value in real_parameters)
        ),
    ]
    lines.extend(real_array("fitted_coefficients", coefficients))
    lines.extend(real_shaped_array("powers", "1561,6", powers))
    lines.extend(real_shaped_array("kk", "230,2", kk_columns[0]+kk_columns[1]))
    lines.extend(real_shaped_array("aa", "132,2", aa_columns[0]+aa_columns[1]))
    lines.append("end module chempotpy_h2co_pip_parameters")
    return stable_format("\n".join(lines)+"\n")


def ch4_module() -> str:
    """Render the released CH4 T8 coefficient expansion."""

    source = PROVENANCE_ROOT / "CH4/ch4pes"
    source_lines = source.read_text().splitlines()
    base_1 = source_lines[1].replace(",", " ").split()
    base_2 = source_lines[2].replace(",", " ").split()
    base_3 = source_lines[3].replace(",", " ").split()
    setup = source_lines[4].replace(",", " ").split()
    if tuple(map(len, (base_1, base_2, base_3, setup))) != (3, 4, 5, 6):
        raise ValueError("CH4 T8 header dimensions are inconsistent")

    def as_float(token: str) -> float:
        return float(token.lower().replace("d", "e"))

    reau = as_float(setup[0])*1.889725989
    ex, exa, beta = setup[1], setup[2], setup[3]
    lm = int(setup[5])
    cs = [0.0]*6
    cb = [0.0]*6
    terms = [[] for _ in range(8)]
    constant = 0.0
    for line in source_lines[5:]:
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 10:
            raise ValueError("CH4 T8 coefficient record has an unexpected width")
        exponents = [int(value) for value in fields[:9]]
        coefficient = as_float(fields[9])
        internal_sum = sum(exponents[2:])
        if internal_sum+exponents[1] == 0 and lm != 0:
            cs[exponents[0]-1] = coefficient*6.0
            coefficient = 0.0
        elif internal_sum+exponents[0] == 0 and lm != 0:
            cb[exponents[1]-1] = coefficient*4.0
            coefficient = 0.0
        shifted = [value+1 for value in exponents]
        nonzero = [(value, index+1) for index, value in enumerate(shifted) if value > 1]
        if not nonzero:
            constant = coefficient
        else:
            terms[len(nonzero)-1].append((coefficient, nonzero))
    nx = [len(items) for items in terms]
    max_terms = max(nx)
    coefficients = []
    exponent_indices = []
    coordinate_indices = []
    term_offsets = []
    for degree in range(8):
        coefficients.extend(
            [repr(item[0]) for item in terms[degree]]
            + ["0.0"]*(max_terms-len(terms[degree]))
        )
    for degree in range(8):
        term_offsets.append(len(exponent_indices))
        for item in terms[degree]:
            exponent_indices.extend(pair[0] for pair in item[1])
            coordinate_indices.extend(pair[1] for pair in item[1])
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store immutable coefficients for the CH4 Schwenke-Partridge T8 surface.",
        "!",
        "! Generated from the released ch4pes model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_ch4_t8_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: a02, a03, a04, a06, a08, beta, beta0, cb, coefficients, constant, &",
        "        & coordinate_indices, cs, ex, ex0, exa, exponent_indices, max_terms, nx, &",
        "        & term_offsets, &",
        "        & r01, r02, r03, r04, re0, reau",
        "    integer, parameter :: max_terms = {}".format(max_terms),
        "    integer, parameter :: nx(8) = [{}]".format(", ".join(map(str, nx))),
        "    integer, parameter :: term_offsets(8) = [{}]".format(
            ", ".join(map(str, term_offsets))
        ),
        "    real(wp), parameter :: re0 = {}, ex0 = {}, beta0 = {}".format(
            *(real_literal(value) for value in base_1)
        ),
        "    real(wp), parameter :: r01 = {}, r02 = {}, r03 = {}, r04 = {}".format(
            *(real_literal(value) for value in base_2)
        ),
        "    real(wp), parameter :: a02 = {}, a03 = {}, a04 = {}, a06 = {}, a08 = {}".format(
            *(real_literal(value) for value in base_3)
        ),
        "    real(wp), parameter :: reau = {}, ex = {}, exa = {}, beta = {}".format(
            real_literal(repr(reau)), real_literal(ex), real_literal(exa), real_literal(beta)
        ),
        "    real(wp), parameter :: constant = {}".format(real_literal(repr(constant))),
    ]
    lines.extend(real_array("cs", [repr(value) for value in cs]))
    lines.extend(real_array("cb", [repr(value) for value in cb]))
    lines.extend(real_shaped_array("coefficients", "max_terms,8", coefficients))
    lines.extend(integer_shaped_array(
        "exponent_indices", str(len(exponent_indices)), exponent_indices
    ))
    lines.extend(integer_shaped_array(
        "coordinate_indices", str(len(coordinate_indices)), coordinate_indices
    ))
    lines.append("end module chempotpy_ch4_t8_parameters")
    return stable_format("\n".join(lines)+"\n")


def n4_bdtc_module() -> str:
    """Render the released N4 BDTC local-fit database as immutable columns."""

    source = PROVENANCE_ROOT / "N4/pes.dat"
    records = source.read_text().splitlines()
    index = 5+20+3
    moltype, pesform = records[index].split(); index += 2
    n, m, eta = records[index].split(); index += 2
    pnu, epsnu = records[index].split(); index += 2
    morse_widths = records[index].split(); index += 2
    equilibrium_distances = records[index].split(); index += 4
    nma, nplmax = records[index].split(); index += 2
    index += 2
    nmincrc, mcrc = records[index].split(); index += 2
    crc = records[index].split(); index += 2
    target = records[index].split(); index += 2
    index += 2
    n_int, m_int, nma_int = int(n), int(m), int(nma)
    geometry_columns = [[] for _ in range(6)]
    for record in records[index:index+n_int+nma_int]:
        fields = record.split()
        if int(fields[-1]) == 0:
            for column in range(6):
                geometry_columns[column].append(fields[2+column])
    if any(len(column) != n_int for column in geometry_columns):
        raise ValueError("N4 unmasked geometry count is inconsistent")
    index += n_int+nma_int+1
    aji_values = [record.strip() for record in records[index:]]
    if len(aji_values) != n_int*m_int:
        raise ValueError("N4 local-fit coefficient count is inconsistent")
    aji_columns = [
        [aji_values[row*m_int+column] for row in range(n_int)]
        for column in range(m_int)
    ]
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    lines = [
        "!> @brief Store the immutable released database for N4_BDTC_2014.",
        "!",
        "! Generated from the released pes.dat model data.",
        "! Source SHA-256: {}".format(digest),
        "module chempotpy_n4_bdtc_parameters",
        "    use chempotpy_kinds, only: wp",
        "    implicit none",
        "    private",
        "    public :: coefficient_value, crc, crcmult, epsnu, eta, equilibrium_distances, &",
        "        & geometry_value, m, mcrc, moltype, morse_widths, n, nitcrc, nma, nmincrc, &",
        "        & nplmax, ntar, ntarmax, ntarmin, pesform, pnu, rwfa, rwfb",
        "    integer, parameter :: moltype = {}, pesform = {}".format(moltype, pesform),
        "    integer, parameter :: n = {}, m = {}, eta = {}".format(n, m, eta),
        "    integer, parameter :: nma = {}, nplmax = {}".format(nma, nplmax),
        "    integer, parameter :: nmincrc = {}, mcrc = {}".format(nmincrc, mcrc),
        "    integer, parameter :: ntarmin = {}, ntarmax = {}, ntar = {}, crcmult = {}, nitcrc = {}".format(
            *target[:5]
        ),
        "    real(wp), parameter :: pnu = {}, epsnu = {}, rwfa = {}, rwfb = {}".format(
            real_literal(pnu), real_literal(epsnu), real_literal(target[5]), real_literal(target[6])
        ),
    ]
    lines.extend(real_array("morse_widths", morse_widths))
    lines.extend(real_array("equilibrium_distances", equilibrium_distances))
    lines.extend(real_array("crc", crc))
    for column, values in enumerate(geometry_columns, 1):
        lines.extend(real_array("geometry_column_{:02d}".format(column), values))
    for column, values in enumerate(aji_columns, 1):
        lines.extend(real_array("coefficient_column_{:02d}".format(column), values))
    lines.extend([
        "",
        "    contains",
        "",
        "    !> @brief Return one released unmasked raw-geometry coordinate.",
        "    pure function geometry_value( point_index, coordinate_index ) result(value)",
        "        integer, intent(in) :: point_index, coordinate_index",
        "        real(wp)            :: value",
        "        select case( coordinate_index )",
    ])
    for column in range(1, 7):
        lines.append("        case( {} ); value = geometry_column_{:02d}(point_index)".format(column, column))
    lines.extend([
        "        case default; value = 0.0_wp",
        "        end select",
        "    end function geometry_value",
        "",
        "    !> @brief Return one released local-fit coefficient.",
        "    pure function coefficient_value( point_index, coefficient_index ) result(value)",
        "        integer, intent(in) :: point_index, coefficient_index",
        "        real(wp)            :: value",
        "        select case( coefficient_index )",
    ])
    for column in range(1, m_int+1):
        lines.append("        case( {} ); value = coefficient_column_{:02d}(point_index)".format(column, column))
    lines.extend([
        "        case default; value = 0.0_wp",
        "        end select",
        "    end function coefficient_value",
        "end module chempotpy_n4_bdtc_parameters",
    ])
    return stable_format("\n".join(lines)+"\n")


def generated_documents() -> dict[Path, str]:
    """Return every generated miscellaneous PES parameter module."""

    return {
        REPOSITORY_ROOT / "chempotpy/ArNO/chempotpy_arno_spline_parameters.f90": (
            arno_module()
        ),
        REPOSITORY_ROOT / "chempotpy/F2H2/chempotpy_f2h2_jbkkl_parameters.f90": (
            f2h2_module()
        ),
        REPOSITORY_ROOT / "chempotpy/HOBr/chempotpy_hobr_3app_parameters.f90": (
            hobr_module()
        ),
        REPOSITORY_ROOT / "chempotpy/COCO/chempotpy_coco_pipnn_parameters.f90": (
            coco_module()
        ),
        REPOSITORY_ROOT / "chempotpy/C3H7NO/chempotpy_c3h7no_pip_parameters.f90": (
            c3h7no_module()
        ),
        REPOSITORY_ROOT / "chempotpy/H4O2/chempotpy_h4o2_pip_parameters.f90": (
            h4o2_module()
        ),
        REPOSITORY_ROOT / "chempotpy/H2CO/chempotpy_h2co_pip_parameters.f90": (
            h2co_module()
        ),
        REPOSITORY_ROOT / "chempotpy/CH4/chempotpy_ch4_t8_parameters.f90": (
            ch4_module()
        ),
        REPOSITORY_ROOT / "chempotpy/N4/chempotpy_n4_bdtc_parameters.f90": (
            n4_bdtc_module()
        ),
    }


def main() -> int:
    """Write or verify the generated modules."""

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
    for path in failures:
        print("outdated: {}".format(path.relative_to(REPOSITORY_ROOT)))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
