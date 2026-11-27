#!/usr/bin/env python3
"""Generate the system-level PES summary table in the project README."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import SurfaceSource, discover_surface_sources


README = REPOSITORY_ROOT / "README.md"
FIXTURE = REPOSITORY_ROOT / "tests" / "data" / "all_surface_references.json"
START_MARKER = "<!-- BEGIN GENERATED PES CATALOG -->"
END_MARKER = "<!-- END GENERATED PES CATALOG -->"
MODEL_FAMILIES = (
    (re.compile(r"3x3_DPEM", re.IGNORECASE), "neural-network coupled"),
    (re.compile(r"sGDML", re.IGNORECASE), "symmetric gradient-domain machine-learning"),
    (re.compile(r"BMKP", re.IGNORECASE), "BMKP global analytic"),
    (re.compile(r"PIPNN", re.IGNORECASE), "PIP neural-network"),
    (re.compile(r"APRP", re.IGNORECASE), "anchor-points reactive"),
    (re.compile(r"VBMM", re.IGNORECASE), "valence-bond molecular-mechanics"),
    (re.compile(r"LEPS|(?:^|_)LEP(?:_|$)", re.IGNORECASE), "LEPS reactive"),
    (re.compile(r"(?:^|_)PIP(?:_|$)", re.IGNORECASE), "permutationally invariant polynomial"),
    (re.compile(r"Surfgen", re.IGNORECASE), "Surfgen multistate"),
    (re.compile(r"(?:^|_)EAM\d*(?:_|$)", re.IGNORECASE), "embedded-atom"),
    (re.compile(r"(?:^|_)VBO\d*(?:_|$)", re.IGNORECASE), "valence-bond-order"),
    (re.compile(r"RKHS", re.IGNORECASE), "RKHS-fitted"),
    (re.compile(r"(?:^|_)MBP(?:_|$)", re.IGNORECASE), "many-body polynomial"),
    (re.compile(r"(?:^|_)ZV(?:_|$)", re.IGNORECASE), "ZV analytic"),
    (re.compile(r"(?:^|_)DIM(?:_|$)", re.IGNORECASE), "diatomics-in-molecules"),
    (re.compile(r"(?:^|_)GEN(?:_|$)", re.IGNORECASE), "general analytic"),
)

SPECIAL_COMPOSITIONS = {
    "AlmHn": "Al<sub>m</sub>H<sub>n</sub>",
    "CH4_H2O_H2O": "CH<sub>4</sub>·(H<sub>2</sub>O)<sub>2</sub>",
    "COCO": "C<sub>2</sub>O<sub>2</sub>",
    "HOCO_Ar": "HOCO·Ar",
    "malonaldehyde": "C<sub>3</sub>H<sub>4</sub>O<sub>2</sub> (malonaldehyde)",
}
GENERIC_SYSTEMS = {"ALC", "MCH", "MXH", "YRH"}
CHARGED_SUFFIXES = {"p": "+", "m": "−"}


def initial_capital(text: str) -> str:
    """Capitalize a sentence without lowercasing scientific acronyms."""

    return text[:1].upper()+text[1:]


def chemical_composition(system_name: str) -> str:
    """Render a system identifier as a compact HTML chemical composition."""

    if system_name in SPECIAL_COMPOSITIONS:
        return SPECIAL_COMPOSITIONS[system_name]
    if system_name in GENERIC_SYSTEMS:
        return system_name
    charge = ""
    formula = system_name
    if formula[-1:] in CHARGED_SUFFIXES:
        charge = "<sup>{}</sup>".format(CHARGED_SUFFIXES[formula[-1]])
        formula = formula[:-1]
    tokens = re.findall(r"[A-Z][a-z]?|\d+|[mn]", formula)
    if "".join(tokens) != formula:
        return system_name+charge
    rendered = "".join(
        "<sub>{}</sub>".format(token)
        if token.isdigit() or token in {"m", "n"}
        else token
        for token in tokens
    )
    return rendered+charge


def surface_description(source: SurfaceSource) -> str:
    """Return a concise description for one callable PES or DPEM entry."""

    family = "analytic"
    for pattern, description in MODEL_FAMILIES:
        if pattern.search(source.surface_name):
            family = description
            break
    representation = (
        "diabatic representation"
        if source.entry_point == "dpem"
        else "potential"
    )
    return "{} {}".format(initial_capital(family), representation)


def generated_catalog() -> str:
    """Render the grouped Markdown surface table."""

    sources = tuple(sorted(
        (
            source for source in discover_surface_sources()
            if source.is_surface
        ),
        key=lambda source: (
            source.system_name.casefold(), source.surface_name.casefold()
        ),
    ))
    records = json.loads(FIXTURE.read_text())["surfaces"]
    keys = {
        "{}/{}".format(source.system_name, source.surface_name)
        for source in sources
    }
    if keys != set(records):
        raise ValueError("README catalog and numerical fixture differ")
    systems = defaultdict(list)
    for source in sources:
        systems[source.system_name].append(source)
    lines = [
        "The catalog currently exposes **{} callable surfaces** across "
        "**{} chemical systems**. The final column gives the number of "
        "electronic surfaces represented by each callable PES or DPEM entry.".format(
            len(sources), len(systems)
        ),
        "",
        "<table>",
        "<thead>",
        "<tr><th>System composition</th><th>Surface name</th>"
        "<th>Surface description</th><th>Electronic surfaces</th></tr>",
        "</thead>",
        "<tbody>",
    ]
    for system_name in sorted(systems, key=str.casefold):
        system_sources = systems[system_name]
        for source_index, source in enumerate(system_sources):
            key = "{}/{}".format(source.system_name, source.surface_name)
            cells = []
            if source_index == 0:
                cells.append(
                    '<td rowspan="{}" valign="top"><strong>{}</strong></td>'.format(
                        len(system_sources), chemical_composition(system_name)
                    )
                )
            cells.append("<td><code>{}</code></td>".format(source.surface_name))
            cells.append("<td>{}</td>".format(surface_description(source)))
            cells.append(
                '<td align="center">{}</td>'.format(records[key]["nstates"])
            )
            lines.append("<tr>{}</tr>".format("".join(cells)))
    lines.extend(("</tbody>", "</table>"))
    return "\n".join(lines)


def updated_readme(original: str) -> str:
    """Return README text with the generated catalog replaced."""

    if original.count(START_MARKER) != 1 or original.count(END_MARKER) != 1:
        raise ValueError("README must contain one PES catalog marker pair")
    start = original.index(START_MARKER)+len(START_MARKER)
    end = original.index(END_MARKER)
    if start > end:
        raise ValueError("README PES catalog markers are reversed")
    return (
        original[:start]
        + "\n\n"
        + generated_catalog()
        + "\n\n"
        + original[end:]
    )


def main() -> int:
    """Write or verify the README PES catalog."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    original = README.read_text()
    updated = updated_readme(original)
    if arguments.check:
        if original != updated:
            print("outdated: README.md PES catalog")
            return 1
    else:
        README.write_text(updated)
        print("README catalog surfaces: {}".format(
            len(json.loads(FIXTURE.read_text())["surfaces"])
        ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
