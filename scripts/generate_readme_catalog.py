#!/usr/bin/env python3
"""Generate the system-level PES summary table in the project README."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import sys
from typing import List


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import SurfaceSource, discover_surface_sources


README = REPOSITORY_ROOT / "README.md"
FIXTURE = REPOSITORY_ROOT / "tests" / "data" / "all_surface_references.json"
START_MARKER = "<!-- BEGIN GENERATED PES CATALOG -->"
END_MARKER = "<!-- END GENERATED PES CATALOG -->"
MODEL_FAMILIES = (
    (re.compile(r"PIPNN", re.IGNORECASE), "PIP neural-network"),
    (re.compile(r"APRP", re.IGNORECASE), "anchor-points reactive potential"),
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

SYSTEM_DESCRIPTIONS = {
    "ALC": "Extended-LEPS model surfaces for atom–linear-molecule collisions.",
    "AlmHn": "Valence-bond-order models for aluminum–hydrogen clusters.",
    "Aln": "Aluminum-cluster embedded-atom, extended-Rydberg, and many-body models.",
    "H3": "Benchmark analytic H₃ surfaces, including BKMP, DMBE, LSTH, and related fits.",
    "N2O": "Ground- and excited-state N₂O fits, including PIP and ZV representations.",
    "N4": "PIP, neural-network, many-body-expansion, and global analytic N₄ surfaces.",
    "O3": "Single- and multistate ozone surfaces in analytic, many-body, and diabatic forms.",
    "O4": "Singlet, triplet, and quintet O₄ many-body and analytic surfaces.",
}


def joined_phrases(phrases: List[str]) -> str:
    """Join short descriptive phrases with natural English punctuation."""

    if len(phrases) == 1:
        return phrases[0]
    if len(phrases) == 2:
        return "{} and {}".format(*phrases)
    return "{}, and {}".format(", ".join(phrases[:-1]), phrases[-1])


def initial_capital(text: str) -> str:
    """Capitalize a sentence without lowercasing scientific acronyms."""

    return text[:1].upper()+text[1:]


def system_description(system_name: str, sources: List[SurfaceSource]) -> str:
    """Return a concise description of one system's surface families."""

    if system_name in SYSTEM_DESCRIPTIONS:
        return SYSTEM_DESCRIPTIONS[system_name]
    families = []
    for source in sources:
        for pattern, description in MODEL_FAMILIES:
            if pattern.search(source.surface_name):
                if description not in families:
                    families.append(description)
                break
        else:
            if "analytic" not in families:
                families.append("analytic")
    description = "{} {}.".format(
        initial_capital(joined_phrases(families)),
        "surface" if len(sources) == 1 else "surfaces",
    )
    if any(source.entry_point == "dpem" for source in sources):
        description = description[:-1] + " with paired adiabatic and diabatic forms."
    return description


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
        "**{} chemical systems**. Paired PES and DPEM entry points are counted "
        "separately because each is independently callable.".format(
            len(sources), len(systems)
        ),
        "",
        "| System | Surface description | Number of surfaces |",
        "| :--- | :--- | ---: |",
    ]
    for system_name in sorted(systems, key=str.casefold):
        system_sources = systems[system_name]
        lines.append(
            "| `{}` | {} | {} |".format(
                system_name,
                system_description(system_name, system_sources),
                len(system_sources),
            )
        )
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
