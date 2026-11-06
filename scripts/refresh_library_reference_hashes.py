#!/usr/bin/env python3
"""Refresh source hashes after proving captured PES outputs are unchanged."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


REFERENCE_RELATIVE_TOLERANCE = 1.0e-13
REFERENCE_ABSOLUTE_TOLERANCE = 1.0e-13


def equivalent_values(expected: list[str], actual: list[str]) -> bool:
    """Return whether captured values satisfy the regression tolerance."""

    if len(expected) != len(actual):
        return False
    return all(
        expected_value == actual_value
        or math.isclose(
            float(expected_value),
            float(actual_value),
            rel_tol=REFERENCE_RELATIVE_TOLERANCE,
            abs_tol=REFERENCE_ABSOLUTE_TOLERANCE,
        )
        for expected_value, actual_value in zip(expected, actual)
    )


def main() -> int:
    """Merge verified source hashes while preserving fixture provenance."""

    parser = argparse.ArgumentParser()
    parser.add_argument("fixture", type=Path)
    parser.add_argument("capture", type=Path)
    arguments = parser.parse_args()
    fixture_document = json.loads(arguments.fixture.read_text())
    capture_document = json.loads(arguments.capture.read_text())
    fixture = fixture_document["surfaces"]
    capture = capture_document["surfaces"]
    if not set(capture).issubset(fixture):
        raise ValueError("capture contains surfaces absent from the fixture")
    compared_fields = ("entry_point", "natoms", "nstates", "signature")
    for key in sorted(capture):
        for field in compared_fields:
            if fixture[key][field] != capture[key][field]:
                raise ValueError("{} changed for {}".format(field, key))
        if not equivalent_values(fixture[key]["values"], capture[key]["values"]):
            raise ValueError("values changed for {}".format(key))
        fixture[key]["source_sha256"] = capture[key]["source_sha256"]
    arguments.fixture.write_text(
        json.dumps(fixture_document, indent=2, sort_keys=True) + "\n"
    )
    print("source_hashes_refreshed: {}".format(len(capture)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
