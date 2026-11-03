#!/usr/bin/env python3
"""Refresh source hashes after proving captured PES outputs are unchanged."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


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
    compared_fields = ("entry_point", "natoms", "nstates", "signature", "values")
    for key in sorted(capture):
        for field in compared_fields:
            if fixture[key][field] != capture[key][field]:
                raise ValueError("{} changed for {}".format(field, key))
        fixture[key]["source_sha256"] = capture[key]["source_sha256"]
    arguments.fixture.write_text(
        json.dumps(fixture_document, indent=2, sort_keys=True) + "\n"
    )
    print("source_hashes_refreshed: {}".format(len(capture)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
