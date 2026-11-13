#!/usr/bin/env python3
"""Reject numerical changes to any released CHEMPOTPY surface."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import discover_surface_sources
from scripts.capture_library_references import capture_source

DEFAULT_FIXTURE = REPOSITORY_ROOT / "tests" / "data" / \
    "all_surface_references.json"
METADATA_FIELDS = ("entry_point", "natoms", "nstates", "signature")
RELATIVE_TOLERANCE = 1.0e-13
ABSOLUTE_TOLERANCE = 1.0e-13
BASELINE_VERSION = 1
MAX_WORKFLOW_ANNOTATIONS = 50


def escape_workflow_command(value: str) -> str:
    """Escape text for the data field of a GitHub workflow command."""

    return value.replace("%", "%25").replace("\r", "%0D").replace(
        "\n", "%0A"
    )


def emit_workflow_annotations(diagnostics: List[str]) -> None:
    """Expose regression details as GitHub check annotations when in Actions."""

    if os.environ.get("GITHUB_ACTIONS") != "true":
        return
    for diagnostic in diagnostics[:MAX_WORKFLOW_ANNOTATIONS]:
        print(
            "::error title=CHEMPOTPY numerical regression::{}".format(
                escape_workflow_command(diagnostic)
            )
        )
    omitted = len(diagnostics) - MAX_WORKFLOW_ANNOTATIONS
    if omitted > 0:
        print(
            "::error title=CHEMPOTPY numerical regression::{}".format(
                "{} additional diagnostics omitted".format(omitted)
            )
        )


def load_fixture(path: Path) -> Tuple[Optional[int], Dict[str, dict]]:
    """Load the baseline version and surface mapping from one fixture."""

    try:
        document = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("cannot read {}: {}".format(path, error)) from error
    surfaces = document.get("surfaces")
    if not isinstance(surfaces, dict) or not surfaces:
        raise ValueError("{} has no non-empty 'surfaces' mapping".format(path))
    version = document.get("baseline_version")
    if version is not None and not isinstance(version, int):
        raise ValueError("{} has an invalid baseline_version".format(path))
    return version, surfaces


def select_baseline(
    baseline_version: Optional[int],
    baseline: Dict[str, dict],
    current_version: Optional[int],
    current: Dict[str, dict],
) -> Tuple[Dict[str, dict], bool]:
    """Select the immutable baseline, allowing only the initial bootstrap."""

    if current_version != BASELINE_VERSION:
        raise ValueError(
            "current fixture must declare baseline_version {}".format(
                BASELINE_VERSION
            )
        )
    if baseline_version is None:
        return current, True
    if baseline_version != BASELINE_VERSION:
        raise ValueError(
            "base fixture uses unsupported baseline_version {}".format(
                baseline_version
            )
        )
    return baseline, False


def compare_records(
    key: str,
    actual: dict,
    expected: dict,
) -> Tuple[List[str], int]:
    """Return diagnostics and the number of numerical values compared."""

    diagnostics = []
    for field in METADATA_FIELDS:
        if actual.get(field) != expected.get(field):
            diagnostics.append(
                "{}: {} changed from {!r} to {!r}".format(
                    key, field, expected.get(field), actual.get(field)
                )
            )
    expected_values = expected.get("values")
    actual_values = actual.get("values")
    if not isinstance(expected_values, list) or not isinstance(
        actual_values, list
    ):
        diagnostics.append("{}: reference values are not arrays".format(key))
        return diagnostics, 0
    if len(actual_values) != len(expected_values):
        diagnostics.append(
            "{}: value count changed from {} to {}".format(
                key, len(expected_values), len(actual_values)
            )
        )
        return diagnostics, min(len(expected_values), len(actual_values))
    for index, (expected_text, actual_text) in enumerate(
        zip(expected_values, actual_values)
    ):
        try:
            expected_number = float(expected_text)
            actual_number = float(actual_text)
        except (TypeError, ValueError):
            diagnostics.append(
                "{}: value {} is not numeric: expected {!r}, got {!r}".format(
                    key, index, expected_text, actual_text
                )
            )
            continue
        if not math.isfinite(expected_number) or not math.isfinite(
            actual_number
        ):
            diagnostics.append(
                "{}: value {} is non-finite: expected {!r}, got {!r}".format(
                    key, index, expected_text, actual_text
                )
            )
            continue
        if expected_text == actual_text:
            continue
        if (
            abs(expected_number) < 1.0e-300
            and abs(actual_number) < 1.0e-300
        ):
            continue
        if not math.isclose(
            actual_number,
            expected_number,
            rel_tol=RELATIVE_TOLERANCE,
            abs_tol=ABSOLUTE_TOLERANCE,
        ):
            diagnostics.append(
                "{}: value {} changed from {} to {}".format(
                    key, index, expected_text, actual_text
                )
            )
    return diagnostics, len(expected_values)


def comparison_record(
    key: str,
    actual: dict,
    baseline: Dict[str, dict],
    current_fixture: Dict[str, dict],
) -> dict:
    """Return the protected record, bootstrapping values for a new PES."""

    if key in baseline:
        return baseline[key]
    expected = dict(current_fixture[key])
    expected["values"] = actual.get("values")
    return expected


def main() -> int:
    """Compile every current PES and compare it with an immutable baseline."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        required=True,
        help="reference fixture extracted from the base Git revision",
    )
    parser.add_argument(
        "--current-fixture",
        type=Path,
        default=DEFAULT_FIXTURE,
        help="current fixture, used only for newly added surfaces",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(4, os.cpu_count() or 1),
    )
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    compiler = shutil.which("gfortran")
    if compiler is None:
        parser.error("gfortran is required")
    try:
        baseline_version, base_fixture = load_fixture(args.baseline.resolve())
        current_version, current_fixture = load_fixture(
            args.current_fixture.resolve()
        )
        baseline, bootstrapped = select_baseline(
            baseline_version,
            base_fixture,
            current_version,
            current_fixture,
        )
    except ValueError as error:
        parser.error(str(error))

    sources = {
        "{}/{}".format(item.system_name, item.surface_name): item.path
        for item in discover_surface_sources()
        if item.is_surface
    }
    diagnostics = []
    missing_current_records = sorted(set(sources) - set(current_fixture))
    stale_current_records = sorted(set(current_fixture) - set(sources))
    removed_baseline_surfaces = sorted(set(baseline) - set(sources))
    for key in missing_current_records:
        diagnostics.append("{}: missing from the current fixture".format(key))
    for key in stale_current_records:
        diagnostics.append("{}: current fixture has no matching PES".format(key))
    for key in removed_baseline_surfaces:
        diagnostics.append("{}: released PES was removed".format(key))
    if diagnostics:
        for diagnostic in diagnostics:
            print("FAIL {}".format(diagnostic))
        emit_workflow_annotations(diagnostics)
        return 1

    ordered_sources = [sources[key] for key in sorted(sources)]
    with tempfile.TemporaryDirectory(
        prefix="chempotpy-numerical-regression-"
    ) as tmp:
        work_root = Path(tmp)
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.jobs
        ) as executor:
            results = list(executor.map(
                lambda source: capture_source(compiler, source, work_root),
                ordered_sources,
            ))

    numerical_values = 0
    new_surfaces = 0
    for result in results:
        if result.error:
            diagnostics.append(result.error)
            continue
        expected = comparison_record(
            result.key, result.record, baseline, current_fixture
        )
        if result.key not in baseline:
            new_surfaces += 1
        record_diagnostics, compared = compare_records(
            result.key, result.record, expected
        )
        diagnostics.extend(record_diagnostics)
        numerical_values += compared

    print("surfaces_checked: {}".format(len(results)))
    print("baseline_version: {}".format(BASELINE_VERSION))
    print("baseline_bootstrap: {}".format("yes" if bootstrapped else "no"))
    print("baseline_surfaces: {}".format(len(baseline)))
    print("new_surfaces: {}".format(new_surfaces))
    print("numerical_values_checked: {}".format(numerical_values))
    print("numerical_regressions: {}".format(len(diagnostics)))
    for diagnostic in diagnostics[:100]:
        print("FAIL {}".format(diagnostic))
    if len(diagnostics) > 100:
        print("FAIL ... {} additional diagnostics omitted".format(
            len(diagnostics) - 100
        ))
    emit_workflow_annotations(diagnostics)
    return 1 if diagnostics else 0


if __name__ == "__main__":
    raise SystemExit(main())
