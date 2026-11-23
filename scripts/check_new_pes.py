#!/usr/bin/env python3
"""Validate every PES source added relative to a selected Git revision."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import ENTRY_POINT_PATTERN, SurfaceSource
from scripts.capture_library_references import (
    BASELINE_VERSION,
    REGRESSION_GEOMETRY_COUNT,
)
from scripts.check_fortran_style import check_source
from scripts.check_modernization import active_code
from scripts.export_standalone import (
    combined_source,
    example_driver,
    surface_interface,
)


DEFAULT_FIXTURE = REPOSITORY_ROOT / "tests" / "data" / \
    "all_surface_references.json"
REQUIRED_METADATA = (
    ("System", r"System"),
    ("Functional form", r"Functional\s+form"),
    ("Common name", r"Common\s+name"),
    ("Reference", r"References?"),
    ("Number of derivatives", r"Number\s+of\s+derivatives"),
    ("Number of bodies", r"Number\s+of\s+bodies"),
    (
        "Number of electronic surfaces",
        r"Number\s+of\s+electronic\s+(?:surfaces|states)",
    ),
    ("Interface", r"Interface"),
    ("Data file", r"Data\s+file"),
    ("Notes", r"Notes"),
)
PLACEHOLDER_PATTERN = re.compile(
    r"CHEMPOTPY_IMPLEMENTATION_REQUIRED|{{[A-Z0-9_]+}}|\b(?:TODO|TBD)\b",
    re.IGNORECASE,
)
FORBIDDEN_SOURCE_PATTERNS = (
    ("external file OPEN", re.compile(r"\bopen\s*\(", re.IGNORECASE)),
    ("external file READ", re.compile(r"\bread\s*\(", re.IGNORECASE)),
    ("COMMON block", re.compile(r"^\s*common\b", re.IGNORECASE)),
    ("DATA statement", re.compile(r"^\s*data\b", re.IGNORECASE)),
    ("ENTRY statement", re.compile(r"^\s*entry\b", re.IGNORECASE)),
    ("saved state", re.compile(r"\bsave\b", re.IGNORECASE)),
    ("GOTO", re.compile(r"\b(?:go\s*to|goto)\b", re.IGNORECASE)),
    (
        "local neural-network activation kernel",
        re.compile(r"^\s*function\s+tranfun\d*\s*\(", re.IGNORECASE),
    ),
    (
        "local D3 dispersion kernel",
        re.compile(r"^\s*subroutine\s+edisp\s*\(", re.IGNORECASE),
    ),
    (
        "local D3 element table",
        re.compile(r"r2r4\s*\(\s*1\s*:\s*94", re.IGNORECASE),
    ),
    (
        "local matrix multiplication kernel",
        re.compile(r"^\s*subroutine\s+mult\s*\(", re.IGNORECASE),
    ),
    (
        "local vector-operation kernel",
        re.compile(
            r"^\s*(?:(?:real\([^)]*\)\s+)?function|subroutine)\s+"
            r"(?:anorm|crosprod|dotprod|sprod|xprod|xprod1)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "local diagonalization kernel",
        re.compile(
            r"^\s*(?:subroutine|call)\s+(?:diagonalize|eign|dsyev)\b",
            re.IGNORECASE,
        ),
    ),
    (
        "local spline kernel",
        re.compile(
            r"^\s*(?:subroutine|call)\s+(?:spline|splint|splnin)\b",
            re.IGNORECASE,
        ),
    ),
)


def metadata_values(source_text: str) -> Tuple[Dict[str, str], List[str]]:
    """Return standardized metadata values and any structural violations."""

    values: Dict[str, str] = {}
    violations = []
    lines = source_text.splitlines()
    first_content = next((line for line in lines if line.strip()), "")
    if re.match(r"^\s*!\s*System\s*:", first_content, re.I) is None:
        violations.append("the metadata block must begin on the first line")
    comment_prefix = []
    for line in lines:
        if not line.strip() or line.lstrip().startswith("!"):
            comment_prefix.append(line)
            continue
        break
    metadata_text = "\n".join(comment_prefix)
    for label, expression in REQUIRED_METADATA:
        match = re.search(
            r"^\s*!\s*{}\s*:\s*(.*?)\s*$".format(expression),
            metadata_text,
            re.IGNORECASE | re.MULTILINE,
        )
        if match is None:
            violations.append("missing '{}' metadata".format(label))
        elif not match.group(1):
            violations.append("'{}' metadata must not be empty".format(label))
        else:
            values[label] = match.group(1)
    return values, violations


def source_contract_violations(path: Path) -> List[str]:
    """Return onboarding violations found directly in one PES source."""

    source_text = path.read_text(errors="replace")
    metadata, violations = metadata_values(source_text)
    if PLACEHOLDER_PATTERN.search(source_text):
        violations.append("remove every scaffold placeholder before review")
    if metadata.get("System") != path.parent.name:
        violations.append(
            "System metadata must match directory '{}'".format(
                path.parent.name
            )
        )
    if "use chempotpy_kinds" not in source_text.lower():
        violations.append("use chempotpy_kinds and real(wp)")
    if re.search(r"^\s*implicit\s+none\b", source_text, re.I | re.M) is None:
        violations.append("declare IMPLICIT NONE")
    if re.search(r"^\s*private\s*$", source_text, re.I | re.M) is None:
        violations.append("make the owning module private by default")
    if re.search(
        r"^\s*public\s*::[^\n]*\bevaluate\b",
        source_text,
        re.IGNORECASE | re.MULTILINE,
    ) is None:
        violations.append("expose only the module evaluator explicitly")
    entry_points = ENTRY_POINT_PATTERN.findall(source_text)
    if len(entry_points) != 1:
        violations.append("define exactly one public PES or DPEM entry point")
    for line_number, line in enumerate(source_text.splitlines(), start=1):
        code = active_code(line)
        for label, pattern in FORBIDDEN_SOURCE_PATTERNS:
            if pattern.search(code):
                violations.append("line {}: {} is forbidden".format(
                    line_number, label
                ))
    return violations


def fixture_violations(
    source: SurfaceSource,
    fixture_document: dict,
) -> List[str]:
    """Return fixture completeness and provenance violations for one PES."""

    key = "{}/{}".format(source.system_name, source.surface_name)
    violations = []
    if fixture_document.get("baseline_version") != BASELINE_VERSION:
        violations.append(
            "the fixture baseline_version must be {}".format(BASELINE_VERSION)
        )
    record = fixture_document.get("surfaces", {}).get(key)
    if not isinstance(record, dict):
        return violations + ["add an independently verified fixture record"]
    entry_point, signature, natoms, nstates = surface_interface(source)
    expected = {
        "entry_point": entry_point,
        "signature": signature,
        "natoms": natoms,
        "nstates": nstates,
        "source_sha256": hashlib.sha256(source.path.read_bytes()).hexdigest(),
    }
    for field, value in expected.items():
        if record.get(field) != value:
            violations.append(
                "fixture field '{}' must equal {!r}".format(field, value)
            )
    values = record.get("values")
    expected_count = REGRESSION_GEOMETRY_COUNT * (
        nstates if entry_point == "pes" else nstates*nstates
    )
    if not isinstance(values, list) or len(values) != expected_count:
        violations.append(
            "fixture must contain {} reference values from {} geometries".format(
                expected_count, REGRESSION_GEOMETRY_COUNT
            )
        )
    else:
        try:
            finite = all(math.isfinite(float(value)) for value in values)
        except (TypeError, ValueError):
            finite = False
        if not finite:
            violations.append("fixture reference values must all be finite")
    return violations


def surface_from_path(path: Path) -> SurfaceSource:
    """Construct the registry record for one newly added PES path."""

    source_text = path.read_text(errors="replace")
    match = ENTRY_POINT_PATTERN.search(source_text)
    return SurfaceSource(
        system_name=path.parent.name,
        surface_name=path.stem,
        path=path,
        source_form="free",
        entry_point=match.group(1).lower() if match else None,
    )


def added_surface_paths(base_ref: str) -> Tuple[Path, ...]:
    """Return callable Fortran sources added since one Git base revision."""

    completed = subprocess.run(
        [
            "git", "diff", "--diff-filter=A", "--name-only",
            "{}...HEAD".format(base_ref), "--", "chempotpy",
        ],
        cwd=str(REPOSITORY_ROOT),
        check=True,
        capture_output=True,
        text=True,
    )
    paths = []
    for relative in completed.stdout.splitlines():
        path = REPOSITORY_ROOT / relative
        if path.suffix.lower() == ".f90" and path.is_file():
            if ENTRY_POINT_PATTERN.search(path.read_text(errors="replace")):
                paths.append(path)
    return tuple(sorted(paths))


def standalone_validation_error(source: SurfaceSource) -> str:
    """Compile and run one dependency-complete standalone PES example."""

    compiler = shutil.which("gfortran")
    if compiler is None:
        return "gfortran is required to validate the standalone export"
    try:
        source_text = combined_source(source, include_validation=True)
        driver_text = example_driver(
            source, geometry_count=1, difference_step=1.0e-3
        )
        with tempfile.TemporaryDirectory(
            prefix="chempotpy-new-pes-standalone-"
        ) as temporary_directory:
            directory = Path(temporary_directory)
            source_path = directory / "surface.f90"
            driver_path = directory / "example.f90"
            executable = directory / "example"
            source_path.write_text(source_text)
            driver_path.write_text(driver_text)
            subprocess.run(
                [
                    compiler,
                    "-ffree-line-length-none",
                    str(source_path),
                    str(driver_path),
                    "-o",
                    str(executable),
                ],
                cwd=str(directory),
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [str(executable)],
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
        detail = getattr(error, "stderr", "") or str(error)
        return "standalone compile/run failed: {}".format(detail.strip())
    return ""


def validate_surface(path: Path, fixture_document: dict) -> List[str]:
    """Return all onboarding violations for one newly added PES."""

    violations = source_contract_violations(path)
    initialization = path.parent / "__init__.py"
    if not initialization.is_file():
        violations.append("the system directory must contain __init__.py")
    elif path.stem not in initialization.read_text(errors="replace"):
        violations.append(
            "add the surface to the system's Python descriptions"
        )
    source = surface_from_path(path)
    if source.entry_point is not None:
        metadata, _ = metadata_values(path.read_text(errors="replace"))
        entry_point, _, _, nstates = surface_interface(source)
        interface = metadata.get("Interface", "").lower()
        if interface != entry_point:
            violations.append(
                "Interface metadata must equal '{}'".format(entry_point.upper())
            )
        state_match = re.match(
            r"\s*(\d+)", metadata.get("Number of electronic surfaces", "")
        )
        if state_match is None or int(state_match.group(1)) != nstates:
            violations.append(
                "electronic-surface metadata must equal nstates ({})".format(
                    nstates
                )
            )
        violations.extend(fixture_violations(source, fixture_document))
        try:
            standalone_error = standalone_validation_error(source)
            if standalone_error:
                violations.append(standalone_error)
        except (OSError, ValueError) as error:
            violations.append("standalone export failed: {}".format(error))
    for violation in check_source(path):
        violations.append(
            "style {} at line {}: {}".format(
                violation.rule, violation.line, violation.message
            )
        )
    return violations


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Validate explicitly selected or newly added PES sources."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--base-ref")
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    args = parser.parse_args(argv)
    if bool(args.paths) == bool(args.base_ref):
        parser.error("provide paths or --base-ref, but not both")
    try:
        fixture_document = json.loads(args.fixture.read_text())
        paths = (
            tuple(path.resolve() for path in args.paths)
            if args.paths
            else added_surface_paths(args.base_ref)
        )
    except (OSError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    diagnostics = []
    for path in paths:
        try:
            relative_path = path.relative_to(REPOSITORY_ROOT)
        except ValueError:
            diagnostics.append("{}: source is outside the repository".format(path))
            continue
        for violation in validate_surface(path, fixture_document):
            diagnostics.append(
                "{}: {}".format(relative_path, violation)
            )
    if diagnostics:
        print("new PES onboarding violations:")
        for diagnostic in diagnostics:
            print("  {}".format(diagnostic))
        print("new_pes_checked: {}".format(len(paths)))
        print("new_pes_violations: {}".format(len(diagnostics)))
        return 1
    print("new_pes_checked: {}".format(len(paths)))
    print("new_pes_violations: 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
