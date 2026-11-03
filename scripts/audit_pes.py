#!/usr/bin/env python3
"""Inventory CHEMPOTPY systems, sources, duplicates, and repository artifacts."""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from chempotpy.core.registry import (  # noqa: E402
    discover_surface_sources,
    duplicate_surface_sources,
)


def audit() -> Dict[str, object]:
    """Collect a deterministic, machine-readable migration baseline."""

    package_root = REPOSITORY_ROOT / "chempotpy"
    sources = discover_surface_sources(package_root)
    duplicates = duplicate_surface_sources(sources)
    python_errors = []
    python_count = 0
    for path in sorted(package_root.rglob("*.py")):
        python_count += 1
        try:
            ast.parse(path.read_text(errors="replace"), filename=str(path))
        except SyntaxError as error:
            python_errors.append(
                {"path": str(path.relative_to(REPOSITORY_ROOT)),
                 "line": error.lineno,
                 "message": error.msg}
            )

    tracked_artifacts = subprocess.run(
        ["git", "ls-files"],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    tracked_artifacts = [
        path for path in tracked_artifacts
        if Path(path).suffix
        in {".log", ".mod", ".o", ".pyd", ".smod", ".so", ".swp", ".swo"}
    ]
    present_tracked_artifacts = [
        path for path in tracked_artifacts
        if (REPOSITORY_ROOT / path).exists()
    ]

    return {
        "system_count": len({source.system_name for source in sources}),
        "fortran_source_count": len(sources),
        "surface_count": sum(source.is_surface for source in sources),
        "adiabatic_surface_count": sum(
            source.entry_point == "pes" for source in sources
        ),
        "diabatic_surface_count": sum(
            source.entry_point == "dpem" for source in sources
        ),
        "support_source_count": sum(
            not source.is_surface for source in sources
        ),
        "fixed_form_count": sum(s.source_form == "fixed" for s in sources),
        "free_form_count": sum(s.source_form == "free" for s in sources),
        "duplicate_surface_count": len(duplicates),
        "duplicate_surfaces": {
            f"{system}/{surface}": [
                str(item.path.relative_to(REPOSITORY_ROOT)) for item in items
            ]
            for (system, surface), items in sorted(duplicates.items())
        },
        "python_file_count": python_count,
        "python_syntax_errors": python_errors,
        "tracked_generated_artifacts": tracked_artifacts,
        "tracked_generated_artifacts_present": present_tracked_artifacts,
    }


def main() -> int:
    """Print the audit as readable text or JSON."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    report = audit()
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        tracked_count = len(report["tracked_generated_artifacts"])
        present_count = len(
            report["tracked_generated_artifacts_present"]
        )
        for key in (
            "system_count", "fortran_source_count", "surface_count",
            "adiabatic_surface_count", "diabatic_surface_count",
            "support_source_count", "fixed_form_count", "free_form_count",
            "duplicate_surface_count", "python_file_count",
        ):
            print(f"{key}: {report[key]}")
        print(f"python_syntax_errors: {len(report['python_syntax_errors'])}")
        print(
            "tracked_generated_artifacts: "
            f"{tracked_count}"
        )
        print(
            "tracked_generated_artifacts_present: "
            f"{present_count}"
        )
        print(
            "tracked_generated_artifacts_pending_removal: "
            f"{tracked_count - present_count}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
