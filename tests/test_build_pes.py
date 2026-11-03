"""Tests for the discovered isolated f2py build commands."""

import unittest
import tempfile
from pathlib import Path

from chempotpy.core.registry import SurfaceSource
from scripts.build_pes import (
    copy_ascii_source,
    core_fortran_dependencies,
    existing_result,
    f2py_command,
    normalized_sources,
    write_f2py_kind_map,
)


class BuildPESTests( unittest.TestCase ):
    def test_build_command_uses_selected_backend( self ) -> None:
        source = SurfaceSource(
            system_name="H3",
            surface_name="example",
            path=Path(__file__),
            source_form="free",
            entry_point="pes",
        )
        command = f2py_command(
            source,
            (Path("example.f90"),),
            "python3",
            "meson",
        )
        self.assertIn("meson", command)
        self.assertIn("example.f90", command)
        self.assertNotIn("-llapack", command)
        self.assertEqual(command[-3:], ["only:", "pes", ":"])

    def test_debug_build_requests_debug_backend_profile( self ) -> None:
        source = SurfaceSource(
            system_name="H3",
            surface_name="example",
            path=Path(__file__),
            source_form="free",
            entry_point="pes",
        )
        command = f2py_command(
            source, (Path("example.f90"),), "python3", "meson", debug=True
        )
        self.assertIn("--debug", command)

    def test_f2py_kind_map_keeps_real64_at_double_precision( self ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            mapping = write_f2py_kind_map(Path(temporary_directory))
            content = mapping.read_text()
        self.assertIn("'real32': 'float'", content)
        self.assertIn("'real64': 'double'", content)

    def test_source_normalization_preserves_fixed_form_columns( self ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = Path(temporary_directory) / "source.f"
            destination = Path(temporary_directory) / "destination.f"
            content = b"      x = 1\nC author: \xc3\xa9\n"
            source.write_bytes(content)
            copy_ascii_source(source, destination)
            normalized = destination.read_bytes()
        self.assertEqual(len(normalized), len(content))
        self.assertTrue(normalized.isascii())
        self.assertEqual(normalized.splitlines()[0], b"      x = 1")

    def test_core_dependencies_are_resolved_before_their_consumer( self ) -> None:
        dependencies = core_fortran_dependencies(
            "use chempotpy_geometry, only: cartesian_distance\n"
        )
        self.assertEqual(
            [path.name for path in dependencies],
            ["chempotpy_kinds.f90", "chempotpy_geometry.f90"],
        )

    def test_normalization_bundles_only_referenced_core_modules( self ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source_path = root / "example.f90"
            source_path.write_text(
                "use chempotpy_legacy_units, "
                "only: legacy_bohr_to_angstrom\n"
                "subroutine pes()\nend subroutine pes\n"
            )
            source = SurfaceSource(
                system_name="H3",
                surface_name="example",
                path=source_path,
                source_form="free",
                entry_point="pes",
            )
            build_directory = root / "build"
            build_directory.mkdir()
            normalized = normalized_sources(source, build_directory)
        self.assertEqual(
            [path.name for path in normalized],
            [
                "chempotpy_kinds.f90",
                "chempotpy_legacy_units.f90",
                "example.f90",
            ],
        )

    def test_resume_recognizes_one_existing_extension( self ) -> None:
        source = SurfaceSource(
            system_name="H3",
            surface_name="example",
            path=Path(__file__),
            source_form="free",
            entry_point="pes",
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            output_directory = root / "output" / "H3"
            output_directory.mkdir(parents=True)
            artifact = output_directory / "example.cpython-test.so"
            artifact.touch()
            result = existing_result(
                source, root / "output", root / "logs"
            )
        self.assertIsNotNone(result)
        self.assertTrue(result.passed)
        self.assertTrue(result.skipped)

    def test_resume_does_not_confuse_surface_name_prefixes( self ) -> None:
        source = SurfaceSource(
            system_name="NaFH",
            surface_name="example",
            path=Path(__file__),
            source_form="free",
            entry_point="pes",
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            output_directory = root / "output" / "NaFH"
            output_directory.mkdir(parents=True)
            (output_directory / "example_DPEM.cpython-test.so").touch()
            result = existing_result(
                source, root / "output", root / "logs"
            )
        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
