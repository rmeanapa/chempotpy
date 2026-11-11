"""Tests for the discovered isolated f2py build commands."""

import unittest
import tempfile
from pathlib import Path

from chempotpy.core.registry import SurfaceSource
from scripts.build_pes import (
    CORE_FORTRAN_ROOT,
    copy_ascii_source,
    core_fortran_dependencies,
    extract_entry_point_source,
    existing_result,
    f2py_command,
    normalized_sources,
    write_f2py_kind_map,
)
from scripts.check_fortran import compiler_command


class BuildPESTests( unittest.TestCase ):
    def test_core_fortran_sources_live_directly_under_core( self ) -> None:
        self.assertEqual(CORE_FORTRAN_ROOT.name, "core")
        self.assertTrue((CORE_FORTRAN_ROOT / "chempotpy_kinds.f90").is_file())
        self.assertFalse((CORE_FORTRAN_ROOT / "fortran").exists())

    def test_entry_point_extraction_omits_internal_modules( self ) -> None:
        source = """\
module owner
contains
    subroutine evaluate( x )
        real(8), intent(in) :: x
    end subroutine evaluate
end module owner

subroutine pes( x )
    use owner, only: evaluate
    real(8), intent(in) :: x
    call evaluate(x)
end subroutine pes
"""
        interface = extract_entry_point_source(source, "pes")
        self.assertTrue(interface.startswith("subroutine pes"))
        self.assertIn("call evaluate(x)", interface)
        self.assertNotIn("module owner", interface)

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
        self.assertTrue(any("-O0" in argument for argument in command))
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
        self.assertTrue(any("-fcheck=all" in argument for argument in command))

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

    def test_neural_network_parameters_precede_their_consumer( self ) -> None:
        source_path = (
            Path(__file__).resolve().parents[1]
            / "chempotpy" / "N4" / "N4_1A_PIPNN.f90"
        )
        source = SurfaceSource(
            system_name="N4",
            surface_name="N4_1A_PIPNN",
            path=source_path,
            source_form="free",
            entry_point="pes",
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            normalized = normalized_sources(
                source, Path(temporary_directory)
            )
        self.assertEqual(
            [path.name for path in normalized[-2:]],
            [
                "chempotpy_n4_pipnn_parameters.f90",
                "N4_1A_PIPNN.f90",
            ],
        )
        command = compiler_command(source, Path("modules"), "gfortran")
        support_index = command.index(str(source_path.parent / (
            "chempotpy_n4_pipnn_parameters.f90"
        )))
        self.assertLess(support_index, command.index(str(source_path)))

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
