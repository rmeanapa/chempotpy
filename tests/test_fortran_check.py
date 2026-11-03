"""Tests for deterministic isolated Fortran syntax commands."""

import unittest
import tempfile
from pathlib import Path

from chempotpy.core.registry import SurfaceSource
from scripts.check_fortran import compiler_command


class FortranCheckTests( unittest.TestCase ):
    def test_fixed_form_command( self ) -> None:
        source = SurfaceSource(
            system_name="H3",
            surface_name="example",
            path=Path("example.f"),
            source_form="fixed",
            entry_point="pes",
        )
        command = compiler_command(source, Path("modules"), "gfortran")
        self.assertIn("-ffixed-form", command)
        self.assertIn("-ffixed-line-length-72", command)
        self.assertNotIn("-ffree-line-length-none", command)

    def test_free_form_command( self ) -> None:
        source = SurfaceSource(
            system_name="H3",
            surface_name="example",
            path=Path("example.f90"),
            source_form="free",
            entry_point="pes",
        )
        command = compiler_command(source, Path("modules"), "gfortran")
        self.assertIn("-ffree-form", command)
        self.assertIn("-ffree-line-length-none", command)
        self.assertNotIn("-ffixed-line-length-none", command)

    def test_free_form_command_includes_core_dependencies( self ) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source_path = Path(temporary_directory) / "example.f90"
            source_path.write_text(
                "use chempotpy_geometry, only: cartesian_distance\n"
            )
            source = SurfaceSource(
                system_name="H3",
                surface_name="example",
                path=source_path,
                source_form="free",
                entry_point="pes",
            )
            command = compiler_command(source, Path("modules"), "gfortran")
        kinds_index = next(
            index for index, value in enumerate(command)
            if value.endswith("chempotpy_kinds.f90")
        )
        geometry_index = next(
            index for index, value in enumerate(command)
            if value.endswith("chempotpy_geometry.f90")
        )
        source_index = command.index(str(source_path))
        self.assertLess(kinds_index, geometry_index)
        self.assertLess(geometry_index, source_index)


if __name__ == "__main__":
    unittest.main()
