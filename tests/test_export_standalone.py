"""Tests for dependency-complete standalone PES exports."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from chempotpy.core.registry import discover_surface_sources
from scripts.export_standalone import (
    combined_source,
    example_driver,
    export_surface,
    select_surface,
)


class ExportStandaloneTests(unittest.TestCase):
    def test_every_callable_surface_can_generate_an_example(self) -> None:
        sources = tuple(
            source for source in discover_surface_sources() if source.is_surface
        )
        self.assertEqual(len(sources), 363)
        for source in sources:
            with self.subTest(
                surface="{}/{}".format(
                    source.system_name, source.surface_name
                )
            ):
                driver = example_driver(
                    source, geometry_count=1, difference_step=1.0e-3
                )
                self.assertIn("program chempotpy_standalone_example", driver)

    def test_combined_source_contains_core_before_surface(self) -> None:
        source = select_surface("H3/H3_GEN_A2_2002")
        combined = combined_source(source)
        self.assertIn("module chempotpy_legacy_coordinates", combined)
        self.assertIn("module chempotpy_h3_2002_basis", combined)
        self.assertIn("module chempotpy_h3_gen_a2_2002_pes", combined)
        self.assertLess(
            combined.index("module chempotpy_legacy_coordinates"),
            combined.index("module chempotpy_h3_2002_basis"),
        )
        self.assertLess(
            combined.index("module chempotpy_h3_2002_basis"),
            combined.index("module chempotpy_h3_gen_a2_2002_pes"),
        )

    def test_example_reports_numerical_gradient_and_hessian(self) -> None:
        source = select_surface("H3_GEN_A2_2002")
        driver = example_driver(source, geometry_count=1, difference_step=1.0e-3)
        self.assertIn("call pes( coordinates, 0, energy, gradient, coupling )", driver)
        self.assertIn("call central_finite_difference", driver)
        self.assertIn("numerical_gradient", driver)
        self.assertIn("numerical_hessian", driver)
        self.assertNotIn("mixed_pp", driver)

    def test_export_writes_source_and_optional_driver(self) -> None:
        source = select_surface("H3/H3_GEN_A2_2002")
        with tempfile.TemporaryDirectory() as temporary_directory:
            outputs = export_surface(
                source,
                Path(temporary_directory),
                with_driver=True,
                geometry_count=1,
                difference_step=1.0e-3,
            )
            self.assertEqual(len(outputs), 2)
            self.assertTrue(all(path.is_file() for path in outputs))
            self.assertIn(
                "module chempotpy_finite_difference",
                outputs[0].read_text(),
            )

    def test_exported_example_compiles_and_runs( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        source = select_surface("H3/H3_GEN_A2_2002")
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            source_output, driver_output = export_surface(
                source,
                directory,
                with_driver=True,
                geometry_count=1,
                difference_step=1.0e-3,
            )
            executable = directory / "example"
            subprocess.run(
                [
                    compiler,
                    "-ffree-line-length-none",
                    str(source_output),
                    str(driver_output),
                    "-o",
                    str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            completed = subprocess.run(
                [str(executable)],
                check=True,
                capture_output=True,
                text=True,
            )
        self.assertIn("numerical_gradient", completed.stdout)
        self.assertIn("numerical_hessian", completed.stdout)


if __name__ == "__main__":
    unittest.main()
