"""Tests for dependency-complete standalone PES exports."""

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
        self.assertIn("module chempotpy_h3_gen_a2_2002_pes", combined)
        self.assertLess(
            combined.index("module chempotpy_legacy_coordinates"),
            combined.index("module chempotpy_h3_gen_a2_2002_pes"),
        )

    def test_example_reports_numerical_gradient_and_hessian(self) -> None:
        source = select_surface("H3_GEN_A2_2002")
        driver = example_driver(source, geometry_count=1, difference_step=1.0e-3)
        self.assertIn("call pes( coordinates, 0, energy, gradient, coupling )", driver)
        self.assertIn("numerical_gradient", driver)
        self.assertIn("numerical_hessian", driver)

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


if __name__ == "__main__":
    unittest.main()
