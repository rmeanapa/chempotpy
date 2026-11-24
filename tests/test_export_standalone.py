"""Tests for dependency-complete standalone PES exports."""

import json
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
        fixture = json.loads(
            (
                Path(__file__).resolve().parent
                / "data"
                / "all_surface_references.json"
            ).read_text()
        )["surfaces"]
        self.assertEqual(len(sources), len(fixture))
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

    def test_co2_export_orders_rkhs_dependencies(self) -> None:
        source = select_surface("CO2/CO2_1Ap_RKHS")
        combined = combined_source(source)
        self.assertLess(
            combined.index("module chempotpy_rkhs"),
            combined.index("module chempotpy_co2_rkhs_runtime"),
        )
        self.assertLess(
            combined.index("module chempotpy_co2_rkhs_runtime"),
            combined.index("module chempotpy_co2_1ap_rkhs_pes"),
        )

    def test_formic_acid_dimer_export_orders_invariant_dependencies(self) -> None:
        source = select_surface("C2H4O4/C2H4O4_FAD_PIP")
        combined = combined_source(source)
        self.assertLess(
            combined.index("module chempotpy_c2h4o4_parameters"),
            combined.index("module chempotpy_c2h4o4_fad_pip_pes"),
        )
        self.assertLess(
            combined.index("module chempotpy_c2h4o4_invariants"),
            combined.index("module chempotpy_c2h4o4_fad_pip_pes"),
        )

    def test_malonaldehyde_export_orders_parameter_dependency(self) -> None:
        source = select_surface(
            "malonaldehyde/malonaldehyde_MHT_2014"
        )
        combined = combined_source(source)
        self.assertLess(
            combined.index(
                "module chempotpy_malonaldehyde_mht_parameters"
            ),
            combined.index(
                "module chempotpy_malonaldehyde_mht_2014_pes"
            ),
        )

    def test_hni_export_orders_fixed_environment_dependency(self) -> None:
        source = select_surface("HNi/HNi_EAM5_1995")
        combined = combined_source(source)
        self.assertLess(
            combined.index("module chempotpy_hni_eam5_parameters"),
            combined.index("module chempotpy_hni_eam5_1995_pes"),
        )
        self.assertLess(
            combined.index("module chempotpy_hni_ni111_parameters"),
            combined.index("module chempotpy_hni_eam5_1995_pes"),
        )

    def test_hni_eam67_exports_order_runtime_dependencies(self) -> None:
        for surface_name, module_name in (
            ("HNi/HNi_EAM6_1996", "chempotpy_hni_eam6_1996_pes"),
            ("HNi/HNi_EAM7_1998", "chempotpy_hni_eam7_1998_pes"),
        ):
            combined = combined_source(select_surface(surface_name))
            self.assertLess(
                combined.index("module chempotpy_hni_eam67_parameters"),
                combined.index("module chempotpy_hni_eam67_runtime"),
            )
            self.assertLess(
                combined.index("module chempotpy_hni_eam67_runtime"),
                combined.index("module {}".format(module_name)),
            )
            self.assertLess(
                combined.index("module chempotpy_hni_ni111_parameters"),
                combined.index("module {}".format(module_name)),
            )

    def test_ch4ni_export_orders_fixed_environment_and_runtime(self) -> None:
        combined = combined_source(select_surface(
            "CH4Ni/CH4Ni_EAM6_J1_1996"
        ))
        self.assertLess(
            combined.index("module chempotpy_ch4ni_ni111_parameters"),
            combined.index("module chempotpy_ch4ni_runtime"),
        )
        self.assertLess(
            combined.index("module chempotpy_ch4ni_runtime"),
            combined.index("module chempotpy_ch4ni_eam6_j1_1996_pes"),
        )
        self.assertLess(
            combined.index("module chempotpy_finite_difference"),
            combined.index("module chempotpy_ch4ni_eam6_j1_1996_pes"),
        )

    def test_example_reports_numerical_gradient_and_hessian(self) -> None:
        source = select_surface("H3_GEN_A2_2002")
        driver = example_driver(source, geometry_count=1, difference_step=1.0e-3)
        self.assertIn("interface", driver)
        self.assertIn("subroutine pes( x, igrad, energy, gradient, coupling )", driver)
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
        source = select_surface("HeH2p/HeH2p_DIM_1976")
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
                    "-std=f2018",
                    "-pedantic",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    "-Wimplicit-interface",
                    "-fimplicit-none",
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
