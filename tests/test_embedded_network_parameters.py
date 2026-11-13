"""Regression tests for source-embedded neural-network parameters."""

import hashlib
import re
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PARAMETER_ROOT = REPOSITORY_ROOT / "tests/data/neural_network_parameters"


class EmbeddedNetworkParameterTests( unittest.TestCase ):
    def test_generated_modules_record_provenance_hashes( self ) -> None:
        single_file_models = {
            "N4": "N4/chempotpy_n4_pipnn_parameters.f90",
            "FH2O": "FH2O/chempotpy_fh2o_pipnn_parameters.f90",
            "H3ClO": "H3ClO/chempotpy_h3clo_pipnn_parameters.f90",
            "CH4OF": "CH4OF/chempotpy_ch4of_pipnn_parameters.f90",
            "CH2O": "CH2O/chempotpy_ch2o_pipnn_parameters.f90",
            "CH2O2": "CH2O2/chempotpy_ch2o2_pipnn_parameters.f90",
            "CH4OCl": "CH4OCl/chempotpy_ch4ocl_pipnn_parameters.f90",
            "O2H3": "O2H3/chempotpy_o2h3_pipnn_parameters.f90",
            "C2H2N2O": "C2H2N2O/chempotpy_c2h2n2o_pipnn_parameters.f90",
            "H3S": "H3S/chempotpy_h3s_pipnn_parameters.f90",
        }
        for system, module_relative_path in single_file_models.items():
            root = PARAMETER_ROOT / system
            weights = next(root.glob("weight*"))
            biases = next(root.glob("bias*"))
            module = REPOSITORY_ROOT / "chempotpy" / module_relative_path
            with self.subTest(module=module.name):
                self.assert_parameter_hash(module, "Weights", 1, weights)
                self.assert_parameter_hash(module, "Biases", 1, biases)
                self.assert_private_module(module)

    def test_ensemble_modules_record_every_network_hash( self ) -> None:
        models = {
            "C2H2O": "C2H2O/chempotpy_c2h2o_pipnn_parameters.f90",
            "HClOH": "HClOH/chempotpy_hcloh_pipnn_parameters.f90",
            "HO3": "HO3/chempotpy_ho3_pipnn_parameters.f90",
            "K2Rb2": "K2Rb2/chempotpy_k2rb2_pipnn_parameters.f90",
            "N2HOCp": "N2HOCp/chempotpy_n2hocp_pipnn_parameters.f90",
            "H3Cl": "H3Cl/chempotpy_h3cl_pipnn_parameters.f90",
            "ClH2O": "ClH2O/chempotpy_clh2o_pipnn_parameters.f90",
        }
        for system, module_relative_path in models.items():
            root = PARAMETER_ROOT / system
            weights = sorted(root.glob("weight*"))
            biases = sorted(root.glob("bias*"))
            module = REPOSITORY_ROOT / "chempotpy" / module_relative_path
            with self.subTest(module=module.name):
                self.assertEqual(len(weights), 3)
                self.assertEqual(len(biases), 3)
                for index, path in enumerate(weights, 1):
                    self.assert_parameter_hash(module, "Weights", index, path)
                for index, path in enumerate(biases, 1):
                    self.assert_parameter_hash(module, "Biases", index, path)
                self.assert_private_module(module)

    def test_no2_modules_record_provenance_hashes( self ) -> None:
        for state in ("2Ap", "4Ap", "6Ap"):
            stem = "NO2_{}_PIPNN".format(state)
            module = REPOSITORY_ROOT / (
                "chempotpy/NO2/chempotpy_no2_{}_pipnn_parameters.f90".format(
                    state.lower()
                )
            )
            self.assert_parameter_hash(
                module, "Weights", 1, PARAMETER_ROOT / "NO2/weights_{}.txt".format(stem)
            )
            self.assert_parameter_hash(
                module, "Biases", 1, PARAMETER_ROOT / "NO2/biases_{}.txt".format(stem)
            )
            self.assert_private_module(module)

    def test_c2n_generated_modules_record_all_network_hashes( self ) -> None:
        cases = (
            ("prime", "chempotpy_c2n_pipnn_ap_parameters.f90"),
            ("dprime", "chempotpy_c2n_pipnn_app_parameters.f90"),
        )
        for state, module_name in cases:
            module = REPOSITORY_ROOT / "chempotpy/C2N" / module_name
            text = module.read_text()
            for index, label in enumerate("abc", 1):
                source = PARAMETER_ROOT / (
                    "C2N/c2n_2a_{}_para_{}.txt".format(state, label)
                )
                digest = hashlib.sha256(source.read_bytes()).hexdigest()
                self.assertIn(
                    "Network {} SHA-256: {}".format(index, digest), text
                )
            self.assert_private_module(module)

    def test_clh2o_module_records_spline_provenance( self ) -> None:
        source = (
            REPOSITORY_ROOT
            / "tests/data/external_parameters/ClH2O/spline_tables.txt"
        )
        module = (
            REPOSITORY_ROOT
            / "chempotpy/ClH2O/chempotpy_clh2o_pipnn_parameters.f90"
        )
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        self.assertIn("Spline tables SHA-256: {}".format(digest), module.read_text())

    def test_runtime_sources_do_not_read_model_files( self ) -> None:
        forbidden = (
            "chempotpy_neural_network_io",
            "pes_init",
            "allocate(",
            "deallocate(",
        )
        provenance_names = {
            path.name.lower() for path in PARAMETER_ROOT.rglob("*") if path.is_file()
        }
        for surface in REPOSITORY_ROOT.glob("chempotpy/*/*PIPNN*.f90"):
            if "parameters" in surface.name.lower():
                continue
            with self.subTest(surface=surface.name):
                text = surface.read_text().lower()
                for token in forbidden:
                    self.assertNotIn(token, text)
                for name in provenance_names:
                    self.assertNotIn(name, text)

    def test_model_files_are_not_installed_as_runtime_data( self ) -> None:
        for source in PARAMETER_ROOT.rglob("*"):
            if source.is_file():
                relative_path = source.relative_to(PARAMETER_ROOT)
                self.assertFalse(
                    (REPOSITORY_ROOT / "chempotpy" / relative_path).exists()
                )

    def assert_parameter_hash(
        self, module: Path, label: str, index: int, source: Path
    ) -> None:
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        self.assertIn(
            "{} {} SHA-256: {}".format(label, index, digest),
            module.read_text(),
        )

    def assert_private_module( self, module: Path ) -> None:
        text = module.read_text()
        self.assertRegex(text, re.compile(r"^\s*implicit none\s*$", re.MULTILINE))
        self.assertRegex(text, re.compile(r"^\s*private\s*$", re.MULTILINE))


if __name__ == "__main__":
    unittest.main()
