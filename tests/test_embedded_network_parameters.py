"""Regression tests for source-embedded neural-network parameters."""

import hashlib
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PARAMETER_ROOT = REPOSITORY_ROOT / "tests/data/neural_network_parameters"


class EmbeddedNetworkParameterTests( unittest.TestCase ):
    def test_generated_modules_record_provenance_hashes( self ) -> None:
        cases = (
            (
                PARAMETER_ROOT / "N4/weights.txt-Noconnected",
                PARAMETER_ROOT / "N4/biases.txt-Noconnected",
                REPOSITORY_ROOT
                / "chempotpy/N4/chempotpy_n4_pipnn_parameters.f90",
            ),
            (
                PARAMETER_ROOT / "FH2O/weights.txt",
                PARAMETER_ROOT / "FH2O/biases.txt",
                REPOSITORY_ROOT
                / "chempotpy/FH2O/chempotpy_fh2o_pipnn_parameters.f90",
            ),
        )
        for weights, biases, module in cases:
            with self.subTest(module=module.name):
                text = module.read_text()
                weights_hash = hashlib.sha256(weights.read_bytes()).hexdigest()
                biases_hash = hashlib.sha256(biases.read_bytes()).hexdigest()
                self.assertIn("Weights SHA-256: "+weights_hash, text)
                self.assertIn("Biases SHA-256: "+biases_hash, text)
                self.assertIn("implicit none", text)
                self.assertIn("private", text)

    def test_runtime_sources_do_not_read_model_files( self ) -> None:
        surfaces = (
            REPOSITORY_ROOT / "chempotpy/N4/N4_1A_PIPNN.f90",
            REPOSITORY_ROOT / "chempotpy/FH2O/FH2O_PIPNN.f90",
            REPOSITORY_ROOT / "chempotpy/FH2O/FH2O_SOC_PIPNN.f90",
        )
        forbidden = (
            "weights.txt",
            "biases.txt",
            "pes_init",
            "chempotpy_neural_network_io",
            "allocate(",
            "deallocate(",
        )
        for surface in surfaces:
            with self.subTest(surface=surface.name):
                text = surface.read_text().lower()
                for token in forbidden:
                    self.assertNotIn(token, text)
                self.assertIn("subroutine evaluate( x, igrad, p, g, d )", text)
                self.assertIn("subroutine pes( x, igrad, path, p, g, d )", text)

    def test_model_files_are_not_installed_as_runtime_data( self ) -> None:
        self.assertFalse((REPOSITORY_ROOT / "chempotpy/FH2O/weights.txt").exists())
        self.assertFalse((REPOSITORY_ROOT / "chempotpy/FH2O/biases.txt").exists())
        self.assertFalse(
            (REPOSITORY_ROOT / "chempotpy/N4/weights.txt-Noconnected").exists()
        )
        self.assertFalse(
            (REPOSITORY_ROOT / "chempotpy/N4/biases.txt-Noconnected").exists()
        )


if __name__ == "__main__":
    unittest.main()
