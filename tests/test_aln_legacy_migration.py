"""Regression tests for Aln kernels migrated from released fixed form."""

import json
import shutil
import tempfile
import unittest
from pathlib import Path

from scripts.capture_pes_references import capture_source


class AlnLegacyMigrationTests( unittest.TestCase ):
    @classmethod
    def setUpClass( cls ) -> None:
        cls.repository = Path(__file__).resolve().parents[1]
        fixture = cls.repository / "tests" / "data" / \
            "aln_legacy_references.json"
        cls.references = json.loads(fixture.read_text())

    def test_modern_sources_match_legacy_references( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        with tempfile.TemporaryDirectory(prefix="chempotpy-aln-migration-") as tmp:
            directory = Path(tmp)
            for surface_name, reference in self.references["surfaces"].items():
                with self.subTest(surface=surface_name):
                    source = self.repository / "chempotpy" / "Aln" / \
                        (surface_name + ".f90")
                    self.assertTrue(source.is_file())
                    actual = capture_source(compiler, source, directory)
                    self.assertEqual(
                        actual["energies_ev"],
                        reference["energies_ev"],
                        "modern source changed a released energy",
                    )

    def test_reference_fixture_has_legacy_provenance( self ) -> None:
        self.assertEqual(len(self.references["surfaces"]), 42)
        for surface_name, reference in self.references["surfaces"].items():
            with self.subTest(surface=surface_name):
                digest = reference["legacy_source_sha256"]
                self.assertEqual(len(digest), 64)
                int(digest, 16)

    def test_all_aln_sources_are_free_form( self ) -> None:
        fixed_sources = sorted(
            (self.repository / "chempotpy" / "Aln").glob("*.f")
        )
        self.assertEqual(fixed_sources, [])


if __name__ == "__main__":
    unittest.main()
