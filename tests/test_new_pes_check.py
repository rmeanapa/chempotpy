"""Tests for enforceable new-PES onboarding requirements."""

import hashlib
import tempfile
import unittest
from pathlib import Path

from chempotpy.core.registry import SurfaceSource
from scripts.add_pes import scaffold_surface
from scripts.capture_library_references import (
    BASELINE_VERSION,
    REGRESSION_GEOMETRY_COUNT,
)
from scripts.check_new_pes import (
    fixture_violations,
    metadata_values,
    source_contract_violations,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_ROOT = REPOSITORY_ROOT / "templates"


class NewPESCheckTests(unittest.TestCase):
    def scaffold(self, directory: Path) -> Path:
        """Return one generated PES source in a temporary package."""

        package_root = directory / "chempotpy"
        system_directory = package_root / "H3"
        system_directory.mkdir(parents=True)
        (system_directory / "__init__.py").write_text('"""H3."""\n')
        return scaffold_surface(
            system="H3",
            surface="H3_TEST_2026",
            natoms=3,
            nstates=1,
            entry_point="pes",
            functional_form="analytic test",
            common_name="test surface",
            reference="Independent reference values",
            package_root=package_root,
            template_root=TEMPLATE_ROOT,
        )[0]

    def test_template_contains_complete_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = self.scaffold(Path(temporary_directory))
            metadata, violations = metadata_values(source.read_text())
            self.assertEqual(violations, [])
            self.assertEqual(metadata["System"], "H3")
            self.assertEqual(metadata["Number of electronic surfaces"], "1")

    def test_unfinished_scaffold_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = self.scaffold(Path(temporary_directory))
            violations = source_contract_violations(source)
            self.assertIn(
                "remove every scaffold placeholder before review", violations
            )

    def test_external_parameter_reads_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source = self.scaffold(Path(temporary_directory))
            source.write_text(
                source.read_text().replace(
                    "p = 0.e0_wp;",
                    "open(newunit=igrad, file='weights.txt')\n"
                    "        read(igrad, *) p\n"
                    "        p = 0.e0_wp;",
                )
            )
            violations = source_contract_violations(source)
            self.assertTrue(any("external file OPEN" in item for item in violations))
            self.assertTrue(any("external file READ" in item for item in violations))

    def test_complete_fixture_record_passes_provenance_checks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            source_path = self.scaffold(Path(temporary_directory))
            source = SurfaceSource(
                system_name="H3",
                surface_name="H3_TEST_2026",
                path=source_path,
                source_form="free",
                entry_point="pes",
            )
            record = {
                "entry_point": "pes",
                "natoms": 3,
                "nstates": 1,
                "signature": "x,igrad,p,g,d",
                "source_sha256": hashlib.sha256(
                    source_path.read_bytes()
                ).hexdigest(),
                "values": ["0.0"]*REGRESSION_GEOMETRY_COUNT,
            }
            fixture = {
                "baseline_version": BASELINE_VERSION,
                "surfaces": {"H3/H3_TEST_2026": record},
            }
            self.assertEqual(fixture_violations(source, fixture), [])


if __name__ == "__main__":
    unittest.main()
