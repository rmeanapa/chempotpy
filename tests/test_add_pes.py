"""Tests for the new-PES scaffolding command."""

import tempfile
import unittest
from pathlib import Path

from scripts.add_pes import scaffold_surface
from scripts.format_fortran import format_source_text


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TEMPLATE_ROOT = REPOSITORY_ROOT / "templates"


class AddPESTests(unittest.TestCase):
    def test_scaffolds_a_pes_in_an_existing_system(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            package_root = Path(temporary_directory) / "chempotpy"
            system_directory = package_root / "H3"
            system_directory.mkdir(parents=True)
            (system_directory / "__init__.py").write_text('"""H3."""\n')
            created = scaffold_surface(
                system="H3",
                surface="H3_TEST_2026",
                natoms=3,
                nstates=1,
                entry_point="pes",
                functional_form="analytic test",
                common_name="test surface",
                reference="Independent test reference",
                package_root=package_root,
                template_root=TEMPLATE_ROOT,
            )
            self.assertEqual(len(created), 1)
            source = created[0].read_text()
            self.assertEqual(source, format_source_text(source))
            self.assertIn("module chempotpy_h3_h3_test_2026_pes", source)
            self.assertIn("subroutine pes( x, igrad, p, g, d )", source)
            self.assertIn("CHEMPOTPY_IMPLEMENTATION_REQUIRED", source)

    def test_new_system_requires_and_records_canonical_atoms(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            package_root = Path(temporary_directory) / "chempotpy"
            created = scaffold_surface(
                system="H2X",
                surface="H2X_MODEL_2026",
                natoms=3,
                nstates=2,
                entry_point="dpem",
                functional_form="diabatic polynomial",
                common_name="model",
                reference="Independent test reference",
                atoms="H,H,X",
                package_root=package_root,
                template_root=TEMPLATE_ROOT,
            )
            self.assertEqual(len(created), 2)
            initialization = (package_root / "H2X" / "__init__.py").read_text()
            source = (package_root / "H2X" / "H2X_MODEL_2026.f90").read_text()
            self.assertIn("ATOM_SYMBOLS = ('H', 'H', 'X')", initialization)
            self.assertIn("subroutine dpem( x, igrad, u, ug )", source)

    def test_refuses_to_overwrite_an_existing_surface(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            package_root = Path(temporary_directory) / "chempotpy"
            system_directory = package_root / "H3"
            system_directory.mkdir(parents=True)
            (system_directory / "__init__.py").write_text('"""H3."""\n')
            (system_directory / "H3_TEST.f90").write_text("preserve me\n")
            with self.assertRaises(FileExistsError):
                scaffold_surface(
                    system="H3",
                    surface="H3_TEST",
                    natoms=3,
                    nstates=1,
                    entry_point="pes",
                    functional_form="test",
                    common_name="test",
                    reference="test reference",
                    package_root=package_root,
                    template_root=TEMPLATE_ROOT,
                )
            self.assertEqual(
                (system_directory / "H3_TEST.f90").read_text(),
                "preserve me\n",
            )

    def test_new_system_rejects_an_incomplete_atom_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            with self.assertRaisesRegex(ValueError, "exactly 3"):
                scaffold_surface(
                    system="H2X",
                    surface="H2X_TEST",
                    natoms=3,
                    nstates=1,
                    entry_point="pes",
                    functional_form="test",
                    common_name="test",
                    reference="test reference",
                    atoms="H,X",
                    package_root=Path(temporary_directory) / "chempotpy",
                    template_root=TEMPLATE_ROOT,
                )


if __name__ == "__main__":
    unittest.main()
