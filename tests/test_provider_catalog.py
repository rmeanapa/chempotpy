"""Tests for the installed catalog used by all external bridges."""

import unittest

from chempotpy.core.catalog import resolve_atom_symbols, surface_catalog


class ProviderCatalogTests(unittest.TestCase):
    def test_catalog_covers_every_callable_surface(self) -> None:
        records = surface_catalog()
        self.assertEqual(len(records), 414)
        self.assertEqual(len({record.key for record in records}), 414)

    def test_every_fixed_system_resolves_its_atom_symbols(self) -> None:
        for record in surface_catalog():
            if not record.variable_atom_count:
                with self.subTest(surface=record.key):
                    self.assertEqual(
                        len(resolve_atom_symbols(record)), record.atom_count
                    )

    def test_variable_system_accepts_explicit_symbols(self) -> None:
        record = next(
            item for item in surface_catalog() if item.system_name == "Aln"
        )
        self.assertEqual(resolve_atom_symbols(record, ["Al"] * 4), ("Al",) * 4)


if __name__ == "__main__":
    unittest.main()
