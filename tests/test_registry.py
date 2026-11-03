"""Tests for automatic CHEMPOTPY surface discovery."""

import unittest
from pathlib import Path

from chempotpy.core.registry import (
    discover_surface_sources,
    duplicate_surface_sources,
)


class RegistryTests( unittest.TestCase ):
    def test_discovers_library_sources( self ) -> None:
        sources = discover_surface_sources()
        self.assertGreater(len(sources), 300)
        self.assertGreater(len({source.system_name for source in sources}), 80)
        self.assertTrue(all(source.path.is_file() for source in sources))
        self.assertEqual(
            sum(source.entry_point == "pes" for source in sources), 323
        )
        self.assertEqual(
            sum(source.entry_point == "dpem" for source in sources), 40
        )
        self.assertEqual(sum(not source.is_surface for source in sources), 1)

    def test_reports_competing_fixed_and_free_form_sources( self ) -> None:
        duplicates = duplicate_surface_sources(discover_surface_sources())
        self.assertEqual(duplicates, {})


if __name__ == "__main__":
    unittest.main()
