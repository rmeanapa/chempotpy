"""Catalog-wide numerical regression tests for every callable PES."""

import concurrent.futures
import hashlib
import json
import math
import shutil
import tempfile
import unittest
from pathlib import Path

from chempotpy.core.registry import discover_surface_sources
from scripts.capture_library_references import capture_source


class AllSurfaceRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repository = Path(__file__).resolve().parents[1]
        fixture = cls.repository / "tests" / "data" / \
            "all_surface_references.json"
        cls.references = json.loads(fixture.read_text())["surfaces"]
        cls.sources = {
            "{}/{}".format(source.system_name, source.surface_name): source.path
            for source in discover_surface_sources()
            if source.is_surface
        }

    def test_fixture_covers_the_callable_catalog(self) -> None:
        self.assertEqual(len(self.references), 363)
        self.assertEqual(set(self.references), set(self.sources))

    def test_fixture_tracks_the_authoritative_sources(self) -> None:
        for key, source in self.sources.items():
            with self.subTest(surface=key):
                digest = hashlib.sha256(source.read_bytes()).hexdigest()
                self.assertEqual(
                    self.references[key]["source_sha256"], digest
                )

    def test_every_surface_matches_its_reference_outputs(self) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        sources = [self.sources[key] for key in sorted(self.references)]
        with tempfile.TemporaryDirectory(
            prefix="chempotpy-all-surfaces-"
        ) as tmp:
            work_root = Path(tmp)
            with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
                results = list(pool.map(
                    lambda source: capture_source(
                        compiler, source, work_root
                    ),
                    sources,
                ))
        for result in results:
            with self.subTest(surface=result.key):
                self.assertFalse(result.error, result.error)
                expected_record = self.references[result.key]
                for field in ("entry_point", "natoms", "nstates", "signature"):
                    self.assertEqual(
                        result.record[field], expected_record[field]
                    )
                expected = expected_record["values"]
                actual = result.record["values"]
                self.assertEqual(len(actual), len(expected))
                for expected_value, actual_value in zip(expected, actual):
                    if expected_value == actual_value:
                        continue
                    expected_number = float(expected_value)
                    actual_number = float(actual_value)
                    if (
                        abs(expected_number) < 1.0e-300
                        and abs(actual_number) < 1.0e-300
                    ):
                        continue
                    self.assertTrue(math.isclose(
                        actual_number,
                        expected_number,
                        rel_tol=1.0e-13,
                        abs_tol=1.0e-13,
                    ))


if __name__ == "__main__":
    unittest.main()
