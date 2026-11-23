"""Tests for tolerance-guarded reference hash refreshes."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts.refresh_library_reference_hashes import equivalent_values


class RefreshLibraryReferenceHashesTests(unittest.TestCase):
    """Verify that hash refreshes accept roundoff, not numerical changes."""

    def test_identical_values_are_equivalent(self) -> None:
        self.assertTrue(equivalent_values(["1.0"], ["1.0"]))

    def test_roundoff_level_values_are_equivalent(self) -> None:
        self.assertTrue(equivalent_values(["1.0"], ["1.000000000005"]))

    def test_materially_changed_values_are_rejected(self) -> None:
        self.assertFalse(equivalent_values(["1.0"], ["1.000001"]))

    def test_different_value_counts_are_rejected(self) -> None:
        self.assertFalse(equivalent_values(["1.0"], ["1.0", "2.0"]))

    def test_add_new_requires_an_explicit_review_flag(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            fixture = directory / "fixture.json"
            capture = directory / "capture.json"
            fixture.write_text(json.dumps({
                "baseline_version": 4,
                "surfaces": {},
            }))
            record = {
                "entry_point": "pes",
                "natoms": 3,
                "nstates": 1,
                "signature": "x,igrad,p,g,d",
                "source_sha256": "abc",
                "values": ["1.0"],
            }
            capture.write_text(json.dumps({
                "baseline_version": 4,
                "surfaces": {"H3/new": record},
            }))
            script = (
                Path(__file__).resolve().parents[1]
                / "scripts"
                / "refresh_library_reference_hashes.py"
            )
            rejected = subprocess.run(
                [sys.executable, str(script), str(fixture), str(capture)],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(rejected.returncode, 0)
            subprocess.run(
                [
                    sys.executable,
                    str(script),
                    str(fixture),
                    str(capture),
                    "--add-new",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            merged = json.loads(fixture.read_text())
            self.assertEqual(merged["surfaces"]["H3/new"], record)


if __name__ == "__main__":
    unittest.main()
