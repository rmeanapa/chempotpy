"""Regression tests for the repository-wide fixed-form PES migration."""

import concurrent.futures
import json
import math
import shutil
import tempfile
import unittest
from pathlib import Path

from scripts.capture_library_references import capture_source


CORRECTED_RELEASE_DEFECTS = {
    # The released wrapper reapplied its unit conversion before every point.
    "ClH2/ClH2_LEPS_ALTG2_1996": [
        "3.85497619023016469E+000",
        "2.50328966056197944E+000",
    ],
    # The released wrappers skipped the basis and damping setup completely.
    "H3/H3_GEN_A2_2002": [
        "5.25459088704119393E-001",
        "4.87261550575933011E-001",
    ],
    "H3/H3_GEN_A3_2002": [
        "5.58475196600794588E-001",
        "5.35016594053392747E-001",
    ],
    "H3/H3_GEN_A4_2002": [
        "5.56554524813657303E-001",
        "5.35268464705464853E-001",
    ],
    "H3/H3_GEN_CCI_2002": [
        "5.55896650189018793E-001",
        "5.35434886933361942E-001",
    ],
    "BrH2/BrH2_LEPS_LTBZSO_1995": [
        "8.25785718003896285E+000",
        "2.16553758038227357E+001",
        "6.08016929813028817E+000",
        "1.85708590815010552E+001",
    ],
    "BrH2/BrH2_LEPS_LTBZSO_1995_DPEM": [
        "8.26131912680582481E+000",
        "-2.15335810324679189E-001",
        "-2.15335810324679189E-001",
        "2.16519138570558702E+001",
        "6.08388272801880348E+000",
        "-2.15335810324679189E-001",
        "-2.15335810324679189E-001",
        "1.85671456516125453E+001",
    ],
}


class LibraryLegacyMigrationTests( unittest.TestCase ):
    @classmethod
    def setUpClass( cls ) -> None:
        cls.repository = Path(__file__).resolve().parents[1]
        fixture = cls.repository / "tests" / "data" / \
            "library_legacy_references.json"
        cls.references = json.loads(fixture.read_text())["surfaces"]

    def test_modern_sources_match_legacy_references( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        sources = [
            self.repository / "chempotpy" / (key + ".f90")
            for key in self.references
        ]
        with tempfile.TemporaryDirectory(
            prefix="chempotpy-library-migration-"
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
                expected = CORRECTED_RELEASE_DEFECTS.get(
                    result.key,
                    self.references[result.key]["values"],
                )
                actual = result.record["values"]
                self.assertEqual(len(actual), len(expected))
                for expected_value, actual_value in zip(expected, actual):
                    if expected_value == actual_value:
                        continue
                    expected_number = float(expected_value)
                    actual_number = float(actual_value)
                    self.assertTrue(math.isclose(
                        actual_number,
                        expected_number,
                        rel_tol=1.0e-13,
                        abs_tol=1.0e-13,
                    ))

    def test_reference_fixture_has_legacy_provenance( self ) -> None:
        self.assertEqual(len(self.references), 205)
        for surface_name, reference in self.references.items():
            with self.subTest(surface=surface_name):
                digest = reference["legacy_source_sha256"]
                self.assertEqual(len(digest), 64)
                int(digest, 16)

    def test_repository_has_no_fixed_form_pes_sources( self ) -> None:
        fixed_sources = sorted(
            (self.repository / "chempotpy").glob("*/*.f")
        )
        self.assertEqual(fixed_sources, [])


if __name__ == "__main__":
    unittest.main()
