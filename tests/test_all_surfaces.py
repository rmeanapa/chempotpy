"""Catalog-wide numerical regression tests for every callable PES."""

import concurrent.futures
import hashlib
import json
import math
import re
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

    def test_legacy_storage_and_setup_constructs_are_absent(self) -> None:
        definitions = 0
        calls = 0
        common_blocks = 0
        block_data_units = 0
        entry_points = 0
        for source in (self.repository / "chempotpy").rglob("*.f90"):
            text = source.read_text(errors="replace")
            definitions += len(re.findall(
                r"(?im)^\s*subroutine\s+prepot[a-z0-9_]*\b",
                text,
            ))
            calls += len(re.findall(
                r"(?im)^\s*call\s+prepot[a-z0-9_]*\b",
                text,
            ))
            common_blocks += len(re.findall(r"(?im)^\s*common\b", text))
            block_data_units += len(re.findall(
                r"(?im)^\s*block\s+data\b",
                text,
            ))
            entry_points += len(re.findall(r"(?im)^\s*entry\b", text))
        self.assertEqual(definitions, 0)
        self.assertEqual(calls, 0)
        self.assertEqual(common_blocks, 0)
        self.assertEqual(block_data_units, 0)
        self.assertEqual(entry_points, 0)

    def test_migrated_fits_have_no_initialization_lifecycle(self) -> None:
        targets = (
            self.repository / "chempotpy" / "core" / "fortran" / "chempotpy_h3_basis.f90",
            self.sources["ClH2/ClH2_LEPS_ALTG2_1996"],
            self.sources["ClH2/ClH2_LEPS_STSBLTG1_1989"],
            self.sources["ClH2/ClH2_LEPS_STSBLTG2_1989"],
            self.sources["ClH2/ClH2_LEPS_GTM4_1982"],
            self.sources["ClH2/ClH2_LEPS_TTGI1_1985"],
            self.sources["ClH2/ClH2_LEPS_TTGI2_1985"],
            *(self.sources["ClH2/{}".format(name)] for name in (
                "ClH2_LEPS_GTM1_1981",
                "ClH2_LEPS_GTM2_1981",
                "ClH2_LEPS_GTM3_1981",
                "ClH2_LEPS_KNPRY_1966",
                "ClH2_LEPS_PK1_1966",
                "ClH2_LEPS_PK2_1966",
                "ClH2_LEPS_SPK_1973",
            )),
            self.sources["OH2/OH2_LEPS_JW1977"],
            self.sources["OH2/OH2_LEPS_JWS1985"],
            self.sources["Cl2H/Cl2H_LEPS_KNPRY_1966"],
            self.sources["Cl2H/Cl2H_LEPS_BCMR_1983"],
            self.sources["Cl2H/Cl2H_LEPS_PK2_1987"],
            self.sources["Cl2H/Cl2H_LEPS_PK3_1987"],
            *(self.sources["FH2/{}".format(name)] for name in (
                "FH2_LEPS_M5SK_1966",
                "FH2_LEPS_M5_1991",
                "FH2_LEPS_NO1_1984",
                "FH2_LEPS_SK_1980",
            )),
            self.sources["F2H/F2H_LEPS_KNPRY_1966"],
            self.sources["BrH2/BrH2_LEPS_W_1959"],
            *(self.sources["BrHCl/BrHCl_LEPS{}".format(index)] for index in range(1, 4)),
            self.sources["IH2/IH2_LEPS_PPW1_1974"],
            self.sources["IH2/IH2_LEPS_KNPRY_1966"],
            self.sources["I2H/I2H_LEPS_PT1_1971"],
            self.sources["I2H/I2H_LEPS_PT2_1971"],
            self.sources["NH3/NH3_GEN_FFW1_2006"],
            self.sources["NH3/NH3_GEN_FFW1_2006_DPEM"],
            self.sources["NH3/NH3_GEN_FFW2_2007"],
            self.sources["NH3/NH3_GEN_FFW2_2007_DPEM"],
            self.sources["H2OBr/H2OBr_PIP"],
            *(self.sources["H3/{}".format(name)] for name in (
                "H3_GEN_A2_2002",
                "H3_GEN_A3_2002",
                "H3_GEN_A4_2002",
                "H3_GEN_CCI_2002",
            )),
            *(source for key, source in self.sources.items()
              if key.startswith("MCH/")),
            *(source for key, source in self.sources.items()
              if key.startswith("YRH/")),
            *(source for key, source in self.sources.items()
              if key.startswith("MXH/")),
            *(source for key, source in self.sources.items()
              if key.startswith("NaH2/")),
            *(source for key, source in self.sources.items()
              if key.startswith("LiFH/")),
            *(source for key, source in self.sources.items()
              if key.startswith("NaFH/")),
            self.sources["H3/H3_GEN_BKMP_1991"],
            self.sources["H3/H3_GEN_BH_2009"],
            self.sources["H3/H3_GEN_PK2_1964"],
            self.sources["BrH2/BrH2_LEPS_LTBZ2_1995"],
            self.sources["BrH2/BrH2_LEPS_LTBZ2_1995_DPEM"],
            self.sources["CH3OH/CH3OH_PIP"],
            self.sources["CH4_H2O_H2O/CH4_H2O_H2O_PIP"],
            self.sources["HOCO_Ar/HOCO_Ar_PIP"],
            self.sources["O2H/O2H_GEN"],
        )
        for source in targets:
            with self.subTest(source=source.name):
                self.assertNotRegex(
                    source.read_text(errors="replace"),
                    r"(?im)\b(?:prepot[a-z0-9_]*|initialize_h3_basis|"
                    r"initialize_only|first_time_data|first_time_diag)\b",
                )

    def test_every_surface_has_one_private_owner_module(self) -> None:
        for key, source in self.sources.items():
            with self.subTest(surface=key):
                owner = "chempotpy_{}_pes".format(
                    re.sub(r"[^a-z0-9_]", "_", source.stem.lower())
                )
                text = source.read_text()
                modules = re.findall(
                    r"^\s*module\s+{}\b".format(re.escape(owner)),
                    text,
                    flags=re.IGNORECASE | re.MULTILINE,
                )
                self.assertEqual(len(modules), 1)
                body = re.search(
                    r"^\s*module\s+{}\b(.*?)"
                    r"^\s*end\s+module\s+{}\b".format(
                        re.escape(owner), re.escape(owner)
                    ),
                    text,
                    flags=re.IGNORECASE | re.MULTILINE | re.DOTALL,
                )
                self.assertIsNotNone(body)
                self.assertRegex(body.group(1), r"(?im)^\s*private\s*$")
                self.assertRegex(
                    body.group(1),
                    r"(?im)^\s*public\s*::\s*evaluate\s*$",
                )

    def test_owner_only_surfaces_have_no_external_helpers(self) -> None:
        for key, source in self.sources.items():
            with self.subTest(surface=key):
                owner = "chempotpy_{}_pes".format(
                    re.sub(r"[^a-z0-9_]", "_", source.stem.lower())
                )
                text = source.read_text()
                module = re.search(
                    r"^\s*module\s+{}\b.*?"
                    r"^\s*end\s+module\s+{}\s*$".format(
                        re.escape(owner), re.escape(owner)
                    ),
                    text,
                    flags=re.IGNORECASE | re.MULTILINE | re.DOTALL,
                )
                self.assertIsNotNone(module)
                outside = text[:module.start()]+text[module.end():]
                adapter = re.search(
                    r"^!>\s*@brief\s+Preserve the historical.*?"
                    r"^\s*end\s+subroutine(?:\s+(?:pes|dpem))?\s*$",
                    outside,
                    flags=re.IGNORECASE | re.MULTILINE | re.DOTALL,
                )
                self.assertIsNotNone(adapter)
                outside = outside[:adapter.start()]+outside[adapter.end():]
                if re.search(
                    r"^\s*(?:module\s+(?!procedure\b)|block\s+data\b)",
                    outside,
                    flags=re.IGNORECASE | re.MULTILINE,
                ):
                    continue
                self.assertNotRegex(
                    outside,
                    r"(?im)^\s*(?:subroutine|"
                    r"(?:(?:double\s+precision|real(?:\([^)]*\))?|"
                    r"integer)\s+)?function)\b",
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
