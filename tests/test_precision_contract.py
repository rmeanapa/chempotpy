"""Tests for the centralized Fortran working-precision contract."""

import re
import unittest
from pathlib import Path


WORKING_PRECISION_CORE_SOURCES = (
    "chempotpy_blas.f90",
    "chempotpy_dsyev.f90",
    "chempotpy_geometry.f90",
    "chempotpy_internal_coordinates.f90",
    "chempotpy_legacy_coordinates.f90",
    "chempotpy_legacy_units.f90",
    "chempotpy_physical_constants.f90",
    "chempotpy_units.f90",
)


class PrecisionContractTests( unittest.TestCase ):
    @classmethod
    def setUpClass( cls ) -> None:
        cls.repository = Path(__file__).resolve().parents[1]
        cls.core = cls.repository / "chempotpy" / "core"

    @staticmethod
    def executable_fortran( source: str ) -> str:
        """Return Fortran code with comments and character contents removed."""

        output = []
        quote = ""
        index = 0
        while index < len(source):
            character = source[index]
            if quote:
                if character == quote:
                    if index+1 < len(source) and source[index+1] == quote:
                        output.extend("  ")
                        index += 2
                        continue
                    quote = ""
                output.append("\n" if character == "\n" else " ")
            elif character in {"'", '"'}:
                quote = character
                output.append(" ")
            elif character == "!":
                while index < len(source) and source[index] != "\n":
                    output.append(" ")
                    index += 1
                continue
            else:
                output.append(character)
            index += 1
        return "".join(output)

    def test_working_precision_defaults_to_portable_double( self ) -> None:
        source = (self.core / "chempotpy_kinds.f90").read_text()
        self.assertRegex(
            source,
            r"integer, parameter, public :: dp\s*=\s*"
            r"selected_real_kind\(\s*15,\s*307\s*\)",
        )
        self.assertRegex(
            source,
            r"integer, parameter, public :: qp\s*=\s*"
            r"selected_real_kind\(\s*33,\s*4931\s*\)",
        )
        self.assertRegex(
            source,
            r"integer, parameter, public :: wp\s*=\s*dp",
        )

    def test_shared_numerical_core_uses_working_precision( self ) -> None:
        hard_wired_kind = re.compile(
            r"\bdouble\s+precision\b|\breal\s*\(\s*8\s*\)|"
            r"\breal\s*\(\s*real64\s*\)|_real64\b",
            flags=re.IGNORECASE,
        )
        for source_name in WORKING_PRECISION_CORE_SOURCES:
            with self.subTest(source=source_name):
                source = self.executable_fortran(
                    (self.core / source_name).read_text()
                )
                self.assertIn("wp", source)
                self.assertIsNone(hard_wired_kind.search(source))

    def test_catalog_uses_only_the_shared_working_precision( self ) -> None:
        hard_wired_precision = re.compile(
            r"\bdouble\s*precision|\breal\s*\*\s*8\b|"
            r"\breal\s*\(\s*(?:kind\s*=\s*)?8\s*\)|"
            r"\breal\s*\(\s*real64\s*\)|_real64\b|"
            r"\breal\s*\([^\n)]*,\s*(?:kind\s*=\s*)?8\s*\)|"
            r"(?<![A-Za-z0-9_])(?:\d+\.\d*|\.\d+|\d+)[dD][+-]?\d+|"
            r"(?<![A-Za-z0-9_])(?:\d+\.\d*|\.\d+|\d+)"
            r"(?:[eE][+-]?\d+)?_dp\b|"
            r"\b(?:dabs|dacos|dasin|datan2|datan|dcosh|dcos|dexp|"
            r"dlog10|dlog|dmax1|dmin1|dsin|dsqrt|dtanh|dtan|dble|"
            r"dcmplx)\s*\(",
            flags=re.IGNORECASE,
        )
        for source in sorted(self.repository.glob("chempotpy/**/*.f90")):
            if source == self.core / "chempotpy_kinds.f90":
                continue
            if source.parent.name == "save":
                continue
            with self.subTest(source=str(source.relative_to(self.repository))):
                self.assertIsNone(
                    hard_wired_precision.search(self.executable_fortran(
                        source.read_text(errors="replace")
                    ))
                )

    def test_internal_linear_algebra_uses_working_precision( self ) -> None:
        blas = (self.core / "chempotpy_blas.f90").read_text().lower()
        dsyev = (self.core / "chempotpy_dsyev.f90").read_text().lower()
        self.assertIn("real(wp)", blas)
        self.assertIn("real(wp)", dsyev)
        self.assertNotIn("double precision", self.executable_fortran(blas))
        self.assertNotIn("real(8)", self.executable_fortran(dsyev))


if __name__ == "__main__":
    unittest.main()
