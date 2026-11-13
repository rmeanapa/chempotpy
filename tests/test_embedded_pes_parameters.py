"""Regression tests for non-network PES parameters embedded in source."""

import hashlib
import re
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PROVENANCE_ROOT = REPOSITORY_ROOT / "tests/data/external_parameters"


class EmbeddedPESParameterTests( unittest.TestCase ):
    def test_small_parameter_loaders_are_removed( self ) -> None:
        surfaces = (
            "chempotpy/ArNO/ArNO_Ap_PIPNN.f90",
            "chempotpy/ArNO/ArNO_App_PIPNN.f90",
            "chempotpy/CH4F/CH4F_GEN_GC_1996.f90",
            "chempotpy/CH4F/CH4F_GEN_GCSO_1996.f90",
            "chempotpy/CH4Br/CH4Br_VBMM_2002.f90",
            "chempotpy/F2H2/F2H2_GEN_JBKKL_1990.f90",
            "chempotpy/HOBr/HOBr_3App_PIP.f90",
            "chempotpy/COCO/COCO_PIPNN.f90",
            "chempotpy/C3H7NO/C3H7NO_PIP.f90",
            "chempotpy/CH4/ch4_t8.f90",
            "chempotpy/ClH2O/ClH2O_PIPNN.f90",
            "chempotpy/H2CO/H2CO_PIP.f90",
            "chempotpy/H4O2/H4O2_PIP.f90",
            "chempotpy/N4/N4_BDTC_2014.f90",
        )
        forbidden = ("trim(path)", "open(", "read(", "close(")
        for relative_path in surfaces:
            with self.subTest(surface=relative_path):
                text = (REPOSITORY_ROOT / relative_path).read_text().lower()
                for token in forbidden:
                    self.assertNotIn(token, text)

    def test_released_files_are_provenance_only( self ) -> None:
        cases = {
            "ArNO/spline.dat": (
                "ec5c8c9faedefb83fb6f24c8716cccdd9279c87f244682063003582edfb37662"
            ),
            "CH4F/CONST-noso": (
                "3ca4951e559524a8954ae9fce036665a9c78abe96fbc269bdece3f477694053e"
            ),
            "CH4F/CONST-so": (
                "3f3be5eee6b95f97ce1a87240def86403527b9fda3a87b4aca9779f712ba29fe"
            ),
            "CH4Br/CONST.dat": (
                "a76be154bf20b50c18c1d0ccd572aa29f19769e16be54662c30cb2b3a36adc43"
            ),
            "F2H2/bjkkl.par": (
                "10e2c911d6975a0b8caae9aa55db283aec5990f287b306d57391d1a58856cb78"
            ),
            "HOBr/o_hbr_a_pp.par": (
                "ee5d8ff9a341a7c4cee6e2b725559711ba44ec72f81785f34b49f200ec476cf6"
            ),
            "COCO/nn_ococ_w20.txt": (
                "cd7436c1138e537b6da9b8a95a55e57badb8b45d4e0ec9b1d5014f3dbab722e5"
            ),
            "C3H7NO/coeff.dat": (
                "ca1ee08b7b8220a42e6a00f6d6dc05def5e9c9847530257fccc6bb98422d6f06"
            ),
            "CH3OH/coeff.dat": (
                "63183735d626d6b047bbfe29f0877487dad0d232f3c0526eaeea9456883286db"
            ),
            "CH4/ch4pes": (
                "c4809fdded194e955072fd2c52248e2117db66b56f479c9e88c016f0f7b303dd"
            ),
            "CH4_H2O_H2O/coeff_CH4-H2O-H2O.dat": (
                "1ef6f392deaade04eacc5dd731ab44359d2e858fa619f5bde0bf2ef19986b204"
            ),
            "CH5p/fullfit.inc": (
                "46ccddfe237392cb71ef940ebf4f3673d55e8ff38e97502bd60e5e0f06e58cfe"
            ),
            "ClH2O/spline_tables.txt": (
                "c58b95e9997b4dc9858bb876beaa594eb4e05b791041eace2512b1d699facf86"
            ),
            "H2CO/coeff8": (
                "cd9467aae6aee8ff319243b88526f629444e352698afdc342ff2219af4186200"
            ),
            "H2CO/index8": (
                "75e4b4c9e9d87f6c155a46e6114ebe75744c80ea636f74a38fe5651d0b1cef42"
            ),
            "H2CO/pothco_hmk.para": (
                "95fd4828349d93f67589972b905dadd793fe9f7aa52a659bdc764a5c1debca1b"
            ),
            "H4O2/h4o2.pes4.coeff.dat": (
                "83c0d3d63e05f390d86a31dd51020bd993f63eef95d249492c7f940da0052fa7"
            ),
            "HOCO_Ar/coeff_HOCO-Ar_under1.5_12432_C.dat": (
                "bf6e117091c98f9b3ba2fb1829b1ab874423aa4c2c5857750a13aa8a5daadcc2"
            ),
            "N4/pes.dat": (
                "5ca040dae649a2afc84e475a54ecf981cee9021a9c6faeb1f9afb57d21a2d1e5"
            ),
        }
        for relative_path, expected_hash in cases.items():
            with self.subTest(data=relative_path):
                path = PROVENANCE_ROOT / relative_path
                actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
                self.assertEqual(actual_hash, expected_hash)
                self.assertFalse(
                    (REPOSITORY_ROOT / "chempotpy" / relative_path).exists()
                )

    def test_catalog_has_no_runtime_file_io( self ) -> None:
        statement = re.compile(r"^\s*(open|read|close)\s*\(", re.IGNORECASE)
        violations = []
        for source in (REPOSITORY_ROOT / "chempotpy").rglob("*.f90"):
            for line_number, line in enumerate(source.read_text().splitlines(), 1):
                if line.lstrip().startswith("!"):
                    continue
                if statement.search(line):
                    violations.append(
                        "{}:{}:{}".format(
                            source.relative_to(REPOSITORY_ROOT),
                            line_number,
                            line.strip(),
                        )
                    )
        self.assertEqual(violations, [])


if __name__ == "__main__":
    unittest.main()
