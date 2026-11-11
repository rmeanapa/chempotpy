"""Regression tests for non-network PES parameters embedded in source."""

import hashlib
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PROVENANCE_ROOT = REPOSITORY_ROOT / "tests/data/external_parameters"


class EmbeddedPESParameterTests( unittest.TestCase ):
    def test_small_parameter_loaders_are_removed( self ) -> None:
        surfaces = (
            "chempotpy/CH4F/CH4F_GEN_GC_1996.f90",
            "chempotpy/CH4F/CH4F_GEN_GCSO_1996.f90",
            "chempotpy/CH4Br/CH4Br_VBMM_2002.f90",
            "chempotpy/F2H2/F2H2_GEN_JBKKL_1990.f90",
            "chempotpy/HOBr/HOBr_3App_PIP.f90",
            "chempotpy/COCO/COCO_PIPNN.f90",
            "chempotpy/C3H7NO/C3H7NO_PIP.f90",
        )
        forbidden = ("trim(path)", "open(", "read(", "close(")
        for relative_path in surfaces:
            with self.subTest(surface=relative_path):
                text = (REPOSITORY_ROOT / relative_path).read_text().lower()
                for token in forbidden:
                    self.assertNotIn(token, text)

    def test_released_files_are_provenance_only( self ) -> None:
        cases = {
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
        }
        for relative_path, expected_hash in cases.items():
            with self.subTest(data=relative_path):
                path = PROVENANCE_ROOT / relative_path
                actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
                self.assertEqual(actual_hash, expected_hash)
                self.assertFalse(
                    (REPOSITORY_ROOT / "chempotpy" / relative_path).exists()
                )


if __name__ == "__main__":
    unittest.main()
