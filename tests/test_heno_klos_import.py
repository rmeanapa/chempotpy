"""Independent release checks for the modern Klos He--NO import."""

import json
import math
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


class HeNOKlosImportTests(unittest.TestCase):
    """Verify the complete angle-and-distance release checkpoint set."""

    def test_python_adapter_orders_he_n_o(self) -> None:
        """Put the rare-gas atom first and retain the N-to-O bond direction."""

        from chempotpy.HeNO import check

        geometry = [
            ["O", 0.0, 0.0, 0.6],
            ["He", -3.0, 0.0, 0.0],
            ["N", 0.0, 0.0, -0.6],
        ]
        coordinates, rank = check("HeNO", "HeNO_Klos_2000", geometry)
        self.assertEqual(rank[:, 0].tolist(), ["He", "N", "O"])
        self.assertEqual(rank[:, 1].tolist(), [1, 2, 0])
        np.testing.assert_allclose(coordinates, np.asarray([geometry[1][1:], geometry[2][1:], geometry[0][1:]]))

    def test_all_release_components(self) -> None:
        """Match all 20 Vsum and Vdif values sampled from the released evaluator."""

        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        fixture = json.loads(
            (repository / "tests" / "data" / "heno_klos_2000_release.json").read_text()
        )
        driver_source = """
program test_heno_klos_import
    use chempotpy_heno_klos_2000_pes, only: evaluate_components
    use chempotpy_kinds,               only: wp
    implicit none
    real(wp) :: cosine_theta, distance_bohr, v_difference, v_sum
    integer  :: input_status
    do
        read(*,*,iostat=input_status) distance_bohr, cosine_theta
        if( input_status /= 0 ) exit
        call evaluate_components(distance_bohr, cosine_theta, v_sum, v_difference)
        write(*,'(2(ES26.17E3,1X))') v_sum, v_difference
    enddo
end program test_heno_klos_import
"""
        core = repository / "chempotpy" / "core"
        source = repository / "chempotpy" / "HeNO" / "HeNO_Klos_2000.f90"
        with tempfile.TemporaryDirectory(prefix="chempotpy-heno-klos-") as temporary:
            directory = Path(temporary)
            driver = directory / "driver.f90"
            executable = directory / "test_heno_klos"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    "-std=f2018",
                    "-O2",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    "-Wno-unused-dummy-argument",
                    "-Wno-compare-reals",
                    "-ffree-line-length-none",
                    str(core / "chempotpy_kinds.f90"),
                    str(core / "chempotpy_constants.f90"),
                    str(core / "chempotpy_units.f90"),
                    str(core / "chempotpy_atomic_elements.f90"),
                    str(core / "chempotpy_geometry.f90"),
                    str(core / "chempotpy_internal_coordinates.f90"),
                    str(core / "chempotpy_dispersion.f90"),
                    str(source),
                    str(driver),
                    "-o",
                    str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            )
            input_rows = [
                "{:.17e} {:.17e}".format(
                    record["distance_bohr"],
                    math.cos(math.radians(record["angle_degree"])),
                )
                for record in fixture["records"]
            ]
            completed = subprocess.run(
                [str(executable)],
                input="\n".join(input_rows) + "\n",
                check=True,
                capture_output=True,
                text=True,
            )
        rows = [np.fromstring(line, sep=" ") for line in completed.stdout.splitlines()]
        self.assertEqual(len(rows), len(fixture["records"]))
        for row, record in zip(rows, fixture["records"]):
            with self.subTest(
                distance=record["distance_bohr"], angle=record["angle_degree"]
            ):
                np.testing.assert_allclose(
                    row,
                    [
                        record["vsum_inverse_centimeter"],
                        record["vdif_inverse_centimeter"],
                    ],
                    atol=2.0e-6,
                    rtol=1.1e-7,
                )


if __name__ == "__main__":
    unittest.main()
