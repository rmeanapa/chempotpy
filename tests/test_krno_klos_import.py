"""Independent release checks for the modern Klos Kr--NO import."""

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


class KrNOKlosImportTests(unittest.TestCase):
    """Verify the complete release checkpoint grid and Python adapter."""

    def test_python_adapter_orders_kr_n_o(self) -> None:
        """Put krypton first and retain the N-to-O bond direction."""

        from chempotpy.KrNO import check

        geometry = [
            ["O", 0.0, 0.0, 0.6],
            ["Kr", -3.0, 0.0, 0.0],
            ["N", 0.0, 0.0, -0.6],
        ]
        coordinates, rank = check("KrNO", "KrNO_Klos_2009", geometry)
        self.assertEqual(rank[:, 0].tolist(), ["Kr", "N", "O"])
        self.assertEqual(rank[:, 1].tolist(), [1, 2, 0])
        np.testing.assert_allclose(
            coordinates,
            np.asarray([geometry[1][1:], geometry[2][1:], geometry[0][1:]]),
        )

    def test_all_release_components(self) -> None:
        """Match all 49 Vsum and Vdif values sampled from the released evaluator."""

        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        fixture = json.loads(
            (repository / "tests" / "data" / "krno_klos_2009_release.json").read_text()
        )
        driver_source = """
program test_krno_klos_import
    use chempotpy_kinds,              only: wp
    use chempotpy_krno_klos_2009_pes, only: evaluate_components
    implicit none
    real(wp) :: angle_degree, cosine_theta, distance_bohr, v_difference, v_sum
    integer  :: input_status
    do
        read(*,*,iostat=input_status) distance_bohr, angle_degree
        if( input_status /= 0 ) exit
        cosine_theta = cos(angle_degree*acos(-1.0_wp)/180.0_wp)
        call evaluate_components(distance_bohr, cosine_theta, v_sum, v_difference)
        write(*,'(2(ES26.17E3,1X))') v_sum, v_difference
    enddo
end program test_krno_klos_import
"""
        core = repository / "chempotpy" / "core"
        source = repository / "chempotpy" / "KrNO" / "KrNO_Klos_2009.f90"
        with tempfile.TemporaryDirectory(prefix="chempotpy-krno-klos-") as temporary:
            directory = Path(temporary)
            driver = directory / "driver.f90"
            executable = directory / "test_krno_klos"
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
                    "-ffree-line-length-none",
                    str(core / "chempotpy_kinds.f90"),
                    str(core / "chempotpy_constants.f90"),
                    str(core / "chempotpy_units.f90"),
                    str(core / "chempotpy_atomic_elements.f90"),
                    str(core / "chempotpy_geometry.f90"),
                    str(core / "chempotpy_internal_coordinates.f90"),
                    str(core / "chempotpy_special_functions.f90"),
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
                "{:.17e} {:.17e}".format(record["distance_bohr"], record["angle_degree"])
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
        expected = [
            [record["vsum_inverse_centimeter"], record["vdif_inverse_centimeter"]]
            for record in fixture["records"]
        ]
        self.assertEqual(len(rows), len(expected))
        np.testing.assert_allclose(np.asarray(rows), np.asarray(expected), atol=1.0e-7, rtol=2.0e-10)


if __name__ == "__main__":
    unittest.main()
