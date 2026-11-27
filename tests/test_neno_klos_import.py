"""Independent release checks for the modern Klos Ne--NO imports."""

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


class NeNOKlosImportTests(unittest.TestCase):
    """Verify both vibrational averages against the complete checkpoint grid."""

    def test_python_adapter_orders_ne_n_o(self) -> None:
        """Put the rare-gas atom first and retain the N-to-O bond direction."""

        from chempotpy.NeNO import check

        geometry = [
            ["O", 0.0, 0.0, 0.6],
            ["Ne", -3.0, 0.0, 0.0],
            ["N", 0.0, 0.0, -0.6],
        ]
        coordinates, rank = check("NeNO", "NeNO_Klos_2016_v0", geometry)
        self.assertEqual(rank[:, 0].tolist(), ["Ne", "N", "O"])
        self.assertEqual(rank[:, 1].tolist(), [1, 2, 0])
        np.testing.assert_allclose(
            coordinates,
            np.asarray([geometry[1][1:], geometry[2][1:], geometry[0][1:]]),
        )

    def test_all_release_components(self) -> None:
        """Match 70 Vsum and Vdif values sampled from both released evaluators."""

        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        fixture = json.loads(
            (repository / "tests" / "data" / "neno_klos_2016_release.json").read_text()
        )
        driver_source = """
program test_neno_klos_import
    use chempotpy_kinds,        only: wp
    use chempotpy_neno_runtime, only: evaluate_neno_vibrational_components
    implicit none
    real(wp) :: angle_degree, distance_bohr, v_difference, v_sum
    integer  :: input_status, vibrational_index
    do
        read(*,*,iostat=input_status) distance_bohr, angle_degree, vibrational_index
        if( input_status /= 0 ) exit
        call evaluate_neno_vibrational_components(distance_bohr, angle_degree, vibrational_index, &
            v_sum, v_difference)
        write(*,'(2(ES26.17E3,1X))') v_sum, v_difference
    enddo
end program test_neno_klos_import
"""
        core = repository / "chempotpy" / "core"
        source = repository / "chempotpy" / "NeNO"
        with tempfile.TemporaryDirectory(prefix="chempotpy-neno-klos-") as temporary:
            directory = Path(temporary)
            driver = directory / "driver.f90"
            executable = directory / "test_neno_klos"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    "-std=f2018",
                    "-O2",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    "-Wno-compare-reals",
                    "-ffree-line-length-none",
                    str(core / "chempotpy_kinds.f90"),
                    str(core / "chempotpy_constants.f90"),
                    str(core / "chempotpy_interpolation.f90"),
                    str(core / "chempotpy_quadrature.f90"),
                    str(core / "chempotpy_special_functions.f90"),
                    str(source / "chempotpy_neno_parameters.f90"),
                    str(source / "chempotpy_neno_runtime.f90"),
                    str(driver),
                    "-o",
                    str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            )
            input_rows = []
            expected_rows = []
            for vibrational_index, state in enumerate(("v0", "v2"), start=1):
                for record in fixture["states"][state]:
                    input_rows.append(
                        "{:.17e} {:.17e} {}".format(
                            record["distance_bohr"],
                            record["angle_degree"],
                            vibrational_index,
                        )
                    )
                    expected_rows.append(
                        [
                            record["vsum_inverse_centimeter"],
                            record["vdif_inverse_centimeter"],
                        ]
                    )
            completed = subprocess.run(
                [str(executable)],
                input="\n".join(input_rows) + "\n",
                check=True,
                capture_output=True,
                text=True,
            )
        rows = [np.fromstring(line, sep=" ") for line in completed.stdout.splitlines()]
        self.assertEqual(len(rows), len(expected_rows))
        np.testing.assert_allclose(
            np.asarray(rows),
            np.asarray(expected_rows),
            atol=1.0e-5,
            rtol=3.0e-8,
        )


if __name__ == "__main__":
    unittest.main()
