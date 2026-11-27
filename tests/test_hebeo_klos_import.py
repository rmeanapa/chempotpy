"""Independent release checks for the modern Klos He--BeO import."""

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


class HeBeOKlosImportTests(unittest.TestCase):
    """Verify release checkpoints, the RKHS core, and Python ordering."""

    def test_python_adapter_orders_he_be_o(self) -> None:
        """Put helium first and retain the Be-to-O bond direction."""

        from chempotpy.HeBeO import check

        geometry = [
            ["O", 0.0, 0.0, 0.7],
            ["He", -3.0, 0.0, 0.0],
            ["Be", 0.0, 0.0, -0.7],
        ]
        coordinates, rank = check("HeBeO", "HeBeO_Klos_2013", geometry)
        self.assertEqual(rank[:, 0].tolist(), ["He", "Be", "O"])
        self.assertEqual(rank[:, 1].tolist(), [1, 2, 0])
        np.testing.assert_allclose(
            coordinates,
            np.asarray([geometry[1][1:], geometry[2][1:], geometry[0][1:]]),
        )

    def test_release_checkpoint_grid(self) -> None:
        """Match 56 points sampled from the release's complete output grid."""

        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        fixture = json.loads(
            (repository / "tests" / "data" / "hebeo_klos_2013_release.json").read_text()
        )
        driver_source = """
program test_hebeo_klos_import
    use chempotpy_hebeo_klos_2013_pes, only: evaluate_fit
    use chempotpy_kinds,                only: wp
    implicit none
    real(wp) :: angle_degree, cosine_theta, distance_bohr, energy_inverse_centimeter
    integer  :: input_status
    do
        read(*,*,iostat=input_status) distance_bohr, angle_degree
        if( input_status /= 0 ) exit
        cosine_theta = cos(angle_degree*acos(-1.0_wp)/180.0_wp)
        call evaluate_fit(distance_bohr, cosine_theta, energy_inverse_centimeter)
        write(*,'(ES26.17E3)') energy_inverse_centimeter
    enddo
end program test_hebeo_klos_import
"""
        core = repository / "chempotpy" / "core"
        source = repository / "chempotpy" / "HeBeO" / "HeBeO_Klos_2013.f90"
        with tempfile.TemporaryDirectory(prefix="chempotpy-hebeo-klos-") as temporary:
            directory = Path(temporary)
            driver = directory / "driver.f90"
            executable = directory / "test_hebeo_klos"
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
                    str(core / "chempotpy_rkhs.f90"),
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
        actual = np.fromstring(completed.stdout, sep=" ")
        expected = np.asarray(
            [record["energy_inverse_centimeter"] for record in fixture["records"]]
        )
        self.assertEqual(actual.size, expected.size)
        np.testing.assert_allclose(actual, expected, atol=5.0e-4, rtol=2.0e-10)


if __name__ == "__main__":
    unittest.main()
