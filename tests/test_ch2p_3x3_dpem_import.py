"""Independent release checks for the modern CH2+ three-state DPEM import."""

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from chempotpy.core.units import (
    EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR,
    EV_TO_HARTREE,
)


class CH2pDPEMImportTests(unittest.TestCase):
    """Verify the modern evaluator against the archived authors' release."""

    def test_python_adapter_orders_hydrogen_hydrogen_carbon(self) -> None:
        """Reorder arbitrary input into the H, H, C convention of the fit."""

        from chempotpy.CH2p import check

        geometry = [
            ["C", 1.5, 0.2, 0.0],
            ["H", -0.4, 0.0, 0.0],
            ["H", 0.4, 0.0, 0.0],
        ]
        coordinates, rank = check("CH2p", "CH2p_3x3_DPEM_2025", geometry)
        np.testing.assert_allclose(
            coordinates,
            np.asarray([geometry[1][1:], geometry[2][1:], geometry[0][1:]]),
        )
        self.assertEqual(rank[:, 1].tolist(), [1, 2, 0])

    def test_all_release_matrices_and_gradients(self) -> None:
        """Match nine independent matrices and all released analytic derivatives."""

        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        fixture = json.loads(
            (repository / "tests" / "data" / "ch2p_3x3_dpem_2025_release.json").read_text()
        )
        driver_source = """
program test_ch2p_3x3_dpem_import
    use chempotpy_ch2p_3x3_dpem_2025_pes, only: evaluate
    use chempotpy_kinds,                       only: wp
    implicit none
    real(wp) :: coordinates(3,3), matrix(3,3), matrix_gradient(3,3,3,3)
    integer  :: atom, direction, input_status, state_i, state_j
    do
        read(*,*,iostat=input_status) ((coordinates(atom,direction), direction=1,3), atom=1,3)
        if( input_status /= 0 ) exit
        call evaluate(coordinates, 1, matrix, matrix_gradient)
        write(*,'(90(ES26.17E3,1X))') ((matrix(state_i,state_j), state_j=1,3), state_i=1,3), &
            ((((matrix_gradient(state_i,state_j,atom,direction), direction=1,3), atom=1,3), &
            state_j=1,3), state_i=1,3)
    enddo
end program test_ch2p_3x3_dpem_import
"""
        core = repository / "chempotpy" / "core"
        system = repository / "chempotpy" / "CH2p"
        with tempfile.TemporaryDirectory(prefix="chempotpy-ch2p-dpem-") as temporary:
            directory = Path(temporary)
            driver = directory / "driver.f90"
            executable = directory / "test_ch2p_dpem"
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
                    str(core / "chempotpy_status.f90"),
                    str(core / "chempotpy_geometry.f90"),
                    str(core / "chempotpy_neural_network.f90"),
                    str(system / "chempotpy_ch2p_long_range.f90"),
                    str(system / "CH2p_3x3_DPEM_2025.f90"),
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
                " ".join(
                    "{:.17e}".format(value)
                    for value in np.asarray(record["geometry_angstrom"]).flat
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
            with self.subTest(geometry=record["name"]):
                matrix_hartree = row[:9].reshape(3, 3) * EV_TO_HARTREE
                gradient_hartree_per_bohr = (
                    row[9:].reshape(3, 3, 3, 3)
                    * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR
                )
                np.testing.assert_allclose(
                    matrix_hartree,
                    np.asarray(record["matrix_hartree"]),
                    atol=3.0e-13,
                    rtol=0.0,
                )
                np.testing.assert_allclose(
                    gradient_hartree_per_bohr,
                    np.asarray(record["matrix_gradient_hartree_per_bohr"]),
                    atol=3.0e-11,
                    rtol=0.0,
                )


if __name__ == "__main__":
    unittest.main()
