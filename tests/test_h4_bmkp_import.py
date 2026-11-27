"""Independent release checks for the modern BMKP H4 import."""

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from chempotpy.core.units import (
    BOHR_TO_ANGSTROM,
    EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR,
    EV_TO_HARTREE,
)


class H4BMKPImportTests(unittest.TestCase):
    """Verify every checkpoint distributed with the authoritative release."""

    def test_python_adapter_accepts_four_equivalent_hydrogens(self) -> None:
        """Keep equivalent hydrogen atoms in the caller's order."""

        from chempotpy.H4 import check

        geometry = [
            ["H", 0.0, 0.0, -0.5],
            ["H", 0.0, 0.0, 0.5],
            ["H", 0.0, 0.0, 3.1],
            ["H", 0.0, 0.0, -3.1],
        ]
        coordinates, rank = check("H4", "H4_BMKP_2002", geometry)
        np.testing.assert_allclose(coordinates, np.asarray([row[1:] for row in geometry]))
        self.assertEqual(rank[:, 1].tolist(), [0, 1, 2, 3])

    def test_all_release_energies_and_gradients(self) -> None:
        """Match all 22 energy and Cartesian-gradient release checkpoints."""

        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        fixture = json.loads(
            (repository / "tests" / "data" / "h4_bmkp_2002_release.json").read_text()
        )
        driver_source = """
program test_h4_bmkp_import
    use chempotpy_h4_bmkp_2002_pes, only: evaluate
    use chempotpy_kinds,             only: wp
    implicit none
    real(wp) :: coupling(1,1,4,3), coordinates(4,3), energy(1), gradient(1,4,3)
    integer  :: atom, direction, input_status
    do
        read(*,*,iostat=input_status) ((coordinates(atom,direction), direction=1,3), atom=1,4)
        if( input_status /= 0 ) exit
        call evaluate(coordinates, 1, energy, gradient, coupling)
        write(*,'(13(ES26.17E3,1X))') energy(1), &
            ((gradient(1,atom,direction), direction=1,3), atom=1,4)
    enddo
end program test_h4_bmkp_import
"""
        core = repository / "chempotpy" / "core"
        source = repository / "chempotpy" / "H4" / "H4_BMKP_2002.f90"
        with tempfile.TemporaryDirectory(prefix="chempotpy-h4-bmkp-") as temporary:
            directory = Path(temporary)
            driver = directory / "driver.f90"
            executable = directory / "test_h4_bmkp"
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
                    str(core / "chempotpy_status.f90"),
                    str(core / "chempotpy_constants.f90"),
                    str(core / "chempotpy_units.f90"),
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
            input_rows = []
            for record in fixture["records"]:
                coordinates = np.asarray(record["geometry_bohr"]) * BOHR_TO_ANGSTROM
                input_rows.append(" ".join("{:.17e}".format(value) for value in coordinates.flat))
            completed = subprocess.run(
                [str(executable)],
                input="\n".join(input_rows) + "\n",
                check=True,
                capture_output=True,
                text=True,
            )
        rows = [np.fromstring(line, sep=" ") for line in completed.stdout.splitlines()]
        self.assertEqual(len(rows), 22)
        for index, (row, record) in enumerate(zip(rows, fixture["records"]), start=1):
            with self.subTest(checkpoint=index):
                energy_hartree = row[0] * EV_TO_HARTREE
                gradient_hartree_per_bohr = (
                    row[1:].reshape(4, 3) * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHR
                )
                self.assertAlmostEqual(
                    energy_hartree,
                    record["energy_hartree"],
                    delta=3.0e-13,
                )
                np.testing.assert_allclose(
                    gradient_hartree_per_bohr,
                    np.asarray(record["gradient_hartree_per_bohr"]),
                    atol=3.0e-11,
                    rtol=0.0,
                )


if __name__ == "__main__":
    unittest.main()
