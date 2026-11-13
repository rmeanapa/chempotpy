"""Regression tests for the first fully formatted Aln PES migration."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.build_pes import core_fortran_dependencies


class AlnHalPModernizationTests( unittest.TestCase ):
    def test_codata_energies_remain_close_to_released_references(
        self,
    ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        source = repository / "chempotpy" / "Aln" / "Aln_HalP.f90"
        driver_source = """
program test_aln_halp
    implicit none
    real(8) :: coordinates(10000,3)
    real(8) :: energy(1), gradient(1,4,3), coupling(1,1,4,3)

    coordinates = 0.0d0
    coordinates(2,1) = 2.8d0
    call pes( coordinates, 2, 0, energy, gradient, coupling )
    if( abs(energy(1)-(-0.34621314477322201d0)) > &
        max(1.0d-8, 1.0d-6*abs(-0.34621314477322201d0)) ) &
        error stop 'two-atom CODATA migration envelope exceeded'

    coordinates(3,2) = 2.8d0
    coordinates(4,3) = 2.8d0
    call pes( coordinates, 4, 0, energy, gradient, coupling )
    if( abs(energy(1)-(-1.4003824541328054d0)) > &
        max(1.0d-8, 1.0d-6*abs(-1.4003824541328054d0)) ) &
        error stop 'four-atom CODATA migration envelope exceeded'
end program test_aln_halp
"""
        with tempfile.TemporaryDirectory(prefix="chempotpy-aln-halp-") as tmp:
            directory = Path(tmp)
            driver = directory / "driver.f90"
            executable = directory / "test_aln_halp"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    *(str(path) for path in core_fortran_dependencies(
                        source.read_text()
                    )),
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
            subprocess.run([str(executable)], check=True)


if __name__ == "__main__":
    unittest.main()
