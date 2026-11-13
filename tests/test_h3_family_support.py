"""Compile and exercise support shared by the H3 2002 fit family."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class H3FamilySupportTests( unittest.TestCase ):
    def test_h3_2002_basis_contract( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        basis = repository / "chempotpy" / "H3" / "chempotpy_h3_2002_basis.f90"
        driver_source = """\
program test_h3_2002_basis
    use chempotpy_h3_2002_basis, only: get_h3_basis_term, primary_term_count, &
        secondary_term_count
    implicit none
    integer :: coefficient_index, first_power, second_power, third_power
    if( primary_term_count /= 418 .or. secondary_term_count /= 330 ) &
        error stop 'H3 2002 basis term counts changed'
    call get_h3_basis_term(1, first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [0,1,1,1]) ) &
        error stop 'first primary H3 2002 basis term changed'
    call get_h3_basis_term(primary_term_count, first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [4,4,4,89]) ) &
        error stop 'last primary H3 2002 basis term changed'
    call get_h3_basis_term(primary_term_count+1, first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [0,1,1,1]) ) &
        error stop 'first secondary H3 2002 basis term changed'
    call get_h3_basis_term(primary_term_count+secondary_term_count, &
        first_power, second_power, third_power, coefficient_index)
    if( any([first_power, second_power, third_power, coefficient_index] /= [4,4,3,71]) ) &
        error stop 'last secondary H3 2002 basis term changed'
end program test_h3_2002_basis
"""
        with tempfile.TemporaryDirectory(prefix="chempotpy-h3-2002-basis-") as tmp:
            directory = Path(tmp)
            driver = directory / "driver.f90"
            executable = directory / "test_h3_2002_basis"
            driver.write_text(driver_source)
            subprocess.run(
                [compiler, str(basis), str(driver), "-o", str(executable)],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            )
            subprocess.run([str(executable)], check=True, cwd=directory)


if __name__ == "__main__":
    unittest.main()
