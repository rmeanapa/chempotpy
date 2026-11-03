"""Compile and exercise the bundled modern Fortran core modules."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class FortranCoreTests( unittest.TestCase ):
    def test_units_geometry_and_status_contract( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        core = repository / "chempotpy" / "core" / "fortran"
        source_names = (
            "chempotpy_kinds.f90",
            "chempotpy_units.f90",
            "chempotpy_legacy_units.f90",
            "chempotpy_geometry.f90",
            "chempotpy_status.f90",
        )
        driver_source = """
program test_fortran_core
    use chempotpy_geometry, only: cartesian_distance, &
        cosine_between_vectors, dot_product_3d, fill_distance_matrix
    use chempotpy_kinds, only: real64
    use chempotpy_legacy_units, only: legacy_bohr_to_angstrom, &
        legacy_coordinates_to_bohr, legacy_energy_to_electronvolt, &
        legacy_hartree_to_electronvolt
    use chempotpy_status, only: status_derivative_unavailable, &
        status_invalid_input, status_success
    use chempotpy_units, only: angstrom_to_bohr, bohr_to_angstrom, &
        electronvolt_to_hartree, hartree_to_electronvolt
    implicit none

    real(real64) :: first_point(3), second_point(3), coordinates(3,3)
    real(real64) :: distances(3,3), converted(1,3)

    if( bohr_to_angstrom /= 0.529177210903_real64 ) &
        error stop 'canonical length conversion changed'
    if( hartree_to_electronvolt /= 27.211386245988_real64 ) &
        error stop 'canonical energy conversion changed'
    if( abs(angstrom_to_bohr*bohr_to_angstrom-1.0_real64) &
        > epsilon(1.0_real64) ) error stop 'length conversion is inconsistent'
    if( abs(electronvolt_to_hartree*hartree_to_electronvolt-1.0_real64) &
        > epsilon(1.0_real64) ) error stop 'energy conversion is inconsistent'

    if( legacy_bohr_to_angstrom /= 0.52917718887329101562_real64 ) &
        error stop 'legacy length rounding changed'
    if( legacy_hartree_to_electronvolt /= &
        27.211385726928710938_real64 ) &
        error stop 'legacy energy rounding changed'
    call legacy_coordinates_to_bohr( &
        reshape([legacy_bohr_to_angstrom, 0.0_real64, 0.0_real64], [1,3]), &
        converted )
    if( converted(1,1) /= 1.0_real64 ) &
        error stop 'legacy coordinate conversion failed'
    if( legacy_energy_to_electronvolt(1.0_real64) /= &
        legacy_hartree_to_electronvolt ) &
        error stop 'legacy energy conversion failed'

    first_point = [0.0_real64, 0.0_real64, 0.0_real64]
    second_point = [3.0_real64, 4.0_real64, 12.0_real64]
    if( cartesian_distance(first_point,second_point) /= 13.0_real64 ) &
        error stop 'Cartesian distance failed'
    if( dot_product_3d([1.0_real64,2.0_real64,3.0_real64], &
        [4.0_real64,5.0_real64,6.0_real64]) /= 32.0_real64 ) &
        error stop 'dot product failed'
    if( abs(cosine_between_vectors([1.0_real64,0.0_real64,0.0_real64], &
        [0.0_real64,1.0_real64,0.0_real64])) > epsilon(1.0_real64) ) &
        error stop 'angle cosine failed'
    coordinates = 0.0_real64
    coordinates(2,1) = 3.0_real64
    coordinates(3,2) = 4.0_real64
    call fill_distance_matrix( coordinates, distances )
    if( distances(1,2) /= 3.0_real64 &
        .or. distances(1,3) /= 4.0_real64 &
        .or. distances(2,3) /= 5.0_real64 ) &
        error stop 'distance matrix failed'
    if( any(distances /= transpose(distances)) ) &
        error stop 'distance matrix is not symmetric'

    if( status_success /= 0 .or. status_derivative_unavailable /= 1 &
        .or. status_invalid_input /= 2 ) error stop 'status contract changed'
end program test_fortran_core
"""
        with tempfile.TemporaryDirectory(prefix="chempotpy-fortran-core-") as tmp:
            directory = Path(tmp)
            driver = directory / "driver.f90"
            executable = directory / "test_fortran_core"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    *(str(core / name) for name in source_names),
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
