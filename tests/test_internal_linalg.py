"""Compile and validate the self-contained DSYEV compatibility routine."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class InternalLinearAlgebraTests( unittest.TestCase ):
    def test_symmetric_eigenproblem( self ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        implementation = (
            repository
            / "chempotpy"
            / "core"
            / "chempotpy_dsyev.f90"
        )
        kinds = (
            repository
            / "chempotpy"
            / "core"
            / "chempotpy_kinds.f90"
        )
        driver_source = """
program test_dsyev
    implicit none
    real(8) :: matrix(5,5), original(5,5), eigenvalues(5), work(32)
    real(8) :: identity(5,5), residual(5,5)
    integer :: info, i
    work = 0.0d0
    call dsyev('V', 'U', 5, matrix, 5, eigenvalues, work, -1, info)
    if( info /= 0 .or. work(1) < 14.0d0 ) &
        error stop 'workspace query failed'

    matrix = reshape([ &
        4.0d0, 1.0d0, 0.5d0, 0.0d0, -0.2d0, &
        1.0d0, 3.0d0, 0.0d0, 0.3d0, 0.0d0, &
        0.5d0, 0.0d0, 2.0d0, 0.4d0, 0.1d0, &
        0.0d0, 0.3d0, 0.4d0, 1.0d0, 0.2d0, &
       -0.2d0, 0.0d0, 0.1d0, 0.2d0, 5.0d0], [5,5])
    original = matrix
    call dsyev('V', 'L', 5, matrix, 5, eigenvalues, work, 32, info)
    if( info /= 0 ) error stop 'five-state DSYEV failed'
    do i = 1, 4
        if( eigenvalues(i) > eigenvalues(i+1) ) &
            error stop 'eigenvalues not sorted'
    enddo
    identity = matmul(transpose(matrix), matrix)
    do i = 1, 5
        identity(i,i) = identity(i,i)-1.0d0
    enddo
    if( maxval(abs(identity)) > 1.0d-12 ) &
        error stop 'five-state eigenvectors not orthogonal'
    residual = matmul(original,matrix)
    do i = 1, 5
        residual(:,i) = residual(:,i)-eigenvalues(i)*matrix(:,i)
    enddo
    if( maxval(abs(residual)) > 1.0d-11 ) &
        error stop 'five-state eigenpair residual too large'

    matrix = 0.0d0
    matrix(1:2,1:2) = reshape( &
        [2.0d0, 1.0d0, 1.0d0, 2.0d0], [2,2])
    call dsyev('V', 'U', 2, matrix, 5, eigenvalues, work, 16, info)
    if( info /= 0 ) error stop 'DSYEV failed'
    if( maxval(abs(eigenvalues(1:2)-[1.0d0,3.0d0])) > 1.0d-12 ) &
        error stop 'eigenvalues differ'
    if( abs(dot_product(matrix(:,1),matrix(:,2))) > 1.0d-12 ) &
        error stop 'eigenvectors not orthogonal'
end program test_dsyev
"""
        with tempfile.TemporaryDirectory(prefix="chempotpy-dsyev-test-") as tmp:
            directory = Path(tmp)
            driver = directory / "driver.f90"
            executable = directory / "test_dsyev"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    str(kinds),
                    str(implementation),
                    str(driver),
                    "-o",
                    str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run([str(executable)], check=True)


if __name__ == "__main__":
    unittest.main()
