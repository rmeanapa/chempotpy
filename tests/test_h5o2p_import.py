"""Independent release and derivative checks for H5O2p_PIP_2005."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class H5O2pImportTests(unittest.TestCase):
    RELEASED_ENERGY_HARTREE = -153.01000221
    RELEASED_GRADIENT_HARTREE_PER_BOHR = (
        0.00316817, -0.00000000, -0.00022061,
        -0.00316817, -0.00000000, 0.00022061,
        0.00000000, 0.00000000, 0.00000000,
        -0.00094431, -0.00015328, 0.00029557,
        -0.00094431, 0.00015328, 0.00029557,
        0.00094431, 0.00015328, -0.00029557,
        0.00094431, -0.00015328, -0.00029557,
    )

    def test_python_adapter_canonicalizes_oxygen_then_hydrogen(self) -> None:
        import numpy as np

        from chempotpy.H5O2p import check, reverse_order

        geometry = [
            ["H", 3.0, 0.0, 0.0],
            ["O", 1.0, 0.0, 0.0],
            ["H", 4.0, 0.0, 0.0],
            ["O", 2.0, 0.0, 0.0],
            ["H", 5.0, 0.0, 0.0],
            ["H", 6.0, 0.0, 0.0],
            ["H", 7.0, 0.0, 0.0],
        ]
        coordinates, rank = check("H5O2p", "H5O2p_PIP_2005", geometry)
        self.assertEqual(coordinates[:,0].tolist(), [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        gradient = np.arange(21, dtype=float).reshape(1,7,3)
        coupling = np.arange(21, dtype=float).reshape(1,1,7,3)
        _, reordered_gradient, reordered_coupling = reverse_order(
            np.zeros(1), gradient, coupling, geometry, rank
        )
        self.assertEqual(reordered_gradient[0,0,:].tolist(), gradient[0,2,:].tolist())
        self.assertEqual(reordered_gradient[0,1,:].tolist(), gradient[0,0,:].tolist())
        self.assertEqual(reordered_coupling[0,0,0,:].tolist(), coupling[0,0,2,:].tolist())

    def test_released_checkpoint_and_analytic_gradient(self) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        core = repository / "chempotpy" / "core"
        system = repository / "chempotpy" / "H5O2p"
        driver_source = """
program test_h5o2p_import
    use chempotpy_kinds, only: wp
    use chempotpy_units, only: hartree_per_bohr_to_ev_per_angstrom, hartree_to_electronvolt
    implicit none
    real(wp), parameter :: step = 1.0e-5_wp
    real(wp) :: coupling(1,1,7,3), energy(1), forward_energy(1), backward_energy(1)
    real(wp) :: gradient(1,7,3), coordinates(7,3), displaced(7,3), numerical_gradient(7,3)
    integer  :: atom, direction
    coordinates(1,:) = [-1.1962131_wp, 0.0_wp, 0.0_wp]
    coordinates(2,:) = [1.1962131_wp, 0.0_wp, 0.0_wp]
    coordinates(3,:) = [0.0_wp, 0.0_wp, 0.0_wp]
    coordinates(4,:) = [-1.6321110_wp, 0.7820072_wp, 0.3711158_wp]
    coordinates(5,:) = [-1.6321110_wp, -0.7820072_wp, 0.3711158_wp]
    coordinates(6,:) = [1.6321110_wp, -0.7820072_wp, -0.3711158_wp]
    coordinates(7,:) = [1.6321110_wp, 0.7820072_wp, -0.3711158_wp]
    call pes(coordinates, 1, energy, gradient, coupling)
    do atom=1, 7
        do direction=1, 3
            displaced                 = coordinates
            displaced(atom,direction) = displaced(atom,direction)+step
            call pes(displaced, 0, forward_energy, gradient, coupling)
            displaced(atom,direction) = coordinates(atom,direction)-step
            call pes(displaced, 0, backward_energy, gradient, coupling)
            numerical_gradient(atom,direction) = (forward_energy(1)-backward_energy(1))/(2.0_wp*step)
        enddo
    enddo
    call pes(coordinates, 1, energy, gradient, coupling)
    write(*,'(ES26.17E3)') energy(1)/hartree_to_electronvolt
    do atom=1, 7
        do direction=1, 3
            write(*,'(2(ES26.17E3,1X))') &
                gradient(1,atom,direction)/hartree_per_bohr_to_ev_per_angstrom, &
                numerical_gradient(atom,direction)
        enddo
    enddo
end program test_h5o2p_import
"""
        with tempfile.TemporaryDirectory(prefix="chempotpy-h5o2p-") as tmp:
            directory = Path(tmp)
            driver = directory / "driver.f90"
            executable = directory / "test_h5o2p"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    "-std=f2018",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    "-ffree-line-length-none",
                    str(core / "chempotpy_kinds.f90"),
                    str(core / "chempotpy_geometry.f90"),
                    str(core / "chempotpy_status.f90"),
                    str(core / "chempotpy_constants.f90"),
                    str(core / "chempotpy_units.f90"),
                    str(system / "chempotpy_h5o2p_pip_parameters.f90"),
                    str(system / "chempotpy_h5o2p_pip_basis.f90"),
                    str(system / "H5O2p_PIP_2005.f90"),
                    str(driver),
                    "-o",
                    str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=directory,
            )
            completed = subprocess.run(
                [str(executable)],
                check=True,
                capture_output=True,
                text=True,
            )
        rows = [tuple(map(float, line.split())) for line in completed.stdout.splitlines()]
        self.assertAlmostEqual(rows[0][0], self.RELEASED_ENERGY_HARTREE, delta=3.0e-9)
        self.assertEqual(len(rows[1:]), 21)
        for row, released_gradient in zip(rows[1:], self.RELEASED_GRADIENT_HARTREE_PER_BOHR):
            self.assertAlmostEqual(row[0], released_gradient, delta=5.1e-9)
            self.assertAlmostEqual(
                row[1],
                row[0]*51.4220674763259,
                delta=3.0e-4,
            )


if __name__ == "__main__":
    unittest.main()
