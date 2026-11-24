"""Independent release-checkpoint tests for the fixed-slab CH4/Ni adapter."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class CH4NiAdapterTests(unittest.TestCase):
    GEOMETRIES_BOHR = (
        (
            (25.89063, 25.82544, 28.06512),
            (24.21066, 26.80046, 25.30616),
            (27.58478, 26.79587, 25.31574),
            (25.89579, 23.87490, 25.31191),
            (25.89579, 25.82436, 26.00000),
        ),
        (
            (25.89579, 25.82436, -1.897653),
            (24.18331, 26.81306, 4.147893),
            (27.60826, 26.81305, 4.147894),
            (25.89579, 23.84697, 4.147894),
            (25.89579, 25.82436, 3.331697),
        ),
        (
            (25.89579, 25.82436, -0.1038303),
            (24.17093, 26.82021, 4.265322),
            (27.62065, 26.82021, 4.265322),
            (25.89579, 23.83266, 4.265322),
            (25.89579, 25.82436, 3.452416),
        ),
        (
            (25.59756, 26.03049, 2.190675),
            (25.75269, 26.50924, 5.435653),
            (28.03949, 28.46122, 3.034229),
            (28.54979, 24.67531, 3.382665),
            (27.11911, 26.45022, 3.785052),
        ),
        (
            (25.91800, 25.84158, 1.779192),
            (25.83421, 25.50006, 4.912542),
            (27.84865, 28.36447, 4.273411),
            (29.25709, 25.10935, 4.085836),
            (27.50351, 26.30543, 3.795846),
        ),
        (
            (25.89226, 25.82436, 1.817841),
            (26.73501, 25.95380, 4.158997),
            (28.45482, 28.93087, 4.167944),
            (30.16104, 25.96076, 4.156673),
            (28.45340, 26.95821, 3.343454),
        ),
    )
    RELEASED_COMPONENTS_HARTREE = (
        (-1.7594223550077004e-5, 1.7646716880232551e-5),
        (1.9421758751136053e-2, 6.5940727504787341e-3),
        (5.1886733836070675e-2, 6.2008948652437440e-3),
        (-1.1013728212375895e-2, 1.1774060521582266e-2),
        (8.6729064300499203e-3, 1.3627036774743332e-2),
        (3.2840793669635104e-3, 6.6481739534724055e-3),
    )

    def test_all_six_released_component_checkpoints(self) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        core = repository / "chempotpy" / "core"
        system = repository / "chempotpy" / "CH4Ni"
        driver_source = """
program test_ch4ni_release
    use chempotpy_ch4ni_runtime, only: evaluate_ch4ni_energy
    use chempotpy_kinds,         only: wp
    use chempotpy_units,         only: bohr_to_angstrom
    implicit none
    real(wp) :: components(2), coordinates_bohr(5,3), energy
    integer  :: atom, geometry, geometry_count
    read(*,*) geometry_count
    do geometry=1, geometry_count
        do atom=1, 5
            read(*,*) coordinates_bohr(atom,:)
        enddo
        call evaluate_ch4ni_energy(coordinates_bohr*bohr_to_angstrom, energy, components)
        write(*,'(3(ES26.17E3,1X))') energy, components
    enddo
end program test_ch4ni_release
"""
        input_lines = [str(len(self.GEOMETRIES_BOHR))]
        for geometry in self.GEOMETRIES_BOHR:
            input_lines.extend(
                " ".join("{:.16e}".format(value) for value in atom)
                for atom in geometry
            )
        with tempfile.TemporaryDirectory(prefix="chempotpy-ch4ni-") as tmp:
            directory = Path(tmp)
            driver = directory / "driver.f90"
            executable = directory / "test_ch4ni"
            driver.write_text(driver_source)
            subprocess.run(
                [
                    compiler,
                    "-std=f2018",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                    str(core / "chempotpy_kinds.f90"),
                    str(core / "chempotpy_constants.f90"),
                    str(core / "chempotpy_units.f90"),
                    str(system / "chempotpy_ch4ni_ni111_parameters.f90"),
                    str(system / "chempotpy_ch4ni_runtime.f90"),
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
                input="\n".join(input_lines) + "\n",
                check=True,
                capture_output=True,
                text=True,
            )
        rows = [tuple(map(float, line.split())) for line in completed.stdout.splitlines()]
        self.assertEqual(len(rows), len(self.RELEASED_COMPONENTS_HARTREE))
        for actual, expected_components in zip(rows, self.RELEASED_COMPONENTS_HARTREE):
            self.assertAlmostEqual(actual[0], actual[1] + actual[2], places=14)
            self.assertLessEqual(abs(actual[1] - expected_components[0]), 3.0e-6)
            self.assertLessEqual(abs(actual[2] - expected_components[1]), 6.0e-7)


if __name__ == "__main__":
    unittest.main()
