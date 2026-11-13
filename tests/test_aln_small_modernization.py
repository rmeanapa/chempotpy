"""Regression tests for the first batch of compact modern Aln kernels."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.build_pes import core_fortran_dependencies


RELEASED_ENERGIES = {
    "Aln_ER": (
        -0.7698150904265624,
        -1.8585626930086288,
        -4.36487998234856,
        -9.899897424732698,
    ),
    "Aln_ER2": (
        -1.5177589431030127,
        -3.610192253857711,
        -9.383220387098891,
        -20.80649466306326,
    ),
    "Aln_Gol": (
        -1.7456424732159792,
        -3.079648292229208,
        -6.4297045279892515,
        -12.147978286279777,
    ),
    "Aln_NP_D": (
        -1.9747498268404122,
        -3.321171799040598,
        -6.854724456533508,
        -12.64465862729535,
    ),
    "Aln_PetW": (
        0.14996324177805817,
        0.31397835334423674,
        0.790047848038414,
        1.6427650991651668,
    ),
    "Aln_PapCEP": (
        -1.815653317463402,
        -3.2305200279029935,
        -6.634934172197606,
        -12.4587695892521,
    ),
    "Aln_PapKEP": (
        -1.530240413027754,
        -2.8393459180844216,
        -5.940632498167081,
        -11.43307465329597,
    ),
    "Aln_deSPH_M": (
        -0.5607281799658681,
        -1.1155275482075038,
        -2.652214003211972,
        -5.73132378183758,
    ),
    "Aln_reGol": (
        -1.6423363501715131,
        -2.9612272553080383,
        -6.401176072458937,
        -12.324116303968902,
    ),
    "Aln_reHalP": (
        -0.028077372087186626,
        -0.06996760220197504,
        -0.16736109217422857,
        -0.40398372220376116,
    ),
    "Aln_StrM": (
        -2.4231260277958606,
        -3.8984931677800816,
        -7.85046726302588,
        -14.127221181464625,
    ),
    "Aln_SutC": (
        -2.8292340903668514,
        -4.433086373337149,
        -8.909824439712663,
        -15.542661782548183,
    ),
    "Aln_rePetW": (
        -0.9567658356544403,
        -1.4058305128968,
        -2.8574101028236,
        -5.400213693836597,
    ),
    "Aln_reStrM": (
        -1.6450083825861372,
        -2.949673148688903,
        -6.216490163089939,
        -11.779637198886185,
    ),
    "Aln_reSutC": (
        -2.1557213713508308,
        -3.6542847461215704,
        -7.444284214968483,
        -13.623390708513137,
    ),
    "Aln_redeSPH_M": (
        -0.6876882059594223,
        -2.344934541892896,
        -5.494214075909317,
        -13.147261528830326,
    ),
}


class AlnSmallModernizationTests( unittest.TestCase ):
    def test_codata_energies_remain_close_to_released_references(
        self,
    ) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]

        with tempfile.TemporaryDirectory(prefix="chempotpy-aln-small-") as tmp:
            directory = Path(tmp)
            for surface_name, expected in RELEASED_ENERGIES.items():
                with self.subTest(surface=surface_name):
                    source = (
                        repository
                        / "chempotpy"
                        / "Aln"
                        / (surface_name + ".f90")
                    )
                    driver = directory / (surface_name + "_driver.f90")
                    executable = directory / (surface_name + "_test")
                    driver.write_text(self._driver_source(expected))
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

    @staticmethod
    def _driver_source( expected: tuple[float, ...] ) -> str:
        checks = "\n".join(
            """    call pes( coordinates, {count}, 0, energy, gradient, coupling )
    if( abs(energy(1)-({energy:.17g}d0)) > &
        max(1.0d-8, 1.0d-6*abs({energy:.17g}d0)) ) &
        error stop 'CODATA migration envelope exceeded for {count} atoms'""".format(
                count=count,
                energy=energy,
            )
            for count, energy in zip((2, 3, 5, 8), expected)
        )
        return """
program test_aln_surface
    implicit none
    real(8) :: coordinates(10000,3)
    real(8) :: energy(1), gradient(1,8,3), coupling(1,1,8,3)

    coordinates = 0.0d0
    coordinates(2,:) = [2.8d0, 0.2d0, 0.1d0]
    coordinates(3,:) = [0.4d0, 3.1d0, 0.3d0]
    coordinates(4,:) = [0.2d0, 0.6d0, 3.3d0]
    coordinates(5,:) = [3.0d0, 3.2d0, 0.4d0]
    coordinates(6,:) = [3.2d0, 0.5d0, 3.4d0]
    coordinates(7,:) = [0.7d0, 3.5d0, 3.6d0]
    coordinates(8,:) = [3.7d0, 3.8d0, 3.9d0]
{checks}
end program test_aln_surface
""".format(checks=checks)


if __name__ == "__main__":
    unittest.main()
