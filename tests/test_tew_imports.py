"""Independent release and derivative checks for the imported Tew surfaces."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class TewSurfaceImportTests(unittest.TestCase):
    CASES = {
        "CH2O2/CH2O2_Tew_2016": {
            "geometry": (
                (0.4135438432, 0.0, 0.0990841188),
                (1.5052953939, 0.0, 0.0329901771),
                (-0.2152558268, 0.0, 1.1196963142),
                (-0.1216759225, 0.0, -1.1319495094),
                (-1.0819397210, 0.0, -1.0183007056),
            ),
            "released_energy_ev": -4.804217027579388e-06,
            "released_gradient_ev_per_angstrom": (
                4.211633799893609e-04, 0.0, 1.022832445640923e-02,
                1.823934577457524e-03, 0.0, -1.459132469891710e-03,
                -1.457397096377681e-04, 0.0, -2.257646007628804e-03,
                -3.596185022708470e-05, 0.0, -7.793612141620695e-03,
                -2.063396397581855e-03, 0.0, 1.282066162731984e-03,
            ),
        },
        "HOPO/HOPO_Tew": {
            "geometry": (
                (0.0, -0.45920589, -0.09184438),
                (0.0, 0.31692522, 1.31251474),
                (0.0, 0.49167473, -1.21110024),
                (0.0, 1.27987369, 1.21307381),
            ),
            "released_energy_ev": -7.670313263014210e-06,
            "released_gradient_ev_per_angstrom": (
                0.0, 1.728640007713108e-07, -1.822647750601132e-07,
                0.0, 5.636398429849412e-08, -3.546178962596138e-08,
                0.0, -1.545594709022505e-07, 2.491886793416624e-07,
                0.0, -7.466851398914837e-08, -3.146211463328709e-08,
            ),
        },
        "H3FO/H3FO_Tew_2022": {
            "support": ("chempotpy_h3fo_tew_parameters.f90",),
            "geometry": (
                (0.0, -0.18003783278112614, 0.71363839708627863),
                (-0.76330172306824107, 0.23992429279980815, 1.1149221682055686),
                (0.76330172306824107, 0.23992429279980815, 1.1149221682055686),
                (0.0, -0.17237330414243346, -1.0056275056315507),
                (0.0, -0.12743744867605664, -1.9378552278658703),
            ),
            "released_energy_ev": -1.5233445954127221e-13,
            "released_gradient_ev_per_angstrom": (
                -1.1315130369843112e-13, 1.5776310025889554e-07, 1.8047676111954424e-07,
                2.6462737883651960e-07, -8.105269910099866e-08, -1.0657507108754492e-07,
                -2.6462726596397552e-07, -8.105263699615042e-08, -1.0657501238077604e-07,
                2.6760920290333790e-16, -1.7474436032562562e-08, 1.1254869968967817e-07,
                1.1150383454305745e-17, 2.1816671879178887e-08, -7.987537759736026e-08,
            ),
        },
    }

    def test_python_adapters_recognize_connectivity(self) -> None:
        from chempotpy.CH2O2 import check as check_formic
        from chempotpy.HOPO import check as check_hopo
        from chempotpy.H3FO import check as check_h3fo

        formic = [
            ["O", -0.1217, 0.0, -1.1319],
            ["H", -1.0819, 0.0, -1.0183],
            ["C", 0.4135, 0.0, 0.0991],
            ["O", -0.2153, 0.0, 1.1197],
            ["H", 1.5053, 0.0, 0.0330],
        ]
        _, formic_rank = check_formic("CH2O2", "CH2O2_Tew_2016", formic)
        self.assertEqual(formic_rank[:, 0].tolist(), ["C", "H", "O", "O", "H"])
        self.assertEqual(formic_rank[:, 1].tolist(), [2, 4, 3, 0, 1])

        hopo = [
            ["O", 0.0, 0.4917, -1.2111],
            ["H", 0.0, 1.2799, 1.2131],
            ["P", 0.0, -0.4592, -0.0918],
            ["O", 0.0, 0.3169, 1.3125],
        ]
        _, hopo_rank = check_hopo("HOPO", "HOPO_Tew", hopo)
        self.assertEqual(hopo_rank[:, 0].tolist(), ["P", "O", "O", "H"])
        self.assertEqual(hopo_rank[:, 1].tolist(), [2, 3, 0, 1])

        h3fo = [
            ["H", 0.0, -0.1724, -1.0056],
            ["F", 0.0, -0.1274, -1.9379],
            ["H", 0.7633, 0.2399, 1.1149],
            ["O", 0.0, -0.1800, 0.7136],
            ["H", -0.7633, 0.2399, 1.1149],
        ]
        _, h3fo_rank = check_h3fo("H3FO", "H3FO_Tew_2022", h3fo)
        self.assertEqual(h3fo_rank[:, 0].tolist(), ["O", "H", "H", "H", "F"])
        self.assertEqual(h3fo_rank[:, 1].tolist(), [3, 2, 4, 0, 1])

    def test_released_checkpoints_and_analytic_gradients(self) -> None:
        compiler = shutil.which("gfortran")
        if compiler is None:
            self.skipTest("gfortran is not available")
        repository = Path(__file__).resolve().parents[1]
        core = repository / "chempotpy" / "core"
        for selector, case in self.CASES.items():
            with self.subTest(surface=selector):
                system, surface = selector.split("/")
                source = repository / "chempotpy" / system / (surface + ".f90")
                support = [
                    repository / "chempotpy" / system / name
                    for name in case.get("support", ())
                ]
                natoms = len(case["geometry"])
                assignments = "\n".join(
                    "    coordinates({0},:) = [{1}]".format(
                        index,
                        ", ".join("{:.17e}_wp".format(value) for value in point),
                    )
                    for index, point in enumerate(case["geometry"], start=1)
                )
                driver_source = """
program test_tew_import
    use chempotpy_kinds, only: wp
    implicit none
    integer, parameter  :: natoms = {natoms}
    real(wp), parameter :: step = 1.0e-5_wp
    real(wp)            :: coupling(1,1,natoms,3), coordinates(natoms,3)
    real(wp)            :: displaced(natoms,3), energy(1), forward(1), backward(1)
    real(wp)            :: gradient(1,natoms,3), numerical
    integer             :: atom, direction
{assignments}
    call pes(coordinates, 1, energy, gradient, coupling)
    write(*,'(ES26.17E3)') energy(1)
    do atom=1, natoms
        do direction=1, 3
            displaced = coordinates
            displaced(atom,direction) = displaced(atom,direction)+step
            call pes(displaced, 0, forward, gradient, coupling)
            displaced(atom,direction) = coordinates(atom,direction)-step
            call pes(displaced, 0, backward, gradient, coupling)
            numerical = (forward(1)-backward(1))/(2.0_wp*step)
            call pes(coordinates, 1, energy, gradient, coupling)
            write(*,'(2(ES26.17E3,1X))') gradient(1,atom,direction), numerical
        enddo
    enddo
end program test_tew_import
""".format(natoms=natoms, assignments=assignments)
                with tempfile.TemporaryDirectory(prefix="chempotpy-tew-") as tmp:
                    directory = Path(tmp)
                    driver = directory / "driver.f90"
                    executable = directory / "test_tew"
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
                            str(core / "chempotpy_status.f90"),
                            str(core / "chempotpy_constants.f90"),
                            str(core / "chempotpy_units.f90"),
                            *(str(path) for path in support),
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
                    completed = subprocess.run(
                        [str(executable)], check=True, capture_output=True, text=True
                    )
                rows = [tuple(map(float, line.split())) for line in completed.stdout.splitlines()]
                self.assertAlmostEqual(
                    rows[0][0], case["released_energy_ev"], delta=1.0e-10
                )
                for row, released in zip(
                    rows[1:], case["released_gradient_ev_per_angstrom"]
                ):
                    self.assertAlmostEqual(row[0], released, delta=8.0e-8)
                    self.assertAlmostEqual(row[0], row[1], delta=2.0e-7)


if __name__ == "__main__":
    unittest.main()
