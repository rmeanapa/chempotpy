"""Cross-project regression for the CHEMPOTPY/POLYRATE CH5 J1 surface."""

import types
import unittest

import numpy as np

from chempotpy.CH5.adapters import create_ch5_j1_surface
from chempotpy.core import PESRequest
from chempotpy.core.units import (
    BOHR_TO_ANGSTROM,
)


class CH5J1AdapterTests( unittest.TestCase ):
    def test_matches_polyrate_atomic_unit_reference( self ) -> None:
        coordinates_bohr = np.array(
            [
                [1.20, 1.10, 1.30],
                [-1.10, -1.25, 1.20],
                [-1.25, 1.15, -1.10],
                [1.10, -1.20, -1.30],
                [0.00, 0.00, 0.00],
                [4.00, 0.30, -0.20],
            ]
        )
        raw_gradient = np.array(
            [
                -3.4058480084848027e-3,
                4.2720304431447631e-1,
                5.8535177503525937e-1,
                -1.7082691647960174e-1,
                6.2098687790046925e-1,
                -1.4593089327621185,
                -1.6672108409294806e-1,
                3.0177965618828911e-1,
                -7.7562127634199851e-1,
                -7.5633866689531107e-2,
                8.2912219501502660e-1,
                -1.1292562407883798e-1,
                3.3345148893240517e-1,
                -1.8139948456802346e-1,
                4.0623932295742948e-1,
                -2.9076594849912044e-2,
                -6.0880674859475903e-1,
                7.9592016122859890e-2,
            ]
        ).reshape((1, 6, 3), order="F")
        expected_gradient = np.array(
            [
                -6.6248287447924058e-5,
                -3.2422235718002892e-3,
                6.4845823614700874e-3,
                8.3077920076466208e-3,
                5.8686971834779576e-3,
                -3.5276744917340465e-3,
                1.1383298104278345e-2,
                -1.5083449361586688e-2,
                7.9001130385225622e-3,
                -3.3220683120423703e-3,
                -1.4708290496616321e-3,
                -5.6543298116664071e-4,
                1.2076275273804822e-2,
                1.6123859181225993e-2,
                -1.1839406698705809e-2,
                -2.8379048786239491e-2,
                -2.1960543816553421e-3,
                1.5478187716138473e-3,
            ]
        ).reshape((1, 6, 3))

        def released_kernel( coordinates_angstrom, derivative_order ):
            np.testing.assert_allclose(
                coordinates_angstrom,
                coordinates_bohr * BOHR_TO_ANGSTROM,
                rtol=0.0,
                atol=0.0,
            )
            self.assertEqual(derivative_order, 1)
            return (
                np.array([6.1879370436879544e-1]),
                raw_gradient,
                np.zeros((1, 1, 6, 3)),
            )

        extension = types.SimpleNamespace(pes=released_kernel)
        surface = create_ch5_j1_surface(extension)
        result = surface.evaluate(coordinates_bohr, PESRequest(gradient=True))

        self.assertEqual(
            surface.metadata.atom_symbols,
            ("H", "H", "H", "H", "C", "H"),
        )
        np.testing.assert_allclose(
            result.energy,
            np.array([2.2740256820584732e-2]),
            rtol=0.0,
            atol=8.0e-9,
        )
        np.testing.assert_allclose(
            result.gradient,
            expected_gradient,
            rtol=0.0,
            atol=2.0e-8,
        )


if __name__ == "__main__":
    unittest.main()
