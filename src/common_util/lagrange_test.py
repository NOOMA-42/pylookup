import random, unittest
from src.common_util.lagrange import lagrange_basis
from src.common_util.curve import Scalar

class LagrangeBasisTest(unittest.TestCase):
    def test_lagrange_basis(self) -> None:
        """
        evaluate following larange polynomial at x = 1, should be 0
        (x-1)(x-3)
        ----------
        (2-1)(2-3)
        """
        l_i = lagrange_basis(1, list(map(Scalar, [1, 2, 3])))
        assert l_i.coeff_eval(Scalar(1)) == 0
        assert l_i.coeff_eval(Scalar(3)) == 0
