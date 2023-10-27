import os
import unittest
import random
import py_ecc.bn128 as b
from common_util.poly import Polynomial, Basis
from common_util.curve import Scalar
from common_util.kzg import Setup
from common_util.fk20_single import fk20
cwd = os.getcwd()

class Fk20SingleTest(unittest.TestCase):
    def test_toeplitz(self) -> None:        
        # compute SRS
        MAX_DEGREE_POLY = b.curve_order - 1 # curve order is modulus
        N_POINTS = 512

        setup = Setup.from_file(cwd + "/src/cq/test/powersOfTau28_hez_final_11.ptau")
        polynomial = Polynomial(
            [random.randint(1, MAX_DEGREE_POLY) for _ in range(N_POINTS)],
            Basis.LAGRANGE,
        )
        commitment = setup.commit_G2(polynomial)
        fk20(polynomial, setup)