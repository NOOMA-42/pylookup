import os
import unittest
import random
import py_ecc.bn128 as b
from common_util.poly import Polynomial, Basis
from common_util.curve import Scalar
from common_util.kzg import Setup
from common_util.fk20_single import test_toeplitz_part1#fk20
cwd = os.getcwd()

class Fk20SingleTest(unittest.TestCase):
    def test_toeplitz(self) -> None:        
        # compute SRS
        MAX_DEGREE_POLY = b.curve_order - 1 # curve order is modulus
        # N_POINTS = 512

        setup = Setup.from_file(cwd + "/src/common_util/setup_file/powersOfTau28_hez_final_08.ptau")
        """ setup.powers_of_x = setup.powers_of_x[:4]
        setup.powers_of_x2 = setup.powers_of_x2[:4]
        setup.length = 4 """
        """ polynomial = Polynomial(
            list(map(Scalar, [random.randint(1, MAX_DEGREE_POLY) for _ in range(setup.length)])),
            Basis.LAGRANGE,
        ) """
        #commitment = setup.commit_G2(polynomial)
        test_toeplitz_part1(setup)
        # fk20(polynomial, setup)