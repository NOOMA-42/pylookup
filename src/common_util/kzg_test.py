import os
import unittest
import random
import py_ecc.bn128 as b
from src.common_util.curve import ec_lincomb, G1Point, G2Point, Scalar
from src.common_util.poly import Polynomial, Basis
from src.common_util.kzg import Setup
#from common_util.fk20_single import fk20
cwd = os.getcwd()

class KzgTest(unittest.TestCase):
    def disable_test_setup(self) -> None:      
        setup = Setup.from_file(cwd + "/src/common_util/setup_file/powersOfTau28_hez_final_11.ptau")
        dummy_values = Polynomial(
            list(map(Scalar, [1, 2, 3, 4, 5, 6, 7, 8])), Basis.LAGRANGE
        )
        commitment = setup.commit_G1(dummy_values)
        assert commitment == G1Point(
            (
                16120260411117808045030798560855586501988622612038310041007562782458075125622,
                3125847109934958347271782137825877642397632921923926105820408033549219695465,
            )
        )
        #print("Pass Setup Test")