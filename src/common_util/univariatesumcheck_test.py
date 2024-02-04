import random, unittest
from src.common_util.poly import Polynomial, Basis, Scalar
from src.common_util.univariatesumcheck import create_col_mapping, construct_Dx
from src.common_util.lagrange import lagrange_basis
from src.cq.setup import Setup # TODO: refactor to common_util


""" 
(i) e (t, [D]2) − e ([1]1u2, [1]2) = e ([R]1, [x]2) + e ([Q2]1, [zI ]2) 
"""

# lookup matrix
M = [
    [0, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
]
# the table we look up on
c = [1, 2, 3, 4]
# the result be looked up
a = [3, 2, 3]
# unique of a
t = [2, 3]

class UnivariateSumcheckTest(unittest.TestCase):
    def test_matrix_encode(self) -> None:
        col_mapping = create_col_mapping(M)
        assert list(col_mapping) == [3, 2, 3]
    
    def test_construct_Dx(self) -> None:        
        m = len(M) # number of rows = number of mu
        N = len(c) # table size

        V = Scalar.roots_of_unity(m) # set for mu
        H = Scalar.roots_of_unity(N) # set for rho
        alpha = Scalar(3)
        construct_Dx(M, V, H, a, alpha)
        # TODO: check this statement: A good formation of H_i requires caulk+ core. t is 2, 3, so H_i corresponding to their index are H[1], H[2]
        

    """ def test_lagrange_basis(self) -> None:
        # public table
        # table = [1...256]
        table = []
        for i in range(1, 257):
            table.append(i)
        print("table: ", table)

        group_order_N = len(table)
        # number of powers of tau
        powers = group_order_N * 2
        # do setup
        setup = Setup.execute(powers, tau, table)

        # Encode the matrix M to 
        M = [
            [0, 1, 0],
            [1, 0, 0],
            [0, 0, 1],
        ]
        col_mapping = create_col_mapping(M)
        print(col_mapping)


        # phi_x is encoding of "a" vector, aka the values be looked up
        phi_x = Polynomial(list(map(Scalar, [3, 2, 3])), Basis.LAGRANGE)
        # t_x is encoding of "t" vector, aka the unique value in a 
        t_x = Polynomial(list(map(Scalar, [2, 3])), Basis.LAGRANGE)

        z_I_x """