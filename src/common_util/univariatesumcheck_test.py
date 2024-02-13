import random, unittest
import numpy as np
from src.common_util.poly import (
    Polynomial, 
    Basis, 
    Scalar, 
    construct_non_multiplicativegroup_vanishing_poly
)
from src.common_util.curve import left_rotate_root_of_unity
from src.common_util.univariatesumcheck import (
    create_col_mapping,
    construct_H_I,
    construct_Dx,
    construct_rho_col_j, 
    prove as univariate_sumcheck_prove
)
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
M2 = [
    [0, 0, 1, 0],
    [0, 0, 1, 0],
    [0, 0, 1, 0],
]
# the table we look up on
c = [1, 2, 3, 4]
# the result be looked up
a = [3, 2, 3]
a2 = [3, 3, 3]
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
        H = left_rotate_root_of_unity(H)
        alpha = Scalar(3)
        Dx = construct_Dx(M, V, H, a, alpha)
        # TODO: check this statement: A proof of formation of H_i requires caulk+ core. t is 2, 3, so H_i corresponding to their index are H[1], H[2]
        
        V = list(map(Scalar, [1, 2, 3]))
        H = list(map(Scalar, [1, 2, 3, 4]))
        alpha = Scalar(1)
        Dx = construct_Dx(M2, V, H, a2, alpha)
        assert Dx == Scalar(1)

    def test_construct_rho_col_j(self) -> None:
        """ 
        M = [
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
        ],
        col_mapping = [3, 2, 3]
        """
        col_mapping = create_col_mapping(M)
        H = Scalar.roots_of_unity(len(c))
        H = left_rotate_root_of_unity(H)
        rho_col_j = construct_rho_col_j(col_mapping, H)
        assert len(rho_col_j) == 3, "3 columns in col_mapping"
        assert np.all(rho_col_j[0].values == rho_col_j[2].values), "1st and 3rd element are rho_3"
        # assert rho_col_j == [Polynomial(list(map(Scalar, [1, 2, 3, 4])), Basis.LAGRANGE), Polynomial(list(map(Scalar, [1, 2, 3, 4])), Basis.LAGRANGE), Polynomial(list(map(Scalar, [1, 2, 3, 4])), Basis.LAGRANGE)]

    def test_prove(self) -> None:
        m = len(M) # number of rows = number of mu
        N = len(c) # table size

        V = Scalar.roots_of_unity(m) # set for mu
        # V = left_rotate_root_of_unity(V)
        H = Scalar.roots_of_unity(N) # set for rho
        # H = left_rotate_root_of_unity(H)
        alpha = Scalar(3)
        
        Dx = construct_Dx(M, V, H, a, alpha)
        # t_x = Polynomial(list(map(Scalar, [2, 3])), Basis.LAGRANGE)
        
        zero_polynomial = Polynomial([Scalar(0)], Basis.MONOMIAL)
        # interpolate through H_I, shouldn't use default lagrange interpolation and fft 
        H_I = construct_H_I(H, M)
        t_x = sum([lagrange_basis(i, H_I) * Scalar(e) for i, e in enumerate([2, 3])], zero_polynomial)


        # phi_x = Polynomial(list(map(Scalar, [3, 2, 3])), Basis.LAGRANGE)
            # Note: not 2^k exponentiation, shouldn't use fft
        phi_x = sum([lagrange_basis(i, V) * Scalar(e) for i, e in enumerate([3, 2, 3])], zero_polynomial)
        zI_x = construct_non_multiplicativegroup_vanishing_poly(H_I)

        setup = Setup.execute(2, 2, [1, 2], False) # TODO: check arguments
        pi2 = univariate_sumcheck_prove(setup, Dx, t_x, phi_x, zI_x, alpha, N, m)


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