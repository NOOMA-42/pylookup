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
    prove as univariate_sumcheck_prove,
    verify as univariate_sumcheck_verify
)
from src.common_util.lagrange import lagrange_basis, multilinear_lagrange_kernel_to_uni
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
        H = Scalar.roots_of_unity(N) # set for rho
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

        # The largest exponent is R_hat with exponent 4 in this test
        setup = Setup.execute(4, 2, [1, 2], False) # FIXME: modify generate_srs and use it instead
        pi2 = univariate_sumcheck_prove(setup, Dx, t_x, phi_x, zI_x, alpha, N, m)

    @unittest.skip(reason="This test is disabled, because it's too slow.")
    def test_verify(self) -> None:
        m = len(M) # number of rows = number of mu
        N = len(c) # table size
        V = Scalar.roots_of_unity(m) # set for mu
        H = Scalar.roots_of_unity(N) # set for rho
        alpha = Scalar(3)
        Dx = construct_Dx(M, V, H, a, alpha)        
        zero_polynomial = Polynomial([Scalar(0)], Basis.MONOMIAL)
        H_I = construct_H_I(H, M)
        t_x = sum([lagrange_basis(i, H_I) * Scalar(e) for i, e in enumerate([2, 3])], zero_polynomial)
        phi_x = sum([lagrange_basis(i, V) * Scalar(e) for i, e in enumerate([3, 2, 3])], zero_polynomial)
        zI_x = construct_non_multiplicativegroup_vanishing_poly(H_I)

        setup = Setup.execute(4, 2, [1, 2], False)
        pi = univariate_sumcheck_prove(setup, Dx, t_x, phi_x, zI_x, alpha, N, m)

        univariate_sumcheck_verify(pi, (t_x, zI_x), setup)

    def test_multilinear(self) -> None:
        v = Polynomial(list(map(Scalar, [1, 2])), Basis.LAGRANGE)
        c = multilinear_lagrange_kernel_to_uni([[1, -1], [1, 1], [-1, 1], [-1, -1]], [6, 3])
        f = [v_i * c_i for v_i, c_i in zip(v.values, c.values)]
        f_x = Polynomial(f, Basis.LAGRANGE)
        alpha = Scalar(3)

        f_x = f_x.ifft()
        v = v.ifft()
        c = c.ifft()
        
        setup = Setup.execute(4, 2, [1, 2], False)
        H = Scalar.roots_of_unity(len([1, 2]))
        z_x = construct_non_multiplicativegroup_vanishing_poly(H)
        pi = univariate_sumcheck_prove(setup, v, c, f_x, z_x, alpha, 1, 1) # N M are dummy rn
        
        univariate_sumcheck_verify(pi, (c, z_x), setup)