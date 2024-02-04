from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar
from src.cq.setup import Setup
from src.common_util.lagrange import lagrange_basis
""" 
This implements generalized univariate sumcheck from Baloo, originally from RZ21

refer to P.16
 Compute D(X ) = M (X, α) = Σmj=1 µj (α)ˆ
τcol(j)(X ) and find R(X ), Q2(X ) such that
D(X )t(X ) − φ(α) = XR(X ) + zI (X )Q2(X )
– Set ˆR = XN−m+2 – Output π2 = [D]2 = [D(x)]2, [R]1 = [R(x)]1, [ ˆR]1 = [ ˆR(x)]1, [Q2]1 = [Q2(x)]1.
"""

#mu_j
#rho_hat_col_j
#D(X )t(X ) − φ(α) = XR(X ) + zI (X )Q2(X )

import numpy as np

# TODO: refactor to baloo
def create_col_mapping(M):
    """ 
    Input:
        M = [
            [0, 1],
            [1, 0],
            [0, 1],
        ]
    Returns:
        col: [m] → [k] is such that col(j) = i if and only if Mj,i is the only non-zero element in row j of M (P.8 in Baloo)
        [2, 1, 2]
    """
    M = np.array(M)
    col_indices = np.argmax(M != 0, axis=1) + 1  # Add 1 for 1-indexing
    return col_indices

# TODO: refactor to baloo
def construct_Dx(M: list[list], V: list[Scalar], H: list[Scalar], vec_a: list, alpha: Scalar): # -> Polynomial
    """
    Constructs D_x based on matrices M, V, H, and vector vec_a.
    
    Parameters:
    - M: Lookup Matrix
    - V: Vector or list of 2nd set of roots of unity
    - H: Vector or list of 1st set of roots of unity
    - vec_a: Vector or list of values to be looked up
    
    Returns:
    - D_x: M(X, Y) evaluated partially at X = alpha 

    Note:
    [a,   b,  a]
     |    |   |
     v    v   v
    [L1,L2]   X
        | |
         v
    [L1, L2, L1]
    """
    assert len(V) == len(M), "V (mu) size must be equal to the number of rows in M"
    assert len(H) == len(M[0]), "H (rho) size must be equal to the number of columns in M"
    assert len(vec_a) == len(M), "size of a, the looked up value, must be equal to the number of rows in M"
    col_mapping = create_col_mapping(M)
    c_index = [(i, H[i]) for i in set(col_mapping)]
    H_I = [e for _, e in c_index]
    rho_col_j = {e[0]: lagrange_basis(i, H_I) for i, e in enumerate(c_index)}
    rho_col_j_to_multiply = [rho_col_j[i] for i in col_mapping if i in rho_col_j]
    mu_j = [lagrange_basis(i, V) for i in range(len(V))]
    
    mu_j = [i if isinstance(i, Scalar) else i.coeff_eval(alpha) for i in mu_j]
    rho_col_j_to_multiply_0 = [i if isinstance(i, Scalar) else i.coeff_eval(0) for i in rho_col_j_to_multiply]
    D_x = [rho_col_j_to_multiply[i] * mu_j[i] / rho_col_j_to_multiply_0[i] for i in range(len(mu_j))]
     
    if (isinstance(rho_col_j_to_multiply_0[0], Scalar) and 
        isinstance(mu_j[0], Scalar) and
        isinstance(rho_col_j_to_multiply[0], Scalar)
    ):
        zero_polynomial = Scalar(0)
    else:
        zero_polynomial = Polynomial([Scalar(0)], Basis.MONOMIAL)
    D_x = sum(D_x, zero_polynomial)
    
    return D_x

def prove(setup: Setup, D_x: Polynomial, t_x: Polynomial, phi_x: Polynomial, zI_x: Polynomial, alpha: Scalar):
    # long division to get R(X) and Q2(X)
    phi_alpha = phi_x(alpha)
    Q2_X, XR_X = (D_x * t_x - phi_alpha).polydiv(zI_x)
    R_X = XR_X / Polynomial(list(map(Scalar, [0, 1])), Basis.MONOMIAL)

    R_hat_x = Polynomial(Scalar(x) ** (N - m + 2))
    D = setup.commit_g2(D_x)
    R = setup.commit_g1(R_x)
    R_hat = setup.commit_g1(R_hat_x)
    Q2 = setup.commit_g1(Q2_x)

""" def verify():
    D, tX, phi_alpha, zI = """

