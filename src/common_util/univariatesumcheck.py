import py_ecc.bn128 as b
from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar, ec_lincomb
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

def construct_H_I(H: list[Scalar], M: list[list]) -> list[Scalar]:
    """
    Input:
        H, e.g. [omega^2, omega^1, omega^2]
        M, e.g. [
            [0, 1],
            [1, 0],
            [0, 1],
        ]
    Returns:
        H_I = [omega^1, omega^2]
    """
    col_mapping = create_col_mapping(M)
    H_I = [H[i - 1] for i in set(col_mapping)]
    return H_I

def construct_rho_col_j(col_mapping: list[int], H) -> list[Polynomial]:
    """
    Input:
        given
        M = [
            [0, 1],
            [1, 0],
            [0, 1],
        ]

        col_mapping = [2, 1, 2]
        H = [omega^2, omega^1, omega^2]
    Returns:
        rho_col_j = [rho_2, rho_1, rho_2]
    """
    index_and_Hi = [(i, H[i - 1]) for i in set(col_mapping)]
    H_I = [e for _, e in index_and_Hi]
    rho_col_j_unique = {e[0]: lagrange_basis(i, H_I) for i, e in enumerate(index_and_Hi)}
    rho_col_j = [rho_col_j_unique[i] for i in col_mapping if i in rho_col_j_unique]
    return rho_col_j

def construct_Dx(M: list[list], V: list[Scalar], H: list[Scalar], vec_a: list, alpha: Scalar) -> Polynomial:
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
    vec a       [a,  b,  a]
                 |   |   |
                 v   v   v
    rho_j       [L1, L2] X
                    | |
                     v
    rho_j mul   [L1, L2, L1]
    """
    assert len(V) == len(M), "V, corresponding with mu, its size must be equal to the number of rows in M"
    assert len(H) == len(M[0]), "H, corresponding with rho, its size must be equal to the number of columns in M"
    assert len(vec_a) == len(M), "size of a, the looked up value, must be equal to the number of rows in M"
    col_mapping = create_col_mapping(M)

    rho_col_j = construct_rho_col_j(col_mapping, H)
    rho_col_j_to_multiply_0 = [i if isinstance(i, Scalar) else i.coeff_eval(0) for i in rho_col_j]

    mu_j = [lagrange_basis(i, V) for i in range(len(V))]
    mu_j = [i if isinstance(i, Scalar) else i.coeff_eval(alpha) for i in mu_j]

    D_x = [rho_col_j[i] * mu_j[i] / rho_col_j_to_multiply_0[i] for i in range(len(mu_j))]

    if (isinstance(rho_col_j_to_multiply_0[0], Scalar) and
        isinstance(mu_j[0], Scalar) and
        isinstance(rho_col_j[0], Scalar)
    ):
        zero_polynomial = Scalar(0)
    else:
        zero_polynomial = Polynomial([Scalar(0)], Basis.MONOMIAL)
    D_x = sum(D_x, zero_polynomial)

    return D_x

def prove(setup: Setup, D_x: Polynomial, t_x: Polynomial, phi_x: Polynomial, zI_x: Polynomial, alpha: Scalar, N: int, m: int):
    """
    Parameters:
    - setup: Setup
    - D_x: M(X, Y) evaluated partially at X = α
    - t_x: t(X), unique value in a
    - phi_x: φ(X), encoding of "a" vector, aka the values be looked up
    - zI_x: zI(X), vanishing polynomial of H_I
    - alpha: α
    - N: table size = # of column of M

    """

    # long division to get XR(X) and Q2(X)

    #phi_alpha = phi_x.barycentric_eval(alpha)
    phi_alpha = phi_x.coeff_eval(alpha)
    #d_t_x = D_x * t_x
    # phi_alpha2 = d_t_x.coeff_eval(Scalar(0))
    #print("d_t_x.values: ", d_t_x.values)
    #print("phi_alpha2: ", phi_alpha2)
    Q2_x, XR_x = (D_x * t_x - phi_alpha).div_with_remainder(zI_x)
    # Q2_x, XR_x = ((D_x * t_x - phi_alpha2)).div_with_remainder(zI_x)
    # Q2_x, XR_x = d_t_x.div_with_remainder(zI_x)
    x = Polynomial(list(map(Scalar, [0, 1])), Basis.MONOMIAL)
    print("Q2_x.values: ", Q2_x.values)
    print("XR_x.values: ", XR_x.values)
    R_x = XR_x / x
    # R_x = XR_x

    # R_hat_x = x ** (N - m + 2)
    D = setup.commit_g2(D_x)
    R = setup.commit_g1(R_x)
    # R_hat = setup.commit_g1(R_hat_x)
    Q2 = setup.commit_g1(Q2_x)
    pi = (D, R, Q2, phi_alpha) # since we only do unisumcheck, we output pi2 and phi_alpha together
    return pi

def verify(pi, common_input, setup: Setup):
    D, R, Q2, phi_alpha = pi
    t_X, Z_I = common_input

    t = setup.commit_g1(t_X)
    scalar_one_g1 = setup.commit_g1(Polynomial([Scalar(1)], Basis.MONOMIAL))
    scalar_one_g2 = setup.commit_g2(Polynomial([Scalar(1)], Basis.MONOMIAL))
    x = setup.commit_g2(Polynomial(list(map(Scalar, [0, 1])), Basis.MONOMIAL))
    Z_I = setup.commit_g2(Z_I)

    """
    original pairing check is as follows
    (b.pairing(t, D) - b.pairing(ec_lincomb([(scalar_one_g1, phi_alpha)]), scalar_one_g2)) \
        == (b.pairing(R, x) + b.pairing(Q2, Z_I))

    1. the pairing input has to be G2, G1 in that order
    2. + has to be * and - has to be moved to the other side
    """

    assert b.pairing(D, t) == \
        (b.pairing(b.G2, R) * b.pairing(Z_I, Q2) *
         b.pairing(scalar_one_g2, ec_lincomb([(scalar_one_g1, phi_alpha)])))
