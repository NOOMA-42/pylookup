from common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar
from src.cq.setup import Setup

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
    M = np.array(M)
    col_indices = np.argmax(matrix_M != 0, axis=1) + 1  # Add 1 for 1-indexing
    col_mapping = {j: col_indices[j-1] for j in range(1, matrix_M.shape[0]+1)}
    return col_mapping

def matrix_encode():
    pass

def prove(Setup: setup, Polynomial: D_x, Polynomial: t_x, Polynomial: phi_x, Polynomial: zI_x, Scalar: alpha):
    # long division to get R(X) and Q2(X)
    phi_alpha = phi_x(alpha)
    Q2_X, XR_X = (D_x * t_x - phi_alpha).polydiv(zI_x)
    R_X = XR_X / Polynomial(list(map(Scalar, [0, 1])), Basis.MONOMIAL)

    R_hat_x = Polynomial(Scalar(x) ** (N - m + 2))
    D = setup.commit_g2(D_x)
    R = setup.commit_g1(R_x)
    R_hat = setup.commit_g1(R_hat_x)
    Q2 = setup.commit_g1(Q2_x)

def verify():
    D, tX, phi_alpha, zI =

