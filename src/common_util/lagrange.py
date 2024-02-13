from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar

def lagrange_basis(i: int, xi: list) -> Polynomial:
    """
    Calculate the i-th Lagrange polynomial at point X.
    
    Parameters:
    - i: Index of the Lagrange polynomial to calculate.
    - xi: Set of interpolation points.
    - X: Point at which to evaluate the polynomial.
    
    Returns:
    The value of the i-th Lagrange polynomial at X.
    
    Note:
    τcol(j) in baloo paper is not compatible with fft from lagrange basis to monomial.
    we really have to calculate the lagrange basis and combine them into polynomial and commit

    TODO: review this
    """
    l_i = Scalar(1)
    for j, xj in enumerate(xi):
        if j != i:
            X = Polynomial(list(map(Scalar, [0, 1])), Basis.MONOMIAL)
            l_i = (X - Scalar(xj)) / (Scalar(xi[i]) - Scalar(xj)) * l_i
    return l_i



