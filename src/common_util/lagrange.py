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


def multilinear_lagrange_kernel(X, x):
    """
    Follow the formula in logup+GKR
    """
    assert len(X) == len(x), "X and x must have the same length"
    l_i = Scalar(1)
    for Xi, xi in zip(X, x):
        print("Xi, xi", Xi, xi)
        l_i = (l_i * (Scalar(1) + Scalar(Xi) * Scalar(xi)))
    l_i = l_i / pow(Scalar(2), len(X))
    return l_i

def multilinear_lagrange_kernel_to_uni(Xs, t):
    """
    Parameters:
    - Xs are a set of multilinear vectors where lagrange kernel is interpolate on
    - t is a multivariate query point and will be fed into the lagrange kernel

    Returns:
    - Polynomial of the lagrange kernel at t
    """
    c = []
    for X in Xs:
        assert len(X) == len(t), "X and t must have the same length"
        c.append(multilinear_lagrange_kernel(X, t))

    return Polynomial(list(map(Scalar, c)), Basis.LAGRANGE)