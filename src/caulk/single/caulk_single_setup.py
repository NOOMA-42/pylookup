from dataclasses import dataclass
import math

from src.common_util.curve_optimized import Scalar, G1Point, G1, ec_mul
from src.common_util.kzg_optimized import KZGSetup
from src.common_util.poly_optimized import Polynomial, Basis


@dataclass
class CaulkSingleSetup:
    kzgSetup: KZGSetup
    c: list[Scalar]
    c_poly: Polynomial
    c_commit: G1Point
    h: Scalar
    h_1: G1Point
    N: int
    logN: int
    roots_N: list[Scalar]
    n: int
    roots_n: list[Scalar]
    poly_prod: Polynomial

    @classmethod
    def example_setup(cls):
        kzgSetup = KZGSetup.manual_setup()
        c = list(map(Scalar, [1, 2, 3, 4]))
        c_poly = Polynomial(c, Basis.LAGRANGE)
        c_commit = kzgSetup.commit_G1(c_poly)
        h = Scalar(123)
        h_1 = ec_mul(G1, h)
        N = len(c)
        roots_N = Scalar.roots_of_unity(N)
        logN = math.ceil(math.log2(N))
        # n = logN + 6
        n = logN + 6
        roots_n = Scalar.roots_of_unity(n)

        # poly_prod = (X - 1) (X - w) (X - w^2) (X - w^3) (X - w^4) (X - w^(5 + logN)) (X - w^(6 + logN))
        poly_prod = Polynomial([Scalar(1)])
        for i in range(n):
            if i < 5 or i >= 5 + logN:
                poly_prod = poly_prod * Polynomial([Scalar(-roots_n[i]), Scalar(1)])

        return cls(kzgSetup, c, c_poly, c_commit, h, h_1
                   , N, logN, roots_N, n, roots_n, poly_prod)


if __name__ == "__main__":
    setup = CaulkSingleSetup.example_setup()
    print(setup.roots_N)
