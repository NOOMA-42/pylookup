from dataclasses import dataclass

from src.common_util.curve_optimized import Scalar, G1Point, G1, ec_mul
from src.common_util.kzg_optimized import KZGSetup
from src.common_util.poly_optimized import Polynomial, Basis


@dataclass
class Setup:
    kzgSetup: KZGSetup
    c: list[Scalar]
    c_poly: Polynomial
    c_commit: G1Point
    h: Scalar
    h_1: G1Point
    N: int
    roots_N: list[Scalar]
    n: int
    roots_n: list[Scalar]

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
        # n = logN + 6
        n = N.bit_length() + 6
        roots_n = Scalar.roots_of_unity(n)

        return cls(kzgSetup, c, c_poly, c_commit, h, h_1
                   , N, roots_N, n, roots_n)


if __name__ == "__main__":
    setup = Setup.example_setup()
    print(setup.roots_N)
