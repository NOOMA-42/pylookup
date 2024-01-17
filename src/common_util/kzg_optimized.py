import py_ecc.optimized_bn128 as b
from src.common_util.curve_optimized import ec_lincomb, G1Point, G2Point, Scalar, G1, G2, ec_pairing, ec_mul, ec_sub
from src.common_util.poly_optimized import Polynomial, Basis
from dataclasses import dataclass

# Recover the trusted setup from a file in the format used in
# https://github.com/iden3/snarkjs#7-prepare-phase-2
SETUP_FILE_G1_STARTPOS = 80
SETUP_FILE_POWERS_POS = 60


def _commit(poly: Polynomial, powers: list):
    if poly.basis == Basis.LAGRANGE:
        # inverse FFT from Lagrange basis to monomial basis
        coeffs = poly.ifft().values
    else:
        coeffs = poly.values

    if len(coeffs) > len(powers):
        raise Exception("Not enough powers in setup")
    return ec_lincomb([(s, x) for s, x in zip(powers, coeffs)])


@dataclass
class KZGSetup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( H,    xH,  ...,  x^{d-1}H ), where H is a generator of G_2
    powers_of_x2: list[G2Point]
    length: int

    @classmethod
    def manual_setup(cls, length: int = 100, tau: Scalar = Scalar(5)):
        powers_of_x = [G1]
        powers_of_x2 = [G2]
        for i in range(1, length):
            powers_of_x.append(b.multiply(powers_of_x[-1], int(tau)))
            powers_of_x2.append(b.multiply(powers_of_x2[-1], int(tau)))

        assert b.pairing(b.G2, powers_of_x[1]) == b.pairing(powers_of_x2[1], b.G1)
        return cls(powers_of_x, powers_of_x2, length)

    # Encodes the KZG commitment that evaluates to the given values in the group
    def commit_G1(self, poly: Polynomial) -> G1Point:
        return _commit(poly, self.powers_of_x)

    def commit_G2(self, poly: Polynomial) -> G2Point:
        return _commit(poly, self.powers_of_x2)

    def open(self, poly: Polynomial, x: Scalar) -> (Scalar, G1Point):
        v = poly.eval(x)
        if poly.basis == Basis.LAGRANGE:
            poly = poly.ifft()
        q = (poly - v) / Polynomial([Scalar(-x), Scalar(1)], Basis.MONOMIAL)
        return v, self.commit_G1(q)

    def verify(self, C: G1Point, pi: G1Point, x: Scalar, v: Scalar):
        LHS = ec_pairing(G2, ec_sub(C, ec_mul(G1, v)))
        comm = self.commit_G2(Polynomial([Scalar(-x), Scalar(1)], Basis.MONOMIAL))
        RHS = ec_pairing(comm, pi)
        assert LHS == RHS


if __name__ == "__main__":
    kzgSetup = KZGSetup.manual_setup()
    poly = Polynomial([Scalar(1), Scalar(2), Scalar(3), Scalar(4)], Basis.LAGRANGE)
    c = kzgSetup.commit_G1(poly)
    v, pi = kzgSetup.open(poly, Scalar(2))
    print(v, pi)
    kzgSetup.verify(c, pi, Scalar(2), v)
