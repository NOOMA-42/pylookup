import py_ecc.bn128 as b
from dataclasses import dataclass
from src.common_util.curve import ec_lincomb, G1Point, G2Point
from src.common_util.poly import Polynomial, Basis

@dataclass
class Setup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    powers_of_x2: list[G2Point]

    def __init__(self, powers: int, tau: int):
        self.generate_srs(powers, tau)

    @classmethod
    # tau: a random number whatever you choose
    def generate_srs(self, powers: int, tau: int):
        print("Start to generate structured reference string")

        # Initialize powers_of_x with 0 values
        powers_of_x = [0] * powers
        # powers_of_x[0] =  b.G1 * tau**0 = b.G1
        # powers_of_x[1] =  b.G1 * tau**1 = powers_of_x[0] * tau
        # powers_of_x[2] =  b.G1 * tau**2 = powers_of_x[1] * tau
        # ...
        # powers_of_x[i] =  b.G1 * tau**i = powers_of_x[i - 1] * tau
        powers_of_x[0] = b.G1

        for i in range(powers):
            if i > 0:
                powers_of_x[i] = b.multiply(powers_of_x[i - 1], tau)

        assert b.is_on_curve(powers_of_x[1], b.b)
        print("Generated G1 side, X^1 point: {}".format(powers_of_x[1]))

        powers_of_x2 = [0] * (powers + 1)
        powers_of_x2[0] = b.G2
        # TODO check paper if this is correct
        for i in range(powers + 1):
            if i > 0:
                powers_of_x2[i] = b.multiply(powers_of_x2[i - 1], tau)

        assert b.is_on_curve(powers_of_x2[1], b.b2)
        print("Generated G2 side, X^1 point: {}".format(powers_of_x2[1]))

        # assert b.pairing(b.G2, powers_of_x[1]) == b.pairing(powers_of_x2[1], b.G1)
        print("X^1 points checked consistent")
        print("Finished to generate structured reference string")
        self.powers_of_x = powers_of_x
        self.powers_of_x2 = powers_of_x2
        return True

    # Encodes the KZG commitment that evaluates to the given values in the group on G1
    @classmethod
    def commit_g1(self, values: Polynomial) -> G1Point:
        if (values.basis == Basis.LAGRANGE):
            # inverse FFT from Lagrange basis to monomial basis
            coeffs = values.ifft().values
        elif (values.basis == Basis.MONOMIAL):
            coeffs = values.values
        if len(coeffs) > len(self.powers_of_x):
            raise Exception("Not enough powers in setup")
        return ec_lincomb([(s, x) for s, x in zip(self.powers_of_x, coeffs)])

    # Encodes the KZG commitment that evaluates to the given values in the group on G2
    @classmethod
    def commit_g2(self, values: Polynomial) -> G2Point:
        if (values.basis == Basis.LAGRANGE):
            # inverse FFT from Lagrange basis to monomial basis
            coeffs = values.ifft().values
        elif (values.basis == Basis.MONOMIAL):
            coeffs = values.values
        if len(coeffs) > len(self.powers_of_x2):
            raise Exception("Not enough powers in setup")
        return ec_lincomb([(s, x) for s, x in zip(self.powers_of_x2, coeffs)])
