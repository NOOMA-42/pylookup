import py_ecc.bn128 as b
from dataclasses import dataclass
from src.common_util.curve import ec_lincomb, G1Point, G2Point, Scalar
from src.common_util.poly import Polynomial, Basis, PolyUtil
from src.baloo.verifier import VerificationKey

@dataclass
class Setup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    powers_of_x2: list[G2Point]
    z_H_comm_1: G2Point
    t_comm_1: G2Point
    # coeffs of quotient commitment polynomial of T(X) in coefficient form
    Q_T_comm_poly_coeffs: list[G1Point]

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

        for i in range(1, powers):
            powers_of_x[i] = b.multiply(powers_of_x[i - 1], tau)

        assert b.is_on_curve(powers_of_x[1], b.b)
        print("Generated G1 side, X^1 point: {}".format(powers_of_x[1]))

        powers_of_x2 = [0] * (powers + 1)
        powers_of_x2[0] = b.G2
        for i in range(1, powers + 1):
            powers_of_x2[i] = b.multiply(powers_of_x2[i - 1], tau)

        assert b.is_on_curve(powers_of_x2[1], b.b2)
        print("Generated G2 side, X^1 point: {}".format(powers_of_x2[1]))

        # assert b.pairing(b.G2, powers_of_x[1]) == b.pairing(powers_of_x2[1], b.G1)
        print("X^1 points checked consistent")
        print("Finished to generate structured reference string")
        self.powers_of_x = powers_of_x
        self.powers_of_x2 = powers_of_x2
        return (powers_of_x, powers_of_x2)

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

    @classmethod
    def execute(self, powers: int, tau: int, public_table: list):
        # 1. generate_srs: will do in the runtime
        self.generate_srs(powers, tau)
        table_len = len(public_table)
        # 2. Compute and output [z_H(x)] * G1
        # vanishing polynomial: X^N - 1, N = group_order_N - 1
        poly_util = PolyUtil()
        z_H_array = poly_util.vanishing_poly(table_len)
        # in coefficient form
        z_H_poly = Polynomial(z_H_array, Basis.MONOMIAL)

        z_H_comm_1 = self.commit_g1(z_H_poly)
        print("Commitment of z_H(X) on G1: ", z_H_comm_1)
        # 3. Compute and output [T(x)] * G1
        t_values = [Scalar(val) for val in public_table]
        t_poly = Polynomial(t_values, Basis.LAGRANGE)
        t_comm_1 = self.commit_g1(t_poly)
        print("Commitment of T(X) on G1: ", t_comm_1)
        self.z_H_comm_1 = z_H_comm_1
        self.t_comm_1 = t_comm_1
        print("setup complete")
        return self

    @classmethod
    def verification_key(self) -> VerificationKey:
        return VerificationKey(
            self.powers_of_x2[1],
            self.t_comm_1,
            self.z_H_comm_1,
        )
