import py_ecc.bn128 as b
from dataclasses import dataclass
from src.common_util.curve import ec_lincomb, G1Point, G2Point, Scalar
from src.common_util.poly import Polynomial, Basis
from src.cq.verifier import VerificationKey
from src.cq.program import CommonPreprocessedInput

@dataclass
class Setup(object):
    #   ([1]₁, [x]₁, ..., [x^{d-1}]₁)
    # = ( G,    xG,  ...,  x^{d-1}G ), where G is a generator of G_1
    powers_of_x: list[G1Point]
    powers_of_x2: list[G2Point]
    Z_V_comm_2: G2Point
    T_comm_2: G2Point

    @classmethod
    # tau: a random number whatever you choose
    def generate_srs(self, powers: int, tau: int, verbose):
        if verbose:
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
        if verbose:
            print("Generated G1 side, X^1 point: {}".format(powers_of_x[1]))

        powers_of_x2 = [0] * (powers + 1)
        powers_of_x2[0] = b.G2
        # TODO check paper if this is correct
        for i in range(powers + 1):
            if i > 0:
                powers_of_x2[i] = b.multiply(powers_of_x2[i - 1], tau)

        assert b.is_on_curve(powers_of_x2[1], b.b2)
        if verbose:
            print("Generated G2 side, X^1 point: {}".format(powers_of_x2[1]))

        # assert b.pairing(b.G2, powers_of_x[1]) == b.pairing(powers_of_x2[1], b.G1)
        if verbose:
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

    @classmethod
    def execute(self, powers: int, tau: int, public_table: list, verbose: bool = True):
        """  
        you only need power and tau if you're not working on cq
        """
        # 1. generate_srs: will do in the runtime
        self.generate_srs(powers, tau, verbose)
        # 2. Compute and output [ZV(x)] * G2
        # vanishing polynomial: X^N - 1, N = group_order_N - 1
        Z_V_array = [Scalar(-1)] + [Scalar(0)] * (len(public_table) - 1) + [Scalar(1)]
        # in coefficient form
        Z_V_poly = Polynomial(Z_V_array, Basis.MONOMIAL)

        Z_V_comm_2 = self.commit_g2(Z_V_poly)
        if verbose:
            print("Commitment of Z_V(X) on G2: ", Z_V_comm_2)
        # 3. Compute and output [T(x)] * G2
        # TODO: optimization
        t_values = [Scalar(val) for val in public_table]
        T_poly = Polynomial(t_values, Basis.LAGRANGE)
        T_comm_2 = self.commit_g2(T_poly)
        if verbose:
            print("Commitment of T(X) on G2: ", T_comm_2)
        # 4. (a): qi = [Qi(x)] * G1
        # 4. (b): [Li(x)] * G1
        # 4. (c): [Li(x)−Li(0) / x] * G1
        if verbose:
            print("setup complete")
        self.Z_V_comm_2 = Z_V_comm_2
        self.T_comm_2 = T_comm_2

        return self

    @classmethod
    def verification_key(self, pk: CommonPreprocessedInput) -> VerificationKey:
        return VerificationKey(
            pk.group_order_N,
            pk.group_order_n,
            Scalar.root_of_unity(pk.group_order_N),
            Scalar.root_of_unity(pk.group_order_n),
            self.powers_of_x2,
            self.T_comm_2,
            self.Z_V_comm_2
        )
