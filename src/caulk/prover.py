from dataclasses import dataclass

from py_ecc.secp256k1.secp256k1 import bytes_to_int

from setup import Setup
from src.caulk.util import hash_ec_points
from src.common_util.curve_optimized import Scalar, ec_mul, G1, ec_add, ec_sub, ec_pairing, G2, G1Point, G2Point
from src.common_util.merlin.keccak import SHA3_256
from src.common_util.poly_optimized import Polynomial, Basis


# utility functions
def vanishing_poly(n: int) -> Polynomial:
    vals = [Scalar(-1)] + [Scalar(0)] * (n - 1) + [Scalar(1)]
    return Polynomial(vals, Basis.MONOMIAL)


def lagrange_polys(n: int):
    polys = []
    for i in range(n):
        poly = Polynomial([Scalar(0)] * i + [Scalar(1)] + [Scalar(0)] * (n - i - 1), Basis.LAGRANGE)
        polys.append(poly.ifft())

    return polys


# Change of variable e.g. f(X) -> f(aX)


@dataclass
class Proof_pederson:
    R: G1Point
    t1: Scalar
    t2: Scalar


@dataclass
class Proof:
    cm: G1Point
    z_comm: G2Point
    T_comm: G1Point
    S_comm: G2Point
    proof_pederson: Proof_pederson


class Prover:
    def __init__(self, setup: Setup):
        self.setup = setup

    def prove(self, v: Scalar) -> Proof:
        setup = self.setup
        i = setup.c.index(v)
        assert setup.c_poly.ifft().eval(setup.roots_N[i]) == v

        # sample prover randomness
        a, s, r = Scalar(666), Scalar(777), Scalar(888)

        # pederson commitment
        cm = ec_mul(G1, v + setup.h * r)

        # z(x) = a(x - omega^(i-1))
        z_mono = Polynomial([a * -setup.roots_N[i], a], Basis.MONOMIAL)

        # T(x) = (c_poly - v) / z(x) + sh (using monomial division)
        T_mono = (setup.c_poly.ifft() - v) / z_mono + s * setup.h

        # S(x) = -r - sz(x)
        S_mono = z_mono * -s - r

        z_comm = setup.kzgSetup.commit_G2(z_mono)
        T_comm = setup.kzgSetup.commit_G1(T_mono)
        S_comm = setup.kzgSetup.commit_G2(S_mono)

        proof_pederson = self.prove_pederson(r, cm, v)
        self.prove_unity(a, i)

        return Proof(cm, z_comm, T_comm, S_comm, proof_pederson)

    def prove_unity(self, a: Scalar, i: int):
        setup = self.setup
        n = setup.n
        b = a * self.setup.roots_N[i]
        r0, r1, r2, r3 = Scalar(11), Scalar(22), Scalar(33), Scalar(44)
        r_poly = Polynomial([r1, r2, r3], Basis.MONOMIAL)
        z_vn = vanishing_poly(n)
        rho_polys = lagrange_polys(n)
        sigma = setup.roots_n[1]

        # setup f(X)
        f_values = [Scalar(0)] * (setup.N.bit_length() + 6)
        f_values[0] = a - b
        f_values[1] = a * sigma - b
        f_values[2] = a
        f_values[3] = b
        f_values[4] = a / b
        for i in range(5, len(f_values)):
            f_values[i] = f_values[i - 1] ** 2
        f_values[setup.N.bit_length() + 4] += r0
        f_poly = Polynomial(f_values, Basis.LAGRANGE).ifft()
        f_poly = f_poly + r_poly * z_vn

        f_shift_1 = f_poly.scale(setup.roots_n[-1])  # f(sigma^-1 * X)
        f_shift_2 = f_poly.scale(setup.roots_n[-2])  # f(sigma^-2 * X)

        # setup p(X)
        p_poly = (f_poly - (Polynomial([-b, a]))) * (rho_polys[0] + rho_polys[1])
        p_poly += ((1 - sigma) * f_poly - f_shift_2 + f_shift_1) * rho_polys[2]
        p_poly += (f_poly + f_shift_2 - sigma * f_shift_1) * rho_polys[3]
        p_poly += (f_poly * f_shift_1 - f_shift_2) * rho_polys[4]

        # poly_prod = (X - 1) (X - w) (X - w^2) (X - w^3) (X - w^4) (X - w^(5 + logN)) (X - w^(6 + logN))
        poly_prod = Polynomial([Scalar(1)])
        for i in range(n):
            if i < 5 or i >= 5 + setup.N.bit_length():
                poly_prod = poly_prod * Polynomial([Scalar(-setup.roots_n[i]), Scalar(1)])




    def prove_pederson(self, r: Scalar, cm: G1, v: Scalar):
        # pederson commitment proof. (s1, s2 are verifier randomness)
        s1 = Scalar(11)
        s2 = Scalar(22)

        R = ec_mul(G1, (s1 + self.setup.h * s2))
        # c = H(cm, R)
        c = hash_ec_points(R, cm)

        t1 = s1 + v * c
        t2 = s2 + r * c

        return Proof_pederson(R, t1, t2)


if __name__ == "__main__":
    n = 8
    print(vanishing_poly(n))
    print(lagrange_polys(n))
    # setup = Setup.example_setup()
    # prover = Prover(setup)
    # proof = prover.prove(Scalar(2))
