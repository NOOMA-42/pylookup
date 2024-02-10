from dataclasses import dataclass

from py_ecc.secp256k1.secp256k1 import bytes_to_int

from caulk_setup import Setup
from src.caulk.util import hash_ec_points
from src.common_util.curve_optimized import Scalar, ec_mul, G1, ec_add, ec_sub, ec_pairing, G2, G1Point, G2Point, \
    ec_lincomb, ec_eq
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


def single_term_poly(degree: int):
    vals = [Scalar(0)] * degree + [Scalar(1)]
    return Polynomial(vals, Basis.MONOMIAL)


# Change of variable e.g. f(X) -> f(aX)


@dataclass
class Proof_pederson:
    R: G1Point
    t1: Scalar
    t2: Scalar


@dataclass
class Proof_unity:
    z_comm: G2Point
    f_comm: G1Point
    h_comm: G1Point
    alpha: Scalar
    v1: Scalar
    v2: Scalar
    pi_f_1: G1Point
    pi_f_2: G1Point
    pi_p_alpha: G1Point


@dataclass
class Proof:
    cm: G1Point
    z_comm: G2Point
    T_comm: G1Point
    S_comm: G2Point
    proof_pederson: Proof_pederson
    proof_unity: Proof_unity


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
        proof_unity = self.prove_unity(a, i)

        return Proof(cm, z_comm, T_comm, S_comm, proof_pederson, proof_unity)

    def prove_unity(self, a: Scalar, i: int):
        setup = self.setup
        logN = setup.logN
        n = setup.n
        b = a * self.setup.roots_N[i]
        r0, r1, r2, r3 = Scalar(11), Scalar(22), Scalar(33), Scalar(44)
        r_poly = Polynomial([r1, r2, r3], Basis.MONOMIAL)
        z_vn = vanishing_poly(n)
        rho_polys = lagrange_polys(n)
        sigma = setup.roots_n[1]
        z_poly = Polynomial([-b, a])
        z_comm = setup.kzgSetup.commit_G2(z_poly)

        # setup f(X)
        f_values = [Scalar(0)] * (logN + 6)
        f_values[0] = a - b
        f_values[1] = a * sigma - b
        f_values[2] = a
        f_values[3] = b
        f_values[4] = a / b
        for i in range(5, logN + 5):
            f_values[i] = f_values[i - 1] ** 2
        f_values[logN + 5] += r0
        f_poly = Polynomial(f_values, Basis.LAGRANGE).ifft()
        f_poly += r_poly * z_vn

        f_shift_1 = f_poly.scale(setup.roots_n[-1])  # f(sigma^-1 * X)
        assert f_shift_1.eval(sigma) == f_poly.eval(Scalar(1))
        f_shift_2 = f_poly.scale(setup.roots_n[-2])  # f(sigma^-2 * X)

        # setup p(X)
        p_poly = (f_poly - z_poly) * (rho_polys[0] + rho_polys[1])
        p_poly += (f_poly * (Scalar(1) - sigma) - f_shift_2 + f_shift_1) * rho_polys[2]
        p_poly += (f_poly + f_shift_2 - f_shift_1 * sigma) * rho_polys[3]
        p_poly += (f_poly * f_shift_1 - f_shift_2) * rho_polys[4]

        p_poly += (f_poly - (f_shift_1 * f_shift_1)) * setup.poly_prod
        p_poly += (f_shift_1 - Scalar(1)) * rho_polys[n - 1]

        h_hat_poly = p_poly / z_vn
        h_poly = h_hat_poly + (single_term_poly(setup.kzgSetup.length - 2) * z_poly)

        f_comm = setup.kzgSetup.commit_G1(f_poly)
        h_comm = setup.kzgSetup.commit_G1(h_poly)

        # alpha is verifier randomness
        alpha = Scalar(999)
        alpha_1 = setup.roots_n[-1] * alpha
        alpha_2 = setup.roots_n[-2] * alpha
        f_alpha_1 = f_poly.eval(alpha_1)
        f_alpha_2 = f_poly.eval(alpha_2)
        v1, pi_f_1 = setup.kzgSetup.open(f_poly, alpha_1)
        v2, pi_f_2 = setup.kzgSetup.open(f_poly, alpha_2)

        p_alpha_poly = h_hat_poly * (-z_vn.eval(alpha))
        p_alpha_poly += (f_poly - z_poly) * (rho_polys[0].eval(alpha) + rho_polys[1].eval(alpha))
        p_alpha_poly += (f_poly * (Scalar(1) - sigma) - f_alpha_2 + f_alpha_1) * rho_polys[2].eval(alpha)
        p_alpha_poly += (f_poly + f_alpha_2 - sigma * f_alpha_1) * rho_polys[3].eval(alpha)
        p_alpha_poly += (f_poly * f_alpha_1 - f_alpha_2) * rho_polys[4].eval(alpha)
        p_alpha_poly += (f_poly - f_alpha_1 * f_alpha_1) * setup.poly_prod.eval(alpha)
        p_alpha_poly += (f_alpha_1 - Scalar(1)) * rho_polys[n - 1].eval(alpha)

        assert p_alpha_poly.eval(alpha) == Scalar(0)

        _, pi_p_alpha = setup.kzgSetup.open(p_alpha_poly, alpha)

        return Proof_unity(z_comm, f_comm, h_comm, alpha,
                           v1, v2, pi_f_1, pi_f_2, pi_p_alpha)

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
    setup = Setup.example_setup()
    prover = Prover(setup)
    proof = prover.prove(Scalar(2))
