from dataclasses import dataclass

from py_ecc.secp256k1.secp256k1 import bytes_to_int

from setup import Setup
from src.caulk.util import hash_ec_points
from src.common_util.curve_optimized import Scalar, ec_mul, G1, ec_add, ec_sub, ec_pairing, G2, G1Point, G2Point
from src.common_util.merlin.keccak import SHA3_256
from src.common_util.poly_optimized import Polynomial, Basis


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

        return Proof(cm, z_comm, T_comm, S_comm, proof_pederson)

    def prove_unity(self, a: Scalar, i: int):
        setup = self.setup
        b = a * self.setup.roots_N[i]

        f_values = []
        f_values.append(a - b)
        f_values.append(a * setup.roots_n[1] - b)
        f_values.append(a)
        f_values.append(b)

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
    setup = Setup.example_setup()
    prover = Prover(setup)
    proof = prover.prove(Scalar(2))
