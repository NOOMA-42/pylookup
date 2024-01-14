from dataclasses import dataclass
from setup import Setup
from src.caulk.prover import Proof, Prover, Proof_pederson
from src.caulk.util import hash_ec_points
from src.common_util.curve_optimized import Scalar, ec_mul, G1, ec_add, ec_sub, ec_pairing, G2, G1Point, ec_eq
from src.common_util.poly_optimized import Polynomial, Basis


class Verifier:
    setup: Setup

    def __init__(self, setup: Setup):
        self.setup = setup

    def verify(self, proof: Proof):
        # LHS = e(c-cm, [1]_2)
        LHS = ec_pairing(G2, ec_sub(setup.c_commit, proof.cm))

        b = ec_pairing(proof.z_comm, proof.T_comm)
        c = ec_pairing(proof.S_comm, setup.h_1)

        assert LHS == b * c
        print("pairing check passed")

        self.verify_pederson(proof.cm, proof.proof_pederson)
        print("pederson check passed")

    def verify_pederson(self, cm: G1Point, proof: Proof_pederson):
        c = hash_ec_points(proof.R, cm)
        LHS = ec_add(proof.R, ec_mul(cm, c))
        RHS = ec_mul(G1, proof.t1 + self.setup.h * proof.t2)
        assert ec_eq(LHS, RHS)


if __name__ == "__main__":
    setup = Setup.example_setup()
    prover = Prover(setup)
    proof = prover.prove(Scalar(2))

    verifier = Verifier(setup)
    verifier.verify(proof)
