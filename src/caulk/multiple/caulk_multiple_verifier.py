from src.caulk.multiple.caulk_multiple_prover import Proof
from src.caulk.multiple.caulk_multiple_setup import CaulkMultipleSetup
from src.caulk.multiple.caulk_multiple_transcript import Transcript
from src.caulk.util import vanishing_poly
from src.common_util.curve_optimized import ec_add, ec_mul, ec_eq, Scalar, ec_lincomb, G1, ec_pairing, G2, ec_sub


class CaulkMultipleVerifier:

    def __init__(self, setup: CaulkMultipleSetup):
        self.setup = setup

    def verify(self, proof: Proof):
        transcript = Transcript(b"caulk-multiple")

        chi = transcript.round1(proof.message1)
        alpha = transcript.round2(proof.message2)

        z_v_m_poly = vanishing_poly(proof.m)
        z_v_m_alpha = z_v_m_poly.eval(alpha)

        g1_P1 = ec_add(proof.message1.g1_Z_I, ec_mul(proof.message1.g1_C_I, chi))
        assert ec_eq(g1_P1, proof.test_g1_P1)

        g1_P2 = ec_lincomb([
            (G1, proof.message3.v2),
            (proof.cm, -chi),
            (proof.message2.g1_H2, -z_v_m_alpha)
        ])
        assert ec_eq(g1_P2, proof.test_g1_P2)

        self.setup.kzgSetup.verify(proof.message1.g1_u, proof.message3.pi_1, alpha, proof.message3.v1)
        self.setup.kzgSetup.verify(g1_P1, proof.message3.pi_2, proof.message3.v1, proof.message3.v2)
        self.setup.kzgSetup.verify(g1_P2, proof.message3.pi_3, alpha, Scalar(0))

        LHS = ec_pairing(G2, ec_sub(self.setup.c_commit, proof.message1.g1_C_I))
        RHS = ec_pairing(proof.message1.g2_H1, proof.message1.g1_Z_I)
        assert LHS == RHS

        print("Verification successful!")
