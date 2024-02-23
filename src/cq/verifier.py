import py_ecc.bn128 as b
from dataclasses import dataclass
from src.cq.transcript import Transcript
from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import ec_lincomb, G2Point, Scalar

@dataclass
class VerificationKey:
    group_order_N: int
    group_order_n: int
    N_w: Scalar
    n_w: Scalar
    powers_of_x2: list[G2Point]
    T_comm_2: G2Point
    Z_V_comm_2: G2Point

    def verify_proof(self, pf, setup) -> bool:
        print("Start to verify proof")
        beta, gamma, eta = self.compute_challenges(pf)
        group_order_N = self.group_order_N
        group_order_n = self.group_order_n
        powers_of_x2 = self.powers_of_x2
        proof = pf.flatten()
        # get commitments
        m_comm_1 = proof["m_comm_1"]
        A_comm_1 = proof["A_comm_1"]
        Q_A_comm_1 = proof["Q_A_comm_1"]
        f_comm_1 = proof["f_comm_1"]
        B_0_comm_1 = proof["B_0_comm_1"]
        Q_B_comm_1 = proof["Q_B_comm_1"]
        P_comm_1 = proof["P_comm_1"]
        b_0_at_gamma = proof["b_0_at_gamma"]
        f_at_gamma = proof["f_at_gamma"]
        a_at_0 = proof["a_at_0"]
        pi_gamma = proof["pi_gamma"]
        a_0_comm_1 = proof["a_0_comm_1"]

        ### Check 1: round 2.11: A encodes the correct values ###
        print("=== Started Check 1: round 2.11: A encodes the correct values ===")
        comb = ec_lincomb([
            (m_comm_1, 1),
            (A_comm_1, -beta)
        ])
        A_check_lhs1 = b.pairing(self.T_comm_2, A_comm_1)
        A_check_rhs1 = b.pairing(self.Z_V_comm_2, Q_A_comm_1)
        A_check_rhs2 = b.pairing(b.G2, comb)
        assert A_check_lhs1 == A_check_rhs1 * A_check_rhs2, "Check 1 failed: A encodes the correct values"
        print("=== Finished Check 1: round 2.11: A encodes the correct values ===")

        ### Check 2: round 2.12: B_0 has the appropriate degree ###
        print("=== Started Check 2: B_0 has the appropriate degree ===")
        # TODO Put it into common preprocessed input?
        x_exponent_order = group_order_N - 1 - (group_order_n - 2)
        x_exponent_values_in_coeff = [Scalar(0)] * (x_exponent_order) + [Scalar(1)]
        x_exponent_poly = Polynomial(x_exponent_values_in_coeff, Basis.MONOMIAL)
        # commit x_exponent_poly
        x_exponent_comm_2 = setup.commit_g2(x_exponent_poly)

        B_0_check_lhs = b.pairing(x_exponent_comm_2, B_0_comm_1)
        B_0_check_rhs = b.pairing(b.G2, P_comm_1)
        assert B_0_check_lhs == B_0_check_rhs, "Check 2 failed: B0 degree check"
        print("=== Finished Check 2: B_0 has the appropriate degree ===")

        ### Check 3: 3.6 (c) ###
        print("=== Started Check 3: batched KZG check for the correctness of b_0_at_gamma, f_at_gamma, Q_b_at_gamma ===")
        # compute c
        b_at_0 = group_order_N * a_at_0 / group_order_n
        Z_H_gamma = gamma ** group_order_n - 1
        b_gamma = b_0_at_gamma * gamma + b_at_0
        Q_b_at_gamma = (b_gamma * (f_at_gamma + beta) - Scalar(1)) / Z_H_gamma
        # (a) both P and V compute v
        v = self.rlc(b_0_at_gamma, f_at_gamma, Q_b_at_gamma, eta)
        # v computes c
        c = ec_lincomb([
            (B_0_comm_1, 1),
            (f_comm_1, eta),
            (Q_B_comm_1, eta * eta)
        ])

        # batched KZG check for the correctness of b_0_at_gamma, f_at_gamma, Q_b_at_gamma
        comb_batch = ec_lincomb([
            (c, 1),
            (b.G1, -v),
            (pi_gamma, gamma)
        ])
        batch_check_lhs = b.pairing(b.G2, comb_batch)
        batch_check_rhs = b.pairing(powers_of_x2[1], pi_gamma)
        assert batch_check_lhs == batch_check_rhs, "Check 3 failed: batched KZG check for the correctness of b_0_at_gamma, f_at_gamma, Q_b_at_gamma"
        print("=== Finished Check 3: batched KZG check for the correctness of b_0_at_gamma, f_at_gamma, Q_b_at_gamma ===")

        ### Check 4: 3.7 (b) ###
        print("=== Started Check 4: KZG check for the correctness of a_at_0 ===")
        a_0_check_comb = ec_lincomb([
            # A_comm_1 - a_at_0
            (A_comm_1, 1),
            (b.G1, -a_at_0)
        ])
        x_poly = Polynomial([Scalar(0), Scalar(1)], Basis.MONOMIAL)
        x_comm_2 = setup.commit_g2(x_poly)
        assert x_comm_2 == powers_of_x2[1], "failed x commitment ==========="
        one_poly = Polynomial([Scalar(1)], Basis.MONOMIAL)
        one_comm_2 = setup.commit_g2(one_poly)
        assert one_comm_2 == b.G2, "failed 1 commitment ==========="
        a_0_check_lhs = b.pairing(one_comm_2, a_0_check_comb)
        a_0_check_rhs = b.pairing(x_comm_2, a_0_comm_1)

        assert a_0_check_lhs == a_0_check_rhs, "Check 4 failed: a_at_0 check"
        print("=== Finished Check 4: KZG check for the correctness of a_at_0 ===")

        print("Finished to verify proof")
        return True

    # Compute challenges (should be same as those computed by prover)
    def compute_challenges(
        self, proof
    ) -> tuple[Scalar, Scalar, Scalar, Scalar, Scalar, Scalar]:
        transcript = Transcript(b"plonk")
        beta = transcript.round_1(proof.msg_1)
        gamma, eta = transcript.round_2(proof.msg_2)

        return beta, gamma, eta

    def verify_commitment(self, proof, W, W_quot_key, eval_key, zeta):
        W_quot = proof[W_quot_key]
        eval = proof[eval_key]
        ec_comb = ec_lincomb(
            [
                (W, 1),
                (W_quot, zeta),
                (b.G1, -eval),
            ]
        )

        assert b.pairing(self.X_2, W_quot) == b.pairing(b.G2, ec_comb)
        print(f"Done KZG10 commitment check for {eval_key} polynomial")

    def rlc(self, term_1, term_2, term_3, eta):
        return term_1 + term_2 * eta + term_3 * eta * eta
