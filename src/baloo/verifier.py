import py_ecc.bn128 as b
from dataclasses import dataclass
from src.baloo.transcript import Transcript
from src.common_util.poly import Polynomial, Basis, PolyUtil
from src.common_util.curve import ec_lincomb, G1Point, G2Point, Scalar

poly_util = PolyUtil()

@dataclass
class VerificationKey:
    x2: G2Point
    # commitment to public table
    t_comm_1: G1Point
    # commitment to vanishing polynomial of public table
    z_H_comm_1: G1Point

    # m: len(lookup)
    def verify_proof(self, proof_instance, setup, lookup_table_len) -> bool:
        print("Start to verify proof")
        alpha, beta, gamma, zeta = self.compute_challenges(proof_instance)
        x2 = self.x2
        t_comm_1 = self.t_comm_1
        z_H_comm_1 = self.z_H_comm_1
        proof = proof_instance.flatten()
        # get commitments
        z_I_comm_2 = proof["z_I_comm_2"]
        v_comm_1 = proof["v_comm_1"]
        t_I_comm_1 = proof["t_I_comm_1"]
        phi_comm_1 = proof["phi_comm_1"]
        D_comm_1 = proof["D_comm_1"]
        R_comm_1 = proof["R_comm_1"]
        Q_D_comm_1 = proof["Q_D_comm_1"]
        E_comm_1 = proof["E_comm_1"]
        Q_E_comm_1 = proof["Q_E_comm_1"]
        v1 = proof["v1"]
        v2 = proof["v2"]
        v3 = proof["v3"]
        v4 = proof["v4"]
        v5 = proof["v5"]
        w1_comm_1 = proof["w1_comm_1"]
        w2_comm_1 = proof["w2_comm_1"]
        w3_comm_1 = proof["w3_comm_1"]
        w4_comm_1 = proof["w4_comm_1"]
        w5_comm_1 = proof["w5_comm_1"]
        w6_comm_1 = proof["w6_comm_1"]
        scalar_one = Scalar(1)
        # length of SRS on G1
        d = len(setup.powers_of_x)
        m = lookup_table_len
        # calculate commitment [P_D(X)]1
        P_D_comm_1 = ec_lincomb([
            (t_I_comm_1, v1),
            (b.G1, -v2),
            (R_comm_1, -scalar_one),
            (Q_D_comm_1, v4),
        ])
        # X^m - 1
        z_V_values = poly_util.vanishing_poly(m)
        z_V_poly = Polynomial(z_V_values, Basis.MONOMIAL)
        z_V_poly_at_zeta = z_V_poly.coeff_eval(Scalar(zeta))
        print("z_V_poly_at_zeta: ", z_V_poly_at_zeta)
        print("d, m: ", d, m)

        # calculate commitment [P_E(X)]1
        P_E_comm_1 = ec_lincomb([
            (b.G1, v5 * beta),
            (E_comm_1, v5),
            (v_comm_1, -v4 / v3),
            (Q_E_comm_1, -z_V_poly_at_zeta),
        ])

        # 1. verify subtable
        subtable_lhs = ec_lincomb([
            (t_comm_1, scalar_one),
            (t_I_comm_1, -scalar_one),
            (z_H_comm_1, gamma),
        ])
        subtable_rhs = ec_lincomb([
            (w5_comm_1, scalar_one),
            (w6_comm_1, gamma),
        ])
        assert b.pairing(b.G2, subtable_lhs) == b.pairing(
            z_I_comm_2, subtable_rhs), "subtable paring check failed"
        print("Finished to verify: subtable")

        # 2. verify w1 for X = α
        # X^(d-m+1)
        x_exponent_poly_1 = poly_util.x_exponent_poly(d - m + 1)
        x_exponent_1_comm_2 = setup.commit_g2(x_exponent_poly_1)
        w1_rhs1 = ec_lincomb([
            (E_comm_1, scalar_one),
            (phi_comm_1, gamma),
            (b.G1, -v1),
            (b.G1, -gamma * v2),
        ])
        w1_rhs2 = ec_lincomb([
            (w1_comm_1, alpha),
        ])
        assert b.pairing(x2, w1_comm_1) == b.pairing(
            x_exponent_1_comm_2, w1_rhs1) * b.pairing(b.G2, w1_rhs2), "w1 paring check failed"
        print("Finished to verify: w1")

        # 3. verify w2 for X = 0
        # X^(d-m+2)
        x_exponent_poly_2 = self.x_exponent_poly(d - m + 2)
        x_exponent_2_comm_1 = setup.commit_g1(x_exponent_poly_2)
        # X^m
        x_exponent_poly_3 = self.x_exponent_poly(m)
        x_exponent_3_comm_1 = setup.commit_g1(x_exponent_poly_3)

        w2_rhs1 = ec_lincomb([
            (b.G1, scalar_one),
            (x_exponent_2_comm_1, gamma ** 2),
        ])
        w2_rhs2 = ec_lincomb([
            (R_comm_1, gamma ** 3),
            (x_exponent_3_comm_1, -gamma ** 3),
        ])
        w2_rhs3 = ec_lincomb([
            (R_comm_1, gamma),
            (b.G1, v3),
        ])
        assert b.pairing(x2, w2_comm_1) == b.pairing(z_I_comm_2, w2_rhs1) * b.pairing(
            x_exponent_2_comm_1, w2_rhs2) * b.pairing(b.G2, w2_rhs3), "w1 paring check failed"
        print("Finished to verify: w2")

        # 4. verify w3 for X = β
        w3_rhs1 = ec_lincomb([
            (D_comm_1, scalar_one),
            (w3_comm_1, beta),
            (b.G1, -gamma * v4),
            (b.G1, -v1),
            (P_D_comm_1, gamma ** 2)
        ])
        w3_rhs2 = ec_lincomb([
            (b.G1, gamma),
        ])
        assert b.pairing(x2, w3_comm_1) == b.pairing(
            b.G2, w3_rhs1) * b.pairing(z_I_comm_2, w3_rhs2), "w3 paring check failed"
        print("Finished to verify: w3")

        # 5. verify w4 for X = ζ
        w4_rhs = ec_lincomb([
            (E_comm_1, scalar_one),
            (w4_comm_1, zeta),
            (P_E_comm_1, gamma)
            (b.G1, -gamma * v5),
        ])
        assert b.pairing(x2, w4_comm_1) == b.pairing(
            b.G2, w4_rhs), "w4 paring check failed"
        print("Finished to verify: w4")

        print("Finished to verify proof")
        return True

    # Compute challenges (should be same as those computed by prover)
    def compute_challenges(
        self, proof
    ) -> tuple[Scalar, Scalar, Scalar, Scalar, Scalar, Scalar]:
        transcript = Transcript(b"plonk")
        alpha, beta = transcript.round_1(proof.msg_1)
        gamma, zeta = transcript.round_2(proof.msg_2)

        return alpha, beta, gamma, zeta

