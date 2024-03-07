from caulk_single_setup import CaulkSingleSetup
from src.caulk.single.caulk_single_prover import Proof, CaulkSingleProver, Proof_pederson, Proof_unity
from src.caulk.util import hash_ec_points, vanishing_poly, lagrange_polys, single_term_poly
from src.common_util.curve_optimized import Scalar, ec_mul, G1, ec_add, ec_sub, ec_pairing, G2, G1Point, ec_eq, \
    ec_lincomb
from src.common_util.poly_optimized import Polynomial


class CaulkSingleVerifier:
    setup: CaulkSingleSetup

    def __init__(self, setup: CaulkSingleSetup):
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

        self.verify_unity(proof.proof_unity)
        print("unity check passed")

    def verify_unity(self, proof: Proof_unity):
        sigma = setup.roots_n[1]
        alpha = proof.alpha
        alpha_1 = setup.roots_n[-1] * alpha
        alpha_2 = setup.roots_n[-2] * alpha
        rho_polys = lagrange_polys(setup.n)
        rho_alphas = [poly.eval(alpha) for poly in rho_polys]
        poly_prod_alpha = setup.poly_prod.eval(alpha)

        z_vn_alpha = vanishing_poly(setup.n).eval(alpha)

        # first verify v1 and v2
        setup.kzgSetup.verify(proof.f_comm, proof.pi_f_1, alpha_1, proof.v1)
        setup.kzgSetup.verify(proof.f_comm, proof.pi_f_2, alpha_2, proof.v2)

        # prepare to calculate commitment of p_alpha
        f_power = rho_alphas[0] + rho_alphas[1]
        f_power += rho_alphas[2] * (Scalar(1) - sigma)
        f_power += rho_alphas[3]
        f_power += rho_alphas[4] * proof.v1
        f_power += poly_prod_alpha

        g1_power = rho_alphas[2] * (proof.v1 - proof.v2)
        g1_power += rho_alphas[3] * (proof.v2 - proof.v1 * sigma)
        g1_power += rho_alphas[4] * -proof.v2
        g1_power += poly_prod_alpha * -proof.v1 ** 2
        g1_power += rho_alphas[setup.n - 1] * (proof.v1 - Scalar(1))

        p_comm = ec_lincomb(
            ((proof.h_comm, -z_vn_alpha),
             (proof.f_comm, f_power),
             (G1, g1_power)))

        x_d_poly = single_term_poly(setup.kzgSetup.length - 2)
        x_d_comm = setup.kzgSetup.commit_G1(x_d_poly)
        g1_term = ec_mul(G1, -rho_alphas[0] - rho_alphas[1])
        g1_term = ec_add(g1_term, ec_mul(x_d_comm, z_vn_alpha))
        LHS = ec_pairing(proof.z_comm, g1_term)
        LHS = LHS * ec_pairing(G2, p_comm)

        x_minus_alpha_comm = setup.kzgSetup.commit_G2(Polynomial([Scalar(-alpha), Scalar(1)]))
        RHS = ec_pairing(x_minus_alpha_comm, proof.pi_p_alpha)

        assert LHS == RHS

    def verify_pederson(self, cm: G1Point, proof: Proof_pederson):
        c = hash_ec_points(proof.R, cm)
        LHS = ec_add(proof.R, ec_mul(cm, c))
        RHS = ec_mul(G1, proof.t1 + self.setup.h * proof.t2)
        assert ec_eq(LHS, RHS)


if __name__ == "__main__":
    setup = CaulkSingleSetup.example_setup()
    prover = CaulkSingleProver(setup)
    proof = prover.prove(Scalar(2))

    verifier = CaulkSingleVerifier(setup)
    verifier.verify(proof)
