from dataclasses import dataclass
from typing import List

from caulk_multiple_setup import CaulkMultipleSetup
from caulk_multiple_transcript import Message1, Message2, Message3, Transcript
from src.caulk.util import vanishing_poly
from src.common_util.curve_optimized import Scalar, G1Point
from src.common_util.poly_optimized import Polynomial, Basis


@dataclass
class Proof:
    cm: G1Point
    m: int
    message1: Message1
    message2: Message2
    message3: Message3

    test_g1_P1: G1Point
    test_g1_P2: G1Point


class CaulkMultipleProver:
    def __init__(self, setup: CaulkMultipleSetup):
        self.setup = setup

    def prove(self, values: List[Scalar]) -> Proof:
        # reference to rust_implementation/src/multi/prover.rs

        transcript = Transcript(b"caulk-multiple")
        N = self.setup.N
        m = len(values)
        # m should be power of 2
        assert m & (m - 1) == 0
        phi_poly = Polynomial(values, Basis.LAGRANGE).ifft()
        cm = self.setup.kzgSetup.commit_G1(phi_poly)
        z_v_m_poly = vanishing_poly(m)

        positions = [self.setup.c.index(v) for v in values]

        blinders = [Scalar.get_random() for _ in range(7)]
        unique_positions = list(set(positions))

        # z_I(X) = r1 prod_{i in I} (X - w^i)
        z_I_poly = Polynomial([blinders[0]], Basis.MONOMIAL)
        for i in unique_positions:
            z_I_poly *= Polynomial([Scalar(-self.setup.roots_N[i]), Scalar(1)], Basis.MONOMIAL)

        # C_I(X) = (r_2+r_3X + r4X^2)*Z_I(X) + sum_j c_j*tau_j(X)
        tau_polys = self.tau_polys(unique_positions)
        C_I_poly = sum([tau_polys[i] * self.setup.c[j] for i, j in enumerate(unique_positions)],
                       start=Polynomial([Scalar(0)]))
        C_I_blinder_poly = Polynomial([blinders[1], blinders[2], blinders[3]], Basis.MONOMIAL)
        C_I_poly += C_I_blinder_poly * z_I_poly

        H1_poly = (self.setup.c_poly - C_I_poly) / z_I_poly
        g2_H1 = self.setup.kzgSetup.commit_G2(H1_poly)

        # u(X) = sum_j w^{i_j} mu_j(X) + (r5 + r6 X + r7 X^2) z_{Vm}(X)
        u_poly = Polynomial([self.setup.roots_N[i_j] for i_j in positions], Basis.LAGRANGE).ifft()
        u_poly += Polynomial([blinders[4], blinders[5], blinders[6]], Basis.MONOMIAL) * z_v_m_poly
        g1_u = self.setup.kzgSetup.commit_G1(u_poly)

        g1_C_I = self.setup.kzgSetup.commit_G1(C_I_poly)
        g1_Z_I = self.setup.kzgSetup.commit_G1(z_I_poly)

        message1 = Message1(g1_C_I, g1_Z_I, g1_u, g2_H1)

        # using Fiat-Shamir heuristic to get challenge chi
        chi = transcript.round1(message1)

        z_I_u_poly = z_I_poly.compose(u_poly)
        c_I_u_poly = C_I_poly.compose(u_poly)

        tmp_poly = z_I_u_poly + (c_I_u_poly - phi_poly) * chi
        H2_poly = tmp_poly / z_v_m_poly
        g1_H2 = self.setup.kzgSetup.commit_G1(H2_poly)

        message2 = Message2(g1_H2)

        # using Fiat-Shamir heuristic to get challenge alpha
        alpha = transcript.round2(message2)

        p1_poly = z_I_poly + C_I_poly * chi
        test_g1_P1 = self.setup.kzgSetup.commit_G1(p1_poly)

        p2_poly = -H2_poly * z_v_m_poly.eval(alpha) + z_I_u_poly.eval(alpha) + chi * (
                -phi_poly + c_I_u_poly.eval(alpha))
        test_g1_P2 = self.setup.kzgSetup.commit_G1(p2_poly)

        v1, pi_1 = self.setup.kzgSetup.open(u_poly, alpha)
        v2, pi_2 = self.setup.kzgSetup.open(p1_poly, v1)
        self.setup.kzgSetup.verify(test_g1_P1, pi_2, v1, v2)
        v3, pi_3 = self.setup.kzgSetup.open(p2_poly, alpha)
        assert v3 == 0
        self.setup.kzgSetup.verify(test_g1_P2, pi_3, alpha, Scalar(0))

        message3 = Message3(v1, v2, pi_1, pi_2, pi_3)

        return Proof(cm, m, message1, message2, message3, test_g1_P1, test_g1_P2)

    def tau_polys(self, unique_positions: List[int]):
        polys = []
        for j in unique_positions:
            tau_j = Polynomial([Scalar(1)])
            for i in unique_positions:
                if i != j:
                    tau_j *= Polynomial([Scalar(-self.setup.roots_N[i]), Scalar(1)], Basis.MONOMIAL)
                    tau_j /= Scalar(self.setup.roots_N[j] - self.setup.roots_N[i])
            polys.append(tau_j)

        return polys

    def prove_unity(self, a: Scalar, i: int):
        pass
