from dataclasses import dataclass
import numpy as np
from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar
from src.baloo.setup import *
from src.baloo.transcript import Transcript, Message1, Message2, Message3, Message4


@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3
    msg_4: Message4

    def flatten(self):
        proof = {}
        # msg_1
        proof["phi_comm_1"] = self.msg_1.phi_comm_1
        # msg_2
        proof["A_comm_1"] = self.msg_2.A_comm_1
        proof["Q_A_comm_1"] = self.msg_2.Q_A_comm_1
        proof["f_comm_1"] = self.msg_2.f_comm_1
        proof["B_0_comm_1"] = self.msg_2.B_0_comm_1
        proof["Q_B_comm_1"] = self.msg_2.Q_B_comm_1
        proof["P_comm_1"] = self.msg_2.P_comm_1
        # msg_3
        proof["b_0_at_gamma"] = self.msg_3.b_0_at_gamma
        proof["f_at_gamma"] = self.msg_3.f_at_gamma
        proof["a_at_0"] = self.msg_3.a_at_0
        proof["pi_gamma"] = self.msg_3.pi_gamma
        proof["a_0_comm_1"] = self.msg_3.a_0_comm_1

        return proof

@dataclass
class Prover:
    group_order_N: int
    group_order_n: int
    setup: Setup
    table: list

    def __init__(self, setup: Setup, table: list, group_order_n: int):
        self.setup = setup
        self.table = table
        self.t_values = [Scalar(val) for val in self.table]
        self.group_order_N = len(table)
        self.group_order_n = group_order_n
        self.roots_of_unity_N = Scalar.roots_of_unity(self.group_order_N)
        self.roots_of_unity_n = Scalar.roots_of_unity(group_order_n)
        self.powers_of_x = setup.powers_of_x
    def prove(self, witness) -> Proof:
        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")
        self.lookup_table = [Scalar(val) for val in witness]

        # Round 1
        msg_1 = self.round_1(witness)
        self.alpha = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.beta = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()
        self.rho, self.gamma = transcript.round_3(msg_3)

        msg_4 = self.round_4()

        return Proof(msg_1, msg_2, msg_3, msg_4)

    """
    table = [1, 2, 3, 4, 5, 6, 7, 8]
    witness = [3, 1, 2, 1]
    t:  [1, 2, 3] # choose non-duplicated elements from witness
    I:  [0, 1, 2] # get indexes from table
    k = len(I) = 3
    vanishing polynomial z_I(X) = (X - I_0)(X - I_1)(X - I_2)
                                = (X - 0)(X - 1)(X - 2)
    M * t = witness
    M = [
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0],
        [1, 0, 0],
    ]
    m = len(witness) # 4
    col[0] = M[0].index(1)
    col[1] = M[1].index(1)
    col[2] = M[2].index(1)
    col[3] = M[3].index(1)
    col: [2, 0, 1, 0]
    v_values[0] = 1 / xi_values[col[0]]
    v_values[1] = 1 / xi_values[col[1]]
    v_values[2] = 1 / xi_values[col[2]]
    v_values[3] = 1 / xi_values[col[3]]
    """
    # Output π1 = [z_I]_2 = [z_I(x)]_2, [v]_1 = [v(x)]_1, t = [t(x)]_1
    def round_1(self, witness) -> Message1:
        setup = self.setup
        # calculate lookup table polynomial φ(X)
        phi_values = [Scalar(val) for val in witness]
        phi_poly = Polynomial(phi_values, Basis.LAGRANGE)
        # coefficient form
        self.phi_poly = phi_poly.ifft()
        # commit phi(X) on G1(by default)
        self.phi_comm_1 = setup.commit_g1(self.phi_poly)
        print("Commitment of phi(X): ", self.phi_comm_1)
        # remove duplicated elements
        t_values = list(set(witness))
        # transform to Scalar
        t_values = [Scalar(elem) for elem in t_values]
        print("t: ", t_values)
        # I: the index of t_values elements in public table
        I_values = [Scalar(self.table.index(elem)) for elem in t_values]
        print("I: ", I_values)

        k = len(I_values)
        # vanishing polynomial in coefficient form
        z_I_poly = Polynomial([Scalar(1)], Basis.MONOMIAL)
        # multiple all root polynomial: (X - I_0)(X - I_1)(X - I_2)...
        for i in range(k):
            z_I_poly = z_I_poly * Polynomial(
                [-I_values[i], Scalar(1)], Basis.MONOMIAL)
        print("z_I_poly.values: ", z_I_poly.values)
        self.z_I_poly = z_I_poly

        root_of_unity = Scalar.root_of_unity(k)
        m = len(witness)  # 4
        col_values = []
        v_values = []
        xi_values = []
        for i in range(m):
            col_i = t_values.index(witness[i])  # find the index of 1 in jth row of M
            col_values.append(col_i)
            # ξ（Xi）
            xi = root_of_unity ** col_i
            xi_values.append(xi)
            v = 1 / xi
            v_values.append(v)
        # ξ（Xi）
        print("xi_values: ", xi_values)
        print("col_values: ", col_values)
        assert np.array_equal([t_values[elem] for elem in col_values], witness)
        print("v_values: ", v_values)
        v_poly = Polynomial(v_values, Basis.LAGRANGE)
        t_poly = Polynomial(t_values, Basis.LAGRANGE)

        self.v_poly = v_poly.ifft()
        self.t_poly = t_poly.ifft()

        # commit
        self.z_I_comm_2 = setup.commit_g2(self.z_I_poly)
        self.v_comm_1 = setup.commit_g1(self.v_poly)
        self.t_comm_1 = setup.commit_g1(self.t_poly)
        print("self.z_I_comm_2: ", self.z_I_comm_2)
        print("self.v_comm_1: ", self.v_comm_1)
        print("self.t_comm_1: ", self.t_comm_1)
        return Message1(self.z_I_comm_2, self.v_comm_1, self.t_comm_1)


    def round_2(self) -> Message2:
        setup = self.setup
        alpha = self.alpha
        tau_col_values = []
        tau_col_poly = Polynomial(tau_col_values, Basis.LAGRANGE)
        tau_col_at_0 = tau_col_poly.barycentric_eval(Scalar(0))
        return Message2(
        )

    def round_3(self) -> Message3:
        # 1. V sends random γ,η ∈ F.: Step 1 in the paper
        setup = self.setup
        beta = self.beta
        gamma = self.gamma
        group_order_N = self.group_order_N
        group_order_n = self.group_order_n

        # 2. compute b_0_at_gamma: Step 2 in the paper
        b_0_at_gamma = self.B_0_poly.coeff_eval(gamma)
        # compute f_at_gamma
        f_at_gamma = self.f_poly.coeff_eval(gamma)
        # 3. compute a_at_0: Step 3 in the paper
        a_at_0 = self.A_poly.coeff_eval(Scalar(0))
        # 4. compute b_at_0: Step 4 in the paper
        b_at_0 = group_order_N * a_at_0 / group_order_n
        # 5. compute b_at_gamma, and Q_b_at_gamma: Step 5 in the paper
        Z_H_at_gamma = gamma ** group_order_n - 1
        b_at_gamma = b_0_at_gamma * gamma + b_at_0
        Q_b_at_gamma = (b_at_gamma * (f_at_gamma + beta) - Scalar(1)) / Z_H_at_gamma

        # 6. batch KZG check: Step 6 in the paper
        # (a) both P and V compute v
        v = self.rlc(b_0_at_gamma, f_at_gamma, Q_b_at_gamma)
        # (b) compute commitment: pi_gamma = [h(X)]_1
        h_poly = (self.rlc(self.B_0_poly, self.f_poly, self.Q_B_poly) - v) / (self.x_poly - gamma)
        pi_gamma = setup.commit_g1(h_poly)

        # 3.7 commit A_0(X): Step 7 in the paper
        # (a) compute a_0_comm_1
        a_0_poly = (self.A_poly - a_at_0) / self.x_poly
        a_0_comm_1 = setup.commit_g1(a_0_poly)
        print("Prover: a_0_comm_1: ", a_0_comm_1)

        return Message3(b_0_at_gamma, f_at_gamma, a_at_0, pi_gamma, a_0_comm_1)

    def round_4(self) -> Message4:
        # todo
        return Message3()

    # random linear combination
    def rlc(self, term_1, term_2, term_3):
        return term_1 + term_2 * self.eta + term_3 * self.eta * self.eta
