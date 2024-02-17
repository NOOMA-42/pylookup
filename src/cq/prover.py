from dataclasses import dataclass
from collections import Counter
from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar
from src.cq.setup import *
from src.cq.transcript import Transcript, Message1, Message2, Message3


@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3

    def flatten(self):
        proof = {}
        # msg_1
        proof["m_comm_1"] = self.msg_1.m_comm_1
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
        self.beta = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.gamma, self.eta = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()

        return Proof(msg_1, msg_2, msg_3)

    # Prover sends commitment of m(X)
    def round_1(self, witness) -> Message1:
        setup = self.setup
        duplicates = dict(Counter(witness))
        self.m_values = [Scalar(duplicates.get(val, 0)) for val in self.table]
        print("m_values: ", self.m_values)
        m_poly = Polynomial(self.m_values, Basis.LAGRANGE)
        self.m_poly = m_poly.ifft()
        # commit A(X) on G1(by default)
        self.m_comm_1 = setup.commit_g1(self.m_poly)
        print("Commitment of m(X): ", self.m_comm_1)

        return Message1(self.m_comm_1)

    # Prover sends commitment of A(X), Q_A(X), B_0(X), Q_B(X), P(X)
    def round_2(self) -> Message2:
        setup = self.setup
        group_order_N = self.group_order_N
        group_order_n = self.group_order_n
        beta = self.beta
        m_values = self.m_values
        t_values = self.t_values
        table_len = group_order_N
        # 1. commit A(X): Step 1-3 in the paper
        # 1.a. compute A_i values
        self.A_values = []
        for i, t_i in enumerate(t_values):
            A_i = m_values[i]/(beta + t_i)
            self.A_values.append(A_i)
            # sanity check
            assert A_i == m_values[i]/(beta + t_i), "A: not equal"
        print("A_values: ", self.A_values)

        # 1.b. compute A(X) from A_i values
        A_poly = Polynomial(self.A_values, Basis.LAGRANGE)
        # A(X) in coefficient form
        self.A_poly = A_poly.ifft()
        assert A_poly.barycentric_eval(Scalar(0))  == self.A_poly.coeff_eval(Scalar(0)), "A value at 0 should be equal"

        # 1.c. commit A(X)
        self.A_comm_1 = setup.commit_g1(self.A_poly)
        print("Commitment of A(X): ", self.A_comm_1)

        # 2. commit Q_A(X): Step 4 in the paper
        # 2.a. T(X) in lagrange form
        T_poly = Polynomial(t_values, Basis.LAGRANGE)
        # T(X) in coefficient form
        self.T_poly = T_poly.ifft()
        # 2.b. vanishing polynomial: X^N - 1, N = group_order_N - 1
        ZV_array = [Scalar(-1)] + [Scalar(0)] * (group_order_N - 1) + [Scalar(1)]
        # vanishing polynomial: X^n - 1, N = group_order_n - 1
        ZH_array = [Scalar(-1)] + [Scalar(0)] * (group_order_n - 1) + [Scalar(1)]
        # vanishing polynomial in coefficient form
        ZH_poly = Polynomial(ZH_array, Basis.MONOMIAL)

        t_poly_coeffs = self.T_poly.values
        print("\nt_poly_coeffs: ", t_poly_coeffs)
        print("\nsetup.powers_of_x: ", setup.powers_of_x)
        srs = setup.powers_of_x[:len(t_poly_coeffs)]
        print("\nsrs: ", srs)

        roots = Scalar.roots_of_unity(len(self.A_values))
        print("roots: ", roots)
        for i in range(len(t_values)):
            eval = self.T_poly.coeff_eval(roots[i])
            print("=====> eval: ", eval)  # 1, 2, 3, 4
            assert eval == t_values[i]

        # Precomputation happens here with FK algorithm
        Q_T_comm_poly_coeffs = setup.Q_T_comm_poly_coeffs
        print("\n------> Q_T_comm_poly_coeffs: \n", Q_T_comm_poly_coeffs)

        Q_A_Comm = b.Z1
        for i in range(table_len):
            K_T_Comm = b.Z1
            root = Scalar(1)
            for j in range(table_len):
                K_T_Comm = b.add(K_T_Comm, b.multiply(
                    Q_T_comm_poly_coeffs[j], root.n))
                root = root * roots[i]
            A_val = self.A_values[i].n
            scale = roots[i]/table_len
            # Compute Quotient polynomial commitment of T(X)
            Q_T_Comm = b.multiply(K_T_Comm, scale.n)
            A_times_Q_T_Comm = b.multiply(Q_T_Comm, A_val)
            # Do the accumulation
            Q_A_Comm = b.add(Q_A_Comm, A_times_Q_T_Comm)


        self.Q_A_comm_1 = Q_A_Comm
        print("\nCommitment of Q_A(X):  \n", self.Q_A_comm_1)

        print(" \n*********** Finish to precomputed commitment of Q_A(X) **********")

        # 3. commit B_0(X): Step 5-7 in the paper
        # 3.a. compute B_0_i values
        self.B_values = []
        f_values = self.lookup_table
        for i, f_i in enumerate(f_values):
            B_i = 1 / (beta + f_i)
            self.B_values.append(B_i)
            # sanity check
            assert B_i == 1 / (beta + f_i), "B: not equal"
        # 3.b. compute B_0(X) from B_0_i values, B_0(X) = (B(X) - B(0)) / X
        B_poly = Polynomial(self.B_values, Basis.LAGRANGE)
        # in coefficient form
        self.B_poly = B_poly.ifft()
        # f(X) = X, coefficient form: [0, 1]
        self.x_poly = Polynomial([Scalar(0), Scalar(1)], Basis.MONOMIAL)
        B_at_0 = self.B_poly.coeff_eval(Scalar(0))
        print("B_at_0: ", B_at_0)
        self.B_0_poly: Polynomial = (self.B_poly - B_at_0) / self.x_poly
        # sanity check
        for i in range(group_order_n):
            point = self.roots_of_unity_n[i]
            b_value = self.B_poly.coeff_eval(point)
            b_0_value = self.B_0_poly.coeff_eval(point)
            assert b_value == self.B_values[i], "B_value and self.B_values[i]: Not equal"
            assert b_0_value == (b_value - B_at_0) / point, "B_0: Not equal"
        # 3.c. commit B_0(X)
        self.B_0_comm_1 = setup.commit_g1(self.B_0_poly)
        print("Commitment of B_0(X): ", self.B_0_comm_1)

        # 4. commit Q_B(X): Step 9 in the paper
        # 4.a. f(X) in coefficient form
        f_poly = Polynomial(f_values, Basis.LAGRANGE)
        # in coefficient form
        self.f_poly = f_poly.ifft()
        # commit f(X)
        self.f_comm_1 = setup.commit_g1(self.f_poly)
        print("Commitment of f(X): ", self.f_comm_1)

        # sanity check
        for i, B_i in enumerate(self.B_values):
            point = self.roots_of_unity_n[i]
            b_value = self.B_poly.coeff_eval(point)
            f_value = self.f_poly.coeff_eval(point)
            assert b_value == 1 / (beta + f_value) , "B quotient: Not equal"
        # 4.b. Q_B(X) in coefficient form
        self.Q_B_poly = (self.B_poly * (self.f_poly + beta) - Scalar(1)) / ZH_poly
        # 4.c. commit Q_B(X): Step 9 in the paper
        self.Q_B_comm_1 = setup.commit_g1(self.Q_B_poly)
        print("Commitment of Q_B(X): ", self.Q_B_comm_1)

        # 5. commit P(X): Step 10 in the paper
        # N - 1 - (n - 2)
        x_exponent_order = group_order_N - 1 - (group_order_n - 2)
        x_exponent_values_in_coeff = [Scalar(0)] * (x_exponent_order) + [Scalar(1)]
        x_exponent_poly = Polynomial(x_exponent_values_in_coeff, Basis.MONOMIAL)
        self.P_poly = self.B_0_poly * x_exponent_poly
        # 5.c. commit P(X)
        self.P_comm_1 = setup.commit_g1(self.P_poly)
        print("Commitment of P(X): ", self.P_comm_1)

        return Message2(
            self.A_comm_1,
            self.Q_A_comm_1,
            self.f_comm_1,
            self.B_0_comm_1,
            self.Q_B_comm_1,
            self.P_comm_1
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

    # random linear combination
    def rlc(self, term_1, term_2, term_3):
        return term_1 + term_2 * self.eta + term_3 * self.eta * self.eta
