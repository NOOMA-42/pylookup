from dataclasses import dataclass
import numpy as np
from src.common_util.poly import Polynomial, Basis, InterpolationPoly, PolyUtil
from src.common_util.curve import Scalar
from src.baloo.setup import *
from src.baloo.transcript import Transcript, Message1, Message2, Message3

poly_util = PolyUtil()

@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3

    def flatten(self):
        proof = {}
        # msg_1
        proof["z_I_comm_2"] = self.msg_1.z_I_comm_2
        proof["v_comm_1"] = self.msg_1.v_comm_1
        proof["t_I_comm_1"] = self.msg_1.t_I_comm_1
        proof["phi_comm_1"] = self.msg_1.phi_comm_1
        # msg_2
        proof["D_comm_1"] = self.msg_2.D_comm_1
        proof["R_comm_1"] = self.msg_2.R_comm_1
        proof["Q_D_comm_1"] = self.msg_2.Q_D_comm_1
        proof["E_comm_1"] = self.msg_2.E_comm_1
        proof["Q_E_comm_1"] = self.msg_2.Q_E_comm_1
        # msg_3
        proof["v1"] = self.msg_3.v1
        proof["v2"] = self.msg_3.v2
        proof["v3"] = self.msg_3.v3
        proof["v4"] = self.msg_3.v4
        proof["v5"] = self.msg_3.v5
        proof["w1_comm_1"] = self.msg_3.w1_comm_1
        proof["w2_comm_1"] = self.msg_3.w2_comm_1
        proof["w3_comm_1"] = self.msg_3.w3_comm_1
        proof["w4_comm_1"] = self.msg_3.w4_comm_1
        proof["w5_comm_1"] = self.msg_3.w5_comm_1
        proof["w6_comm_1"] = self.msg_3.w6_comm_1

        return proof

@dataclass
class Prover:
    setup: Setup
    table: list

    def __init__(self, setup: Setup, table: list):
        self.setup = setup
        self.table = table
        self.t_values = [Scalar(val) for val in table]
        self.t_poly = Polynomial(self.t_values, Basis.LAGRANGE).ifft()
        self.roots_of_unity_N = Scalar.roots_of_unity(len(table))
        # vanishing polynomial: X^N - 1, N = len(table)
        z_H_array = poly_util.vanishing_poly(len(table))
        # in coefficient form
        self.z_H_poly = Polynomial(z_H_array, Basis.MONOMIAL)

    def prove(self, lookup) -> Proof:
        self.m = len(lookup)

        # Initialise Fiat-Shamir transcript
        transcript = Transcript(b"plonk")
        self.lookup_table = [Scalar(val) for val in lookup]

        # Round 1
        msg_1 = self.round_1(lookup)
        self.alpha, self.beta = transcript.round_1(msg_1)

        # Round 2
        msg_2 = self.round_2()
        self.gamma, self.zeta = transcript.round_2(msg_2)

        # Round 3
        msg_3 = self.round_3()

        return Proof(msg_1, msg_2, msg_3)

    """
    H = [1, ω, ω^2, ω^3, ω^4, ω^5, ω^6, ω^7]
    table = [1, 2, 3, 4, 5, 6, 7, 8]
    lookup = [3, 7, 3, 4]
    t_I:  [3, 4, 7] # choose non-duplicated elements from lookup
    I:  [2, 3, 6] # get indexes from table
    s ∈ I = [2, 3, 6]
    H_I = {ω^s} = [ω^2, ω^3, ω^6]
    k = len(I) = 3
    vanishing polynomial z_I(X) = (X - H_I_0)(X - H_I_1)(X - H_I_2)
                                = (X - ω^2)(X - ω^3)(X - ω^6)

    M * t = lookup
    M = [
        [1, 0, 0],
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0],
    ]
    m = len(lookup) # 4
    col[0] = M[0].index(1)
    col[1] = M[1].index(1)
    col[2] = M[2].index(1)
    col[3] = M[3].index(1)
    col: [0, 2, 0, 1]
    ξ(xi): xi = root_of_unity ** col_i
    get v(X)) values
    v_values[0] = 1 / H_I[col[0]]
    v_values[1] = 1 / H_I[col[1]]
    v_values[2] = 1 / H_I[col[2]]
    v_values[3] = 1 / H_I[col[3]]
    """
    def round_1(self, lookup) -> Message1:
        setup = self.setup
        m = self.m
        # calculate lookup table polynomial φ(X)
        phi_values = [Scalar(val) for val in lookup]
        phi_poly = Polynomial(phi_values, Basis.LAGRANGE)
        # coefficient form
        self.phi_poly = phi_poly.ifft()
        # commit phi(X) on G1
        self.phi_comm_1 = setup.commit_g1(self.phi_poly)
        print("Commitment of phi(X): ", self.phi_comm_1)
        # remove duplicated elements
        t_values = list(set(lookup))
        # transform to Scalar
        t_values = [Scalar(elem) for elem in t_values]
        print("t: ", t_values)
        # I: the index of t_values elements in sub table t_I
        I_values = [Scalar(self.table.index(elem)) for elem in t_values]
        print("I: ", I_values)
        # H_I = {ξ_i} , i = [1, k], ξ(Xi)
        H_I = [self.roots_of_unity_N[i.n] for i in I_values]  # len(H_I) = k
        self.H_I = H_I
        print("H_I: ", H_I)
        H_I_interp_poly = InterpolationPoly(H_I, t_values)

        # multiple all root polynomial: (X - I_0)(X - I_1)(X - I_2)...
        col_values = []
        v_values = []
        for i in range(m):
            # find the index of 1 in jth row of M
            col_i = t_values.index(lookup[i])
            col_values.append(col_i)
            col_i_root = H_I[col_i]
            # Note: v = 1 / col_i_root in paper
            # Here we use different construction that does not affect the verification
            v = col_i_root
            v_values.append(v)

        # vanishing polynomial in coefficient form
        z_I_poly = H_I_interp_poly.vanishing_poly()

        print("z_I_poly.values: ", z_I_poly.values)
        self.z_I_poly = z_I_poly

        print("col_values: ", col_values)
        self.col = col_values
        assert np.array_equal([t_values[elem] for elem in col_values], lookup)
        print("v_values: ", v_values)
        v_poly = Polynomial(v_values, Basis.LAGRANGE)
        # refer to section 5. can not use FFT due to it's not in multiplicative subgroup
        t_interp_poly = InterpolationPoly(H_I, t_values)
        self.t_I_poly = t_interp_poly.poly()

        self.v_poly = v_poly.ifft()

        # commit
        # π1 = ([z_I]_2 = [z_I(x)]_2, [v]_1 = [v(x)]_1, t = [t(x)]_1)
        self.z_I_comm_2 = setup.commit_g2(self.z_I_poly)
        self.v_comm_1 = setup.commit_g1(self.v_poly)
        self.t_I_comm_1 = setup.commit_g1(self.t_I_poly)

        return Message1(
            self.z_I_comm_2,
            self.v_comm_1,
            self.t_I_comm_1,
            self.phi_comm_1
        )

    """
    Calculate μ_i(X), i = [0, m - 1], V is a multiplicative subgroup
    col_i = col[i]
    v_i = V[i], v_col_i = V[col_i]
    μ_i(X) = z_V(X) / (z_V'(v_i) * (X - v_i)) = v_i / m * (X^m - 1) / (X - v_i)

    Calculate Normalized Lagrange Polynomial: τ_col(i)(X) / τ_col(i)(0):
    ξ_{col(i)}: h_i = H_I[col_i]
    normalized_lag_poly = τ_col(i)(X) / τ_col(i)(0)
                        = z_I(X) / z_I(0) * (-ξ_{col(i)}) / (X - ξ_{col(i)})


    Calculate R(X): Theorem 5 (Inner Product Polynomial Relation) on Baloo paper
    root = H_I[col_i]
    a_i = μ_i(α)
    b_i = t_I(root)
    R(X) = Σ_i(a_i * b_i * normalized_lag_poly(root)) - Σ_i(a_i * b_i)

    σ = Σ_i(a_i * b_i)
    φ(X): Polynomial(lookup, Basis.MONOMIAL), lookup table interpolation polynomial
    assert σ == φ(α)

    Calculate Q_D(X), Q_E(X)
    D(X) = Σ_i(μ_i(α) * normalized_lag_poly)
    E(X) = Σ_i(μ_i(X) * normalized_lag_poly(β))

    Calculate Q_D(X)
    D(X) * t_I(X) - φ(α) - R(X) = z_I(X) * Q_D(X)
    Q_D(X) = (D(X) * t_I(X) - φ(α) - R(X)) / z_I(X)

    Calculate Q_E(X)
    1) Baloo paper uses this construction
    v_i = 1 / H_I[col_i]
    v(X): interpolate polynomial with v_i values
    E(X) * (βv(X) - 1) + z_I(β) / z_I(0) = z_V(X) * Q_E(X)
    Q_E(X) = (E(X) * (βv(X) - 1) + z_I(β) / z_I(0)) / z_V(X)

    2) Our code uses this optimized construction:
    v_i = H_I[col_i]
    v(X): interpolate polynomial with v_i values
    E(X) * (β - v(X)) + v(X) * z_I(β) / z_I(0) = z_V(X) * Q_E(X)
    Q_E(X) = (E(X) * (β - v(X)) + v(X) * z_I(β) / z_I(0)) / z_V(X)
    """
    def round_2(self) -> Message2:
        setup = self.setup
        alpha = self.alpha
        beta = self.beta
        z_I_poly = self.z_I_poly
        t_I_poly = self.t_I_poly
        v_poly = self.v_poly
        phi_poly = self.phi_poly
        col = self.col
        H_I = self.H_I
        m = self.m
        zero_poly = Polynomial([Scalar(0)], Basis.MONOMIAL)

        V = Scalar.roots_of_unity(m)
        # X^m - 1
        z_V_values = poly_util.vanishing_poly(m)
        z_V_poly = Polynomial(z_V_values, Basis.MONOMIAL)
        # z_I(0)
        z_I_at_0 = z_I_poly.coeff_eval(Scalar(0))
        print("z_I_at_0: ", z_I_at_0, -z_I_at_0)
        # calculate D(X) = Σ_{0, m-1} μ_i(α) * τ^_{col(i)}(X)
        D_poly = zero_poly
        E_poly = zero_poly

        d_t_sum_poly = zero_poly
        sum = Scalar(0)
        for i in range(m):
            col_i = col[i]
            v_root = V[i]
            v_root_poly = Polynomial([-Scalar(v_root), Scalar(1)], Basis.MONOMIAL)
            # ξ_i
            root = H_I[col_i]
            # X - ξ_i
            x_root_poly = poly_util.root_poly(root)
            # Lagrange polynomial on V: μ_i(X)
            mu_poly = z_V_poly / v_root_poly * v_root / Scalar(m)
            # Normalized Lagrange Polynomial: τ_col(i)(X) / τ_col(i)(0)
            normalized_lag_poly: Polynomial = z_I_poly / \
                x_root_poly * (-root) / z_I_at_0
            # μ_i(α)
            mu_poly_at_alpha = mu_poly.coeff_eval(alpha)
            normalized_lag_poly_at_beta = normalized_lag_poly.coeff_eval(beta)
            print("mu_poly_at_alpha: ", alpha, mu_poly_at_alpha)
            # D(X) = Σ_i(μ_i(α) * normalized_lag_poly)
            D_poly += normalized_lag_poly * mu_poly_at_alpha
            # E(X) = Σ_i(μ_i(X) * normalized_lag_poly(β))
            E_poly += mu_poly * normalized_lag_poly_at_beta

            a_i = mu_poly_at_alpha
            b_i = t_I_poly.coeff_eval(root)
            print("a_i: ", a_i, -a_i)
            print("b_i: ", b_i, -b_i)
            sum_accu = a_i * b_i
            sum += sum_accu
            poly_accu = normalized_lag_poly * sum_accu
            d_t_sum_poly += poly_accu
        print("D_poly.values: ", D_poly.values)
        D_t_poly = D_poly * t_I_poly
        print("D_t_poly.values: ", D_t_poly.values)
        print("sum: ", sum)
        pha_poly_at_alpha = phi_poly.coeff_eval(alpha)
        assert sum == pha_poly_at_alpha, "should equal for sum == pha()"
        R_poly = d_t_sum_poly - pha_poly_at_alpha

        # Compute commitment:
        # Q_D(X) = (D(X) * t_I(X) - φ(α) - R(X)) / z_I(X)
        Q_D_poly = (D_t_poly - pha_poly_at_alpha - R_poly) / z_I_poly
        # Q_E(X) = (E(X) * (β - v(X)) + v(X) * z_I(β) / z_I(0)) / z_V(X)
        z_I_at_beta = z_I_poly.coeff_eval(Scalar(beta))
        Q_E_poly = (E_poly * (v_poly * Scalar(-1) + beta) +
                    v_poly * z_I_at_beta / z_I_at_0) / z_V_poly

        self.z_V_poly = z_V_poly
        self.D_poly = D_poly
        self.E_poly = E_poly
        self.R_poly = R_poly
        self.Q_D_poly = Q_D_poly
        self.Q_E_poly = Q_E_poly

        # π2 = ([D]1 = [D(x)]1, [R]1 = [R(x)]1, [Q2]1 = [Q2(x)]1)
        # π3 = ([E]1 = [E(x)]1, [Q1]1 = [Q1(x)]1)
        R_comm_1 = setup.commit_g1(R_poly)
        D_comm_1 = setup.commit_g1(D_poly)
        Q_D_comm_1 = setup.commit_g1(Q_D_poly)
        E_comm_1 = setup.commit_g1(E_poly)
        Q_E_comm_1 = setup.commit_g1(Q_E_poly)

        return Message2(
            D_comm_1,
            R_comm_1,
            Q_D_comm_1,
            E_comm_1,
            Q_E_comm_1
        )

    def round_3(self) -> Message3:
        setup = self.setup
        alpha = self.alpha
        beta = self.beta
        gamma = self.gamma
        zeta = self.zeta
        phi_poly = self.phi_poly
        z_I_poly = self.z_I_poly
        t_I_poly = self.t_I_poly
        v_poly = self.v_poly
        z_V_poly = self.z_V_poly
        D_poly = self.D_poly
        E_poly = self.E_poly
        R_poly = self.R_poly
        Q_D_poly = self.Q_D_poly
        Q_E_poly = self.Q_E_poly
        t_poly = self.t_poly
        z_H_poly = self.z_H_poly
        d = len(setup.powers_of_x)
        m = self.m

        # calculate v1, v2, v3, v4, v5
        # v1 = e(α)
        v1 = E_poly.coeff_eval(alpha)
        # v2 = a(α)
        v2 = phi_poly.coeff_eval(alpha)
        # v3 = z_I(0)
        v3 = z_I_poly.coeff_eval(Scalar(0))
        # v4 = z_I(β)
        v4 = z_I_poly.coeff_eval(beta)
        # v5 = e(ζ)
        v5 = E_poly.coeff_eval(zeta)
        # calculate P_D(X), P_E(X)
        # P_D(X) = D(β) * t_I(X) - φ(α) - R(X) - z_I(β) * Q_D(X)
        # P_E(X) = E(ζ) * (β - v(X)) + v(X) * z_I(β) / z_I(0) - z_V(ζ) * Q_E(X)
        D_poly_at_beta = D_poly.coeff_eval(beta)
        P_D_poly = t_I_poly * D_poly_at_beta - v2 - R_poly - Q_D_poly * v4
        E_poly_at_zeta = E_poly.coeff_eval(zeta)
        z_V_poly_at_zeta = z_V_poly.coeff_eval(zeta)
        P_E_poly = (v_poly * Scalar(-1) + beta) * E_poly_at_zeta + \
            v_poly * v4 / v3 - Q_E_poly * z_V_poly_at_zeta

        # X^(d-m+1)
        x_exponent_poly = poly_util.x_exponent_poly(d - m + 1)
        # calculate [w1]1, [w2]1, [w2]1, [w4]1
        # X - α
        x_alpha_poly = poly_util.root_poly(Scalar(alpha))
        # calculate w1
        w1_poly = x_exponent_poly * \
            (E_poly - v1 + (phi_poly - v2) * gamma) / x_alpha_poly
        x_poly = poly_util.root_poly(Scalar(0))
        # x^m
        x_m_exponent_poly = poly_util.x_exponent_poly(m)
        # calculate w2
        w2_poly = (z_I_poly - v3) / x_poly + R_poly * gamma / x_poly + \
            x_exponent_poly * (
                (z_I_poly - x_m_exponent_poly) * gamma ** 2 +
                R_poly * gamma ** 3
            )
        # X - β
        x_beta_poly = poly_util.root_poly(Scalar(beta))
        # calculate w3
        w3_poly = (D_poly - v1) / x_beta_poly + \
            (z_I_poly - v4) * gamma / x_beta_poly + \
            P_D_poly * gamma ** 2 / x_beta_poly
        # X - ζ
        x_zeta_poly = poly_util.root_poly(Scalar(zeta))
        # calculate w4
        w4_poly = (E_poly - v5) / x_zeta_poly + \
            P_E_poly * gamma / x_zeta_poly

        w1_comm_1 = setup.commit_g1(w1_poly)
        w2_comm_1 = setup.commit_g1(w2_poly)
        w3_comm_1 = setup.commit_g1(w3_poly)
        w4_comm_1 = setup.commit_g1(w4_poly)

        # caulk+
        w5_poly = (t_poly - t_I_poly) / z_I_poly
        w6_poly = z_H_poly / z_I_poly
        w5_comm_1 = setup.commit_g1(w5_poly)
        w6_comm_1 = setup.commit_g1(w6_poly)

        return Message3(
            v1,
            v2,
            v3,
            v4,
            v5,
            w1_comm_1,
            w2_comm_1,
            w3_comm_1,
            w4_comm_1,
            w5_comm_1,
            w6_comm_1
        )
