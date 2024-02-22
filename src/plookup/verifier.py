import py_ecc.bn128 as b
from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point, ec_lincomb
from src.common_util.poly import Polynomial, Basis
from src.plookup.program import Params, aggregate, aggregate_comm
from src.plookup.prover import Proof
from src.plookup.setup import Setup
from src.plookup.transcript import Transcript


@dataclass
class Verifier:
    # the notations follow the plookup paper (https://eprint.iacr.org/2020/315.pdf)
    table: list[int]
    n: int
    H: list[Scalar] # [1, ..., g^n]

    def __init__(self, setup: Setup, params: Params):
        self.setup = setup
        self.table = params.table
        self.n = params.order - 1
        self.H = params.roots
        # vanish = x^(n+1) - 1
        self.vanish = Polynomial([Scalar(-1)]+[Scalar(0)]*self.n+[Scalar(1)], Basis.MONOMIAL)
    
    def verify(self, pf: Proof) -> bool:
        print("Start to verify proof")
        transcript = Transcript(b"plonk_plookup")
        beta, gamma = transcript.round_1(pf.msg_1)
        delta = transcript.round_2(pf.msg_2)
        zeta = transcript.round_3(pf.msg_3, self.n+1)
        eps = transcript.round_4(pf.msg_4)

        # get proofs
        proof = pf.flatten()
        f_comm = proof["f_comm"]
        h1_comm = proof["h1_comm"]
        h2_comm = proof["h2_comm"]
        z_comm = proof["z_comm"]
        q_comm = proof["q_comm"]
        f_eval = proof["f_eval"]
        h1_eval = proof["h1_eval"]
        h2_eval = proof["h2_eval"]
        z_eval = proof["z_eval"]
        h1_g_eval = proof["h1_g_eval"]
        h2_g_eval = proof["h2_g_eval"]
        z_g_eval = proof["z_g_eval"]
        agg_witness_comm = proof["agg_witness_comm"]
        agg_g_witness_comm = proof["agg_g_witness_comm"]

        # verify
        L0 = Polynomial([Scalar(1) if i == 0 else Scalar(0) for i in range(self.n+1)], Basis.LAGRANGE).ifft()
        Ln = Polynomial([Scalar(1) if i == self.n else Scalar(0) for i in range(self.n+1)], Basis.LAGRANGE).ifft()
        L0_zeta = L0.eval(zeta)
        Ln_zeta = Ln.eval(zeta)
        g = self.H[1]

        def get_b():
            t_poly = Polynomial(list(map(Scalar, self.table)), Basis.LAGRANGE)
            t_eval = t_poly.eval(zeta)
            t_g_eval = t_poly.eval(g*zeta)
            lhs = (zeta - self.H[self.n]) * z_eval * (Scalar(1) + beta) * \
                (gamma + f_eval) * (gamma * (Scalar(1) + beta) + t_eval + beta * t_g_eval)
            rhs = (zeta - self.H[self.n]) * z_g_eval * (gamma * (1+beta) + h1_eval + beta * h1_g_eval) * \
                (gamma * (Scalar(1) + beta) + h2_eval + beta * h2_g_eval)
            return lhs - rhs
        
        a = L0_zeta * (z_eval - Scalar(1))
        b = get_b()
        c = Ln_zeta * (h1_eval - h2_g_eval)
        d = Ln_zeta * (z_eval - Scalar(1))
        q_eval = aggregate(delta, [a, b, c, d]) / self.vanish.eval(zeta)
        
        print("=== Started Check 1: evaluation on zeta ===")
        # agg(x) = aggregate(eps, [f, h1, h2, z, q])
        agg_comm = aggregate_comm(eps, [f_comm, h1_comm, h2_comm, z_comm, q_comm]) # agg(tau)
        agg_eval = aggregate(eps, [f_eval, h1_eval, h2_eval, z_eval, q_eval])      # agg(zeta)
        self.check_pairing(agg_comm, agg_eval, agg_witness_comm, zeta)
        print("=== Finished Check 1: evaluation on zeta ===")

        print("=== Started Check 2: evaluation on g*zeta ===")
        # agg_g(x) = aggregate(eps, [h1, h2, z])
        agg_g_comm = aggregate_comm(eps, [h1_comm, h2_comm, z_comm])   # agg_g(tau)
        agg_g_eval = aggregate(eps, [h1_g_eval, h2_g_eval, z_g_eval])  # agg_g(g*zeta)
        self.check_pairing(agg_g_comm, agg_g_eval, agg_g_witness_comm, g*zeta)
        print("=== Finished Check 2: evaluation on g*zeta ===")

        print("Finished to verify proof")
        return True
    
    def check_pairing(self, comm: G1Point, eval: Scalar, witness_comm: G1Point, point: Scalar):
        # comm: P(tau)
        # eval: P(point)
        # witness_com: Q(tau), Q(x) = P(x)-P(point)/(x-point)
        # e(P(tau)-P(point), 1) = e(Q(tau), tau-point)
        comb_lhs = ec_lincomb([(comm, Scalar(1)), (b.G1, -eval)])
        comb_rhs = self.setup.commit_g2(Polynomial([-point, Scalar(1)], Basis.MONOMIAL))
        pairing_lhs = b.pairing(b.G2, comb_lhs)
        pairing_rhs = b.pairing(comb_rhs, witness_comm)
        assert(pairing_lhs == pairing_rhs)
