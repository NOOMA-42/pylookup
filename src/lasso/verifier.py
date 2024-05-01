import py_ecc.bn128 as b
from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point, ec_lincomb
from src.common_util.poly import Polynomial, Basis
from src.lasso.program import Params, SOSTable, GrandProductData, Hash
from src.lasso.prover import Proof
from src.lasso.setup import Setup
from src.lasso.transcript import Transcript


@dataclass
class Verifier:
    # the notations follow the lasso paper (https://eprint.iacr.org/2023/1216.pdf)
    table: SOSTable
    l: int
    c: int
    k: int
    alpha: int

    def __init__(self, setup: Setup, params: Params):
        self.setup = setup
        self.table = params.table
        self.l = self.table.l
        self.c = self.table.c
        self.k = self.table.k
        self.alpha = self.table.alpha
    
    def verify(self, pf: Proof) -> bool:
        print("Start to verify proof")
        transcript = Transcript(b"lasso")
        r = transcript.round_1(pf.msg_1)
        rz = transcript.round_2(pf.msg_2)
        tau, gamma, r_prime2, r_prime3 = transcript.round_3(pf.msg_3)
        r_prime2, r_prime3, r_prime4, r_prime5 = transcript.round_4(pf.msg_4)

        # get proofs
        proof = pf.flatten()
        a_comm = proof["a_comm"]
        logm = proof["logm"]
        dim_comm = proof["dim_comm"]
        a_eval = proof["a_eval"]
        a_PIOP = proof["a_PIOP"]
        E_comm = proof["E_comm"]
        read_ts_comm = proof["read_ts_comm"]
        final_cts_comm = proof["final_cts_comm"]
        sumcheck_h_data = proof["sumcheck_h_data"]
        E_eval = proof["E_eval"]
        E_PIOP = proof["E_PIOP"]
        S_comm = proof["S_comm"]
        RS_comm = proof["RS_comm"]
        WS1_comm = proof["WS1_comm"]
        WS2_comm = proof["WS2_comm"]
        sumcheck_S_data = proof["sumcheck_S_data"]
        sumcheck_RS_data = proof["sumcheck_RS_data"]
        sumcheck_WS1_data = proof["sumcheck_WS1_data"]
        sumcheck_WS2_data = proof["sumcheck_WS2_data"]
        S_data = proof["S_data"]
        RS_data = proof["RS_data"]
        WS1_data = proof["WS1_data"]
        WS2_data = proof["WS2_data"]
        E_eval2 = proof["E_eval2"]
        dim_eval = proof["dim_eval"]
        read_ts_eval = proof["read_ts_eval"]
        final_cts_eval = proof["final_cts_eval"]
        E_PIOP2 = proof["E_PIOP2"]
        dim_PIOP = proof["dim_PIOP"]
        read_ts_PIOP = proof["read_ts_PIOP"]
        final_cts_PIOP = proof["final_cts_PIOP"]

        print("=== Started Check 1: check value of a(r) ===")
        assert(self.setup.PIOP_verify(a_comm, r, a_eval, a_PIOP))
        print("=== Finished Check 1: check value of a(r) ===")
        
        print("=== Started Check 1: sum check protocol of h ===")
        h_eval = self.verify_sumcheck_data(sumcheck_h_data, a_eval)
        assert(h_eval == self.setup.eq_mle(r, rz) * self.table.g_func(E_eval))
        print("=== Finished Check 1: sum check protocol of h ===")

        print("=== Started Check 2: check values of E(rz) ===")
        for eval, comm, PIOP in range(zip(E_eval, E_comm, E_PIOP)):
            assert(self.setup.PIOP_verify(comm, rz, eval, PIOP))
        print("=== Finished Check 2: check values of E(rz) ===")
        
        print("=== Started Check 3: sum check protocol of grand product ===")
        for i in range(self.alpha):
            self.verify_sumcheck_data(sumcheck_S_data[i], Scalar(0))
            self.verify_sumcheck_data(sumcheck_RS_data[i], Scalar(0))
            self.verify_sumcheck_data(sumcheck_WS1_data[i], Scalar(0))
            self.verify_sumcheck_data(sumcheck_WS2_data[i], Scalar(0))
        print("=== Finished Check 3: sum check protocol of grand product ===")

        print("=== Started Check 4: check value of E, dim, read_ts, final_cts ===")
        for i in range(self.alpha):
            self.verify_grand_product_PIOP(S_data[i], S_comm[i], r_prime2[i])
            self.verify_grand_product_PIOP(RS_data[i], RS_comm[i], r_prime3[i])
            self.verify_grand_product_PIOP(WS1_data[i], WS1_comm[i], r_prime4[i])
            self.verify_grand_product_PIOP(WS2_data[i], WS2_comm[i], r_prime5[i])
            assert(S_data[i].product * RS_data[i].product == WS1_data[i].product * WS2_data[i].product)
            assert(self.setup.PIOP_verify(E_comm[i], r_prime3, E_eval2[i], E_PIOP2[i]))
            assert(self.setup.PIOP_verify(dim_comm[i//self.k], r_prime3, dim_eval, dim_PIOP))
            assert(self.setup.PIOP_verify(read_ts_comm, r_prime3, read_ts_eval, read_ts_PIOP))
            assert(self.setup.PIOP_verify(final_cts_comm, r_prime2, final_cts_eval, final_cts_PIOP))
            assert(RS_data[i].f_0_r == Hash((dim_eval[i], E_eval2[i], read_ts_eval[i]), tau, gamma))
            identity_poly = Polynomial([Scalar(i) for i in range(2**self.l)], Basis.LAGRANGE)
            identity_eval = self.setup.multivar_eval(identity_poly, r_prime2[i])
            table_poly = Polynomial(map(self.table.tables[i], Scalar), Basis.LAGRANGE)
            table_eval = self.setup.multivar_eval(table_poly, r_prime2[i])
            assert(S_data[i].f_0_r == Hash((identity_eval, table_eval, final_cts_eval[i]), tau, gamma))
        print("=== Finished Check 4: check value of E, dim, read_ts, final_cts ===")

        print("Finished to verify proof")
        return True

    def verify_sumcheck_data(self, data: list, sum: Scalar):
        now = sum
        # For each round, get evaluation on [0,1,2] for a degree-2 polynomial p
        # assert (p(0)+p(1) == now)
        # get the next bit of r (Fiat-Shamir?)
        # recover the degree-2 polynomial p and set now=p(r)
        return now

    def verify_grand_product_PIOP(self, data: GrandProductData, comm: G1Point, point: list):
        assert(self.setup.PIOP_verify(comm, point, data.f_0_r, data.f_0_r_PIOP))
        assert(self.setup.PIOP_verify(comm, point, data.f_1_r, data.f_1_r_PIOP))
        assert(self.setup.PIOP_verify(comm, point, data.f_r_0, data.f_r_0_PIOP))
        assert(self.setup.PIOP_verify(comm, point, data.f_r_1, data.f_r_1_PIOP))
        assert(self.setup.PIOP_verify(comm, point, data.product, data.product_PIOP))