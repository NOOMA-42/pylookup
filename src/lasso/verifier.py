import py_ecc.bn128 as b
from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point
from src.common_util.mle_poly import get_multi_poly_lagrange
from src.common_util.sumcheck import verify_sumcheck
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
        transcript.round_2(pf.msg_2)
        tau, gamma = transcript.round_3(pf.msg_3)
        transcript.round_4(pf.msg_4)

        # get proofs
        proof = pf.flatten()
        a_comm = proof["a_comm"]
        logm = proof["logm"]
        dim_comm = proof["dim_comm"]
        a_eval = proof["a_eval"]
        a_eval_proof = proof["a_eval_proof"]
        E_comm = proof["E_comm"]
        read_ts_comm = proof["read_ts_comm"]
        final_cts_comm = proof["final_cts_comm"]
        h_sumcheck_proof = proof["h_sumcheck_proof"]
        rz = proof["rz"]
        E_eval = proof["E_eval"]
        E_eval_proof = proof["E_eval_proof"]
        S_comm = proof["S_comm"]
        RS_comm = proof["RS_comm"]
        WS1_comm = proof["WS1_comm"]
        WS2_comm = proof["WS2_comm"]
        S_sumcheck_proof = proof["S_sumcheck_proof"]
        RS_sumcheck_proof = proof["RS_sumcheck_proof"]
        WS1_sumcheck_proof = proof["WS1_sumcheck_proof"]
        WS2_sumcheck_proof = proof["WS2_sumcheck_proof"]
        r_prime2 = proof["r_prime2"]
        r_prime3 = proof["r_prime3"]
        r_prime4 = proof["r_prime4"]
        r_prime5 = proof["r_prime5"]
        S_data = proof["S_data"]
        RS_data = proof["RS_data"]
        WS1_data = proof["WS1_data"]
        WS2_data = proof["WS2_data"]
        E_eval2 = proof["E_eval2"]
        dim_eval = proof["dim_eval"]
        read_ts_eval = proof["read_ts_eval"]
        final_cts_eval = proof["final_cts_eval"]
        E_eval2_proof = proof["E_eval2_proof"]
        dim_eval_proof = proof["dim_eval_proof"]
        read_ts_eval_proof = proof["read_ts_eval_proof"]
        final_cts_eval_proof = proof["final_cts_eval_proof"]

        print("=== Started Check 1: check value of a(r) ===")
        assert(self.setup.verify(a_comm, r, a_eval, a_eval_proof))
        print("=== Finished Check 1: check value of a(r) ===")
        
        print("=== Started Check 2: sum check protocol of h ===")
        assert(verify_sumcheck(a_eval, h_sumcheck_proof, rz, logm))
        # Todo: assert(h(rz) == chi(r, rz) * self.table.g_func(E_eval))
        print("=== Finished Check 2: sum check protocol of h ===")

        print("=== Started Check 3: check values of E(rz) ===")
        for eval, comm, proof in zip(E_eval, E_comm, E_eval_proof):
            assert(self.setup.verify(comm, rz, eval, proof))
        print("=== Finished Check 3: check values of E(rz) ===")
        
        print("=== Started Check 4: sum check protocol of grand product ===")
        for i in range(self.alpha):
            verify_sumcheck(Scalar(0), S_sumcheck_proof[i], r_prime2[i], self.l)
            verify_sumcheck(Scalar(0), RS_sumcheck_proof[i], r_prime3[i], logm)
            verify_sumcheck(Scalar(0), WS1_sumcheck_proof[i], r_prime4[i], logm)
            verify_sumcheck(Scalar(0), WS2_sumcheck_proof[i], r_prime5[i], self.l)
        print("=== Finished Check 4: sum check protocol of grand product ===")

        print("=== Started Check 5: check value of E, dim, read_ts, final_cts ===")
        for i in range(self.alpha):
            self.verify_grand_product(S_data[i], S_comm[i], r_prime2[i])
            self.verify_grand_product(RS_data[i], RS_comm[i], r_prime3[i])
            self.verify_grand_product(WS1_data[i], WS1_comm[i], r_prime4[i])
            self.verify_grand_product(WS2_data[i], WS2_comm[i], r_prime5[i])
            assert(S_data[i].product * RS_data[i].product == WS1_data[i].product * WS2_data[i].product)
            assert(self.setup.verify(E_comm[i], r_prime3[i], E_eval2[i], E_eval2_proof[i]))
            assert(self.setup.verify(dim_comm[i//self.k], r_prime3[i], dim_eval[i], dim_eval_proof[i]))
            assert(self.setup.verify(read_ts_comm[i], r_prime3[i], read_ts_eval[i], read_ts_eval_proof[i]))
            assert(self.setup.verify(final_cts_comm[i], r_prime2[i], final_cts_eval[i], final_cts_eval_proof[i]))
            assert(RS_data[i].base == Hash((dim_eval[i], E_eval2[i], read_ts_eval[i]), tau, gamma))
            identity_poly = get_multi_poly_lagrange([Scalar(i) for i in range(2**self.l)], self.l)
            identity_eval = identity_poly.eval(r_prime2[i])
            table_poly = get_multi_poly_lagrange(list(map(Scalar, self.table.tables[i])), self.l)
            table_eval = table_poly.eval(r_prime2[i])
            assert(S_data[i].base == Hash((identity_eval, table_eval, final_cts_eval[i]), tau, gamma))
        print("=== Finished Check 5: check value of E, dim, read_ts, final_cts ===")

        print("Finished to verify proof")
        return True

    def verify_grand_product(self, data: GrandProductData, comm: G1Point, point: list[Scalar]):
        assert(self.setup.verify(comm, [Scalar(0)]+point, data.base, data.base_proof))
        assert(self.setup.verify(comm, [Scalar(1)]*len(point)+[Scalar(0)], data.product, data.product_proof))
