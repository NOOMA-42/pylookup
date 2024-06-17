from dataclasses import dataclass
from src.common_util.curve import Scalar
from src.common_util.mle_poly import polynomial, get_multi_poly_lagrange, eq_mle_poly
from src.common_util.sumcheck import prove_sumcheck
from src.lasso.program import Params, SOSTable, GrandProductData, log_ceil, Hash
from src.lasso.setup import Setup
from src.lasso.transcript import Transcript, Message1, Message2, Message3, Message4, Message5

@dataclass
class Proof:
    msg_1: Message1
    msg_2: Message2
    msg_3: Message3
    msg_4: Message4
    msg_5: Message5

    def flatten(self):
        proof = {}
        # msg_1
        proof["a_comm"] = self.msg_1.a_comm
        proof["logm"] = self.msg_1.logm
        proof["dim_comm"] = self.msg_1.dim_comm
        # msg_2
        proof["a_eval"] = self.msg_2.a_eval
        proof["a_eval_proof"] = self.msg_2.a_eval_proof
        proof["E_comm"] = self.msg_2.E_comm
        proof["read_ts_comm"] = self.msg_2.read_ts_comm
        proof["final_cts_comm"] = self.msg_2.final_cts_comm
        # msg_3
        proof["h_sumcheck_proof"] = self.msg_3.h_sumcheck_proof
        proof["rz"] = self.msg_3.rz
        proof["E_eval"] = self.msg_3.E_eval
        proof["E_eval_proof"] = self.msg_3.E_eval_proof
        # msg_4
        proof["S_comm"] = self.msg_4.S_comm
        proof["RS_comm"] = self.msg_4.RS_comm
        proof["WS1_comm"] = self.msg_4.WS1_comm
        proof["WS2_comm"] = self.msg_4.WS2_comm
        # msg_5
        proof["S_sumcheck_proof"] = self.msg_5.S_sumcheck_proof
        proof["RS_sumcheck_proof"] = self.msg_5.RS_sumcheck_proof
        proof["WS1_sumcheck_proof"] = self.msg_5.WS1_sumcheck_proof
        proof["WS2_sumcheck_proof"] = self.msg_5.WS2_sumcheck_proof
        proof["r_prime2"] = self.msg_5.r_prime2
        proof["r_prime3"] = self.msg_5.r_prime3
        proof["r_prime4"] = self.msg_5.r_prime4
        proof["r_prime5"] = self.msg_5.r_prime5
        proof["S_data"] = self.msg_5.S_data
        proof["RS_data"] = self.msg_5.RS_data
        proof["WS1_data"] = self.msg_5.WS1_data
        proof["WS2_data"] = self.msg_5.WS2_data
        proof["E_eval2"] = self.msg_5.E_eval2
        proof["dim_eval"] = self.msg_5.dim_eval
        proof["read_ts_eval"] = self.msg_5.read_ts_eval
        proof["final_cts_eval"] = self.msg_5.final_cts_eval
        proof["E_eval2_proof"] = self.msg_5.E_eval2_proof
        proof["dim_eval_proof"] = self.msg_5.dim_eval_proof
        proof["read_ts_eval_proof"] = self.msg_5.read_ts_eval_proof
        proof["final_cts_eval_proof"] = self.msg_5.final_cts_eval_proof
        return proof

@dataclass
class Prover:
    # the notations follow the lasso paper (https://eprint.iacr.org/2023/1216.pdf)
    table: SOSTable
    l: int
    c: int
    k: int
    alpha: int
    m: int
    logm: int

    def __init__(self, setup: Setup, params: Params):
        self.setup = setup
        self.table = params.table
        self.l = self.table.l
        self.c = self.table.c
        self.k = self.table.k
        self.alpha = self.table.alpha

    def prove(self, witness: list[int]):
        transcript = Transcript(b"lasso")
        # Round 1
        msg_1 = self.round_1(witness)

        self.r = transcript.round_1(msg_1)
        # Round 2
        msg_2 = self.round_2()
        
        transcript.round_2(msg_2)
        # Round 3
        msg_3 = self.round_3()

        self.tau, self.gamma = transcript.round_3(msg_3)
        # Round 4
        msg_4 = self.round_4()

        transcript.round_4(msg_4)
        # Round 5
        msg_5 = self.round_5()

        return Proof(msg_1, msg_2, msg_3, msg_4, msg_5)

    def round_1(self, witness: list[int]) -> Message1:
        self.m = len(witness)
        self.logm = log_ceil(self.m)
        assert(self.logm <= self.l)

        self.a = [Scalar(0) for _ in range(2**self.logm)]
        for i in range(self.m):
            self.a[i] = Scalar(witness[i])
        self.a_poly = get_multi_poly_lagrange(self.a, self.logm)
        self.a_comm = self.setup.commit_g1(self.a_poly, self.logm)
        self.indexes = []
        for w in witness:
            self.indexes.append(self.table.get_index(w))

        self.dim_poly = []
        for i in range(self.c):
            values = [0 for _ in range(2**self.logm)]
            for j in range(self.m):
                values[j] = self.indexes[j][i]
            self.dim_poly.append(get_multi_poly_lagrange(list(map(Scalar, values)), self.logm))

        self.dim_comm = [self.setup.commit_g1(poly, self.logm) for poly in self.dim_poly]
        return Message1(self.a_comm, self.logm, self.dim_comm)
    
    def round_2(self):
        self.a_eval = self.a_poly.eval(self.r)
        self.a_eval_proof = self.setup.prove(self.a_poly, self.r, self.a_eval)

        self.E_poly, self.read_poly, self.write_poly, self.final_poly = [], [], [], []
        for i in range(self.alpha):
            # Offline memory checking, or "Memory in the head"
            # See Spartan(https://eprint.iacr.org/2019/550.pdf), Section 7.2.1
            E_val = [0 for _ in range(2**self.logm)]
            values = [0 for _ in range(2**self.logm)]
            read_ts = [0 for _ in range(2**self.logm)]
            write_cts = [0 for _ in range(2**self.logm)]
            final_cts = [0 for _ in range(2**self.l)]
            ts = 0
            for j in range(self.m):
                values[j] = self.indexes[j][i//self.k]
                E_val[j] = self.table.tables[i][values[j]]
                ts = final_cts[values[j]]
                read_ts[j] = ts
                write_cts[j] = ts+1
                final_cts[values[j]] = ts+1

            self.E_poly.append(get_multi_poly_lagrange(list(map(Scalar, E_val)), self.logm))
            self.read_poly.append(get_multi_poly_lagrange(list(map(Scalar, read_ts)), self.logm))
            self.write_poly.append(get_multi_poly_lagrange(list(map(Scalar, write_cts)), self.logm))
            self.final_poly.append(get_multi_poly_lagrange(list(map(Scalar, final_cts)), self.l))
        
        self.E_comm = [self.setup.commit_g1(poly, self.logm) for poly in self.E_poly]
        self.read_comm = [self.setup.commit_g1(poly, self.logm) for poly in self.read_poly]
        self.final_comm = [self.setup.commit_g1(poly, self.l) for poly in self.final_poly]
        return Message2(self.a_eval, self.a_eval_proof, self.E_comm, self.read_comm, self.final_comm)
    
    def round_3(self):
        # sumcheck protocol on h(k) = eq(r, k) * g(...E_i(k))
        h_poly = eq_mle_poly(self.r) * self.table.g_func(self.E_poly)
        self.h_sumcheck_proof, self.rz = prove_sumcheck(h_poly, self.logm, 1)
        self.E_eval = [E.eval(self.rz) for E in self.E_poly]
        self.E_eval_proof = [self.setup.prove(e_poly, self.rz, e_eval) for (e_poly, e_eval) in zip(self.E_poly, self.E_eval)]
        return Message3(self.h_sumcheck_proof, self.rz, self.E_eval, self.E_eval_proof)
    
    def round_4(self):
        self.S_poly, self.RS_poly, self.WS1_poly, self.WS2_poly = [], [], [], []
        self.S_comm, self.RS_comm, self.WS1_comm, self.WS2_comm = [], [], [], []
        for i in range(self.alpha):
            S, RS, WS1, WS2 = [], [], [], []
            for j in range(2**self.logm):
                RS.append((self.dim_poly[i].values[j], self.E_poly[i].values[j], self.read_poly[i].values[j]))
                WS1.append((self.dim_poly[i].values[j], self.E_poly[i].values[j], self.write_poly[i].values[j]))
            for j in range(2**self.l):
                S.append((Scalar(j), Scalar(self.table.tables[i][j]), self.final_poly[i].values[j]))
                WS2.append((Scalar(j), Scalar(self.table.tables[i][j]), Scalar(0)))

            S_poly = self.grand_product_poly(S, self.l)
            RS_poly = self.grand_product_poly(RS, self.logm)
            WS1_poly = self.grand_product_poly(WS1, self.logm)
            WS2_poly = self.grand_product_poly(WS2, self.l)
            self.S_poly.append(S_poly)
            self.RS_poly.append(RS_poly)
            self.WS1_poly.append(WS1_poly)
            self.WS2_poly.append(WS2_poly)
            self.S_comm.append(self.setup.commit_g1(S_poly, self.l+1))
            self.RS_comm.append(self.setup.commit_g1(RS_poly, self.logm+1))
            self.WS1_comm.append(self.setup.commit_g1(WS1_poly, self.logm+1))
            self.WS2_comm.append(self.setup.commit_g1(WS2_poly, self.l+1))

        return Message4(self.S_comm, self.RS_comm, self.WS1_comm, self.WS2_comm)
    
    def round_5(self):
        self.S_sumcheck_proof, self.RS_sumcheck_proof, self.WS1_sumcheck_proof, self.WS2_sumcheck_proof = [], [], [], []
        self.r_prime2, self.r_prime3, self.r_prime4, self.r_prime5 = [], [], [], []
        self.S_data, self.RS_data, self.WS1_data, self.WS2_data = [], [], [], []
        self.E_eval2, self.dim_eval, self.read_eval, self.final_eval = [], [], [], []
        self.E_eval2_proof, self.dim_eval_proof, self.read_eval_proof, self.final_eval_proof = [], [], [], []
        for i in range(self.alpha):
            self.handle_grand_product_sumcheck(self.S_sumcheck_proof, self.r_prime2, self.S_poly[i], self.l)
            self.handle_grand_product_sumcheck(self.RS_sumcheck_proof, self.r_prime3, self.RS_poly[i], self.logm)
            self.handle_grand_product_sumcheck(self.WS1_sumcheck_proof, self.r_prime4, self.WS1_poly[i], self.logm)
            self.handle_grand_product_sumcheck(self.WS2_sumcheck_proof, self.r_prime5, self.WS2_poly[i], self.l)
            self.S_data.append(self.generate_grand_product_data(self.S_poly[i], self.r_prime2[i]))
            self.RS_data.append(self.generate_grand_product_data(self.RS_poly[i], self.r_prime3[i]))
            self.WS1_data.append(self.generate_grand_product_data(self.WS1_poly[i], self.r_prime4[i]))
            self.WS2_data.append(self.generate_grand_product_data(self.WS2_poly[i], self.r_prime5[i]))
            self.E_eval2.append(self.E_poly[i].eval(self.r_prime3[i]))
            self.dim_eval.append(self.dim_poly[i//self.k].eval(self.r_prime3[i]))
            self.read_eval.append(self.read_poly[i].eval(self.r_prime3[i]))
            self.final_eval.append(self.final_poly[i].eval(self.r_prime2[i]))
            self.E_eval2_proof.append(self.setup.prove(self.E_poly[i], self.r_prime3[i], self.E_eval2[i]))
            self.dim_eval_proof.append(self.setup.prove(self.dim_poly[i//self.k], self.r_prime3[i], self.dim_eval[i]))
            self.read_eval_proof.append(self.setup.prove(self.read_poly[i], self.r_prime3[i], self.read_eval[i]))
            self.final_eval_proof.append(self.setup.prove(self.final_poly[i], self.r_prime2[i], self.final_eval[i]))

        return Message5(self.S_sumcheck_proof, self.RS_sumcheck_proof, self.WS1_sumcheck_proof, self.WS2_sumcheck_proof,
                        self.r_prime2, self.r_prime3, self.r_prime4, self.r_prime5,
                        self.S_data, self.RS_data, self.WS1_data, self.WS2_data,
                        self.E_eval2, self.dim_eval, self.read_eval, self.final_eval,
                        self.E_eval2_proof, self.dim_eval_proof, self.read_eval_proof, self.final_eval_proof)
    
    def grand_product_poly(self, multiset: list[tuple[Scalar, Scalar, Scalar]], length: int):
        # see https://eprint.iacr.org/2020/1275.pdf, section 5
        f = [Hash(s, self.tau, self.gamma) for s in multiset] # f(0,x) = v(x)
        for i in range(len(f)-1):
            f.append(f[2*i]*f[2*i+1]) # f(1,x) = f(x,0) * f(x,1)
        f.append(Scalar(0))
        poly_f = get_multi_poly_lagrange(f, length+1)
        return poly_f
    
    def grand_product_sumcheck(self, poly: polynomial, length: int):
        # run sumcheck protocol on g(k) = eq(r,k) * (f(1,k)-f(k,0)*f(k,1))
        diff_poly = poly.eval_i(Scalar(1), 1) - poly.eval_i(Scalar(0), length+1) * poly.eval_i(Scalar(1), length+1)
        g_poly = eq_mle_poly(self.r) * diff_poly
        return prove_sumcheck(g_poly, length, 1)
    
    def handle_grand_product_sumcheck(self, data_list, r_list, poly: polynomial, length: int):
        data, r = self.grand_product_sumcheck(poly, length)
        data_list.append(data)
        r_list.append(r)

    def generate_grand_product_data(self, f: polynomial, r: list[Scalar]):
        base = f.eval([Scalar(0)]+r)
        point = [Scalar(1)]*len(r)+[Scalar(0)]
        product = f.eval(point)
        base_proof = self.setup.prove(f, [Scalar(0)]+r, base)
        product_proof = self.setup.prove(f, point, product)
        return GrandProductData(base, product, base_proof, product_proof)
