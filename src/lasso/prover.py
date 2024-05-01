from dataclasses import dataclass
from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar
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
        proof["a_PIOP"] = self.msg_2.a_PIOP
        proof["E_comm"] = self.msg_2.E_comm
        proof["read_ts_comm"] = self.msg_2.read_ts_comm
        proof["final_cts_comm"] = self.msg_2.final_cts_comm
        # msg_3
        proof["sumcheck_h_data"] = self.msg_3.sumcheck_h_data
        proof["E_eval"] = self.msg_3.E_eval
        proof["E_PIOP"] = self.msg_3.E_PIOP
        # msg_4
        proof["S_comm"] = self.msg_4.S_comm
        proof["RS_comm"] = self.msg_4.RS_comm
        proof["WS1_comm"] = self.msg_4.WS1_comm
        proof["WS2_comm"] = self.msg_4.WS2_comm
        # msg_5
        proof["sumcheck_S_data"] = self.msg_5.sumcheck_S_data
        proof["sumcheck_RS_data"] = self.msg_5.sumcheck_RS_data
        proof["sumcheck_WS1_data"] = self.msg_5.sumcheck_WS1_data
        proof["sumcheck_WS2_data"] = self.msg_5.sumcheck_WS2_data
        proof["S_data"] = self.msg_5.S_data
        proof["RS_data"] = self.msg_5.RS_data
        proof["WS1_data"] = self.msg_5.WS1_data
        proof["WS2_data"] = self.msg_5.WS2_data
        proof["E_eval2"] = self.msg_5.E_eval2
        proof["dim_eval"] = self.msg_5.dim_eval
        proof["read_ts_eval"] = self.msg_5.read_ts_eval
        proof["final_cts_eval"] = self.msg_5.final_cts_eval
        proof["E_PIOP2"] = self.msg_5.E_PIOP2
        proof["dim_PIOP"] = self.msg_5.dim_PIOP
        proof["read_ts_PIOP"] = self.msg_5.read_ts_PIOP
        proof["final_cts_PIOP"] = self.msg_5.final_cts_PIOP
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
        
        # Todo: rz should be gotten from sumcheck protocol
        # To be non-iteractive, we can get a seed here and use Fiat-Shamir to generate the whole rz
        self.rz = transcript.round_2(msg_2)
        # Round 3
        msg_3 = self.round_3()

        self.tau, self.gamma = transcript.round_3(msg_3)
        # Round 4
        msg_4 = self.round_4()

        self.r_prime2, self.r_prime3, self.r_prime4, self.r_prime5 = transcript.round_4(msg_4)
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
        self.a_poly = Polynomial(self.a, Basis.LAGRANGE)
        self.a_comm = self.setup.commit(self.a_poly)
        self.indexes = []
        for w in witness:
            self.indexes.append(self.table.get_index(w))

        self.dim_poly = []
        for i in range(self.c):
            values = [0 for _ in range(2**self.logm)]
            for j in range(self.m):
                values[j] = self.indexes[j][i]
            self.dim_poly.append(Polynomial(list(map(values, Scalar)), Basis.LAGRANGE))

        # Todo: use multivariate poly comm
        self.dim_comm = [self.setup.commit(poly) for poly in self.dim_poly]
        return Message1(self.a_comm, self.logm, self.dim_comm)
    
    def round_2(self):
        self.a_eval = self.setup.multivar_eval(self.a_poly, self.r)
        self.a_PIOP = self.setup.PIOP_prove(self.a, self.r, self.a_eval)

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

            # Todo: use multilinear poly
            self.E_poly.append(Polynomial(list(map(E_val, Scalar)), Basis.LAGRANGE))
            self.read_poly.append(Polynomial(list(map(read_ts, Scalar)), Basis.LAGRANGE))
            self.write_poly.append(Polynomial(list(map(write_cts, Scalar)), Basis.LAGRANGE))
            self.final_poly.append(Polynomial(list(map(final_cts, Scalar)), Basis.LAGRANGE))
        
        self.E_comm = [self.setup.commit(poly) for poly in self.E_poly]
        self.read_comm = [self.setup.commit(poly) for poly in self.read_poly]
        self.final_comm = [self.setup.commit(poly) for poly in self.final_poly]
        return Message2(self.a_eval, self.a_PIOP, self.E_comm, self.read_comm, self.final_comm)
    
    def round_3(self):
        # sum-check protocol on h(k) = eq(r, k) * g(...E_i(k))
        self.sumcheck_h_data = []
        k = [Scalar(0) for _ in range(self.logm)]
        for i in range(self.logm):
            # For each round, evaluate the ith param on 0, 1, 2
            # eq_mle * g has degree 2 for each variable, so we need three points
            # or evaluate at g, g^2, g^3, g^4=1?
            vals = []
            for j in range(3):
                k[0] = Scalar(j)
                # Calculate eq(r,k), E_i(k) for all possible k
                # Todo: maybe can reduce duplicate?
                val = Scalar(0)
                for bits in range(2**(self.logm-i-1)):
                    now = bits
                    for pos in range(i+1, self.logm):
                        k[pos] = Scalar(now%2)
                        now /= 2
                    E_eval = [self.setup.multivar_eval(E, k) for E in self.E_poly]
                    val += self.setup.eq_mle(self.r, k) * self.table.g_func(E_eval)
                vals.append(val)
            self.sumcheck_h_data.append(vals)
            k[i] = self.rz[i]
        self.E_eval = [self.setup.multivar_eval(E, self.rz) for E in self.E_poly]
        self.E_PIOP = [self.setup.PIOP_prove(e_poly, self.rz, e_eval) for (e_poly, e_eval) in zip(self.E_poly, self.E_eval)]
        return Message3(self.sumcheck_h_data, self.E_eval, self.E_PIOP)
    
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
            S_poly = self.grand_product_poly(S)
            RS_poly = self.grand_product_poly(RS)
            WS1_poly = self.grand_product_poly(WS1)
            WS2_poly = self.grand_product_poly(WS2)
            self.S_poly.append(S_poly)
            self.RS_poly.append(RS_poly)
            self.WS1_poly.append(WS1_poly)
            self.WS2_poly.append(WS2_poly)
            self.S_comm.append(self.setup.commit(S_poly))
            self.RS_comm.append(self.setup.commit(RS_poly))
            self.WS1_comm.append(self.setup.commit(WS1_poly))
            self.WS2_comm.append(self.setup.commit(WS2_poly))

        return Message4(self.S_comm, self.RS_comm, self.WS1_comm, self.WS2_comm)
    
    def round_5(self):
        self.sumcheck_S_data, self.sumcheck_RS_data, self.sumcheck_WS1_data, self.sumcheck_WS2_data = [], [], [], []
        self.S_data, self.RS_data, self.WS1_data, self.WS2_data = [], [], [], []
        self.E_eval2, self.dim_eval, self.read_eval, self.final_eval = [], [], [], []
        self.E_PIOP2, self.dim_PIOP, self.read_PIOP, self.final_PIOP = [], [], [], []
        for i in range(self.alpha):
            # Todo: 
            # For RS and WS1
            # run sumcheck protocol on g(k) = eq(r,k) * (fi(1,k)-fi(k,0)fi(k,1))
            # will get g(r') where len(r') = logm
            # warning! fi(1,k)-fi(k,0)fi(k,1) is not linear
            # For S and WS2
            # run sumcheck protocol with len(r') = l
            self.sumcheck_S_data.append([])
            self.sumcheck_RS_data.append([])
            self.sumcheck_WS1_data.append([])
            self.sumcheck_WS2_data.append([])
            self.S_data.append(GrandProductData(self.S_poly[i], self.r_prime2[i]))
            self.RS_data.append(GrandProductData(self.RS_poly[i], self.r_prime3[i]))
            self.WS1_data.append(GrandProductData(self.WS1_poly[i], self.r_prime4[i]))
            self.WS2_data.append(GrandProductData(self.WS2_poly[i], self.r_prime5[i]))
            self.E_eval2.append(self.setup.multivar_eval(self.E_poly[i], self.r_prime3[i]))
            self.dim_eval.append(self.setup.multivar_eval(self.dim_poly[i//self.k], self.r_prime3[i]))
            self.read_eval.append(self.setup.multivar_eval(self.read_poly[i], self.r_prime3[i]))
            self.final_eval.append(self.setup.multivar_eval(self.final_poly[i], self.r_prime2[i]))
            self.E_PIOP2.append(self.setup.PIOP_prove(self.E_poly[i], self.r_prime3[i], self.E_eval2[i]))
            self.dim_PIOP.append(self.setup.PIOP_prove(self.dim_poly[i//self.k], self.r_prime3[i], self.dim_eval[i]))
            self.read_PIOP.append(self.setup.PIOP_prove(self.read_poly[i], self.r_prime3[i], self.read_eval[i]))
            self.final_PIOP.append(self.setup.PIOP_prove(self.final_poly[i], self.r_prime2[i], self.final_eval[i]))

        return Message5(self.sumcheck_S_data, self.sumcheck_RS_data,
                        self.sumcheck_WS1_data, self.sumcheck_WS2_data,
                        self.S_data, self.RS_data, self.WS1_data, self.WS2_data,
                        self.E_eval2, self.dim_eval, self.read_eval, self.final_eval,
                        self.E_PIOP2, self.dim_PIOP, self.read_PIOP, self.final_PIOP)
    
    def grand_product_poly(self, multiset: list[tuple[Scalar, Scalar, Scalar]]):
        # see https://eprint.iacr.org/2020/1275.pdf, section 5
        f = [Hash(s, self.gamma, self.tau) for s in multiset] # f(0,x) = v(x)
        for i in range(len(f)-1):
            f.append(f[2*i]*f[2*i+1]) # f(1,x) = f(x,0) * f(x,1)
        f.append(Scalar(0))
        poly_f = Polynomial(f, Basis.LAGRANGE)
        return poly_f
    
    def generateGrandProductData(self, f: Polynomial, r: list[Scalar]):
        f_0_r = f.eval([Scalar(0)]+r)
        f_1_r = f.eval([Scalar(1)]+r)
        f_r_0 = f.eval(r+[Scalar(0)])
        f_r_1 = f.eval(r+[Scalar(1)])
        product = f.values[2**len(r)-2]
        f_0_r_PIOP = self.setup.PIOP_prove(f, r, f_0_r)
        f_1_r_PIOP = self.setup.PIOP_prove(f, r, f_1_r)
        f_r_0_PIOP = self.setup.PIOP_prove(f, r, f_r_0)
        f_r_1_PIOP = self.setup.PIOP_prove(f, r, f_r_1)
        product_PIOP = self.setup.PIOP_prove(f, r, product)
        return GrandProductData(f_0_r, f_1_r, f_r_0, f_r_1, product,
                f_0_r_PIOP, f_1_r_PIOP, f_r_0_PIOP, f_r_1_PIOP, product_PIOP)
