from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point
from src.common_util.transcript import CommonTranscript
from src.lasso.program import GrandProductData

@dataclass
class Message1:
    a_comm: G1Point
    logm: int
    l: int
    dim_comm: list[G1Point] # multivariate polynomial commitment of f

@dataclass
class Message2:
    a_eval: Scalar
    a_PIOP: G1Point
    E_comm: list[G1Point]
    read_ts_comm: list[G1Point]
    final_cts_comm: list[G1Point]

@dataclass
class Message3:
    sumcheck_h_data: list
    E_eval: list[Scalar]
    E_PIOP: list[G1Point]

@dataclass
class Message4:
    S_comm: G1Point
    RS_comm: G1Point
    WS1_comm: G1Point
    WS2_comm: G1Point

@dataclass
class Message5:
    sumcheck_S_data: list
    sumcheck_RS_data: list
    sumcheck_WS1_data: list
    sumcheck_WS2_data: list
    S_data: list[GrandProductData]
    RS_data: list[GrandProductData]
    WS1_data: list[GrandProductData]
    WS2_data: list[GrandProductData]
    E_eval2: list[Scalar]
    dim_eval: list[Scalar]
    read_ts_eval: list[Scalar]
    final_cts_eval: list[Scalar]
    E_PIOP2: list[G1Point]
    dim_PIOP: list[G1Point]
    read_ts_PIOP: list[G1Point]
    final_cts_PIOP: list[G1Point]

class Transcript(CommonTranscript):
    def round_1(self, message: Message1) -> list[Scalar]:
        self.append_point(b"a_comm", message.a_comm)
        self.c = len(message.dim_comm)
        for i in range(self.c):
            self.append_point(bytes("dim"+str(i)+"_comm"), message.dim_comm[i])
        self.logm = message.logm
        self.l = message.l
        self.r = [self.get_and_append_challenge(b"r") for _ in range(self.logm)]
        # self.r_prime = [[self.get_and_append_challenge(b"r_prime") for _ in range(self.l)] for __ in range(self.c)]
        return self.r
    
    def round_2(self, message: Message2) -> tuple[Scalar]:
        self.append_scalar(b"a_eval", message.a_eval)
        self.append_point(b"a_PIOP", message.a_PIOP)
        self.alpha = len(message.E_comm)
        for i in range(self.alpha):
            self.append_point(bytes("E"+str(i)+"_comm"), message.E_comm[i])
            self.append_point(bytes("read_ts"+str(i)+"_comm"), message.read_ts_comm[i])
            self.append_point(bytes("final_cts"+str(i)+"_comm"), message.final_cts_comm[i])
        self.rz = [self.get_and_append_challenge(b"rz") for _ in range(self.logm)]
        return self.rz
    
    def round_3(self, message: Message3) -> tuple[Scalar, Scalar]:
        for i in range(self.alpha):
            self.append_scalar(bytes("E"+str(i)+"_eval"), message.E_eval[i])
            self.append_point(bytes("E"+str(i)+"_PIOP"), message.E_PIOP[i])
        self.tau = self.get_and_append_challenge(b"tau")
        self.gamma = self.get_and_append_challenge(b"gamma")
        return self.tau, self.gamma

    def round_4(self, message: Message4) -> tuple[list[Scalar], list[Scalar]]:
        self.append_point(b"S_comm", message.S_comm)
        self.append_point(b"RS_comm", message.RS_comm)
        self.append_point(b"WS1_comm", message.WS1_comm)
        self.append_point(b"WS2_comm", message.WS2_comm)
        self.r_prime2 = [[self.get_and_append_challenge(b"r_prime2") for _ in range(self.l)] for __ in range(self.alpha)]
        self.r_prime3 = [[self.get_and_append_challenge(b"r_prime3") for _ in range(self.logm)] for __ in range(self.alpha)]
        self.r_prime4 = [[self.get_and_append_challenge(b"r_prime3") for _ in range(self.logm)] for __ in range(self.alpha)]
        self.r_prime5 = [[self.get_and_append_challenge(b"r_prime2") for _ in range(self.l)] for __ in range(self.alpha)]
        return self.r_prime2, self.r_prime3, self.r_prime4, self.r_prime5
