from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point
from src.common_util.transcript import CommonTranscript
from src.lasso.program import GrandProductData

@dataclass
class Message1:
    a_comm: G1Point
    logm: int
    dim_comm: list[G1Point] # multivariate polynomial commitment of f

@dataclass
class Message2:
    a_eval: Scalar
    a_eval_proof: G1Point
    E_comm: list[G1Point]
    read_ts_comm: list[G1Point]
    final_cts_comm: list[G1Point]

@dataclass
class Message3:
    h_sumcheck_proof: list[list[Scalar]]
    rz: list[Scalar]
    E_eval: list[Scalar]
    E_eval_proof: list[G1Point]

@dataclass
class Message4:
    S_comm: list[G1Point]
    RS_comm: list[G1Point]
    WS1_comm: list[G1Point]
    WS2_comm: list[G1Point]

@dataclass
class Message5:
    S_sumcheck_proof: list[list[list[Scalar]]]
    RS_sumcheck_proof: list[list[list[Scalar]]]
    WS1_sumcheck_proof: list[list[list[Scalar]]]
    WS2_sumcheck_proof: list[list[list[Scalar]]]
    r_prime2: list[list[Scalar]]
    r_prime3: list[list[Scalar]]
    r_prime4: list[list[Scalar]]
    r_prime5: list[list[Scalar]]
    S_data: list[GrandProductData]
    RS_data: list[GrandProductData]
    WS1_data: list[GrandProductData]
    WS2_data: list[GrandProductData]
    E_eval2: list[Scalar]
    dim_eval: list[Scalar]
    read_ts_eval: list[Scalar]
    final_cts_eval: list[Scalar]
    E_eval2_proof: list[G1Point]
    dim_eval_proof: list[G1Point]
    read_ts_eval_proof: list[G1Point]
    final_cts_eval_proof: list[G1Point]

class Transcript(CommonTranscript):
    def round_1(self, message: Message1) -> list[Scalar]:
        self.append_point(b"a_comm", message.a_comm)
        self.c = len(message.dim_comm)
        for i in range(self.c):
            self.append_point(b"dim_comm", message.dim_comm[i])
        self.logm = message.logm
        self.r = [self.get_and_append_challenge(b"r") for _ in range(self.logm)]
        return self.r
    
    def round_2(self, message: Message2):
        self.append_scalar(b"a_eval", message.a_eval)
        self.append_point(b"a_eval_proof", message.a_eval_proof)
        self.alpha = len(message.E_comm)
        for i in range(self.alpha):
            self.append_point(b"E_comm", message.E_comm[i])
            self.append_point(b"read_ts_comm", message.read_ts_comm[i])
            self.append_point(b"final_cts_comm", message.final_cts_comm[i])
    
    def round_3(self, message: Message3) -> tuple[Scalar, Scalar]:
        for i in range(self.alpha):
            self.append_scalar(b"E_eval", message.E_eval[i])
            self.append_point(b"E_eval_proof", message.E_eval_proof[i])
        self.tau = self.get_and_append_challenge(b"tau")
        self.gamma = self.get_and_append_challenge(b"gamma")
        return self.tau, self.gamma

    def round_4(self, message: Message4):
        for i in range(self.alpha):
            self.append_point(b"S_comm", message.S_comm[i])
            self.append_point(b"RS_comm", message.RS_comm[i])
            self.append_point(b"WS1_comm", message.WS1_comm[i])
            self.append_point(b"WS2_comm", message.WS2_comm[i])
