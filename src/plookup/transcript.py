from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point
from src.common_util.transcript import CommonTranscript

@dataclass
class Message1:
    f_comm: G1Point   # kzg commitment of f
    h1_comm: G1Point  # kzg commitment of h1
    h2_comm: G1Point  # kzg commitment of h2

@dataclass
class Message2:
    z_comm: G1Point   # kzg commitment of z

@dataclass
class Message3:
    # we want to prove agg(x)=0 for all x in H
    # q = agg(x)/vanish_polynomial
    q_comm: G1Point   # kzg commitment of q

@dataclass
class Message4:
    f_eval: Scalar    # f(zeta)
    h1_eval: Scalar   # h1(zeta)
    h2_eval: Scalar   # h2(zeta)
    z_eval: Scalar    # z(zeta)
    h1_g_eval: Scalar # h1(g*zeta)
    h2_g_eval: Scalar # h2(g*zeta)
    z_g_eval: Scalar  # z(g*zeta)

@dataclass
class Message5:
    agg_witness_comm: G1Point
    agg_g_witness_comm: G1Point

class Transcript(CommonTranscript):
    def round_1(self, message: Message1) -> tuple[Scalar, Scalar]:
        self.append_point(b"f_comm", message.f_comm)
        self.append_point(b"h1_comm", message.h1_comm)
        self.append_point(b"h2_comm", message.h2_comm)
        beta = self.get_and_append_challenge(b"beta")
        gamma = self.get_and_append_challenge(b"gamma")
        return beta, gamma
    
    def round_2(self, message: Message2) -> tuple[Scalar]:
        self.append_point(b"z_comm", message.z_comm)
        self.delta = self.get_and_append_challenge(b"delta")
        return self.delta
    
    def round_3(self, message: Message3, order: int) -> tuple[Scalar]:
        self.append_point(b"q_comm", message.q_comm)
        self.zeta = self.get_and_append_point(b"zeta", order)
        return self.zeta
    
    def round_4(self, message: Message4) -> tuple[Scalar]:
        self.append_scalar(b"f_eval", message.f_eval)
        self.append_scalar(b"h1_eval", message.h1_eval)
        self.append_scalar(b"h2_eval", message.h2_eval)
        self.append_scalar(b"z_eval", message.z_eval)
        self.append_scalar(b"h1_g_eval", message.h1_g_eval)
        self.append_scalar(b"h2_g_eval", message.h2_g_eval)
        self.append_scalar(b"z_g_eval", message.z_g_eval)
        self.eps = self.get_and_append_challenge(b"eps")
        return self.eps
