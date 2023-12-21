from src.common_util.curve import Scalar, G1Point
from src.common_util.merlin.merlin_transcript import MerlinTranscript
from py_ecc.secp256k1.secp256k1 import bytes_to_int
from dataclasses import dataclass


@dataclass
class Message1:
    # Commitments in G1
    m_comm_1: G1Point


@dataclass
class Message2:
    # Commitments in G1
    A_comm_1: G1Point
    Q_A_comm_1: G1Point
    f_comm_1: G1Point
    B_0_comm_1: G1Point
    Q_B_comm_1: G1Point
    P_comm_1: G1Point

@dataclass
class Message3:
    # Commitments in G1
    b_0_at_gamma: G1Point
    f_at_gamma: G1Point
    a_at_0: G1Point
    pi_gamma: G1Point
    a_0_comm_1: G1Point

# https://merlin.cool/
class Transcript(MerlinTranscript):
    def append(self, label: bytes, item: bytes) -> None:
        self.append_message(label, item)

    def append_scalar(self, label: bytes, item: Scalar):
        self.append_message(label, item.n.to_bytes(32, "big"))

    def append_point(self, label: bytes, item: G1Point):
        self.append_message(label, item[0].n.to_bytes(32, "big"))
        self.append_message(label, item[1].n.to_bytes(32, "big"))

    def get_and_append_challenge(self, label: bytes) -> Scalar:
        while True:
            challenge_bytes = self.challenge_bytes(label, 255)
            f = Scalar(bytes_to_int(challenge_bytes))
            if f != Scalar.zero():  # Enforce challenge != 0
                self.append(label, challenge_bytes)
                return f

    def round_1(self, message: Message1) -> tuple[Scalar]:
        self.append_point(b"m_comm_1", message.m_comm_1)
        beta = self.get_and_append_challenge(b"beta")

        return beta

    def round_2(self, message: Message2) -> tuple[Scalar, Scalar]:
        self.append_point(b"A_comm_1", message.A_comm_1)
        self.append_point(b"Q_A_comm_1", message.Q_A_comm_1)
        self.append_point(b"f_comm_1", message.f_comm_1)
        self.append_point(b"B_0_comm_1", message.B_0_comm_1)
        self.append_point(b"Q_B_comm_1", message.Q_B_comm_1)
        self.append_point(b"P_comm_1", message.P_comm_1)

        gamma = self.get_and_append_challenge(b"gamma")
        eta = self.get_and_append_challenge(b"eta")

        return gamma, eta
