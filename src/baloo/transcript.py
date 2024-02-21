from src.common_util.curve import Scalar, G1Point, G2Point
from src.common_util.merlin.merlin_transcript import MerlinTranscript
from py_ecc.secp256k1.secp256k1 import bytes_to_int
from dataclasses import dataclass


@dataclass
class Message1:
    # Commitments in G1
    z_I_comm_2: G2Point
    v_comm_1: G1Point
    t_comm_1: G1Point

@dataclass
class Message2:
    # Commitments in G1
    D_comm_1: G1Point
    R_comm_1: G1Point
    Q_D_comm_1: G1Point
    E_comm_1: G1Point
    Q_E_comm_1: G1Point

@dataclass
class Message3:
    # Commitments in G1
    a_comm_1: G1Point

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
        self.append_point(b"v_comm_1", message.v_comm_1)
        self.append_point(b"t_comm_1", message.t_comm_1)

        alpha = self.get_and_append_challenge(b"alpha")
        beta = self.get_and_append_challenge(b"beta")

        return alpha, beta

    def round_2(self, message: Message2) -> tuple[Scalar, Scalar]:
        self.append_point(b"R_comm_1", message.R_comm_1)
        self.append_point(b"D_comm_1", message.D_comm_1)
        self.append_point(b"Q_D_comm_1", message.Q_D_comm_1)
        self.append_point(b"E_comm_1", message.E_comm_1)
        self.append_point(b"Q_E_comm_1", message.Q_E_comm_1)

        rho = self.get_and_append_challenge(b"rho")
        gamma = self.get_and_append_challenge(b"gamma")
        zeta = self.get_and_append_challenge(b"zeta")

        return rho, gamma, zeta
