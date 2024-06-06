from dataclasses import dataclass

from src.common_util.curve_optimized import Scalar, G1Point, G2Point
from src.common_util.transcript_optimized import CommonTranscript


@dataclass
class Message1:
    g1_C_I: G1Point
    g1_Z_I: G1Point
    g1_u: G1Point
    g2_H1: G2Point


@dataclass
class Message2:
    g1_H2: G1Point


@dataclass
class Message3:
    v1: Scalar
    v2: Scalar
    pi_1: G1Point
    pi_2: G1Point
    pi_3: G1Point


class Transcript(CommonTranscript):

    def round1(self, message: Message1) -> Scalar:
        self.append_g1_point(b"g1_C_I", message.g1_C_I)
        self.append_g1_point(b"g1_Z_I", message.g1_Z_I)
        self.append_g1_point(b"g1_u", message.g1_u)
        self.append_g2_point(b"g2_H1", message.g2_H1)
        return self.get_and_append_challenge(b"round1")

    def round2(self, message: Message2) -> Scalar:
        self.append_g1_point(b"g1_H2", message.g1_H2)
        return self.get_and_append_challenge(b"round2")
