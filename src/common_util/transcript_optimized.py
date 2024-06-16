from py_ecc.secp256k1.secp256k1 import bytes_to_int

from src.caulk.util import g1_point_to_bytes, g2_point_to_bytes
from src.common_util.curve_optimized import Scalar, G1Point, G2Point
from src.common_util.merlin.merlin_transcript import MerlinTranscript


class CommonTranscript(MerlinTranscript):
    def append(self, label: bytes, item: bytes) -> None:
        self.append_message(label, item)

    def append_scalar(self, label: bytes, item: Scalar):
        self.append_message(label, item.n.to_bytes(32, "big"))

    def append_g1_point(self, label: bytes, item: G1Point):
        self.append_message(label, g1_point_to_bytes(item))

    def append_g2_point(self, label: bytes, item: G2Point):
        self.append_message(label, g2_point_to_bytes(item))

    def get_and_append_challenge(self, label: bytes) -> Scalar:
        while True:
            challenge_bytes = self.challenge_bytes(label, 255)
            f = Scalar(bytes_to_int(challenge_bytes))
            if f != Scalar.zero():  # Enforce challenge != 0
                self.append_scalar(label, f)
                return f

    def get_and_append_point(self, label: bytes, order: int) -> Scalar:
        while True:
            challenge_bytes = self.challenge_bytes(label, 255)
            f = Scalar(bytes_to_int(challenge_bytes))
            if f ** order != Scalar.one():  # Enforce point not a root of unity
                self.append_scalar(label, f)
                return f
