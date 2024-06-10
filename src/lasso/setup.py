import py_ecc.bn128 as b
from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point
from src.common_util.mle_poly import polynomial

@dataclass
class Setup(object):

    def __init__(self, powers: int, tau: int):
        self.generate_srs(powers, tau)

    @classmethod
    # tau: a random number whatever you choose
    def generate_srs(self, powers: int, tau: int):
        print("Start to generate structured reference string")

        print("Finished to generate structured reference string")
        return True

    @classmethod
    def commit(self, values: polynomial) -> G1Point:
        # Todo
        return b.G1
    
    @classmethod
    def multivar_eval(self, poly: polynomial, point: list[Scalar]) -> Scalar:
        return poly.eval(point)
    
    @classmethod
    def PIOP_prove(self, poly: polynomial, point: list[Scalar], eval: Scalar) -> G1Point:
        # Todo
        return b.G1
    
    @classmethod
    def PIOP_verify(self, commitment: G1Point, point: list[Scalar], eval: Scalar, proof: G1Point) -> bool:
        # Todo
        return True
