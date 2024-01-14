from py_ecc.secp256k1.secp256k1 import bytes_to_int

from src.common_util.curve_optimized import G1Point, Scalar
from src.common_util.merlin.keccak import SHA3_256
from src.common_util.kzg_optimized import KZGSetup
from src.common_util.poly_optimized import Polynomial


def hash_ec_points(*points: G1Point) -> Scalar:
    c = SHA3_256(b''.join(ec_point_to_bytes(pt) for pt in points))
    return Scalar(bytes_to_int(c))


def ec_point_to_bytes(pt: G1Point) -> bytes:
    return pt[0].n.to_bytes(32, "big") + pt[1].n.to_bytes(32, "big")


def kzg_open(setup: KZGSetup, poly: Polynomial, x: Scalar):
    v = poly.eval(x)
