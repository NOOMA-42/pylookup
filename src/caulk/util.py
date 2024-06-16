from py_ecc.secp256k1.secp256k1 import bytes_to_int

from src.common_util.curve_optimized import G1Point, Scalar, G2Point
from src.common_util.kzg_optimized import KZGSetup
from src.common_util.merlin.keccak import SHA3_256
from src.common_util.poly_optimized import Polynomial, Basis


def hash_ec_points(*points: G1Point) -> Scalar:
    c = SHA3_256(b''.join(g1_point_to_bytes(pt) for pt in points))
    return Scalar(bytes_to_int(c))


def g1_point_to_bytes(pt: G1Point) -> bytes:
    return pt[0].n.to_bytes(32, "big") + pt[1].n.to_bytes(32, "big")


def g2_point_to_bytes(pt: G2Point) -> bytes:
    return pt[0].coeffs[0].to_bytes(32, "big") + pt[1].coeffs[0].to_bytes(32, "big")


def kzg_open(setup: KZGSetup, poly: Polynomial, x: Scalar):
    v = poly.eval(x)


def vanishing_poly(n: int) -> Polynomial:
    vals = [Scalar(-1)] + [Scalar(0)] * (n - 1) + [Scalar(1)]
    return Polynomial(vals, Basis.MONOMIAL)


def lagrange_polys(n: int):
    polys = []
    for i in range(n):
        poly = Polynomial([Scalar(0)] * i + [Scalar(1)] + [Scalar(0)] * (n - i - 1), Basis.LAGRANGE)
        polys.append(poly.ifft())

    return polys


def single_term_poly(degree: int):
    vals = [Scalar(0)] * degree + [Scalar(1)]
    return Polynomial(vals, Basis.MONOMIAL)
