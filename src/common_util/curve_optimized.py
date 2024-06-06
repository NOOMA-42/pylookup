import random
from typing import NewType

import py_ecc.optimized_bn128 as b
from py_ecc.fields.optimized_field_elements import FQ as Field

PRIMITIVE_ROOT = 5
G1Point = NewType("G1Point", tuple[b.FQ, b.FQ, b.FQ])
G2Point = NewType("G2Point", tuple[b.FQ2, b.FQ2, b.FQ2])
G1: G1Point = b.G1
G2: G2Point = b.G2


class Scalar(Field):
    field_modulus = b.curve_order

    # Gets the first root of unity of a given group order
    @classmethod
    def root_of_unity(cls, group_order: int):
        assert (cls.field_modulus - 1) % group_order == 0
        return Scalar(PRIMITIVE_ROOT) ** ((cls.field_modulus - 1) // group_order)

    # Gets the full list of roots of unity of a given group order
    @classmethod
    def roots_of_unity(cls, group_order: int):
        o = [Scalar(1), cls.root_of_unity(group_order)]
        while len(o) < group_order:
            o.append(o[-1] * o[1])
        return o

    @classmethod
    def get_random(cls):
        value = random.randint(0, cls.field_modulus - 1)
        return Scalar(value)

    def inv(self):
        return Scalar(pow(self.n, -1, self.field_modulus))

    def __mul__(self, other):
        from src.common_util.poly_optimized import Polynomial
        if isinstance(other, Polynomial):
            return other * self
        return Scalar(super().__mul__(other))

    def to_g1(self):
        return ec_mul(G1, self)

    def to_g2(self):
        return ec_mul(G2, self)


Base = NewType("Base", b.FQ)


def ec_eq(p1, p2):
    return b.eq(p1, p2)


def ec_add(p1, p2):
    return b.add(p1, p2)


def ec_neg(pt):
    return b.neg(pt)


def ec_sub(p1, p2):
    return b.add(p1, ec_neg(p2))


def ec_pairing(p1: G2Point, p2: G1Point):
    return b.pairing(p1, p2)


def ec_mul(pt, coeff):
    if hasattr(coeff, "n"):
        coeff = coeff.n
    return b.multiply(pt, coeff % b.curve_order)


def ec_lincomb(pairs):
    if isinstance(pairs[0][0][0], b.FQ):
        result = b.Z1
    else:
        result = b.Z2

    for pt, coeff in pairs:
        result = b.add(result, ec_mul(pt, coeff))
    return result
