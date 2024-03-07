from py_ecc.fields.optimized_field_elements import FQ as Field
import py_ecc.optimized_bn128 as b
from typing import NewType

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
