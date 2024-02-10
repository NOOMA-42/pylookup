from py_ecc.fields.field_elements import FQ as Field
import py_ecc.bn128 as b
from typing import NewType

PRIMITIVE_ROOT = 5
G1Point = NewType("G1Point", tuple[b.FQ, b.FQ])
G2Point = NewType("G2Point", tuple[b.FQ2, b.FQ2])


class Scalar(Field):
    field_modulus = b.curve_order

    # Gets the first root of unity of a given group order
    @classmethod
    def root_of_unity(cls, group_order: int):
        assert (cls.field_modulus - 1) % group_order == 0
        return Scalar(PRIMITIVE_ROOT) ** ((cls.field_modulus - 1) // group_order)

    """  
    Gets the full list of roots of unity of a given group order
    
    Returns: 
    [1, omega, omega^2, ..., omega^(group_order-1)]
    
    Note: 
    Some papers use [omega^1, omega^2, ..., omega^(group_order)], where omega^(group_order) = 1
    """
    @classmethod
    def roots_of_unity(cls, group_order: int):
        o = [Scalar(1), cls.root_of_unity(group_order)]
        while len(o) < group_order:
            o.append(o[-1] * o[1])
        return o

def left_rotate_root_of_unity(roots_of_unity: list[Scalar]) -> list[Scalar]:
    return roots_of_unity[1:] + roots_of_unity[:1]


Base = NewType("Base", b.FQ)


def ec_mul(pt, coeff):
    if hasattr(coeff, "n"):
        coeff = coeff.n
    return b.multiply(pt, coeff % b.curve_order)


def ec_lincomb(pairs):
    o = b.Z1
    for pt, coeff in pairs:
        o = b.add(o, ec_mul(pt, coeff))
    return o
