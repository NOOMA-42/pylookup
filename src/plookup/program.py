from src.common_util.curve import Scalar, G1Point, ec_lincomb
from src.common_util.poly import Polynomial

class Params:
    table: list[int]
    order: int
    roots: list[Scalar]

    def __init__(self, table: list[int]):
        self.table = sorted(table)
        self.order = len(table)
        self.roots = Scalar.roots_of_unity(self.order)

def aggregate(power: Scalar, items: list):
    assert(len(items) > 0)
    weight = power
    result = items[0]
    for item in items[1:]:
        result += item * weight
        weight *= power
    return result

def aggregate_poly(power: Scalar, polys: list[Polynomial]):
    assert(len(polys) > 0)
    weight = power
    result = polys[0]
    for poly in polys[1:]:
        result = result.force_add(poly * weight)
        weight *= power
    return result

def aggregate_comm(power: Scalar, comms: list[G1Point]) -> G1Point:
    assert(len(comms) > 0)
    weight = power
    result = comms[0]
    for comm in comms[1:]:
        result = ec_lincomb([(result, 1), (comm, weight)])
        weight *= power
    return result
