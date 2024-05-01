from dataclasses import dataclass
from src.common_util.curve import Scalar, G1Point, ec_lincomb
from src.common_util.poly import Polynomial

class SOSTable:
    l: int
    c: int
    k: int
    alpha: int
    tables: list[list[int]]

    def __init__(self, l: int, c: int, k: int, tables: list[list[int]]):
        assert(len(tables) == k * c)
        self.l = l
        self.c = c
        self.k = k
        self.alpha = k * c
        self.tables = tables

    def g_func(self, r: list[int]):
        assert len(r) == self.c
        # the g function
        # here we have g(r_1, ..., r_c) = 1 + r_1 + r_2 * 2**l + ... + r_c * 2**(l*(c-1))
        # represent the table [1, 2, ..., 2**(l*c)]
        # need to override this function
        ret = 1
        mul = 1
        for i in range(self.c):
            ret += self.tables[i][r[i]] * mul
            mul *= (2**self.l)
        return ret
    
    def get_index(self, value: int):
        # need to override this function
        val = value - 1
        index = []
        for i in range(self.c):
            index.append(val % (2**self.l))
            val //= (2**self.l)
        return index

class Params:
    table: SOSTable
    order: int
    roots: list[Scalar]

    def __init__(self, table: SOSTable):
        self.table = table
        self.order = len(table)
        self.roots = Scalar.roots_of_unity(self.order)

@dataclass
class GrandProductData:
    f_0_r: Scalar
    f_1_r: Scalar
    f_r_0: Scalar
    f_r_1: Scalar
    product: Scalar
    f_0_r_PIOP: G1Point
    f_1_r_PIOP: G1Point
    f_r_0_PIOP: G1Point
    f_r_1_PIOP: G1Point
    product_PIOP: G1Point

def log_ceil(n: int):
    # return smallest int that >= log(n)
    return len(bin(n-1))-2

def Hash(element: tuple[Scalar, Scalar, Scalar], tau: Scalar, gamma: Scalar) -> Scalar:
    return element[0]*gamma*gamma + element[1]*gamma + element[2] - tau

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
