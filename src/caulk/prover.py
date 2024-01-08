from dataclasses import dataclass
from setup import Setup
from src.common_util.poly import Polynomial, Basis


@dataclass
class PublicInput:


@dataclass
class Prove:
    pi_pederson: object
    pi_unity: object
    z_comm: object

@dataclass
class Witness:

    # blinders (uniformly random data)
    a: object
    s: object

class Prover:
    setup: Setup
    table: list
    poly_c: Polynomial

    def __init__(self, setup: Setup, table: list):
        self.setup = setup
        self.table = table
        poly_c = Polynomial(table, Basis.LAGRANGE)

    def prove_single(self):

