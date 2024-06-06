import random
from src.plookup.setup import Setup
from src.plookup.program import Params
from src.plookup.prover import Prover, Proof
from src.plookup.verifier import Verifier


# setup: public setup includes srs
# public_table: public table
# witness: values to lookup
def prover(setup: Setup, params: Params, witness: list[int]):
    print("Beginning prover test")

    print("table: ", params.table)
    print("witness: ", witness)

    prover = Prover(setup, params)
    proof = prover.prove(witness)
    print("Prover test success")

    return proof

def prover_simple_array(setup: Setup, params: Params):
    print("Beginning prover_simple_array test")
    # values to lookup
    witness = [1, 1, 5, 5, 6, 6, 5] # twinkle twinkle little star
    proof = prover(setup, params, witness)
    return proof

def prover_random_lookup(setup: Setup, params: Params):
    print("Beginning prover_random_lookup test")
    # values to lookup
    witness = []
    for _ in range(32):
        witness.append(random.randint(1, 256))

    proof = prover(setup, params, witness)
    return proof

def verifier(setup: Setup, params: Params, proof: Proof):
    print("Beginning verifier test")
    verifier = Verifier(setup, params)
    assert verifier.verify(proof)
    print("Verifier test success")

def simple_test():
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100
    # public table
    table = [1, 2, 3, 4, 5, 6, 7, 8]

    group_order_N = len(table)
    # number of powers of tau
    powers = group_order_N * 3
    # do setup
    setup = Setup(powers, tau)
    # set public params
    params = Params(table)
    # run prover
    proof = prover_simple_array(setup, params)
    # run verifier
    verifier(setup, params, proof)

def random_test():
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100

    # public table
    # table = [1...256]
    table = []
    for i in range(1, 257):
        table.append(i)
    
    group_order_N = len(table)
    # number of powers of tau
    powers = group_order_N * 3
    # do setup
    setup = Setup(powers, tau)
    # set public params
    params = Params(table)
    # run prover
    proof = prover_random_lookup(setup, params)
    # run verifier
    verifier(setup, params, proof)


if __name__ == "__main__":
    simple_test()
    random_test()
