import random
from src.baloo.setup import Setup
from src.baloo.program import CommonPreprocessedInput
from src.baloo.prover import Prover

# setup: public setup includes srs
# public_table: public table
# lookup: values to lookup
def prover(setup: Setup, public_table: [any], lookup: [any]):
    print("Beginning prover test")

    print("table: ", public_table)
    print("lookup: ", lookup)
    group_order_n = len(lookup)

    prover = Prover(setup, public_table)
    proof = prover.prove(lookup)
    print("Prover test success")

    return proof

def verifier(setup: Setup, proof, m: int):
    print("Beginning verifier test")
    vk = setup.verification_key()
    assert vk.verify_proof(proof, setup, m)
    print("Verifier test success")

def simple_test():
    print("===========> Beginning simple test ===========> ")
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100
    # public table
    table = [1, 2, 3, 4, 5, 6, 7, 8]
    # values to lookup
    lookup = [3, 7, 3, 4]

    # number of powers of tau
    powers = len(table) * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof = prover(setup, table, lookup)
    # run verifier
    m = len(lookup)
    verifier(setup, proof, m)
    print("===========> End simple test ===========> ")

def random_test():
    print("===========> Beginning random test ===========> ")
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100

    # public table
    # table = [1...32]
    table_len = 32
    table = []
    for i in range(1, table_len + 1):
        table.append(i)
    print("table: ", table)

    # values to lookup
    lookup = []
    for _ in range(table_len):
        lookup.append(random.randint(1, table_len))

    # number of powers of tau
    powers = len(table) * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof = prover(setup, table, lookup)
    # run verifier
    m = len(lookup)
    verifier(setup, proof, m)


if __name__ == "__main__":
    simple_test()
    random_test()