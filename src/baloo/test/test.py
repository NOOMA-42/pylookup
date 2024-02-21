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

    prover = Prover(setup, public_table, group_order_n)
    proof = prover.prove(lookup)
    print("Prover test success")

    return proof

def verifier(setup: Setup, proof, group_order_N, group_order_n):
    print("Beginning verifier test")
    common_preprocessed_input = CommonPreprocessedInput(group_order_N, group_order_n)
    vk = setup.verification_key(common_preprocessed_input)
    assert vk.verify_proof(proof, setup)
    print("Verifier test success")

def simple_test():
    print("===========> Beginning simple test ===========> ")
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100
    # public table
    table = [1, 2, 3, 4, 5, 6, 7, 8]
    # values to lookup
    lookup = [3, 7, 3, 4]

    group_order_N = len(table)
    group_order_n = len(lookup)
    # number of powers of tau
    powers = group_order_N * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof = prover(setup, table, lookup)
    # run verifier
    verifier(setup, proof, group_order_N, group_order_n)
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

    group_order_N = len(table)
    group_order_n = len(lookup)
    # number of powers of tau
    powers = group_order_N * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof = prover(setup, table, lookup)
    # run verifier
    verifier(setup, proof, group_order_N, group_order_n)


if __name__ == "__main__":
    simple_test()
    random_test()