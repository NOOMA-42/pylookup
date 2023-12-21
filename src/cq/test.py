from setup import Setup
from program import CommonPreprocessedInput
from prover import Prover
import random

# setup: public setup includes srs
# public_table: public table
# witness: values to lookup
def prover(setup: Setup, public_table: [any], witness: [any]):
    print("Beginning prover test")

    print("table: ", public_table)
    print("witness: ", witness)
    group_order_N = len(public_table)
    group_order_n = len(witness)

    prover = Prover(setup, public_table, group_order_N, group_order_n)
    proof = prover.prove(witness)
    print("Prover test success")

    return proof, group_order_N, group_order_n

def prover_simple_array(setup: Setup, public_table):
    print("Beginning prover_simple_array test")
    # values to lookup
    witness = [1, 2, 1, 7, 3, 6, 3, 2]
    proof, group_order_N, group_order_n = prover(setup, public_table, witness)
    return proof, group_order_N, group_order_n

def prover_random_lookup(setup: Setup, public_table):
    print("Beginning prover_random_lookup test")
    # values to lookup
    witness = []
    for _ in range(32):
        witness.append(random.randint(1, 256))

    proof, group_order_N, group_order_n = prover(setup, public_table, witness)
    return proof, group_order_N, group_order_n

def verifier(setup, proof, group_order_N, group_order_n):
    print("Beginning verifier test")
    common_preprocessed_input = CommonPreprocessedInput(group_order_N, group_order_n)
    vk = setup.verification_key(common_preprocessed_input)
    assert vk.verify_proof(proof, setup)
    print("Verifier test success")

def simple_test():
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100
    # public table
    table = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    group_order_N = len(table)
    # number of powers of tau
    powers = group_order_N * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof, group_order_N, group_order_n = prover_simple_array(setup, table)
    # run verifier
    verifier(setup, proof, group_order_N, group_order_n)

def random_test():
    # random number, normally comes from MPC(Multi-Party Computation)
    tau = 100

    # public table
    # table = [1...256]
    table = []
    for i in range(1, 257):
        table.append(i)
    print("table: ", table)

    group_order_N = len(table)
    # number of powers of tau
    powers = group_order_N * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof, group_order_N, group_order_n = prover_random_lookup(setup, table)
    # run verifier
    verifier(setup, proof, group_order_N, group_order_n)


if __name__ == "__main__":
    simple_test()
    random_test()