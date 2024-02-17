import random
from src.baloo.setup import Setup
from src.baloo.program import CommonPreprocessedInput
from src.baloo.prover import Prover

# setup: public setup includes srs
# public_table: public table
# witness: values to lookup
def prover(setup: Setup, public_table: [any], witness: [any]):
    print("Beginning prover test")

    print("table: ", public_table)
    print("witness: ", witness)
    group_order_n = len(witness)

    prover = Prover(setup, public_table, group_order_n)
    proof = prover.prove(witness)
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
    table = [1, 2, 3, 4]
    # values to lookup
    witness = [1, 2, 1, 3]

    group_order_N = len(table)
    group_order_n = len(witness)
    # number of powers of tau
    powers = group_order_N * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof = prover(setup, table, witness)
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
    witness = []
    for _ in range(table_len):
        witness.append(random.randint(1, table_len))

    group_order_N = len(table)
    group_order_n = len(witness)
    # number of powers of tau
    powers = group_order_N * 2
    # do setup
    setup = Setup.execute(powers, tau, table)
    # run prover
    proof = prover(setup, table, witness)
    # run verifier
    verifier(setup, proof, group_order_N, group_order_n)


if __name__ == "__main__":
    simple_test()
    random_test()