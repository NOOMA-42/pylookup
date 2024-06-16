from src.caulk.multiple.caulk_multiple_prover import CaulkMultipleProver
from src.caulk.multiple.caulk_multiple_setup import CaulkMultipleSetup
from src.caulk.multiple.caulk_multiple_verifier import CaulkMultipleVerifier
from src.common_util.curve_optimized import Scalar

if __name__ == "__main__":
    setup = CaulkMultipleSetup.example_setup()
    prover = CaulkMultipleProver(setup)
    proof = prover.prove([Scalar(1), Scalar(3)])

    verifier = CaulkMultipleVerifier(setup)
    verifier.verify(proof)
