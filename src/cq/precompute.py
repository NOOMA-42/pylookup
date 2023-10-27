import os
from common_util.kzg import Setup
from common_util.poly import Polynomial, Basis
from common_util.curve import Scalar
import py_ecc.bn128 as b

cwd = os.getcwd()

def gen(N, t):
    # compute SRS
    setup = Setup.from_file(cwd + "/src/cq/test/powersOfTau28_hez_final_11.ptau")
    # TODO: consider commit_G2 with a monomial basis, current implementation only supports Lagrange basis. That implementation need to handle condition such as X^n - 1 tho.
    z_v = b.add(setup.powers_of_x2[setup.n - 1], setup.powers_of_x2[0])
    t = Polynomial(
        [Scalar(t_i) for t_i in t],
        Basis.LAGRANGE,
    )
    t_commitment = setup.commit_G2(t)
    
    # FIXME move this to unit test
    assert b.is_on_curve(t_commitment, b.b2)

    
