import common_util.polycommit.kzg as kzg
import os
cwd = os.getcwd()

def gen(N, t):
    # compute SRS
    setup = kzg.Setup.from_file(cwd + "/src/cq/test/powersOfTau28_hez_final_11.ptau")
    setup.powers_of_x
    setup.powers_of_x2
    # Z_v = 