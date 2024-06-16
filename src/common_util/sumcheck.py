# Code modified from https://github.com/jeong0982/gkr
#mle sumcheck instead of binary
from src.common_util.mle_poly import polynomial, generate_binary, eval_univariate
from src.common_util.curve import Scalar
from src.common_util.util import *

def prove_sumcheck(g: polynomial, v: int, offset=0) -> tuple[list[list[Scalar]], list[Scalar]]:
    # TODO: think about removing offset
    """
    params:
    g: the polynomial to prove
    v: number of variables
    offset: index to start from e.g. x_5: evaluate g at x_5, x_6, ... x_v
    
    returns:
    proof: containing all the coefficients of each round
    r: randomness generated in each round
    
    NOTE: 
    1. the offset begins from 0, meaning index begins from x_2, because i + 2 + start, because we make x_1 a variable in the first round
    2. this notation follows Proof Argument and Zero Knowledge by Justin Thaler
    """
    proof: list[list[Scalar]] = [] # containing all the coefficients of each round
    r: list[Scalar] = [] # hash of the coefficients of each round as a randomness
    # first round
    # g1(X1)=∑(x2,⋯,xv)∈{0,1}^v g(X_1,x_2,⋯,x_v)    
    if v > 1:
        g_1 = polynomial([])
        assignments = generate_binary(v - 1)
        for assignment in assignments:
            g_1_sub = polynomial(g.terms[:], g.constant)
            
            # Loop through every bit of the assignment
            for i, x_i in enumerate(assignment):
                idx = i + 2 + offset # the offset begins from 0, meaning index begins from x_2, because i + 2 + start
                g_1_sub = g_1_sub.eval_i(x_i, idx)
            g_1 += g_1_sub
        proof.append(g_1.get_all_coefficients()) # TODO: not sure if the proof should contain the coefficient of the polynomial

        r_1 = Scalar(sum(list(map(lambda x : int(x), g_1.get_all_coefficients())))) # FIXME: sum in this line should be hash
        r.append(r_1)

    # 1 < j < v round
    for j in range(1, v - 1):
        g_j = polynomial(g.terms[:], g.constant)
        assignments = generate_binary(v - j - 1)
        for i, r_i in enumerate(r):
            idx = i + 1 + offset
            g_j = g_j.eval_i(r_i, idx)
        
        res_g_j = polynomial([])
        for assignment in assignments:
            g_j_sub = polynomial(g_j.terms[:], g_j.constant)
            for k, x_i in enumerate(assignment):
                idx = j + k + offset + 1
                g_j_sub = g_j_sub.eval_i(x_i, idx)
            res_g_j += g_j_sub
        proof.append(res_g_j.get_all_coefficients())

        r_n = Scalar(sum(list(map(lambda x : int(x), proof[len(proof) - 1]))))
        r.append(r_n)

    # last round
    g_v = polynomial(g.terms[:], g.constant)
    for i, r_i in enumerate(r):
        idx = i + 1 + offset
        g_v = g_v.eval_i(r_i, idx)
    proof.append(g_v.get_all_coefficients())

    if v == 1:
        r_v = Scalar(sum(list(map(lambda x : int(x), g.get_all_coefficients())))) # TODO is this correct way to treat univariate polynomial?
    else:
        r_v = Scalar(sum(list(map(lambda x : int(x), proof[len(proof) - 1])))) #FIXME sum in this line should be hash
    r.append(r_v)

    return proof, r

# TODO accommodate +1 -1 case 
def verify_sumcheck(claim: Scalar, proof: list[list[Scalar]], r: list[Scalar], v: int, config="DEFAULT", g=None, p_q_plus_one_dict=None) -> bool:
    """  
    params:
    v: number of variables
    config: "DEFAULT" or "FRACTIONAL_GKR"
    
    With DEFAULT config:
    g: 
        optional
        the polynomial to prove, verifier evaluate herself
    
    With FRACTIONAL_GKR config:
    p_q_plus_one_dict: 
        optional
        dict[str, Scalar]
        given as univariate polynomial, q'(x), p'(x)

    """
    bn = len(proof)
    # Univariate case
    if(v == 1 and (eval_univariate(proof[0], Scalar.zero()) + eval_univariate(proof[0], Scalar.one())) == claim):
        return True
    # round 1 to round v: g_j-1(r_j-1) ?= g_j(0) + g_j(1)
    expected = claim
    for i in range(bn):
        q_zero = eval_univariate(proof[i], Scalar.zero())
        q_one = eval_univariate(proof[i], Scalar.one())

        if q_zero + q_one != expected:
            return False
        if Scalar(sum(list(map(lambda x : int(x), proof[i])))) != r[i]:
            return False
        expected = eval_univariate(proof[i], r[i])
    
    # Final check: g_v(r_v) ?= g(r1, r2, ..., rv)
    if config == "DEFAULT":
        if g is None:
            raise ValueError("g must be provided in default sumcheck")
        g_result = polynomial(g.terms[:], g.constant)
        g_result_value = Scalar(1)
        for j, x in enumerate(r, 1):
            if j == len(r):
                g_result_value = g_result.eval_univariate(x)
            g_result: polynomial = g_result.eval_i(x, j)
        if g_result_value == eval_univariate(proof[bn - 1], r[bn - 1]):
            return True
    elif config == "FRACTIONAL_GKR":
        if p_q_plus_one_dict is None:
            raise ValueError("p_q_plus_one_dict must be provided in fractional gkr sumcheck")
        if p_q_plus_one_dict["p_k_plus_one_one"] is None or p_q_plus_one_dict["p_k_plus_one_zero"] is None or p_q_plus_one_dict["q_k_plus_one_one"] is None or p_q_plus_one_dict["q_k_plus_one_zero"] is None:
            raise ValueError("Invalid p_q_plus_one_dict")
        # TODO make this random linear combined
        f_r_k: Scalar = (p_q_plus_one_dict["p_k_plus_one_one"] * p_q_plus_one_dict["q_k_plus_one_zero"] + p_q_plus_one_dict["p_k_plus_one_zero"] * p_q_plus_one_dict["q_k_plus_one_one"]) + p_q_plus_one_dict["q_k_plus_one_one"] * p_q_plus_one_dict["q_k_plus_one_zero"]
        if claim == f_r_k:
            return True
    else:
        raise ValueError("Invalid config")
    return False
