import unittest
from src.common_util.sumcheck import prove_sumcheck, verify_sumcheck
from src.common_util.curve import Scalar
from src.common_util.mle_poly import (
    polynomial, term, monomial, 
    generate_binary, eval_ext, get_multi_ext, eval_expansion
)

def test_polynomial():
    a = Scalar(2)
    b = Scalar(3)
    c = Scalar(4)
    d = Scalar(5)
    e = Scalar(1)

    term1 = term(coeff=a, i=1, const=e)  # (2 * x_1 + 1)
    term2 = term(coeff=b, i=2, const=c)  # (3 * x_2 + 4)
    term3 = term(coeff=a, i=3, const=a)  # (2 * x_3 + 2)

    monomial1 = monomial(coeff=d, terms=[term1, term2, term3])  # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4) * (2 * x_3 + 2))
    poly = polynomial(terms=[monomial1]) # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4) * (2 * x_3 + 2)) + 6
    return poly

def test_polynomial_function(values):
    x_1, x_2 , x_3 = values
    return 5 * ((2 * x_1 + 1) * (3 * x_2 + 4) * (2 * x_3 + 2))

def test_polynomial2():
    a = Scalar(2)
    b = Scalar(3)
    c = Scalar(4)
    d = Scalar(5)
    e = Scalar(1)

    term1 = term(coeff=a, i=1, const=e)  # (2 * x_1 + 1)
    term2 = term(coeff=b, i=2, const=c)  # (3 * x_2 + 4)
    
    monomial1 = monomial(coeff=d, terms=[term1, term2])  # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4))
    # you can't have extra const term in this polynomial
    poly = polynomial(terms=[monomial1]) # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4))
    return poly

def test_polynomial_function2(values):
    x_1, x_2  = values
    return 5 * (2 * x_1 + 1) * (3 * x_2 + 4)

zero = Scalar.zero()
one = Scalar.one()

def mult_layer_zero(arr: list[Scalar]) -> Scalar:
    if len(arr) == 3:
        if arr == [zero, zero, one]:
            return one
        elif arr == [one, one, zero]:
            return one
        elif arr == [zero, one, zero]: # [0, 1, 0]
            return one
        else:
            return zero
    else:
        raise ValueError("Invalid input length")

def W0Func(bitstring):
  if bitstring == [zero]:
    return Scalar(36)
  elif bitstring == [one]:
    return Scalar(6)

def W1Func(bitstring):
  if bitstring == [zero, zero]:
    return Scalar(9)
  elif bitstring == [zero, one]:
    return Scalar(4)
  elif bitstring == [one, zero]:
    return Scalar(6)
  elif bitstring == [one, one]:
    return Scalar(1)

class TestSumcheck(unittest.TestCase):
    def test_prove_sumcheck(self):
        try:
            f = test_polynomial()
            prove_sumcheck(f, 3, 1)
        except Exception as e:
            self.fail(e)
    
    def test_t(self):
        def D_func(arr: list[Scalar]) -> Scalar:
            if arr == [Scalar.zero()]:
                return Scalar(36)
            elif arr == [Scalar.one()]:
                return Scalar(6)
            else:
                raise ValueError("Invalid input")
        print(get_multi_ext(D_func, 1))

    def test_verify_sumcheck(self):
        try:
            v = 2

            # calculate g with all x_i evaluated across all possible assignments
            claim = Scalar.zero()
            f = test_polynomial2()

            assignments = generate_binary(v)
            for assignment in assignments:
                # Note: test_polynomial_function2 
                multi_expansion = get_multi_ext(test_polynomial_function2, v) 
                claim += eval_expansion(multi_expansion, assignment)
            
            proof, r = prove_sumcheck(f, v, 0)
            self.assertTrue(verify_sumcheck(claim, proof, r, v), "Verification failed")
        except Exception as e:
            self.fail(e)

    def test_prove_sumcheck_start_idx(self):
        v = 2
        f = test_polynomial2()
        # start from x_3
        prove_sumcheck(f, v, 1)