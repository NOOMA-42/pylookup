import unittest
from src.common_util.mle_poly import (
    generate_binary, eval_ext, eval_expansion, get_multi_ext,
    term, monomial, polynomial
)
from src.common_util.curve import Scalar

def mult_layer_zero(arr: list[Scalar]) -> Scalar:
    zero = Scalar.zero()
    one = Scalar.one()
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

def test_polynomial():
    a = Scalar(2)
    b = Scalar(3)
    c = Scalar(4)
    d = Scalar(5)
    e = Scalar(1)

    term1 = term(coeff=a, i=1, const=e)  # (2 * x_1 + 1)
    term2 = term(coeff=b, i=2, const=c)  # (3 * x_2 + 4)
    
    monomial1 = monomial(coeff=d, terms=[term1, term2])  # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4))
    poly = polynomial(terms=[monomial1], c=Scalar(6)) # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4)) + 6
    return poly

def test_polynomial2():
    # 15 * ((3 * x_2 + 4)) + 6
    coeff1 = Scalar(3)
    const1 = Scalar(4)
    coeff_monomial = Scalar(15)
    const_poly = Scalar(6)

    # Create term for (3 * x_2 + 4)
    term1 = term(coeff=coeff1, i=2, const=const1)

    # Create monomial for 15 * (term1)
    monomial1 = monomial(coeff=coeff_monomial, terms=[term1])

    # Create polynomial with monomial1 and constant term 6
    poly = polynomial(terms=[monomial1], c=const_poly)

def test_polynomial_function(values):
    # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4)) + 6
    x_1, x_2 = values
    return 5 * ((2 * x_1 + 1) * (3 * x_2 + 4)) + 6

class TestMlePoly(unittest.TestCase):
    def test_generate_binary(self):
        bit_count = 3
        binary = generate_binary(bit_count)
        self.assertEqual(len(binary), 2 ** bit_count)
        for i in range(2 ** bit_count):
            self.assertEqual(len(binary[i]), bit_count)
    
    def test_gen_bin_to_func(self):
        f = mult_layer_zero
        w = generate_binary(3)[6] # [1, 1, 0]
        self.assertEqual(f(w), Scalar.one())
    
    def test_eval_i(self):
        poly = test_polynomial()
        poly.eval_i(Scalar(1), 1)

    def test_eval_univariate(self):
        # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4)) + 6
        poly = test_polynomial()
        self.assertEqual(poly.eval_univariate(Scalar(1)), 111)

    def test_apply(self):
        """  
        p1: 2 * ((3 * x_1 + 4)) + 1 * ((0 * x_1 + 5)) + 0
        p2: 6 * ((3 * x_1 + 4) * (1 * x_1 + 2)) + 3 * ((0 * x_1 + 5) * (1 * x_1 + 2)) + 0
        p3: 6 * ((3 * x_1 + 4) * (1 * x_1 + 2)) + 15 * ((1 * x_1 + 2)) + 0
        """
        p1 = polynomial([monomial(Scalar(2), [term(Scalar(3), 1, Scalar(4))]), monomial(Scalar(1), [term(Scalar(0), 1, Scalar(5))])])
        # Performing some operations that might result in non-simplified terms
        p2 = p1 * polynomial([monomial(Scalar(3), [term(Scalar(1), 1, Scalar(2))])])
        p3 = p2.apply_all()

    def test_get_expansion(self):
        # 5 * ((2 * x_1 + 1) * (3 * x_2 + 4)) + 6
        poly = test_polynomial()
        expansion = poly.get_expansion()

    def test_get_multi_ext(self):
        multi_ext: list[list[Scalar]] = get_multi_ext(test_polynomial_function, 2)
        result = [[30, 1, 1], [40, 1, 0], [15, 0, 1], [26, 0, 0]]
        scalar_result = [[Scalar(term[0]), Scalar(term[1]), Scalar(term[2])] for term in result]
        self.assertEqual(multi_ext, scalar_result)

    def test_eval_expansion(self):
        # Define the polynomial: 3 + 2(x_1) + 4(x_2)^2
        poly = [[Scalar(3), Scalar(0), Scalar(0)], 
                [Scalar(2), Scalar(1), Scalar(0)], 
                [Scalar(4), Scalar(0), Scalar(2)]]
        
        # Define the point at which to evaluate the polynomial: (x, y) = (2, 3)
        r = [Scalar(2), Scalar(3)]
        
        # Evaluate the polynomial at the point (2, 3)
        result = eval_expansion(poly, r)
        
        # Expected result: 3 + 2*2 + 4*3^2 = 3 + 4 + 36 = 43
        expected = Scalar(43)
        
        # Assert the result is as expected
        self.assertEqual(result, expected)
        
    def test_eval_ext(self):
        # Define a simple polynomial function f(x, y) = x + y
        def f(w):
            return w[0] + w[1]

        # Define the point at which to evaluate the polynomial: (x, y) = (2, 3)
        r = [Scalar(2), Scalar(3)]

        # Evaluate the polynomial using eval_ext
        eval_ext(f, r)

