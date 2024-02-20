import numpy as np
import py_ecc.bn128 as b
from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar

def poly_test():
    vals = [1, 2, 3, 4]
    vals_scalar = [Scalar(int(x)) for x in vals]
    roots_of_unity = Scalar.roots_of_unity(4)

    poly_lag = Polynomial(vals_scalar, Basis.LAGRANGE)
    poly_coeff = poly_lag.ifft()
    points = roots_of_unity + [Scalar(2), Scalar(3), Scalar(4)]
    for i in range(len(points)):
      point = points[i]
      eval_lag = poly_lag.barycentric_eval(point)
      coeff_eval = poly_coeff.coeff_eval(point)
      assert eval_lag == coeff_eval

    quo = poly_coeff / Scalar(2)
    print("quo: ", quo.values)
if __name__ == "__main__":
    print("===========> Beginning test <===========")
    poly_test()
    print("===========> Test success <===========")
