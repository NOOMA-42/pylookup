from src.common_util.poly import Polynomial, Basis
from src.common_util.univariatesumcheck import create_col_mapping
from src.cq.setup import Setup # TODO: refactor to common_util


""" 
(i) e (t, [D]2) − e ([1]1u2, [1]2) = e ([R]1, [x]2) + e ([Q2]1, [zI ]2) 
"""

def test():
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

    # Encode the matrix M to 
    M = [
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 1],
    ]
    col_mapping = create_col_mapping(M)
    print(col_mapping)


    # phi_x is encoding of "a" vector, aka the values be looked up
    phi_x = Polynomial(list(map(Scalar, [3, 2, 3])), Basis.LAGRANGE)
    # t_x is encoding of "t" vector, aka the unique value in a 
    t_x = Polynomial(list(map(Scalar, [2, 3])), Basis.LAGRANGE)

    z_I_x


if __name__ == "__main__":
    test()