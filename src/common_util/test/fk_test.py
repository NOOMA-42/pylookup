import random
import py_ecc.bn128 as b
from collections import Counter
from src.cq.setup import Setup
from src.common_util.poly import Polynomial, Basis
from src.common_util.curve import Scalar
from src.common_util.fk import fk
import time


def fk_test(table, powers_of_x):
    print("\n ***************** start fk_test ****************")
    t_values = [Scalar(val) for val in table]
    # get T(X) with lagrange interpolation
    T_poly_lag = Polynomial(t_values, Basis.LAGRANGE)
    # get T(X) in coefficient form with ifft
    T_poly = T_poly_lag.ifft()
    t_poly_coeffs = T_poly.values

    # compute h values with fk
    return fk(t_poly_coeffs, powers_of_x)


def compute_commitment_with_fk(table, witness, Q_T_comm_poly_coeffs):
    print("\n ***************** start compute_commitment_with_fk ****************")
    table_len = len(table)
    t_values = [Scalar(val) for val in table]
    duplicates = dict(Counter(witness))
    m_values = [Scalar(duplicates.get(val, 0)) for val in table]

    beta = Scalar(5)
    A_values = []
    for i, t_i in enumerate(t_values):
        A_i = m_values[i]/(beta + t_i)
        A_values.append(A_i)
        # sanity check
        assert A_i == m_values[i]/(beta + t_i), "A: not equal"

    roots = Scalar.roots_of_unity(len(A_values))

    # Compute Quotient polynomial commitment of A(X)
    # Q_A_Comm_Accu: Accumulator for computing quotient polynomial commitment of A(X)
    Q_A_Comm = b.Z1
    for i in range(table_len):
        K_T_Comm = b.Z1
        root = Scalar(1)
        for j in range(table_len):
            K_T_Comm = b.add(K_T_Comm, b.multiply(Q_T_comm_poly_coeffs[j], root.n))
            root = root * roots[i]
        A_val = A_values[i].n
        scale = roots[i]/table_len
        # Compute Quotient polynomial commitment of T(X)
        Q_T_Comm = b.multiply(K_T_Comm, scale.n)
        A_times_Q_T_Comm = b.multiply(Q_T_Comm, A_val)
        # Do the accumulation
        Q_A_Comm = b.add(Q_A_Comm, A_times_Q_T_Comm)

    print("\n Commitment of A(X) with FK: \n", Q_A_Comm)

    return Q_A_Comm



def compute_commitment_without_fk(table, setup, witness):
    print("\n ***************** start compute_commitment_without_fk ****************")
    beta = Scalar(5)
    t_values = [Scalar(val) for val in table]
    duplicates = dict(Counter(witness))
    m_values = [Scalar(duplicates.get(val, 0)) for val in table]

    A_values = []
    for i, t_i in enumerate(t_values):
        A_i = m_values[i]/(beta + t_i)
        A_values.append(A_i)
        # sanity check
        assert A_i == m_values[i]/(beta + t_i), "A: not equal"

    T_poly_lag = Polynomial(t_values, Basis.LAGRANGE)
    # T(X) in coefficient form
    T_poly = T_poly_lag.ifft()

    # vanishing polynomial: X^N - 1, N = group_order_N - 1
    group_order_N = len(table)
    ZV_array = [Scalar(-1)] + [Scalar(0)] * (group_order_N - 1) + [Scalar(1)]

    # vanishing polynomial in coefficient form
    ZV_poly = Polynomial(ZV_array, Basis.MONOMIAL)
    m_poly_lag = Polynomial(m_values, Basis.LAGRANGE)
    m_poly = m_poly_lag.ifft()
    # 2.c. Q_A(X) in coefficient form
    A_poly_lag = Polynomial(A_values, Basis.LAGRANGE)
    A_poly = A_poly_lag.ifft()
    Q_A_poly = (A_poly * (T_poly + beta) - m_poly) / ZV_poly
    Q_A_comm = setup.commit_g1(Q_A_poly)
    print("\n Commitment of A(X) without FK:  \n", Q_A_comm)

    return Q_A_comm


if __name__ == "__main__":
    # random number for setup SRS
    tau = 100
    table_len = 64
    table = [x for x in range(1, table_len + 1)]
    witness = []
    for _ in range(table_len + 1):
        witness.append(random.randint(1, table_len + 1))

    setup = Setup.execute(table_len * 2, tau, table)

    # only FK
    powers_of_x_in_table_len = setup.powers_of_x[:table_len]
    start_time0 = time.perf_counter()
    Q_T_comm_poly_coeffs = fk_test(table, powers_of_x_in_table_len)
    end_time0 = time.perf_counter()
    elapsed_time0 = end_time0 - start_time0
    print(f"fk only: {elapsed_time0} seconds")

    # compute commitment with FK
    start_time1 = time.perf_counter()
    comm_with_fk = compute_commitment_with_fk(table, witness, Q_T_comm_poly_coeffs)
    end_time1 = time.perf_counter()

    elapsed_time1 = end_time1 - start_time1
    print(f"with fk: {elapsed_time1} seconds")

    # compute commitment without FK
    start_time2 = time.perf_counter()
    comm_without_fk = compute_commitment_without_fk(table, setup, witness)
    end_time2 = time.perf_counter()
    elapsed_time2 = end_time2 - start_time2

    assert comm_with_fk == comm_without_fk

    print(f"fk only: {elapsed_time0} seconds")
    print(f"with fk: {elapsed_time1} seconds")
    print(f"without fk: {elapsed_time2} seconds")



