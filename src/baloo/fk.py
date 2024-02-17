import numpy as np
import py_ecc.bn128 as b
from src.common_util.curve import Scalar
from src.cq.fft import fft, ec_fft, ec_ifft, is_power_of_two

# https://eprint.iacr.org/2023/033
def fk(coeffs, powers_of_x):
    print("\n ***************** Start fk() ****************")
    assert len(coeffs) == len(powers_of_x), "length should be equal"
    n = len(coeffs)
    assert is_power_of_two(n), "length should be power of 2"
    # Get first column of circulant matrix in length of 2 * len(coeffs)
    # For example: coeffs is [1, 2, 3, 4]
    # The first column of circulant matrix should be: [4, 0, 0, 0, 0, 0, 2, 3]
    first_col = coeffs.copy()
    # first coefficient is unused, so set it to Scalar(0)
    first_col[0] = Scalar(0)

    # get first column of circulant matrix in 2n size
    # 1. padding 0
    first_col = np.pad(first_col, (n, 0), 'constant', constant_values=(Scalar(0),))
    # 2. roll by 1 to right
    first_col = np.roll(first_col, 1)

    # inverse srs: delete last one then inverse
    inv_powers_of_x = powers_of_x[:-1][::-1]
    inv_powers_of_x.append(b.Z1)
    # padding n 0s to the end
    ec_neutral_vals = [b.Z1] * n
    padded_x = inv_powers_of_x + ec_neutral_vals

    # We have circulant matrix C, C = F_inv * diag(F * first_col) * F
    # F: DFT matrix, F_inv: inverse DFT matrix
    # We want to get Q_T_comm_poly_coeffs = C * x = F_inv * diag(F * first_col) * F * x
    # 1. right hand side: F * x
    rhs = ec_fft(padded_x)

    # 2. middle hand side: F * first_col
    mhs = fft(first_col)

    # middle * right (element wise) to get diagonal: diag(F * first_col) * F * x
    m_r_hs = [b.multiply(rhs[i], mhs[i].n) for i in range(len(rhs))]

    # 3. ifft
    result = ec_ifft(m_r_hs)

    # 4. return firt n values
    Q_comm_poly_coeffs = result[:n]
    print("\n ***************** End fk() ****************")
    return Q_comm_poly_coeffs

