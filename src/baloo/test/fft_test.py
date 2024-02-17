import numpy as np
import py_ecc.bn128 as b
from src.common_util.curve import Scalar
from src.cq.fft import fft, ifft, ec_fft, ec_ifft, next_power_of_2, is_power_of_two


def ec_fft_test():
    # public table
    arr = [1, 2, 3, 4, 5]
    n = len(arr)
    # padding with 0
    next_power = next_power_of_2(len(arr))
    padding_len = next_power - n
    ext_arr = np.pad(arr, (0, padding_len))
    ext_arr = [b.multiply(b.G1, x) for x in ext_arr]
    print("ec_fft_test ext_arr: ", ext_arr)
    orig_fft = ec_fft(ext_arr)
    print("ec_fft_test: fft result: ", orig_fft)
    orig_from_ifft = ec_ifft(orig_fft)
    print("ec_fft_test: ifft result: ", orig_from_ifft)

    for i in range(len(ext_arr)):
        assert orig_from_ifft[i] == ext_arr[i]



def fft_test():
    arr = [1, 2, 3, 4, 5]
    n = len(arr)
    # padding with 0
    next_power = next_power_of_2(len(arr))
    padding_len = next_power - n
    ext_arr = np.pad(arr, (0, padding_len))
    print("fft_test ext_arr: ", ext_arr)
    ext_arr = [Scalar(int(x)) for x in ext_arr]
    orig_fft = fft(ext_arr)
    print("fft_test: fft result: ", orig_fft)
    orig_from_ifft = ifft(orig_fft)
    print("fft_test: ifft result: ", orig_from_ifft)
    assert np.array_equal(ext_arr, orig_from_ifft)


def fft_test_should_raise_exception():
    # This is not a power of 2, so it should raise the exception
    test_values = [1, 2, 3]
    try:
        ext_arr = [Scalar(int(x)) for x in test_values]
        fft(ext_arr)
        return "Test Failed: No Exception Raised"
    except AssertionError as e:
        if str(e) == "fft: values length should be powers of 2":
            return f"Test Failed: Exception Raised - {e}"
        else:
            return f"Test Failed: Different Exception Raised - {e}"


def ec_fft_test_should_raise_exception():
    # This is not a power of 2, so it should raise the exception
    test_values = [1, 2, 3]
    try:
        ext_arr = [Scalar(int(x)) for x in test_values]
        ec_fft(ext_arr)
        return "Test Failed: No Exception Raised"
    except AssertionError as e:
        if str(e) == "ec_fft: values length should be powers of 2":
            return f"Test Failed: Exception Raised - {e}"
        else:
            return f"Test Failed: Different Exception Raised - {e}"


def next_power_of_2_test():
    assert next_power_of_2(0) == 1
    assert next_power_of_2(1) == 1
    assert next_power_of_2(2) == 2
    assert next_power_of_2(3) == 4
    assert next_power_of_2(4) == 4
    assert next_power_of_2(5) == 8
    assert next_power_of_2(6) == 8
    assert next_power_of_2(7) == 8
    assert next_power_of_2(8) == 8
    assert next_power_of_2(9) == 16

def is_power_of_two_test():
    assert is_power_of_two(0) == False
    assert is_power_of_two(1) == True
    assert is_power_of_two(2) == True
    assert is_power_of_two(3) == False
    assert is_power_of_two(4) == True
    assert is_power_of_two(5) == False
    assert is_power_of_two(6) == False
    assert is_power_of_two(7) == False
    assert is_power_of_two(8) == True


if __name__ == "__main__":
    print("===========> Beginning test <===========")
    is_power_of_two_test()
    next_power_of_2_test()
    fft_test()
    ec_fft_test()
    fft_test_should_raise_exception()
    ec_fft_test_should_raise_exception()
    print("===========> End test <===========")
