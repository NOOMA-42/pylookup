from src.common_util.curve import Scalar
from py_ecc.fields.field_elements import FQ as Field
import py_ecc.bn128 as b
from src.common_util.curve import G1Point, Scalar

PRIMITIVE_ROOT = 5

class Scalar(Field):
    field_modulus = b.curve_order

    # Gets the first root of unity of a given group order
    @classmethod
    def root_of_unity(cls, group_order: int):
        assert (cls.field_modulus - 1) % group_order == 0
        return Scalar(PRIMITIVE_ROOT) ** ((cls.field_modulus - 1) // group_order)

    # Gets the full list of roots of unity of a given group order
    @classmethod
    def roots_of_unity(cls, group_order: int):
        o = [Scalar(1), cls.root_of_unity(group_order)]
        while len(o) < group_order:
            o.append(o[-1] * o[1])
        return o

def next_power_of_2(n):
    """
    Get the next power of 2, for example if n = 5, then it will return 8

    :param n: The number to be processed.
    """
    return 1 if n == 0 else 2**(n - 1).bit_length()


def is_power_of_two(n):
    """
    Check if a given number is a power of two.

    :param n: The number to be checked.
    :return: True if n is a power of two, False otherwise.
    """
    if n <= 0:
        return False
    else:
        return (n & (n - 1)) == 0

def fft(values: list[Scalar], inv=False):
    def _fft(vals, modulus, roots_of_unity):
        if len(vals) == 1:
            return vals
        L = _fft(vals[::2], modulus, roots_of_unity[::2])
        R = _fft(vals[1::2], modulus, roots_of_unity[::2])
        o = [0] * len(vals)
        for i, (x, y) in enumerate(zip(L, R)):
            y_times_root = y * roots_of_unity[i]
            o[i] = (x + y_times_root) % modulus
            o[i + len(L)] = (x - y_times_root) % modulus
        return o

    assert is_power_of_two(len(values)), "fft: values length should be powers of 2"
    roots = [x.n for x in Scalar.roots_of_unity(len(values))]
    o, nvals = Scalar.field_modulus, [x.n for x in values]
    if inv:
        # Inverse FFT
        invlen = Scalar(1) / len(values)
        reversed_roots = [roots[0]] + roots[1:][::-1]
        return [Scalar(x) * invlen for x in _fft(nvals, o, reversed_roots)]
    else:
        # Regular FFT
        return [Scalar(x) for x in _fft(nvals, o, roots)]


def ifft(values: list[Scalar]):
    return fft(values, True)


def ec_fft(values: list[G1Point], inv=False):
    def _fft(vals: list[G1Point], modulus, roots_of_unity):
        if len(vals) == 1:
            return vals
        L = _fft(vals[::2], modulus, roots_of_unity[::2])
        R = _fft(vals[1::2], modulus, roots_of_unity[::2])
        o = [0] * len(vals)
        for i, (x, y) in enumerate(zip(L, R)):
            y_times_root = b.multiply(y, roots_of_unity[i])
            o[i] = b.add(x, y_times_root)
            o[i + len(L)] = b.add(x, b.neg(y_times_root))
        return o

    assert is_power_of_two(
        len(values)), "ec_fft: values length should be powers of 2"
    roots = [x.n for x in Scalar.roots_of_unity(len(values))]
    o, nvals = Scalar.field_modulus, values
    if inv:
        # Inverse FFT
        invlen = (Scalar(1) / len(values)).n
        reversed_roots = [roots[0]] + roots[1:][::-1]
        return [b.multiply(x, invlen) for x in _fft(nvals, o, reversed_roots)]
    else:
        # Regular FFT
        return _fft(nvals, o, roots)


def ec_ifft(values: list[G1Point]):
    return ec_fft(values, True)
