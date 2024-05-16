from src.common_util.curve import Scalar

# generate input {0, 1}^(bit_count)
def generate_binary(bit_count) -> list[list[Scalar]]:
    binary = []

    def genbin(n, bs=[]):
        if len(bs) == n:
            binary.append(bs)
        else:
            b_zero = bs + [Scalar.zero()]
            b_one = bs + [Scalar.one()]
            genbin(n, b_zero)
            genbin(n, b_one)

    genbin(bit_count)
    return binary