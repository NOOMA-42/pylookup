from src.common_util.curve import Scalar

def length_expansion(l: list[Scalar], v: int):
    if len(l) == v:
        return l
    elif len(l) < v:
        k = [Scalar.zero()] * (v - len(l))
        return l + k
    else:
        raise IndexError
