import unittest
from src.common_util.curve import Scalar

class TestScalar(unittest.TestCase):
    def test_pow(self):
        Scalar(3)
        assert Scalar(3).exp(2) == Scalar(9)
        assert Scalar(3) ** 2 == Scalar(9)
