import unittest
from .precompute import gen

class CqTest(unittest.TestCase):
    def test_precompute(self) -> None:        
        # t for table
        dummy_t = [1, 2, 3, 4]
        gen(1, dummy_t)