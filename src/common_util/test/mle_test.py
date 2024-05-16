import unittest
from src.common_util.mle import generate_binary

class TestMLE(unittest.TestCase):
    def test_generate_binary(self):
        bit_count = 3
        binary = generate_binary(bit_count)
        self.assertEqual(len(binary), 2 ** bit_count)
        for i in range(2 ** bit_count):
            self.assertEqual(len(binary[i]), bit_count)
    
