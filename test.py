from src.cq.test import test as cq_test
import unittest

if __name__ == "__main__":
    loader = unittest.TestLoader()
    suite = loader.discover(start_dir='./common_util', pattern='*_test.py')    
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
    # cq_test()