import unittest
import sys

if __name__ == "__main__":
    sys.path.append('/Users/paulyu/opensource_project/pylookup/src/')    

    # common util
    loader = unittest.TestLoader()
    # lagrange
    # util_suite = loader.discover(start_dir='./src/common_util', pattern='lagrange_test.py')
    # lagrange
    util_suite = loader.discover(start_dir='./src/common_util', pattern='univariatesumcheck_test.py')
    runner = unittest.TextTestRunner(verbosity=3)
    runner.run(util_suite)

    # protocol    
