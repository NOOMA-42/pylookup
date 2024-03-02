import unittest
import sys

if __name__ == "__main__":
    sys.path.append('/Users/paulyu/opensource_project/pylookup/src/')    

    # common util
    loader = unittest.TestLoader()
    runner = unittest.TextTestRunner(verbosity=3)
    
    # univariatesumcheck
    util_suite = loader.discover(start_dir='./src/common_util', pattern='univariatesumcheck_test.py')
    runner.run(util_suite)

    # lagrange
    #util_suite = loader.discover(start_dir='./src/common_util', pattern='lagrange_test.py')
    #runner.run(util_suite)
    
    # multilinear

