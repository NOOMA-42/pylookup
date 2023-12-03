import unittest
import sys

if __name__ == "__main__":
    sys.path.append('/Users/paulyu/opensource_project/pylookup/src/')    

    # util
    loader = unittest.TestLoader()
    util_suite = loader.discover(start_dir='./src/common_util', pattern='*_test.py')
    runner = unittest.TextTestRunner(verbosity=3)
    runner.run(util_suite)

    # protocol    
